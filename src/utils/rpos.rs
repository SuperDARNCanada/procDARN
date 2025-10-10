use crate::error::ProcdarnError;
use crate::gridding::grid::GridError;
use crate::gridding::grid_table::RADIUS_EARTH;
use crate::utils::hdw::HdwInfo;
use geodesy::prelude::*;
use igrf::declination;
use std::f64::consts::PI;
use time::Date;

/// Normalize a vector.
fn norm_vector(v: &Coor3D) -> Coor3D {
    let len = v.hypot3(&Coor3D::origin());
    Coor3D::raw(v[0] / len, v[1] / len, v[2] / len)
}

/// Convert geocentric coordinates with radian angles to cartesian
pub fn geocentric_to_cartesian(coord: &Coor3D) -> Coor3D {
    let x = coord[2] * coord[1].cos() * coord[0].cos();
    let y = coord[2] * coord[1].cos() * coord[0].sin();
    let z = coord[2] * coord[1].sin();
    Coor3D::raw(x, y, z)
}

/// Converts from geodetic coordinates to geocentric spherical coordinates.
/// The distance of the point in km from the center of the Earth is also calculated.
/// The WGS84 model of Earth is used.
pub fn geodetic_to_geocentric(geodetic_coords: &Coor2D) -> Coor3D {
    let semi_major_axis: f64 = 6378.137;
    let flattening: f64 = 1.0 / 298.257223563;
    let semi_minor_axis: f64 = semi_major_axis * (1.0 - flattening);
    let second_eccentricity_squared: f64 =
        (semi_major_axis * semi_major_axis) / (semi_minor_axis * semi_minor_axis) - 1.0;

    let gclat = ((semi_minor_axis * semi_minor_axis) / (semi_major_axis * semi_major_axis)
        * geodetic_coords[1].tan())
    .atan();
    let mut gclon = geodetic_coords[0];
    if gclon.to_degrees() > 180.0 {
        gclon -= 360.0_f64.to_radians();
    }
    let rho =
        semi_major_axis / (1.0 + second_eccentricity_squared * gclat.sin() * gclat.sin()).sqrt();

    Coor3D::raw(gclon, gclat, rho)
}

/// Converts from geocentric coordinates to geodetic coordinates.
/// The distance of the point in km to the center of the Earth is also calculated.
/// The WGS84 model of Earth is used.
pub fn geocentric_to_geodetic(geocentric_coords: Coor2D) -> Coor3D {
    let semi_major_axis: f64 = 6378.137;
    let flattening: f64 = 1.0 / 298.257223563;
    let semi_minor_axis: f64 = semi_major_axis * (1.0 - flattening);
    let second_eccentricity_squared: f64 =
        (semi_major_axis * semi_major_axis) / (semi_minor_axis * semi_minor_axis) - 1.0;

    let gdlat = ((semi_major_axis * semi_major_axis) / (semi_minor_axis * semi_minor_axis)
        * geocentric_coords[1].tan())
    .atan();
    let gdlon = geocentric_coords[0];

    let rho = semi_major_axis
        / (1.0
            + second_eccentricity_squared
                * geocentric_coords[1].sin()
                * geocentric_coords[1].sin())
        .sqrt();

    Coor3D::raw(gdlon, gdlat, rho)
}

/// Convert a vector v from radar-to-range/beam cell into local south/east/vertical
/// (horizontal) coordinates at location loc in geocentric coordinates
fn cartesian_to_local(loc: &Coor3D, v: &Coor3D) -> Coor3D {
    // Rotate v about the z-axis by the longitude
    let sx = loc[0].cos() * v[0] + loc[0].sin() * v[1];
    let sy = -loc[0].sin() * v[0] + loc[0].cos() * v[1];
    let sz = v[2];

    // Calculate the colatitude
    let lax = PI / 2.0 - loc[1];

    // Rotate the vector about the east-axis by the colatitude
    let tx = lax.cos() * sx - lax.sin() * sz;
    let ty = sy;
    let tz = lax.sin() * sx + lax.cos() * sz;

    Coor3D::raw(tx, ty, tz)
}

/// Calculates the slant range to a range gate in km.
/// Called slant_range in cnvtcoord.c of RST
pub fn slant_range(
    first_range: i32,
    range_sep: i32,
    rx_rise: i32,
    range_edge: i32,
    range_gate: i32,
) -> f32 {
    // The next two lines truncate to integers, for some reason
    let lag_to_first_range = first_range * 20 / 3; // microseconds
    let sample_separation = range_sep * 20 / 3; // microseconds

    (lag_to_first_range - rx_rise + (range_gate * sample_separation) + range_edge) as f32 * 0.15
}

/// Adjusts a point in geodetic coordinates to account for the oblateness of the Earth.
/// Called geocnvrt in cnvtcoord.c of RST
fn geocnvrt(point: &Coor2D, xal: f64, xel: f64) -> Coor2D {
    let kxg = xel.cos() * xal.sin();
    let kyg = xel.cos() * xal.cos();
    let kzg = xel.sin();

    let point_gc = geodetic_to_geocentric(point);
    let del = point[1] - point_gc[1];

    let kxr = kxg;
    let kyr = kyg * del.cos() + kzg * del.sin();
    let kzr = -kyg * del.sin() + kzg * del.cos();

    let ral = kxr.atan2(kyr);
    let rel = (kzr / (kxr * kxr + kyr * kyr).sqrt()).atan();

    Coor2D::raw(ral, rel)
}

/// Calculate a destination point (lat, lon) from a start point, distance, and bearing in degrees
/// East of North using the Haversine formula.
/// Called fldpnt_sph in invmag.c of RST
fn fieldpoint_sphere(start: Coor3D, bearing: f64, range: f64) -> Coor2D {
    // start: lon, lat, alt
    let start_lon = start[0];
    let start_lat = start[1];
    let start_alt = start[2];

    // Solving spherical triangle
    let c_side = PI / 2.0 - start_lat;
    let a_angle: f64;
    if bearing > 180.0 {
        a_angle = (bearing - 360.0).to_radians();
    } else {
        a_angle = bearing.to_radians();
    }

    let b_side = range / start_alt;
    let mut arg = b_side.cos() * c_side.cos() + b_side.sin() * c_side.sin() * a_angle.cos();

    if arg <= -1.0 {
        arg = -1.0;
    } else if arg >= 1.0 {
        arg = 1.0;
    }

    let a_side = arg.acos();
    arg = (b_side.cos() - a_side.cos() * c_side.cos()) / (a_side.sin() * c_side.sin());

    if arg <= -1.0 {
        arg = -1.0;
    } else if arg >= 1.0 {
        arg = 1.0;
    }

    let mut b_angle = arg.acos();
    if a_angle < 0.0 {
        b_angle = -b_angle;
    }

    let end_lat = PI / 2.0 - a_side;
    let mut end_lon = start_lon + b_angle;
    if end_lon < 0.0 {
        end_lon += 2.0 * PI;
    } else if end_lon > 2.0 * PI {
        end_lon -= 2.0 * PI;
    }

    Coor2D::raw(end_lon, end_lat)
}

/// Uses the Haversine formula to calculate bearing from a start point to an end point,
/// assuming a spherical Earth. Inputs in degrees, output in radians.
/// Called fldpnt_azm in invmag.c of RST
fn fieldpoint_azimuth(start_lat: f64, start_lon: f64, end_lat: f64, end_lon: f64) -> f64 {
    let a_side = (90.0 - end_lat).to_radians();
    let c_side = (90.0 - start_lat).to_radians();
    let b_angle = (end_lon - start_lon).to_radians();

    let mut arg = a_side.cos() * c_side.cos() + a_side.sin() * c_side.sin() * b_angle.cos();
    let b_side = arg.acos();

    arg = (a_side.cos() - b_side.cos() * c_side.cos()) / (b_side.sin() * c_side.sin());
    let mut a_angle = arg.acos();

    if b_angle < 0.0 {
        a_angle = -a_angle;
    }

    let mut bearing = a_angle;
    if bearing.is_nan() {
        bearing = 0.0;
    }

    bearing
}

/// Calculates the geocentric coordinates of a point located `direction` from `radar_location`,
/// with `radar_location` given in geocentric coordinates [lon, lat, rho] and `direction` given
/// in local azimuth, elevation, and slant range.
/// Called fldpnt in cnvtcoord.c of RST.
fn fieldpoint(radar_location: &Coor3D, direction: &Coor3D) -> Coor3D {
    /* Convert from global spherical [lon, lat, rho] to global Cartesian [x, y, z]
     * (rx,ry,rz: Earth centered) */
    let sin_colat = (PI / 2.0 - radar_location[1]).sin();
    let rx = radar_location[2] * sin_colat * radar_location[0].cos();
    let ry = radar_location[2] * sin_colat * radar_location[0].sin();
    let rz = radar_location[2] * (PI / 2.0 - radar_location[1]).cos();

    /* Convert from local spherical (ral, rel, r) to local Cartesian
     * (sx,sy,sz: south,east,up) */
    let mut sx = -direction[2] * direction[1].cos() * direction[0].cos();
    let mut sy = direction[2] * direction[1].cos() * direction[0].sin();
    let mut sz = direction[2] * direction[1].sin();

    /* Convert from local Cartesian to global Cartesian */
    let mut tx =
        (PI / 2.0 - radar_location[1]).cos() * sx + (PI / 2.0 - radar_location[1]).sin() * sz;
    let mut ty = sy;
    let mut tz =
        -(PI / 2.0 - radar_location[1]).sin() * sx + (PI / 2.0 - radar_location[1]).cos() * sz;
    sx = radar_location[0].cos() * tx - radar_location[0].sin() * ty;
    sy = radar_location[0].sin() * tx + radar_location[0].cos() * ty;
    sz = tz;

    /* Find global Cartesian coordinates of new point by vector addition */
    tx = rx + sx;
    ty = ry + sy;
    tz = rz + sz;

    /* Convert from global Cartesian to global spherical */
    let frho = ((tx * tx) + (ty * ty) + (tz * tz)).sqrt();
    let flat = PI / 2.0 - (tz / frho).acos();
    let flon = {
        if (tx == 0.0) && (ty == 0.0) {
            0.0
        } else {
            ty.atan2(tx)
        }
    };

    Coor3D::raw(flon, flat, frho)
}

/// Calculate the geocentric coordinates of a radar field point using either the standard or
/// Chisham virtual height model.
/// Called fldpnth in cnvtcoord.c of RST
fn fieldpoint_height(
    point: Coor3D,
    bearing_off_boresight: f32,
    boresight_bearing: f32,
    height: f32,
    slant_range: f32,
    chisham: bool,
) -> Result<Coor3D, GridError> {
    let mut virtual_height: f32;
    if chisham {
        if slant_range < 787.5 {
            virtual_height =
                108.974 + 0.0191271 * slant_range + 6.68283e-5 * slant_range * slant_range;
        } else if slant_range < 2137.5 {
            virtual_height =
                384.416 - 0.17864 * slant_range + 1.81405e-4 * slant_range * slant_range;
        } else {
            virtual_height =
                1098.28 - 0.354557 * slant_range + 9.39961e-5 * slant_range * slant_range;
        }
        if slant_range < 115.0 {
            virtual_height = slant_range / 115.0 * 112.0;
        }
    } else {
        if height <= 150.0 {
            virtual_height = height;
        } else {
            if slant_range <= 600.0 {
                virtual_height = 115.0;
            } else if slant_range < 800.0 {
                virtual_height = (slant_range - 600.0) / 200.0 * (height - 115.0) + 115.0;
            } else {
                virtual_height = height;
            }
        }
        if slant_range < 150.0 {
            virtual_height = (slant_range / 150.0) * 115.0;
        }
    }

    // let ellipse = Ellipsoid::named("WGS84")?;
    // let radar_geo = ellipse.cartesian(&point); // geodetic
    let radar_geo = geodetic_to_geocentric(&Coor2D::raw(point[0], point[1]));

    let radar_radius = radar_geo[2]; // Radius of Earth beneath point
    let mut earth_rad_under_point = radar_radius; // Will update with calculations
    let mut point_rho;

    let mut point_sph = Coor3D::default();

    // This will prevent elevation angle from being NaN later on
    let range = if slant_range == 0.0 {
        0.1
    } else {
        slant_range as f64
    };

    let mut point_height = virtual_height + 1.0; // Initialize to make the below loop a do-while loop
    while (point_height - virtual_height).abs() > 0.5 {
        point_rho = earth_rad_under_point + virtual_height as f64;

        // Elevation angle relative to horizon [radians]
        let angle_above_horizon =
            ((point_rho * point_rho - radar_radius * radar_radius - range * range)
                / (2.0 * radar_radius * range))
                .asin();

        // Need to calculate actual elevation angle for 1.5-hop propagation when using Chisham model
        // for coning angle correction
        let xel: f64;
        if chisham && range > 2137.5 {
            let gamma = ((radar_radius * radar_radius + point_rho * point_rho - range * range)
                / (2.0 * radar_radius * point_rho))
                .acos();
            let beta = (radar_radius * (gamma / 3.0).sin() / (range / 3.0)).asin();
            xel = PI / 2.0 - beta - (gamma / 3.0);
        } else {
            xel = angle_above_horizon;
        }

        // Estimate the off-array-normal azimuth
        let off_boresight_rad = bearing_off_boresight.to_radians() as f64;
        let boresight_bearing_rad = boresight_bearing.to_radians() as f64;
        let tan_azimuth: f64;
        if off_boresight_rad.cos() * off_boresight_rad.cos() - xel.sin() * xel.sin() < 0.0 {
            tan_azimuth = 1e32;
        } else {
            tan_azimuth = (off_boresight_rad.sin() * off_boresight_rad.sin()
                / (off_boresight_rad.cos() * off_boresight_rad.cos() - xel.sin() * xel.sin()))
            .sqrt();
        }
        let azimuth: f64;
        if off_boresight_rad > 0.0 {
            azimuth = tan_azimuth.atan();
        } else {
            azimuth = -tan_azimuth.atan();
        }

        // Pointing azimuth in radians east of north
        let xal = azimuth + boresight_bearing_rad;

        // Adjust azimuth and elevation for oblateness of the Earth
        let adjusted_point = geocnvrt(&Coor2D::raw(point[0], point[1]), xal, xel);

        // Obtain the global spherical coordinates of the field point
        let point_sph_new = fieldpoint(
            &radar_geo,
            &Coor3D::raw(adjusted_point[0], angle_above_horizon, range),
        );
        if point_sph_new == point_sph {
            panic!("stagnation!!")
        } else {
            point_sph = point_sph_new;
        }

        // Recalculate the radius of the Earth beneath the field point
        let geodetic = geocentric_to_geodetic(Coor2D::raw(point_sph[0], point_sph[1]));
        earth_rad_under_point = geodetic[2];

        point_height = (point_sph[2] - earth_rad_under_point) as f32;
    }

    Ok(point_sph)
}

/// This function converts a gate/beam coordinate to geographic position. The height of the
/// transformation is given by height - if this value is less than 90 then it is assumed to be the
/// elevation angle from the radar. If center is not equal to zero, then the calculation is assumed
/// to be for the center of the cell, not the edge. The calculated values are returned in geocentric
/// coordinates.
/// Called RPosGeo in cnvtcoord.c of RST
fn rpos_geo(
    center: bool,
    beam_num: i32,
    range_gate: i32,
    hdw: &HdwInfo,
    first_range: f32,
    range_sep: f32,
    rx_rise_time: f32,
    altitude: f32,
    chisham: bool,
) -> Result<Coor4D, GridError> {
    let mut beam_edge: f32 = 0.0;
    let mut range_edge: f32 = 0.0;

    if !center {
        beam_edge = -0.5 * hdw.beam_separation;
        range_edge = -0.5 * range_sep * 20.0 / 3.0;
    }

    let rx_rise = match rx_rise_time {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise_time,
    };

    let offset = hdw.max_num_beams as f32 / 2.0 - 0.5;

    // Calculate deviation from boresight in degrees
    let psi = hdw.beam_separation * (beam_num as f32 - offset) + beam_edge + hdw.boresight_shift;

    // Calculate the slant range to the range gate in km
    let distance = slant_range(
        first_range as i32,
        range_sep as i32,
        rx_rise as i32,
        range_edge as i32,
        range_gate,
    );
    // If the input altitude is below 90, then it is actually an input elevation angle in degrees.
    // If so, we calculate the field point height
    let field_point_height = if altitude < 90.0 {
        -RADIUS_EARTH
            + ((RADIUS_EARTH * RADIUS_EARTH)
                + 2.0 * distance * RADIUS_EARTH * altitude.to_radians().sin()
                + distance * distance)
                .sqrt()
    } else {
        altitude
    };

    // Calculate the geocentric coordinates of the field point
    let result = fieldpoint_height(
        Coor3D::geo(
            hdw.latitude as f64,
            hdw.longitude as f64,
            hdw.altitude as f64,
        ),
        psi,
        hdw.boresight,
        field_point_height,
        distance,
        chisham,
    )?;

    Ok(Coor4D::raw(
        result[0],
        result[1],
        result[2],
        distance as f64,
    ))
}

pub fn rpos_range_beam_azimuth_elevation(
    beam: i32,
    range: i32,
    year: i32,
    hdw: &HdwInfo,
    first_range: f32,
    range_sep: f32,
    rx_rise: f32,
    altitude: f32,
    chisham: bool,
) -> Result<(f32, f32, f32), GridError> {
    let site_location_geo = Coor3D::geo(
        hdw.latitude as f64,
        hdw.longitude as f64,
        hdw.altitude as f64,
    );

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise,
    };

    // let ellipse = Ellipsoid::named("WGS84")?;

    // Convert center of range/beam cell to geocentric latitude/longitude/altitude
    let cell_geoc = rpos_geo(
        true,
        beam,
        range,
        hdw,
        first_range,
        range_sep,
        rx_rise_time,
        altitude,
        chisham,
    )?;
    // Convert range/beam position from geocentric coordinates to global Cartesian coordinates
    // let cell_cartesian = ellipse.cartesian(&cell_geoc);
    let slant_range = cell_geoc[3];
    let cell_wout_srng = Coor3D::raw(cell_geoc[0], cell_geoc[1], cell_geoc[2]);
    let cell_cartesian = geocentric_to_cartesian(&cell_wout_srng);

    // Convert radar geocentric coordinates to global Cartesian coordinates
    // let site_location_cartesian = ellipse.cartesian(&site_location_geo);
    let site_location_cartesian = geocentric_to_cartesian(&site_location_geo);

    // Calculate vector from site to center of range/beam cell
    let del = cell_cartesian - site_location_cartesian;

    // Normalize the vector
    let normed_del = norm_vector(&del);

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let local_del = cartesian_to_local(&cell_wout_srng, &normed_del);

    // Normalize the local horizontal vector
    let mut normed_local_del = norm_vector(&local_del);

    // Calculate the magnetic field vector in nT at the geocentric spherical range/beam position
    let igrf_field = declination(
        cell_geoc[1],
        cell_geoc[0],
        cell_geoc[2] as u32,
        Date::from_calendar_date(year, time::Month::January, 1)
            .map_err(|_| ProcdarnError::Timestamp(format!("bad year: {year}")))?,
    )?;

    // Convert from north/east/down coordinates to south/east/up
    let b_field = Coor3D::raw(-igrf_field.x, igrf_field.y, -igrf_field.z);

    // Normalize the magnetic field vector
    let normed_b = norm_vector(&b_field);

    // Calculate a new local vertical component such that the radar-to-range/beam vector becomes
    // orthogonal to the magnetic field at the range/beam position
    normed_local_del[2] =
        -(normed_b[0] * normed_local_del[0] + normed_b[1] * normed_local_del[1]) / normed_b[2];

    // Normalize the new radar-to-range/beam vector
    normed_local_del = norm_vector(&normed_local_del);

    // Calculate the azimuth and elevation angles of the orthogonal radar-to-range/beam vector
    let elevation = normed_local_del[2].atan2(
        (normed_local_del[0] * normed_local_del[0] + normed_local_del[1] * normed_local_del[1])
            .sqrt(),
    );
    let azimuth = normed_local_del[1].atan2(-normed_local_del[0]);

    Ok((azimuth as f32, elevation as f32, slant_range as f32))
}

pub fn rpos_inv_mag(
    beam: i32,
    range: i32,
    year: i32,
    hdw: &HdwInfo,
    first_range: f32,
    range_sep: f32,
    rx_rise: f32,
    altitude: f32,
    chisham: bool,
    _old_aacgm: bool,
) -> Result<(Coor2D, f32, f32), GridError> {
    let site_location_geod = Coor3D::geo(
        hdw.latitude as f64,
        hdw.longitude as f64,
        hdw.altitude as f64 + RADIUS_EARTH as f64,
    );

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise,
    };

    // let ellipse = Ellipsoid::named("WGS84")?;

    // Convert center of range/beam cell to geocentric latitude/longitude/altitude
    let cell_geoc = rpos_geo(
        true,
        beam,
        range,
        hdw,
        first_range,
        range_sep,
        rx_rise_time,
        altitude,
        chisham,
    )?;

    // Convert range/beam position from geocentric coordinates to global Cartesian coordinates
    // let cell_cartesian = ellipse.cartesian(&cell_geoc);
    let cell_wout_srng = Coor3D::raw(cell_geoc[0], cell_geoc[1], cell_geoc[2]);
    let slant_range = cell_geoc[3];
    let cell_cartesian = geocentric_to_cartesian(&cell_wout_srng);

    let site_location_geoc =
        geodetic_to_geocentric(&Coor2D::raw(site_location_geod[0], site_location_geod[1]));

    // Convert radar geocentric coordinates to global Cartesian coordinates
    // let site_location_cartesian = ellipse.cartesian(&site_location_geo);
    let site_location_cartesian = geocentric_to_cartesian(&Coor3D::raw(
        site_location_geoc[0],
        site_location_geoc[1],
        site_location_geoc[2],
    ));

    // Calculate vector from site to center of range/beam cell
    let del = cell_cartesian - site_location_cartesian;

    // Normalize the vector
    let normed_del = norm_vector(&del);

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let local_del = cartesian_to_local(&cell_wout_srng, &normed_del);

    // Normalize the local horizontal vector
    let mut normed_local_del = norm_vector(&local_del);

    // Calculate the magnetic field vector in nT at the geocentric spherical range/beam position
    let igrf_field = declination(
        cell_geoc[1],
        cell_geoc[0],
        cell_geoc[2] as u32,
        Date::from_calendar_date(year, time::Month::January, 1)
            .map_err(|_| ProcdarnError::Timestamp(format!("invalid year: {year}")))?,
    )?;

    // Convert from north/east/down coordinates to south/east/up
    let b_field = Coor3D::raw(-igrf_field.x, igrf_field.y, -igrf_field.z);

    // Normalize the magnetic field vector
    let normed_b = norm_vector(&b_field);

    // Calculate a new local vertical component such that the radar-to-range/beam vector becomes
    // orthogonal to the magnetic field at the range/beam position
    normed_local_del[2] =
        -(normed_b[0] * normed_local_del[0] + normed_b[1] * normed_local_del[1]) / normed_b[2];

    // Normalize the new radar-to-range/beam vector
    normed_local_del = norm_vector(&normed_local_del);

    // Calculate the azimuth angle of the orthogonal radar-to-range/beam vector
    let azimuth = normed_local_del[1].atan2(-normed_local_del[0]);

    // Get geodetic coordinates of cell location
    let cell_geod = geocentric_to_geodetic(Coor2D::raw(cell_geoc[0], cell_geoc[1]));

    // Calculate virtual height of range/beam position
    let virtual_height = cell_geoc[2] - cell_geod[2];

    // TODO: Accept old_aacgm option
    // Convert range/beam position from geocentric lat/lon at virtual height to AACGM magnetic
    // lat/lon
    let mut mag_lat = 0.0;
    let mut mag_lon = 0.0;
    let mut mag_rad = 0.0;
    unsafe {
        aacgmv2_rs::AACGM_v2_Convert(
            cell_geoc[1].to_degrees(),
            cell_geoc[0].to_degrees(),
            virtual_height,
            &mut mag_lat,
            &mut mag_lon,
            &mut mag_rad,
            0,
        );
    }

    // Calculate pointing direction lat/lon given distance and bearing from the radar position
    // at the field point radius
    let pointing_loc = fieldpoint_sphere(cell_wout_srng, azimuth, range_sep as f64);

    // TODO: Accept old_aacgm option
    // Convert pointing direction position from geocentric lat/lon at virtual height to AACGM
    // magnetic coordinates
    let mut pointing_mag_lat: f64 = 0.0;
    let mut pointing_mag_lon: f64 = 0.0;
    let mut pointing_mag_rad: f64 = 0.0;
    unsafe {
        aacgmv2_rs::AACGM_v2_Convert(
            pointing_loc[1].to_degrees(),
            pointing_loc[0].to_degrees(),
            virtual_height,
            &mut pointing_mag_lat,
            &mut pointing_mag_lon,
            &mut pointing_mag_rad,
            0,
        );
    }

    // Make sure pointing_mag_lon lies between +/- 180 degrees
    if pointing_mag_lon - mag_lon > 180.0 {
        pointing_mag_lon -= 360.0;
    } else if pointing_mag_lon - mag_lon < -180.0 {
        pointing_mag_lon += 360.0;
    }

    // Calculate bearing (azimuth) to pointing direction lat/lon from the radar position in magnetic
    // coordinates
    let azimuth = fieldpoint_azimuth(mag_lat, mag_lon, pointing_mag_lat, pointing_mag_lon);

    Ok((
        Coor2D::geo(mag_lat, mag_lon),
        azimuth as f32,
        slant_range as f32,
    ))
}
