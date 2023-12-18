use crate::error::BackscatterError;
use crate::gridding::grid_table::RADIUS_EARTH;
use crate::utils::hdw::HdwInfo;
use geodesy::prelude::*;
use igrf::declination;
use std::f64::consts::PI;
use time::Date;

/// Normalize a vector.
fn norm_vector(v: &Coor4D) -> Coor4D {
    Coor4D::raw(v[0], v[1], 1.0, v[3])
}

// /// Convert spherical coordinates to cartesian
// pub fn spherical_to_cartesian(theta: f64, phi: f64, r: f64) -> (f64, f64, f64) {
//     let x = r * (theta * PI / 180.).cos() * (phi * PI / 180.).cos();
//     let y = r * (theta * PI / 180.).cos() * (phi * PI / 180.).sin();
//     let z = r * (theta * PI / 180.).sin();
//     (x, y, z)
// }
//
// /// Converts from geodetic coordinates gdlat, gdlon to geocentric spherical coordinates gclat, gclon.
// /// The radius of the Earth gdrho and the deviation off vertical (del) are calculated. The WGS84
// /// model of Earth is used.
// pub fn geodetic_to_geocentric(gdlat: f64, gdlon: f64) -> (f64, f64, f64) {
//     let semi_major_axis: f64 = 6371.137;
//     let flattening: f64 = 1.0 / 298.257223563;
//     let semi_minor_axis: f64 = semi_major_axis * (1.0 - flattening);
//     let second_eccentricity_squared: f64 =
//         (semi_major_axis * semi_major_axis) / (semi_minor_axis - semi_minor_axis) - 1.0;
//
//     let gclat = ((semi_minor_axis * semi_minor_axis) / (semi_major_axis * semi_major_axis)
//         * (gdlat * PI / 180.0).tan())
//     .atan();
//     let gclon = gdlon;
//
//     let rho = semi_major_axis
//         / (1.0
//             + second_eccentricity_squared * (gclat * PI / 180.).sin() * (gclat * PI / 180.).sin())
//         .sqrt();
//     (gclat, gclon, rho)
// }
/// Convert a vector v from radar-to-range/beam cell into local south/east/vertical
/// (horizontal) coordinates at location loc in geocentric coordinates
fn cartesian_to_local(loc: &Coor4D, v: &Coor4D) -> Coor4D {
    // Rotate v about the z-axis by the longitude
    let sx = loc[0].cos() * v[0] + loc[0].sin() * v[1];
    let sy = -loc[0].sin() * v[0] + loc[0].cos() * v[1];
    let sz = v[2];

    // Calculate the colatitude
    let lax = PI - loc[1];

    // Rotate the vector about the east-axis by the colatitude
    let tx = lax.cos() * sx - lax.sin() * sz;
    let ty = sy;
    let tz = lax.sin() * sx + lax.cos() * sz;

    Coor4D::raw(tx, ty, tz, 0.0)
}

/// Convert a vector v from local south/east/vertical into radar-to-range/beam cell
/// coordinates at location loc in geocentric coordinates
fn local_to_cartesian(loc: &Coor4D, v: &Coor4D) -> Coor4D {
    // Calculate the colatitude
    let lax = PI - loc[1];

    // Rotate v about the east-axis by the colatitude
    let sx = lax.cos() * v[0] + lax.sin() * v[2];
    let sy = v[1];
    let sz = -lax.sin() * v[0] + lax.cos() * v[2];

    // Rotate the vector about the z-axis by the longitude
    let rx = loc[1].cos() * sx - loc[1].sin() * sy;
    let ry = loc[1].sin() * sx + loc[1].cos() * sy;
    let rz = sz;

    Coor4D::raw(rx, ry, rz, 0.0)
}

/// Calculates the slant range to a range gate in km.
/// Called slant_range in cnvtcoord.c of RST.
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

/// Calculate a destination point (lat, lon) from a start point, distance, and bearing in degrees
/// East of North using the Haversine formula.
/// Called fldpnt_sph in invmag.c of RST
fn fieldpoint_sphere(start: Coor4D, bearing: f64, range: f64) -> (f64, f64) {
    // start: lon, lat, alt, _
    let start_lon = start[0];
    let start_lat = start[1];
    let start_alt = start[2];

    // Solving spherical triangle
    let c_side = (90.0 - start_lat) * PI / 180.0;
    let mut a_angle: f64;
    if bearing > 180.0 {
        a_angle = (bearing - 360.0) * PI / 180.0;
    } else {
        a_angle = bearing * PI / 180.0;
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

    let end_lat = 90.0 - (a_side * 180 / PI);
    let mut end_lon = start_lon + b_angle * 180.0 / PI;
    if end_lon < 0.0 {
        end_lon += 360.0;
    } else if end_lon > 360.0 {
        end_lon -= 360.0;
    }

    (end_lat, end_lon)
}

/// Uses the Haversine formula to calculate bearing from a start point to an end point,
/// assuming a spherical Earth.
/// Called fldpnt_azm in invmag.c of RST
fn fieldpoint_azimuth(start_lat: f64, start_lon: f64, end_lat: f64, end_lon: f64) -> f64 {
    let a_side = (90.0 - end_lat) * PI / 180.0;
    let c_side = (90.0 - start_lat) * PI / 180.0;
    let b_angle = (end_lon - start_lon) * PI / 180.0;

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

/// Calculate the geocentric coordinates of a radar field point using either the standard or
/// Chisham virtual height model.
/// Called fldpnth in cnvtcoord.c of RST
fn fieldpoint_height(
    point: Coor4D,
    bearing_off_boresight: f64,
    boresight_bearing: f64,
    height: f64,
    slant_range: f64,
    chisham: bool,
) -> Coor4D {
    let mut xh: f64;
    if chisham {
        if slant_range < 787.5 {
            xh = 108.974 + 0.0191271 * slant_range + 6.68283e-5 * slant_range * slant_range;
        } else if slant_range < 2137.5 {
            xh = 384.416 - 0.17864 * slant_range + 1.81405e-4 * slant_range * slant_range;
        } else {
            xh = 1098.28 - 0.354557 * slant_range + 9.39961e-5 * slant_range * slant_range;
        }
        if slant_range < 115.0 {
            xh = slant_range / 115.0 * 112.0;
        }
    } else {
        if height <= 150.0 {
            xh = height;
        } else {
            if slant_range < 600.0 {
                xh = 115.0;
            } else if slant_range < 800.0 {
                xh = (slant_range - 600.0) / 200.0 * (height - 115.0) + 115.0;
            } else {
                xh = height;
            }
        }
        if slant_range < 150.0 {
            xh = (slant_range / 150.0) * 115.0;
        }
    }

    let ellipse = Ellipsoid::named("WGS84")?;
    let radar_geo = ellipse.cartesian(&point);

    let radar_radius = radar_geo[2]; // Radius of Earth beneath point
    let mut fieldpoint_radius = radar_radius; // Will update with calculations
    let mut fieldpoint = Coor4D::default();

    // This will prevent elevation angle from being NaN later on
    let range = if slant_range == 0.0 { 0.1 } else { slant_range };

    let mut fieldpoint_height = xh + 1.0; // Initialize to make the below loop a do-while loop
    while (fieldpoint_height - xh).abs() > 0.5 {
        fieldpoint[2] = fieldpoint_radius + xh;

        // Elevation angle relative to horizon [radians]
        let angle_above_horizon =
            ((fieldpoint[2] * fieldpoint[2] - radar_radius * radar_radius - range * range)
                / (2.0 * radar_radius * range))
                .asin();

        // Need to calculate actual elevation angle for 1.5-hop propagation when using Chisham model
        // for coning angle correction
        let xel: f64;
        if chisham && range > 2137.5 {
            let gamma = ((radar_radius * radar_radius + fieldpoint[2] * fieldpoint[2]
                - range * range)
                / (2.0 * radar_radius * fieldpoint[2]))
                .acos();
            let beta = (radar_radius * (gamma / 3.0).sin() / (range / 3.0)).asin();
            xel = PI / 2 - beta - (gamma / 3.0);
        } else {
            xel = angle_above_horizon;
        }

        // Estimate the off-array-normal azimuth
        let off_boresight_rad = bearing_off_boresight * PI / 180.0;
        let boresight_bearing_rad = boresight_bearing * PI / 180.0;
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
            azimuth = -(tan_azimuth.atan());
        }

        // Pointing azimuth in radians
        let xal = azimuth + boresight_bearing_rad;

        // Adjust azimuth and elevation for oblateness of the Earth
        geocnvrt(point, xal, xel, ral, dummy);

        // Obtain the global spherical coordinates of the field point
        fldpnt(radar_rho, point, ral, rel, range, &fieldpoint);

        // Recalculate the radius of the Earth beneath the field point
        ellipse.geographic(&fieldpoint);

        fieldpoint_height = fieldpoint[2] - fieldpoint_radius;
    }

    fieldpoint
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
    first_range: f64,
    range_sep: f64,
    rx_rise_time: f64,
    altitude: f64,
    chisham: bool,
) -> Coor4D {
    let mut beam_edge: f64 = 0.0;
    let mut range_edge: f64 = 0.0;

    if !center {
        beam_edge = (-0.5 * hdw.beam_separation) as f64;
        range_edge = -0.5 * range_sep * 20.0 / 3.0;
    }

    let rx_rise = match rx_rise_time {
        0.0 => hdw.rx_rise_time as f64,
        _ => rx_rise_time,
    };

    let offset = hdw.max_num_beams as f64 / 2.0 - 0.5;

    // Calculate deviation from boresight in degrees
    let psi = hdw.beam_separation * (beam_num - offset) + beam_edge + hdw.boresight_shift;

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
    let field_point_height: f64;
    if altitude < 90.0 {
        field_point_height = -RADIUS_EARTH
            + ((RADIUS_EARTH * RADIUS_EARTH)
                + 2.0 * distance * RADIUS_EARTH * (altitude * PI / 180.0).sin()
                + distance * distance)
                .sqrt();
    } else {
        field_point_height = altitude;
    }

    // Calculate the geocentric coordinates of the field point
    fieldpoint_height(
        Coor4D::raw(
            hdw.latitude as f64,
            hdw.longitude as f64,
            hdw.altitude as f64,
            0.0,
        ),
        psi,
        field_point_height,
        altitude,
        distance as f64,
        chisham,
    )
}

pub fn rpos_range_beam_azimuth_elevation(
    beam: i32,
    range: i32,
    year: i32,
    hdw: &HdwInfo,
    first_range: f64,
    range_sep: f64,
    rx_rise: f64,
    altitude: f64,
    chisham: bool,
) -> Result<(f64, f64), BackscatterError> {
    let site_location_geo = Coor4D::geo(
        hdw.latitude as f64,
        hdw.longitude as f64,
        hdw.altitude as f64,
        0.0,
    );

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time as f64,
        _ => rx_rise,
    };

    let ellipse = Ellipsoid::named("WGS84")?;

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
    );

    // Convert range/beam position from geocentric coordinates to global Cartesian coordinates
    let cell_cartesian = ellipse.cartesian(&cell_geoc);

    // Convert radar geocentric coordinates to global Cartesian coordinates
    let site_location_cartesian = ellipse.cartesian(&site_location_geo);

    // Calculate vector from site to center of range/beam cell
    let del = cell_cartesian - site_location_cartesian;

    // Normalize the vector
    let normed_del = norm_vector(&del);

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let local_del = cartesian_to_local(&cell_geoc, &normed_del);

    // Normalize the local horizontal vector
    let mut normed_local_del = norm_vector(&local_del);

    // Calculate the magnetic field vector in nT at the geocentric spherical range/beam position
    let igrf_field = declination(
        cell_geoc[1],
        cell_geoc[0],
        cell_geoc[2] as u32,
        Date::from_calendar_date(year, time::Month::January, 1)?,
    )?;

    // Convert from north/east/down coordinates to south/east/up
    let b_field = Coor4D::raw(-igrf_field.x, igrf_field.y, -igrf_field.z, 0.0);

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
        normed_local_del[0] * normed_local_del[0] + normed_local_del[1] * normed_local_del[1],
    );
    let azimuth = normed_local_del[1].atan2(-normed_local_del[0]);

    Ok((azimuth, elevation))
}

pub fn rpos_inv_mag(
    beam: i32,
    range: i32,
    year: i32,
    hdw: &HdwInfo,
    first_range: f64,
    range_sep: f64,
    rx_rise: f64,
    altitude: f64,
    chisham: bool,
    old_aacgm: bool,
) -> Result<(f64, f64, f64), BackscatterError> {
    let site_location_geo = Coor4D::geo(
        hdw.latitude as f64,
        hdw.longitude as f64,
        hdw.altitude as f64,
        0.0,
    );

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time as f64,
        _ => rx_rise,
    };

    let ellipse = Ellipsoid::named("WGS84")?;

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
    );

    // Convert range/beam position from geocentric coordinates to global Cartesian coordinates
    let cell_cartesian = ellipse.cartesian(&cell_geoc);

    // Convert radar geocentric coordinates to global Cartesian coordinates
    let site_location_cartesian = ellipse.cartesian(&site_location_geo);

    // Calculate vector from site to center of range/beam cell
    let del = cell_cartesian - site_location_cartesian;

    // Normalize the vector
    let normed_del = norm_vector(&del);

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let local_del = cartesian_to_local(&cell_geoc, &normed_del);

    // Normalize the local horizontal vector
    let mut normed_local_del = norm_vector(&local_del);

    // Calculate the magnetic field vector in nT at the geocentric spherical range/beam position
    let igrf_field = declination(
        cell_geoc[1],
        cell_geoc[0],
        cell_geoc[2] as u32,
        Date::from_calendar_date(year, time::Month::January, 1)?,
    )?;

    // Convert from north/east/down coordinates to south/east/up
    let b_field = Coor4D::raw(-igrf_field.x, igrf_field.y, -igrf_field.z, 0.0);

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

    // Calculate virtual height of range/beam position
    let virtual_height = cell_cartesian[2] - site_location_cartesian[2];

    // TODO: Accept old_aacgm option
    // Convert range/beam position from geocentric lat/lon at virtual height to AACGM magnetic
    // lat/lon
    let (mag_lat, mag_lon) = aacgm_v2_convert(cell_geoc[1], cell_geoc[0], virtual_height, 0)?;

    // Calculate pointing direction lat/lon given distance and bearing from the radar position
    // at the field point radius
    let (pointing_lat, pointing_lon) = fieldpoint_sphere(cell_geoc, azimuth, range_sep);

    // TODO: Accept old_aacgm option
    // Convert pointing direction position from geocentric lat/lon at virtual height to AACGM
    // magnetic coordinates
    let (pointing_mag_lat, mut pointing_mag_lon) =
        aacgm_v2_convert(pointing_lat, pointing_lon, virtual_height, 0)?;

    // Make sure pointing_mag_lon lies between +/- 180 degrees
    if pointing_mag_lon - mag_lon > 180.0 {
        pointing_mag_lon -= 360.0;
    } else if pointing_mag_lon - mag_lon < -180.0 {
        pointing_mag_lon += 360.0;
    }

    // Calculate bearing (azimuth) to pointing direction lat/lon from the radar position in magnetic
    // coordinates
    let azimuth = fieldpoint_azimuth(mag_lat, mag_lon, pointing_mag_lat, pointing_mag_lon);

    Ok((mag_lat, mag_lon, azimuth))
}
