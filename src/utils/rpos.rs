use crate::error::ProcdarnError;
use crate::gridding::grid::GridError;
use crate::gridding::grid_table::RADIUS_EARTH;
use crate::utils::hdw::HdwInfo;
use igrf::declination;
use std::f64::consts::PI;
use time::Date;
use crate::utils::coords::{CartesianCoords, GeocentricCoords, GeodeticCoords, LocalAngularCoords, MagneticCoords};

/// Calculates the slant range to a range gate in km.
/// Called slant_range in cnvtcoord.c of RST
pub fn slant_range(
    first_range: i32,
    range_sep: i32,
    rx_rise: f64,
    range_edge: f64,
    range_gate: i32,
) -> f64 {
    // The next two lines truncate to integers, for some reason
    let lag_to_first_range = (first_range * 20 / 3) as f64; // microseconds
    let sample_separation = (range_sep * 20 / 3) as f64; // microseconds

    (lag_to_first_range - rx_rise + ((range_gate - 1) as f64 * sample_separation) + range_edge) * 0.15
}



/// Calculate a destination point (lat, lon) from a start point, distance, and bearing in degrees
/// East of North using the Haversine formula.
/// Called fldpnt_sph in invmag.c of RST
fn fieldpoint_sphere(start: GeocentricCoords, bearing: f64, range: f64) -> GeocentricCoords {
    // Solving spherical triangle
    let c_side = PI / 2.0 - start.lat;
    let a_angle: f64;
    if bearing > 180.0 {
        a_angle = (bearing - 360.0).to_radians();
    } else {
        a_angle = bearing.to_radians();
    }

    let b_side = range / start.rad;
    let mut arg = b_side.cos() * c_side.cos() + b_side.sin() * c_side.sin() * a_angle.cos();

    arg = arg.max(-1.0).min(1.0);

    let a_side = arg.acos();
    arg = (b_side.cos() - a_side.cos() * c_side.cos()) / (a_side.sin() * c_side.sin());

    arg = arg.max(-1.0).min(1.0);

    let mut b_angle = arg.acos();
    if a_angle < 0.0 {
        b_angle = -b_angle;
    }

    let end_lat = PI / 2.0 - a_side;
    let mut end_lon = start.lon + b_angle;
    if end_lon < 0.0 {
        end_lon += 2.0 * PI;
    } else if end_lon > 2.0 * PI {
        end_lon -= 2.0 * PI;
    }

    GeocentricCoords::new(end_lat, end_lon, 0.0)
}

/// Uses the Haversine formula to calculate bearing from a start point to an end point,
/// assuming a spherical Earth. Inputs in degrees, output in radians.
/// Called fldpnt_azm in invmag.c of RST
fn fieldpoint_azimuth(start: &GeocentricCoords, end: &GeocentricCoords) -> f64 {
    let a_side = PI / 2.0 - end.lat;
    let c_side = PI / 2.0 - start.lat;
    let b_angle = end.lon - start.lon;

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
fn fieldpoint(radar_location: &GeocentricCoords, direction: &LocalAngularCoords) -> GeocentricCoords {
    /* Convert from global spherical [lon, lat, rho] to global Cartesian [x, y, z]
     * (rx,ry,rz: Earth centered) */
    let sin_colat = (PI / 2.0 - radar_location.lat).sin();
    let rx = radar_location.rad * sin_colat * radar_location.lon.cos();
    let ry = radar_location.rad * sin_colat * radar_location.lon.sin();
    let rz = radar_location.rad * (PI / 2.0 - radar_location.lat).cos();

    /* Convert from local spherical (ral, rel, r) to local Cartesian
     * (sx,sy,sz: south,east,up) */
    let mut sx = -direction.range * direction.el.cos() * direction.az.cos();
    let mut sy = direction.range * direction.el.cos() * direction.az.sin();
    let mut sz = direction.range * direction.el.sin();

    /* Convert from local Cartesian to global Cartesian */
    let mut tx =
        (PI / 2.0 - radar_location.lat).cos() * sx + (PI / 2.0 - radar_location.lat).sin() * sz;
    let mut ty = sy;
    let mut tz =
        -(PI / 2.0 - radar_location.lat).sin() * sx + (PI / 2.0 - radar_location.lat).cos() * sz;
    sx = radar_location.lon.cos() * tx - radar_location.lon.sin() * ty;
    sy = radar_location.lon.sin() * tx + radar_location.lon.cos() * ty;
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

    GeocentricCoords::new(flat, flon, frho)
}

/// Calculate the geocentric coordinates of a radar field point using either the standard or
/// Chisham virtual height model.
/// Called fldpnth in cnvtcoord.c of RST
fn fieldpoint_height(
    point: &GeodeticCoords,
    bearing_off_boresight: f64,
    boresight_bearing: f64,
    height: f64,
    slant_range: f64,
    chisham: bool,
) -> Result<GeocentricCoords, GridError> {
    let mut virtual_height: f64;
    if chisham {
        if slant_range < 787.5 {
            virtual_height =
                108.974 + 0.0191271 * slant_range + 6.68283e-5 * slant_range * slant_range;
        } else if slant_range <= 2137.5 {
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

    let radar_geo = point.to_geocentric();
    let radar_radius = radar_geo.rad; // Radius of Earth beneath point
    let mut earth_rad_under_point = radar_radius; // Will update with calculations
    let mut point_rho;

    let mut point_sph = GeocentricCoords::default();

    // This will prevent elevation angle from being NaN later on
    let range = if slant_range == 0.0 {
        0.1
    } else {
        slant_range
    };

    let mut point_height = virtual_height + 1.0; // Initialize to make the below loop a do-while loop
    while (point_height - virtual_height).abs() > 0.5 {
        point_rho = earth_rad_under_point + virtual_height;

        // Elevation angle relative to horizon [radians]
        let angle_above_horizon =
            ((point_rho * point_rho - radar_geo.rad * radar_geo.rad - range * range)
                / (2.0 * radar_geo.rad * range))
                .asin();

        // Need to calculate actual elevation angle for 1.5-hop propagation when using Chisham model
        // for coning angle correction
        let xel: f64;
        if chisham && range > 2137.5 {
            let gamma = ((radar_geo.rad * radar_geo.rad + point_rho * point_rho - range * range)
                / (2.0 * radar_geo.rad * point_rho))
                .acos();
            let beta = (radar_geo.rad * (gamma / 3.0).sin() / (range / 3.0)).asin();
            xel = PI / 2.0 - beta - (gamma / 3.0);
        } else {
            xel = angle_above_horizon;
        }

        // Estimate the off-array-normal azimuth
        let off_boresight_rad = bearing_off_boresight.to_radians();
        let boresight_bearing_rad = boresight_bearing.to_radians();
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
        let (ral, _) = point.correct_look_dir(xal, xel);

        // Obtain the global spherical coordinates of the field point
        let point_sph_new = fieldpoint(
            &radar_geo,
            &LocalAngularCoords::new(ral, angle_above_horizon, range),
        );
        if point_sph_new == point_sph {
            panic!("stagnation!!")
        } else {
            point_sph = point_sph_new;
        }

        // Recalculate the radius of the Earth beneath the field point
        let geodetic = point_sph.to_geodetic();
        earth_rad_under_point = geodetic.rad;

        point_height = point_sph.rad - earth_rad_under_point;
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
) -> Result<(GeocentricCoords, f64), GridError> {
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
        rx_rise as f64,
        range_edge as f64,
        range_gate + 1,
    ) as f32;
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
        &GeodeticCoords::new(
            hdw.latitude.to_radians() as f64,
            hdw.longitude.to_radians() as f64,
            hdw.altitude as f64,
        ),
        psi as f64,
        hdw.boresight as f64,
        field_point_height as f64,
        distance as f64,
        chisham,
    )?;

    Ok((result, distance as f64))
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
) -> Result<LocalAngularCoords, GridError> {
    let site_location_geod = GeodeticCoords::new(
        hdw.latitude.to_radians() as f64,
        hdw.longitude.to_radians() as f64,
        0.0
    );

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise,
    };

    // Convert center of range/beam cell to geocentric latitude/longitude/altitude
    let (cell_geoc, slant_range) = rpos_geo(
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
    let cell_cartesian = cell_geoc.to_cartesian();

    // Convert radar geocentric coordinates to global Cartesian coordinates
    // let site_location_cartesian = ellipse.cartesian(&site_location_geo);
    let site_location_geoc = site_location_geod.to_geocentric();
    let site_location_cartesian = site_location_geoc.to_cartesian();

    // Calculate vector from site to center of range/beam cell
    let mut del = cell_cartesian - site_location_cartesian;

    // Normalize the vector
    del.norm();

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let mut local_del = cell_geoc.cartesian_to_local(&del);

    // Normalize the local horizontal vector
    local_del.norm();

    // Calculate the magnetic field vector in nT at the geocentric spherical range/beam position
    let igrf_field = declination(
        cell_geoc.lat.to_degrees(),
        cell_geoc.lon.to_degrees(),
        cell_geoc.rad,
        Date::from_calendar_date(year, time::Month::January, 1)
            .map_err(|_| ProcdarnError::Timestamp(format!("bad year: {year}")))?,
    )?;

    // Convert from north/east/down coordinates to south/east/up
    let mut b_field = CartesianCoords { x: igrf_field.x, y: igrf_field.y, z: igrf_field.z };

    // Normalize the magnetic field vector
    b_field.norm();

    // Calculate a new local vertical component such that the radar-to-range/beam vector becomes
    // orthogonal to the magnetic field at the range/beam position
    local_del.up =
        -(b_field.x * local_del.south + b_field.y * local_del.east) / b_field.z;

    // Normalize the new radar-to-range/beam vector
    local_del.norm();

    // Calculate the azimuth and elevation angles of the orthogonal radar-to-range/beam vector
    let elevation = local_del.up.atan2(
        (local_del.south * local_del.south + local_del.east * local_del.east)
            .sqrt(),
    );
    let azimuth = local_del.east.atan2(-local_del.south);

    Ok(LocalAngularCoords::new(azimuth, elevation, slant_range))
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
) -> Result<(MagneticCoords, f32, f32), GridError> {
    let site_location_geod = GeodeticCoords::new(
        hdw.latitude.to_radians() as f64,
        hdw.longitude.to_radians() as f64,
        hdw.altitude as f64 + RADIUS_EARTH as f64,
    );

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise,
    };

    // Convert center of range/beam cell to geocentric latitude/longitude/altitude
    let (cell_geoc, slant_range) = rpos_geo(
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
    let cell_cartesian = cell_geoc.to_cartesian();

    let site_location_geoc = site_location_geod.to_geocentric();

    // Convert radar geocentric coordinates to global Cartesian coordinates
    let site_location_cartesian = site_location_geoc.to_cartesian();

    // Calculate vector from site to center of range/beam cell
    let mut del = cell_cartesian - site_location_cartesian;

    // Normalize the vector
    del.norm();

    // Convert the normalized vector from cartesian into local south/east/vertical coordinates
    let mut local_del = cell_geoc.cartesian_to_local(&del);

    // Normalize the local horizontal vector
    local_del.norm();

    // Calculate the magnetic field vector in nT at the geocentric spherical range/beam position
    let igrf_field = declination(
        cell_geoc.lat.to_degrees(),
        cell_geoc.lon.to_degrees(),
        cell_geoc.rad,
        Date::from_calendar_date(year, time::Month::January, 1)
            .map_err(|_| ProcdarnError::Timestamp(format!("invalid year: {year}")))?,
    )?;

    let mut b_field = CartesianCoords::new(igrf_field.x, igrf_field.y, igrf_field.z);

    // Normalize the magnetic field vector
    b_field.norm();

    // Calculate a new local vertical component such that the radar-to-range/beam vector becomes
    // orthogonal to the magnetic field at the range/beam position
    local_del.up =
        -(b_field.x * local_del.south + b_field.y * local_del.east) / b_field.z;

    // Normalize the new radar-to-range/beam vector
    local_del.norm();

    // Calculate the azimuth angle of the orthogonal radar-to-range/beam vector
    let azimuth = local_del.east.atan2(-local_del.south);

    // Get geodetic coordinates of cell location
    let cell_geod = cell_geoc.to_geodetic();

    // Calculate virtual height of range/beam position
    let virtual_height = cell_geoc.rad - cell_geod.rad;

    // TODO: Accept old_aacgm option
    // Convert range/beam position from geocentric lat/lon at virtual height to AACGM magnetic
    // lat/lon
    let mut geoc_with_virtual_height = cell_geoc.clone();
    geoc_with_virtual_height.rad = virtual_height;
    let mag_coords: GeocentricCoords;
    unsafe { mag_coords = geoc_with_virtual_height.aacgmv2_convert(); }

    // Calculate pointing direction lat/lon given distance and bearing from the radar position
    // at the field point radius
    let mut pointing_loc = fieldpoint_sphere(cell_geoc, azimuth.to_degrees(), range_sep as f64);
    pointing_loc.rad = virtual_height;

    // TODO: Accept old_aacgm option
    // Convert pointing direction position from geocentric lat/lon at virtual height to AACGM
    // magnetic coordinates
    let mut pointing_mag: GeocentricCoords;
    unsafe { pointing_mag = pointing_loc.aacgmv2_convert(); }

    // Make sure pointing_mag_lon lies between +/- 180 degrees
    if pointing_mag.lon - mag_coords.lon > PI {
        pointing_mag.lon -= 2.0 * PI;
    } else if pointing_mag.lon - mag_coords.lon < -PI {
        pointing_mag.lon += 2.0 * PI;
    }

    // Calculate bearing (azimuth) to pointing direction lat/lon from the radar position in magnetic
    // coordinates
    let azimuth = fieldpoint_azimuth(&mag_coords, &pointing_mag);
    Ok((
        MagneticCoords::new(mag_coords.lat, mag_coords.lon),
        azimuth as f32,
        slant_range as f32,
    ))
}


#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use chrono::TimeZone;
    use super::*;

    #[test]
    fn test_fieldpoint_sphere() {
        let rel = 1e-8;
        let start = GeocentricCoords::geo(69.51941199, -133.91889036, 6474.25014983);
        let bearing = -2.56489722;
        let range = 45.0;
        let res = fieldpoint_sphere(start, bearing, range);
        assert_relative_eq!(res.lat.to_degrees(), 69.91724618, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), 226.02920893, max_relative = rel);

        let start = GeocentricCoords::geo(88.01854931, -3.04211092, 7205.35616109);
        let bearing = 134.27412781;
        let res = fieldpoint_sphere(start, bearing, range);
        assert_relative_eq!(res.lat.to_degrees(), 87.75409320, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), 3.51003980, max_relative = rel);

        let start = GeocentricCoords::geo(70.31822085, -70.16621865, 7178.67916230);
        let bearing = 115.48199583;
        let res = fieldpoint_sphere(start, bearing, range);
        assert_relative_eq!(res.lat.to_degrees(), 70.16115496, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), 290.78917074, max_relative = rel);
    }

    #[test]
    fn test_fieldpoint_azimuth() {
        let rel = 1e-7;
        let start = GeocentricCoords::geo(72.13155188, -80.96482800, 0.0);
        let end = GeocentricCoords::geo(72.50725310, -81.35993126, 0.0);
        let res = fieldpoint_azimuth(&start, &end).to_degrees();
        assert_relative_eq!(res, -17.52508198, max_relative = rel);

        let start = GeocentricCoords::geo(80.31958289, 130.46130811, 0.0);
        let end = GeocentricCoords::geo(80.01375090, 129.79161762, 0.0);
        let res = fieldpoint_azimuth(&start, &end).to_degrees();
        assert_relative_eq!(res, -159.16561423, max_relative = rel);

        let start = GeocentricCoords::geo(87.46967400, -174.50635030, 0.0);
        let end = GeocentricCoords::geo(87.37635121, -181.88776787, 0.0);
        let res = fieldpoint_azimuth(&start, &end).to_degrees();
        assert_relative_eq!(res, -101.99796960, max_relative = rel);
    }

    #[test]
    fn test_fieldpoint_height() {
        let rel = 1e-9;
        let point = GeodeticCoords::geo(68.413, -133.769, 0.0);
        let bearing_off_boresight = -24.3;
        let boresight_bearing = 29.5;
        let height = 300.0;
        let slant_range = 180.0;
        let chisham = true;
        let res = fieldpoint_height(&point, bearing_off_boresight, boresight_bearing, height, slant_range, chisham).unwrap();
        assert_relative_eq!(res.lat.to_degrees(), 69.519411986, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), -133.918890359, max_relative = rel);
        assert_relative_eq!(res.rad, 6474.250149832, max_relative = rel);

        let bearing_off_boresight = 21.06;
        let slant_range = 810.0;
        let res = fieldpoint_height(&point, bearing_off_boresight, boresight_bearing, height, slant_range, chisham).unwrap();
        assert_relative_eq!(res.lat.to_degrees(), 71.497889906, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), -117.673408068, max_relative = rel);
        assert_relative_eq!(res.rad, 6717.634104790, max_relative = rel);
    }

    #[test]
    fn test_fieldpoint() {
        let rel = 1e-9;
        let radar_location = GeocentricCoords::geo(68.281017392, -133.769, 6359.668034912);
        let direction = LocalAngularCoords::new(-2.425048105_f64.to_radians(), 38.913774588_f64.to_radians(), 180.0);
        let res = fieldpoint(&radar_location, &direction);
        assert_relative_eq!(res.lat.to_degrees(), 69.519411986, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), -133.918890359, max_relative = rel);
        assert_relative_eq!(res.rad, 6474.250149832, max_relative = rel);
    }

    #[test]
    fn test_rpos_range_beam_az_el() {
        let rel = 1e-6;
        let beam: i32 = 0;
        let range: i32 = 0;
        let year: i32 = 2025;
        let hdw: &HdwInfo = &HdwInfo::new(64, chrono::Utc.with_ymd_and_hms(2025, 7, 12, 0, 0, 0).unwrap()).unwrap();
        let first_range: f32 = 180.0;
        let range_sep: f32 = 45.0;
        let rx_rise: f32 = 0.0;
        let altitude: f32 = 300.0;
        let chisham: bool = true;
        let res = rpos_range_beam_azimuth_elevation(beam, range, year, hdw, first_range, range_sep, rx_rise, altitude, chisham).unwrap();
        assert_relative_eq!(res.az.to_degrees(), -2.564897223, max_relative = rel);
        assert_relative_eq!(res.el.to_degrees(), 7.618535464, max_relative = 1e-2);
        assert_relative_eq!(res.range, 180.0, max_relative = rel);
    }
    #[test]
    fn test_rpos_geo() {
        let rel = 1e-7;
        let center = true;
        let beam_num = 0;
        let range_gate = 0;
        let hdw = HdwInfo::new(64, chrono::Utc.with_ymd_and_hms(2025, 7, 12, 0, 0, 0).unwrap()).unwrap();
        let first_range = 180.0;
        let range_sep = 45.0;
        let rx_rise_time = 0.0;
        let altitude = 300.0;
        let chisham = true;
        let (end, srng) = rpos_geo(
            center,
            beam_num,
            range_gate,
            &hdw,
            first_range,
            range_sep,
            rx_rise_time,
            altitude,
            chisham,
        ).unwrap();
        assert_relative_eq!(end.lat.to_degrees(), 69.51941199, max_relative = rel);
        assert_relative_eq!(end.lon.to_degrees(), -133.91889036, max_relative = rel);
        assert_relative_eq!(end.rad, 6474.25014983, max_relative = rel);
        assert_relative_eq!(srng, 180.0, max_relative = rel);
    }

    #[test]
    fn test_rpos_inv_mag() {
        let rel = 1e-5;

        unsafe {
            aacgmv2_rs::AACGM_v2_SetDateTime(2025, 7, 12, 0, 0, 0);
        }

        let beam_num = 0;
        let range_gate = 0;
        let year = 2025;
        let hdw = HdwInfo::new(64, chrono::Utc.with_ymd_and_hms(2025, 7, 12, 0, 0, 0).unwrap()).unwrap();
        let first_range = 180.0;
        let range_sep = 45.0;
        let rx_rise_time = 0.0;
        let altitude = 300.0;
        let chisham = true;
        let (coords, azimuth, slant_range) = rpos_inv_mag(
            beam_num,
            range_gate,
            year,
            &hdw,
            first_range,
            range_sep,
            rx_rise_time,
            altitude,
            chisham,
            false,
        ).unwrap();
        assert_relative_eq!(coords.lon.to_degrees(), -80.964827995, max_relative = rel);
        assert_relative_eq!(coords.lat.to_degrees(), 72.131551884, max_relative = rel);
        assert_relative_eq!(azimuth.to_degrees() as f64, -17.525081983, max_relative = rel);
        assert_relative_eq!(slant_range, 180.0);
    }

    #[test]
    fn test_igrf_field() {
        let rel = 1e-2;
        let (lat, lon, alt) = (69.51941199, -133.91889036, 6474.25014983);
        let igrf_field = declination(
            lat,
            lon,
            alt,
            Date::from_calendar_date(2025, time::Month::January, 1).unwrap()
        ).unwrap();
        assert_relative_eq!(igrf_field.x, -7334.09740294, max_relative = rel);
        assert_relative_eq!(igrf_field.y, 2496.73900915, max_relative = rel);
        assert_relative_eq!(igrf_field.z, -53940.93134632, max_relative = rel);
    }
}