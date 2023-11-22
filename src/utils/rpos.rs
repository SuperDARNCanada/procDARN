use crate::error::BackscatterError;
use crate::gridding::grid_table::RADIUS_EARTH;
use crate::utils::hdw::HdwInfo;
use geodesy::prelude::*;
use igrf::declination;
use time::Date;
use std::f64::consts::PI;

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

/// Calculates the slant range to a range gate
fn slant_range(
    first_range: f64,
    range_sep: f64,
    rx_rise: f64,
    range_edge: f64,
    range_gate: i32,
) -> f64 {
    let lag_to_first_range = first_range * 20.0 / 3.0;
    let sample_separation = range_sep * 20.0 / 3.0;
    (lag_to_first_range - rx_rise + (range_gate - 1) * sample_separation + range_edge) * 0.15
}

/// This function converts a gate/beam coordinate to geographic position. The height of the
/// transformation is given by height - if this value is less than 90 then it is assumed to be the
/// elevation angle from the radar. If center is not equal to zero, then the calculation is assumed
/// to be for the center of the cell, not the edge. The calculated values are returned in geocentric
/// coordinates.
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
        range_edge = -0.5 * range_sep * 20 / 3;
    }

    let rx_rise = match rx_rise_time {
        0.0 => hdw.rx_rise_time as f64,
        _ => rx_rise_time,
    };

    let offset = hdw.max_num_beams / 2.0 - 0.5;

    // Calculate deviation from boresight in degrees
    let psi = hdw.beam_separation * (beam_num - offset) + beam_edge + hdw.boresight_shift;

    // Calculate the slant range to the range gate in km
    let distance = slant_range(first_range, range_sep, rx_rise, range_edge, range_gate + 1);

    // If the input altitude is below 90, then it is actually an input elevation angle in degrees.
    // If so, we calculate the field point height
    let field_point_height: f64;
    if altitude < 90.0 {
        field_point_height = -RADIUS_EARTH
            + ((RADIUS_EARTH * RADIUS_EARTH)
                + 2 * distance * RADIUS_EARTH * (altitude * PI / 180.0).sin()
                + distance * distance)
                .sqrt();
    } else {
        field_point_height = altitude;
    }

    // Calculate the geocentric coordinates of the field point
    field_point_calculation(hdw, psi, field_point_height, distance, chisham)
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
    let igrf_field = declination(cell_geoc[1], cell_geoc[0], cell_geoc[2] as u32, Date::from_calendar_date(year, time::Month::January, 1)?)?;

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
    let igrf_field = declination(cell_geoc[1], cell_geoc[0], cell_geoc[2] as u32, Date::from_calendar_date(year, time::Month::January, 1)?)?;

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
    let (pointing_lat, pointing_lon) = fieldpoint_sphere(
        cell_geoc,
        azimuth,
        range_sep,
    );

    // TODO: Accept old_aacgm option
    // Convert pointing direction position from geocentric lat/lon at virtual height to AACGM
    // magnetic coordinates
    let (pointing_mag_lat, pointing_mag_lon) =
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
