use crate::error::BackscatterError;
use crate::utils::hdw::HdwInfo;

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
) -> Result<(f32, f32), BackscatterError> {
    let geodetic_lat = hdw.latitude;   // degrees
    let geodetic_lon: hdw.longitude;       // degrees

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise,
    };

    // TODO: Convert all these function calls into library calls. No sense reinventing the wheel
    // Convert center of range/beam cell to geocentric spherical latitude/longitude and distance
    // from the center of the oblate spheroid plus virtual height
    let (spherical_lat, spherical_lon, spherical_alt) = rpos_geo(1, beam, range, hdw, first_range, range_sep, rx_rise_time, altitude, chisham);

    // Convert range/beam position from geocentric spherical coordinates to global Cartesian
    // coordinates
    let (x, y, z) = spherical_to_cartesian(spherical_lat, spherical_lon, spherical_alt);

    // Convert radar site geodetic latitude/longitude to geocentric spherical coordinates and
    // distance from the center to the surface of the oblate spheroid
    let (site_lat, site_lon, site_alt) = geodetic_to_geocentric(1, geodetic_lat, geodetic_lon);

    // Convert radar geocentric coordinates to global Cartesian coordinates
    let (site_x, site_y, site_z) = spherical_to_cartesian(site_lat, site_lon, site_alt);

    // Calculate vector from site to center of range/beam cell
    let dx = x - site_x;
    let dy = y - site_y;
    let dz = z - site_z;

    // Normalize the vector
    let (normed_x, normed_y, normed_z) = norm_vector(dx, dy, dz);

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let (mut local_x, mut local_y, mut local_z) = cartesian_to_local(normed_x, normed_y, normed_z);

    // Normalize the local horizontal vector
    (local_x, local_y, local_z) = norm_vector(local_x, local_y, local_z);

    // Calculate the magnetic field vector at the geocentric spherical range/beam position in
    // local south/east/vertical coordinates
    let (mut bx, mut by, mut bz) = igrf_magnetic_cmp(year, spherical_lat, spherical_lon, spherical_alt)?;

    // Normalize the magnetic field vector
    (bx, by, bz) = norm_vector(bx, by, bz);

    // Calculate a new local vertical component such that the radar-to-range/beam vector becomes
    // orthogonal to the magnetic field at the range/beam position
    local_z = -(bx * local_x + by * local_y) / bz;

    // Normalize the new radar-to-range/beam vector
    (local_x, local_y, local_z) = norm_vector(local_x, local_y, local_z);

    // Calculate the elevation angle of the orthogonal radar-to-range/beam vector
    let elevation = local_z.atan2(local_x * local_x + local_y * local_y);
    let azimuth = local_y.atan2(-local_x);

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
) -> Result<f32, BackscatterError> {
    let geodetic_lat = hdw.latitude;   // degrees
    let geodetic_lon: hdw.longitude;       // degrees

    let rx_rise_time = match rx_rise {
        0.0 => hdw.rx_rise_time,
        _ => rx_rise,
    };

    // TODO: Convert all these function calls into library calls. No sense reinventing the wheel
    // Convert center of range/beam cell to geocentric spherical latitude/longitude and distance
    // from the center of the oblate spheroid plus virtual height
    let (spherical_lat, spherical_lon, spherical_alt) = rpos_geo(1, beam, range, hdw, first_range, range_sep, rx_rise_time, altitude, chisham);

    // Convert range/beam position from geocentric spherical coordinates to global Cartesian
    // coordinates
    let (x, y, z) = spherical_to_cartesian(spherical_lat, spherical_lon, spherical_alt);

    // Convert radar site geodetic latitude/longitude to geocentric spherical coordinates and
    // distance from the center to the surface of the oblate spheroid
    let (site_lat, site_lon, site_alt) = geodetic_to_geocentric(1, geodetic_lat, geodetic_lon);

    // Convert radar geocentric coordinates to global Cartesian coordinates
    let (site_x, site_y, site_z) = spherical_to_cartesian(site_lat, site_lon, site_alt);

    // Calculate vector from site to center of range/beam cell
    let dx = x - site_x;
    let dy = y - site_y;
    let dz = z - site_z;

    // Normalize the vector
    let (normed_x, normed_y, normed_z) = norm_vector(dx, dy, dz);

    // Convert the normalized vector from radar-to-range/beam cell into local south/east/vertical
    // (horizontal) coordinates
    let (mut local_x, mut local_y, mut local_z) = cartesian_to_local(normed_x, normed_y, normed_z);

    // Normalize the local horizontal vector
    (local_x, local_y, local_z) = norm_vector(local_x, local_y, local_z);

    // Calculate the magnetic field vector at the geocentric spherical range/beam position in
    // local south/east/vertical coordinates
    let (mut bx, mut by, mut bz) = igrf_magnetic_cmp(year, spherical_lat, spherical_lon, spherical_alt)?;

    // Normalize the magnetic field vector
    (bx, by, bz) = norm_vector(bx, by, bz);

    // Calculate a new local vertical component such that the radar-to-range/beam vector becomes
    // orthogonal to the magnetic field at the range/beam position
    local_z = -(bx * local_x + by * local_y) / bz;

    // Normalize the new radar-to-range/beam vector
    (local_x, local_y, local_z) = norm_vector(local_x, local_y, local_z);

    // Calculate the elevation angle of the orthogonal radar-to-range/beam vector
    let azimuth = local_y.atan2(-local_x);

    // Calculate virtual height of range/beam position
    let virtual_height = spherical_alt - site_alt;

    // TODO: Accept old_aacgm option
    // Convert range/beam position from geocentric lat/lon at virtual height to AACGM magnetic
    // lat/lon
    let (mag_lat, mag_lon) = aacgm_v2_convert(spherical_lat, spherical_lon, virtual_height, 0)?;

    // Calculate pointing direction lat/lon given distance and bearing from the radar position
    // at the field point radius
    let (pointing_lat, pointing_lon) = fieldpoint_sphere(spherical_lat, spherical_lon, spherical_alt, azimuth, range_sep);

    // TODO: Accept old_aacgm option
    // Convert pointing direction position from geocentric lat/lon at virtual height to AACGM
    // magnetic coordinates
    let (pointing_mag_lat, pointing_mag_lon) = aacgm_v2_convert(pointing_lat, pointing_lon, virtual_height, 0)?;

    // Make sure pointing_mag_lon lies between +/- 180 degrees
    if pointing_mag_lon - mag_lon > 180.0 {
        pointing_mag_lon -= 360.0;
    }
    else if pointing_mag_lon - mag_lon < -180.0 {
        pointing_mag_lon += 360.0;
    }

    // Calculate bearing (azimuth) to pointing direction lat/lon from the radar position in magnetic
    // coordinates
    let azimuth = fieldpoint_azimuth(mag_lat, mag_lon, pointing_mag_lat, pointing_mag_lon);

    Ok(azimuth)
}