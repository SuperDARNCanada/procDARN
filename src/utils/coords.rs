use std::f64::consts::PI;
use std::fmt::Display;

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct GeocentricCoords {
    /// Latitude (radians)
    pub lat: f64,
    /// Longitude (radians)
    pub lon: f64,
    /// Distance from the center of the Earth (km)
    pub rad: f64,
}
impl GeocentricCoords {
    pub fn new(lat: f64, lon: f64, rad: f64) -> GeocentricCoords {
        GeocentricCoords { lat, lon, rad }
    }
    /// Constructor with input `lat`, `lon` in degrees.
    pub(crate) fn geo(lat: f64, lon: f64, rad: f64) -> GeocentricCoords {
        GeocentricCoords::new(lat.to_radians(), lon.to_radians(), rad)
    }

    /// Converts to geodetic coordinates. The WGS84 Earth model is used.
    pub fn to_geodetic(&self) -> GeodeticCoords {
        let semi_major_axis: f64 = 6378.137;
        let flattening: f64 = 1.0 / 298.257223563;
        let semi_minor_axis: f64 = semi_major_axis * (1.0 - flattening);
        let second_eccentricity_squared: f64 =
            (semi_major_axis * semi_major_axis) / (semi_minor_axis * semi_minor_axis) - 1.0;

        let gdlat = ((semi_major_axis * semi_major_axis) / (semi_minor_axis * semi_minor_axis)
            * self.lat.tan())
        .atan();
        let gdlon = self.lon;

        let rho = semi_major_axis
            / (1.0 + second_eccentricity_squared * self.lat.sin() * self.lat.sin()).sqrt();

        GeodeticCoords::new(gdlat, gdlon, rho)
    }

    pub fn to_cartesian(&self) -> CartesianCoords {
        let x = self.rad * self.lat.cos() * self.lon.cos();
        let y = self.rad * self.lat.cos() * self.lon.sin();
        let z = self.rad * self.lat.sin();
        CartesianCoords { x, y, z }
    }

    /// Convert a Cartesian vector `v` centered at `self` into local south/east/vertical coordinates.
    pub(crate) fn cartesian_to_local(&self, v: &CartesianCoords) -> LocalCartesianCoords {
        // Rotate v about the z-axis by the longitude
        let sx = self.lon.cos() * v.x + self.lon.sin() * v.y;
        let sy = -self.lon.sin() * v.x + self.lon.cos() * v.y;
        let sz = v.z;

        // Calculate the colatitude
        let colat = PI / 2.0 - self.lat;

        // Rotate the vector about the east-axis by the colatitude
        let tx = colat.cos() * sx - colat.sin() * sz;
        let ty = sy;
        let tz = colat.sin() * sx + colat.cos() * sz;

        LocalCartesianCoords::new(tx, ty, tz)
    }

    pub(crate) unsafe fn aacgmv2_convert(&self) -> GeocentricCoords {
        let mut mag_coords = GeocentricCoords::default();
        unsafe {
            aacgmv2_rs::AACGM_v2_Convert(
                self.lat.to_degrees(),
                self.lon.to_degrees(),
                self.rad,
                &mut mag_coords.lat,
                &mut mag_coords.lon,
                &mut mag_coords.rad,
                0,
            );
        }
        mag_coords.lat = mag_coords.lat.to_radians();
        mag_coords.lon = mag_coords.lon.to_radians();

        mag_coords
    }
}
impl Display for GeocentricCoords {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "({}, {}, {})",
            self.lat.to_degrees(),
            self.lon.to_degrees(),
            self.rad
        ))
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct GeodeticCoords {
    /// Latitude (degrees)
    pub lat: f64,
    /// Longitude (degrees)
    pub lon: f64,
    /// Distance from the center of the Earth (km)
    pub rad: f64,
}

impl GeodeticCoords {
    pub fn new(lat: f64, lon: f64, rad: f64) -> GeodeticCoords {
        GeodeticCoords { lat, lon, rad }
    }
    /// Constructor with input `lat`, `lon` in degrees.
    pub fn geo(lat: f64, lon: f64, rad: f64) -> GeodeticCoords {
        GeodeticCoords::new(lat.to_radians(), lon.to_radians(), rad)
    }

    /// Converts to geocentric spherical coordinates. The WGS84 Earth model is used.
    pub fn to_geocentric(&self) -> GeocentricCoords {
        let semi_major_axis: f64 = 6378.137;
        let flattening: f64 = 1.0 / 298.257223563;
        let semi_minor_axis: f64 = semi_major_axis * (1.0 - flattening);
        let second_eccentricity_squared: f64 =
            (semi_major_axis * semi_major_axis) / (semi_minor_axis * semi_minor_axis) - 1.0;

        let gclat = ((semi_minor_axis * semi_minor_axis) / (semi_major_axis * semi_major_axis)
            * self.lat.tan())
        .atan();
        let mut gclon = self.lon;
        if gclon.to_degrees() > 180.0 {
            gclon -= 360.0_f64.to_radians();
        }
        let rho = semi_major_axis
            / (1.0 + second_eccentricity_squared * gclat.sin() * gclat.sin()).sqrt();

        GeocentricCoords::new(gclat, gclon, rho)
    }

    /// Corrects a look direction at geodetic coordinates to account for the oblateness of the Earth.
    ///
    /// `az` and `el` are given in radians.
    ///
    /// Called geocnvrt in cnvtcoord.c of RST
    pub(crate) fn correct_look_dir(&self, az: f64, el: f64) -> (f64, f64) {
        let kxg = el.cos() * az.sin();
        let kyg = el.cos() * az.cos();
        let kzg = el.sin();

        let point_gc = self.to_geocentric();
        let del = self.lat - point_gc.lat;

        let kxr = kxg;
        let kyr = kyg * del.cos() + kzg * del.sin();
        let kzr = -kyg * del.sin() + kzg * del.cos();

        let ral = kxr.atan2(kyr);
        let rel = (kzr / (kxr * kxr + kyr * kyr).sqrt()).atan();

        (ral, rel)
    }
}
impl Display for GeodeticCoords {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "({}, {}, {})",
            self.lat.to_degrees(),
            self.lon.to_degrees(),
            self.rad
        ))
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct CartesianCoords {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl CartesianCoords {
    pub fn new(x: f64, y: f64, z: f64) -> CartesianCoords {
        CartesianCoords { x, y, z }
    }
    /// Scale to unit length.
    pub fn norm(&mut self) {
        let len = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        self.x /= len;
        self.y /= len;
        self.z /= len;
    }
}
impl std::ops::Sub for CartesianCoords {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}
impl Display for CartesianCoords {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_fmt(format_args!("({}, {}, {})", self.x, self.y, self.z))
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct MagneticCoords {
    pub lat: f64,
    pub lon: f64,
}
impl MagneticCoords {
    pub fn new(lat: f64, lon: f64) -> MagneticCoords {
        MagneticCoords { lat, lon }
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct LocalCartesianCoords {
    /// Distance in km
    pub south: f64,
    /// Distance in km
    pub east: f64,
    /// Distance in km
    pub up: f64,
}
impl LocalCartesianCoords {
    pub fn new(south: f64, east: f64, up: f64) -> LocalCartesianCoords {
        LocalCartesianCoords { south, east, up }
    }
    /// Scale to unit length.
    pub fn norm(&mut self) {
        let len = (self.south * self.south + self.east * self.east + self.up * self.up).sqrt();
        self.south /= len;
        self.east /= len;
        self.up /= len;
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct LocalAngularCoords {
    /// East of North (radians)
    pub az: f64,
    /// Up from horizon (radians)
    pub el: f64,
    /// Distance in km
    pub range: f64,
}
impl LocalAngularCoords {
    pub fn new(az: f64, el: f64, range: f64) -> LocalAngularCoords {
        LocalAngularCoords { az, el, range }
    }
    /// Constructor with inputs in degrees.
    pub fn from_degrees(az: f64, el: f64, range: f64) -> LocalAngularCoords {
        LocalAngularCoords::new(az.to_radians(), el.to_radians(), range)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_geocentric_to_cartesian() {
        let rel = 1e-7;
        let p = GeocentricCoords::geo(69.519412, -133.91889, 6474.25015);
        let res = p.to_cartesian();
        assert_relative_eq!(res.x, -1571.284220, max_relative = rel);
        assert_relative_eq!(res.y, -1631.728799, max_relative = rel);
        assert_relative_eq!(res.z, 6065.017892, max_relative = rel);
    }

    #[test]
    fn test_geodetic_to_geocentric() {
        let rel = 1e-7;
        let p = GeodeticCoords::geo(68.413, -133.769, 0.0);
        let res = p.to_geocentric();
        assert_relative_eq!(res.lat.to_degrees(), 68.281017, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), -133.769, max_relative = rel);
        assert_relative_eq!(res.rad, 6359.668035, max_relative = rel);
    }

    #[test]
    fn test_geocentric_to_geodetic() {
        let rel = 1e-9;
        let p = GeocentricCoords::geo(69.519411986, -133.918890359, 0.0);
        let res = p.to_geodetic();
        assert_relative_eq!(res.lat.to_degrees(), 69.645235706, max_relative = rel);
        assert_relative_eq!(res.lon.to_degrees(), -133.918890359, max_relative = rel);
        assert_relative_eq!(res.rad, 6359.358742609, max_relative = rel);
    }

    #[test]
    fn test_cartesian_to_local() {
        let rel = 1e-9;
        let loc = GeocentricCoords::geo(69.519411986, -133.918890359, 6474.25015);
        let v = CartesianCoords::new(0.315016678, 0.376445908, 0.871236461);
        let res = loc.cartesian_to_local(&v);
        assert_relative_eq!(res.south, -0.763555664, max_relative = rel);
        assert_relative_eq!(res.east, -0.034204109, max_relative = rel);
        assert_relative_eq!(res.up, 0.644835504, max_relative = rel);
    }

    #[test]
    fn test_correct_look_dir() {
        let rel = 1e-9;
        let point = GeodeticCoords::geo(68.413, -133.769, 0.0);
        let az = -2.429550020_f64;
        let el = 38.913774588_f64;
        let (az, el) = point.correct_look_dir(az.to_radians(), el.to_radians());
        assert_relative_eq!(az.to_degrees(), -2.425048105, max_relative = rel);
        assert_relative_eq!(el.to_degrees(), 38.781910399, max_relative = rel);
    }

    #[test]
    fn test_aacgmv2_convert() {
        let rel = 1e-7;

        let point = GeocentricCoords::geo(69.917246, 226.029209, 114.891407);
        let mag_point: GeocentricCoords;
        unsafe {
            mag_point = point.aacgmv2_convert();
        }
        assert_relative_eq!(mag_point.lat.to_degrees(), 72.507253, max_relative = rel);
        assert_relative_eq!(mag_point.lon.to_degrees(), -81.359931, max_relative = rel);
        assert_relative_eq!(mag_point.rad, 1.016164, max_relative = rel);
    }
}
