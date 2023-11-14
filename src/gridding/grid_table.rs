use crate::gridding::grid::GridError;
use crate::utils::scan::{RadarBeam, RadarScan};
use std::f64::consts::PI;
use crate::utils::hdw::HdwInfo;

pub const VELOCITY_ERROR_MIN: f64 = 100.0;  // m/s
pub const POWER_LIN_ERROR_MIN: f64 = 1.0;   // a.u. in linear scale
pub const WIDTH_LIN_ERROR_MIN: f64 = 1.0;   // m/s

pub const RADIUS_EARTH: f64 = 6378.0;   // km TODO: Confirm value

#[derive(Debug)]
pub struct GridBeam {
    pub beam: i32,         // bm in RST
    pub first_range: i32,  // frang in RST, km
    pub range_sep: i32,    // rsep in RST, km
    pub rx_rise: i32,      // rxrise in RST, microseconds?
    pub num_ranges: i32,   // nrang in RST
    pub azimuth: Vec<f64>, // azm in RST, degrees?
    pub ival: Vec<f64>,    // ival in RST
    pub index: Vec<i32>,   // inx in RST
}

#[derive(Debug)]
pub struct GridPoint {
    pub max: i32,                   // max in RST
    pub count: i32,                 // cnt in RST
    pub reference: i32,             // ref in RST
    pub magnetic_lat: f64,          // mlat in RST
    pub magnetic_lon: f64,          // mlon in RST
    pub azimuth: f64,               // azm in RST, degrees?
    pub velocity_median: f64,       // vel.median in RST, m/s
    pub velocity_median_north: f64, // vel.median_n in RST, m/s
    pub velocity_median_east: f64,  // vel.median_e in RST, m/s
    pub velocity_stddev: f64,       // vel.sd in RST, m/s
    pub power_median: f64,          // pwr.median in RST, a.u. in linear scale
    pub power_stddev: f64,          // pwr.sd in RST, a.u. in linear scale
    pub spectral_width_median: f64, // wdt.median in RST, m/s
    pub spectral_width_stddev: f64, // wdt.sd in RST, m/s
}
impl GridPoint {
    pub fn clear(&mut self) {
        self.azimuth = 0.0;
        self.velocity_median_north = 0.0;
        self.velocity_median_east = 0.0;
        self.velocity_stddev = 0.0;
        self.power_median = 0.0;
        self.power_stddev = 0.0;
        self.spectral_width_median = 0.0;
        self.spectral_width_stddev = 0.0;
        self.count = 0;
    }
}

#[derive(Debug)]
pub struct GridTable {
    pub start_time: f64,         // st_time in RST
    pub end_time: f64,           // ed_time in RST
    pub channel: i32,            // chn in RST
    pub status: i32,             // status in RST
    pub station_id: i32,         // st_id in RST
    pub program_id: i32,         // prog_id in RST
    pub num_scans: i32,          // nscan in RST
    pub num_points_npnt: i32,    // npnt in RST, number of grid points
    pub freq: f64,               // freq in RST
    pub noise_mean: f64,         // noise.mean in RST
    pub noise_stddev: f64,       // noise.sd in RST
    pub groundscatter: i32,      // gsct in RST
    pub min_power: f64,          // min[0] in RST, a.u. in linear scale
    pub min_velocity: f64,       // min[1] in RST, m/s
    pub min_spectral_width: f64, // min[2] in RST, m/s
    pub min_velocity_error: f64, // min[3] in RST, m/s
    pub max_power: f64,          // max[0] in RST, a.u. in linear scale
    pub max_velocity: f64,       // max[1] in RST, m/s
    pub max_spectral_width: f64, // max[2] in RST, m/s
    pub max_velocity_error: f64, // max[3] in RST, m/s
    pub num_beams: i32,          // bnum in RST
    pub beams: Vec<GridBeam>,    // bm in RST
    pub num_points_pnum: i32,    // pnum in RST
    pub points: Vec<GridPoint>,  // pnt in RST
}
impl GridTable {
    /// Called GridTableZero in RST
    pub fn clear(&mut self) {
        for p in self.points.iter_mut() {
            p.clear()
        }
    }

    /// Tests whether gridded data should be written to a file.
    /// Called GridTableTest in RST
    pub fn test(mut self, scan: &RadarScan) -> bool {
        let time = (&scan.start_time + &scan.end_time) / 2.0;

        if self.start_time == -1.0 {
            return false;
        }

        if time <= self.end_time {
            return false;
        }

        self.num_points_npnt = 0;

        // Average values across all scans included in the grid table
        let num_scans: &f64 = &(self.num_scans as f64);
        self.freq /= num_scans;
        self.noise_mean /= num_scans;
        self.noise_stddev /= num_scans;

        for point in self.points.iter_mut() {
            if point.count != 0 {
                if point.count <= &self.num_scans * &point.max / 4 {
                    point.count = 0;
                } else {
                    // Update the total number of grid points in the grid table
                    self.num_points_npnt += 1;

                    // Calculate weighted mean of north/east velocity components
                    point.velocity_median_north /= &point.velocity_stddev;
                    point.velocity_median_east /= &point.velocity_stddev;

                    // Calculate the magnitude of weighted mean velocity error
                    point.velocity_median = (&point.velocity_median_north
                        * &point.velocity_median_north
                        + &point.velocity_median_east * &point.velocity_median_east)
                        .sqrt();

                    // Calculate azimuth of weighted mean velocity vector
                    point.azimuth = &point
                        .velocity_median_east
                        .atan2(point.velocity_median_north.clone())
                        * 180.0
                        / PI;

                    // Calculate weighted mean of spectral width and power
                    point.spectral_width_median /= &point.spectral_width_stddev;
                    point.power_median /= &point.power_stddev;

                    // Calculate standard deviation of velocity, power, and spectral width
                    point.velocity_stddev = 1.0 / &point.velocity_stddev.sqrt();
                    point.spectral_width_stddev = 1.0 / &point.spectral_width_stddev.sqrt();
                    point.power_stddev = 1.0 / &point.power_stddev.sqrt();
                }
            }
        }
        self.status = 0;
        true
    }

    /// Returns the index of the pointer to a newly added grid cell in the structure
    /// storing gridded radar data.
    /// Called GridTableAddPoint in RST
    // pub fn add_point

    /// Returns the index of the point in the table whose reference number matches the input.
    /// Called GridTableFindPoint in RST
    pub fn find_point(&self, reference: i32) -> Result<i32, GridError> {
        self.points
            .iter()
            .position(|x| x.reference == reference)
            .map(|x| x as i32)
            .ok_or(GridError::Message(format!(
                "Point {} not in grid table",
                reference
            )))
    }

    /// Adds a grid beam to the grid table.
    /// Called GridTableAddBeam in RST
    pub fn add_beam(&mut self, hdw: &HdwInfo, altitude: f64, time: f64, beam: RadarBeam, chisham: bool, old_aacgm: bool) {
        let velocity_correction: f64 = (2.0 * PI / 86400.0)*RADIUS_EARTH*1000.0*(PI*&hdw.latitude as f64/180.0).cos();
        self.num_beams += 1;

        // TODO: Convert tval to year, month, day, hour, minute, seconds
        for r in 0..beam.num_ranges {
            // TODO: Calculate geographic azimuth and elevation
            // TODO: Calculate magnetic latitude and longitude
            // TODO: Ensure magnetic azimuth and longitude between 0-360 degrees
            // TODO: Calculate magnetic grid cell latitude (eg, 72.1->72.5, 57.8->57.5, etc)
            // TODO: Calculate magnetic grid longitude spacing at grid latitude
            // TODO: Calculate magnetic grid cell longitude
            // TODO: Calculate reference number for cell
            // TODO: Find index of GridPoint corresponding to reference number for cell, make new GridPoint if none found
            // TODO: Update the total number of range gates that map to GridPoint (GridPoint.max)
            // TODO: Set magnetic lat/lon for GridPoint
            // TODO: Set index, magnetic azimuth, inertial velocity correction factor of beam
            // TODO: Return index of beam number added to self
        }
    }

    /// Find the index of the beam in the grid table whose beam number and operating parameters
    /// match those of the input.
    /// Called GridTableFindBeam in RST
    pub fn find_beam(&self, beam: &RadarBeam) -> Result<i32, GridError> {
        self.beams
            .iter()
            .position(|x| {
                x.beam == beam.beam
                    && x.first_range == beam.first_range
                    && x.range_sep == beam.range_sep
                    && x.num_ranges == beam.num_ranges
            })
            .map(|x| x as i32)
            .ok_or(GridError::Message(format!("Beam not found in grid table",)))
    }

    /// Maps radar scan data to an equal-area grid in magnetic coordinates.
    /// Called GridTableMap in RST
    pub fn map (&mut self, scan: &RadarScan, hdw: &HdwInfo, tlen: i32, iflg: i32, altitude: f64, chisham: bool, old_aacgm: bool) {
        let time = (&scan.start_time + &scan.end_time) / 2.0;
        if self.status == 0 {
            self.status = 1;
            self.noise_mean = 0.0;
            self.noise_stddev = 0.0;
            self.freq = 0.0;
            self.num_scans = 0;
            self.start_time = scan.start_time.clone();
            self.end_time = scan.end_time.clone();
            self.station_id = scan.station_id.clone();
        }

        for beam in scan.beams.iter() {
            if beam.beam != -1 {
                match self.find_beam(beam) {
                    Ok(_) => {},
                    Err(_) => self.add_beam(hdw, altitude, time, beam, chisham, old_aacgm),
                };
            }

            for range in 0..beam.num_ranges as usize {
                if beam.scatter[range as usize] != 0 {
                    let velocity_error = beam.cells[range].
                }
            }
            // TODO: A whole lot more...
        }


    }
}
