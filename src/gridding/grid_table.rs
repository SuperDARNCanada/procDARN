use crate::gridding::grid::GridError;
use crate::utils::hdw::HdwInfo;
use crate::utils::rpos::{rpos_inv_mag, rpos_range_beam_azimuth_elevation};
use crate::utils::scan::{RadarBeam, RadarScan};
use chrono::{DateTime, Datelike, TimeDelta, Timelike, Utc};
use dmap::formats::grid::GridRecord;
use dmap::record::Record;
use dmap::types::DmapField;
use indexmap::IndexMap;
use numpy::ndarray::array;
use numpy::ndarray::Array;
use std::f32::consts::PI;
use std::iter;

pub const GRID_REVISION_MAJOR: i32 = 2;
pub const GRID_REVISION_MINOR: i32 = 0;
pub const VELOCITY_ERROR_MIN: f32 = 100.0; // m/s
pub const POWER_LIN_ERROR_MIN: f32 = 1.0; // a.u. in linear scale
pub const WIDTH_LIN_ERROR_MIN: f32 = 1.0; // m/s

pub const RADIUS_EARTH: f32 = 6371.2; // km

#[derive(Debug, Default)]
pub struct GridBeam {
    pub beam: i32,             // bm in RST
    pub first_range: i32,      // frang in RST, km
    pub range_sep: i32,        // rsep in RST, km
    pub rx_rise: i32,          // rxrise in RST, microseconds?
    pub num_ranges: i32,       // nrang in RST
    pub azimuth: Vec<f32>,     // azm in RST, radians
    pub slant_range: Vec<f32>, // srng in RST, km
    pub ival: Vec<f32>,        // ival in RST
    pub index: Vec<i32>,       // inx in RST
}

#[derive(Debug, Default)]
pub struct GridPoint {
    pub max: i32,                   // max in RST
    pub count: i32,                 // cnt in RST
    pub reference: i32,             // ref in RST
    pub magnetic_lat: f32,          // mlat in RST
    pub magnetic_lon: f32,          // mlon in RST
    pub azimuth: f32,               // azm in RST, degrees
    pub slant_range: f32,           // srng in RST, km
    pub velocity_median: f32,       // vel.median in RST, m/s
    pub velocity_median_north: f32, // vel.median_n in RST, m/s
    pub velocity_median_east: f32,  // vel.median_e in RST, m/s
    pub velocity_stddev: f32,       // vel.sd in RST, m/s
    pub power_median: f32,          // pwr.median in RST, a.u. in linear scale
    pub power_stddev: f32,          // pwr.sd in RST, a.u. in linear scale
    pub spectral_width_median: f32, // wdt.median in RST, m/s
    pub spectral_width_stddev: f32, // wdt.sd in RST, m/s
}
impl GridPoint {
    pub fn clear(&mut self) {
        self.azimuth = 0.0;
        self.slant_range = 0.0;
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

#[derive(Debug, Default)]
pub struct GridTable {
    pub start_time: DateTime<Utc>, // st_time in RST
    pub end_time: DateTime<Utc>,   // ed_time in RST
    pub channel: i16,              // chn in RST
    pub status: i32,               // status in RST
    pub station_id: i16,           // st_id in RST
    pub program_id: i32,           // prog_id in RST
    pub num_scans: i32,            // nscan in RST
    pub num_points_npnt: i32,      // npnt in RST, number of grid points
    pub freq: f32,                 // freq in RST
    pub noise_mean: f32,           // noise.mean in RST
    pub noise_stddev: f32,         // noise.sd in RST
    pub groundscatter: i32,        // gsct in RST
    pub min_power: f32,            // min[0] in RST, a.u. in linear scale
    pub min_velocity: f32,         // min[1] in RST, m/s
    pub min_spectral_width: f32,   // min[2] in RST, m/s
    pub min_velocity_error: f32,   // min[3] in RST, m/s
    pub max_power: f32,            // max[0] in RST, a.u. in linear scale
    pub max_velocity: f32,         // max[1] in RST, m/s
    pub max_spectral_width: f32,   // max[2] in RST, m/s
    pub max_velocity_error: f32,   // max[3] in RST, m/s
    pub num_beams: i32,            // bnum in RST
    pub beams: Vec<GridBeam>,      // bm in RST
    pub num_points_pnum: i32,      // pnum in RST
    pub points: Vec<GridPoint>,    // pnt in RST
}
impl GridTable {
    /// Called GridTableZero in RST
    pub fn clear(&mut self) {
        for p in self.points.iter_mut() {
            p.clear()
        }
    }

    /// Tests whether gridded data should be written to a file.
    /// Modifies self, averaging measurements within grid cells if returning true.
    ///
    /// Called `GridTableTest` in RST.
    pub fn test(&mut self, scan: &RadarScan) -> bool {
        let time_micros =
            (scan.start_time.timestamp_micros() + scan.end_time.timestamp_micros()) / 2;
        let time: DateTime<Utc>;
        match DateTime::from_timestamp_micros(time_micros) {
            Some(x) => {
                time = x;
            }
            None => return false,
        }

        if self.start_time == DateTime::<Utc>::default() {
            return false;
        }

        if time <= self.end_time {
            return false;
        }

        self.num_points_npnt = 0;

        // Average values across all scans included in the grid table
        let num_scans: &f32 = &(self.num_scans as f32);
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
                    point.azimuth = point
                        .velocity_median_east
                        .atan2(point.velocity_median_north.clone())
                        .to_degrees();

                    // Calculate average slant range of velocity vector
                    point.slant_range = point.slant_range / point.count as f32;

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
    pub fn add_point(&mut self) -> usize {
        self.points.push(GridPoint::default());
        self.num_points_pnum += 1;
        self.points.len() - 1
    }

    /// Returns the index of the point in the table whose reference number matches the input.
    /// Called GridTableFindPoint in RST
    pub fn find_point(&self, reference: i32) -> Option<usize> {
        self.points.iter().position(|x| x.reference == reference)
    }

    /// Adds a grid beam to the grid table.
    /// Called GridTableAddBeam in RST
    pub fn add_beam(
        &mut self,
        hdw: &HdwInfo,
        altitude: f32,
        time: DateTime<Utc>,
        scan_beam: &RadarBeam,
        chisham: bool,
        old_aacgm: bool,
    ) -> Result<usize, GridError> {
        let velocity_correction: f32 =
            (2.0 * PI / 86400.0) * RADIUS_EARTH * 1000.0 * hdw.latitude.to_radians().cos();
        self.num_beams += 1;

        let mut grid_beam = GridBeam {
            beam: scan_beam.beam,
            first_range: scan_beam.first_range,
            range_sep: scan_beam.range_sep,
            rx_rise: scan_beam.rx_rise,
            num_ranges: scan_beam.num_ranges,
            ..Default::default()
        };

        for range in 0..grid_beam.num_ranges {
            // Calculate geographic azimuth and elevation to scatter point
            let result = rpos_range_beam_azimuth_elevation(
                grid_beam.beam,
                range,
                time.year(),
                hdw,
                grid_beam.first_range as f32,
                grid_beam.range_sep as f32,
                grid_beam.rx_rise as f32,
                altitude,
                chisham,
            )?;
            let azimuth_geo = result.az;

            // Calculate magnetic latitude, longitude, azimuth, and slant range of scatter point
            let (mut mag_loc, mut azimuth_mag, srng_mag) = rpos_inv_mag(
                grid_beam.beam,
                range,
                time.year(),
                hdw,
                grid_beam.first_range as f32,
                grid_beam.range_sep as f32,
                grid_beam.rx_rise as f32,
                altitude,
                chisham,
                old_aacgm,
            )?;

            // Ensure magnetic azimuth and longitude between 0-360 degrees
            if azimuth_mag < 0.0 {
                azimuth_mag += 2.0 * PI;
            }
            if mag_loc.lon < 0.0 {
                mag_loc.lon += 2.0 * PI as f64;
            }

            // Calculate magnetic grid cell latitude, (e.g. 72.1->72.5, 57.8->57.5, etc)
            let grid_lat: f32;
            if mag_loc.lat > 0.0 {
                grid_lat = mag_loc.lat.to_degrees().floor() as f32 + 0.5;
            } else {
                grid_lat = mag_loc.lat.to_degrees().floor() as f32 - 0.5;
            }

            // Calculate magnetic grid longitude spacing at grid latitude
            let lon_spacing = (360.0 * grid_lat.abs().to_radians().cos() + 0.5).floor() / 360.0;

            // Calculate magnetic grid cell longitude
            let grid_lon =
                ((mag_loc.lon.to_degrees() as f32 * lon_spacing).floor() + 0.5) / lon_spacing;

            // Calculate reference number for cell
            let reference: i32;
            if mag_loc.lat > 0.0 {
                reference = (1000.0 * mag_loc.lat.to_degrees().floor() as f32
                    + (mag_loc.lon.to_degrees() as f32 * lon_spacing).floor())
                    as i32;
            } else {
                reference = (-1000.0 * (-1.0 * mag_loc.lat.to_degrees()).floor() as f32
                    - (mag_loc.lon.to_degrees() as f32 * lon_spacing).floor())
                    as i32;
            }

            // Find GridPoint corresponding to reference number for cell, make new GridPoint if none found
            let index = match self.find_point(reference) {
                Some(x) => x,
                None => self.add_point(),
            };
            let point = &mut self.points[index];

            // Update the total number of range gates that map to GridPoint (GridPoint.max)
            point.reference = reference;
            point.max += 1;

            // Set magnetic lat/lon for GridPoint
            point.magnetic_lat = grid_lat;
            point.magnetic_lon = grid_lon;

            // Set index, magnetic azimuth, slant range, and inertial velocity correction factor of beam
            grid_beam.index.push(index as i32);
            grid_beam.azimuth.push(azimuth_mag);
            grid_beam.slant_range.push(srng_mag);
            grid_beam
                .ival
                .push(velocity_correction * (-azimuth_geo.sin() as f32));
        }
        self.beams.push(grid_beam);
        // Return index of beam number added to self
        Ok((self.num_beams - 1) as usize)
    }

    /// Find the index of the beam in the grid table whose beam number and operating parameters
    /// match those of the input.
    /// Called GridTableFindBeam in RST
    pub fn find_beam(&self, beam: &RadarBeam) -> Option<usize> {
        self.beams.iter().position(|x| {
            x.beam == beam.beam
                && x.first_range == beam.first_range
                && x.range_sep == beam.range_sep
                && x.num_ranges == beam.num_ranges
        })
    }

    /// Maps radar scan data to an equal-area grid in magnetic coordinates.
    /// Called GridTableMap in RST
    pub fn map(
        &mut self,
        scan: &RadarScan,
        hdw: &HdwInfo,
        tlen: i32,
        iflg: bool,
        altitude: f32,
        chisham: bool,
        old_aacgm: bool,
    ) -> Result<(), GridError> {
        let time_micros =
            (scan.start_time.timestamp_micros() + scan.end_time.timestamp_micros()) / 2;
        let time = DateTime::from_timestamp_micros(time_micros).ok_or_else(|| {
            GridError::InvalidFitacf("Invalid datetime for GridTable".to_string())
        })?;
        if self.status == 0 {
            // set to zero by self.test()
            self.status = 1;
            self.noise_mean = 0.0;
            self.noise_stddev = 0.0;
            self.freq = 0.0;
            self.num_scans = 0;
            self.clear();
            self.start_time = scan.start_time.clone();
            self.end_time = scan.start_time.clone() + TimeDelta::seconds(tlen as i64);
            self.station_id = scan.station_id.clone();
        }

        for scan_beam in scan.beams.iter() {
            if scan_beam.beam == -1 {
                continue;
            }
            let beam_index = match self.find_beam(scan_beam) {
                Some(i) => i,
                None => self.add_beam(hdw, altitude, time, scan_beam, chisham, old_aacgm)?,
            };
            let grid_beam = &self.beams[beam_index];

            for range in 0..scan_beam.num_ranges.clone() as usize {
                if scan_beam.scatter[range] == 0 {
                    continue;
                }

                let velocity_error = scan_beam.cells[range]
                    .velocity_error
                    .max(VELOCITY_ERROR_MIN);
                let power_lin_error = scan_beam.cells[range]
                    .power_lin_error
                    .max(POWER_LIN_ERROR_MIN);
                let width_lin_error = scan_beam.cells[range]
                    .spectral_width_lin_error
                    .max(WIDTH_LIN_ERROR_MIN);

                // Get grid cell of radar beam/gate measurement
                let grid_cell = &mut self.points[grid_beam.index[range] as usize];

                // Add slant range of gate measurement
                grid_cell.slant_range += grid_beam.slant_range[range];

                if iflg {
                    grid_cell.velocity_median_north -= (scan_beam.cells[range].velocity
                        + grid_beam.ival[range])
                        * grid_beam.azimuth[range].cos()
                        / (velocity_error * velocity_error);
                    grid_cell.velocity_median_east -= (scan_beam.cells[range].velocity
                        + grid_beam.ival[range])
                        * grid_beam.azimuth[range].sin()
                        / (velocity_error * velocity_error);
                } else {
                    grid_cell.velocity_median_north -= scan_beam.cells[range].velocity
                        * grid_beam.azimuth[range].cos()
                        / (velocity_error * velocity_error);
                    grid_cell.velocity_median_east -= scan_beam.cells[range].velocity
                        * grid_beam.azimuth[range].sin()
                        / (velocity_error * velocity_error);
                }

                grid_cell.power_median +=
                    scan_beam.cells[range].power_lin / (power_lin_error * power_lin_error);
                grid_cell.spectral_width_median +=
                    scan_beam.cells[range].spectral_width_lin / (width_lin_error * width_lin_error);

                grid_cell.velocity_stddev += 1.0 / (velocity_error * velocity_error);
                grid_cell.power_stddev += 1.0 / (power_lin_error * power_lin_error);
                grid_cell.spectral_width_stddev += 1.0 / (width_lin_error * width_lin_error);
                grid_cell.count += 1;
            }
        }

        // TODO: Check if somehow all beams in scan not considered?

        let mut freq = 0.0;
        let mut noise: f64 = 0.0;
        let mut variance: f64 = 0.0;
        let mut count: f64 = 0.0;

        for scan_beam in scan.beams.iter().filter(|beam| beam.beam != -1) {
            self.program_id = scan_beam.program_id;

            // Sum the frequency and noise values
            freq += scan_beam.freq as f32;
            noise += scan_beam.noise as f64;
            count += 1.0;
        }

        // Average frequency and noise over all beams in scan
        freq = freq / count as f32;
        noise = noise / count;

        for scan_beam in scan.beams.iter().filter(|beam| beam.beam != -1) {
            variance += (scan_beam.noise as f64 - noise) * (scan_beam.noise as f64 - noise);
        }
        self.noise_mean += noise as f32;
        self.noise_stddev += (variance / count).sqrt() as f32;
        self.freq += freq;
        self.num_scans += 1;

        Ok(())
    }

    /// Converts the GridTable to a GridRecord for writing to file.
    /// Equivalent to GridTableWrite in RST.
    pub fn to_dmap_record(&self, extended_flag: bool) -> Result<GridRecord, GridError> {
        let mut grid_rec: IndexMap<String, DmapField> = IndexMap::new();

        // Find the valid points in the grid
        let valid_points: Vec<&GridPoint> = self.points.iter().filter(|&p| p.count > 0).collect();
        let num_points = valid_points.len();

        // These vector fields require accessing the points of grid_table
        let magnetic_lat: Vec<f32> = valid_points.iter().map(|&p| p.magnetic_lat).collect();
        let magnetic_lon = valid_points.iter().map(|&p| p.magnetic_lon).collect();
        let azimuth = valid_points.iter().map(|&p| p.azimuth).collect();
        let slant_range = valid_points.iter().map(|&p| p.slant_range).collect();
        let index: Vec<i32> = valid_points.iter().map(|&p| p.reference).collect();
        let velocity_median = valid_points.iter().map(|&p| p.velocity_median).collect();
        let velocity_stddev = valid_points.iter().map(|&p| p.velocity_stddev).collect();
        let power_median = valid_points.iter().map(|&p| p.power_median).collect();
        let power_stddev = valid_points.iter().map(|&p| p.power_stddev).collect();
        let spectral_width_median = valid_points
            .iter()
            .map(|&p| p.spectral_width_median)
            .collect();
        let spectral_width_stddev: Vec<f32> = valid_points
            .iter()
            .map(|&p| p.spectral_width_stddev)
            .collect();
        let station_ids: Vec<i16> = iter::repeat(self.station_id)
            .take(valid_points.len())
            .collect();
        let channels: Vec<i16> = iter::repeat(self.channel)
            .take(valid_points.len())
            .collect();

        grid_rec.insert(
            "start.year".to_string(),
            (self.start_time.year() as i16).into(),
        );
        grid_rec.insert(
            "start.month".to_string(),
            (self.start_time.month() as i16).into(),
        );
        grid_rec.insert(
            "start.day".to_string(),
            (self.start_time.day() as i16).into(),
        );
        grid_rec.insert(
            "start.hour".to_string(),
            (self.start_time.hour() as i16).into(),
        );
        grid_rec.insert(
            "start.minute".to_string(),
            (self.start_time.minute() as i16).into(),
        );
        grid_rec.insert(
            "start.second".to_string(),
            (self.start_time.second() as f64 + (self.start_time.nanosecond() as f64) * 1e-9).into(),
        );
        grid_rec.insert("end.year".to_string(), (self.end_time.year() as i16).into());
        grid_rec.insert(
            "end.month".to_string(),
            (self.end_time.month() as i16).into(),
        );
        grid_rec.insert("end.day".to_string(), (self.end_time.day() as i16).into());
        grid_rec.insert("end.hour".to_string(), (self.end_time.hour() as i16).into());
        grid_rec.insert(
            "end.minute".to_string(),
            (self.end_time.minute() as i16).into(),
        );
        grid_rec.insert(
            "end.second".to_string(),
            (self.end_time.second() as f64 + (self.end_time.nanosecond() as f64) * 1e-9).into(),
        );
        grid_rec.insert(
            "stid".to_string(),
            array![self.station_id].into_dyn().into(),
        );
        grid_rec.insert(
            "channel".to_string(),
            array![self.channel].into_dyn().into(),
        );
        grid_rec.insert(
            "nvec".to_string(),
            array![num_points as i16].into_dyn().into(),
        );
        grid_rec.insert("freq".to_string(), array![self.freq].into_dyn().into());
        grid_rec.insert(
            "major.revision".to_string(),
            array![GRID_REVISION_MAJOR as i16].into_dyn().into(),
        );
        grid_rec.insert(
            "minor.revision".to_string(),
            array![GRID_REVISION_MINOR as i16].into_dyn().into(),
        );
        grid_rec.insert(
            "program.id".to_string(),
            array![self.program_id as i16].into_dyn().into(),
        );
        grid_rec.insert(
            "noise.mean".to_string(),
            array![self.noise_mean].into_dyn().into(),
        );
        grid_rec.insert(
            "noise.sd".to_string(),
            array![self.noise_stddev].into_dyn().into(),
        );
        grid_rec.insert(
            "gsct".to_string(),
            array![self.groundscatter as i16].into_dyn().into(),
        );
        grid_rec.insert(
            "v.min".to_string(),
            array![self.min_velocity].into_dyn().into(),
        );
        grid_rec.insert(
            "v.max".to_string(),
            array![self.max_velocity].into_dyn().into(),
        );
        grid_rec.insert(
            "p.min".to_string(),
            array![self.min_power].into_dyn().into(),
        );
        grid_rec.insert(
            "p.max".to_string(),
            array![self.max_power].into_dyn().into(),
        );
        grid_rec.insert(
            "w.min".to_string(),
            array![self.min_spectral_width].into_dyn().into(),
        );
        grid_rec.insert(
            "w.max".to_string(),
            array![self.max_spectral_width].into_dyn().into(),
        );
        grid_rec.insert(
            "ve.min".to_string(),
            array![self.min_velocity_error].into_dyn().into(),
        );
        grid_rec.insert(
            "ve.max".to_string(),
            array![self.max_velocity_error].into_dyn().into(),
        );
        if !magnetic_lat.is_empty() {
            grid_rec.insert(
                "vector.mlat".to_string(),
                Array::from_vec(magnetic_lat).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.mlon".to_string(),
                Array::from_vec(magnetic_lon).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.kvect".to_string(),
                Array::from_vec(azimuth).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.srng".to_string(),
                Array::from_vec(slant_range).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.stid".to_string(),
                Array::from_vec(station_ids).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.channel".to_string(),
                Array::from_vec(channels).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.index".to_string(),
                Array::from_vec(index).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.vel.median".to_string(),
                Array::from_vec(velocity_median).into_dyn().into(),
            );
            grid_rec.insert(
                "vector.vel.sd".to_string(),
                Array::from_vec(velocity_stddev).into_dyn().into(),
            );
            if extended_flag {
                grid_rec.insert(
                    "vector.pwr.median".to_string(),
                    Array::from_vec(power_median).into_dyn().into(),
                );
                grid_rec.insert(
                    "vector.pwr.sd".to_string(),
                    Array::from_vec(power_stddev).into_dyn().into(),
                );
                grid_rec.insert(
                    "vector.wdt.median".to_string(),
                    Array::from_vec(spectral_width_median).into_dyn().into(),
                );
                grid_rec.insert(
                    "vector.wdt.sd".to_string(),
                    Array::from_vec(spectral_width_stddev).into_dyn().into(),
                );
            }
        }
        GridRecord::new(&mut grid_rec).map_err(|e| e.into())
    }
}
