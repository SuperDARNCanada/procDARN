use crate::error::BackscatterError;
use crate::gridding::grid::GridError;
use crate::utils::dmap::convert_to_dmapvec;
use crate::utils::hdw::HdwInfo;
use crate::utils::rpos::{rpos_inv_mag, rpos_range_beam_azimuth_elevation};
use crate::utils::scan::{RadarBeam, RadarScan};
use chrono::NaiveDateTime;
use dmap::formats::GridRecord;
use dmap::DmapType;
use std::f64::consts::PI;
use std::iter;

pub const GRID_REVISION_MAJOR: i32 = 2;
pub const GRID_REVISION_MINOR: i32 = 0;
pub const VELOCITY_ERROR_MIN: f64 = 100.0; // m/s
pub const POWER_LIN_ERROR_MIN: f64 = 1.0; // a.u. in linear scale
pub const WIDTH_LIN_ERROR_MIN: f64 = 1.0; // m/s

pub const RADIUS_EARTH: f64 = 6371.2; // km

#[derive(Debug, Default)]
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

#[derive(Debug, Default)]
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

#[derive(Debug, Default)]
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
    pub fn find_point(&self, reference: i32) -> Result<usize, GridError> {
        self.points
            .iter()
            .position(|x| x.reference == reference)
            .ok_or(GridError::Message(format!(
                "Point {} not in grid table",
                reference
            )))
    }

    /// Adds a grid beam to the grid table.
    /// Called GridTableAddBeam in RST
    pub fn add_beam(
        &mut self,
        hdw: &HdwInfo,
        altitude: f64,
        time: f64,
        scan_beam: &RadarBeam,
        chisham: bool,
        old_aacgm: bool,
    ) -> Result<usize, BackscatterError> {
        let velocity_correction: f64 = (2.0 * PI / 86400.0)
            * RADIUS_EARTH
            * 1000.0
            * (PI * hdw.latitude.clone() as f64 / 180.0).cos();
        self.num_beams += 1;

        let mut grid_beam = GridBeam {
            beam: scan_beam.beam,
            first_range: scan_beam.first_range,
            range_sep: scan_beam.range_sep,
            rx_rise: scan_beam.rx_rise,
            num_ranges: scan_beam.num_ranges,
            ..Default::default()
        };

        // TODO: Convert tval to year, month, day, hour, minute, seconds

        for range in 0..grid_beam.num_ranges {
            // Calculate geographic azimuth and elevation to scatter point
            let (azimuth_geo, elevation_geo) = rpos_range_beam_azimuth_elevation(
                grid_beam.beam,
                range,
                year,
                hdw,
                first_range,
                range_sep,
                rx_rise,
                altitude,
                chisham,
            )?;

            // Calculate magnetic latitude, longitude, and azimuth of scatter point
            let (mag_lat, mut mag_lon, mut azimuth_mag) = rpos_inv_mag(
                grid_beam.beam,
                range,
                year,
                hdw,
                first_range,
                range_sep,
                rx_rise,
                altitude,
                chisham,
                old_aacgm,
            )?;

            // Ensure magnetic azimuth and longitude between 0-360 degrees
            if azimuth_mag < 0.0 {
                azimuth_mag += 360.0;
            }
            if mag_lon < 0.0 {
                mag_lon += 360.0;
            }

            // Calculate magnetic grid cell latitude, (e.g. 72.1->72.5, 57.8->57.5, etc)
            let grid_lat: f64;
            if mag_lat > 0.0 {
                grid_lat = mag_lat.floor() + 0.5;
            } else {
                grid_lat = mag_lat.floor() - 0.5;
            }

            // Calculate magnetic grid longitude spacing at grid latitude
            let lon_spacing = (360.0 * (grid_lat.abs() * PI / 180.0).cos() + 0.5).floor() / 360.0;

            // Calculate magnetic grid cell longitude
            let _grid_lon = (mag_lon * lon_spacing + 0.5) / lon_spacing;

            // Calculate reference number for cell
            let reference: i32;
            if mag_lat > 0.0 {
                reference = (1000.0 * mag_lat.floor() + (mag_lon * lon_spacing).floor()) as i32;
            } else {
                reference =
                    (-1000.0 * (-1.0 * mag_lat).floor() - (mag_lon * lon_spacing).floor()) as i32;
            }

            // Find GridPoint corresponding to reference number for cell, make new GridPoint if none found
            let index = self.find_point(reference)?;
            let mut point = &mut self.points[index];

            // Update the total number of range gates that map to GridPoint (GridPoint.max)
            point.reference = reference;

            // Set magnetic lat/lon for GridPoint
            point.magnetic_lat = mag_lat;
            point.magnetic_lon = mag_lon;

            // Set index, magnetic azimuth, inertial velocity correction factor of beam
            grid_beam.index[range as usize] = index as i32;
            grid_beam.azimuth[range as usize] = azimuth_mag;
            grid_beam.ival[range as usize] =
                velocity_correction * (PI * (azimuth_geo + 90.0) / 180.0).cos();
        }
        // Return index of beam number added to self
        Ok((self.num_beams - 1) as usize)
    }

    /// Find the index of the beam in the grid table whose beam number and operating parameters
    /// match those of the input.
    /// Called GridTableFindBeam in RST
    pub fn find_beam(&self, beam: &RadarBeam) -> Result<usize, GridError> {
        self.beams
            .iter()
            .position(|x| {
                x.beam == beam.beam
                    && x.first_range == beam.first_range
                    && x.range_sep == beam.range_sep
                    && x.num_ranges == beam.num_ranges
            })
            .ok_or(GridError::Message(format!("Beam not found in grid table",)))
    }

    /// Maps radar scan data to an equal-area grid in magnetic coordinates.
    /// Called GridTableMap in RST
    pub fn map(
        &mut self,
        scan: &RadarScan,
        hdw: &HdwInfo,
        tlen: i32,
        iflg: bool,
        altitude: f64,
        chisham: bool,
        old_aacgm: bool,
    ) -> Result<(), GridError> {
        let time = (&scan.start_time + &scan.end_time) / 2.0;
        if self.status == 0 {
            self.status = 1;
            self.noise_mean = 0.0;
            self.noise_stddev = 0.0;
            self.freq = 0.0;
            self.num_scans = 0;
            self.start_time = scan.start_time.clone();
            self.end_time = scan.start_time.clone() + tlen;
            self.station_id = scan.station_id.clone();
        }

        for scan_beam in scan.beams.iter() {
            let mut beam_index: usize = 0;
            if scan_beam.beam != -1 {
                beam_index = match self.find_beam(scan_beam) {
                    Ok(i) => i,
                    Err(_) => self.add_beam(hdw, altitude, time, scan_beam, chisham, old_aacgm)?,
                };
            }

            let grid_beam = &self.beams[beam_index];

            for range in 0..scan_beam.num_ranges.clone() as usize {
                if scan_beam.scatter[range] == 0 {
                    continue;
                }

                let mut velocity_error = scan_beam.cells[range].velocity_error;
                let mut power_lin_error = scan_beam.cells[range].power_lin_error;
                let mut width_lin_error = scan_beam.cells[range].spectral_width_lin_error;

                if velocity_error < VELOCITY_ERROR_MIN {
                    velocity_error = VELOCITY_ERROR_MIN;
                }
                if power_lin_error < POWER_LIN_ERROR_MIN {
                    power_lin_error = POWER_LIN_ERROR_MIN;
                }
                if width_lin_error < WIDTH_LIN_ERROR_MIN {
                    width_lin_error = WIDTH_LIN_ERROR_MIN;
                }

                // Get grid cell of radar beam/gate measurement
                let mut grid_cell = &mut self.points[grid_beam.index[range] as usize];

                // Add magnetic azimuth of radar beam/gate measurement
                grid_cell.azimuth += grid_beam.azimuth[range];

                if iflg {
                    grid_cell.velocity_median_north -= (scan_beam.cells[range].velocity
                        + grid_beam.ival[range])
                        * (grid_beam.azimuth[range] * PI / 180.).cos()
                        / (velocity_error * velocity_error);
                    grid_cell.velocity_median_east -= (scan_beam.cells[range].velocity
                        + grid_beam.ival[range])
                        * (grid_beam.azimuth[range] * PI / 180.).sin()
                        / (velocity_error * velocity_error);
                } else {
                    grid_cell.velocity_median_north -= scan_beam.cells[range].velocity
                        * (grid_beam.azimuth[range] * PI / 180.).cos()
                        / (velocity_error * velocity_error);
                    grid_cell.velocity_median_east -= scan_beam.cells[range].velocity
                        * (grid_beam.azimuth[range] * PI / 180.).sin()
                        / (velocity_error * velocity_error);
                }

                grid_cell.power_median +=
                    scan_beam.cells[range].power_lin / (power_lin_error * power_lin_error);
                grid_cell.spectral_width_median +=
                    scan_beam.cells[range].spectral_width_lin / (width_lin_error * width_lin_error);

                grid_cell.velocity_stddev /= velocity_error * velocity_error;
                grid_cell.power_stddev /= power_lin_error * power_lin_error;
                grid_cell.spectral_width_stddev /= width_lin_error * width_lin_error;
                grid_cell.count += 1;
            }
        }

        // TODO: Check if somehow all beams in scan not considered?

        let mut freq: f64 = 0.0;
        let mut noise: f64 = 0.0;
        let mut variance: f64 = 0.0;
        let mut count: f64 = 0.0;

        for scan_beam in scan.beams.iter().filter(|beam| beam.beam != -1) {
            self.program_id = scan_beam.program_id;

            // Sum the frequency and noise values
            freq += scan_beam.freq as f64;
            noise += scan_beam.noise as f64;
            count += 1.0;
        }

        // Average frequency and noise over all beams in scan
        freq = freq / count;
        noise = noise / count;

        for scan_beam in scan.beams.iter().filter(|beam| beam.beam != -1) {
            variance += (scan_beam.noise as f64 - noise) * (scan_beam.noise as f64 - noise);
        }
        self.noise_mean += noise;
        self.noise_stddev += (variance / count).sqrt();
        self.freq += freq;
        self.num_scans += 1;

        Ok(())
    }

    /// Converts the GridTable to a GridRecord for writing to file.
    /// Equivalent to GridTableWrite in RST.
    pub fn to_dmap_record(&self) -> Result<GridRecord, GridError> {
        let start_time = NaiveDateTime::from_timestamp_micros(self.start_time * 1000.0 as i64)?;

        // Find the valid points in the grid
        let valid_points: Vec<&GridPoint> = self.points.iter().filter(|&p| p.count > 0).collect();
        let num_points = valid_points.len();

        // These vector fields require accessing the points of grid_table
        let magnetic_lat: Vec<DmapType::DOUBLE> =
            valid_points.iter().map(|&p| p.magnetic_lat).collect();
        let magnetic_lon: Vec<DmapType::DOUBLE> =
            valid_points.iter().map(|&p| p.magnetic_lon).collect();
        let azimuth: Vec<DmapType::DOUBLE> = valid_points.iter().map(|&p| p.azimuth).collect();
        let index: Vec<DmapType::INT> = valid_points.iter().map(|&p| p.reference).collect();
        let velocity_median: Vec<DmapType::DOUBLE> =
            valid_points.iter().map(|&p| p.velocity_median).collect();
        let velocity_stddev: Vec<DmapType::DOUBLE> =
            valid_points.iter().map(|&p| p.velocity_stddev).collect();
        let power_median: Vec<DmapType::DOUBLE> =
            valid_points.iter().map(|&p| p.power_median).collect();
        let power_stddev: Vec<DmapType::DOUBLE> =
            valid_points.iter().map(|&p| p.power_stddev).collect();
        let spectral_width_median: Vec<DmapType::DOUBLE> = valid_points
            .iter()
            .map(|&p| p.spectral_width_median)
            .collect();
        let spectral_width_stddev: Vec<DmapType::DOUBLE> = valid_points
            .iter()
            .map(|&p| p.spectral_width_stddev)
            .collect();
        let station_ids: Vec<DmapType::SHORT> = iter::repeat(self.station_id)
            .take(valid_points.len())
            .collect();
        let channels: Vec<DmapType::SHORT> = iter::repeat(self.channel)
            .take(valid_points.len())
            .collect();

        Ok(GridRecord {
            start_year: start_time.format("%Y").to_string().parse::<i16>()?,
            start_month: start_time.format("%m").to_string().parse::<i16>()?,
            start_day: start_time.format("%d").to_string().parse::<i16>()?,
            start_hour: start_time.format("%H").to_string().parse::<i16>()?,
            start_minute: start_time.format("%M").to_string().parse::<i16>()?,
            start_second: start_time.format("%S.%.6f").to_string().parse::<f64>()?,
            end_year: start_time.format("%Y").to_string().parse::<i16>()?,
            end_month: start_time.format("%m").to_string().parse::<i16>()?,
            end_day: start_time.format("%d").to_string().parse::<i16>()?,
            end_hour: start_time.format("%H").to_string().parse::<i16>()?,
            end_minute: start_time.format("%M").to_string().parse::<i16>()?,
            end_second: start_time.format("%S.%.6f").to_string().parse::<f64>()?,
            station_ids: convert_to_dmapvec(vec![self.station_id as i16]),
            channels: convert_to_dmapvec(vec![self.channel as i16]),
            num_vectors: convert_to_dmapvec(vec![num_points as i16]),
            freq: convert_to_dmapvec(vec![self.freq as f32]),
            grid_major_revision: convert_to_dmapvec(vec![GRID_REVISION_MAJOR as i16]),
            grid_minor_revision: convert_to_dmapvec(vec![GRID_REVISION_MINOR as i16]),
            program_ids: convert_to_dmapvec(vec![self.program_id as i16]),
            noise_mean: convert_to_dmapvec(vec![self.noise_mean as f32]),
            noise_stddev: convert_to_dmapvec(vec![self.noise_stddev as f32]),
            groundscatter: convert_to_dmapvec(vec![self.groundscatter as i16]),
            velocity_min: convert_to_dmapvec(vec![self.min_velocity as f32]),
            velocity_max: convert_to_dmapvec(vec![self.max_velocity as f32]),
            power_min: convert_to_dmapvec(vec![self.min_power as f32]),
            power_max: convert_to_dmapvec(vec![self.min_power as f32]),
            spectral_width_min: convert_to_dmapvec(vec![self.min_spectral_width as f32]),
            spectral_width_max: convert_to_dmapvec(vec![self.max_spectral_width as f32]),
            velocity_error_min: convert_to_dmapvec(vec![self.min_velocity_error as f32]),
            velocity_error_max: convert_to_dmapvec(vec![self.max_velocity_error as f32]),
            magnetic_lat: convert_to_dmapvec(magnetic_lat),
            magnetic_lon: convert_to_dmapvec(magnetic_lon),
            magnetic_azi: convert_to_dmapvec(azimuth),
            station_id_vector: convert_to_dmapvec(station_ids),
            channel_vector: convert_to_dmapvec(channels),
            grid_cell_index: convert_to_dmapvec(index),
            velocity_median: convert_to_dmapvec(velocity_median),
            velocity_stddev: convert_to_dmapvec(velocity_stddev),
            power_median: convert_to_dmapvec(power_median),
            power_stddev: convert_to_dmapvec(power_stddev),
            spectral_width_median: convert_to_dmapvec(spectral_width_median),
            spectral_width_stddev: convert_to_dmapvec(spectral_width_stddev),
        })
    }
}
