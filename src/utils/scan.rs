use crate::gridding::grid_table::GridTable;
use crate::utils::rpos::slant_range;
use crate::utils::sugar::get_datetime;
use dmap::formats::fitacf::FitacfRecord;
use numpy::ndarray::ArrayD;
use crate::error::ProcdarnError;
use crate::gridding::grid::GridError;

#[derive(Copy, Clone, Default, PartialEq)]
pub struct RadarCell {
    pub groundscatter: i8,             // gsct in RST
    pub power_lag_zero: f32,           // pwr0 in RST
    pub power_error_lag_zero: f32,     // pwr0_e in RST
    pub velocity: f32,                 // v in RST
    pub velocity_error: f32,           // v_e in RST
    pub spectral_width_lin: f32,       // w_l in RST
    pub spectral_width_lin_error: f32, // w_l_e in RST
    pub power_lin: f32,                // p_l in RST
    pub power_lin_error: f32,          // p_l_e in RST
    pub phi_zero: f32,                 // phi0 in RST
    pub elevation: f32,                // elv in RST
}

#[derive(Clone, Default, PartialEq)]
pub struct RadarBeam {
    pub scan: i32,                // scan in RST
    pub beam: i32,                // bm in RST
    pub beam_azimuth: f32,        // bmazm in RST
    pub time: f64,                // time in RST
    pub program_id: i32,          // cpid in RST
    pub integration_time_s: i32,  // intt.sc in RST
    pub integration_time_us: i32, // intt.us in RST
    pub num_averages: i32,        // nave in RST
    pub first_range: i32,         // frang in RST
    pub range_sep: i32,           // rsep in RST
    pub rx_rise: i32,             // rxrise in RST
    pub freq: i32,                // freq in RST
    pub noise: i32,               // noise in RST
    pub attenuation: i32,         // atten in RST
    pub channel: i32,             // channel in RST
    pub num_ranges: i32,          // nrang in RST
    pub scatter: Vec<i8>,         // sct in RST
    pub cells: Vec<RadarCell>,    // rng in RST
}
impl RadarBeam {
    pub fn reset(&mut self) {
        self.scatter.clear();
        self.cells.clear();
    }
}

#[derive(Clone, Default, PartialEq)]
pub struct RadarScan {
    pub station_id: i32,       // stid in RST
    pub version_major: i32,    // version.major in RST
    pub version_minor: i32,    // version.minor in RST
    pub start_time: f64,       // st_time in RST
    pub end_time: f64,         // ed_time in RST
    pub beams: Vec<RadarBeam>, // bm in RST
}
impl RadarScan {
    /// Clears the beams
    /// Called RadarScanReset in RST
    pub fn reset(&mut self) {
        self.beams.clear();
    }

    /// Remove beams whose beam number is in beam_list
    /// Called RadarScanResetBeam in RST
    pub fn reset_beams(&mut self, beam_list: &Vec<i32>) -> Result<(), GridError> {
        // remove beams from self.beams that are in beam_list
        self.beams = self
            .beams
            .clone()
            .into_iter()
            .filter(|beam| !beam_list.contains(&beam.beam))
            .collect();
        Ok(())
    }

    /// Called RadarScanAddBeam in RST
    pub fn add_beam(&mut self, num_ranges: i32) {
        self.beams.push(RadarBeam {
            num_ranges,
            ..Default::default()
        })
    }

    /// Exclude beams that are not part of a scan.
    /// Called exclude_outofscan in radarscan.c of RST.
    pub fn exclude_outofscan(&mut self) {
        self.beams = self
            .beams
            .clone()
            .into_iter()
            .filter(|beam| beam.scan >= 0)
            .collect();
    }

    /// Read a full scan of data from a vector of FitacfRecords. If scan_length is `Some(x)`, will
    /// grab the first records spanning `x` seconds. Otherwise, uses the scan flag in the FitacfRecords
    /// to determine the end of the scan.
    /// Called FitReadRadarScan in fitscan.c of RST.
    pub fn get_first_scan(
        fit_records: &[FitacfRecord],
        scan_length: Option<u32>,
    ) -> Result<RadarScan, ProcdarnError> {
        if fit_records.len() == 0 {
            return Err(ProcdarnError::ZeroRecords("in `get_first_scan()`"));
        }
        let mut rec = &fit_records[0];
        let mut scan: RadarScan = RadarScan {
            station_id: i32::try_from(
                rec.get(&"stid".to_string())
                    .ok_or(ProcdarnError::MissingField("stid"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("stid"))?,
            version_major: i32::try_from(
                rec.get(&"radar.revision.major".to_string())
                    .ok_or(ProcdarnError::MissingField("radar.revision.major"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("radar.revision.major"))?,
            version_minor: i32::try_from(
                rec.get(&"radar.revision.minor".to_string())
                    .ok_or(ProcdarnError::MissingField("radar.revision.minor"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("radar.revision.minor"))?,
            start_time: get_datetime(rec)?.timestamp_micros() as f64,
            ..Default::default()
        };

        for i in 0..fit_records.len() {
            rec = &fit_records[i];

            let mut beam = RadarBeam {
                time: get_datetime(rec)?.timestamp_micros() as f64,
                scan: i32::try_from(
                    rec.get(&"scan".to_string())
                        .ok_or(ProcdarnError::MissingField("scan"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("scan"))?,
                beam: i32::try_from(
                    rec.get(&"bmnum".to_string())
                        .ok_or(ProcdarnError::MissingField("bmnum"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("bmnum"))?,
                beam_azimuth: f32::try_from(
                    rec.get(&"bmazm".to_string())
                        .ok_or(ProcdarnError::MissingField("bmazm"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("bmazm"))?,
                program_id: i32::try_from(
                    rec.get(&"cp".to_string())
                        .ok_or(ProcdarnError::MissingField("cp"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("cp"))?,
                integration_time_s: i32::try_from(
                    rec.get(&"intt.sc".to_string())
                        .ok_or(ProcdarnError::MissingField("intt.sc"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("intt.sc"))?,
                integration_time_us: i32::try_from(
                    rec.get(&"intt.us".to_string())
                        .ok_or(ProcdarnError::MissingField("intt.us"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("intt.us"))?,
                num_averages: i32::try_from(
                    rec.get(&"nave".to_string())
                        .ok_or(ProcdarnError::MissingField("nave"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("nave"))?,
                first_range: i32::try_from(
                    rec.get(&"frang".to_string())
                        .ok_or(ProcdarnError::MissingField("frang"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("frang"))?,
                range_sep: i32::try_from(
                    rec.get(&"rsep".to_string())
                        .ok_or(ProcdarnError::MissingField("rsep"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("rsep"))?,
                rx_rise: i32::try_from(
                    rec.get(&"rxrise".to_string())
                        .ok_or(ProcdarnError::MissingField("rxrise"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("rxrise"))?,
                freq: i32::try_from(
                    rec.get(&"tfreq".to_string())
                        .ok_or(ProcdarnError::MissingField("tfreq"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("tfreq"))?,
                noise: i32::try_from(
                    rec.get(&"noise.search".to_string())
                        .ok_or(ProcdarnError::MissingField("noise.search"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("noise.search"))?,
                attenuation: i32::try_from(
                    rec.get(&"atten".to_string())
                        .ok_or(ProcdarnError::MissingField("atten"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("atten"))?,
                channel: i32::try_from(
                    rec.get(&"channel".to_string())
                        .ok_or(ProcdarnError::MissingField("channel"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("channel"))?,
                num_ranges: i32::try_from(
                    rec.get(&"nrang".to_string())
                        .ok_or(ProcdarnError::MissingField("nrang"))?
                        .clone(),
                )
                .map_err(|_| ProcdarnError::WrongType("nrang"))?,
                ..Default::default()
            };
            let groundscatter = ArrayD::try_from(
                rec.get(&"gflg".to_string())
                    .ok_or(ProcdarnError::MissingField("gflg"))?
                    .clone(),
            )
            .map_err(|e| ProcdarnError::Channel(format!("gflg - {e}")))?;
            let power_lag_zero = ArrayD::try_from(
                rec.get(&"pwr0".to_string())
                    .ok_or(ProcdarnError::MissingField("pwr0"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("pwr0"))?;
            let power_error_lag_zero = 0.0;
            let velocity = ArrayD::try_from(
                rec.get(&"v".to_string())
                    .ok_or(ProcdarnError::MissingField("v"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("v"))?;
            let power_lin = ArrayD::try_from(
                rec.get(&"p_l".to_string())
                    .ok_or(ProcdarnError::MissingField("p_l"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("p_l"))?;
            let spectral_width_lin = ArrayD::try_from(
                rec.get(&"w_l".to_string())
                    .ok_or(ProcdarnError::MissingField("w_l"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("w_l"))?;
            let velocity_error = ArrayD::try_from(
                rec.get(&"v_e".to_string())
                    .ok_or(ProcdarnError::MissingField("v_e"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("v_e"))?;
            let quality_flag = ArrayD::try_from(
                rec.get(&"qflg".to_string())
                    .ok_or(ProcdarnError::MissingField("qflg"))?
                    .clone(),
            )
            .map_err(|_| ProcdarnError::WrongType("qflg"))?;
            let slist = ArrayD::try_from(
                rec.get(&"slist".to_string())
                    .ok_or(ProcdarnError::MissingField("slist"))?
                    .clone(),
            )
                .map_err(|_| ProcdarnError::WrongType("slist"))?;
            let lag_zero_phi = rec.get(&"phi0".to_string());
            let elevation = rec.get(&"elv".to_string());
            for r in 0..beam.num_ranges as usize {
                let slist_idx = slist.iter().position(|&x: &i16| x as usize == r);
                let cell = match slist_idx {
                    Some(idx) => {
                        beam.scatter.push(quality_flag[idx]);
                        let mut cell_tmp = RadarCell {
                            groundscatter: groundscatter[idx],
                            power_lag_zero: power_lag_zero[idx],
                            power_error_lag_zero: power_error_lag_zero.clone(),
                            velocity: velocity[idx],
                            power_lin: power_lin[idx],
                            spectral_width_lin: spectral_width_lin[idx],
                            velocity_error: velocity_error[idx],
                            ..Default::default()
                        };
                        if let Some(x) = lag_zero_phi {
                            cell_tmp.phi_zero = ArrayD::try_from(x.clone())
                                .map_err(|_| ProcdarnError::WrongType("phi0"))?[idx]
                        }
                        if let Some(x) = elevation {
                            cell_tmp.elevation = ArrayD::try_from(x.clone())
                                .map_err(|_| ProcdarnError::WrongType("elv"))?[idx]
                        }
                        cell_tmp
                    },
                    None => {
                        beam.scatter.push(0);
                        RadarCell::default()
                    },
                };

                // Add the measurement (RadarCell) to the beam
                beam.cells.push(cell);
            }

            // Add the beam to the scan
            scan.beams.push(beam);

            // Update the end time of the scan
            scan.end_time = get_datetime(rec)?.timestamp_micros() as f64;

            // Conditions for finding the end of the scan
            match scan_length {
                // If the scan has spanned longer than scan_length
                Some(x) => {
                    if scan.end_time - scan.start_time >= x as f64 {
                        break;
                    }
                }
                // If the next record is the start of a new scan
                None => {
                    if i < fit_records.len() - 1
                        && i8::try_from(
                            fit_records[i + 1]
                                .get(&"scan".to_string())
                                .ok_or(ProcdarnError::MissingField("scan"))?
                                .clone(),
                        )
                        .map_err(|_| ProcdarnError::WrongType("scan"))?
                        .abs()
                            == 1
                    {
                        break;
                    }
                }
            }
        }
        Ok(scan)
    }

    /// Filters data in the scan based on optional min and max range gates or slant ranges.
    /// Called exclude_range in make_grid.c of RST.
    pub fn exclude_range(
        &mut self,
        min_range_gate: Option<usize>,
        max_range_gate: Option<usize>,
        min_slant_range: Option<f32>,
        max_slant_range: Option<f32>,
    ) {
        let range_edge = 0;
        for beam in self.beams.iter_mut().filter(|b| b.beam != -1) {
            // If either min or max slant range given, then exclude data using slant range filters
            if min_slant_range.is_some() || max_slant_range.is_some() {
                for rg in 0..beam.num_ranges {
                    let range_slant = slant_range(
                        beam.first_range,
                        beam.range_sep,
                        beam.rx_rise,
                        range_edge,
                        rg,
                    );
                    match (min_slant_range, max_slant_range) {
                        (Some(min), Some(max)) => {
                            if min > range_slant || range_slant > max {
                                beam.scatter[rg as usize] = 0;
                            }
                        }
                        (Some(min), None) => {
                            if min > range_slant {
                                beam.scatter[rg as usize] = 0;
                            }
                        }
                        (None, Some(max)) => {
                            if range_slant > max {
                                beam.scatter[rg as usize] = 0;
                            }
                        }
                        (None, None) => {}
                    }
                }
            } else {
                // Exclude data using range gate filters
                match (min_range_gate, max_range_gate) {
                    (Some(min), Some(max)) => {
                        for scat in beam.scatter[..min].iter_mut() {
                            *scat = 0;
                        }
                        for scat in beam.scatter[max..].iter_mut() {
                            *scat = 0;
                        }
                    }
                    (Some(min), None) => {
                        for scat in beam.scatter[..min].iter_mut() {
                            *scat = 0;
                        }
                    }
                    (None, Some(max)) => {
                        for scat in beam.scatter[max..].iter_mut() {
                            *scat = 0;
                        }
                    }
                    (None, None) => {}
                }
            }
        }
    }

    /// Excludes ground scatter from the radar scan.
    /// Called FilterBoundType in bound.c of RST
    pub fn exclude_groundscatter(&mut self) {
        for beam in self.beams.iter_mut() {
            for rg in 0..beam.num_ranges as usize {
                if beam.scatter[rg] == 0 {
                    continue;
                }
                if beam.cells[rg].groundscatter == 1 {
                    beam.scatter[rg] = 0;
                }
            }
        }
    }

    /// Excludes ionospheric scatter from the radar scan.
    /// Called FilterBoundType in bound.c of RST
    pub fn exclude_ionospheric_scatter(&mut self) {
        for beam in self.beams.iter_mut() {
            for rg in 0..beam.num_ranges as usize {
                if beam.scatter[rg] == 0 {
                    continue;
                }
                if beam.cells[rg].groundscatter == 0 {
                    beam.scatter[rg] = 0;
                }
            }
        }
    }

    pub fn exclude_outofbounds(&mut self, grid_table: &GridTable) {
        for beam in self.beams.iter_mut() {
            for rg in 0..beam.num_ranges as usize {
                if beam.scatter[rg] == 0 {
                    continue;
                }
                let cell = beam.cells[rg];
                let discard_cell = cell.velocity.abs() < grid_table.min_velocity
                    || cell.velocity.abs() > grid_table.max_velocity
                    || cell.power_lin < grid_table.min_power
                    || cell.power_lin > grid_table.max_power
                    || cell.spectral_width_lin < grid_table.min_spectral_width
                    || cell.spectral_width_lin > grid_table.max_spectral_width
                    || cell.velocity_error < grid_table.min_velocity_error
                    || cell.velocity_error > grid_table.max_velocity_error;
                if discard_cell {
                    beam.scatter[rg] = 0;
                }
            }
        }
    }
}
