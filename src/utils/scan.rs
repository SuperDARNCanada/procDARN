use crate::error::BackscatterError;
use chrono::NaiveDate;
use dmap::formats::FitacfRecord;
use time::Month::December;

#[derive(Copy, Clone, Default)]
pub struct RadarCell {
    pub groundscatter: i32,            // gsct in RST
    pub power_lag_zero: f64,           // p_0 in RST
    pub power_error_lag_zero: f64,     // p_0_e in RST
    pub velocity: f64,                 // v in RST
    pub velocity_error: f64,           // v_e in RST
    pub spectral_width_lin: f64,       // w_l in RST
    pub spectral_width_lin_error: f64, // w_l_e in RST
    pub power_lin: f64,                // p_l in RST
    pub power_lin_error: f64,          // p_l_e in RST
    pub phi_zero: f64,                 // phi0 in RST
    pub elevation: f64,                // elv in RST
}

#[derive(Clone, Default)]
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
    pub scatter: Vec<u8>,         // sct in RST
    pub cells: Vec<RadarCell>,    // rng in RST
}
impl RadarBeam {
    pub fn reset(&mut self) {
        self.scatter.clear();
        self.cells.clear();
    }
}

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

    /// Keeps beams whose number have index of beam_num in beam_list.
    /// Called RadarScanResetBeam in RST
    pub fn reset_beams(
        &mut self,
        beam_list: Vec<i32>,
        beam_num: i32,
    ) -> Result<(), BackscatterError> {
        if beam_num == 0 {
            Err(BackscatterError {
                details: format!("Bad beam_num {}, must be greater than zero", beam_num),
            })
        } else {
            // keep beams from self.beams that are in beam_list at position beam_num
            self.beams = self
                .beams
                .clone()
                .into_iter()
                .filter(|beam| {
                    beam_list
                        .iter()
                        .position(|&x| x == beam.beam)
                        .is_some_and(|x| x as i32 == beam_num)
                })
                .collect();
            Ok(())
        }
    }

    /// Called RadarScanAddBeam in RST
    pub fn add_beam(&mut self, num_ranges: i32) {
        self.beams.push(RadarBeam {
            num_ranges,
            ..Default::default()
        })
    }

    /// Exclude beams that are not part of a scan
    pub fn exclude_outofscan(&mut self) {
        self.beams = self
            .beams
            .clone()
            .into_iter()
            .filter(|beam| beam.scan >= 0)
            .collect();
    }

    /// Read a full scan of data from a vector of FitacfRecords. If scan_length is Some(x), will
    /// grab the first records spanning x seconds. Otherwise, uses the scan flag in the FitacfRecords
    /// to determine the end of the scan.
    /// Called FitReadRadarScan in RST.
    pub fn get_first_scan(
        fit_records: &Vec<FitacfRecord>,
        scan_length: Option<u32>,
    ) -> Result<RadarScan, BackscatterError> {
        for i in 0..fit_records.len() {
            let rec = &fit_records[i];
            let mut scan_rec = RadarScan {
                station_id: rec.station_id as i32,
                version_major: rec.radar_revision_major as i32,
                version_minor: rec.radar_revision_minor as i32,
                start_time: NaiveDate::from_ymd_opt(
                    rec.year as i32,
                    rec.month as u32,
                    rec.day as u32,
                )?
                .and_hms_opt(rec.hour as u32, rec.minute as u32, rec.second as u32)?
                .timestamp()
                    + (rec.microsecond as f64) / 1e6,
                ..Default::default()
            };
            let mut beam = RadarBeam {
                time: NaiveDate::from_ymd_opt(rec.year as i32, rec.month as u32, rec.day as u32)?
                    .and_hms_opt(rec.hour as u32, rec.minute as u32, rec.second as u32)?
                    .timestamp()
                    + (rec.microsecond as f64) / 1e6,
                scan: rec.scan_flag as i32,
                beam: rec.beam_num as i32,
                beam_azimuth: rec.beam_azimuth,
                program_id: rec.control_program as i32,
                integration_time_s: rec.intt_second as i32,
                integration_time_us: rec.intt_microsecond,
                num_averages: rec.num_averages as i32,
                first_range: rec.first_range as i32,
                range_sep: rec.range_sep as i32,
                rx_rise: rec.rx_rise_time as i32,
                freq: rec.tx_freq as i32,
                noise: rec.search_noise as i32,
                attenuation: rec.attenuation as i32,
                channel: rec.channel as i32,
                num_ranges: rec.num_ranges as i32,
                ..Default::default()
            };
            for r in 0..beam.num_ranges {
                beam.scatter.push(rec.quality_flag.clone().collect());

                // Create a new measurement (RadarCell) and populate it
                let mut cell = RadarCell {
                    groundscatter: rec.ground_flag[r],
                    power_lag_zero: rec.lag_zero_power[r],
                    power_error_lag_zero: 0.0,
                    velocity: rec.velocity[r],
                    power_lin: rec.lambda_power[r],
                    spectral_width_lin: rec.lambda_spectral_width[r],
                    velocity_error: rec.velocity_error[r],
                    ..Default::default()
                };
                if let Some(x) = rec.lag_zero_phi.clone() {
                    cell.phi_zero = x[r]
                } else {
                    cell.phi_zero = 0.0
                }
                if let Some(x) = rec.elevation.clone() {
                    cell.elevation = x[r]
                } else {
                    cell.elevation = 0.0
                }

                // Add the measurement (RadarCell) to the beam
                beam.cells.push(cell);
            }

            // Update the end time of the scan
            scan_rec.end_time =
                NaiveDate::from_ymd_opt(rec.year as i32, rec.month as u32, rec.day as u32)?
                    .and_hms_opt(rec.hour as u32, rec.minute as u32, rec.second as u32)?
                    .timestamp()
                    + (rec.microsecond as f64) / 1e6;

            // Conditions for finding the end of the scan
            match scan_length {
                // If the scan has spanned longer than scan_length
                Some(x) => {
                    if scan_rec.end_time - scan_rec.start_time >= x as f64 {
                        break;
                    }
                }
                // If the next record is the start of a new scan
                None => {
                    if i < fit_records.len() - 1 && fit_records[i + 1].scan_flag.abs() == 1 {
                        break;
                    }
                }
            }
        }
        Ok(scan_rec)
    }
}
