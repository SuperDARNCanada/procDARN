use crate::fitting::fitacf3::fitstruct::{FitType, LagNode, RangeNode};
use crate::fitting::fitacf3::least_squares::LeastSquares;
use crate::hdw::hdw::HdwInfo;
use chrono::NaiveDateTime;
use dmap::formats::{FitacfRecord, RawacfRecord};
use dmap::{DmapVec, InDmap};
use std::error::Error;
use std::f32::consts::PI as PI_f32;
use std::f64::consts::PI;
use std::fmt;
use std::fmt::Display;
use std::iter::zip;
use crate::fitting::fitacf3::filtering;

type Result<T> = std::result::Result<T, Fitacf3Error>;

pub const FITACF_REVISION_MAJOR: i32 = 3;
pub const FITACF_REVISION_MINOR: i32 = 0;
pub const V_MAX: f32 = 30.0;
pub const W_MAX: f32 = 90.0;
pub const FLUCTUATION_CUTOFF_COEFFICIENT: f32 = 2.0;
pub const ALPHA_CUTOFF: f32 = 2.0;
pub const ACF_SNR_CUTOFF: f32 = 1.0;
pub const MIN_LAGS: i16 = 3;

#[derive(Debug, Clone)]
pub enum Fitacf3Error {
    Message(String),
    Lookup(String),
    Mismatch { msg: String },
}

impl Error for Fitacf3Error {}

impl Display for Fitacf3Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Fitacf3Error::Message(msg) => write!(f, "{}", msg),
            Fitacf3Error::Lookup(msg) => write!(f, "{}", msg),
            Fitacf3Error::Mismatch{msg} => write!(f, "{}", msg)
        }
    }
}

/// Creates the lag table based on the data.
fn create_lag_list(record: &RawacfRecord) -> Vec<LagNode> {
    let lag_table = &record.lag_table;
    let pulse_table = &record.pulse_table;
    let multi_pulse_increment = record.multi_pulse_increment;
    let sample_separation = record.sample_separation;

    let mut lags = vec![];
    for i in 0..record.num_lags as usize {
        let mut pulse_1_idx = 0;
        let mut pulse_2_idx = 0;
        let number = lag_table.data[2 * i + 1] - lag_table.data[2 * i]; // flattened, we want row i, cols 1 and 0
        for j in 0..record.num_pulses as usize {
            if lag_table.data[2 * i] == pulse_table.data[j] {
                pulse_1_idx = j;
            }
            if lag_table.data[2 * i + 1] == pulse_table.data[j] {
                pulse_2_idx = j;
            }
        }
        let sample_base_1 =
            (lag_table.data[2 * i] * (multi_pulse_increment / sample_separation)) as i32;
        let sample_base_2 =
            (lag_table.data[2 * i + 1] * (multi_pulse_increment / sample_separation)) as i32;
        lags.push(LagNode {
            lag_num: number as i32,
            pulses: [pulse_1_idx, pulse_2_idx],
            lag_idx: 0,
            sample_base_1,
            sample_base_2,
        });
    }
    lags
}

pub fn fit_rawacf_record(record: &RawacfRecord) -> Result<FitacfRecord> {
    let lags = create_lag_list(record);

    let noise_power;
    if record.num_averages <= 0 {
        noise_power = 1.0;
    } else {
        noise_power = acf_cutoff_power(record);
    }
    let mut range_list = vec![];
    for i in 0..record.range_list.data.len() {
        let range_num = record.range_list.data[i];
        if record.lag_zero_power.data[range_num as usize] != 0.0 {
            range_list.push(RangeNode::new(i, range_num as usize, record, &lags)?)
        }
    }
    filtering::filter_tx_overlapped_lags(record, lags, &mut range_list);
    filtering::filter_infinite_lags(&mut range_list);
    filtering::filter_low_power_lags(record, &mut range_list);
    filtering::filter_bad_acfs(record, &mut range_list, noise_power);
    acf_power_fitting(&mut range_list)?;
    calculate_phase_and_elev_sigmas(&mut range_list, record)?;
    // something wrong with phi values (sigmas are good, but phi is not)
    acf_phase_unwrap(&mut range_list);
    acf_phase_fitting(&mut range_list)?;
    filtering::filter_bad_fits(&mut range_list)?;
    xcf_phase_unwrap(&mut range_list)?;
    xcf_phase_fitting(&mut range_list)?;

    let dets = determinations(record, range_list, noise_power);
    dets
}

/// Passing
fn cutoff_power_correction(rec: &RawacfRecord) -> f64 {
    let std_dev = 1.0 / (rec.num_averages as f64).sqrt();

    let mut i = 0.0;
    let mut cumulative_pdf = 0.0;
    let mut cumulative_pdf_x_norm_power = 0.0;
    let mut normalized_power;
    while cumulative_pdf < (10.0 / rec.num_ranges as f64) {
        // Normalized power for calculating model PDF (Gaussian)
        normalized_power = i / 1000.0;
        let x = -(normalized_power - 1.0) * (normalized_power - 1.0) / (2.0 * std_dev * std_dev);
        let pdf = x.exp() / std_dev / (2.0 * PI).sqrt() / 1000.0;
        cumulative_pdf += pdf;

        // Cumulative value of PDF * x  -> needed for calculating the mean
        cumulative_pdf_x_norm_power += pdf * normalized_power;
        i += 1.0;
    }
    // Correcting factor as the inverse of a normalized mean
    cumulative_pdf / cumulative_pdf_x_norm_power
}

/// Calculates the minimum power value for ACFs in the record (passing)
fn acf_cutoff_power(rec: &RawacfRecord) -> f32 {
    let mut sorted_power_levels = rec.lag_zero_power.data.clone();
    sorted_power_levels.sort_by(|a, b| a.total_cmp(&b)); // sort floats
    let mut i: usize = 0;
    let mut j: f64 = 0.0;
    let mut min_power: f64 = 0.0;
    while j < 10.0 && i < rec.num_ranges as usize / 3 {
        if sorted_power_levels[i] > 0.0 {
            j += 1.0;
        }
        min_power += sorted_power_levels[i] as f64;
        i += 1;
    }
    if j <= 0.0 {
        j = 1.0;
    }
    min_power *= cutoff_power_correction(rec) / j;
    if min_power < 1.0 && rec.search_noise > 0.0 {
        min_power = rec.search_noise as f64;
    }
    min_power as f32
}

/// passing
fn acf_power_fitting(ranges: &mut Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);

    for mut range in ranges {
        let log_powers = &range.powers.ln_power;
        let sigmas = &range.powers.std_dev;
        let t = &range.powers.t;
        let num_points = range.powers.ln_power.len();
        if t.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::Message(format!(
                "Cannot perform acf power fitting - dimension mismatch"
            )))?
        }
        range.lin_pwr_fit =
            Some(lsq.two_parameter_line_fit(t, log_powers, sigmas, FitType::Linear));
        range.quad_pwr_fit =
            Some(lsq.two_parameter_line_fit(t, log_powers, sigmas, FitType::Quadratic));

        let log_corrected_sigmas: Vec<f64> = zip(sigmas.iter(), log_powers.iter())
            .map(|(s, l)| s / l.exp())
            .collect();

        range.lin_pwr_fit_err =
            Some(lsq.two_parameter_line_fit(t, log_powers, &log_corrected_sigmas, FitType::Linear));
        range.quad_pwr_fit_err = Some(lsq.two_parameter_line_fit(
            t,
            log_powers,
            &log_corrected_sigmas,
            FitType::Quadratic,
        ));
    }
    Ok(())
}

/// passing
fn acf_phase_fitting(ranges: &mut Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);
    for mut range in ranges {
        let phases = &range.phases.phases;
        let sigmas = &range.phases.std_dev;
        let t = &range.phases.t;

        let num_points = t.len();
        if phases.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::Message(format!(
                "Cannot perform acf phase fitting - dimension mismatch"
            )))?
        }
        range.phase_fit = Some(lsq.one_parameter_line_fit(t, phases, sigmas));
    }
    Ok(())
}

/// passing
fn xcf_phase_fitting(ranges: &mut Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);
    for mut range in ranges {
        let phases = &range.elev.phases;
        let sigmas = &range.elev.std_dev;
        let t = &range.elev.t;

        let num_points = t.len();
        if phases.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::Message(format!(
                "Cannot perform xcf phase fitting - dimension mismatch"
            )))?
        }
        range.elev_fit = Some(lsq.two_parameter_line_fit(t, phases, sigmas, FitType::Linear));
    }
    Ok(())
}

/// passing
fn calculate_phase_and_elev_sigmas(ranges: &mut Vec<RangeNode>, rec: &RawacfRecord) -> Result<()> {
    for mut range in ranges {
        let inverse_alpha_2: Vec<f64> = range.phase_alpha_2.iter().map(|x| 1.0 / x).collect();
        // let elevs_inverse_alpha_2: Vec<f64> = range.alpha_2.iter().map(|x| 1.0 / x).collect();
        let pwr_values: Vec<f64> = range
            .phases
            .t
            .iter()
            .map(|t| (-1.0 * range.lin_pwr_fit.as_ref().unwrap().slope.abs() * t).exp())
            .collect();
        let inverse_pwr_squared: Vec<f64> = pwr_values.iter().map(|x| 1.0 / (x * x)).collect();
        let phase_numerator: Vec<f64> = zip(inverse_alpha_2.iter(), inverse_pwr_squared.iter())
            .map(|(x, y)| x * y - 1.0)
            .collect();
        // let elev_numerator: Vec<f64> = zip(inverse_alpha_2.iter(), inverse_pwr_squared.iter())
        //     .map(|(x, y)| x * y - 1.0)
        //     .collect();
        let denominator = 2.0 * rec.num_averages as f64;
        let mut phase_sigmas: Vec<f64> = phase_numerator
            .iter()
            .map(|x| (x / denominator).sqrt())
            .collect();
        // let elev_sigmas: Vec<f64> = elev_numerator
        //     .iter()
        //     .map(|x| (x/denominator).sqrt())
        //     .collect();
        let _check: Vec<&f64> = phase_sigmas
            .iter()
            .filter(|&x| !x.is_finite())
            .collect();
        if _check.len() > 0 {
            Err(Fitacf3Error::Message(format!(
                "Phase sigmas bad at range {}",
                range.range_idx
            )))?
        }
        // let _check: Vec<&f64> = elev_sigmas
        //     .iter()
        //     .filter(|&x| (x.is_infinite() || x.is_nan()))
        //     .collect();
        // if _check.len() > 0 {
        //     Err(Fitacf3Error::Message(format!("Elevation sigmas bad at range {}", range.range_idx)))?
        // }

        // elev_sigmas[0] = elev_sigmas[1];
        range.phases.std_dev = phase_sigmas.clone();
        // Since lag 0 phase is included for elevation fit, set lag 0 sigma the same as lag 1 sigma
        phase_sigmas[0] = phase_sigmas[1];
        range.elev.std_dev = phase_sigmas; // = elev_sigmas;
    }
    Ok(())
}

/// passing
fn acf_phase_unwrap(ranges: &mut Vec<RangeNode>) {
    for mut range in ranges {
        let (mut slope_numerator, mut slope_denominator) = (0.0, 0.0);

        let phases = &range.phases.phases;
        let sigmas = &range.phases.std_dev;
        let t = &range.phases.t;

        // This is to skip the first element
        let mut phase_prev = phases[0];
        let mut sigma_prev = sigmas[0];
        let mut t_prev = t[0];

        let mut first_time = true;
        for (p, (s, t)) in zip(phases.iter(), zip(sigmas.iter(), t.iter())) {
            if !first_time {
                let phase_diff = p - phase_prev;
                let sigma_bar = (s + sigma_prev) / 2.0;
                let t_diff = t - t_prev;
                if phase_diff.abs() > PI {
                    slope_numerator += phase_diff / (sigma_bar * sigma_bar * t_diff);
                    slope_denominator += 1.0 / (sigma_bar * sigma_bar)
                }
                phase_prev = *p;
                sigma_prev = *s;
                t_prev = *t;
            } else {
                first_time = false;
            }
        }

        let piecewise_slope_estimate = slope_numerator / slope_denominator;
        let (new_phases, num_phase_jumps) = phase_correction(piecewise_slope_estimate, &phases, &t);
        if num_phase_jumps > 0 {
            let (mut sum_xx, mut sum_xy) = (0.0, 0.0);
            for (p, (s, t)) in zip(new_phases.iter(), zip(sigmas.iter(), t.iter())) {
                if *s > 0.0 {
                    sum_xy += p * t / (s * s);
                    sum_xx += (t * t) / (s * s);
                }
            }
            let corr_slope_estimate = sum_xy / sum_xx;
            let mut corr_slope_error = 0.0;
            for (p, (s, t)) in zip(new_phases.iter(), zip(sigmas.iter(), t.iter())) {
                if *s > 0.0 {
                    let temp = corr_slope_estimate * t - p;
                    corr_slope_error += temp * temp / (s * s);
                }
            }
            (sum_xx, sum_xy) = (0.0, 0.0);
            for (p, (s, t)) in zip(phases.iter(), zip(sigmas.iter(), t.iter())) {
                if *s > 0.0 {
                    sum_xy += p * t / (s * s);
                    sum_xx += t * t / (s * s);
                }
            }
            let orig_slope_estimate = sum_xy / sum_xx;
            let mut orig_slope_error = 0.0;
            for (p, (s, t)) in zip(new_phases.iter(), zip(sigmas.iter(), t.iter())) {
                if *s > 0.0 {
                    let temp = orig_slope_estimate * t - p;
                    orig_slope_error += temp * temp / (s * s);
                }
            }
            if orig_slope_error > corr_slope_error {
                range.phases.phases = new_phases;
            }
        }
    }
}

/// passing
fn xcf_phase_unwrap(ranges: &mut Vec<RangeNode>) -> Result<()> {
    for mut range in ranges {
        let (mut sum_xy, mut sum_xx) = (0.0, 0.0);

        let phases = &range.elev.phases;
        let sigmas = &range.elev.std_dev;
        let t = &range.elev.t;

        match range.phase_fit.as_ref() {
            None => Err(Fitacf3Error::Message(format!(
                "Phase fit must be defined to unwrap XCF phase"
            )))?,
            Some(fit) => {
                let mut new_phases = phase_correction(fit.slope, &phases, &t).0;
                for (p, (s, t)) in zip(new_phases.iter(), zip(sigmas.iter(), t.iter())) {
                    if *s > 0.0 {
                        sum_xy += p * t / (s * s);
                        sum_xx += t * t / (s * s);
                    }
                }
                let slope_estimate = sum_xy / sum_xx;
                new_phases = phase_correction(slope_estimate, &new_phases, &t).0;
                range.elev.phases = new_phases;
            }
        }
    }
    Ok(())
}

/// passing
fn phase_correction(slope_estimate: f64, phases: &Vec<f64>, times: &Vec<f64>) -> (Vec<f64>, i32) {
    let phase_predicted: Vec<f64> = times.iter().map(|t| t * slope_estimate).collect();

    // Round to 4 decimals places, so that 0.4999 rounds to 0.5, then up to 1.0
    let phase_diff: Vec<i32> = zip(phases.iter(), phase_predicted.iter())
        .map(|(p, pred)| ((((pred - p) / (2.0 * PI)) * 10000.0).round() / 10000.0).round() as i32)
        .collect();
    let corrected_phase: Vec<f64> = zip(phases.iter(), phase_diff.iter())
        .map(|(p, &corr)| p + corr as f64 * 2.0 * PI)
        .collect();
    let total_corrections: i32 = phase_diff.iter().map(|x| x.abs()).sum();
    (corrected_phase, total_corrections)
}

fn determinations(
    rec: &RawacfRecord,
    ranges: Vec<RangeNode>,
    noise_power: f32,
) -> Result<FitacfRecord> {
    let file_datetime = NaiveDateTime::parse_from_str(
        format!(
            "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
            rec.year, rec.month, rec.day, rec.hour, rec.minute, rec.second
        )
        .as_str(),
        "%Y%m%d %H:%M:%S",
    )
    .map_err(|_| Fitacf3Error::Message("Unable to interpret record timestamp".to_string()))?;
    let hdw = HdwInfo::new(rec.station_id, file_datetime)
        .map_err(|e| Fitacf3Error::Message(e.details))?;
    let range_list: Vec<i16> = ranges.iter().map(|r| r.range_num as i16).collect();
    let lag_0_power_db: Vec<f32> = rec
        .lag_zero_power
        .data
        .iter()
        .map(|p| {
            if p - noise_power > 0.0 {
                10.0 * ((p - noise_power) / noise_power).log10()
            } else {
                -50.0
            }
        })
        .collect();
    if range_list.len() == 0 {
        Ok(FitacfRecord {
            radar_revision_major: rec.radar_revision_major,
            radar_revision_minor: rec.radar_revision_minor,
            origin_code: rec.origin_code,
            origin_time: "".to_string(),    // TODO: Get current time
            origin_command: "".to_string(), // TODO: Get this
            control_program: rec.control_program,
            station_id: rec.station_id,
            year: rec.year,
            month: rec.month,
            day: rec.day,
            hour: rec.hour,
            minute: rec.minute,
            second: rec.second,
            microsecond: rec.microsecond,
            tx_power: rec.tx_power,
            num_averages: rec.num_averages,
            attenuation: rec.attenuation,
            lag_to_first_range: rec.lag_to_first_range,
            sample_separation: rec.sample_separation,
            error_code: rec.error_code,
            agc_status: rec.agc_status,
            low_power_status: rec.low_power_status,
            search_noise: rec.search_noise,
            mean_noise: rec.mean_noise,
            channel: rec.channel,
            beam_num: rec.beam_num,
            beam_azimuth: rec.beam_azimuth,
            scan_flag: rec.scan_flag,
            offset: rec.offset,
            rx_rise_time: rec.rx_rise_time,
            intt_second: rec.intt_second,
            intt_microsecond: rec.intt_microsecond,
            tx_pulse_length: rec.tx_pulse_length,
            multi_pulse_increment: rec.multi_pulse_increment,
            num_pulses: rec.num_pulses,
            num_lags: rec.num_lags,
            num_lags_extras: rec.num_lags_extras,
            if_mode: rec.if_mode,
            num_ranges: rec.num_ranges,
            first_range: rec.first_range,
            range_sep: rec.range_sep,
            xcf_flag: rec.xcf_flag,
            tx_freq: rec.tx_freq,
            max_power: rec.max_power,
            max_noise_level: rec.max_noise_level,
            comment: rec.comment.clone(),
            algorithm: None,
            fitacf_revision_major: FITACF_REVISION_MAJOR,
            fitacf_revision_minor: FITACF_REVISION_MINOR,
            sky_noise: noise_power,
            lag_zero_noise: 0.0,
            velocity_noise: 0.0,
            tdiff: None,
            pulse_table: rec.pulse_table.clone(),
            lag_table: rec.lag_table.clone(),
            lag_zero_power: convert_to_dmapvec(lag_0_power_db),
            range_list: convert_to_dmapvec(vec![]),
            fitted_points: convert_to_dmapvec(vec![]),
            quality_flag: convert_to_dmapvec(vec![]),
            ground_flag: convert_to_dmapvec(vec![]),
            lambda_power: convert_to_dmapvec(vec![]),
            lambda_power_error: convert_to_dmapvec(vec![]),
            sigma_power: convert_to_dmapvec(vec![]),
            sigma_power_error: convert_to_dmapvec(vec![]),
            velocity: convert_to_dmapvec(vec![]),
            velocity_error: convert_to_dmapvec(vec![]),
            lambda_spectral_width: convert_to_dmapvec(vec![]),
            lambda_spectral_width_error: convert_to_dmapvec(vec![]),
            sigma_spectral_width: convert_to_dmapvec(vec![]),
            sigma_spectral_width_error: convert_to_dmapvec(vec![]),
            lambda_std_dev: convert_to_dmapvec(vec![]),
            sigma_std_dev: convert_to_dmapvec(vec![]),
            phi_std_dev: convert_to_dmapvec(vec![]),
            xcf_quality_flag: None,
            xcf_ground_flag: None,
            lambda_xcf_power: None,
            lambda_xcf_power_error: None,
            sigma_xcf_power: None,
            sigma_xcf_power_error: None,
            xcf_velocity: None,
            xcf_velocity_error: None,
            lambda_xcf_spectral_width: None,
            lambda_xcf_spectral_width_error: None,
            sigma_xcf_spectral_width: None,
            sigma_xcf_spectral_width_error: None,
            lag_zero_phi: None,
            lag_zero_phi_error: None,
            elevation: None,
            elevation_fitted: None,
            elevation_error: None,
            elevation_low: None,
            elevation_high: None,
            lambda_xcf_std_dev: None,
            sigma_xcf_std_dev: None,
            phi_xcf_std_dev: None,
        })
    } else {
        let num_lags: Vec<i16> = ranges
            .iter()
            .map(|r| r.powers.ln_power.len() as i16)
            .collect();
        let quality_flag: Vec<i8> = range_list.iter().map(|_| 1).collect();
        let noise_db: f32 = 10.0 * noise_power.log10();
        let power_linear: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * r
                    .lin_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf without linear fitted power")
                    .intercept as f32
                    / (10.0 as f32).ln()
                    - noise_db
            })
            .collect();
        let power_linear_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * (r
                    .lin_pwr_fit_err
                    .as_ref()
                    .expect("Unable to make fitacf without linear fitted power error")
                    .variance_intercept as f32)
                    .sqrt()
                    / (10.0 as f32).ln()
            })
            .collect();
        let power_quadratic: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * (r
                    .quad_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf without quadratic fitted power")
                    .intercept as f32)
                    / (10.0 as f32).ln()
                    - noise_db
            })
            .collect();
        let power_quadratic_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * (r
                    .quad_pwr_fit_err
                    .as_ref()
                    .expect("Unable to make fitacf without quadratic fitted power error")
                    .variance_intercept as f32)
                    .sqrt()
                    / (10.0 as f32).ln()
            })
            .collect();
        let velocity_conversion: f32 =
            299792458.0 * hdw.velocity_sign / (4.0 * PI_f32 * rec.tx_freq as f32 * 1000.0);
        let velocity: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.phase_fit
                    .as_ref()
                    .expect("Unable to make fitacf without fitted velocity")
                    .slope as f32)
                    * velocity_conversion
            })
            .collect();
        let velocity_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.phase_fit
                    .as_ref()
                    .expect("Unable to make fitacf without fitted velocity")
                    .variance_slope as f32)
                    .sqrt()
                    * velocity_conversion
            })
            .collect();
        let width_conversion: f32 =
            299792458.0 * 2.0 / (4.0 * PI_f32 * rec.tx_freq as f32 * 1000.0);
        let spectral_width_linear: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.lin_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf spectral width without fitted power")
                    .slope as f32)
                    .abs()
                    * width_conversion
            })
            .collect();
        let spectral_width_linear_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.lin_pwr_fit_err
                    .as_ref()
                    .expect("Unable to make fitacf spectral width error without fitted power error")
                    .variance_slope as f32)
                    .sqrt()
                    * width_conversion
            })
            .collect();
        let quadratic_width_conversion: f32 =
            299792458.0 * (2.0 as f32).ln().sqrt() / (PI_f32 * rec.tx_freq as f32 * 1000.0);
        let spectral_width_quadratic: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.quad_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf quadratic spectral width without fitted power")
                    .slope as f32)
                    .abs()
                    .sqrt()
                    * quadratic_width_conversion
            })
            .collect();
        let spectral_width_quadratic_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.quad_pwr_fit_err.as_ref().expect("Unable to make fitacf quadratic spectral width error without fitted power error")
                    .variance_slope as f32).sqrt() * quadratic_width_conversion /
                    ((r.quad_pwr_fit.as_ref().expect("Unable to make fitacf quadratic spectral width error without fitted power error")
                        .slope as f32).abs().sqrt() * 2.0)
            })
            .collect();
        let std_dev_linear: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.lin_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf linear std deviation without fitted power")
                    .chi_squared as f32
            })
            .collect();
        let std_dev_quadratic: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.quad_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf quadratic std deviation without fitted power")
                    .chi_squared as f32
            })
            .collect();
        let std_dev_phi: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.phase_fit
                    .as_ref()
                    .expect("Unable to make fitacf phi std deviation")
                    .chi_squared as f32
            })
            .collect();
        let groundscatter_flag: Vec<i8> = zip(velocity.iter(), spectral_width_linear.iter())
            .map(|(v, w)| {
                if v.abs() - (V_MAX - w * (V_MAX / W_MAX)) < 1.0 {
                    1
                } else {
                    0
                }
            })
            .collect();
        let xcfs = &rec
            .xcfs
            .as_ref()
            .expect("Unable to make fitacf xcf_phi0")
            .data;
        let xcf_phi0: Vec<f32> = ranges
            .iter()
            .map(|r| {
                xcfs[r.range_idx * rec.num_lags as usize * 2 + 1]
                    .atan2(xcfs[r.range_idx * rec.num_lags as usize * 2])
                    * hdw.phase_sign
            })
            .collect();
        let xcf_phi0_err: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.elev_fit
                    .as_ref()
                    .expect("Unable to make fitacf xcf_phi0_err")
                    .variance_intercept as f32)
                    .sqrt()
            })
            .collect();
        let xcf_phi_std_dev: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.elev_fit
                    .as_ref()
                    .expect("Unable to make fitacf xcf_phi_std_dev")
                    .chi_squared as f32
            })
            .collect();
        let (elevation_low, elevation_normal, elevation_high) =
            calculate_elevation(&ranges, rec, &xcf_phi0, &hdw);

        let float_zeros = DmapVec {
            data: quality_flag.iter().map(|_| 0.0 as f32).collect(),
            dimensions: vec![quality_flag.len() as i32],
        };
        let i8_zeros = DmapVec {
            data: quality_flag.iter().map(|_| 0 as i8).collect(),
            dimensions: vec![quality_flag.len() as i32],
        };

        Ok(FitacfRecord {
            radar_revision_major: rec.radar_revision_major,
            radar_revision_minor: rec.radar_revision_minor,
            origin_code: rec.origin_code,
            origin_time: "".to_string(),    // TODO: Get current time
            origin_command: "".to_string(), // TODO: Get this
            control_program: rec.control_program,
            station_id: rec.station_id,
            year: rec.year,
            month: rec.month,
            day: rec.day,
            hour: rec.hour,
            minute: rec.minute,
            second: rec.second,
            microsecond: rec.microsecond,
            tx_power: rec.tx_power,
            num_averages: rec.num_averages,
            attenuation: rec.attenuation,
            lag_to_first_range: rec.lag_to_first_range,
            sample_separation: rec.sample_separation,
            error_code: rec.error_code,
            agc_status: rec.agc_status,
            low_power_status: rec.low_power_status,
            search_noise: rec.search_noise,
            mean_noise: rec.mean_noise,
            channel: rec.channel,
            beam_num: rec.beam_num,
            beam_azimuth: rec.beam_azimuth,
            scan_flag: rec.scan_flag,
            offset: rec.offset,
            rx_rise_time: rec.rx_rise_time,
            intt_second: rec.intt_second,
            intt_microsecond: rec.intt_microsecond,
            tx_pulse_length: rec.tx_pulse_length,
            multi_pulse_increment: rec.multi_pulse_increment,
            num_pulses: rec.num_pulses,
            num_lags: rec.num_lags,
            num_lags_extras: rec.num_lags_extras,
            if_mode: rec.if_mode,
            num_ranges: rec.num_ranges,
            first_range: rec.first_range,
            range_sep: rec.range_sep,
            xcf_flag: rec.xcf_flag,
            tx_freq: rec.tx_freq,
            max_power: rec.max_power,
            max_noise_level: rec.max_noise_level,
            comment: rec.comment.clone(),
            algorithm: None,
            fitacf_revision_major: FITACF_REVISION_MAJOR,
            fitacf_revision_minor: FITACF_REVISION_MINOR,
            sky_noise: noise_power,
            lag_zero_noise: 0.0,
            velocity_noise: 0.0,
            tdiff: None,
            pulse_table: rec.pulse_table.clone(),
            lag_table: rec.lag_table.clone(),
            lag_zero_power: convert_to_dmapvec(lag_0_power_db),
            range_list: convert_to_dmapvec(range_list),
            fitted_points: convert_to_dmapvec(num_lags),
            quality_flag: convert_to_dmapvec(quality_flag),
            ground_flag: convert_to_dmapvec(groundscatter_flag),
            lambda_power: convert_to_dmapvec(power_linear),
            lambda_power_error: convert_to_dmapvec(power_linear_error),
            sigma_power: convert_to_dmapvec(power_quadratic),
            sigma_power_error: convert_to_dmapvec(power_quadratic_error),
            velocity: convert_to_dmapvec(velocity),
            velocity_error: convert_to_dmapvec(velocity_error),
            lambda_spectral_width: convert_to_dmapvec(spectral_width_linear),
            lambda_spectral_width_error: convert_to_dmapvec(spectral_width_linear_error),
            sigma_spectral_width: convert_to_dmapvec(spectral_width_quadratic),
            sigma_spectral_width_error: convert_to_dmapvec(spectral_width_quadratic_error),
            lambda_std_dev: convert_to_dmapvec(std_dev_linear),
            sigma_std_dev: convert_to_dmapvec(std_dev_quadratic),
            phi_std_dev: convert_to_dmapvec(std_dev_phi),
            xcf_quality_flag: Some(i8_zeros.clone()),
            xcf_ground_flag: Some(i8_zeros.clone()),
            lambda_xcf_power: Some(float_zeros.clone()),
            lambda_xcf_power_error: Some(float_zeros.clone()),
            sigma_xcf_power: Some(float_zeros.clone()),
            sigma_xcf_power_error: Some(float_zeros.clone()),
            xcf_velocity: Some(float_zeros.clone()),
            xcf_velocity_error: Some(float_zeros.clone()),
            lambda_xcf_spectral_width: Some(float_zeros.clone()),
            lambda_xcf_spectral_width_error: Some(float_zeros.clone()),
            sigma_xcf_spectral_width: Some(float_zeros.clone()),
            sigma_xcf_spectral_width_error: Some(float_zeros.clone()),
            lag_zero_phi: Some(convert_to_dmapvec(xcf_phi0)),
            lag_zero_phi_error: Some(convert_to_dmapvec(xcf_phi0_err)),
            elevation: Some(convert_to_dmapvec(elevation_normal)),
            elevation_fitted: None,
            elevation_error: None,
            elevation_low: Some(convert_to_dmapvec(elevation_low)),
            elevation_high: Some(convert_to_dmapvec(elevation_high)),
            lambda_xcf_std_dev: Some(float_zeros.clone()),
            sigma_xcf_std_dev: Some(float_zeros.clone()),
            phi_xcf_std_dev: Some(convert_to_dmapvec(xcf_phi_std_dev)),
        })
    }
}

fn convert_to_dmapvec<T: InDmap>(vals: Vec<T>) -> DmapVec<T> {
    DmapVec {
        dimensions: vec![vals.len() as i32],
        data: vals,
    }
}

fn calculate_elevation(
    ranges: &Vec<RangeNode>,
    rec: &RawacfRecord,
    xcf_phi0: &Vec<f32>,
    hdw: &HdwInfo,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let x = hdw.intf_offset_x;
    let y = hdw.intf_offset_y;
    let z = hdw.intf_offset_z;

    let array_separation: f32 = (x * x + y * y + z * z).sqrt();
    let mut elevation_corr = (z / array_separation).asin();
    let phi_sign: f32;
    if y > 0.0 {
        phi_sign = 1.0;
    } else {
        phi_sign = -1.0;
        elevation_corr *= -1.0;
    }
    let azimuth_offset = hdw.max_num_beams as f32 / 2.0 - 0.5;
    let phi_0 =
        (hdw.beam_separation * (rec.beam_num as f32 - azimuth_offset) * PI_f32 / 180.0).cos();
    let wave_num = 2.0 * PI_f32 * rec.tx_freq as f32 * 1000.0 / 299792458.0;
    let cable_offset = -2.0 * PI_f32 * rec.tx_freq as f32 * 1000.0 * hdw.tdiff_a * 1.0e-6;
    let phase_diff_max = phi_sign * wave_num * array_separation * phi_0 + cable_offset;
    let mut psi: Vec<f32> = ranges
        .iter()
        .map(|r| {
            let x = r
                .elev_fit
                .as_ref()
                .expect("Unable to find elevation without fitted elevation")
                .intercept as f32;
            let mut y =
                x + 2.0 * PI_f32 * ((phase_diff_max - x) / (2.0 * PI_f32)).floor() - cable_offset;
            if phi_sign < 0.0 {
                y += 2.0 * PI_f32;
            }
            y
        })
        .collect();
    let mut psi_kd: Vec<f32> = psi
        .iter()
        .map(|p| p / (wave_num * array_separation))
        .collect();
    let mut theta: Vec<f32> = psi_kd.iter().map(|p| phi_0 * phi_0 - p * p).collect();
    let elevation: Vec<f32> = theta
        .iter()
        .map(|&t| {
            if t < 0.0 || t.abs() > 1.0 {
                -elevation_corr
            } else {
                t.sqrt().asin()
            }
        })
        .collect();
    let elevation_high: Vec<f32> = elevation
        .iter()
        .map(|e| (e + elevation_corr) * 180.0 / PI_f32)
        .collect();
    let psi_k2d2: Vec<f32> = psi
        .iter()
        .map(|p| p / (wave_num * wave_num * array_separation * array_separation))
        .collect();
    let df_by_dy: Vec<f32> = zip(psi_k2d2.iter(), theta.iter())
        .map(|(p, t)| p / (t * (1.0 - t)).sqrt())
        .collect();
    let errors: Vec<f32> = ranges
        .iter()
        .map(|r| {
            r.elev_fit
                .as_ref()
                .expect("Unable to calculate elevation errors")
                .variance_intercept as f32
        })
        .collect();
    let elevations_low: Vec<f32> = zip(errors.iter(), df_by_dy.iter())
        .map(|(e, d)| e.sqrt() * d.abs() * 180.0 / PI_f32)
        .collect();

    // This time, use the xcf lag0 phase
    psi = xcf_phi0
        .iter()
        .map(|&x| {
            let mut y = x as f32
                + 2.0 * PI_f32 * ((phase_diff_max - x as f32) / (2.0 * PI_f32)).floor()
                - cable_offset;
            if phi_sign < 0.0 {
                y += 2.0 * PI_f32;
            }
            y
        })
        .collect();
    psi_kd = psi
        .iter()
        .map(|p| p / (wave_num * array_separation))
        .collect();
    theta = psi_kd.iter().map(|p| phi_0 * phi_0 - p * p).collect();
    let elevation_normal: Vec<f32> = theta
        .iter()
        .map(|&t| {
            if t < 0.0 || t.abs() > 1.0 {
                -180.0 / PI_f32 * elevation_corr
            } else {
                (t + elevation_corr).sqrt().asin() * 180.0 / PI_f32
            }
        })
        .collect();
    (elevations_low, elevation_normal, elevation_high)
}
