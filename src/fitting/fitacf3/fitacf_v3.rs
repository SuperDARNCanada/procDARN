use crate::fitting::fitacf3::fitstruct::{Alpha, FitType, LagNode, RangeNode};
use dmap::formats::{RawacfRecord, FitacfRecord};
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fmt::Display;
use std::f64::consts::PI;
use std::iter::zip;
use crate::fitting::fitacf3::filtering::filter_bad_fits;
use crate::fitting::fitacf3::least_squares::LeastSquares;

type Result<T> = std::result::Result<T, Fitacf3Error>;

#[derive(Debug, Clone)]
pub enum Fitacf3Error {
    // Parse(String, Vec<u8>),
    // BadVal(String, DmapType),
    Message(String),
    Lookup(String),
    Mismatch { msg: String },
    // CastError(String, PodCastError),
}

impl Error for Fitacf3Error {}

impl Display for Fitacf3Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Fitacf3Error::Message(msg) => write!(f, "{}", msg),
            Fitacf3Error::Lookup(msg) => write!(f, "{}", msg),
            Fitacf3Error::Mismatch{msg} => write!(f, "{}", msg)

            // DmapError::BadVal(msg, val) => write!(f, "{}: {:?}", msg, val),
            // DmapError::Parse(msg, val) => write!(f, "{}: {:?}", msg, val),
            // DmapError::CastError(msg, err) => write!(f, "{}: {}", msg, err.to_string()),
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
        let number = lag_table.data[2*i + 1] - lag_table.data[2*i];   // flattened, we want row i, cols 1 and 0
        for j in 0..record.num_pulses as usize {
            if lag_table.data[2*i] == pulse_table.data[j] {
                pulse_1_idx = j;
            }
            if lag_table.data[2*i + 1] == pulse_table.data[j] {
                pulse_2_idx = j;
            }
        }
        let sample_base_1= (lag_table.data[2*i] * (multi_pulse_increment / sample_separation)) as i32;
        let sample_base_2 = (lag_table.data[2*i + 1] * (multi_pulse_increment / sample_separation)) as i32;
        let pulses = lag_table.data[i];
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

fn fit_rawacf_record(record: &RawacfRecord) -> Result<FitacfRecord> {
    let lags = create_lag_list(record);

    let mut noise_power = 0.0;
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
    acf_power_fitting(&range_list)?;
    calculate_phase_and_elev_sigmas(&range_list, record)?;
    acf_phase_unwrap(&range_list);
    acf_phase_fitting(&range_list);
    filter_bad_fits(&range_list);
    xcf_phase_unwrap(&range_list);
    xcf_phase_fitting(&range_list);

    let determined_parameters = Determinations::new(record, range_list, noise_power, tdiff);
    Err(Fitacf3Error::Message(format!("Unable to fit record")))
}

// fn calculate_alpha_at_lags(lag: LagNode, range: RangeNode, lag_zero_power: f64) -> Alpha {
//     let pulse_i_cri = range.cross_range_interference[lag.pulses[0] as usize];
//     let pulse_j_cri = range.cross_range_interference[lag.pulses[1] as usize];
//
//     let lag_idx = lag.lag_idx;
//     let alpha_2 = lag_zero_power*lag_zero_power / ((lag_zero_power + pulse_i_cri) * (lag_zero_power + pulse_j_cri));
//     // let range.alpha_2 = alpha_2;
//
// }

fn acf_power_fitting(ranges: &Vec<RangeNode>) -> Result<()> {
    let mut lsq = LeastSquares::new(1, 1);

    for mut range in *ranges {
        let log_powers = &range.powers.ln_power;
        let sigmas = &range.powers.std_dev;
        let t = &range.powers.t;
        let num_points = range.powers.ln_power.len();
        if t.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::Message(format!("Cannot perform acf power fitting - dimension mismatch")))?
        }
        range.lin_pwr_fit = Some(lsq.two_parameter_line_fit(t, log_powers, sigmas, FitType::Linear));
        range.quad_pwr_fit = Some(lsq.two_parameter_line_fit(t, log_powers, sigmas, FitType::Quadratic));

        let log_corrected_sigmas: Vec<f64> = zip(sigmas.iter(), log_powers.iter())
            .map(|(s, l)| s / l.exp())
            .collect();

        range.lin_pwr_fit_err = Some(lsq.two_parameter_line_fit(t, log_powers, &log_corrected_sigmas, FitType::Linear));
        range.quad_pwr_fit_err = Some(lsq.two_parameter_line_fit(t, log_powers, &log_corrected_sigmas, FitType::Quadratic));
    }
    Ok(())
}

fn acf_phase_fitting(ranges: &Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);
    for mut range in *ranges {
        let phases = &range.phases.phases;
        let sigmas = &range.phases.std_dev;
        let t = &range.phases.t;

        let num_points = t.len();
        if phases.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::Message(format!("Cannot perform acf phase fitting - dimension mismatch")))?
        }
        range.phase_fit = Some(lsq.one_parameter_line_fit(t, phases, sigmas));
    }
    Ok(())
}

fn xcf_phase_fitting(ranges: &Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);
    for mut range in *ranges {
        let phases = &range.elev.phases;
        let sigmas = &range.elev.std_dev;
        let t = &range.elev.t;

        let num_points = t.len();
        if phases.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::Message(format!("Cannot perform xcf phase fitting - dimension mismatch")))?
        }
        range.elev_fit = Some(lsq.two_parameter_line_fit(t, phases, sigmas, FitType::Linear));
    }
    Ok(())
}

fn calculate_phase_and_elev_sigmas(ranges: &Vec<RangeNode>, rec: &RawacfRecord) -> Result<()>{
    for mut range in *ranges {
        let phase_inverse_alpha_2: Vec<f64> = range.phases.alpha_2.iter().map(|x| 1.0 / x).collect();
        let elevs_inverse_alpha_2: Vec<f64> = range.elev.alpha_2.iter().map(|x| 1.0 / x).collect();

        let pwr_values: Vec<f64> = range.phases.t
            .iter()
            .map(|t| (-1.0 * range.lin_pwr_fit.unwrap().slope.abs() * t).exp())
            .collect();
        let inverse_pwr_squared: Vec<f64> = pwr_values
            .iter()
            .map(|x| 1.0 / (x*x))
            .collect();
        let phase_numerator: Vec<f64> = zip(phase_inverse_alpha_2.iter(), inverse_pwr_squared.iter())
            .map(|(x, y)| x * y - 1.0)
            .collect();
        let elev_numerator: Vec<f64> = zip(elevs_inverse_alpha_2.iter(), inverse_pwr_squared.iter())
            .map(|(x, y)| x * y - 1.0)
            .collect();
        let denominator = 2.0 * rec.num_averages as f64;
        let phase_sigmas: Vec<f64> = phase_numerator
            .iter()
            .map(|x| (x/denominator).sqrt())
            .collect();
        let elev_sigmas: Vec<f64> = elev_numerator
            .iter()
            .map(|x| (x/denominator).sqrt())
            .collect();
        let _check: Vec<&f64> = phase_sigmas
            .iter()
            .filter(|&x| (x.is_infinite() || x.is_nan()))
            .collect();
        if _check.len() > 0 {
            Err(Fitacf3Error::Message(format!("Phase sigmas bad at range {}", range.range_idx)))?
        }
        let _check: Vec<&f64> = elev_sigmas
            .iter()
            .filter(|&x| (x.is_infinite() || x.is_nan()))
            .collect();
        if _check.len() > 0 {
            Err(Fitacf3Error::Message(format!("Elevation sigmas bad at range {}", range.range_idx)))?
        }
        // Since lag 0 phase is included for elevation fit, set lag 0 sigma the same as lag 1 sigma
        elev_sigmas[0] = elev_sigmas[1];
        range.phases.std_dev = phase_sigmas;
        range.elev.std_dev = elev_sigmas;
    }
    Ok(())

}

fn acf_phase_unwrap(ranges: &Vec<RangeNode>) {
    for mut range in *ranges {
        let (mut slope_numerator, mut slope_denominator) = (0.0, 0.0);

        let phases = range.phases.phases;
        let sigmas = range.phases.std_dev;
        let t = range.phases.t;

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

fn xcf_phase_unwrap(ranges: &Vec<RangeNode>) -> Result<()> {
    for mut range in *ranges {
        let (mut sum_xy, mut sum_xx) = (0.0, 0.0);

        let phases = range.elev.phases;
        let sigmas = range.elev.std_dev;
        let t = range.elev.t;

        match range.phase_fit {
            None => Err(Fitacf3Error::Message(format!("Phase fit must be defined to unwrap XCF phase")))?,
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

fn phase_correction(slope_estimate: f64, phases: &Vec<f64>, times: &Vec<f64>) -> (Vec<f64>, i32) {
    let phase_predicted: Vec<f64> = times.iter().map(|t| t * slope_estimate).collect();

    // Round to 4 decimals places, so that 0.4999 rounds to 0.5, then up to 1.0
    let mut phase_diff: Vec<i32> = zip(phases.iter(), phase_predicted.iter())
        .map(|(p, pred)| ((((pred - p) / (2.0 * PI)) * 10000.0).round() / 10000.0).round() as i32)
        .collect();
    let corrected_phase: Vec<f64> = zip(phases.iter(), phase_diff.iter())
        .map(|(p, &corr)| p + corr as f64 * 2.0 * PI)
        .collect();
    let total_corrections: i32 = phase_diff.iter().map(|x| x.abs()).sum();
    (corrected_phase, total_corrections)
}