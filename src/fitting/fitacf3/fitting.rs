use crate::fitting::fitacf3::fitacf_v3::Fitacf3Error;
use crate::fitting::fitacf3::fitstruct::{PowerFitType, RangeNode};
use crate::fitting::fitacf3::least_squares::LeastSquares;
use crate::utils::rawacf::Rawacf;
use std::f64::consts::PI;
use std::iter::zip;

type Result<T> = std::result::Result<T, Fitacf3Error>;

/// Fits the power of ACF data.
pub(crate) fn acf_power_fitting(ranges: &mut Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);

    for range in ranges {
        let log_powers = &range.powers.ln_power;
        let sigmas = &range.powers.std_dev;
        let t = &range.powers.t;
        let num_points = range.powers.ln_power.len();
        if t.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::BadFit(
                "Cannot perform acf power fitting - dimension mismatch".to_string(),
            ))?;
        }
        range.lin_pwr_fit =
            Some(lsq.two_parameter_line_fit(t, log_powers, sigmas, &PowerFitType::Linear));
        range.quad_pwr_fit =
            Some(lsq.two_parameter_line_fit(t, log_powers, sigmas, &PowerFitType::Quadratic));

        let log_corrected_sigmas: Vec<f64> = zip(sigmas.iter(), log_powers.iter())
            .map(|(s, l)| s / l.exp())
            .collect();

        range.lin_pwr_fit_err = Some(lsq.two_parameter_line_fit(
            t,
            log_powers,
            &log_corrected_sigmas,
            &PowerFitType::Linear,
        ));
        range.quad_pwr_fit_err = Some(lsq.two_parameter_line_fit(
            t,
            log_powers,
            &log_corrected_sigmas,
            &PowerFitType::Quadratic,
        ));
    }
    Ok(())
}

/// Fits the phase of ACF data.
pub(crate) fn acf_phase_fitting(ranges: &mut Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);
    for range in ranges {
        let phases = &range.phases.phases;
        let sigmas = &range.phases.std_dev;
        let t = &range.phases.t;

        let num_points = t.len();
        if phases.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::BadFit(
                "Cannot perform acf phase fitting - dimension mismatch".to_string(),
            ))?;
        }
        range.phase_fit = Some(lsq.one_parameter_line_fit(t, phases, sigmas));
    }
    Ok(())
}

/// Fits the phase of XCF data.
pub(crate) fn xcf_phase_fitting(ranges: &mut Vec<RangeNode>) -> Result<()> {
    let lsq = LeastSquares::new(1, 1);
    for range in ranges {
        let phases = &range.elev.phases;
        let sigmas = &range.elev.std_dev;
        let t = &range.elev.t;

        let num_points = t.len();
        if phases.len() != num_points || sigmas.len() != num_points {
            Err(Fitacf3Error::BadFit(
                "Cannot perform xcf phase fitting - dimension mismatch".to_string(),
            ))?;
        }
        range.elev_fit = Some(lsq.two_parameter_line_fit(t, phases, sigmas, &PowerFitType::Linear));
    }
    Ok(())
}

/// Calculates standard deviations for phase and elevation fits.
pub(crate) fn calculate_phase_and_elev_sigmas(
    ranges: &mut Vec<RangeNode>,
    rec: &Rawacf,
) -> Result<()> {
    for range in ranges {
        let inverse_alpha_2: Vec<f64> = range.phase_alpha_2.iter().map(|x| 1.0 / x).collect();
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
        let denominator = 2.0 * rec.nave as f64;
        let mut phase_sigmas: Vec<f64> = phase_numerator
            .iter()
            .map(|x| (x / denominator).sqrt())
            .collect();
        if phase_sigmas.iter().filter(|&x| !x.is_finite()).count() > 0 {
            Err(Fitacf3Error::BadFit(format!(
                "Phase sigmas infinite at range {}",
                range.range_idx
            )))?;
        }
        range.phases.std_dev = phase_sigmas.clone();
        // Since lag 0 phase is included for elevation fit, set lag 0 sigma the same as lag 1 sigma
        phase_sigmas[0] = phase_sigmas[1];
        range.elev.std_dev = phase_sigmas; // = elev_sigmas;
    }
    Ok(())
}

/// Applies 2π phase unwrapping to ACF phases.
pub(crate) fn acf_phase_unwrap(ranges: &mut Vec<RangeNode>) {
    for range in ranges {
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
            if first_time {
                first_time = false;
            } else {
                let phase_diff = p - phase_prev;
                let sigma_bar = (s + sigma_prev) / 2.0;
                let t_diff = t - t_prev;
                if phase_diff.abs() < PI {
                    slope_numerator += phase_diff / (sigma_bar * sigma_bar * t_diff);
                    slope_denominator += 1.0 / (sigma_bar * sigma_bar);
                }
                phase_prev = *p;
                sigma_prev = *s;
                t_prev = *t;
            }
        }

        let piecewise_slope_estimate = slope_numerator / slope_denominator;
        let (new_phases, num_phase_jumps) = phase_correction(piecewise_slope_estimate, phases, t);
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
            for (p, (s, t)) in zip(phases.iter(), zip(sigmas.iter(), t.iter())) {
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

/// Applies 2π phase unwrapping to XCF phases
pub(crate) fn xcf_phase_unwrap(ranges: &mut Vec<RangeNode>) -> Result<()> {
    for range in ranges {
        let (mut sum_xy, mut sum_xx) = (0.0, 0.0);

        let phases = &range.elev.phases;
        let sigmas = &range.elev.std_dev;
        let t = &range.elev.t;

        match range.phase_fit.as_ref() {
            None => Err(Fitacf3Error::BadFit(
                "Phase fit must be defined to unwrap XCF phase".to_string(),
            ))?,
            Some(fit) => {
                let mut new_phases = phase_correction(fit.slope, phases, t).0;
                for (p, (s, t)) in zip(new_phases.iter(), zip(sigmas.iter(), t.iter())) {
                    if *s > 0.0 {
                        sum_xy += p * t / (s * s);
                        sum_xx += t * t / (s * s);
                    }
                }
                let slope_estimate = sum_xy / sum_xx;
                new_phases = phase_correction(slope_estimate, &new_phases, t).0;
                range.elev.phases = new_phases;
            }
        }
    }
    Ok(())
}

/// Determines which points need a 2π phase correction applied
fn phase_correction(slope_estimate: f64, phases: &[f64], times: &[f64]) -> (Vec<f64>, i32) {
    let phase_predicted: Vec<f64> = times.iter().map(|t| t * slope_estimate).collect();

    // Round to 5 decimals places, so that 0.49999 rounds to 0.5, then up to 1.0
    let phase_diff: Vec<i32> = zip(phases.iter(), phase_predicted.iter())
        .map(|(p, pred)| {
            ((((pred - p) / (2.0 * PI)) * 100_000.0).round() / 100_000.0).round() as i32
        })
        .collect();
    let corrected_phase: Vec<f64> = zip(phases.iter(), phase_diff.iter())
        .map(|(p, &corr)| p + corr as f64 * 2.0 * PI)
        .collect();
    let total_corrections: i32 = phase_diff
        .iter()
        .map(|x| x.abs())
        .max()
        .map_or_else(|| 0, |x| x);
    (corrected_phase, total_corrections)
}
