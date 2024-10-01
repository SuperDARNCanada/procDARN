use crate::fitting::common::error::FittingError;
use crate::fitting::common::fitstruct::RangeNode;
use crate::fitting::fitacf3::fitacf_v3::{ALPHA_CUTOFF, FLUCTUATION_CUTOFF_COEFFICIENT, MIN_LAGS};
use crate::utils::rawacf::Rawacf;
use is_close::is_close;

/// Removes all lags that have infinite power values.
pub(crate) fn filter_infinite_lags(ranges: &mut Vec<RangeNode>) {
    for range in ranges {
        let mut infinite_indices = vec![];
        for i in 0..range.powers.ln_power.len() {
            if !range.powers.ln_power[i].is_finite() {
                infinite_indices.push(i);
            }
        }
        for i in infinite_indices.iter().rev() {
            range.powers.remove(*i);
            range.power_alpha_2.remove(*i);
        }
    }
}

/// Removes all lags after a lag with low power
pub(crate) fn filter_low_power_lags(rec: &Rawacf, ranges: &mut Vec<RangeNode>) {
    if rec.nave <= 0 {
        return;
    }
    for range in ranges {
        let range_num = range.range_num;
        if range.powers.ln_power.is_empty() {
            continue;
        }
        let log_sigma_fluc = (FLUCTUATION_CUTOFF_COEFFICIENT * rec.pwr0[range_num as usize]
            / f32::from(2 * rec.nave).sqrt())
        .ln();
        let mut bad_indices = vec![];
        let mut cutoff_lag = rec.mplgs as usize + 1;

        for idx in 0..range.powers.ln_power.len() {
            if idx > cutoff_lag {
                bad_indices.push(idx);
            } else {
                let log_power = range.powers.ln_power[idx];
                let alpha_2 = range.power_alpha_2[idx];
                if ((1_f64 / alpha_2.sqrt()) <= ALPHA_CUTOFF as f64)
                    && ((log_power < log_sigma_fluc as f64)
                        || is_close!(log_power, log_sigma_fluc as f64))
                {
                    cutoff_lag = idx;
                    bad_indices.push(idx);
                }
            }
        }
        for i in bad_indices.iter().rev() {
            range.powers.remove(*i);
            // range.phases.remove(*i);
            // range.elev.remove(*i);
            range.power_alpha_2.remove(*i);
        }
    }
}

/// Removes range gates that either contain too weak of a fit, or used too few lags when fitting.
pub(crate) fn filter_bad_acfs(rec: &Rawacf, ranges: &mut Vec<RangeNode>, noise_power: f32) {
    if rec.nave <= 0 {
        return;
    }
    let cutoff_power = noise_power * 2.0;
    let mut bad_indices = vec![];
    for (idx, range) in ranges.iter().enumerate() {
        let power = rec.pwr0[range.range_num as usize];
        let num_powers = range.powers.ln_power.len();
        if (power <= cutoff_power) || (num_powers < MIN_LAGS as usize) {
            bad_indices.push(idx);
        } else {
            let power_value = range.powers.ln_power[0];
            let mut all_equal = true;
            for pwr in &range.powers.ln_power {
                if !is_close!(*pwr, power_value) {
                    all_equal = false;
                }
            }
            if all_equal {
                bad_indices.push(idx);
            }
        }
    }
    for idx in bad_indices.iter().rev() {
        ranges.remove(*idx);
    }
}

/// Removes all ranges that have not had phase, lambda power, or quadratic power fitted
pub(crate) fn filter_bad_fits(ranges: &mut Vec<RangeNode>) -> Result<(), FittingError> {
    let mut bad_indices = vec![];
    for (idx, range) in ranges.iter().enumerate() {
        if (range
            .phase_fit
            .as_ref()
            .ok_or_else(|| {
                FittingError::BadFit("Cannot filter fits since phase not fit".to_string())
            })?
            .slope
            == 0.0)
            || (range
                .lin_pwr_fit
                .as_ref()
                .ok_or_else(|| {
                    FittingError::BadFit(
                        "Cannot filter fits since power not linearly fit".to_string(),
                    )
                })?
                .slope
                == 0.0)
            || (range
                .quad_pwr_fit
                .as_ref()
                .ok_or_else(|| {
                    FittingError::BadFit(
                        "Cannot filter fits since power not quadratically fit".to_string(),
                    )
                })?
                .slope
                == 0.0)
        {
            bad_indices.push(idx);
        }
    }
    bad_indices.sort();
    for idx in bad_indices.iter().rev() {
        ranges.remove(*idx);
    }
    Ok(())
}
