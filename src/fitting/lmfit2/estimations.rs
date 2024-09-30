use std::f64::consts::PI;
use std::iter::zip;
use itertools::enumerate;
use ndarray::Array1;
use crate::fitting::common::error::FittingError;
use crate::fitting::common::fitstruct::RangeNode;
use crate::utils::rawacf::Rawacf;

const LIGHTSPEED: f64 = 299_792_458.0;
const KHZ_TO_HZ: f64 = 1000.0;

/// Estimate the self clutter for each range gate, from each sample in each lag
pub(crate) fn estimate_self_clutter(range_list: &mut Vec<RangeNode>, rawacf: &Rawacf) {
    for range in range_list {
        range.self_clutter = Some(estimate_maximum_self_clutter(range.range_num, rawacf, &rawacf.pwr0));
    }
}

/// Maximal lag0 Power Based Self-Clutter Estimator
fn estimate_maximum_self_clutter(range_gate: u16, rawacf: &Rawacf, lag0_power: &Array1<f32>) -> Vec<f64> {
    let pulse_width = rawacf.mpinc / rawacf.smsep;
    let first_range_sample = rawacf.lagfr / rawacf.smsep;

    let bad_range = rawacf.nrang;
    let mut self_clutter: Vec<f64> = vec![];
    let mut r1 = Array1::ones(rawacf.mppul) * -1000;
    let mut r2 = r1.clone();
    for lag in 0..rawacf.mplgs {
        self_clutter.push(0.0);

        let sample_1 = pulse_width * rawacf.ltab[[lag, 0]] + range_gate as i16 + first_range_sample;
        let sample_2 = pulse_width * rawacf.ltab[[lag, 1]] + range_gate as i16 + first_range_sample;

        for pulse in 0..rawacf.mppul {
            // Find the pulses that were transmitted before the samples were recorded,
            // then save which range gates each pulse is coming from
            if rawacf.ptab[pulse] * pulse_width <= sample_1 {
                let temp = sample_1 - rawacf.ptab[pulse] * pulse_width - first_range_sample;
                // Also we need to check and make sure we only save interfering range
                // gates where we have valid lag0 power
                if (temp != range_gate as i16) && (temp >= 0) && (temp < rawacf.nrang) && (temp < bad_range) {
                    r1[pulse] = temp;
                }
            }
            // Do the same for the second sample comprising the lag
            if rawacf.ptab[pulse] * pulse_width <= sample_2 {
                let temp = sample_2 - rawacf.ptab[pulse] * pulse_width - first_range_sample;
                if (temp != range_gate as i16) && (temp >= 0) && (temp < rawacf.nrang) && (temp < bad_range) {
                    r2[pulse] = temp;
                }
            }
        }

        let (mut term1, mut term2, mut term3) = (0.0_f64, 0.0_f64, 0.0_f64);

        for pulse in 0..rawacf.mppul {
            // First term in the summation for the self-clutter estimate (P_r*P_j^*)
            if r2[pulse] != -1000 {
                term1 += (lag0_power[range_gate] * lag0_power[r2[pulse]]).sqrt() as f64;
            }
            // Second term in the summation for the self-clutter estimate (P_i*P_r^*)
            if r1[pulse] != -1000 {
                term2 += (lag0_power[range_gate] * lag0_power[r1[pulse]]).sqrt() as f64;
            }
            for pulse2 in 0..rawacf.mppul {
                // Third term in the summation for the self-clutter estimate (P_i*P_j^*)
                if (r1[pulse] != -1000) && (r2[pulse2] != -1000) {
                    term3 += (lag0_power[r1[pulse]] * lag0_power[r2[pulse2]]).sqrt() as f64;
                }
            }
        }
        self_clutter.push(term1 + term2 + term3);
    }
    self_clutter
}

/// Estimate the first-order error for each lag at each range
pub(crate) fn estimate_first_order_error(range_list: &mut Vec<RangeNode>, rawacf: &Rawacf, noise_power: f64) {
    for range in range_list {
        if let Some(clutter) = &range.self_clutter {
            for (lag, sc) in enumerate(clutter.iter()) {
                let error = (rawacf.pwr0[range.range_num] as f64 + noise_power + sc) / (rawacf.nave as f64).sqrt();
                range.powers.std_dev[lag] = error;
                range.phases.std_dev_real[lag] = error;
                range.phases.std_dev_imag[lag] = error;
            }
        }
        // todo: What if clutter is not present?
    }
}

/// Estimate the error for real and imaginary components of each lag at each range
pub(crate) fn estimate_real_imag_error(range_list: &mut Vec<RangeNode>, rawacf: &Rawacf, noise_power: f64) -> Result<(), FittingError> {
    let wavelength: f64 = LIGHTSPEED / (rawacf.tfreq as f64 * KHZ_TO_HZ);

    for range in range_list.iter_mut() {
        let power = range.lin_pwr_fit.ok_or_else(|| FittingError::BadFit("No linear fit present".to_string()))?.intercept;
        // todo: Get the right values for width and vel
        let width = range.lin_pwr_fit.ok_or_else(|| FittingError::BadFit("No linear fit present".to_string()))?.intercept;
        let vel = range.lin_pwr_fit.ok_or_else(|| FittingError::BadFit("No linear fit present".to_string()))?.intercept;

        if let Some(clutter) = &range.self_clutter {
            for (lag, (sc, t)) in enumerate(zip(clutter.iter(), range.phases.t.iter())) {
                let mut rho = (-2.0_f64 * PI *  * t / wavelength).exp();

                if rho > 0.999 { rho = 0.999; }
                rho *= power / (power + noise_power + sc);
                let rho_real = rho * (4.0 * PI * vel * t / wavelength).cos();
                let rho_imag = rho * (4.0 * PI * vel * t / wavelength).sin();

                let real_error = (power + noise_power + sc) * (((1 - rho*rho) / 2.0 + (rho_real*rho_real)) / (rawacf.nave as f64)).sqrt();
                let imag_error = (power + noise_power + sc) * (((1 - rho*rho) / 2.0 + (rho_imag*rho_imag)) / (rawacf.nave as f64)).sqrt();

                range.phases.std_dev_real[lag] = real_error;
                range.phases.std_dev_imag[lag] = imag_error;
            }
        }
        // todo: What if clutter is not present?
    }
    Ok(())
}