use crate::fitting::common::error::FittingError;
use crate::fitting::lmfit2::fitstruct::RangeNode;
use crate::utils::constants::{KHZ_TO_HZ, LIGHTSPEED};
use crate::utils::rawacf::Rawacf;
use numpy::ndarray::Array1;
use std::f64::consts::PI;
use std::iter::zip;

/// Estimate the self clutter for each range gate, from each sample in each lag
pub(crate) fn estimate_self_clutter(range_list: &mut Vec<RangeNode>, rawacf: &Rawacf) {
    for range in range_list {
        range.self_clutter = Some(estimate_maximum_self_clutter(
            range,
            rawacf,
            &rawacf.pwr0,
        ));
    }
}

/// Maximal lag0 power-based Self-Clutter Estimator
fn estimate_maximum_self_clutter(
    range: &RangeNode,
    rawacf: &Rawacf,
    lag0_power: &Array1<f32>,
) -> Vec<f64> {
    let pulse_width = rawacf.mpinc / rawacf.smsep;
    let first_range_sample = rawacf.lagfr / rawacf.smsep;

    let bad_range = rawacf.nrang;
    let mut self_clutter: Vec<f64> = vec![];
    let mut r1: Array1<i16> = Array1::ones(rawacf.mppul as usize) * -1000;
    let mut r2 = r1.clone();
    for lag in 0..rawacf.mplgs as usize {

        let sample_1 = pulse_width * rawacf.ltab[[lag, 0]] + range.range_num as i16 + first_range_sample;
        let sample_2 = pulse_width * rawacf.ltab[[lag, 1]] + range.range_num as i16 + first_range_sample;

        for pulse in 0..rawacf.mppul as usize {
            // Find the pulses that were transmitted before the samples were recorded,
            // then save which range gates each pulse is coming from
            if rawacf.ptab[pulse] * pulse_width <= sample_1 {
                let temp = sample_1 - rawacf.ptab[pulse] * pulse_width - first_range_sample;
                // Also we need to check and make sure we only save interfering range
                // gates where we have valid lag0 power
                if (temp != range.range_num as i16)
                    && (temp >= 0)
                    && (temp < rawacf.nrang)
                    && (temp < bad_range)
                {
                    r1[pulse] = temp;
                }
            }
            // Do the same for the second sample comprising the lag
            if rawacf.ptab[pulse] * pulse_width <= sample_2 {
                let temp = sample_2 - rawacf.ptab[pulse] * pulse_width - first_range_sample;
                if (temp != range.range_num as i16)
                    && (temp >= 0)
                    && (temp < rawacf.nrang)
                    && (temp < bad_range)
                {
                    r2[pulse] = temp;
                }
            }
        }

        let (mut term1, mut term2, mut term3) = (0.0_f64, 0.0_f64, 0.0_f64);

        for pulse in 0..rawacf.mppul as usize {
            // First term in the summation for the self-clutter estimate (P_r*P_j^*)
            if r2[pulse] != -1000 {
                term1 += (lag0_power[range.range_num as usize] * lag0_power[r2[pulse] as usize]).sqrt()
                    as f64;
            }
            // Second term in the summation for the self-clutter estimate (P_i*P_r^*)
            if r1[pulse] != -1000 {
                term2 += (lag0_power[range.range_num as usize] * lag0_power[r1[pulse] as usize]).sqrt()
                    as f64;
            }
            for pulse2 in 0..rawacf.mppul as usize {
                // Third term in the summation for the self-clutter estimate (P_i*P_j^*)
                if (r1[pulse] != -1000) && (r2[pulse2] != -1000) {
                    term3 += (lag0_power[r1[pulse] as usize] * lag0_power[r2[pulse2] as usize])
                        .sqrt() as f64;
                }
            }
        }
        self_clutter.push(term1 + term2 + term3);
    }
    self_clutter
}

/// Estimate the first-order error for each lag at each range
pub(crate) fn estimate_first_order_error(
    range_list: &mut Vec<RangeNode>,
    rawacf: &Rawacf,
    noise_power: f64,
) {
    for range in range_list {
        let mut sigma_real = vec![];
        let mut sigma_imag = vec![];
        if let Some(clutter) = &range.self_clutter {
            for sc in clutter.iter() {
                let error = (rawacf.pwr0[range.range_num as usize] as f64 + noise_power + sc)
                    / (rawacf.nave as f64).sqrt();
                sigma_real.push(error);
                sigma_imag.push(error);
            }
        }
        range.sigma_real = Some(sigma_real);
        range.sigma_imag = Some(sigma_imag);
        // todo: What if clutter is not present?
    }
}

/// Estimate the error for real and imaginary components of each lag at each range
pub(crate) fn estimate_real_imag_error(
    range_list: &mut [RangeNode],
    rawacf: &Rawacf,
    noise_power: f64,
) -> Result<(), FittingError> {
    let wavelength: f64 = LIGHTSPEED as f64 / (rawacf.tfreq as f64 * KHZ_TO_HZ as f64);

    for range in range_list.iter_mut() {
        let power = range
            .lin_fit
            .as_ref()
            .ok_or_else(|| FittingError::BadFit("No linear fit present".to_string()))?
            .pwr;
        let width = range
            .lin_fit
            .as_ref()
            .ok_or_else(|| FittingError::BadFit("No linear fit present".to_string()))?
            .wid;
        let vel = range
            .lin_fit
            .as_ref()
            .ok_or_else(|| FittingError::BadFit("No linear fit present".to_string()))?
            .vel;

        let mut sigma_real = vec![];
        let mut sigma_imag = vec![];

        if let Some(clutter) = &range.self_clutter {
            for (sc, t) in zip(clutter.iter(), range.t.iter()) {
                let mut rho = (-2.0_f64 * PI * width * t / wavelength).exp();

                if rho > 0.999 {
                    rho = 0.999;
                }
                rho *= power / (power + noise_power + sc);
                let rho_real = rho * (4.0 * PI * vel * t / wavelength).cos();
                let rho_imag = rho * (4.0 * PI * vel * t / wavelength).sin();

                let real_error = (power + noise_power + sc)
                    * (((1.0 - rho * rho) / 2.0 + (rho_real * rho_real)) / (rawacf.nave as f64))
                        .sqrt();
                let imag_error = (power + noise_power + sc)
                    * (((1.0 - rho * rho) / 2.0 + (rho_imag * rho_imag)) / (rawacf.nave as f64))
                        .sqrt();
                sigma_real.push(real_error);
                sigma_imag.push(imag_error);
            }
        }
        range.sigma_real = Some(sigma_real);
        range.sigma_imag = Some(sigma_imag);
        // todo: What if clutter is not present?
    }
    Ok(())
}
