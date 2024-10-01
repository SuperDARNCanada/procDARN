use crate::fitting::lmfit2::fitstruct::{LagNode, RangeNode};
use crate::utils::rawacf::Rawacf;
use std::f64::consts::PI;

pub const ACF_SNR_CUTOFF: f64 = 1.0;

/// Creates the lag table based on the data.
pub(crate) fn create_lag_list(record: &Rawacf) -> Vec<LagNode> {
    let lag_table = &record.ltab;
    let pulse_table = &record.ptab;
    let multi_pulse_increment = record.mpinc;
    let sample_separation = record.smsep;

    let mut lags = vec![];
    for i in 0..record.mplgs as usize {
        let mut pulse_1_idx = 0;
        let mut pulse_2_idx = 0;
        let number = lag_table[[i, 1]] - lag_table[[i, 0]];
        for j in 0..record.mppul as usize {
            if lag_table[[i, 0]] == pulse_table[j] {
                pulse_1_idx = j;
            }
            if lag_table[[i, 1]] == pulse_table[j] {
                pulse_2_idx = j;
            }
        }
        let sample_base_1 =
            i32::from(lag_table[[i, 0]] * (multi_pulse_increment / sample_separation));
        let sample_base_2 =
            i32::from(lag_table[[i, 1]] * (multi_pulse_increment / sample_separation));
        lags.push(LagNode {
            lag_num: i32::from(number),
            pulses: [pulse_1_idx, pulse_2_idx],
            sample_base_1,
            sample_base_2,
        });
    }
    lags
}

/// Calculates the minimum power value for ACFs in the record
pub(crate) fn acf_cutoff_power(rec: &Rawacf) -> f32 {
    let mut sorted_power_levels = rec.pwr0.clone().to_vec();
    sorted_power_levels.sort_by(f32::total_cmp); // sort floats
    let mut i: usize = 0;
    let mut j: f64 = 0.0;
    let mut min_power: f64 = 0.0;
    while j < 10.0 && i < rec.nrang as usize / 3 {
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
    let search_noise = rec.noise_search;
    if min_power < ACF_SNR_CUTOFF && search_noise != 0.0 {
        min_power = search_noise as f64;
    }
    min_power as f32
}

/// Applies a correction to the noise power estimate to account for selecting least-powerful ranges
pub(crate) fn cutoff_power_correction(rec: &Rawacf) -> f64 {
    let std_dev = 1.0 / (rec.nave as f64).sqrt();

    let mut i = 0.0_f64;
    let mut cumulative_pdf = 0.0_f64;
    let mut cumulative_pdf_x_norm_power = 0.0_f64;
    let mut normalized_power: f64;
    while cumulative_pdf < (10.0 / rec.nrang as f64) {
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

/// Finds all samples that were collected during transmission of a pulse.
pub(crate) fn mark_bad_samples(rec: &Rawacf) -> Vec<i32> {
    let mut pulses_in_us: Vec<i32> = rec
        .ptab
        .iter()
        .map(|&p| i32::from(p) * i32::from(rec.mpinc))
        .collect();

    if rec.offset != 0 {
        if rec.channel == 1 {
            let pulses_stereo: Vec<i32> = pulses_in_us
                .iter()
                .map(|&p| p - i32::from(rec.offset))
                .collect();
            pulses_in_us.extend(pulses_stereo);
        } else if rec.channel == 2 {
            let pulses_stereo: Vec<i32> = pulses_in_us
                .iter()
                .map(|&p| p + i32::from(rec.offset))
                .collect();
            pulses_in_us.extend(pulses_stereo);
        }
    }
    pulses_in_us.sort();

    let mut ts = i32::from(rec.lagfr);
    let mut t1;
    let mut t2;
    let mut sample = 0;
    let mut bad_samples = vec![];

    for pulse_us in pulses_in_us {
        t1 = pulse_us - i32::from(rec.txpl) / 2;
        t2 = t1 + 3 * i32::from(rec.txpl) / 2 + 100;

        // Start incrementing the sample until we find a sample that lies within a pulse
        while ts < t1 {
            sample += 1;
            ts += i32::from(rec.smsep);
        }

        // Blank all samples within the pulse duration
        while (ts >= t1) && (ts <= t2) {
            bad_samples.push(sample);
            sample += 1;
            ts += i32::from(rec.smsep);
        }
    }
    bad_samples
}

/// Removes all lags that contain samples collected during transmission of a pulse.
pub(crate) fn remove_tx_overlapped_lags(
    rec: &Rawacf,
    lags: &[LagNode],
    ranges: &mut Vec<RangeNode>,
) {
    let bad_samples = mark_bad_samples(rec);
    for range_node in ranges.iter_mut() {
        let mut bad_indices = vec![];
        for (idx, lag) in lags.iter().enumerate() {
            let sample_1 = lag.sample_base_1 + range_node.range_num as i32;
            let sample_2 = lag.sample_base_2 + range_node.range_num as i32;
            if bad_samples.contains(&sample_1) || bad_samples.contains(&sample_2) {
                bad_indices.push(idx);
            }
        }
        for i in bad_indices.iter().rev() {
            range_node.acf_real.remove(*i);
            range_node.acf_imag.remove(*i);
            range_node.t.remove(*i);
            if let Some(ref mut x) = range_node.sigma_real {
                x.remove(*i);
            }
            if let Some(ref mut x) = range_node.sigma_imag {
                x.remove(*i);
            }
        }
    }
}
