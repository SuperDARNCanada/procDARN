use dmap::formats::{RawacfRecord};
use crate::fitting::fitacf3::fitstruct::{LagNode, RangeNode};

pub fn mark_bad_samples(rec: &RawacfRecord) -> Vec<i32> {
    let mut pulses_in_us: Vec<i16> = rec.pulse_table.data
        .iter()
        .map(|p| p * rec.multi_pulse_increment)
        .collect();

    if rec.offset != 0 {
        if rec.channel == 1 {
            let pulses_stereo: Vec<i16> = pulses_in_us
                .iter()
                .map(|p| p - rec.offset)
                .collect();
            pulses_in_us.extend(pulses_stereo);
        } else if rec.channel == 2 {
            let pulses_stereo: Vec<i16> = pulses_in_us
                .iter()
                .map(|p| p + rec.offset)
                .collect();
            pulses_in_us.extend(pulses_stereo);
        }
    }
    pulses_in_us.sort();

    let mut ts = rec.lag_to_first_range;
    let mut t1 = 0;
    let mut t2 = 0;
    let mut sample = 0;
    let mut bad_samples = vec![];

    for pulse_us in pulses_in_us {
        t1 = pulse_us - rec.tx_pulse_length / 2;
        t2 = t1 + 3 * rec.tx_pulse_length / 2 + 100;

        // Start incrementing the sample until we find a sample that lies within a pulse
        while ts < t1 {
            sample += 1;
            ts = ts + rec.sample_separation;
        }

        // Blank all samples within the pulse duration
        while (ts >= t1) && (ts <= t2) {
            bad_samples.push(sample);
            sample += 1;
            ts += rec.sample_separation;
        }
    }
    bad_samples
}

pub fn filter_tx_overlapped_lags(rec: &RawacfRecord, lags: Vec<LagNode>, ranges: Vec<RangeNode>) {
    let bad_samples = mark_bad_samples(rec);
    for mut range_node in ranges {
        let mut bad_indices = vec![];
        for idx in 0..lags.len() {
            let lag = &lags[idx];
            let sample_1 = lag.sample_base_1 + range_node.range_num;
            let sample_2 = lag.sample_base_2 + range_node.range_num;
            if bad_samples.contains(&sample_1) || bad_samples.contains(&sample_2) {
                bad_indices.push(idx);
            }
        }
        bad_indices.iter()
            .map(|i| {
                range_node.powers.remove(*i);
                range_node.phases.remove(*i);
                range_node.elev.remove(*i);
                range_node.alpha_2.remove(*i);
            });
    }
}

pub fn filter_infinite_lags(ranges: Vec<RangeNode>) {
    for range in ranges {
        let mut infinite_indices = vec![];
        for i in 0..range.powers.len() {
            if !range.powers[i].ln_power.is_finite() {
                infinite_indices.push(i);
            }
        }
    }
}

pub fn filter_low_power_lags(rec: &RawacfRecord, ranges: Vec<RangeNode>) {
    if rec.num_averages <= 0 {
        return
    }
    for mut range in ranges {
        let range_num = range.range_num;
        if range.powers.len() == 0 { continue }
        let log_sigma_fluc = (FLUCTUATION_CUTOFF_COEFFICIENT * &rec.lag_zero_power / ((2 * rec.num_averages) as f32).sqrt()).ln();
        let mut bad_indices = vec![];
        let mut cutoff_lag = rec.num_lags as usize + 1;

        for idx in 0..range.powers.len() {
            if idx > cutoff_lag as usize {
                bad_indices.push(idx);
            } else {
                let log_power = range.powers[idx].ln_power;
                let alpha_2 = range.alpha_2[idx].alpha_2;
                if ((1 as f64 / alpha_2.sqrt()) <= ALPHA_CUTOFF) &&
                    ((log_power < log_sigma_fluc) || is_close!(log_power, log_sigma_fluc)) {
                    cutoff_lag = idx;
                    bad_indices.push(idx);
                }
            }
        }
        for idx in bad_indices {
            range.powers.remove(idx);
        }
    }
}

pub fn filter_bad_acfs(rec: &RawacfRecord, mut ranges: Vec<RangeNode>, noise_power: f32) {
    if rec.num_averages <= 0 {
        return
    }
    let cutoff_power = noise_power * 2.0;
    let mut bad_indices = vec![];
    for idx in 0..ranges.len() {
        let range = &ranges[idx];
        let range_num = range.range_num as usize;
        let power = rec.lag_zero_power.data[range_num];
        let num_powers = range.powers.len();
        if (power <= cutoff_power) || (num_powers < MIN_LAGS) {
            bad_indices.push(idx);
        } else {
            let power_value = range.powers[0].ln_power;
            let mut all_equal = true;
            for pwr in range.powers {
                if !isclose!(pwr.ln_power, power_value) { all_equal = false; }
            }
            if all_equal { bad_indices.push(idx); }
        }
    }
    bad_indices.sort();
    for idx in bad_indices.iter().rev() {
        ranges.remove(*idx);
    }
}

pub fn filter_bad_fits(mut ranges: Vec<RangeNode>) {
    let mut bad_indices = vec![];
    for idx in 0..ranges.len() {
        let range = &ranges[idx];
        if (range.phase_fit.b == 0.0) ||
            (range.lin_pwr_fit.b == 0.0) ||
            (range.quad_pwr_fit.b == 0.0) {
            bad_indices.push(idx);
        }
    }
    bad_indices.sort();
    for idx in bad_indices.iter().rev() {
        ranges.remove(*idx);
    }
}