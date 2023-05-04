use crate::fitting::fitacf3::fitacf_v3::{
    Fitacf3Error, ALPHA_CUTOFF, FLUCTUATION_CUTOFF_COEFFICIENT, MIN_LAGS,
};
use crate::fitting::fitacf3::fitstruct::{LagNode, RangeNode};
use dmap::formats::RawacfRecord;
use is_close::is_close;

/// passing
pub fn mark_bad_samples(rec: &RawacfRecord) -> Vec<i32> {
    let mut pulses_in_us: Vec<i32> = rec
        .pulse_table
        .data
        .iter()
        .map(|&p| p as i32 * rec.multi_pulse_increment as i32)
        .collect();

    if rec.offset != 0 {
        if rec.channel == 1 {
            let pulses_stereo: Vec<i32> = pulses_in_us.iter().map(|&p| p - rec.offset as i32).collect();
            pulses_in_us.extend(pulses_stereo);
        } else if rec.channel == 2 {
            let pulses_stereo: Vec<i32> = pulses_in_us.iter().map(|&p| p + rec.offset as i32).collect();
            pulses_in_us.extend(pulses_stereo);
        }
    }
    pulses_in_us.sort();

    let mut ts = rec.lag_to_first_range as i32;
    let mut t1;
    let mut t2;
    let mut sample = 0;
    let mut bad_samples = vec![];

    for pulse_us in pulses_in_us {
        t1 = pulse_us - rec.tx_pulse_length as i32 / 2;
        t2 = t1 + 3 * rec.tx_pulse_length as i32 / 2 + 100;

        // Start incrementing the sample until we find a sample that lies within a pulse
        while ts < t1 {
            sample += 1;
            ts = ts + rec.sample_separation as i32;
        }

        // Blank all samples within the pulse duration
        while (ts >= t1) && (ts <= t2) {
            bad_samples.push(sample);
            sample += 1;
            ts += rec.sample_separation as i32;
        }
    }
    bad_samples
}

/// passing
pub fn filter_tx_overlapped_lags(
    rec: &RawacfRecord,
    lags: Vec<LagNode>,
    ranges: &mut Vec<RangeNode>,
) {
    let bad_samples = mark_bad_samples(rec);
    for range_node in ranges {
        let mut bad_indices = vec![];
        for idx in 0..lags.len() {
            let lag = &lags[idx];
            let sample_1 = lag.sample_base_1 + range_node.range_num as i32;
            let sample_2 = lag.sample_base_2 + range_node.range_num as i32;
            if bad_samples.contains(&sample_1) || bad_samples.contains(&sample_2) {
                bad_indices.push(idx);
            }
        }
        for i in bad_indices.iter().rev() {
            range_node.powers.remove(*i);
            range_node.phases.remove(*i);
            range_node.elev.remove(*i);
            range_node.alpha_2.remove(*i);
        }
    }
}

/// passing
pub fn filter_infinite_lags(ranges: &mut Vec<RangeNode>) {
    for range in ranges {
        let mut infinite_indices = vec![];
        for i in 0..range.powers.ln_power.len() {
            if !range.powers.ln_power[i].is_finite() {
                infinite_indices.push(i);
            }
        }
        for i in infinite_indices.iter().rev() {
            range.powers.remove(*i);
            range.phases.remove(*i);
            range.elev.remove(*i);
            range.alpha_2.remove(*i);
        }
    }
}

/// passing
pub fn filter_low_power_lags(rec: &RawacfRecord, ranges: &mut Vec<RangeNode>) {
    if rec.num_averages <= 0 {
        return;
    }
    for range in ranges {
        let range_num = range.range_num;
        if range.powers.ln_power.len() == 0 {
            continue;
        }
        let log_sigma_fluc = (FLUCTUATION_CUTOFF_COEFFICIENT as f32
            * rec.lag_zero_power.data[range_num]
            / ((2 * rec.num_averages) as f32).sqrt())
        .ln();
        let mut bad_indices = vec![];
        let mut cutoff_lag = rec.num_lags as usize + 1;

        for idx in 0..range.powers.ln_power.len() {
            if idx > cutoff_lag as usize {
                bad_indices.push(idx);
            } else {
                let log_power = range.powers.ln_power[idx];
                let alpha_2 = range.alpha_2[idx];
                if ((1 as f64 / alpha_2.sqrt()) <= ALPHA_CUTOFF as f64)
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
            range.phases.remove(*i);
            range.elev.remove(*i);
            range.alpha_2.remove(*i);
        }
    }
}

/// passing
pub fn filter_bad_acfs(rec: &RawacfRecord, ranges: &mut Vec<RangeNode>, noise_power: f32) {
    if rec.num_averages <= 0 {
        return;
    }
    let cutoff_power = noise_power * 2.0;
    let mut bad_indices = vec![];
    for idx in 0..ranges.len() {
        let range = &ranges[idx];
        let range_num = range.range_num as usize;
        let power = rec.lag_zero_power.data[range_num];
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

pub fn filter_bad_fits(ranges: &mut Vec<RangeNode>) -> Result<(), Fitacf3Error> {
    let mut bad_indices = vec![];
    for idx in 0..ranges.len() {
        let range = &ranges[idx];
        if (range
            .phase_fit
            .as_ref()
            .ok_or_else(|| {
                Fitacf3Error::Message(format!("Cannot filter fits since phase not fit"))
            })?
            .slope
            == 0.0)
            || (range
                .lin_pwr_fit
                .as_ref()
                .ok_or_else(|| {
                    Fitacf3Error::Message(format!(
                        "Cannot filter fits since power not linearly fit"
                    ))
                })?
                .slope
                == 0.0)
            || (range
                .quad_pwr_fit
                .as_ref()
                .ok_or_else(|| {
                    Fitacf3Error::Message(format!(
                        "Cannot filter fits since power not quadratically fit"
                    ))
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
