use std::slice::range;
use dmap::formats::RawacfRecord;

pub struct RangeNode {
    pub range_num: i32,
    pub range_idx: i32,
    pub cross_range_interference: Vec<f64>,
    pub refractive_idx: f32,
    pub alpha_2: Vec<Alpha>,
    pub phases: Vec<PhaseNode>,
    pub powers: Vec<PowerNode>,
    pub elev: Vec<f32>,
    pub lin_pwr_fit: FittedData,
    pub quad_pwr_fit: FittedData,
    pub lin_pwr_fit_err: FittedData,
    pub quad_pwr_fit_err: FittedData,
    pub phase_fit: FittedData,
    pub elev_fit: FittedData,
}
impl RangeNode {
    pub fn new(index: i32, range_num: i32, record: &RawacfRecord, lags: Vec<[i32; 2]>) -> RangeNode {
        let cross_range_interference = RangeNode::calculate_cross_range_interference(range_num as i16, record);
        let alpha_2 = calculate_alpha_2();
        let phases = Phase
        RangeNode {
            range_idx: index,
            range_num,

        }
    }
    fn calculate_cross_range_interference(range_num: i32, rec: &RawacfRecord) -> Vec<f32> {
        let tau: i16;
        if rec.sample_separation != 0 {
            tau = rec.multi_pulse_increment / rec.sample_separation;
        } else {
            // TODO: Log warning?
            tau = rec.multi_pulse_increment / rec.tx_pulse_length;
        }

        let interference_for_pulses = vec![];
        for pulse_to_check in 0..rec.num_pulses as usize {
            let total_interference = 0.0;
            for pulse in 0..rec.num_pulses as usize {
                let pulse_diff = rec.pulse_table.data[pulse_to_check] - rec.pulse_table.data[pulse];
                let range_to_check = (pulse_diff * tau + range_num) as usize;
                if (pulse != pulse_to_check) &&
                    (0 <= range_to_check) &&
                    (range_to_check < rec.num_ranges as usize) {
                    total_interference += rec.lag_zero_power.data[range_to_check];
                }
            }
            interference_for_pulses[pulse_to_check] = total_interference;
        }
        interference_for_pulses
    }
    fn calculate_alphas(range_num: i32, cross_range_interference: Vec<f32>, rec: &RawacfRecord, lags: Vec<LagNode>) {
        let alpha_2 = vec![];
        for idx in 0..lags.len() {
            let lag = lags[idx];
            let pulse_1_interference = cross_range_interference[lag.pulses[0] as usize];
            let pulse_2_interference = cross_range_interference[lag.pulses[1] as usize];
            let lag_zero_power = rec.lag_zero_power.data[range_num as usize];
            alpha_2.push(lag_zero_power*lag_zero_power / (()))
        }
    }
}

struct PhaseNode {
    pub phi: f64,
    pub t: f64,
    pub std_dev: f64,
    pub lag_idx: i32, // TODO: Is this redundant with Alpha?
    pub alpha_2: f64, // TODO: Is this redundant with Alpha?
}

struct PowerNode {
    pub ln_power: f64,
    pub t: f64,
    pub std_dev: f64,
    pub lag_idx: i32, // TODO: Is this redundant with Alpha?
    pub alpha_2: f64, // TODO: Is this redundant with Alpha?
}

pub struct LagNode {
    pub lag_num: i32,
    pub pulses: [i32; 2],
    pub lag_idx: i32,
    pub sample_base_1: i32,
    pub sample_base_2: i32,
}

pub struct Alpha {
    lag_idx: i32,
    pub alpha_2: f64,
}

enum IntfPosition {
    Forward,
    Behind,
}

pub struct FitData {
    pub channel: i32,
    pub offset: i32, /* used for stereo badlags */
    pub cp: i32,
    pub xcf_flag: i32,
    pub transmit_freq: i32,
    pub noise: f32,
    pub num_ranges: i32,
    pub sample_sep: i32,
    pub num_avg: i32,
    pub num_lags: i32,
    pub multi_pulse_increment: i32,
    pub tx_pulse_len: i32,
    pub lag_to_first_range: i32,
    pub num_pulses: i32,
    pub beam_num: i32,
    pub old: i32,
    pub lag: [i32; 2],
    pub pulse: i32,
    pub pwr0: f64,
    pub acfd: Vec<f64>,
    pub xcfd: Vec<f64>,
    pub maxbeam: i32,
    pub beam_offset: f64,
    pub beam_sep: f64,
    pub interferometer_offset: [f64; 3],
    pub phi_diff: f64,
    pub time_diff: f64,
    pub vel_dir: f64,
    pub time: Time,
}

struct Time {
    yr: i16,
    mo: i16,
    dy: i16,
    hr: i16,
    mt: i16,
    sc: i16,
    us: i32,
}

struct FittedData {
    pub delta: f64,
    pub a: f64,
    pub b: f64,
    pub variance_a: f64,
    pub variance_b: f64,
    pub delta_a: f64,
    pub delta_b: f64,
    pub covariance_ab: f64,
    pub residual_ab: f64,
    pub quality: f64,
    pub chi_squared: f64,
    pub sums: Sums,
}

struct Sums {
    sum: f64,
    sum_x: f64,
    sum_y: f64,
    sum_x_squared: f64,
    sum_y_squared: f64,
}

enum FitType {
    Linear,
    Quadratic,
}
