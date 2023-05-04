use crate::fitting::fitacf3::fitacf_v3::Fitacf3Error;
use dmap::formats::RawacfRecord;
use std::iter::zip;

#[derive(Debug)]
pub struct RangeNode {
    pub range_num: usize,
    pub range_idx: usize,
    pub cross_range_interference: Vec<f64>,
    pub refractive_idx: f32,
    pub power_alpha_2: Vec<f64>,
    pub phase_alpha_2: Vec<f64>,
    pub phases: PhaseNode,
    pub powers: PowerNode,
    pub elev: PhaseNode,
    pub lin_pwr_fit: Option<FittedData>,
    pub quad_pwr_fit: Option<FittedData>,
    pub lin_pwr_fit_err: Option<FittedData>,
    pub quad_pwr_fit_err: Option<FittedData>,
    pub phase_fit: Option<FittedData>,
    pub elev_fit: Option<FittedData>,
}
impl RangeNode {
    pub fn new(
        index: usize,
        range_num: usize,
        record: &RawacfRecord,
        lags: &Vec<LagNode>,
    ) -> Result<RangeNode, Fitacf3Error> {
        let cross_range_interference =
            RangeNode::calculate_cross_range_interference(range_num, record);
        let alpha_2 =
            RangeNode::calculate_alphas(range_num, &cross_range_interference, record, &lags);
        let phases = PhaseNode::new(record, "acfd", &lags, index)?;
        let elevations = PhaseNode::new(record, "xcfd", &lags, index)?;
        let powers = PowerNode::new(record, &lags, index, range_num, &alpha_2);
        Ok(RangeNode {
            range_idx: index,
            range_num,
            cross_range_interference,
            refractive_idx: 1.0,
            power_alpha_2: alpha_2.clone(),
            phase_alpha_2: alpha_2,
            phases,
            powers,
            elev: elevations,
            lin_pwr_fit: None,
            quad_pwr_fit: None,
            lin_pwr_fit_err: None,
            quad_pwr_fit_err: None,
            phase_fit: None,
            elev_fit: None,
        })
    }
    fn calculate_cross_range_interference(range_num: usize, rec: &RawacfRecord) -> Vec<f64> {
        let tau: i16;
        if rec.sample_separation != 0 {
            tau = rec.multi_pulse_increment / rec.sample_separation;
        } else {
            // TODO: Log warning?
            tau = rec.multi_pulse_increment / rec.tx_pulse_length;
        }

        let mut interference_for_pulses: Vec<f64> = vec![];
        for pulse_to_check in 0..rec.num_pulses as usize {
            let mut total_interference: f64 = 0.0;
            for pulse in 0..rec.num_pulses as usize {
                let pulse_diff = rec.pulse_table.data[pulse_to_check] - rec.pulse_table.data[pulse];
                let range_to_check = (pulse_diff * tau + range_num as i16) as usize;
                if (pulse != pulse_to_check) && (range_to_check < rec.num_ranges as usize) {
                    total_interference += rec.lag_zero_power.data[range_to_check] as f64;
                }
            }
            interference_for_pulses.push(total_interference);
        }
        interference_for_pulses
    }
    fn calculate_alphas(
        range_num: usize,
        cross_range_interference: &Vec<f64>,
        rec: &RawacfRecord,
        lags: &Vec<LagNode>,
    ) -> Vec<f64> {
        let mut alpha_2: Vec<f64> = vec![];
        for idx in 0..lags.len() {
            let lag = &lags[idx];
            let pulse_1_interference = cross_range_interference[lag.pulses[0] as usize];
            let pulse_2_interference = cross_range_interference[lag.pulses[1] as usize];
            let lag_zero_power = rec.lag_zero_power.data[range_num] as f64;
            alpha_2.push(
                lag_zero_power * lag_zero_power
                    / ((lag_zero_power + pulse_1_interference)
                        * (lag_zero_power + pulse_2_interference)),
            );
        }
        alpha_2
    }
}

#[derive(Debug)]
pub struct PhaseNode {
    pub phases: Vec<f64>,
    pub t: Vec<f64>,
    pub std_dev: Vec<f64>,
}
impl PhaseNode {
    pub fn new(
        rec: &RawacfRecord,
        phase_type: &str,
        lags: &Vec<LagNode>,
        range_idx: usize,
    ) -> Result<PhaseNode, Fitacf3Error> {
        let acfd = match phase_type {
            "acfd" => &rec.acfs.data,
            "xcfd" => match &rec.xcfs {
                Some(x) => &x.data,
                None => Err(Fitacf3Error::Message(format!("Cannot find xcfs in data")))?,
            },
            _ => Err(Fitacf3Error::Message(format!(
                "Unknown type for PhaseNode: {}",
                phase_type
            )))?,
        };
        let start_idx = range_idx * 2 * rec.num_lags as usize;
        let end_idx = start_idx + 2 * rec.num_lags as usize;
        let phases = acfd[start_idx..end_idx]
            .chunks_exact(2)
            .map(|x| (x[1] as f64).atan2(x[0] as f64))
            .collect();
        let t = lags
            .iter()
            .map(|x| (x.lag_num * rec.multi_pulse_increment as i32) as f64 * 1.0e-6)
            .collect();
        let std_dev = (0..rec.num_lags).map(|_| 0.0).collect();
        Ok(PhaseNode { phases, t, std_dev })
    }
    pub fn remove(&mut self, idx: usize) {
        self.phases.remove(idx);
        self.t.remove(idx);
        self.std_dev.remove(idx);
    }
}

#[derive(Debug)]
pub struct PowerNode {
    pub ln_power: Vec<f64>,
    pub t: Vec<f64>,
    pub std_dev: Vec<f64>,
}
impl PowerNode {
    pub fn new(
        rec: &RawacfRecord,
        lags: &Vec<LagNode>,
        range_idx: usize,
        range_num: usize,
        alpha_2: &Vec<f64>,
    ) -> PowerNode {
        let pwr_0 = rec.lag_zero_power.data[range_num] as f64;
        // acfs stores as [num_ranges, num_lags, 2] in memory, with 2 corresponding to real, imag
        let start_idx = range_idx * 2 * rec.num_lags as usize;
        let end_idx = start_idx + 2 * rec.num_lags as usize;
        let powers: Vec<f64> = rec.acfs.data[start_idx..end_idx]
            .chunks_exact(2)
            .map(|x| {
                let real = x[0] as f64;
                let imag = x[1] as f64;
                (real * real + imag * imag).sqrt()
            })
            .collect();
        let normalized_power: Vec<f64> = powers.iter().map(|x| x * x / (pwr_0 * pwr_0)).collect();

        let sigmas: Vec<f64> = zip(normalized_power.iter(), alpha_2.iter())
            .map(|(pwr_norm, alpha)| {
                pwr_0 * ((pwr_norm + 1.0 / alpha) / (2.0 * rec.num_averages as f64)).sqrt()
            })
            .collect();
        let t = lags
            .iter()
            .map(|x| (x.lag_num * rec.multi_pulse_increment as i32) as f64 * 1.0e-6)
            .collect();
        PowerNode {
            ln_power: powers.iter().map(|x| x.ln()).collect(),
            t,
            std_dev: sigmas,
        }
    }
    pub fn remove(&mut self, idx: usize) {
        self.ln_power.remove(idx);
        self.t.remove(idx);
        self.std_dev.remove(idx);
    }
}

#[derive(Debug)]
pub struct LagNode {
    pub lag_num: i32,
    pub pulses: [usize; 2],
    pub lag_idx: i32,
    pub sample_base_1: i32,
    pub sample_base_2: i32,
}

#[derive(Default, Debug)]
pub struct FittedData {
    pub delta: f64,
    pub intercept: f64,
    pub slope: f64,
    pub variance_intercept: f64,
    pub variance_slope: f64,
    pub delta_intercept: f64,
    pub delta_slope: f64,
    pub covariance_intercept_slope: f64,
    pub residual_intercept_slope: f64,
    pub quality: f64,
    pub chi_squared: f64,
}

#[derive(Default, Debug)]
pub struct Sums {
    pub sum: f64,
    pub sum_x: f64,
    pub sum_y: f64,
    pub sum_xx: f64,
    pub sum_xy: f64,
}

pub enum FitType {
    Linear,
    Quadratic,
}
