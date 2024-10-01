use crate::fitting::common::error::FittingError;
use crate::utils::rawacf::Rawacf;
use numpy::ndarray::prelude::*;
use std::iter::zip;

#[derive(Debug)]
pub(crate) struct RangeNode {
    pub range_num: u16,
    pub range_idx: usize,
    // pub cross_range_interference: Vec<f64>,
    // pub refractive_idx: f32,
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
    pub self_clutter: Option<Vec<f64>>,
}
impl RangeNode {
    pub(crate) fn new(
        index: usize,
        range_num: usize,
        record: &Rawacf,
        lags: &[LagNode],
    ) -> Result<RangeNode, FittingError> {
        let cross_range_interference =
            RangeNode::calculate_cross_range_interference(range_num, record);
        let alpha_2 =
            RangeNode::calculate_alphas(range_num, &cross_range_interference, record, lags);
        let phases = PhaseNode::new(record, &PhaseFitType::Acf, lags, index)?;
        let elevations = PhaseNode::new(record, &PhaseFitType::Xcf, lags, index)?;
        let powers = PowerNode::new(record, lags, index, range_num, &alpha_2);
        Ok(RangeNode {
            range_idx: index,
            range_num: range_num as u16,
            // cross_range_interference,
            // refractive_idx: 1.0,
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
            self_clutter: None,
        })
    }
    fn calculate_cross_range_interference(range_num: usize, rec: &Rawacf) -> Vec<f64> {
        let tau: i16 = if rec.smsep != 0 {
            rec.mpinc / rec.smsep
        } else {
            // TODO: Log warning?
            rec.mpinc / rec.txpl
        };

        let mut interference_for_pulses: Vec<f64> = vec![];
        for pulse_to_check in 0..rec.mppul as usize {
            let mut total_interference: f64 = 0.0;
            for pulse in 0..rec.mppul as usize {
                let pulse_diff = rec.ptab[pulse_to_check] - rec.ptab[pulse];
                let range_to_check = (pulse_diff * tau + range_num as i16) as usize;
                if (pulse != pulse_to_check) && (range_to_check < rec.nrang as usize) {
                    total_interference += rec.pwr0[range_to_check] as f64;
                }
            }
            interference_for_pulses.push(total_interference);
        }
        interference_for_pulses
    }
    fn calculate_alphas(
        range_num: usize,
        cross_range_interference: &[f64],
        rec: &Rawacf,
        lags: &[LagNode],
    ) -> Vec<f64> {
        let mut alpha_2: Vec<f64> = vec![];
        for lag in lags {
            let pulse_1_interference = cross_range_interference[lag.pulses[0]];
            let pulse_2_interference = cross_range_interference[lag.pulses[1]];
            let lag_zero_power = rec.pwr0[range_num] as f64;
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
pub(crate) struct PhaseNode {
    pub phases: Vec<f64>,
    pub t: Vec<f64>,
    pub std_dev: Vec<f64>,
    pub std_dev_real: Vec<f64>,
    pub std_dev_imag: Vec<f64>,
}
impl PhaseNode {
    pub(crate) fn new(
        rec: &Rawacf,
        phase_type: &PhaseFitType,
        lags: &[LagNode],
        range_idx: usize,
    ) -> Result<PhaseNode, FittingError> {
        let acfd = match phase_type {
            PhaseFitType::Acf => &rec.acfd,
            PhaseFitType::Xcf => match &rec.xcfd {
                Some(ref x) => x,
                None => Err(FittingError::InvalidRawacf(
                    "Cannot find xcfs in data".to_string(),
                ))?,
            },
        };
        let phases = zip(
            acfd.slice(s![range_idx, .., 0]),
            acfd.slice(s![range_idx, .., 1]),
        )
        .map(|(&x, &y)| {
            let real = x as f64;
            let imag = y as f64;
            imag.atan2(real)
        })
        .collect();
        let t = lags
            .iter()
            .map(|x| (x.lag_num * rec.mpinc as i32) as f64 * 1.0e-6)
            .collect();
        let std_dev: Vec<f64> = (0..rec.mplgs).map(|_| 0.0).collect();
        let std_dev_real = std_dev.clone();
        let std_dev_imag = std_dev.clone();
        Ok(PhaseNode {
            phases,
            t,
            std_dev,
            std_dev_real,
            std_dev_imag,
        })
    }
    pub fn remove(&mut self, idx: usize) {
        self.phases.remove(idx);
        self.t.remove(idx);
        self.std_dev.remove(idx);
        self.std_dev_real.remove(idx);
        self.std_dev_imag.remove(idx);
    }
}

#[derive(Debug)]
pub(crate) struct PowerNode {
    pub ln_power: Vec<f64>,
    pub t: Vec<f64>,
    pub std_dev: Vec<f64>,
}
impl PowerNode {
    pub(crate) fn new(
        rec: &Rawacf,
        lags: &[LagNode],
        range_idx: usize,
        range_num: usize,
        alpha_2: &[f64],
    ) -> PowerNode {
        let pwr_0 = rec.pwr0[range_num] as f64;
        // acfs stores as [num_ranges, num_lags, 2] in memory, with 2 corresponding to real, imag
        let powers: Vec<f64> = zip(
            rec.acfd.slice(s![range_idx, .., 0]),
            rec.acfd.slice(s![range_idx, .., 1]),
        )
        .map(|(&x, &y)| {
            let real = x as f64;
            let imag = y as f64;
            (real * real + imag * imag).sqrt()
        })
        .collect();
        let normalized_power: Vec<f64> = powers.iter().map(|x| x * x / (pwr_0 * pwr_0)).collect();

        let sigmas: Vec<f64> = zip(normalized_power.iter(), alpha_2.iter())
            .map(|(pwr_norm, alpha)| {
                pwr_0 * ((pwr_norm + 1.0 / alpha) / (2.0 * rec.nave as f64)).sqrt()
            })
            .collect();
        let t = lags
            .iter()
            .map(|x| (x.lag_num * rec.mpinc as i32) as f64 * 1.0e-6)
            .collect();
        PowerNode {
            ln_power: powers.iter().map(|x| x.ln()).collect(),
            t,
            std_dev: sigmas,
        }
    }
    pub(crate) fn remove(&mut self, idx: usize) {
        self.ln_power.remove(idx);
        self.t.remove(idx);
        self.std_dev.remove(idx);
    }
}

#[derive(Debug)]
pub(crate) struct LagNode {
    pub lag_num: i32,
    pub pulses: [usize; 2],
    pub sample_base_1: i32,
    pub sample_base_2: i32,
}

#[derive(Default, Debug)]
pub(crate) struct FittedData {
    pub delta: f64,
    pub intercept: f64,
    pub slope: f64,
    pub variance_intercept: f64,
    pub variance_slope: f64,
    pub delta_intercept: f64,
    pub delta_slope: f64,
    pub covariance_intercept_slope: f64,
    pub residual_intercept_slope: f64,
    pub chi_squared: f64,
}

#[derive(Default, Debug)]
pub(crate) struct Sums {
    pub sum: f64,
    pub sum_x: f64,
    pub sum_y: f64,
    pub sum_xx: f64,
    pub sum_xy: f64,
}

#[derive(Copy, Clone)]
pub(crate) enum PowerFitType {
    Linear,
    Quadratic,
}

#[derive(Copy, Clone)]
pub(crate) enum PhaseFitType {
    Acf,
    Xcf,
}
