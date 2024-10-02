use crate::fitting::common::error::FittingError;
use crate::utils::constants::US_TO_S;
use crate::utils::rawacf::Rawacf;
use numpy::ndarray::prelude::*;

#[derive(Debug)]
pub(crate) struct RangeNode {
    pub range_num: u16,
    pub t: Vec<f64>,
    pub lags: Vec<usize>,
    pub acf_real: Vec<f32>,
    pub acf_imag: Vec<f32>,
    pub sigma_real: Option<Vec<f64>>,
    pub sigma_imag: Option<Vec<f64>>,
    pub lin_fit: Option<FittedData>,
    pub quad_fit: Option<FittedData>,
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
        Ok(RangeNode {
            range_num: range_num as u16,
            t: lags
                .iter()
                .map(|x| (x.lag_num * record.mpinc as i32) as f64 * US_TO_S as f64)
                .collect(),
            lags: (0..lags.len()).collect(),
            acf_real: record
                .acfd
                .slice(s![index, .., 0])
                .iter()
                .copied()
                .collect(),
            acf_imag: record
                .acfd
                .slice(s![index, .., 1])
                .iter()
                .copied()
                .collect(),
            sigma_real: None,
            sigma_imag: None,
            lin_fit: None,
            quad_fit: None,
            elev_fit: None,
            self_clutter: None,
        })
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
    pub pwr: f64,
    pub wid: f64,
    pub vel: f64,
    pub phi: f64,
    pub sigma_2_pwr: f64,
    pub sigma_2_wid: f64,
    pub sigma_2_vel: f64,
    pub sigma_2_phi: f64,
    pub chi_squared: f64,
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
