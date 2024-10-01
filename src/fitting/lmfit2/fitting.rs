use crate::fitting::common::error::FittingError;
use crate::fitting::lmfit2::fitstruct::{FittedData, RangeNode};
use crate::utils::constants::{LIGHTSPEED_f64, US_TO_S_f64};
use crate::utils::rawacf::Rawacf;
use itertools::enumerate;
use rmpfit::{MPConfig, MPFitter, MPPar, MPResult};
use std::f64::consts::PI;

pub const NUM_VEL_MODELS: u32 = 30;
const CONFIDENCE: i32 = 1;

pub(crate) fn acf_fit(range_list: &mut Vec<RangeNode>, raw: &Rawacf) {
    for range in range_list {
        range.lin_fit = Some(lmfit(range, raw)?);
    }
}

fn lmfit(range_node: &mut RangeNode, raw: &Rawacf) -> Result<FittedData, FittingError> {
    let wavelength: f64 = LIGHTSPEED_f64 / raw.tfreq as f64;
    let nyquist_vel: f64 = wavelength / (4.0 * raw.mpinc as f64 * US_TO_S_f64);
    let vel_step: f64 = (nyquist_vel + 1.0) / (NUM_VEL_MODELS as f64 - 1.0);
    let delta_chi: i32 = CONFIDENCE * CONFIDENCE;

    // independent variable for our data
    let t: Vec<f64> = range_node.t.clone().extend(range_node.t.clone()).collect(); // repeat since data goes real then imaginary
    let y: Vec<f64> = range_node
        .acf_real
        .clone()
        .extend(range_node.acf_imag.clone())
        .collect();
    let ye: Vec<f64> =
        range_node
            .sigma_real
            .as_ref()
            .ok_or_else(|| FittingError::BadFit("Cannot fit without error estimate".to_string()))?
            .extend(range_node.sigma_imag.as_ref().ok_or_else(|| {
                FittingError::BadFit("Cannot fit without error estimate".to_string())
            })?)
            .collect();

    let mut problem = LevMarProblem {
        x: t,
        y,
        ye,
        wavelength,
        nyquist_vel,
    };

    let mut fit: FittedData = FittedData::default();
    fit.chi_squared = 10e200; // arbitrary large number
    let mut chi_squared: Vec<f64> = vec![];
    let mut powers: Vec<f64> = vec![];
    let mut widths: Vec<f64> = vec![];
    let mut velocities: Vec<f64> = vec![];
    let mut power_err: Vec<f64> = vec![];
    let mut width_err: Vec<f64> = vec![];
    let mut velocity_err: Vec<f64> = vec![];

    for i in 0..NUM_VEL_MODELS {
        let mut params: &[f64] =
            &mut *vec![10_000.0, 200.0, -nyquist_vel / 2.0 + i as f64 * vel_step];
        let result = problem.mpfit(&mut params)?;

        chi_squared.push(result.best_norm);
        powers.push(params[0]);
        widths.push(params[1]);
        velocities.push(params[2]);
        power_err.push(result.xerror[0]);
        width_err.push(result.xerror[1]);
        velocity_err.push(result.xerror[2]);

        if result.best_norm < fit.chi_squared {
            fit.chi_squared = result.best_norm;
            fit.pwr = params[0];
            fit.wid = params[1];
            fit.vel = params[2];
            fit.sigma_2_pwr = CONFIDENCE as f64 * result.xerror[0];
            fit.sigma_2_wid = CONFIDENCE as f64 * result.xerror[1];
            fit.sigma_2_vel = CONFIDENCE as f64 * result.xerror[2];
        }
    }

    for i in 0..NUM_VEL_MODELS {
        if chi_squared[i] <= fit.chi_squared + delta_chi as f64 {
            if fit.sigma_2_pwr < (fit.pwr - powers[i]).abs() {
                fit.sigma_2_pwr = (fit.pwr - powers[i]).abs()
            }
            if fit.sigma_2_wid < (fit.wid - widths[i]).abs() {
                fit.sigma_2_wid = (fit.wid - widths[i]).abs()
            }
            if fit.sigma_2_vel < (fit.vel - velocities[i]).abs() {
                fit.sigma_2_vel = (fit.vel - velocities[i]).abs()
            }
        }
    }
    fit
}

/// Levenberg-Marquardt solver using the rmpfit crate
pub(crate) struct LevMarProblem {
    /// Independent variable of the ACF data
    x: Vec<f64>,

    /// flattened ACF, all real followed by all imaginary
    y: Vec<f64>,

    /// Uncertainty in the ACF components
    ye: Vec<f64>,

    /// The radio wavelength
    wavelength: f64,

    /// The upper limit on observable velocity given by the sampling rate
    nyquist_vel: f64,
}

impl MPFitter for LevMarProblem {
    fn eval(&mut self, params: &[f64], deviates: &mut [f64]) -> MPResult<()> {
        let exponential = (-2.0 * PI * params[1] * &self.x / self.wavelength).exp();
        let coeff = 4.0 * PI * params[2] / self.wavelength;

        for (i, dev) in enumerate(deviates.iter_mut()) {
            if i < deviates.len() / 2 {
                *dev = (self.y[i] - params[0] * exponential[i] * (coeff * self.x[i]).cos())
                    / self.ye[i];
            } else {
                *dev = (self.y[i] - params[0] * exponential[i] * (coeff * self.x[i]).sin())
                    / self.ye[i]
            }
        }
    }

    fn number_of_points(&self) -> usize {
        self.x.len()
    }

    fn config(&self) -> MPConfig {
        let mut default_config = MPConfig::default();
        default_config.ftol = 0.0001;
        default_config.gtol = 0.0001;
        default_config.no_finite_check = false;
        default_config.max_fev = 200;

        default_config
    }

    fn parameters(&self) -> Option<&[MPPar]> {
        let mut pwr_param = MPPar::default();
        pwr_param.limited_low = true;
        pwr_param.limit_low = 0.0;

        let mut wid_param = MPPar::default();
        wid_param.limited_low = true;
        wid_param.limit_low = -100.0;

        let mut vel_param = MPPar::default();
        vel_param.limited_low = true;
        vel_param.limit_low = -self.nyquist_vel / 2.0;
        vel_param.limited_up = true;
        vel_param.limit_up = self.nyquist_vel / 2.0;

        Some(&[pwr_param, wid_param, vel_param])
    }
}
