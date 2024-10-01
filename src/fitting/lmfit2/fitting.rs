use crate::fitting::common::error::FittingError;
use crate::fitting::lmfit2::fitstruct::{FittedData, RangeNode};
use crate::utils::constants::{LIGHTSPEED_f64, US_TO_S_f64};
use crate::utils::rawacf::Rawacf;
use itertools::enumerate;
use rmpfit::{MPConfig, MPFitter, MPPar, MPResult};
use std::f64::consts::PI;

pub const NUM_VEL_MODELS: u32 = 30;
const CONFIDENCE: i32 = 1;

pub(crate) fn acf_fit(range_list: &mut Vec<RangeNode>, raw: &Rawacf) -> Result<(), FittingError> {
    for range in range_list {
        range.lin_fit = Some(lmfit(range, raw)?);
    }
    Ok(())
}

fn lmfit(range_node: &mut RangeNode, raw: &Rawacf) -> Result<FittedData, FittingError> {
    let wavelength: f64 = LIGHTSPEED_f64 / raw.tfreq as f64;
    let nyquist_vel: f64 = wavelength / (4.0 * raw.mpinc as f64 * US_TO_S_f64);
    let vel_step: f64 = (nyquist_vel + 1.0) / (NUM_VEL_MODELS as f64 - 1.0);
    let delta_chi: i32 = CONFIDENCE * CONFIDENCE;

    // independent variable for our data
    let t: Vec<f64> = [range_node.t.clone(), range_node.t.clone()].concat(); // repeat since data goes real then imaginary
    let real_acf: Vec<f64> = range_node.acf_real.iter().map(|&x| x as f64).collect();
    let imag_acf: Vec<f64> = range_node.acf_imag.iter().map(|&x| x as f64).collect();
    let y: Vec<f64> = [real_acf, imag_acf].concat();
    let ye: Vec<f64> = [
        range_node
            .sigma_real
            .as_ref()
            .ok_or_else(|| FittingError::BadFit("Cannot fit without error estimate".to_string()))?
            .clone(),
        range_node
            .sigma_imag
            .as_ref()
            .ok_or_else(|| FittingError::BadFit("Cannot fit without error estimate".to_string()))?
            .clone(),
    ]
    .concat();

    let mut problem = LevMarProblem::new(
        t,
        y,
        ye,
        wavelength,
        nyquist_vel,
    );

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
        let mut params =
            vec![10_000.0, 200.0, -nyquist_vel / 2.0 + i as f64 * vel_step];
        let result = problem
            .mpfit(&mut params)
            .map_err(|e| FittingError::BadFit(format!("Error with MPFit: {e}")))?;

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

    for i in 0..NUM_VEL_MODELS as usize {
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
    Ok(fit)
}

/// Levenberg-Marquardt solver using the rmpfit crate
#[derive(Default)]
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

    /// The actual parameters being optimized
    params: Vec<MPPar>
}

impl LevMarProblem {
    pub fn new(
        t: Vec<f64>,
        y: Vec<f64>,
        ye: Vec<f64>,
        wavelength: f64,
        nyquist_vel: f64,
    ) -> LevMarProblem {
        let mut params: Vec<MPPar> = vec![];

        let mut pwr_param = MPPar {
            limited_low: true,
            limit_low: 0.0,
            ..Default::default()
        };
        params.push(pwr_param);

        let mut wid_param = MPPar {
            limited_low: true,
            limit_low: -100.0,
            ..Default::default()
        };
        params.push(wid_param);

        let mut vel_param = MPPar {
            limited_low: true,
            limit_low: -nyquist_vel / 2.0,
            limited_up: true,
            limit_up: nyquist_vel / 2.0,
            ..Default::default()
        };
        params.push(vel_param);

        LevMarProblem {x: t, y, ye, wavelength, nyquist_vel, params}
    }
}
impl MPFitter for LevMarProblem {
    fn eval(&mut self, params: &[f64], deviates: &mut [f64]) -> MPResult<()> {
        let exponential: Vec<f64> = self
            .x
            .iter()
            .map(|x| (-2.0 * PI * params[1] * x / self.wavelength).exp())
            .collect();
        let coeff = 4.0 * PI * params[2] / self.wavelength;

        let num_points = deviates.len();
        for (i, dev) in enumerate(deviates.iter_mut()) {
            if i < num_points / 2 {
                *dev = (self.y[i] - params[0] * exponential[i] * (coeff * self.x[i]).cos())
                    / self.ye[i];
            } else {
                *dev = (self.y[i] - params[0] * exponential[i] * (coeff * self.x[i]).sin())
                    / self.ye[i]
            }
        }
        Ok(())
    }

    fn number_of_points(&self) -> usize {
        self.x.len()
    }

    fn config(&self) -> MPConfig {
        MPConfig {
            ftol: 0.0001,
            gtol: 0.0001,
            no_finite_check: false,
            max_fev: 200,
            ..Default::default()
        }
    }

    fn parameters(&self) -> Option<&[MPPar]> {
        Some(&*self.params)
    }
}
