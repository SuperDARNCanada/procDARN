use std::f64::consts::PI;
use itertools::enumerate;
use levenberg_marquardt::LeastSquaresProblem;
use nalgebra::{Matrix, Matrix3, Owned, U3, Vector, Vector3, VectorN};
use ndarray::{Array1, Array2, s};
use rmpfit::{MPFitter, MPResult};

pub const NUM_VEL_MODELS: u32 = 30;
const US_TO_S: f64 = 1e-6;

pub(crate) fn acf_fit(acf: Vec<f64>, wavelength: f64, mpinc: u16, confidence: i32, model: i32) {
    let nyquist_vel: f64 = wavelength / (4.0 * mpinc as f64 * US_TO_S);
    let vel_step: f64 = (nyquist_vel + 1.0) / (NUM_VEL_MODELS as f64 - 1.0);

    let delta_chi: i32 = confidence * confidence;

    for i in 0..NUM_VEL_MODELS {

    }
}

/// Levenberg-Marquardt solver using the varpro crate

/// Levenberg-Marquardt solver using the levenberg-marquardt crate
pub(crate) struct LevMarProblem {
    /// Independent variable of the problem (time)
    x: Vec<f64>,

    /// Real component of the ACF
    y: Vec<f64>,

    /// Uncertainty in real component of the ACF
    ye: Vec<f64>,

    /// Imaginary component of the ACF
    z: Vec<f64>,

    /// Uncertainty in real component of the ACF
    ze: Vec<f64>

    /// The current value of the 3 parameters being fitted
    pwr: Vec<f64>,
    wid: Vec<f64>,
    vel: Vec<f64>,


    /// The radio wavelength
    wavelength: f64,

    /// The time of the data (s)
    t: Array1<f64>,

    /// The real and imaginary components of the ACF
    acf: Array2<f64>,

    /// The uncertainty for the real and imaginary ACF components
    sigma: Array2<f64>,

    /// Exponential decay quantity from model
    exponential: Array1<f64>,

    /// Real component of oscillation from model
    cos: Array1<f64>,

    /// Imaginary component of oscillation from model
    sin: Array1<f64>,
}

impl MPFitter for LevMarProblem {

    fn set_params(&mut self, p: &Vector<f64, U3, Self::ParameterStorage>) {
        self.p.copy_from(p);

        let [pwr, wid, vel] = [p.x, p.y, p.z];

        // calculate these values for new parameters, for easily updating residuals and jacobian



    }

    fn params(&self) -> Vector<f64, U3, Self::ParameterStorage> {
        self.p
    }

    fn residuals(&self) -> Option<Vector<f64, U3, Self::ResidualStorage>> {
        let [pwr, wid, vel] = [self.p.x, self.p.y, self.p.z];

        let deviates_real = (&self.acf.slice(s![.., 0]) - pwr * &self.exponential * &self.cos) / &self.sigma.slice(s![.., 0]);
        let deviates_imag = (&self.acf.slice(s![.., 1]) - pwr * &self.exponential * &self.sin) / &self.sigma.slice(s![.., 1]);


        Some(Vector::new())
    }


    fn eval(&mut self, params: &[f64], deviates: &mut [f64]) -> MPResult<()> {
        let exponential = (-2.0 * PI * params[1] * &self.t / self.wavelength).exp();
        let cos = (4.0 * PI * params[2] * &self.t / self.wavelength).cos();
        let sin = (4.0 * PI * params[2] * &self.t / self.wavelength).sin();

        let f_real = params[0] * &exponential * cos;
        let f_imag = params[0] * exponential * sin;

        for (i, dev) in enumerate(deviates
            .iter_mut()) {
            *dev =
        }
            .zip(self.sigma.slice(s![.., 0]))
            .zip(self.sigma.slice(s![.., 1])) {

            *dev = (&self.acf.slice(s![.., 0]) - pwr * &self.exponential * &self.cos) / &self.sigma.slice(s![.., 0]);
        }
    }

    fn number_of_points(&self) -> usize {
        todo!()
    }
}
