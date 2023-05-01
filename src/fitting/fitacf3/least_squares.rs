use crate::fitting::fitacf3::fitstruct::FitType;

#[derive(Default)]
pub struct LeastSquaresValues {
    pub sum: f64,
    pub sum_x: f64,
    pub sum_y: f64,
    pub sum_xx: f64,
    pub sum_xy: f64,
    pub delta: f64,
    pub intercept: f64,
    pub slope: f64,
    pub sigma_squared_intercept: f64,
    pub sigma_squared_slope: f64,
    pub delta_intercept: f64,
    pub delta_slope: f64,
    pub covariance_slope_intercept: f64,
    pub residual_slope_intercept: f64,
    pub q_factor: f64,
    pub chi_squared: f64
}

pub struct LeastSquares {
    pub delta_chi_2: [[f64; 2]; 6],
    pub confidence: usize,
    pub degrees_of_freedom: usize
}
impl LeastSquares {
    pub fn new(confidence: usize, degrees_of_freedom: usize) -> LeastSquares {
        let delta_chi_2 = [
            [1.00, 2.30],
            [2.71, 4.61],
            [4.00, 6.17],
            [6.63, 9.21],
            [9.00, 11.8],
            [15.1, 18.4]];
        LeastSquares {
            delta_chi_2,
            confidence: confidence - 1,
            degrees_of_freedom: degrees_of_freedom - 1
        }
    }
    pub fn two_parameter_line_fit(&self, x_vals: Vec<f64>, y_vals: Vec<f64>, sigmas: Vec<f64>, fit_type: FitType) -> LeastSquaresValues {
        let mut lsq: LeastSquaresValues = Default::default();
        Self::find_sums(&mut lsq, x_vals, y_vals, sigmas, fit_type);

        lsq.delta = lsq.sum * lsq.sum_xx - lsq.sum_x * lsq.sum_x;
        lsq.intercept = (lsq.sum_xx * lsq.sum_y - lsq.sum_x * lsq.sum_xy) / lsq.delta;
        lsq.slope = (lsq.sum * lsq.sum_xy - lsq.sum_x * lsq.sum_y) / lsq.delta;
        lsq.sigma_squared_intercept = lsq.sum_xx / lsq.delta;
        lsq.sigma_squared_slope = lsq.sum / lsq.delta;
        lsq.covariance_slope_intercept = (-1.0 * lsq.sum_x) / lsq.delta;
        lsq.residual_slope_intercept = (-1.0 * lsq.sum_x) / (lsq.sum * lsq.sum_xx).sqrt();

        let delta_chi_2 = self.delta_chi_2[self.confidence][self.degrees_of_freedom];
        lsq.delta_intercept = delta_chi_2.sqrt() * lsq.sigma_squared_intercept.sqrt();
        lsq.delta_slope = delta_chi_2.sqrt() * lsq.sigma_squared_slope.sqrt();
        Self::calculate_chi_2(&mut lsq, x_vals, y_vals, sigmas, fit_type);
        lsq
    }
    pub fn one_parameter_line_fit(&self, x_vals: Vec<f64>, y_vals: Vec<f64>, sigmas: Vec<f64>) -> LeastSquaresValues {
        let mut lsq: LeastSquaresValues = Default::default();
        Self::find_sums(&mut lsq, x_vals, y_vals, sigmas, FitType::Linear);

        lsq.slope = lsq.sum_xy / lsq.sum_xx;
        lsq.sigma_squared_slope = 1.0 / lsq.sum_xx;

        let delta_chi_2 = self.delta_chi_2[self.confidence][self.degrees_of_freedom];
        lsq.delta_slope = delta_chi_2.sqrt() * lsq.sigma_squared_slope.sqrt();
        Self::calculate_chi_2(&mut lsq, x_vals, y_vals, sigmas, FitType::Linear);
        lsq
    }
    fn find_sums(least_squares: &mut LeastSquaresValues, x_vals: Vec<f64>, y_vals: Vec<f64>, sigmas: Vec<f64>, fit_type: FitType) {
        let nonzero_sigma: Vec<usize> = sigmas.iter().enumerate()
            .map(|(i, &x)| {
                if x != 0.0 { Some(i) }
                else { None }
            }).filter(|&x| x.is_some())
            .map(|x| x.unwrap())
            .collect();
        let sigma_squared: Vec<f64> = nonzero_sigma.iter().map(|&x| sigmas[x]*sigmas[x]).collect();

        let mut sum: f64 = sigma_squared.iter().map(|x| 1.0 / x).sum();
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut sum_xx: f64 = 0.0;
        let mut sum_xy: f64 = 0.0;

        match fit_type {
            FitType::Linear => {
                for (new, &orig) in nonzero_sigma.iter().enumerate() {
                    sum_x += x_vals[orig] / sigma_squared[new];
                    sum_y += y_vals[orig] / sigma_squared[new];
                    sum_xx += x_vals[orig]*x_vals[orig] / sigma_squared[new];
                    sum_xy += y_vals[orig]*y_vals[orig] / sigma_squared[new];
                }
            },
            FitType::Quadratic => {
                for (new, &orig) in nonzero_sigma.iter().enumerate() {
                    sum_x += x_vals[orig]*x_vals[orig] / sigma_squared[new];
                    sum_y += y_vals[orig] / sigma_squared[new];
                    sum_xx += x_vals[orig]*x_vals[orig]*x_vals[orig]*x_vals[orig] / sigma_squared[new];
                    sum_xy += x_vals[orig]*x_vals[orig]*y_vals[orig] / sigma_squared[new];
                }
            }
        }
        least_squares.sum = sum;
        least_squares.sum_x = sum_x;
        least_squares.sum_y = sum_y;
        least_squares.sum_xx = sum_xx;
        least_squares.sum_xy = sum_xy;
    }
    fn calculate_chi_2(lsq: &mut LeastSquaresValues, x_vals: Vec<f64>, y_vals: Vec<f64>, sigmas: Vec<f64>, fit_type: FitType) {
        let nonzero_sigma: Vec<usize> = sigmas.iter().enumerate()
            .map(|(i, &x)| {
                if x != 0.0 { Some(i) }
                else { None }
            }).filter(|&x| x.is_some())
            .map(|x| x.unwrap())
            .collect();
        let chi: Vec<f64> = vec![];
        match fit_type {
            FitType::Linear => {
                for &i in nonzero_sigma.iter() {
                    chi.push(((y_vals[i] - lsq.intercept) - (lsq.slope * x_vals[i])) / sigmas[i]);
                }
                lsq.chi_squared = chi.iter().map(|x| x*x).sum();
            }
            FitType::Quadratic => {
                for &i in nonzero_sigma.iter() {
                    chi.push(((y_vals[i] - lsq.intercept) - (lsq.slope * x_vals[i] * x_vals[i])) / sigmas[i]);
                }
                lsq.chi_squared = chi.iter().map(|x| x*x).sum();
            }
        }
    }
}