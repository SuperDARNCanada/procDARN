use crate::fitting::fitacf3::fitstruct::{FittedData, PowerFitType, Sums};

#[derive(Debug)]
pub(crate) struct LeastSquares {
    pub delta_chi_2: [[f64; 2]; 6],
    pub confidence: usize,
    pub degrees_of_freedom: usize,
}
impl LeastSquares {
    pub(crate) fn new(confidence: usize, degrees_of_freedom: usize) -> LeastSquares {
        let delta_chi_2 = [
            [1.00, 2.30],
            [2.71, 4.61],
            [4.00, 6.17],
            [6.63, 9.21],
            [9.00, 11.8],
            [15.1, 18.4],
        ];
        LeastSquares {
            delta_chi_2,
            confidence: confidence - 1,
            degrees_of_freedom: degrees_of_freedom - 1,
        }
    }
    pub(crate) fn two_parameter_line_fit(
        &self,
        x_vals: &[f64],
        y_vals: &[f64],
        sigmas: &[f64],
        fit_type: &PowerFitType,
    ) -> FittedData {
        let mut fitted: FittedData = FittedData::default();
        let sums = Self::find_sums(x_vals, y_vals, sigmas, fit_type);
        fitted.delta = sums.sum * sums.sum_xx - sums.sum_x * sums.sum_x;
        fitted.intercept = (sums.sum_xx * sums.sum_y - sums.sum_x * sums.sum_xy) / fitted.delta;
        fitted.slope = (sums.sum * sums.sum_xy - sums.sum_x * sums.sum_y) / fitted.delta;
        fitted.variance_intercept = sums.sum_xx / fitted.delta;
        fitted.variance_slope = sums.sum / fitted.delta;
        fitted.covariance_intercept_slope = (-1.0 * sums.sum_x) / fitted.delta;
        fitted.residual_intercept_slope = (-1.0 * sums.sum_x) / (sums.sum * sums.sum_xx).sqrt();

        let delta_chi_2 = self.delta_chi_2[self.confidence][self.degrees_of_freedom];
        fitted.delta_intercept = delta_chi_2.sqrt() * fitted.variance_intercept.sqrt();
        fitted.delta_slope = delta_chi_2.sqrt() * fitted.variance_slope.sqrt();
        fitted.chi_squared = Self::calculate_chi_2(&fitted, x_vals, y_vals, sigmas, fit_type);
        fitted
    }
    pub(crate) fn one_parameter_line_fit(
        &self,
        x_vals: &[f64],
        y_vals: &[f64],
        sigmas: &[f64],
    ) -> FittedData {
        let mut fitted: FittedData = FittedData::default();
        let sums = Self::find_sums(x_vals, y_vals, sigmas, &PowerFitType::Linear);

        fitted.slope = sums.sum_xy / sums.sum_xx;
        fitted.variance_slope = 1.0 / sums.sum_xx;

        let delta_chi_2 = self.delta_chi_2[self.confidence][self.degrees_of_freedom];
        fitted.delta_slope = delta_chi_2.sqrt() * fitted.variance_slope.sqrt();
        fitted.delta_intercept = delta_chi_2.sqrt() * fitted.variance_intercept.sqrt();
        fitted.chi_squared =
            Self::calculate_chi_2(&fitted, x_vals, y_vals, sigmas, &PowerFitType::Linear);
        fitted
    }
    /// passing
    fn find_sums(x_vals: &[f64], y_vals: &[f64], sigmas: &[f64], fit_type: &PowerFitType) -> Sums {
        let nonzero_sigma: Vec<usize> = sigmas
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| if x == 0.0 { None } else { Some(i) })
            .collect();
        let sigma_squared: Vec<f64> = nonzero_sigma
            .iter()
            .map(|&x| sigmas[x] * sigmas[x])
            .collect();

        let sum: f64 = sigma_squared.iter().map(|x| 1.0 / x).sum();
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut sum_xx: f64 = 0.0;
        let mut sum_xy: f64 = 0.0;

        match fit_type {
            PowerFitType::Linear => {
                for (new, &orig) in nonzero_sigma.iter().enumerate() {
                    sum_x += x_vals[orig] / sigma_squared[new];
                    sum_y += y_vals[orig] / sigma_squared[new];
                    sum_xx += x_vals[orig] * x_vals[orig] / sigma_squared[new];
                    sum_xy += x_vals[orig] * y_vals[orig] / sigma_squared[new];
                }
            }
            PowerFitType::Quadratic => {
                for (new, &orig) in nonzero_sigma.iter().enumerate() {
                    sum_x += x_vals[orig] * x_vals[orig] / sigma_squared[new];
                    sum_y += y_vals[orig] / sigma_squared[new];
                    sum_xx += x_vals[orig] * x_vals[orig] * x_vals[orig] * x_vals[orig]
                        / sigma_squared[new];
                    sum_xy += x_vals[orig] * x_vals[orig] * y_vals[orig] / sigma_squared[new];
                }
            }
        }
        Sums {
            sum,
            sum_x,
            sum_y,
            sum_xx,
            sum_xy,
        }
    }
    fn calculate_chi_2(
        fitted: &FittedData,
        x_vals: &[f64],
        y_vals: &[f64],
        sigmas: &[f64],
        fit_type: &PowerFitType,
    ) -> f64 {
        let nonzero_sigma: Vec<usize> = sigmas
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| if x == 0.0 { None } else { Some(i) })
            .collect();
        let mut chi: Vec<f64> = vec![];
        match fit_type {
            PowerFitType::Linear => {
                for &i in &nonzero_sigma {
                    chi.push(
                        ((y_vals[i] - fitted.intercept) - (fitted.slope * x_vals[i])) / sigmas[i],
                    );
                }
                chi.iter().map(|x| x * x).sum()
            }
            PowerFitType::Quadratic => {
                for &i in &nonzero_sigma {
                    chi.push(
                        ((y_vals[i] - fitted.intercept) - (fitted.slope * x_vals[i] * x_vals[i]))
                            / sigmas[i],
                    );
                }
                chi.iter().map(|x| x * x).sum()
            }
        }
    }
}
