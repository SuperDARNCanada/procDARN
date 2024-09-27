//! Error type for Fitacfv3 algorithm
use crate::error::ProcdarnError;
use crate::fitting::fitacf3::determinations::determinations;
use crate::fitting::fitacf3::filtering;
use crate::fitting::fitacf3::fitstruct::{LagNode, RangeNode};
use crate::fitting::fitacf3::fitting;
use crate::utils::hdw::HdwInfo;
use crate::utils::rawacf::{get_hdw, Rawacf};
use dmap::formats::{fitacf::FitacfRecord, rawacf::RawacfRecord};
use pyo3::exceptions::PyValueError;
use pyo3::PyErr;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::f64::consts::PI;
use thiserror::Error;

type Result<T> = std::result::Result<T, Fitacf3Error>;

pub const FLUCTUATION_CUTOFF_COEFFICIENT: f32 = 2.0;
pub const ALPHA_CUTOFF: f32 = 2.0;
pub const ACF_SNR_CUTOFF: f64 = 1.0;
pub const MIN_LAGS: i16 = 3;

/// Enum of the possible error variants that may be encountered
#[derive(Error, Debug)]
pub enum Fitacf3Error {
    /// Represents an error in the Rawacf record that is attempting to be fitted
    #[error("{0}")]
    InvalidRawacf(String),

    /// Represents a bad fit of the record, for any reason
    #[error("{0}")]
    BadFit(String),

    /// Unable to get hardware file information
    #[error("{0}")]
    Hdw(#[from] ProcdarnError),
}

impl From<Fitacf3Error> for PyErr {
    fn from(value: Fitacf3Error) -> Self {
        let msg = value.to_string();
        PyValueError::new_err(msg)
    }
}

/// Fits a single `RawacfRecord` into a `FitacfRecord`
///
/// # Errors
/// Will return `Err` if the `RawacfRecord` does not have all required fields for fitting,
/// or if the data within the `RawacfRecord` is unsuitable for fitting for any reason.
fn fit_rawacf_record(record: &RawacfRecord, hdw: &HdwInfo) -> Result<FitacfRecord> {
    let raw: Rawacf = Rawacf::try_from(record).map_err(|e| {
        Fitacf3Error::InvalidRawacf(format!(
            "Could not extract all required fields from rawacf record: {e}"
        ))
    })?;
    let lags = create_lag_list(&raw);

    let noise_power = if raw.nave <= 0 {
        1.0
    } else {
        acf_cutoff_power(&raw)
    };
    let mut range_list = vec![];
    for i in 0..raw.slist.len() {
        let range_num = raw.slist[i];
        if raw.pwr0[range_num as usize] != 0.0 {
            range_list.push(RangeNode::new(i, range_num as usize, &raw, &lags)?);
        }
    }
    filtering::filter_tx_overlapped_lags(&raw, &lags, &mut range_list);
    filtering::filter_infinite_lags(&mut range_list);
    filtering::filter_low_power_lags(&raw, &mut range_list);
    filtering::filter_bad_acfs(&raw, &mut range_list, noise_power);
    fitting::acf_power_fitting(&mut range_list)?;
    fitting::calculate_phase_and_elev_sigmas(&mut range_list, &raw)?;
    fitting::acf_phase_unwrap(&mut range_list);
    fitting::acf_phase_fitting(&mut range_list)?;
    filtering::filter_bad_fits(&mut range_list)?;
    fitting::xcf_phase_unwrap(&mut range_list)?;
    fitting::xcf_phase_fitting(&mut range_list)?;

    determinations(&raw, &range_list, noise_power, hdw)
}

/// Fits a collection of `RawacfRecord`s into `FitacfRecord`s.
///
/// # Errors
/// Will return `Err` if the `RawacfRecord`s do not have all required fields for fitting,
/// or if the data within the `RawacfRecord`s are unsuitable for fitting for any reason.
pub fn fitacf3(raw_recs: Vec<RawacfRecord>) -> Result<Vec<FitacfRecord>> {
    let hdw = get_hdw(&raw_recs[0])?;

    let mut fitacf_records = vec![];
    for rec in raw_recs {
        fitacf_records.push(fit_rawacf_record(&rec, &hdw)?);
    }
    Ok(fitacf_records)
}

/// Fits a collection of `RawacfRecord`s into `FitacfRecord`s in parallel.
///
/// # Errors
/// Will return `Err` if the `RawacfRecord`s do not have all required fields for fitting,
/// or if the data within the `RawacfRecord`s are unsuitable for fitting for any reason.
pub fn par_fitacf3(raw_recs: Vec<RawacfRecord>) -> Result<Vec<FitacfRecord>> {
    let hdw = get_hdw(&raw_recs[0])?;

    // Fit the records!
    let fitacf_results: Vec<Result<FitacfRecord>> = raw_recs
        .par_iter()
        .map(|rec| fit_rawacf_record(rec, &hdw))
        .collect();

    let mut fitacf_records = vec![];
    for res in fitacf_results {
        match res {
            Ok(x) => fitacf_records.push(x),
            Err(e) => Err(e)?,
        }
    }
    Ok(fitacf_records)
}

/// Creates the lag table based on the data.
fn create_lag_list(record: &Rawacf) -> Vec<LagNode> {
    let lag_table = &record.ltab;
    let pulse_table = &record.ptab;
    let multi_pulse_increment = record.mpinc;
    let sample_separation = record.smsep;

    let mut lags = vec![];
    for i in 0..record.mplgs as usize {
        let mut pulse_1_idx = 0;
        let mut pulse_2_idx = 0;
        let number = lag_table[[i, 1]] - lag_table[[i, 0]];
        for j in 0..record.mppul as usize {
            if lag_table[[i, 0]] == pulse_table[j] {
                pulse_1_idx = j;
            }
            if lag_table[[i, 1]] == pulse_table[j] {
                pulse_2_idx = j;
            }
        }
        let sample_base_1 =
            i32::from(lag_table[[i, 0]] * (multi_pulse_increment / sample_separation));
        let sample_base_2 =
            i32::from(lag_table[[i, 1]] * (multi_pulse_increment / sample_separation));
        lags.push(LagNode {
            lag_num: i32::from(number),
            pulses: [pulse_1_idx, pulse_2_idx],
            sample_base_1,
            sample_base_2,
        });
    }
    lags
}

/// Calculates the minimum power value for ACFs in the record
fn acf_cutoff_power(rec: &Rawacf) -> f32 {
    let mut sorted_power_levels = rec.pwr0.clone().to_vec();
    sorted_power_levels.sort_by(f32::total_cmp); // sort floats
    let mut i: usize = 0;
    let mut j: f64 = 0.0;
    let mut min_power: f64 = 0.0;
    while j < 10.0 && i < rec.nrang as usize / 3 {
        if sorted_power_levels[i] > 0.0 {
            j += 1.0;
        }
        min_power += sorted_power_levels[i] as f64;
        i += 1;
    }
    if j <= 0.0 {
        j = 1.0;
    }
    min_power *= cutoff_power_correction(rec) / j;
    let search_noise = rec.noise_search;
    if min_power < ACF_SNR_CUTOFF && search_noise > 0.0 {
        min_power = search_noise as f64;
    }
    min_power as f32
}

/// Applies a correction to the noise power estimate to account for selecting least-powerful ranges
fn cutoff_power_correction(rec: &Rawacf) -> f64 {
    let std_dev = 1.0 / (rec.nave as f64).sqrt();

    let mut i = 0.0_f64;
    let mut cumulative_pdf = 0.0_f64;
    let mut cumulative_pdf_x_norm_power = 0.0_f64;
    let mut normalized_power: f64;
    while cumulative_pdf < (10.0 / rec.nrang as f64) {
        // Normalized power for calculating model PDF (Gaussian)
        normalized_power = i / 1000.0;
        let x = -(normalized_power - 1.0) * (normalized_power - 1.0) / (2.0 * std_dev * std_dev);
        let pdf = x.exp() / std_dev / (2.0 * PI).sqrt() / 1000.0;
        cumulative_pdf += pdf;

        // Cumulative value of PDF * x  -> needed for calculating the mean
        cumulative_pdf_x_norm_power += pdf * normalized_power;
        i += 1.0;
    }
    // Correcting factor as the inverse of a normalized mean
    cumulative_pdf / cumulative_pdf_x_norm_power
}
