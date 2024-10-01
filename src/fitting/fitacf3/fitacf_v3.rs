//! Error type for Fitacfv3 algorithm
use crate::fitting::common::error::FittingError;
use crate::fitting::common::fitstruct::RangeNode;
use crate::fitting::common::preprocessing;
use crate::fitting::fitacf3::determinations::determinations;
use crate::fitting::fitacf3::{filtering, fitting};
use crate::utils::hdw::HdwInfo;
use crate::utils::rawacf::{get_hdw, Rawacf};
use dmap::formats::{fitacf::FitacfRecord, rawacf::RawacfRecord};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

type Result<T> = std::result::Result<T, FittingError>;

pub const FLUCTUATION_CUTOFF_COEFFICIENT: f32 = 2.0;
pub const ALPHA_CUTOFF: f32 = 2.0;
pub const MIN_LAGS: i16 = 3;

/// Fits a single `RawacfRecord` into a `FitacfRecord`
///
/// # Errors
/// Will return `Err` if the `RawacfRecord` does not have all required fields for fitting,
/// or if the data within the `RawacfRecord` is unsuitable for fitting for any reason.
fn fit_rawacf_record(record: &RawacfRecord, hdw: &HdwInfo) -> Result<FitacfRecord> {
    let raw: Rawacf = Rawacf::try_from(record).map_err(|e| {
        FittingError::InvalidRawacf(format!(
            "Could not extract all required fields from rawacf record: {e}"
        ))
    })?;
    let lags = preprocessing::create_lag_list(&raw);

    let noise_power = if raw.nave <= 0 {
        1.0
    } else {
        preprocessing::acf_cutoff_power(&raw)
    };
    let mut range_list = vec![];
    for i in 0..raw.slist.len() {
        let range_num = raw.slist[i];
        if raw.pwr0[range_num as usize] != 0.0 {
            range_list.push(RangeNode::new(i, range_num as usize, &raw, &lags)?);
        }
    }
    preprocessing::remove_tx_overlapped_lags(&raw, &lags, &mut range_list);
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
