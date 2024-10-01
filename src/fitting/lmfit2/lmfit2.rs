//! Error type for Lmfitv2 algorithm
use crate::fitting::common::error::FittingError;
use crate::fitting::lmfit2::determinations::determinations;
use crate::fitting::lmfit2::estimations::{
    estimate_first_order_error, estimate_real_imag_error, estimate_self_clutter,
};
use crate::fitting::lmfit2::fitstruct::RangeNode;
use crate::fitting::lmfit2::fitting::acf_fit;
use crate::fitting::lmfit2::preprocessing;
use crate::utils::hdw::HdwInfo;
use crate::utils::rawacf::{get_hdw, Rawacf};
use dmap::formats::{fitacf::FitacfRecord, rawacf::RawacfRecord};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

type Result<T> = std::result::Result<T, FittingError>;

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

    //filtering::check_range_nodes(&mut range_list);
    estimate_self_clutter(&mut range_list, &raw);
    estimate_first_order_error(&mut range_list, &raw, noise_power as f64);
    acf_fit(&mut range_list, &raw)?;
    estimate_real_imag_error(&mut range_list, &raw, noise_power as f64)?;
    acf_fit(&mut range_list, &raw)?;
    // xcf_fit(&mut range_list, &raw);

    determinations(&raw, &range_list, noise_power, hdw)
}

/// Fits a collection of `RawacfRecord`s into `FitacfRecord`s.
///
/// # Errors
/// Will return `Err` if the `RawacfRecord`s do not have all required fields for fitting,
/// or if the data within the `RawacfRecord`s are unsuitable for fitting for any reason.
pub fn lmfit2(raw_recs: Vec<RawacfRecord>) -> Result<Vec<FitacfRecord>> {
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
pub fn par_lmfit2(raw_recs: Vec<RawacfRecord>) -> Result<Vec<FitacfRecord>> {
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
