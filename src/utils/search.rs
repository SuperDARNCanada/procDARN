use crate::error::BackscatterError;
use chrono::{NaiveDate, NaiveDateTime};
use dmap::formats::FitacfRecord;

/// Finds the first FitacfRecord in fitacf_records which occurs at or after date_time.
/// Called FitSeek/FitFSeek in RST
pub fn fit_seek(
    fitacf_records: &Vec<FitacfRecord>,
    date_time: NaiveDateTime,
) -> Option<(&FitacfRecord, usize)> {
    let record_times = fitacf_records
        .iter()
        .map(|rec| {
            NaiveDate::from_ymd_opt(rec.year as i32, rec.month as u32, rec.day as u32)?
                .and_hms_opt(rec.hour as u32, rec.minute as u32, rec.second as u32)?
        })
        .collect();

    match record_times.into_iter().position(|t| t >= date_time) {
        Some(i) => Some((fitacf_records[i], i)),
        None => None,
    }
}
