use chrono::{NaiveDate, DateTime, Utc};
use dmap::formats::fitacf::FitacfRecord;

/// Finds the first FitacfRecord in fitacf_records which occurs at or after date_time.
/// Called FitSeek/FitFSeek in RST
pub fn fit_seek(
    fitacf_records: &Vec<FitacfRecord>,
    date_time: DateTime<Utc>,
) -> Option<(&FitacfRecord, usize)> {
    let record_times = fitacf_records
        .iter()
        .map(|rec| {
            NaiveDate::from_ymd_opt(
                rec.get(&"year".to_string())?.into(),
                rec.get(&"month".to_string())?.into(),
                rec.get(&"day".to_string())?.into()
            )?
                .and_hms_opt(
                    rec.get(&"hour".to_string())?.into(),
                    rec.get(&"minute".to_string())?.into(),
                    rec.get(&"second".to_string())?.into()
                )?
                .and_utc()
        })
        .collect();

    match record_times.into_iter().position(|t| t >= date_time) {
        Some(i) => Some((&fitacf_records[i], i)),
        None => None,
    }
}
