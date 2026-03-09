use crate::error::ProcdarnError;
use chrono::{DateTime, NaiveDate, Utc};
use dmap::formats::fitacf::FitacfRecord;

/// Finds the first FitacfRecord in fitacf_records which occurs at or after date_time.
/// Called FitSeek/FitFSeek in RST
pub fn fit_seek(
    fitacf_records: &[FitacfRecord],
    date_time: DateTime<Utc>,
) -> Result<Option<(&FitacfRecord, usize)>, ProcdarnError> {
    let mut record_times: Vec<DateTime<Utc>> = vec![];
    for rec in fitacf_records.iter() {
        let tstamp = NaiveDate::from_ymd_opt(
            i32::try_from(
                rec.get(&"time.yr".to_string())
                    .ok_or(ProcdarnError::MissingField("time.yr"))?
                    .clone(),
            )?,
            u32::try_from(
                rec.get(&"time.mo".to_string())
                    .ok_or(ProcdarnError::MissingField("time.mo"))?
                    .clone(),
            )?,
            u32::try_from(
                rec.get(&"time.dy".to_string())
                    .ok_or(ProcdarnError::MissingField("time.dy"))?
                    .clone(),
            )?,
        )
        .ok_or(ProcdarnError::Timestamp("could not parse date".to_string()))?
        .and_hms_opt(
            u32::try_from(
                rec.get(&"time.hr".to_string())
                    .ok_or(ProcdarnError::MissingField("time.hr"))?
                    .clone(),
            )?,
            u32::try_from(
                rec.get(&"time.mt".to_string())
                    .ok_or(ProcdarnError::MissingField("time.mt"))?
                    .clone(),
            )?,
            u32::try_from(
                rec.get(&"time.sc".to_string())
                    .ok_or(ProcdarnError::MissingField("time.sc"))?
                    .clone(),
            )?,
        )
        .ok_or(ProcdarnError::Timestamp("could not parse time".to_string()))?
        .and_utc();
        record_times.push(tstamp);
    }

    match record_times.into_iter().position(|t| t >= date_time) {
        Some(i) => Ok(Some((&fitacf_records[i], i))),
        None => Ok(None),
    }
}
