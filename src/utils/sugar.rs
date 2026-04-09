use crate::error::ProcdarnError;
use chrono::{DateTime, NaiveDate, Utc};
use dmap::formats::fitacf::FitacfRecord;
use dmap::Record;

/// Gets the timestamp from a `FitacfRecord`
pub fn get_datetime(rec: &FitacfRecord) -> Result<DateTime<Utc>, ProcdarnError> {
    let date = NaiveDate::from_ymd_opt(
        rec.get("time.yr")
            .ok_or(ProcdarnError::Timestamp("missing `time.yr`".to_string()))?
            .clone()
            .try_into()?,
        rec.get("time.mo")
            .ok_or(ProcdarnError::Timestamp("missing `time.mo`".to_string()))?
            .clone()
            .try_into()?,
        rec.get("time.dy")
            .ok_or(ProcdarnError::Timestamp("missing `time.dy`".to_string()))?
            .clone()
            .try_into()?,
    )
    .ok_or(ProcdarnError::Timestamp(
        "invalid ymd timestamp in record".to_string(),
    ))?;
    let dt = date
        .and_hms_micro_opt(
            rec.get("time.hr")
                .ok_or(ProcdarnError::Timestamp("missing `time.hr`".to_string()))?
                .clone()
                .try_into()?,
            rec.get("time.mt")
                .ok_or(ProcdarnError::Timestamp("missing `time.mt`".to_string()))?
                .clone()
                .try_into()?,
            rec.get("time.sc")
                .ok_or(ProcdarnError::Timestamp("missing `time.sc`".to_string()))?
                .clone()
                .try_into()?,
            rec.get("time.us")
                .ok_or(ProcdarnError::Timestamp("missing `time.us`".to_string()))?
                .clone()
                .try_into()?,
        )
        .ok_or(ProcdarnError::Timestamp(
            "invalid hms_micro timestamp in record".to_string(),
        ))?
        .and_utc();

    Ok(dt)
}
