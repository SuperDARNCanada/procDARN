use crate::error::ProcdarnError;
use chrono::{DateTime, NaiveDate, Utc};
use dmap::formats::fitacf::FitacfRecord;

/// Gets the timestamp from a `FitacfRecord`
pub fn get_datetime(rec: &FitacfRecord) -> Result<DateTime<Utc>, ProcdarnError> {
    let date = NaiveDate::from_ymd_opt(
        rec.get(&"year".to_string())
            .ok_or(ProcdarnError::Timestamp("missing `year`"))?
            .clone()
            .try_into()?,
        rec.get(&"month".to_string())
            .ok_or(ProcdarnError::Timestamp("missing `month`"))?
            .clone()
            .try_into()?,
        rec.get(&"day".to_string())
            .ok_or(ProcdarnError::Timestamp("missing `day`"))?
            .clone()
            .try_into()?,
    )
    .ok_or(ProcdarnError::Timestamp("invalid ymd timestamp in record"))?;
    let dt = date
        .and_hms_micro_opt(
            rec.get(&"hour".to_string())
                .ok_or(ProcdarnError::Timestamp("missing `hour`"))?
                .clone()
                .try_into()?,
            rec.get(&"minute".to_string())
                .ok_or(ProcdarnError::Timestamp("missing `minute`"))?
                .clone()
                .try_into()?,
            rec.get(&"second".to_string())
                .ok_or(ProcdarnError::Timestamp("missing `second`"))?
                .clone()
                .try_into()?,
            rec.get(&"microsecond".to_string())
                .ok_or(ProcdarnError::Timestamp("missing `microsecond`"))?
                .clone()
                .try_into()?,
        )
        .ok_or(ProcdarnError::Timestamp(
            "invalid hms_micro timestamp in record",
        ))?
        .and_utc();

    Ok(dt)
}
