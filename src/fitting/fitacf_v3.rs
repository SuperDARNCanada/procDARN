use crate::dmap::DmapData;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fmt::Display;

type Result<T> = std::result::Result<T, Fitacf3Error>;

#[derive(Debug, Clone)]
pub enum Fitacf3Error {
    // Parse(String, Vec<u8>),
    // BadVal(String, DmapType),
    Message(String),
    Lookup(String),
    Mismatch { msg: String },
    // CastError(String, PodCastError),
}

impl Error for Fitacf3Error {}

impl Display for Fitacf3Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Fitacf3Error::Message(msg) => write!(f, "{}", msg),
            Fitacf3Error::Lookup(msg) => write!(f, "{}", msg),
            Fitacf3Error::Mismatch{msg} => write!(f, "{}", msg)

            // DmapError::BadVal(msg, val) => write!(f, "{}: {:?}", msg, val),
            // DmapError::Parse(msg, val) => write!(f, "{}: {:?}", msg, val),
            // DmapError::CastError(msg, err) => write!(f, "{}: {}", msg, err.to_string()),
        }
    }
}

/// Creates the lag table based on the data.
fn create_lag_list(record_hash: HashMap<String, DmapData>) -> Result<Vec<[i32; 2]>> {
    let num_lags = match record_hash.get("mplgs") {
        Some(x) => match x {
            DmapData::Scalar(y) => Ok(y),
            DmapData::Array(..) => Err(Fitacf3Error::Mismatch {
                msg: "Scalar type expected for mplgs, got array instead".to_string(),
            }),
        },
        None => Err(Fitacf3Error::Lookup(
            "mplgs not found in record".to_string(),
        )),
    }?;
    let pulse_table = match record_hash.get("ptab") {
        Some(x) => match x {
            DmapData::Array(y) => Ok(y),
            DmapData::Scalar(..) => Err(Fitacf3Error::Mismatch {
                msg: "Array type expected for ptab, got scalar instead".to_string(),
            }),
        },
        None => Err(Fitacf3Error::Lookup("ptab not found in record".to_string())),
    }?;
    let num_pulses = match record_hash.get("mppul") {
        Some(x) => match x {
            DmapData::Scalar(y) => Ok(y),
            DmapData::Array(..) => Err(Fitacf3Error::Mismatch {
                msg: "Scalar type expected for mppul, got array instead".to_string(),
            }),
        },
        None => Err(Fitacf3Error::Lookup(
            "mppul not found in record".to_string(),
        )),
    }?;
    let lag_table = match record_hash.get("ltab") {
        Some(x) => match x {
            DmapData::Scalar(y) => Ok(y),
            DmapData::Array(..) => Err(Fitacf3Error::Mismatch {
                msg: "Array type expected for ltab, got scalar instead".to_string(),
            }),
        },
        None => Err(Fitacf3Error::Lookup("ltab not found in record".to_string())),
    }?;
    let multi_pulse_increment = match record_hash.get("mpinc") {
        Some(x) => match x {
            DmapData::Scalar(y) => Ok(y),
            DmapData::Array(..) => Err(Fitacf3Error::Mismatch {
                msg: "Scalar type expected for mpinc, got array instead".to_string(),
            }),
        },
        None => Err(Fitacf3Error::Lookup(
            "mpinc not found in record".to_string(),
        )),
    }?;
    let sample_separation = match record_hash.get("smsep") {
        Some(x) => match x {
            DmapData::Scalar(y) => Ok(y),
            DmapData::Array(..) => Err(Fitacf3Error::Mismatch {
                msg: "Scalar type expected for smsep, got array instead".to_string(),
            }),
        },
        None => Err(Fitacf3Error::Lookup(
            "smsep not found in record".to_string(),
        )),
    }?;

    let mut lags = vec![];
    for i in 0..num_lags.data {
        let number = lag_table[i][1] - lag_table[i][0];
    }
    Ok(vec![])
}
