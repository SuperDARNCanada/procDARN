use crate::formats::{RawacfRecord, FitacfRecord};
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
fn create_lag_list(record: RawacfRecord) -> Result<Vec<[i32; 2]>> {
    let lag_table = record.lag_table;
    let pulse_table = record.pulse_table;
    let multi_pulse_increment = record.multi_pulse_increment;
    let sample_separation = record.sample_separation;

    let mut lags = vec![];
    for i in 0..record.num_lags as usize {
        let number = lag_table[2*i + 1] - lag_table[2*i];   // flattened, we want row i, cols 1 and 0
        for j in 0..record.num_pulses as usize {
            if lag_table[2*i] == pulse_table[j] {
                let pulse_1_idx = j;
            }
            if lag_table[2*i + 1] == pulse_table[j] {
                let pulse_2_idx = j;
            }
        }
        let sample_base_1 = lag_table[2*i] * (multi_pulse_increment / sample_separation);
        let sample_base_2 = lag_table[2*i + 1] * (multi_pulse_increment / sample_separation);
        let pulses = lag_table[i];
    }
    Ok(lags)
}
