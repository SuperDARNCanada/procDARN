use crate::fitting::fitstruct::{Alpha, LagNode, RangeNode};
use dmap::formats::{RawacfRecord, FitacfRecord};
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
fn create_lag_list(record: &RawacfRecord) -> Result<Vec<[i32; 2]>> {
    let lag_table = record.lag_table;
    let pulse_table = record.pulse_table;
    let multi_pulse_increment = record.multi_pulse_increment;
    let sample_separation = record.sample_separation;

    let mut lags = vec![];
    for i in 0..record.num_lags as usize {
        let number = lag_table.data[2*i + 1] - lag_table.data[2*i];   // flattened, we want row i, cols 1 and 0
        for j in 0..record.num_pulses as usize {
            if lag_table.data[2*i] == pulse_table.data[j] {
                let pulse_1_idx = j;
            }
            if lag_table.data[2*i + 1] == pulse_table.data[j] {
                let pulse_2_idx = j;
            }
        }
        let sample_base_1 = lag_table.data[2*i] * (multi_pulse_increment / sample_separation);
        let sample_base_2 = lag_table.data[2*i + 1] * (multi_pulse_increment / sample_separation);
        let pulses = lag_table.data[i];
    }
    Ok(lags)
}

fn fit_rawacf_record(record: &RawacfRecord) -> Result<FitacfRecord> {
    let lags = create_lag_list(record);

    if record.num_averages <= 0 {
        let noise_power = 1.0;
    } else {
        let noise_power = 0.0;  // = acf_cutoff_power(record)
    }

    let mut range_list = vec![];
    for i in 0..record.range_list.data.len() {
        let range_num = record.range_list.data[i];
        if record.lag_zero_power.data[range_num as usize] != 0.0 {
            range_list.push(RangeNode::new(i, range_num, record, lags))
        }
    }
    Err(Fitacf3Error::Message(format!("Unable to fit record")))
}

fn calculate_alpha_at_lags(lag: LagNode, range: RangeNode, lag_zero_power: f64) -> Alpha {
    let pulse_i_cri = range.cross_range_interference[lag.pulses[0] as usize];
    let pulse_j_cri = range.cross_range_interference[lag.pulses[1] as usize];

    let lag_idx = lag.lag_idx;
    let alpha_2 = lag_zero_power*lag_zero_power / ((lag_zero_power + pulse_i_cri) * (lag_zero_power + pulse_j_cri));
    // let range.alpha_2 = alpha_2;

}