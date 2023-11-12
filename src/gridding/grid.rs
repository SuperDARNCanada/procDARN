use crate::utils::hdw::HdwInfo;
use dmap::formats::{FitacfRecord, GridRecord};
use std::error::Error;
use std::f64::consts::PI;
use std::fmt;
use std::fmt::Display;

type Result<T> = std::result::Result<T, GridError>;

#[derive(Debug, Clone)]
pub enum GridError {
    Message(String),
    Lookup(String),
    Mismatch { msg: String },
}

impl Error for GridError {}

impl Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            GridError::Message(msg) => write!(f, "{}", msg),
            GridError::Lookup(msg) => write!(f, "{}", msg),
            GridError::Mismatch { msg } => write!(f, "{}", msg),
        }
    }
}

// pub fn grid_fitacf_record(record: &FitacfRecord, hdw: &HdwInfo) -> Result<GridRecord> {
//
// }
