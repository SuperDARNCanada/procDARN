use dmap::error::DmapError;
use pyo3::exceptions::PyValueError;
use pyo3::PyErr;
use thiserror::Error;
use crate::error::ProcdarnError;

/// Enum of the possible error variants that may be encountered
#[derive(Error, Debug)]
pub enum FittingError {
    /// Represents an error in the Rawacf record that is attempting to be fitted
    #[error("{0}")]
    InvalidRawacf(String),

    /// Represents a bad fit of the record, for any reason
    #[error("{0}")]
    BadFit(String),

    /// Unable to get hardware file information
    #[error("{0}")]
    Hdw(#[from] ProcdarnError),

    /// Invalid DMAP file
    #[error("{0}")]
    Dmap(#[from] DmapError)
}

impl From<FittingError> for PyErr {
    fn from(value: FittingError) -> Self {
        let msg = value.to_string();
        PyValueError::new_err(msg)
    }
}