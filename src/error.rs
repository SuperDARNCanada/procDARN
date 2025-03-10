use crate::utils::hdw::HdwError;
use dmap::error::DmapError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ProcdarnError {
    /// Represents a bad DMAP record
    #[error("{0}")]
    Dmap(#[from] DmapError),

    /// Unable to get hdw file information
    #[error("{0}")]
    Hdw(#[from] HdwError),

    /// Error in igrf crate usage
    #[error("{0}")]
    Igrf(#[from] igrf::Error),

    /// Error in geodesy crate usage
    #[error("{0}")]
    Geodesy(#[from] geodesy::Error),

    /// Invalid timestamp in a record
    #[error("invalid timestamp `{0}`")]
    Timestamp(String),

    /// Field missing from a record
    #[error("missing `{0}`")]
    MissingField(&'static str),

    /// Field from a record has the wrong type
    #[error("wrong type for `{0}`")]
    WrongType(&'static str),

    /// Zero records available
    #[error("{0}")]
    ZeroRecords(&'static str),

    /// Invalid channel specifier
    #[error("{0}")]
    Channel(String),
}
