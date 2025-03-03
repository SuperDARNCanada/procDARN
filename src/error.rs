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
}
