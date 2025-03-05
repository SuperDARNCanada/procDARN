use crate::error::ProcdarnError;
use crate::utils::scan::RadarScan;
use dmap::error::DmapError;
use thiserror::Error;

type Result<T> = std::result::Result<T, GridError>;

/// Enum of the possible error variants that may be encountered
#[derive(Error, Debug)]
pub enum GridError {
    /// Represents an error in the Rawacf record that is attempting to be fitted
    #[error("{0}")]
    InvalidFitacf(String),

    /// Represents an error in processing of the record, for any reason
    #[error("{0}")]
    ProcessingError(String),

    /// Unable to get hardware file information
    #[error("{0}")]
    Hdw(#[from] ProcdarnError),

    /// Invalid DMAP file
    #[error("{0}")]
    Dmap(#[from] DmapError),

    /// Error in `igrf` crate
    #[error("{0}")]
    Igrf(#[from] igrf::Error),
}

/// Checks to make sure the radar operating parameters do not change significantly between scans.
/// If the frequency, distance to first range, or range separation change between scans, then the
/// scattering location for a range gate will also change, so median filtering the data is
/// nonsensical.
/// Called FilterCheckOps in checkops.c of RST.
pub fn check_operational_params(scans: &Vec<&RadarScan>, max_frequency_var: i32) -> bool {
    // Choose the middle scan of scans being median filtered
    let ref_scan = scans[scans.len() / 2];

    // Loop through other scans that are being median filtered
    for &scan in scans.iter().filter(|&s| *s != ref_scan) {
        // Loop through beams of the reference scan
        for ref_beam in ref_scan.beams.iter() {
            // Loop through beams of the scan under consideration
            for check_beam in scan.beams.iter().filter(|&b| b.beam == ref_beam.beam) {
                // Check if the relevant operating parameters are equal or close to equal
                if ref_beam.first_range != check_beam.first_range {
                    return false;
                }
                if ref_beam.range_sep != check_beam.range_sep {
                    return false;
                }
                if (ref_beam.freq - check_beam.freq).abs() > max_frequency_var {
                    return false;
                }
            }
        }
    }
    // If relevant operating parameters match or are close enough, then return true
    true
}
