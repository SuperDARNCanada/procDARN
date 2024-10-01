use crate::fitting::fitacf3::fitacf_v3::par_fitacf3;
use crate::fitting::common::error::FittingError;
use clap::Parser;
use dmap::error::DmapError;
use dmap::formats::dmap::Record;
use dmap::formats::rawacf::RawacfRecord;
use dmap::types::DmapField;
use indexmap::IndexMap;
use itertools::{Either, Itertools};
use pyo3::prelude::{PyAnyMethods, PyModule};
use pyo3::{pyfunction, pymodule, wrap_pyfunction, Bound, PyErr, PyResult, Python};
use std::path::PathBuf;

pub mod error;
pub mod fitting;
pub mod utils;

/// Fits a list of RAWACF records into FITACF records using the FITACFv3 algorithm.
#[pyfunction]
#[pyo3(name = "fitacf3")]
#[pyo3(text_signature = "(recs: list[dict], /)")]
fn fitacf3_py(
    mut recs: Vec<IndexMap<String, DmapField>>,
) -> PyResult<Vec<IndexMap<String, DmapField>>> {
    let (errors, formatted_recs): (Vec<_>, Vec<_>) =
        recs.iter_mut()
            .enumerate()
            .partition_map(|(i, rec)| match RawacfRecord::try_from(rec) {
                Err(e) => Either::Left((i, e)),
                Ok(x) => Either::Right(x),
            });
    if !errors.is_empty() {
        Err(PyErr::from(DmapError::InvalidRecord(format!(
            "Corrupted records: {errors:?}"
        ))))?
    }
    let fitacf_recs = par_fitacf3(formatted_recs)
        .map_err(PyErr::from)?
        .into_iter()
        .map(|rec| rec.inner())
        .collect();
    Ok(fitacf_recs)
}

/// Fits a RAWACF file into a FITACF record using the FITACFv3 algorithm.
fn file_fitacf3(raw_file: PathBuf, fit_file: PathBuf) -> Result<(), FittingError> {
    let rawacf_records = dmap::read_rawacf(raw_file)?;
    let fitacf_records = par_fitacf3(rawacf_records)?;
    dmap::write_fitacf(fitacf_records, &fit_file)?;
    Ok(())
}

/// Fits a RAWACF file into a FITACF record using the FITACFv3 algorithm.
#[pyfunction]
#[pyo3(name = "file_fitacf3")]
#[pyo3(text_signature = "(rawacf_file: str, fitacf_file: str, /)")]
fn file_fitacf3_py(raw_file: PathBuf, fit_file: PathBuf) -> PyResult<()> {
    file_fitacf3(raw_file, fit_file)?;
    Ok(())
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Fitacf3Args {
    /// Rawacf file to fit
    #[arg()]
    infile: PathBuf,

    /// Output fitacf file path
    #[arg()]
    outfile: PathBuf,
}

/// Fits a RAWACF file into a FITACF file using the FITACFv3 algorithm.
#[pyfunction]
#[pyo3(name = "fit_fitacf3")]
fn fitacf3_cli(py: Python) -> PyResult<()> {
    let argv = py
        .import_bound("sys")?
        .getattr("argv")?
        .extract::<Vec<String>>()?;
    let args = Fitacf3Args::parse_from(argv);

    let rawacf_records = dmap::read_rawacf(args.infile)?;
    let fitacf_records = par_fitacf3(rawacf_records)?;

    // Write to file
    dmap::write_fitacf(fitacf_records, &args.outfile)?;
    Ok(())
}

/// Functions for SuperDARN data processing.
#[pymodule]
fn procdarn(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fitacf3_py, m)?)?;
    m.add_function(wrap_pyfunction!(file_fitacf3_py, m)?)?;
    m.add_wrapped(wrap_pyfunction!(fitacf3_cli))?;

    Ok(())
}
