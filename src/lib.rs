//! Core SuperDARN Processing Tools
//!
//! [![github]](https://github.com/SuperDARNCanada/procdarn)
//!
//! [github]: https://img.shields.io/badge/github-8da0cb?style=for-the-badge&labelColor=555555&logo=github
//!
//! <br>
//!
//! This library also has a Python API using pyo3.
//!
//! This library is a re-implementation of the core tools from [SuperDARN's Radar Software Toolkit
//! (RST)](https://github.com/SuperDARN/rst). Currently, only two binaries from RST are implemented.
//! The goal is for the entire RAWACF -> FITACF -> GRID -> MAP
//! pipeline to be implemented.
//!
//! | RST binaries        | `procdarn` function |
//! | ------------------- | ------------------- |
//! | `make_fit -fitacf3` | [`fitacf3`]         |
//! | `make_grid`         | [`fit2grid`]        |
//!
//! The `procdarn` algorithms can be called with directly with data, or can be used to read data
//! from files and run it through the algorithms.

use clap::Parser;
use dmap::error::DmapError;
use dmap::formats::rawacf::RawacfRecord;
use dmap::record::Record;
use dmap::types::DmapField;
use indexmap::IndexMap;
use itertools::{Either, Itertools};
use pyo3::prelude::{PyAnyMethods, PyModule, PyModuleMethods};
use pyo3::{pyfunction, pymodule, wrap_pyfunction, Bound, PyErr, PyResult, Python};
use std::path::PathBuf;
use dmap::formats::fitacf::FitacfRecord;
use dmap::GridRecord;
use pyo3::types::PyDict;

pub mod error;
pub mod fitting;
pub mod gridding;
mod utils;

pub use crate::fitting::fitacf3::fitacf_v3::{fitacf3, Fitacf3Error};
pub use crate::gridding::grid::{fit2grid, fit2grid_file, GridArgs, GridError};

/// Fits a list of RAWACF records into FITACF records using the FITACFv3 algorithm.
#[pyfunction]
#[pyo3(name = "fitacf3_recs")]
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
    let fitacf_recs = fitacf3(formatted_recs)
        .map_err(PyErr::from)?
        .into_iter()
        .map(|rec| rec.inner())
        .collect();
    Ok(fitacf_recs)
}

/// Fits a RAWACF file into a FITACF record using the FITACFv3 algorithm.
fn fitacf3_file(raw_file: PathBuf, fit_file: PathBuf) -> Result<(), Fitacf3Error> {
    let rawacf_records = RawacfRecord::read_file(raw_file)?;
    let fitacf_records = fitacf3(rawacf_records)?;
    FitacfRecord::write_to_file(&fitacf_records, &fit_file, false)?;
    Ok(())
}

/// Fits a RAWACF file into a FITACF record using the FITACFv3 algorithm.
#[pyfunction]
#[pyo3(name = "fitacf3_file")]
#[pyo3(text_signature = "(rawacf_file: str, fitacf_file: str, /)")]
fn fitacf3_file_py(raw_file: PathBuf, fit_file: PathBuf) -> PyResult<()> {
    fitacf3_file(raw_file, fit_file)?;
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
#[pyo3(name = "raw2fit_cli")]
fn fitacf3_cli(py: Python) -> PyResult<()> {
    let argv = py
        .import("sys")?
        .getattr("argv")?
        .extract::<Vec<String>>()?;
    let args = Fitacf3Args::parse_from(argv);

    let rawacf_records = RawacfRecord::read_file(args.infile)?;
    let fitacf_records = fitacf3(rawacf_records)?;

    // Write to file
    FitacfRecord::write_to_file(&fitacf_records, &args.outfile, false)?;
    Ok(())
}


/// Converts a list of FITACF records into GRID records.
#[pyfunction]
#[pyo3(name = "fit2grid_recs")]
#[pyo3(signature = (recs, /, **py_kwargs))]
#[pyo3(text_signature = "(recs: list[dict], /, **)")]
fn fit2grid_py(
    mut recs: Vec<IndexMap<String, DmapField>>,
    py_kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<IndexMap<String, DmapField>>> {
    let args = match py_kwargs {
        Some(kwargs) => kwargs.extract()?,
        None => GridArgs::parse_from(vec!["fit2grid"]),
    };

    let (errors, formatted_recs): (Vec<_>, Vec<_>) =
        recs.iter_mut()
            .enumerate()
            .partition_map(|(i, rec)| match FitacfRecord::try_from(rec) {
                Err(e) => Either::Left((i, e)),
                Ok(x) => Either::Right(x),
            });
    if !errors.is_empty() {
        Err(PyErr::from(DmapError::InvalidRecord(format!(
            "Corrupted records: {errors:?}"
        ))))?
    }
    let grid_recs = fit2grid(&args, &formatted_recs)
        .map_err(PyErr::from)?
        .into_iter()
        .map(|rec| rec.inner())
        .collect();
    Ok(grid_recs)
}

/// Fits a list of FITACF files into a GRID file.
#[pyfunction]
#[pyo3(name = "fit2grid_file")]
#[pyo3(signature = (fitacf_files, grid_file, /, **py_kwargs))]
#[pyo3(text_signature = "(fitacf_files: list[str], grid_file: str, /, **kwargs)")]
fn fit2grid_file_py(fitacf_files: Vec<PathBuf>, grid_file: PathBuf, py_kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let args = match py_kwargs {
        Some(kwargs) => kwargs.extract()?,
        None => GridArgs::parse_from(vec!["fit2grid"]),
    };
    let grid_recs = fit2grid_file(&fitacf_files, &args)?;
    GridRecord::write_to_file(&grid_recs, &grid_file, false)?;
    Ok(())
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct GridArgsCLI {
    /// Output grid file path
    pub outfile: PathBuf,

    /// Fitacf file(s) to grid
    #[arg(num_args = 1.., last = true)]
    pub infiles: Vec<PathBuf>,

    #[command(flatten)]
    pub grid_args: GridArgs,
}

/// Converts a set of FITACF files into a GRID file.
#[pyfunction]
#[pyo3(name = "fit2grid_cli")]
fn fit2grid_cli(py: Python) -> PyResult<()> {
    let argv = py
        .import("sys")?
        .getattr("argv")?
        .extract::<Vec<String>>()?;
    let args = GridArgsCLI::parse_from(argv);
    let grid_records = fit2grid_file(&args.infiles, &args.grid_args)?;
    GridRecord::write_to_file(&grid_records, &args.outfile, false)?;
    Ok(())
}

/// Functions for SuperDARN data processing.
#[pymodule]
fn procdarn_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fitacf3_py, m)?)?;
    m.add_function(wrap_pyfunction!(fitacf3_file_py, m)?)?;
    m.add_wrapped(wrap_pyfunction!(fitacf3_cli))?;
    m.add_function(wrap_pyfunction!(fit2grid_py, m)?)?;
    m.add_function(wrap_pyfunction!(fit2grid_file_py, m)?)?;
    m.add_wrapped(wrap_pyfunction!(fit2grid_cli))?;

    Ok(())
}
