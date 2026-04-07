use crate::fitting::fitacf3::fitacf_v3::{par_fitacf3, Fitacf3Error};
use crate::gridding::grid::{fit2grid, fit2grid_file, GridArgs};
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
use pyo3::types::PyDict;

pub mod error;
pub mod fitting;
pub mod gridding;
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
fn file_fitacf3(raw_file: PathBuf, fit_file: PathBuf) -> Result<(), Fitacf3Error> {
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


/// Converts a list of FITACF records into GRID records.
#[pyfunction]
#[pyo3(name = "fit2grid")]
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
    dmap::write_grid(grid_recs, &grid_file)?;
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
        .import_bound("sys")?
        .getattr("argv")?
        .extract::<Vec<String>>()?;
    let args = GridArgsCLI::parse_from(argv);
    let grid_records = fit2grid_file(&args.infiles, &args.grid_args)?;
    dmap::write_grid(grid_records, &args.outfile)?;
    Ok(())
}

/// Functions for SuperDARN data processing.
#[pymodule]
fn procdarn(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fitacf3_py, m)?)?;
    m.add_function(wrap_pyfunction!(file_fitacf3_py, m)?)?;
    m.add_wrapped(wrap_pyfunction!(fitacf3_cli))?;
    m.add_function(wrap_pyfunction!(fit2grid_py, m)?)?;
    m.add_function(wrap_pyfunction!(fit2grid_file_py, m)?)?;
    m.add_wrapped(wrap_pyfunction!(fit2grid_cli))?;

    Ok(())
}
