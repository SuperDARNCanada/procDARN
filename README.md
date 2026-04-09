Core SuperDARN Processing Tools
===============================

[![github]](https://github.com/SuperDARNCanada/procdarn)

[github]: https://img.shields.io/badge/github-8da0cb?style=for-the-badge&labelColor=555555&logo=github

## Installation

### Rust
Add the crate to your dependencies in your `Cargo.toml` file, specifying as a git dependency.

### Python
`pip install procdarn`

### From source
If you want to build from source, you first need to have Rust installed on your machine. Then:
1. Clone the repository: `git clone https://github.com/SuperDARNCanada/procdarn`
2. Run `cargo build` in the repository directory
3. If wanting to install the Python API, create a virtual environment and source it, then install `maturin`
4. In the project directory, run `maturin develop` to build and install the Python bindings. This will make a wheel file based on your operating system and architecture that you can install directly on any compatible machine.

## Python API

### In-memory fitting
Use FITACF3 algorithm to fit a list of RAWACF records
```python
fitacf3(recs: list[dict]) -> list[dict]
```

#### Example
```python
import procdarn
import dmap
infile = "path/to/rawacf"
rawacf_records = dmap.read_rawacf(infile)
fitacf_records = procdarn.fitacf3(rawacf_records)
dmap.write_fitacf(fitacf_records, "path/to/fitacf")
```

### File-to-file fitting

Fit a RAWACF file into a FITACF file
```python
fitacf3_file(rawacf_file: str, fitacf_file: str)
```

#### Example
```python
import procdarn
procdarn.fitacf3_file("path/to/rawacf", "path/to/fitacf")
```

## Rust API

* `procdarn::fitacf3_file(raw_file: PathBuf, fit_file: PathBuf) -> Result<(), Fitacf3Error>`: Fits a RAWACF file into a FITACF file
* `procdarn::fitacf3(Vec<RawacfRecord>) -> Result<Vec<FitacfRecord>, Fitacf3Error>`: parallelized FITACFv3 on collection of `RawacfRecord`s (from `dmap` crate)
* `procdarn::fitacf3_single_threaded(Vec<RawacfRecord>) -> Result<Vec<FitacfRecord>, Fitacf3Error>`: single-threaded FITACFv3 implementation

## Binary
`raw2fit`: Command-line tool for fitting a RAWACF file using the FITACF3 algorithm.

```
Usage: raw2fit <INFILE> <OUTFILE>

Arguments:
  <INFILE>   Rawacf file to fit
  <OUTFILE>  Output fitacf file path

Options:
  -h, --help     Print help
  -V, --version  Print version
```