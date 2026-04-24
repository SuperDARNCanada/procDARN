A library for the SuperDARN data processing chain
=================================================

[<img alt="github" src="https://img.shields.io/badge/github-SuperDARNCanada/procdarn-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/SuperDARNCanada/procdarn)

## Features

* `raw2fit`: RAWACF fitting routine using the FITACF3 algorithm, for converting to FITACF data.
* `fit2grid`: FITACF gridding routing paralleling the RST's `make_grid` binary.

## Installation

### Package manager
This package is registered on PyPI as `procdarn`, you can install with your package manager, e.g. `uv pip install procdarn`.

### From source
If you want to build from source, you first need to have Rust installed on your machine. Then:
1. Clone the repository: `git clone https://github.com/SuperDARNCanada/procdarn`
2. Create a virtual environment and source it, then install `maturin`
3. In the project directory, run `maturin develop` to build and install the Python bindings. This will make a wheel file based on your operating system and architecture that you can install directly on any compatible machine.

## Usage

The API is very simple, with only two functions exposed currently: `raw2fit` and `fit2grid`. Each of these functions
takes a `source` and optional `dest`, with behaviour determined by the type of the input. 

#### `raw2fit`
The above arguments are all that are currently supported for `raw2fit`. The FITACF3 algorithm is used; future plans include
support for additional fitting algorithms, e.g. LMFIT2.

```python
"""
Convert RAWACF records to FITACF records using the FITACF3 algorithm.

Parameters
----------
source: Union[str, list[dict]]
    Where to read data from. If input is of type `str`, this is interpreted as the path to a file.
    If input is of type `list[dict]`, this is interpreted as the raw data itself.
dest: Optional[str]
    Where to put fitted data. If `dest` is of type `str`, it is interpreted as the path to a file where the fitted
    data will be written to. If `dest` is `None`, the data is returned as `[list[dict]]`.

Returns
-------
If `dest` is `None`, the fitted data is returned as `[list[dict]]`.
"""
```

#### `fit2grid`
```python
"""
Convert FITACF records to GRID records.

Parameters
----------
source: Union[list[str], list[dict]]
    Where to read data from. If input is of type `list[str]`, this is interpreted as file paths.
    If input is of type `list[dict]`, this is interpreted as the raw data itself.
dest: Optional[str]
    Where to put fitted data. If `dest` is of type `str`, it is interpreted as the path to a file where the fitted
    data will be written to. If `dest` is `None`, the data is returned as `[list[dict]]`.
"""
```
There is also an extensive list of optional arguments to `fit2grid` that can be used to tweak the processing. 