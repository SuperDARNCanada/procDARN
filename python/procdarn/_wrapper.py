"""
Wrappers around the `procdarn` Python API.
"""

from typing import Union, Optional
import dmap
from . import procdarn_rs


def raw2fit(source: Union[str, list[dict]], dest: Optional[str] = None) -> Optional[list[dict]]:
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
    if not isinstance(source, list) and not isinstance(source, str):
        raise TypeError(
            f"invalid type for `source` {type(source)}: expected `str` or `list[dist]`"
        )
    if dest is not None and not isinstance(dest, str):
        raise TypeError(
            f"invalid type for `dest` {type(dest)}: expected `Optional[str]`"
        )
    if isinstance(source, str):
        if isinstance(dest, str):
            return procdarn_rs.fitacf3_file(source, dest)
        else:
            source = dmap.read_rawacf(source, mode="strict")
            return procdarn_rs.fitacf3_recs(source)
    else:
        if isinstance(dest, str):
            recs = procdarn_rs.fitacf3_recs(source)
            return dmap.write_fitacf(recs, dest)
        else:
            return procdarn_rs.fitacf3_recs(source)

def fit2grid(source: Union[list[str], list[dict]], dest: Optional[str] = None, **kwargs) -> Optional[list[dict]]:
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

    Keyword Arguments
    -----------------
    start_time: str
        Start time in HH:MM format
    end_time: str
        End time in HH:MM format
    start_date: str
        Start date in YYYYMMDD format
    end_date: str
        End date in YYYYMMDD format
    interval: str
        Use interval of length HH:MM
    scan_length: int
        Scan length specification in whole seconds, overriding the scan flag
    record_interval: int
        Time interval to store in each grid record, in whole seconds
    channel: str
        Stereo channel identifier, either 'a' or 'b'
    channel_fix: str
        User-defined channel identifier for the output file only
    exclude_beams: list[int]
        Beams to exclude
    min_range_gate: int
        Minimum range gate
    max_range_gate: int
        Maximum range gate
    min_slant_range: float
        Minimum slant range in km
    max_slant_range: float
        Maximum slant range in km
    filter_weighting: int
        Filter weighting mode [default: 0]
    max_power: float
        Maximum power (linear scale) [default: 60]
    max_velocity: float
        Maximum velocity in m/s [default: 2500]
    max_spectral_width: float
        Maximum spectral width in m/s [default: 1000]
    max_velocity_error: float
        Maximum velocity error in m/s [default: 200]
    min_power: float
        Minimum power (linear scale) [default: 3]
    min_velocity: float
        Minimum velocity in m/s [default: 35]
    min_spectral_width: float
        Minimum spectral width in m/s [default: 10]
    min_velocity_error: float
        Minimum velocity error in m/s [default: 0]
    altitude: float
        Altitude at which mapping is done in km [default: 300]
    max_frequency_var: int
        Maximum allowed frequency variation in Hz [default: 500000]
    boxcar_filter_flag: bool
        Flag to disable boxcar median filtering
    no_limits_flag: bool
        Flag to include data that exceeds limits
    op_param_flag: bool
        Flag to exclude data that doesn't match operating parameter requirements
    exclude_neg_scan_flag: bool
        Flag to exclude data with scan flag of -1
    extended_mode_flag: bool
        Extended output, include power and width in output file
    sort_params_flag: bool
        If using a median filter, sort parameters independent of the velocity
    ionosphere_only_flag: bool
        Exclude data marked as ground scatter
    groundscatter_only_flag: bool
        Exclude data not marked as ground scatter
    all_data_flag: bool
        Do not exclude data based on scatter flag
    inertial_frame_flag: bool
        Use inertial reference frame
    chisham_flag: bool
        Map data using Chisham virtual height model
    verbose: bool
        Verbose mode

    Returns
    -------
    If `dest` is `None`, the fitted data is returned as `[list[dict]]`.
    """
    if not isinstance(source, list) or len(source) == 0:
        raise TypeError(
            f"invalid type for `source` {type(source)}: expected non-empty `list[str]` or `list[dist]`"
        )
    if dest is not None and not isinstance(dest, str):
        raise TypeError(
            f"invalid type for `dest` {type(dest)}: expected `Optional[str]`"
        )

    if isinstance(source[0], str):
        if isinstance(dest, str):
            return procdarn_rs.fit2grid_file(source, dest, **kwargs)
        else:
            all_recs = []
            for f in source:
                all_recs.extend(dmap.read_fitacf(f, mode="strict"))
            return procdarn_rs.fit2grid_recs(all_recs, **kwargs)
    else:
        if isinstance(dest, str):
            recs = procdarn_rs.fit2grid_recs(source, **kwargs)
            return dmap.write_grid(recs, dest)
        else:
            return procdarn_rs.fit2grid_recs(source, **kwargs)
