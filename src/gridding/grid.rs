use std::os::raw::c_int;
use crate::error::ProcdarnError;
use crate::gridding::filter::{check_operational_params, median_filter};
use crate::gridding::grid_table::GridTable;
use crate::utils::channel::{set_fix_channel, set_stereo_channel};
use crate::utils::hdw::{HdwError, HdwInfo};
use crate::utils::scan::RadarScan;
use crate::utils::search::fit_seek;
use chrono::{DateTime, Datelike, NaiveDateTime, NaiveTime, TimeDelta, Utc};
use clap::Parser;
use dmap::error::DmapError;
use dmap::formats::grid::GridRecord;
use std::path::PathBuf;
use pyo3::exceptions::PyValueError;
use pyo3::PyErr;
use thiserror::Error;

/// Enum of the possible error variants that may be encountered
#[derive(Error, Debug)]
pub enum GridError {
    /// Represents an error in the Fitacf record that is attempting to be gridded
    #[error("{0}")]
    InvalidFitacf(String),

    /// Represents a field having the wrong type
    #[error("{0}")]
    WrongType(&'static str),

    /// Represents an error in processing of the record, for any reason
    #[error("{0}")]
    ProcessingError(#[from] ProcdarnError),

    /// Represents an error parsing a date, time, or datetime
    #[error("{0}")]
    Datetime(#[from] chrono::ParseError),

    /// Unable to get hardware file information
    #[error("{0}")]
    Hdw(#[from] HdwError),

    /// Invalid DMAP file
    #[error("{0}")]
    Dmap(#[from] DmapError),

    /// Error in `igrf` crate
    #[error("{0}")]
    Igrf(#[from] igrf::Error),

    /// Error in `geodesy` crate
    #[error("{0}")]
    Geodesy(#[from] geodesy::Error),

    /// Error in argument specification
    #[error("{0}")]
    BadArgs(String),
}

impl From<GridError> for PyErr {
    fn from(value: GridError) -> Self {
        let msg = value.to_string();
        PyValueError::new_err(msg)
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct GridArgs {
    /// Output grid file path
    #[arg()]
    pub outfile: PathBuf,

    /// Fitacf file(s) to grid
    #[arg(num_args = 1..)]
    pub infiles: Vec<PathBuf>,

    /// Start time in HH:MM format
    #[arg(long, visible_alias = "st")]
    start_time: Option<String>,

    /// End time in HH:MM format
    #[arg(long, visible_alias = "et")]
    end_time: Option<String>,

    /// Start date in YYYYMMDD format
    #[arg(long, visible_alias = "sd")]
    start_date: Option<String>,

    /// End date in YYYYMMDD format
    #[arg(long, visible_alias = "ed")]
    end_date: Option<String>,

    /// Use interval of length HH:MM
    #[arg(long, visible_alias = "ex", conflicts_with_all = &["end_time", "end_date"])]
    interval: Option<String>,

    /// Scan length specification in whole seconds, overriding the scan flag
    #[arg(long, visible_alias = "tl")]
    scan_length: Option<u32>,

    /// Time interval to store in each grid record, in whole seconds
    #[arg(short = 'i', long, value_parser, default_value = "120")]
    record_interval: u32,

    /// Stereo channel identifier, either 'a' or 'b'
    #[arg(long, visible_alias = "cn", value_parser)]
    channel: Option<char>,

    /// User-defined channel identifier for the output file only
    #[arg(long, visible_alias = "cn_fix", conflicts_with = "channel")]
    channel_fix: Option<char>,

    /// Beams to exclude, as a comma-separated list
    #[arg(long, visible_alias = "ebm", value_delimiter = ',', value_parser)]
    exclude_beams: Option<Vec<i32>>,

    /// Minimum range gate
    #[arg(long, visible_alias = "minrng")]
    min_range_gate: Option<usize>,

    /// Maximum range gate
    #[arg(long, visible_alias = "maxrng")]
    max_range_gate: Option<usize>,

    /// Minimum slant range in km
    #[arg(long, visible_alias = "minsrng")]
    min_slant_range: Option<f32>,

    /// Maximum slant range in km
    #[arg(long, visible_alias = "maxsrng")]
    max_slant_range: Option<f32>,

    /// Filter weighting mode
    #[arg(long, visible_alias = "fwgt", value_parser, default_value = "0")]
    filter_weighting: i32,

    /// Maximum power (linear scale)
    #[arg(
        long,
        visible_alias = "pmax",
        value_parser,
        default_value = "2500",
        requires = "op_param_flag"
    )]
    max_power: f32,

    /// Maximum velocity in m/s
    #[arg(
        long,
        visible_alias = "vmax",
        value_parser,
        default_value = "60",
        requires = "op_param_flag"
    )]
    max_velocity: f32,

    /// Maximum spectral width in m/s
    #[arg(
        long,
        visible_alias = "wmax",
        value_parser,
        default_value = "1000",
        requires = "op_param_flag"
    )]
    max_spectral_width: f32,

    /// Maximum velocity error in m/s
    #[arg(
        long,
        visible_alias = "vemax",
        value_parser,
        default_value = "200",
        requires = "op_param_flag"
    )]
    max_velocity_error: f32,

    /// Minimum power (linear scale)
    #[arg(
        long,
        visible_alias = "pmin",
        value_parser,
        default_value = "35",
        requires = "op_param_flag"
    )]
    min_power: f32,

    /// Minimum velocity in m/s
    #[arg(
        long,
        visible_alias = "vmin",
        value_parser,
        default_value = "3",
        requires = "op_param_flag"
    )]
    min_velocity: f32,

    /// Minimum spectral width in m/s
    #[arg(
        long,
        visible_alias = "wmin",
        value_parser,
        default_value = "10",
        requires = "op_param_flag"
    )]
    min_spectral_width: f32,

    /// Minimum velocity error in m/s
    #[arg(
        long,
        visible_alias = "vemin",
        value_parser,
        default_value = "0",
        requires = "op_param_flag"
    )]
    min_velocity_error: f32,

    /// Altitude at which mapping is done in km
    #[arg(long, visible_alias = "alt", value_parser, default_value = "300")]
    altitude: f32,

    /// Maximum allowed frequency variation in Hz
    #[arg(long, visible_alias = "fmax", value_parser, default_value = "500000")]
    max_frequency_var: i32,

    /// Flag to disable boxcar median filtering
    #[arg(long, visible_alias = "nav", action = clap::ArgAction::SetFalse)]
    boxcar_filter_flag: bool,

    /// Flag to include data that exceeds limits
    #[arg(long, visible_alias = "nlm", action = clap::ArgAction::SetTrue)]
    no_limits_flag: bool,

    /// Flag to exclude data that doesn't match operating parameter requirements
    #[arg(long, visible_alias = "nb", action = clap::ArgAction::SetTrue)]
    op_param_flag: bool,

    /// Flag to exclude data with scan flag of -1
    #[arg(long, visible_alias = "ns", action = clap::ArgAction::SetTrue)]
    exclude_neg_scan_flag: bool,

    /// Extended output, include power and width in output file
    #[arg(long, visible_alias = "xtd", action = clap::ArgAction::SetTrue)]
    extended_mode_flag: bool,

    /// If using a median filter, sort parameters independent of the velocity
    #[arg(long, visible_alias = "isort", action = clap::ArgAction::SetTrue)]
    sort_params_flag: bool,

    /// Exclude data marked as ground scatter
    #[arg(long, visible_alias = "ion", default_value = "true", action = clap::ArgAction::SetTrue)]
    ionosphere_only_flag: bool,

    /// Exclude data not marked as ground scatter
    #[arg(long, visible_alias = "gs", action = clap::ArgAction::SetTrue,
    conflicts_with = "ionosphere_only_flag")]
    groundscatter_only_flag: bool,

    /// Do not exclude data based on scatter flag
    #[arg(long, visible_alias = "both", action = clap::ArgAction::SetTrue,
    conflicts_with_all = &["ionosphere_only_flag", "groundscatter_only_flag"])]
    all_data_flag: bool,

    /// Use inertial reference frame
    #[arg(long, visible_alias = "inertial", action = clap::ArgAction::SetTrue)]
    inertial_frame_flag: bool,

    /// Map data using Chisham virtual height model
    #[arg(long, visible_alias = "chisham", action = clap::ArgAction::SetTrue)]
    chisham_flag: bool,

    /// Map data using old AACGM coefficients, rather than v2
    #[arg(long, visible_alias = "old_aacgm", action = clap::ArgAction::SetTrue)]
    old_aacgm_flag: bool,

    /// Verbose mode
    #[arg(short, long, action = clap::ArgAction::SetTrue)]
    verbose: bool,
}

/// Takes a list of fitacf files and converts them into grid files.
pub fn fit2grid(args: &GridArgs) -> Result<Vec<GridRecord>, GridError> {
    let mut grid_table = GridTable {
        ..Default::default()
    };

    // If "filter_weighting_mode" greater than 0, decrement it by one
    let mut filter_weighting_mode = args.filter_weighting;
    if filter_weighting_mode > 0 {
        filter_weighting_mode -= 1;
    };

    // Set GridTable groundscatter flag
    grid_table.groundscatter = {
        if args.groundscatter_only_flag == true {
            0
        } else if args.ionosphere_only_flag == true {
            1
        } else if args.all_data_flag == true {
            2
        } else {
            return Err(GridError::BadArgs(
                "Cannot interpret data exclusion flags [--ion, --gs, --both]".to_string(),
            ));
        }
    };

    // Set GridTable channel number
    grid_table.channel = {
        if let Some(c) = args.channel {
            set_stereo_channel(c).unwrap_or(-1) // Determine stereo channel, either 'a' or 'b'
        } else if let Some(c) = args.channel_fix {
            set_fix_channel(c).unwrap_or(-1) // Determine appropriate channel for output file
        } else {
            0
        }
    };

    // Store bounding thresholds for power, velocity, velocity error, and spectral width in GridTable
    if !args.no_limits_flag {
        grid_table.min_power = args.min_power;
        grid_table.min_velocity = args.min_velocity;
        grid_table.min_spectral_width = args.min_spectral_width;
        grid_table.min_velocity_error = args.min_velocity_error;
        grid_table.max_power = args.max_power;
        grid_table.max_velocity = args.max_velocity;
        grid_table.max_spectral_width = args.max_spectral_width;
        grid_table.max_velocity_error = args.max_velocity_error;
    }

    // Initialize the size of the boxcar. Default 3 if median filtering being applied, 1 otherwise
    let num_averages: i32;
    if args.boxcar_filter_flag {
        num_averages = 3;
    } else {
        num_averages = 1;
    }
    // Preallocate memory for a vector of records that will be boxcar filtered
    let mut current_scans: Vec<RadarScan> = Vec::with_capacity(num_averages as usize);
    for _ in 0..num_averages {
        current_scans.push(RadarScan::default());
    }
    let mut found_record = false;
    let mut index = 0;
    let mut num_scans = 0;
    let mut record_idx: Option<usize> = None;
    let mut end_time: DateTime<Utc> = DateTime::default();
    let hdw_info: Option<HdwInfo> = None;
    let mut records_for_file: Vec<GridRecord> = vec![];
    let mut found_scan: bool;

    for infile in args.infiles.clone().into_iter() {
        println!("Fitting file {}", infile.display());
        let fitacf_records = dmap::read_fitacf(infile.clone())?;

        // Get the first scan from the file
        match RadarScan::get_first_scan(&fitacf_records, args.scan_length) {
            Ok((x, _)) => {
                println!("current_scans: len({}), index={index}", current_scans.len());
                current_scans[index] = x;
                found_scan = true;
            }
            Err(e) => {
                eprintln!("Unable to get first scan from {} - {e}", infile.display());
                continue;
            }
        };

        let file_datetime = current_scans[index].start_time;

        println!("file starts at {file_datetime}"); // todo: remove

        // Determine the starting time for gridding based on the record and input options
        let mut start_time = file_datetime;
        if !found_record {
            if let (None, None) = (&args.start_date, &args.start_time) {
                start_time = current_scans[0].start_time;
                found_record = true;
            } else {
                let date_string = match &args.start_date {
                    Some(d) => d.clone(),
                    None => current_scans[0].start_time.format("%Y%m%d").to_string()
                };

                let time_string = match &args.start_time {
                    Some(t) => format!("{}", t),
                    // The None branch truncates back to the start of the minute
                    None => current_scans[0].start_time.format("%H:%M").to_string(),
                };

                start_time = NaiveDateTime::parse_from_str(
                    format!("{} {}", date_string, time_string).as_str(),
                    "%Y%m%d %H:%M",
                )
                .map_err(|_| {
                    ProcdarnError::Timestamp(
                        "Unable to parse date and/or time from options or file".to_string(),
                    )
                })?
                .and_utc();
                println!("start_time: {start_time}");
                // If applying boxcar median filter then we need to load data prior to the usual start
                // time, so start_time needs to be adjusted
                if num_averages > 1 {
                    match args.scan_length {
                        Some(x) => {
                            start_time -= TimeDelta::new(x as i64, 0).ok_or_else(|| {
                                ProcdarnError::Timestamp(
                                    "Out of bounds duration when adjusting start_time".to_string(),
                                )
                            })?
                        }
                        None => {
                            let td = current_scans[0].end_time - current_scans[0].start_time + TimeDelta::seconds(15);
                            start_time -= td;
                        }
                    }
                }

                // Find the first record which occurs after the grid start time, if any
                if let Ok(Some((_, idx))) = fit_seek(&fitacf_records, start_time) {
                    record_idx = Some(idx);
                } else {
                    eprintln!(
                        "Ignoring file {} as it ends before requested start time",
                        infile.display()
                    );
                    continue;
                }
                found_record = true;

                // If using scan flag, go to the next beginning of the next scan
                if let None = args.scan_length {
                    if let Some(x) = record_idx {
                        let mut scan_flags: Vec<i16> = vec![];
                        for rec in fitacf_records[x..].iter() {
                            let scn_flg = i16::try_from(
                                rec.get(&"scan".to_string())
                                    .ok_or_else(|| {
                                        GridError::InvalidFitacf(format!(
                                            "missing `scan` flag in {}", infile.display()
                                        ))
                                    })?
                                    .clone(),
                            ).map_err(|_| {
                                GridError::InvalidFitacf(format!(
                                    "bad `scan` flag in {}", infile.display()
                                ))
                            })?;
                            scan_flags.push(scn_flg);
                        }
                        record_idx = Some(scan_flags.iter().position(|&flg| flg == 1).ok_or_else(|| {
                            GridError::InvalidFitacf(format!(
                                "No records with set `scan` flag in {}", infile.display()
                            ))
                        })?);
                    } else {
                        return Err(GridError::BadArgs(
                            "No records match requested scan time".to_string(),
                        ));
                    }
                }

                // Read the first full scan of data corresponding to grid start datetime
                let first_idx = record_idx.unwrap_or_else(|| 0);
                let this_scan =
                    RadarScan::get_first_scan(&fitacf_records[first_idx..], args.scan_length);
                found_scan = this_scan.is_ok();
                current_scans[0] = this_scan?.0;
            }
        }

        if found_record {
            end_time = match &args.end_time {
                Some(t) => {
                    let time_string = format!("{}", t);
                    let date_string = match &args.end_date {
                        Some(d) => format!("{}", d),
                        None => current_scans[0].start_time.format("%Y%m%d").to_string(),
                    };
                    NaiveDateTime::parse_from_str(
                        format!("{} {}", date_string, time_string).as_str(),
                        "%Y%m%d %H:%M",
                    )
                    .map_err(|_| {
                        GridError::BadArgs(
                            "Unable to parse end date and/or time from options".to_string(),
                        )
                    })?
                    .and_utc()
                }
                None => match &args.interval {
                    Some(x) => {
                        let dt = NaiveTime::parse_from_str(x, "%H:%M")?;
                        let dur = dt - NaiveTime::from_hms_opt(0, 0, 0).ok_or_else(|| GridError::BadArgs("This should never happen, trying to make NaiveTime::from_hms(0, 0, 0)".to_string()))?;
                        start_time + dur
                    }
                    None => {
                        return Err(GridError::BadArgs(
                            "No end time or interval specified for grid".to_string(),
                        ))
                    }
                },
            };
            found_record = false
        }

        let year = start_time.year() as c_int;
        let month = start_time.month() as c_int;
        let day = start_time.day() as c_int;
        unsafe {
            aacgmv2_rs::AACGM_v2_SetDateTime(year, month, day, 0, 0, 0);
        }
        num_scans += 1;

        // Grid all data until end of gridding time or end of file
        while found_scan {
            // Exclude scatter in beams listed in args.exclude_beams
            if let Some(b) = &args.exclude_beams {
                eprintln!("excluding beams {:?}", &args.exclude_beams);
                current_scans[index].reset_beams(b)?;
            }

            // Exclude data with scan flag == -1 if args.exclude_neg_scan_flag given
            if args.exclude_neg_scan_flag {
                eprintln!("excluding data with negative scan flag");
                current_scans[index].exclude_outofscan();
            }

            // Exclude scatter in range gates below args.min_range_gate or above args.max_range_gate
            current_scans[index].exclude_range(
                args.min_range_gate,
                args.max_range_gate,
                args.min_slant_range,
                args.max_slant_range,
            );

            // Exclude groundscatter or ionospheric scatter, depending on the args given
            if args.groundscatter_only_flag {
                eprintln!("excluding ionospheric scatter");
                current_scans[index].exclude_ionospheric_scatter();
            } else if args.ionosphere_only_flag {
                eprintln!("excluding ground scatter");
                current_scans[index].exclude_groundscatter();
            }

            // Exclude scatter outside power, velocity, spectral width, and velocity error bounds
            if !args.no_limits_flag {
                eprintln!("Excluding out of bounds");
                current_scans[index].exclude_outofbounds(&grid_table);
            }

            // If enough scans have been loaded and args.no_limit not given, check to make sure the
            // first range, range separation, and transmit frequency have not changed significantly
            let mut passed_check = true;
            if num_scans >= current_scans.capacity()
                && !args.no_limits_flag
                && filter_weighting_mode != -1
            {
                passed_check = check_operational_params(&current_scans, args.max_frequency_var);
            }

            // If enough scans have been loaded, proceed with filtering and gridding
            if passed_check && num_scans >= current_scans.capacity() {
                let grid_record = match filter_weighting_mode {
                    -1 => current_scans[index].clone(),
                    _ => median_filter(
                        filter_weighting_mode,
                        current_scans.capacity() as u32,
                        index as i32,
                        15,
                        args.sort_params_flag,
                        &current_scans,
                    )?,
                };

                // If not already done, load HdwInfo for radar
                let hdw_params = match hdw_info {
                    Some(ref x) => x,
                    None => &HdwInfo::new(grid_record.station_id as i16, start_time)?,
                };

                // Test whether the grid table should be written to file
                if grid_table.test(&grid_record) {
                    // If GridTable good and grid record starts at or after start_time, write to file
                    if grid_table.start_time >= start_time {
                        records_for_file.push(grid_table.to_dmap_record()?);
                    }
                }

                // Map GridTable to equal-area grid in magnetic coordinates
                grid_table.map(
                    &grid_record,
                    &hdw_params,
                    args.scan_length.unwrap_or(0) as i32,
                    args.inertial_frame_flag,
                    args.altitude,
                    args.chisham_flag,
                    args.old_aacgm_flag,
                )?;
            }

            // Update index
            index += 1;
            if index >= current_scans.capacity() {
                index = 0;
            }

            // Get the next scan
            let start_idx = record_idx.unwrap_or_else(|| 0);
            let scan_res = RadarScan::get_first_scan(&fitacf_records[start_idx..], args.scan_length);
            match scan_res {
                Ok((new_scan, num_read)) => {
                    found_scan = true;
                    current_scans[index] = new_scan;
                    record_idx = Some(start_idx + num_read);
                },
                Err(ProcdarnError::ZeroRecords(_)) => {
                    found_scan = false;
                    record_idx = None;
                },
                Err(e) => Err(e)?
            };

            // If scan starts after end_time, this file is done being gridded
            if current_scans[index].start_time > end_time {
                break;
            }
            num_scans += 1;
        }
    }

    // Write to file
    Ok(records_for_file)
}
