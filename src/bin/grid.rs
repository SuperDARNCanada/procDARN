use backscatter_rs::gridding::filter::median_filter;
use backscatter_rs::gridding::grid::check_operational_params;
use backscatter_rs::gridding::grid_table::GridTable;
use backscatter_rs::utils::channel::{set_fix_channel, set_stereo_channel};
use backscatter_rs::utils::hdw::HdwInfo;
use backscatter_rs::utils::scan::RadarScan;
use backscatter_rs::utils::search::fit_seek;
use chrono::{Duration, NaiveDateTime};
use clap::{value_parser, Parser};
use dmap::formats::{to_file, DmapRecord, FitacfRecord, GridRecord};
use rayon::prelude::*;
use std::fs::File;
use std::path::PathBuf;

pub type BinResult<T, E = Box<dyn std::error::Error + Send + Sync>> = Result<T, E>;

/// Status of finding a record to grid. Done indicates certain initializations are done after Found
enum Status {
    Found,
    NotFound,
    Done,
}

fn main() {
    if let Err(e) = bin_main() {
        eprintln!("error: {e}");
        if let Some(e) = e.source() {
            eprintln!("error: {e}")
        }
        std::process::exit(1);
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Output grid file path
    #[arg()]
    outfile: PathBuf,

    /// Fitacf file(s) to grid
    #[arg(num_args = 1..)]
    infiles: Vec<PathBuf>,

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
    min_range_gate: Option<i32>,

    /// Maximum range gate
    #[arg(long, visible_alias = "maxrng")]
    max_range_gate: Option<i32>,

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
    max_power: f64,

    /// Maximum velocity in m/s
    #[arg(
        long,
        visible_alias = "vmax",
        value_parser,
        default_value = "60",
        requires = "op_param_flag"
    )]
    max_velocity: f64,

    /// Maximum spectral width in m/s
    #[arg(
        long,
        visible_alias = "wmax",
        value_parser,
        default_value = "1000",
        requires = "op_param_flag"
    )]
    max_spectral_width: f64,

    /// Maximum velocity error in m/s
    #[arg(
        long,
        visible_alias = "vemax",
        value_parser,
        default_value = "200",
        requires = "op_param_flag"
    )]
    max_velocity_error: f64,

    /// Minimum power (linear scale)
    #[arg(
        long,
        visible_alias = "pmin",
        value_parser,
        default_value = "35",
        requires = "op_param_flag"
    )]
    min_power: f64,

    /// Minimum velocity in m/s
    #[arg(
        long,
        visible_alias = "vmin",
        value_parser,
        default_value = "3",
        requires = "op_param_flag"
    )]
    min_velocity: f64,

    /// Minimum spectral width in m/s
    #[arg(
        long,
        visible_alias = "wmin",
        value_parser,
        default_value = "10",
        requires = "op_param_flag"
    )]
    min_spectral_width: f64,

    /// Minimum velocity error in m/s
    #[arg(
        long,
        visible_alias = "vemin",
        value_parser,
        default_value = "0",
        requires = "op_param_flag"
    )]
    min_velocity_error: f64,

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

fn bin_main() -> BinResult<()> {
    let args = Args::parse();

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
            panic!("Cannot interpret data exclusion flags [--ion, --gs, --both]")
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
    let mut current_scans: Vec<&RadarScan> = Vec::with_capacity(num_averages as usize);
    let mut grid_record: &RadarScan;
    let mut found_record = Status::NotFound;
    let mut index = 0;
    let mut num_scans = 0;
    let mut record_idx: Option<usize> = None;
    let mut end_time: NaiveDateTime;
    let mut hdw_info: Option<HdwInfo>;
    let mut records_for_file: Vec<GridRecord>;

    for infile in args.infiles.clone().into_iter() {
        let fitacf = File::open(infile)?;
        let fitacf_records = FitacfRecord::read_records(fitacf)?;

        // Get the first scan from the file
        let mut this_scan = RadarScan::get_first_scan(&fitacf_records, args.scan_length);
        if this_scan.is_err() {
            eprintln!(format!("Unable to get first scan from {:?}", infile));
            continue;
        }

        current_scans[index] = &this_scan.unwrap();
        let file_datetime =
            NaiveDateTime::from_timestamp_micros(current_scans[index].start_time * 1000.0 as i64);

        // Determine the starting time for gridding based on the record and input options
        let mut start_time = file_datetime?;
        if found_record == Status::NotFound {
            if let (None, None) = (&args.start_date, &args.start_time) {
                match NaiveDateTime::from_timestamp_micros(
                    current_scans[0].start_time * 1000.0 as i64,
                ) {
                    Some(t) => {
                        start_time = t;
                    }
                    None => panic!("Invalid timestamp in first scan"),
                }
                found_record = Status::Found;
            } else {
                let date_string = match &args.start_date {
                    Some(d) => d,
                    None => NaiveDateTime::from_timestamp_micros(
                        current_scans[0].start_time * 1000.0 as i64,
                    )?
                    .format("%Y%m%d")
                    .to_string(),
                };

                let time_string = match &args.start_time {
                    Some(t) => format!("{}", t),
                    // The None branch truncates back to the start of the minute
                    None => NaiveDateTime::from_timestamp_micros(
                        current_scans[0].start_time * 1000.0 as i64,
                    )?
                    .format("%H:%M")
                    .to_string(),
                };

                start_time = NaiveDateTime::parse_from_str(
                    format!("{} {}", date_string, time_string).as_str(),
                    "%Y%m%d %H:%M",
                )
                .map_err(|_| {
                    GridError::Message("Unable to parse date and/or time from options".to_string())
                })?;

                // If applying boxcar median filter then we need to load data prior to the usual start
                // time, so start_time needs to be adjusted
                if num_averages > 1 {
                    match args.scan_length {
                        Some(x) => start_time -= Duration::from_secs(x as u64),
                        None => {
                            start_time -= Duration::from_secs(
                                15 + current_scans[0].end_time - current_scans[0].start_time,
                            )
                        }
                    }
                }

                // Find the first record which occurs after the grid start time, if any
                if let Some((rec, idx)) = fit_seek(&fitacf_records, start_time) {
                    record_idx = Some(idx);
                } else {
                    eprintln!(
                        "Ignoring file {:?} as it ends before requested start time",
                        infile
                    );
                    continue;
                }
                found_record = Status::Found;

                // If using scan flag, go to the next beginning of the next scan
                if let None = args.scan_length {
                    record_idx = fitacf_records[record_idx..]
                        .iter()
                        .position(|rec| rec.scan_flag == 1);
                }

                // Read the first full scan of data corresponding to grid start datetime
                current_scans[0] = match record_idx {
                    Some(i) => &RadarScan::get_first_scan(&fitacf_records[i..], args.scan_length)?,
                    None => &RadarScan::get_first_scan(&fitacf_records, args.scan_length)?,
                };
            }
        }

        if found_record == Status::Found {
            end_time = match &args.end_time {
                Some(t) => {
                    let time_string = format!("{}", t);
                    let date_string = match &args.end_date {
                        Some(d) => format!("{}", d),
                        None => NaiveDateTime::from_timestamp_micros(
                            current_scans[0].start_time * 1000.0 as i64,
                        )?
                        .format("%Y%m%d")
                        .to_string(),
                    };
                    NaiveDateTime::parse_from_str(
                        format!("{} {}", date_string, time_string).as_str(),
                        "%Y%m%d %H:%M",
                    )
                    .map_err(|_| {
                        GridError::Message(
                            "Unable to parse end date and/or time from options".to_string(),
                        )
                    })?
                }
                None => match &args.interval {
                    Some(x) => {
                        start_time
                            + Duration::from_secs(
                                NaiveDateTime::parse_from_str(x?, "%H:%M")?.timestamp(),
                            )
                    }
                    None => panic!("No end time or interval specified for grid"),
                },
            };
            found_record = Status::NotFound
        }

        // TODO: Load AACGM coefficients

        num_scans += 1;

        // Grid all data until end of gridding time or end of file
        while this_scan.is_ok() {
            // Exclude scatter in beams listed in args.exclude_beams
            if let Some(b) = &args.exclude_beams {
                current_scans[index].reset_beams(b)?;
            }

            // Exclude data with scan flag == -1 if args.exclude_neg_scan_flag given
            if args.exclude_neg_scan_flag {
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
                current_scans[index].exclude_ionospheric_scatter();
            } else if args.ionosphere_only_flag {
                current_scans[index].exclude_groundscatter();
            }

            // Exclude scatter outside power, velocity, spectral width, and velocity error bounds
            if !args.no_limits_flag {
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
                if filter_weighting_mode != -1 {
                    match median_filter(
                        filter_weighting_mode,
                        current_scans.capacity() as i32,
                        index as i32,
                        15,
                        args.sort_params_flag,
                        &current_scans,
                    ) {
                        Ok(s) => {
                            grid_record = &s;
                        }
                        Err(e) => GridRecord(e),
                    }
                }
                grid_record = current_scans[index];

                // If not already done, load HdwInfo for radar
                if let None = hdw_info {
                    hdw_info = Some(HdwInfo::new(grid_record.station_id as i16, start_time)?);
                }

                // Test whether the grid table should be written to file
                if grid_table.test(grid_record) {
                    // If GridTable good and grid record starts at or after start_time, write to file
                    if grid_table.start_time >= start_time.timestamp() as f64 {
                        records_for_file.push(grid_table.to_dmap_record()?);
                    }
                }

                // Map GridTable to equal-area grid in magnetic coordinates
                grid_table.map(
                    grid_record,
                    &hdw_info?,
                    args.scan_length? as i32,
                    args.inertial_frame_flag,
                    args.altitude as f64,
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
            this_scan = match record_idx {
                Some(i) => RadarScan::get_first_scan(&fitacf_records[i..], args.scan_length),
                None => RadarScan::get_first_scan(&fitacf_records, args.scan_length),
            };

            if let Ok(scan) = this_scan {
                current_scans[index] = &scan;
            }

            // If scan starts after end_time, this file is done being gridded
            if current_scans[index].start_time > end_time.timestamp() as f64 {
                break;
            }

            num_scans += 1;
        }
    }

    // Write to file
    to_file(args.outfile, &records_for_file)?;
    Ok(())
}
