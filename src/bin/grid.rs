// use backscatter_rs::gridding::grid::{grid_fitacf_record, GridError};
use backscatter_rs::gridding::grid_table::GridTable;
use backscatter_rs::utils::channel::{set_fix_channel, set_stereo_channel};
use backscatter_rs::utils::hdw::HdwInfo;
use chrono::{Duration, NaiveDateTime};
use clap::{value_parser, Parser};
use dmap::formats::{to_file, DmapRecord, FitacfRecord, RawacfRecord};
use rayon::prelude::*;
use std::fs::File;
use std::path::PathBuf;

pub type BinResult<T, E = Box<dyn std::error::Error + Send + Sync>> = Result<T, E>;

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
    #[arg(long, visible_alias = "ex")]
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

    // If "interval" set convert to seconds
    // If "start_time" set convert to seconds
    // If "end_time" set convert to seconds
    // If "start_date" set convert to seconds since epoch
    // If "end_date" set convert to seconds since epoch

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
    if args.op_param_flag {
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
    let mut current_records: Vec<&FitacfRecord> = Vec::with_capacity(num_averages as usize);

    let mut found_record = false;
    let mut index = 0;

    for infile in args.infiles.clone().into_iter() {
        let fitacf = File::open(infile)?;
        let fitacf_records = FitacfRecord::read_records(fitacf)?;

        current_records[index] = &fitacf_records[0];
        let file_datetime = NaiveDateTime::parse_from_str(
            format!(
                "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
                current_records[0].year,
                current_records[0].month,
                current_records[0].day,
                current_records[0].hour,
                current_records[0].minute,
                current_records[0].second
            )
            .as_str(),
            "%Y%m%d %H:%M:%S",
        )
        .map_err(|_| GridError::Message("Unable to interpret record timestamp".to_string()))?;

        // Determine the starting time for gridding based on the record and input options
        let mut start_time;
        if let (None, None) = (args.start_date, args.start_time) {
            start_time = NaiveDateTime::parse_from_str(
                format!(
                    "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
                    current_records[0].year,
                    current_records[0].month,
                    current_records[0].day,
                    current_records[0].hour,
                    current_records[0].minute,
                    current_records[0].second
                )
                .as_str(),
                "%Y%m%d %H:%M",
            )
            .map_err(|_| GridError::Message("Unable to interpret record timestamp".to_string()))?;
        } else {
            let date_string = match args.start_date {
                Some(d) => format!("{}", d),
                None => format!(
                    "{:4}{:0>2}{:0>2}",
                    current_records[0].year, current_records[0].month, current_records[0].day
                ),
            };

            let time_string = match args.start_time {
                Some(t) => format!("{}", t),
                None => format!(
                    "{:0>2}:{:0>2}:{:0>2}",
                    current_records[0].hour, current_records[0].minute, current_records[0].second
                ),
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
                    Some(x) => { start_time -= Duration::from_secs(x as u64) },
                    None => { start_time -= Duration::from_secs(15 + current_records[0].end_time - current_records[0].start_time)}
                }
            }
        }

        let hdw = HdwInfo::new(rec.station_id, file_datetime)
            .map_err(|e| GridError::Message(e.details))?;

        // Grid the records!
        let grid_records: Vec<GridRecord> = fitacf_records
            .par_iter()
            .map(|rec| grid_fitacf_record(rec, &hdw).expect("Unable to fit record"))
            .collect();

        // Write to file
        to_file(args.outfile, &grid_records)?;
    }

    Ok(())
}
