// use backscatter_rs::gridding::grid::{grid_fitacf_record, GridError};
use backscatter_rs::utils::hdw::HdwInfo;
use chrono::NaiveDateTime;
use clap::Parser;
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
    scan_length: Option<i32>,

    /// Time interval to store in each grid record, in whole seconds
    #[arg(short = 'i', long)]
    record_interval: Option<i32>,

    /// Stereo channel identifier, either 'a' or 'b'
    #[arg(long, visible_alias = "cn")]
    channel: Option<char>,

    /// User-defined channel identifier for the output file only
    #[arg(long, visible_alias = "cn_fix")]
    channel_fix: Option<char>,

    /// Beams to exclude, as a comma-separated list
    #[arg(long, visible_alias = "ebm")]
    exclude_beams: Option<String>,

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

    /// Flag to use filter weighting mode
    #[arg(long, visible_alias = "fwgt", action = clap::ArgAction::SetTrue)]
    filter_weighting: bool,

    /// Maximum power in dB
    #[arg(long, visible_alias = "pmax")]
    max_power: Option<f32>,

    /// Maximum velocity in m/s
    #[arg(long, visible_alias = "vmax")]
    max_velocity: Option<f32>,

    /// Maximum spectral width in m/s
    #[arg(long, visible_alias = "wmax")]
    max_spectral_width: Option<f32>,

    /// Maximum velocity error in m/s
    #[arg(long, visible_alias = "vemax")]
    max_velocity_error: Option<f32>,

    /// Minimum power in dB
    #[arg(long, visible_alias = "pmin")]
    min_power: Option<f32>,

    /// Minimum velocity in m/s
    #[arg(long, visible_alias = "vmin")]
    min_velocity: Option<f32>,

    /// Minimum spectral width in m/s
    #[arg(long, visible_alias = "wmin")]
    min_spectral_width: Option<f32>,

    /// Minimum velocity error in m/s
    #[arg(long, visible_alias = "vemin")]
    min_velocity_error: Option<f32>,

    /// Altitude at which mapping is done in km
    #[arg(long, visible_alias = "alt")]
    altitude: Option<f32>,

    /// Maximum allowed frequency variation in Hz
    #[arg(long, visible_alias = "fmax")]
    max_frequency_var: Option<i32>,

    /// Flag to disable boxcar median filtering
    #[arg(long, visible_alias = "nav", action = clap::ArgAction::SetFalse)]
    boxcar_filter_flag: bool,

    /// Flag to include data that exceeds limits
    #[arg(long, visible_alias = "nlm", action = clap::ArgAction::SetTrue)]
    no_limits_flag: bool,

    /// Flag to include data that doesn't match operating parameter requirements
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
    #[arg(long, visible_alias = "ion", action = clap::ArgAction::SetTrue)]
    ionosphere_only_flag: bool,

    /// Exclude data not marked as ground scatter
    #[arg(long, visible_alias = "gs", action = clap::ArgAction::SetTrue)]
    groundscatter_only_flag: bool,

    /// Do not exclude data based on scatter flag
    #[arg(long, visible_alias = "both", action = clap::ArgAction::SetTrue)]
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
}

fn bin_main() -> BinResult<()> {
    let args = Args::parse();

    // If "channel" set then determine stereo channel, either 'A' or 'B'
    // If "channel_fix" set then determine appropriate channel for output file
    // If "exclude_beams" set then parse the beam list
    // If "interval" set convert to seconds
    // If "start_time" set convert to seconds
    // If "end_time" set convert to seconds
    // If "start_date" set convert to seconds since epoch
    // If "end_date" set convert to seconds since epoch
    // If "filter_weighting_mode" greater than 0, decrement it by one
    // Set GridTable groundscatter flag
    // Set GridTable channel number
    // Store bounding thresholds for power, velocity, velocity error, and spectral width in GridTable
    // Initialize the size of the boxcar. Default 3 if median filtering being applied, 1 otherwise

    // let fitacf = File::open(args.infile)?;
    // let fitacf_records = FitacfRecord::read_records(fitacf)?;
    //
    // let rec = &fitacf_records[0];
    // let file_datetime = NaiveDateTime::parse_from_str(
    //     format!(
    //         "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
    //         rec.year, rec.month, rec.day, rec.hour, rec.minute, rec.second
    //     )
    //     .as_str(),
    //     "%Y%m%d %H:%M:%S",
    // )
    // .map_err(|_| GridError::Message("Unable to interpret record timestamp".to_string()))?;
    // let hdw = HdwInfo::new(rec.station_id, file_datetime)
    //     .map_err(|e| GridError::Message(e.details))?;
    //
    // // Fit the records!
    // let grid_records: Vec<GridRecord> = fitacf_records
    //     .par_iter()
    //     .map(|rec| grid_fitacf_record(rec, &hdw).expect("Unable to fit record"))
    //     .collect();
    //
    // // Write to file
    // to_file(args.outfile, &grid_records)?;
    Ok(())
}
