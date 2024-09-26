use backscatter_rs::fitting::fitacf3::fitacf_v3::{fit_rawacf_record, Fitacf3Error};
use backscatter_rs::utils::hdw::HdwInfo;
use chrono::NaiveDate;
use clap::Parser;
use dmap::write_fitacf;
use dmap::formats::{fitacf::FitacfRecord, rawacf::RawacfRecord, dmap::Record};
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
    /// Rawacf file to fit
    #[arg(short, long)]
    infile: PathBuf,

    /// Output fitacf file path
    #[arg(short, long)]
    outfile: PathBuf,
}

fn bin_main() -> BinResult<()> {
    let args = Args::parse();

    let rawacf = File::open(args.infile)?;
    let mut rawacf_records = RawacfRecord::read_records(rawacf)?;

    let rec = &rawacf_records[0];
    let file_datetime = NaiveDate::from_ymd_opt(
        rec.get(&"time.yr".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get time.yr".to_string()))?.clone().try_into()?, 
        rec.get(&"time.mo".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get time.mo".to_string()))?.clone().try_into()?, 
        rec.get(&"time.dy".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get time.dy".to_string()))?.clone().try_into()?
        ).unwrap()
        .and_hms_opt(
            rec.get(&"time.hr".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get time.hr".to_string()))?.clone().try_into()?, 
            rec.get(&"time.mt".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get time.mt".to_string()))?.clone().try_into()?, 
            rec.get(&"time.sc".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get time.sc".to_string()))?.clone().try_into()?
        ).ok_or_else(|| Fitacf3Error::Message("Unable to interpret record timestamp".to_string()))?;
    let hdw = HdwInfo::new(rec.get(&"stid".to_string()).ok_or_else(|| Fitacf3Error::Message("Could not get station ID".to_string()))?.clone().try_into()?, file_datetime)
        .map_err(|e| Fitacf3Error::Message(e.details))?;

    // Fit the records!
    let fitacf_records: Vec<FitacfRecord> = rawacf_records
        .par_iter_mut()
        .map(|rec| fit_rawacf_record(rec, &hdw).expect("Unable to fit record"))
        .collect();

    // Write to file
    write_fitacf(fitacf_records, &args.outfile)?;
    Ok(())
}
