use backscatter_rs::fitting::fitacf3::fitacf_v3::{fit_rawacf_record, Fitacf3Error};
use backscatter_rs::hdw::hdw::HdwInfo;
use clap::Parser;
use dmap::formats::{to_file, DmapRecord, FitacfRecord, RawacfRecord};
use std::fs::File;
use std::path::PathBuf;
use chrono::NaiveDateTime;

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
    let rawacf_records = RawacfRecord::read_records(rawacf)?;
    let mut fitacf_records: Vec<FitacfRecord> = vec![];

    let rec = &rawacf_records[0];
    let file_datetime = NaiveDateTime::parse_from_str(
        format!(
            "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
            rec.year, rec.month, rec.day, rec.hour, rec.minute, rec.second
        )
        .as_str(),
        "%Y%m%d %H:%M:%S",
    )
    .map_err(|_| Fitacf3Error::Message("Unable to interpret record timestamp".to_string()))?;
    let hdw = HdwInfo::new(rec.station_id, file_datetime)
        .map_err(|e| Fitacf3Error::Message(e.details))?;

    // Fit the records!
    let mut i = 0;
    for rec in rawacf_records {
        fitacf_records.push(fit_rawacf_record(&rec, &hdw)?);
        i += 1;
    }

    // Write to file
    to_file(args.outfile, &fitacf_records)?;
    Ok(())
}
