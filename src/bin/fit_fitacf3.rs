use clap::Parser;
use dmap::formats::fitacf::FitacfRecord;
use procdarn::fitting::fitacf3::fitacf_v3::{fit_rawacf_record, Fitacf3Error};
use procdarn::utils::rawacf::get_hdw;
use rayon::prelude::*;
use std::path::PathBuf;

pub type BinResult<T, E = Box<dyn std::error::Error + Send + Sync>> = Result<T, E>;

fn main() {
    if let Err(e) = bin_main() {
        eprintln!("error: {e}");
        if let Some(e) = e.source() {
            eprintln!("error: {e}");
        }
        std::process::exit(1);
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Rawacf file to fit
    #[arg()]
    infile: PathBuf,

    /// Output fitacf file path
    #[arg()]
    outfile: PathBuf,
}

fn bin_main() -> BinResult<()> {
    let args = Args::parse();

    let mut rawacf_records = dmap::read_rawacf(args.infile)?;
    let hdw = get_hdw(&rawacf_records[0])?;

    // Fit the records!
    let fitacf_results: Vec<Result<FitacfRecord, Fitacf3Error>> = rawacf_records
        .par_iter_mut()
        .map(|rec| fit_rawacf_record(rec, &hdw))
        .collect();

    let mut fitacf_records = vec![];
    for res in fitacf_results {
        match res {
            Ok(x) => fitacf_records.push(x),
            Err(e) => Err(e)?,
        }
    }

    // Write to file
    dmap::write_fitacf(fitacf_records, &args.outfile)?;
    Ok(())
}
