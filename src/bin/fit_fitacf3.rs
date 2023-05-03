use backscatter_rs::fitting::fitacf3::fitacf_v3::fit_rawacf_record;
use clap::Parser;
use dmap::formats::{to_file, DmapRecord, FitacfRecord, RawacfRecord};
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
    // let matches = Command::new(clap::crate_name!())
    //     .version(clap::crate_version!())
    //     .arg_required_else_help(true)
    //     .arg(Arg::new("infile")
    //         .long("infile")
    //         .short('i')
    //         .help("Rawacf file to fit")
    //         .num_args(1)
    //         .value_parser(value_parser!(PathBuf))
    //         .required(true))
    //     .arg(Arg::new("outfile")
    //         .long("outfile")
    //         .short('o')
    //         .help("Output fitacf file")
    //         .num_args(1)
    //         .value_parser(value_parser!(PathBuf))
    //         .required(true))
    //     .get_matches_from(wild::args_os());
    let args = Args::parse();

    let rawacf = File::open(args.infile)?;
    let rawacf_records = RawacfRecord::read_records(rawacf)?;
    let mut fitacf_records: Vec<FitacfRecord> = vec![];

    // Fit the records!
    let mut i = 0;
    for rec in rawacf_records {
        println!("Fitting record {}", i);
        fitacf_records.push(fit_rawacf_record(&rec)?);
        i += 1;
    }

    // Write to file
    to_file(args.outfile, &fitacf_records)?;
    Ok(())
}
