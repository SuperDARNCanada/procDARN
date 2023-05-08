use backscatter_rs::fitting::fitacf3::fitacf_v3::{fit_rawacf_record, Fitacf3Error};
use backscatter_rs::utils::hdw::HdwInfo;
use chrono::NaiveDateTime;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use dmap::formats::{DmapRecord, FitacfRecord, RawacfRecord};

use dmap;
use rayon::prelude::*;
use std::fs::File;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Fitacf3", |b| b.iter(|| fitacf3()));
    c.bench_function("Parallel Fitacf3", |b| b.iter(|| rayon_fitacf3()));
}

fn fitacf3() {
    let file =
        File::open("tests/test_files/20210607.1801.00.cly.a.rawacf").expect("Test file not found");
    let rawacf = RawacfRecord::read_records(file).expect("Could not read records");
    let mut fitacf_records = vec![];

    let rec = &rawacf[0];
    let file_datetime = NaiveDateTime::parse_from_str(
        format!(
            "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
            rec.year, rec.month, rec.day, rec.hour, rec.minute, rec.second
        )
        .as_str(),
        "%Y%m%d %H:%M:%S",
    )
    .expect("Unable to interpret record timestamp");
    let hdw = HdwInfo::new(rec.station_id, file_datetime).expect("Unable to read hdw file");
    // fitacf_records.push(fit_rawacf_record(&rawacf[0]).expect("Could not fit rawacf record"));
    for rec in rawacf {
        fitacf_records.push(fit_rawacf_record(&rec, &hdw).expect("Could not fit record"));
    }
    dmap::formats::to_file("tests/test_files/temp.fitacf", &fitacf_records)
        .expect("Unable to write to file");
}

fn rayon_fitacf3() {
    let file =
        File::open("tests/test_files/20210607.1801.00.cly.a.rawacf").expect("Test file not found");
    let rawacf = RawacfRecord::read_records(file).expect("Could not read records");
    let fitacf_records: Vec<FitacfRecord>;

    let rec = &rawacf[0];
    let file_datetime = NaiveDateTime::parse_from_str(
        format!(
            "{:4}{:0>2}{:0>2} {:0>2}:{:0>2}:{:0>2}",
            rec.year, rec.month, rec.day, rec.hour, rec.minute, rec.second
        )
        .as_str(),
        "%Y%m%d %H:%M:%S",
    )
    .expect("Unable to interpret record timestamp");
    let hdw = HdwInfo::new(rec.station_id, file_datetime).expect("Unable to read utils file");

    fitacf_records = rawacf
        .par_iter()
        .map(|rec| fit_rawacf_record(&rec, &hdw).expect("Could not fit record"))
        .collect();
    dmap::formats::to_file("tests/test_files/temp.fitacf", &fitacf_records)
        .expect("Unable to write to file");
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
