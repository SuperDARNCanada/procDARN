use backscatter_rs::fitting::fitacf3::fitacf_v3::fit_rawacf_record;
use backscatter_rs::utils::hdw::HdwInfo;
use chrono::NaiveDate;
use criterion::{criterion_group, criterion_main, Criterion};
use dmap;
use dmap::formats::{dmap::Record, rawacf::RawacfRecord};
use rayon::prelude::*;
use std::fs::File;

const TEST_FILE: &str = "tests/test_files/large.rawacf";

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Fitacf3", |b| b.iter(|| fitacf3()));
    c.bench_function("Parallel Fitacf3", |b| b.iter(|| rayon_fitacf3()));
}

fn fitacf3() {
    let file = File::open(TEST_FILE).expect("Test file not found");
    let rawacf = RawacfRecord::read_records(file).expect("Could not read records");
    let mut fitacf_records = vec![];

    let rec = &rawacf[0];
    let file_datetime = NaiveDate::from_ymd_opt(
        rec.get(&"time.yr".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.yr"),
        rec.get(&"time.mo".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.mo"),
        rec.get(&"time.dy".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.dy"),
    )
    .unwrap()
    .and_hms_opt(
        rec.get(&"time.hr".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.hr"),
        rec.get(&"time.mt".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.mt"),
        rec.get(&"time.sc".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.sc"),
    )
    .unwrap();
    let hdw = HdwInfo::new(
        rec.get(&"stid".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get stid"),
        file_datetime,
    )
    .expect("Unable to get hdw info");

    for rec in rawacf {
        fitacf_records.push(fit_rawacf_record(&rec, &hdw).expect("Could not fit record"));
    }
}

fn rayon_fitacf3() {
    let file = File::open(TEST_FILE).expect("Test file not found");
    let rawacf = RawacfRecord::read_records(file).expect("Could not read records");

    let rec = &rawacf[0];
    let file_datetime = NaiveDate::from_ymd_opt(
        rec.get(&"time.yr".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.yr"),
        rec.get(&"time.mo".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.mo"),
        rec.get(&"time.dy".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.dy"),
    )
    .unwrap()
    .and_hms_opt(
        rec.get(&"time.hr".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.hr"),
        rec.get(&"time.mt".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.mt"),
        rec.get(&"time.sc".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get time.sc"),
    )
    .unwrap();
    let hdw = HdwInfo::new(
        rec.get(&"stid".to_string())
            .unwrap()
            .clone()
            .try_into()
            .expect("Unable to get stid"),
        file_datetime,
    )
    .expect("Unable to get hdw info");

    let _ = rawacf
        .par_iter()
        .map(|rec| fit_rawacf_record(&rec, &hdw).expect("Could not fit record"));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
