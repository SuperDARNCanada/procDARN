use backscatter_rs::fitting::fitacf3::fitacf_v3::fit_rawacf_record;
use backscatter_rs::utils::rawacf::get_hdw;
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
    let hdw = get_hdw(&rawacf[0]).expect("Unable to get hdw info");

    let mut fitacf_records = vec![];
    for rec in rawacf {
        fitacf_records.push(fit_rawacf_record(&rec, &hdw).expect("Could not fit record"));
    }
}

fn rayon_fitacf3() {
    let file = File::open(TEST_FILE).expect("Test file not found");
    let rawacf = RawacfRecord::read_records(file).expect("Could not read records");
    let hdw = get_hdw(&rawacf[0]).expect("Unable to get hdw info");

    let _ = rawacf
        .par_iter()
        .map(|rec| fit_rawacf_record(&rec, &hdw).expect("Could not fit record"));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
