use criterion::{criterion_group, criterion_main, Criterion};
use dmap;
use procdarn::fitting::fitacf3::fitacf_v3::{fitacf3, par_fitacf3};

const TEST_FILE: &str = "tests/test_files/large.rawacf";

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Fitacf3", |b| b.iter(|| single_threaded()));
    c.bench_function("Parallel Fitacf3", |b| b.iter(|| multi_threaded()));
}

fn single_threaded() {
    let rawacf = dmap::read_rawacf(TEST_FILE.to_string().into()).expect("Could not read records");
    fitacf3(rawacf).expect("Error with single-threaded fitacf3");
}

fn multi_threaded() {
    let rawacf = dmap::read_rawacf(TEST_FILE.to_string().into()).expect("Could not read records");
    par_fitacf3(rawacf).expect("Error with multi-threaded fitacf3");
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
