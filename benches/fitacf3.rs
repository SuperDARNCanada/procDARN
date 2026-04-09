use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use dmap::{RawacfRecord, Record};
use procdarn::fitting::fitacf3::fitacf_v3::{fitacf3, par_fitacf3};

const TEST_FILE: &str = "tests/test_files/large.rawacf";
// const TEST_FILE: &str = "/data/dmap_files/20221107.2200.00.rkn.a.rawacf"; // widebeam file

fn criterion_benchmark(c: &mut Criterion) {
    let rawacf = RawacfRecord::read_file(TEST_FILE.to_string()).expect("Could not read records");

    c.bench_function("Fitacf3", |b| {
        b.iter_batched(
            || rawacf.clone(),
            |rawacf| fitacf3(rawacf),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("Parallel Fitacf3", |b| {
        b.iter_batched(
            || rawacf.clone(),
            |rawacf| par_fitacf3(rawacf),
            BatchSize::SmallInput,
        )
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark
}
criterion_main!(benches);
