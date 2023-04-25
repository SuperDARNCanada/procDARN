use criterion::{BenchmarkId, criterion_group, criterion_main, Criterion};
use backscatter_rs::dmap;
use std::fs::File;
use backscatter_rs::dmap::RawDmapRecord;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Read IQDAT", |b| b.iter(|| read_iqdat()));
    c.bench_function("Read MAP", |b| b.iter(|| read_map()));

    let records = read_iqdat();
    c.bench_with_input(BenchmarkId::new("Write IQDAT", "IQDAT Records"), &records,
                       |b, s| b.iter(|| write_iqdat(s)));
}

fn read_iqdat() -> Vec<RawDmapRecord> {
    let file = File::open("tests/test_files/20160316.1945.01.rkn.iqdat")
        .expect("Test file not found");
    dmap::read_records(file).unwrap()
}

fn write_iqdat(records: &Vec<RawDmapRecord>) {
    let file = File::open("tests/test_files/20160316.1945.01.rkn.iqdat")
        .expect("Test file not found");
    dmap::read_records(file).unwrap();
    dmap::to_file("tests/test_files/test.iqdat", records).unwrap();
}

fn read_map() {
    let file = File::open("tests/test_files/20110214.map").expect("Test file not found");
    dmap::read_records(file).unwrap();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);



