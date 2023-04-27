use backscatter_rs::dmap;
use backscatter_rs::dmap::RawDmapRecord;
use backscatter_rs::formats::Rawacf;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::fs::File;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Read IQDAT", |b| b.iter(|| read_iqdat()));
    c.bench_function("Read MAP", |b| b.iter(|| read_map()));
    c.bench_function("Read RAWACF", |b| b.iter(|| read_rawacf()));
    c.bench_function("Read and Parse RAWACF", |b| b.iter(|| read_and_parse_rawacf()));

    let records = read_iqdat();
    c.bench_with_input(
        BenchmarkId::new("Write IQDAT", "IQDAT Records"),
        &records,
        |b, s| b.iter(|| write_iqdat(s)),
    );
}

fn read_rawacf() -> Vec<RawDmapRecord> {
    let file =
        File::open("tests/test_files/20210607.1801.00.cly.a.rawacf.mean").expect("Test file not found");
    dmap::read_records(file).unwrap()
}

fn read_and_parse_rawacf() -> Vec<Rawacf> {
    let file =
        File::open("tests/test_files/20210607.1801.00.cly.a.rawacf.mean").expect("Test file not found");
    let recs = dmap::read_records(file).unwrap();
    let mut rawacf_recs = vec![];
    for rec in recs {
        rawacf_recs.push(Rawacf::new(&rec).unwrap());
    }
    rawacf_recs
}

fn read_iqdat() -> Vec<RawDmapRecord> {
    let file =
        File::open("tests/test_files/20160316.1945.01.rkn.iqdat").expect("Test file not found");
    dmap::read_records(file).unwrap()
}

fn write_iqdat(records: &Vec<RawDmapRecord>) {
    let file =
        File::open("tests/test_files/20160316.1945.01.rkn.iqdat").expect("Test file not found");
    dmap::read_records(file).unwrap();
    dmap::to_file("tests/test_files/test.iqdat", records).unwrap();
}

fn read_map() {
    let file = File::open("tests/test_files/20110214.map").expect("Test file not found");
    dmap::read_records(file).unwrap();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
