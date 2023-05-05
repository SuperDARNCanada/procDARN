use dmap::formats::{FitacfRecord, RawacfRecord, DmapRecord};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use backscatter_rs::fitting::fitacf3::fitacf_v3::fit_rawacf_record;
use dmap;
use std::fs::File;

fn criterion_benchmark(c: &mut Criterion) {
    // c.bench_function("Read IQDAT", |b| b.iter(|| read_iqdat()));
    // c.bench_function("Read MAP", |b| b.iter(|| read_map()));
    // c.bench_function("Read RAWACF", |b| b.iter(|| read_rawacf()));
    c.bench_function("Fitacf3", |b| b.iter(|| fitacf3()));

    // c.bench_function("Read and Parse RAWACF", |b| {
    //     b.iter(|| read_and_parse_rawacf())
    // });
    // c.bench_function("Read and Parse FITACF", |b| {
    //     b.iter(|| read_and_parse_fitacf())
    // });
    //
    // let records = read_iqdat();
    // c.bench_with_input(
    //     BenchmarkId::new("Write IQDAT", "IQDAT Records"),
    //     &records,
    //     |b, s| b.iter(|| write_iqdat(s)),
    // );
}

fn fitacf3() {
    let file =
        File::open("tests/test_files/20210607.1801.00.cly.a.rawacf").expect("Test file not found");
    let rawacf = RawacfRecord::read_records(file).expect("Could not read records");
    let mut fitacf_records = vec![];
    // fitacf_records.push(fit_rawacf_record(&rawacf[0]).expect("Could not fit rawacf record"));
    for rec in rawacf {
        fitacf_records.push(fit_rawacf_record(&rec).expect("Could not fit record"));
    }
    dmap::formats::to_file("tests/test_files/temp.fitacf", &fitacf_records);
}

// fn read_rawacf() -> Vec<RawDmapRecord> {
//     let file = File::open("tests/test_files/20210607.1801.00.cly.a.rawacf.mean")
//         .expect("Test file not found");
//     dmap::read_records(file).unwrap()
// }
//
// fn read_and_parse_rawacf() -> Vec<RawacfRecord> {
//     let file = File::open("tests/test_files/20210607.1801.00.cly.a.rawacf.mean")
//         .expect("Test file not found");
//     let recs = dmap::read_records(file).unwrap();
//     let mut rawacf_recs = vec![];
//     for rec in recs {
//         rawacf_recs.push(RawacfRecord::new(&rec).unwrap());
//     }
//     rawacf_recs
// }
//
// fn read_iqdat() -> Vec<RawDmapRecord> {
//     let file =
//         File::open("tests/test_files/20160316.1945.01.rkn.iqdat").expect("Test file not found");
//     dmap::read_records(file).unwrap()
// }
//
// fn write_iqdat(records: &Vec<RawDmapRecord>) {
//     let file =
//         File::open("tests/test_files/20160316.1945.01.rkn.iqdat").expect("Test file not found");
//     dmap::read_records(file).unwrap();
//     dmap::to_file("tests/test_files/test.iqdat", records).unwrap();
// }
//
// fn read_map() {
//     let file = File::open("tests/test_files/20110214.map").expect("Test file not found");
//     dmap::read_records(file).unwrap();
// }

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
