use dmap::{Record, RawacfRecord, write_rawacf};
use std::fs::File;
use std::path::PathBuf;

fn main() {
    let file =
        File::open("tests/test_files/20210607.1801.00.cly.a.rawacf").expect("Test file not found");
    // let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
    //     .expect("Test file not found");
    // let file = File::open(Path::new("tests/test_files/20110214.map"))
    //     .expect("Test file not found");
    let contents = RawacfRecord::read_records(file).unwrap();

    write_rawacf(contents, &PathBuf::from("tests/test_files/temp.rawacf")).unwrap();
}
