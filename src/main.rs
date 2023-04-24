use backscatter_rs;
use backscatter_rs::DmapError;
use std::fs::File;
use std::io::Read;
use std::path::Path;

fn main() {
    let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    // let file = File::open(Path::new("tests/test_files/20110214.map"))
    //     .expect("Test file not found");
    let contents = backscatter_rs::read_records(file).unwrap();

    backscatter_rs::to_file("tests/test_files/test.iqdat", contents.clone()).unwrap();
}
