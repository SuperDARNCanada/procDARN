use std::fs::File;
use std::path::Path;
use backscatter_rs::RawDmapRead;

#[test]
fn read_iqdat() {
    let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    let _contents = RawDmapRead::new(file);
}