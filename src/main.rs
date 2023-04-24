use backscatter_rs::RawDmapRead;
use std::fs::File;
use std::path::Path;

fn main() {
    let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    println!("{:?}", RawDmapRead::new(file).unwrap().dmap_records[0]);
}
