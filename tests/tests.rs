use backscatter_rs::dmap;
use itertools::izip;
use std::fs::File;
use std::path::Path;

#[test]
fn read_iqdat() {
    let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    let _contents = dmap::read_records(file);
}

#[test]
fn test_read_write_iqdat() {
    let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    let contents = dmap::read_records(file).unwrap();

    dmap::to_file("tests/test_files/test.iqdat", contents.clone()).unwrap();
    let test_file = File::open("tests/test_files/test.iqdat").expect("Test file unwritten");
    let test_contents = dmap::read_records(test_file).unwrap();
    for (read_rec, written_rec) in izip!(contents.iter(), test_contents.iter()) {
        assert_eq!(read_rec, written_rec)
    }
}

#[test]
fn test_read_write_map() {
    let file = File::open(Path::new("tests/test_files/20110214.map")).expect("Test file not found");
    let contents = dmap::read_records(file).unwrap();

    dmap::to_file("tests/test_files/test.map", contents.clone()).unwrap();
    let test_file = File::open("tests/test_files/test.map").expect("Test file unwritten");
    let test_contents = dmap::read_records(test_file).unwrap();
    for (read_rec, written_rec) in izip!(contents.iter(), test_contents.iter()) {
        assert_eq!(read_rec, written_rec)
    }
}
