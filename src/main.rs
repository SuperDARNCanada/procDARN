use backscatter_rs;
use backscatter_rs::DmapError;
use std::fs::File;
use std::path::Path;
use std::io::Read;

fn main() {
    let file = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    // let file = File::open(Path::new("tests/test_files/20110214.map"))
    //     .expect("Test file not found");
    let contents = backscatter_rs::read_records(file).unwrap();
    // println!("num records: {}", contents.len());
    // println!("  num scalars: {:?}", contents[1].num_scalars);
    // for s in contents[1].scalars.iter() {
    //     println!("  {}: {}", s.name, s.data);
    // };
    // println!("  num arrays: {:?}", contents[1].num_arrays);
    // for s in contents[1].arrays.iter() {
    //     println!("  {:?}: {:?}", s.name, s.data);
    // };

    backscatter_rs::to_file("tests/test_files/test.iqdat", contents.clone()).unwrap();
    // let test_file = File::open("tests/test_files/test.iqdat").expect("Test file unwritten");
    // let _test_contents = backscatter_rs::read_records(test_file).unwrap();
    //
    //
    let mut bytes = vec![];
    for rec in contents {
        let mut rec_bytes = rec.to_bytes();
        println!("rec: {} bytes", rec_bytes.len());
        bytes.append(&mut rec_bytes);
    }

    let mut buffer: Vec<u8> = vec![];
    let mut file1 = File::open(Path::new("tests/test_files/20160316.1945.01.rkn.iqdat"))
        .expect("Test file not found");
    // let mut file1 = File::open(Path::new("tests/test_files/20110214.map"))
    //     .expect("Test file not found");
    file1
        .read_to_end(&mut buffer)
        .map_err(|_| DmapError::Message("Could not read data".to_string())).unwrap();

    println!("{}", bytes.len());
    println!("{}", buffer.len());
    for i in 0..50 {
        println!("{:?} -> {:?}", buffer[i], bytes[i]);
    }

}
