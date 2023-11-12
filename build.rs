use git2::Repository;
use std::path::Path;

fn main() {
    // clone the utils repo
    let out_dir = "target/hdw/";
    let url = "https://github.com/SuperDARN/hdw";
    println!("Installing {url} to {out_dir}");
    if !Path::new(&out_dir).is_dir() {
        match Repository::clone(url, out_dir) {
            Ok(r) => r,
            Err(err) => panic!("failed to clone: {}", err),
        };
    } else {
        println!("{out_dir} already exists");
    }
}
