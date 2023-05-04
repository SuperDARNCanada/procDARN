use git2::Repository;
use std::env;
use std::path::Path;

fn main() {
    // clone the hdw repo
    let out_dir = env::var("HDW_DIR").unwrap();
    let url = "https://github.com/SuperDARN/hdw";
    if !Path::new(&out_dir).is_dir() {
        match Repository::clone(url, out_dir) {
            Ok(r) => r,
            Err(err) => panic!("failed to clone: {}", err),
        };
    }
}
