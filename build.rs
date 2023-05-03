use git2::Repository;
use std::env;

fn main() {
    // clone the hdw repo
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let url = "https://github.com/SuperDARN/hdw";
    match Repository::open(format!("{:?}/hdw", out_dir)) {
        Ok(repo) => repo,
        Err(_) => match Repository::clone(url, format!("{:?}/hdw", out_dir)) {
            Ok(r) => r,
            Err(err) => panic!("failed to clone: {}", err),
        },
    };
}
