use std::env;
use std::fs;
use std::path::Path;
use git2::Repository;

fn main() {
    // clone the hdw repo
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let url = "https://github.com/SuperDARN/hdw";
    let repo = match Repository::clone(url, out_dir) {
        Ok(repo) => repo,
        Err(e) => panic!("failed to clone: {}", e),
    };
}