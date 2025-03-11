use git2::Repository;
use std::path::Path;

fn main() {
    // clone the utils repo
    let out_dir = "target/hdw/";
    let url = "https://github.com/SuperDARN/hdw";
    println!("Installing {url} to {out_dir}");
    if Path::new(&out_dir).is_dir() {
        println!("{out_dir} already exists");
    } else {
        match Repository::clone(url, out_dir) {
            Ok(r) => r,
            Err(err) => panic!("failed to clone: {err}"),
        };
    }

    for (key, value) in std::env::vars() {
        println!("{key}: {value}");
    }
    let aacgm_coeffs_path = std::env::var("DEP_AACGMV2_COEFFS_DAT_PREFIX").expect("DEP_AACGMV2_RS_COEFFS_DAT_PREFIX missing");
    let igrf_coeffs_path = std::env::var("DEP_AACGMV2_IGRF_COEFFS").expect("DEP_AACGMV2_RS_IGRF_COEFFS");

    println!("cargo:rustc-env=AACGM_v2_DAT_PREFIX={aacgm_coeffs_path}");
    println!("cargo:rustc-env=IGRF_COEFFS={igrf_coeffs_path}");
}
