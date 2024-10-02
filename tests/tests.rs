use assert_unordered::assert_eq_unordered;
use dmap::types::{DmapField, DmapVec};
use ndarray::Array;
use procdarn::fitting::fitacf3::fitacf_v3::fitacf3;
use std::iter::zip;
use procdarn::fitting::lmfit2::lmfit2::lmfit2;

#[test]
fn test_fitacf3() {
    // Create fitacf file from rawacf file
    let rawacf = dmap::read_rawacf("tests/test_files/test.rawacf".to_string().into())
        .expect("Could not read records");
    let fitacf_records = fitacf3(rawacf).expect("Unable to fit records");

    // Compare to fitacf file generated by RST
    let rst_records = dmap::read_fitacf("tests/test_files/test.fitacf".to_string().into())
        .expect("Could not read test.fitacf records");
    let variable_fields = vec!["origin.time", "origin.command"];
    for (test_rec, rst_rec) in zip(fitacf_records.iter(), rst_records.iter()) {
        assert_eq_unordered!(test_rec.keys(), rst_rec.keys());
        for k in test_rec.keys() {
            if variable_fields.contains(&&**k) {
            } else {
                match test_rec.get(k) {
                    Some(DmapField::Vector(DmapVec::Float(x))) => {
                        assert!(rst_rec.get(k).is_some(), "Testing {k}");
                        if let Some(DmapField::Vector(DmapVec::Float(y))) = rst_rec.get(k) {
                            assert!(
                                Array::<f32, _>::zeros(x.raw_dim())
                                    .abs_diff_eq(&((x - y) / x), 1e-5),
                                "Testing {k}\nleft: {x:?}\nright: {y:?}",
                            )
                        }
                    }
                    Some(_) => {
                        assert_eq!(test_rec.get(k), rst_rec.get(k), "Testing {k}")
                    }
                    None => {}
                }
            }
        }
    }
}

#[test]
fn test_lmfit2() {
    // Create fitacf file from rawacf file
    let rawacf = dmap::read_rawacf("tests/test_files/test.rawacf".to_string().into())
        .expect("Could not read records");
    let fitacf_records = lmfit2(rawacf).expect("Unable to fit records");

    // Compare to fitacf file generated by RST
    let rst_records = dmap::read_fitacf("tests/test_files/test.lmfit2".to_string().into())
        .expect("Could not read test.lmfit2 records");
    let variable_fields = vec!["origin.time", "origin.command", "p_s", "p_s_e", "w_s", "w_s_e", "sd_s", "phi0", "phi0_e", "elv", "elv_fitted", "elv_error", "x_sd_phi"];

    for (test_rec, rst_rec) in zip(fitacf_records.iter(), rst_records.iter()) {
        // assert_eq_unordered!(test_rec.keys(), rst_rec.keys());
        for k in test_rec.keys() {
            if variable_fields.contains(&&**k) {
            } else {
                match test_rec.get(k) {
                    Some(DmapField::Vector(DmapVec::Float(x))) => {
                        assert!(rst_rec.get(k).is_some(), "Testing {k}");
                        if let Some(DmapField::Vector(DmapVec::Float(y))) = rst_rec.get(k) {
                            assert!(
                                Array::<f32, _>::zeros(x.raw_dim())
                                    .abs_diff_eq(&((x - y) / x), 1e-4),
                                "Testing {k}\nleft: {x:?}\nright: {y:?}",
                            )
                        }
                    }
                    Some(_) => {
                        assert_eq!(test_rec.get(k), rst_rec.get(k), "Testing {k}")
                    }
                    None => {}
                }
            }
        }
    }
}