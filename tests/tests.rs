use assert_unordered::assert_eq_unordered;
use backscatter_rs::fitting::fitacf3::fitacf_v3::fit_rawacf_record;
use backscatter_rs::utils::rawacf::get_hdw;
use dmap::types::{DmapField, DmapVec};
use ndarray::Array;
use std::iter::zip;

#[test]
fn test_fitacf3() {
    // Create fitacf file from rawacf file
    let rawacf = dmap::read_rawacf("tests/test_files/test.rawacf".to_string().into())
        .expect("Could not read records");
    let hdw = get_hdw(&rawacf[0]).expect("Unable to get hdw info");

    // Fit the rawacf data
    let mut fitacf_records = vec![];
    for rec in rawacf {
        fitacf_records.push(fit_rawacf_record(&rec, &hdw).expect("Could not fit record"));
    }

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
