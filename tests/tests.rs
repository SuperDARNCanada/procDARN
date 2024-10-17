use assert_unordered::assert_eq_unordered;
use dmap::types::{DmapField, DmapVec};
use itertools::enumerate;
use procdarn::fitting::fitacf3::fitacf_v3::fitacf3;
use std::iter::zip;

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
    for (i, (test_rec, rst_rec)) in enumerate(zip(fitacf_records.iter(), rst_records.iter())) {
        assert_eq_unordered!(test_rec.keys(), rst_rec.keys());
        for k in test_rec.keys() {
            if variable_fields.contains(&&**k) {
            } else {
                match test_rec.get(k) {
                    Some(DmapField::Vector(DmapVec::Float(x))) => {
                        assert!(rst_rec.get(k).is_some(), "Testing rec {i} {k}");
                        if let Some(DmapField::Vector(DmapVec::Float(y))) = rst_rec.get(k) {
                            assert!(
                                x.map(|v| if v.is_nan() { -1_000_000.0 } else { *v })
                                    .relative_eq(
                                        &y.map(|v| if v.is_nan() { -1_000_000.0 } else { *v }),
                                        1e-5,
                                        1e-5
                                    ),
                                "Testing rec {i} {k}: left == right\nleft: {x}\nright: {y}"
                            );
                        }
                    }
                    Some(_) => {
                        assert_eq!(test_rec.get(k), rst_rec.get(k), "Testing rec {i} {k}")
                    }
                    None => {}
                }
            }
        }
    }
}
