use approx::RelativeEq;
use assert_unordered::assert_eq_unordered;
use dmap::formats::grid::GridRecord;
use dmap::record::Record;
use dmap::types::{DmapField, DmapScalar, DmapVec};
use itertools::enumerate;
use procdarn::gridding::grid::{fit2grid, GridArgs};
use std::iter::zip;
use std::path::PathBuf;

fn compare_grid_recs(left_recs: Vec<GridRecord>, right_recs: Vec<GridRecord>) {
    let variable_fields = vec!["origin.time", "origin.command"];
    println!(
        "left_recs: {:?}\nright_recs: {:?}",
        left_recs[0], right_recs[0]
    );
    assert_eq!(left_recs.len(), right_recs.len());
    for (i, (test_rec, rst_rec)) in enumerate(zip(left_recs.iter(), right_recs.iter())) {
        assert_eq_unordered!(test_rec.keys(), rst_rec.keys());
        for k in test_rec.keys() {
            if variable_fields.contains(&&**k) {
            } else {
                eprintln!("testing rec {i} field {k}");
                match test_rec.get(k) {
                    Some(DmapField::Vector(DmapVec::Float(x))) => {
                        assert!(rst_rec.get(k).is_some(), "Testing rec {i} {k}");
                        if let Some(DmapField::Vector(DmapVec::Float(y))) = rst_rec.get(k) {
                            assert!(
                                x.map(|v| if v.is_nan() { -1_000_000.0 } else { *v })
                                    .relative_eq(
                                        &y.map(|v| if v.is_nan() { -1_000_000.0 } else { *v }),
                                        1e-4,
                                        1e-4
                                    ),
                                "Testing rec {i} {k}: left == right\n\tleft: {x}\n\tright: {y}\n\tDiff: {}",
                                (x - y) / x,
                            );
                        }
                    }
                    Some(DmapField::Scalar(DmapScalar::Float(x))) => {
                        assert!(rst_rec.get(k).is_some(), "Testing rec {i} {k}");
                        if let Some(DmapField::Scalar(DmapScalar::Float(y))) = rst_rec.get(k) {
                            assert!(x.relative_eq(y, 1e-5, 1e-5),
                                "Testing rec {i} {k}: left == right\n\nleft: {x}\n\nright: {y}\n\nDiff: {}",
                                (x - y) / x
                            );
                        }
                    }
                    Some(DmapField::Scalar(DmapScalar::Double(x))) => {
                        assert!(rst_rec.get(k).is_some(), "Testing rec {i} {k}");
                        if let Some(DmapField::Scalar(DmapScalar::Double(y))) = rst_rec.get(k) {
                            assert!(x.relative_eq(y, 1e-5, 1e-5),
                                    "Testing rec {i} {k}: left == right\n\nleft: {x}\n\nright: {y}\n\nDiff: {}",
                                    (x - y) / x
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

fn test_grid_with_args(args: &GridArgs, rst_args: Vec<String>) {
    // Create grid file
    let grid_recs = fit2grid(args).expect("Unable to make grid from fitacf");

    // Create grid file using same arguments in RST
    let output = std::process::Command::new("make_grid")
        .args(rst_args)
        .output()
        .expect("Failed to run make_grid");
    print!("{}", String::from_utf8_lossy(&output.stderr[..]));
    let rst_records =
        GridRecord::read_records(&output.stdout[..]).expect("Unable to read from stdout");

    // Compare the two new grid files
    compare_grid_recs(grid_recs, rst_records);
}

fn init_test_env() -> (GridArgs, Vec<String>) {
    let fitacf_files: Vec<PathBuf> = glob::glob("tests/test_files/fitacfs/*inv*")
        .expect("Could not find fitacf files")
        .map(|x| x.expect("Could not read fitacf file"))
        .collect();
    let fitacf_files = vec![fitacf_files[0].clone()];

    let args = GridArgs {
        outfile: "tests/test_files/out.grid".to_string().into(),
        infiles: fitacf_files.clone(),
        start_time: None,
        end_time: None,
        start_date: None,
        end_date: None,
        interval: None,
        scan_length: None,
        record_interval: 120,
        channel: None,
        channel_fix: None,
        exclude_beams: None,
        min_range_gate: None,
        max_range_gate: None,
        min_slant_range: None,
        max_slant_range: None,
        filter_weighting: 0,
        max_power: 60.0,
        max_velocity: 2500.0,
        max_spectral_width: 1000.0,
        max_velocity_error: 200.0,
        min_power: 3.0,
        min_velocity: 35.0,
        min_spectral_width: 10.0,
        min_velocity_error: 0.0,
        altitude: 300.0,
        max_frequency_var: 500000,
        boxcar_filter_flag: true,
        no_limits_flag: false,
        op_param_flag: false,
        exclude_neg_scan_flag: false,
        extended_mode_flag: false,
        sort_params_flag: false,
        ionosphere_only_flag: true,
        groundscatter_only_flag: false,
        all_data_flag: false,
        inertial_frame_flag: false,
        chisham_flag: false,
        verbose: false,
    };

    let mut rst_args = vec![];
    for infile in fitacf_files.iter() {
        rst_args.push(infile.display().to_string());
    }

    (args, rst_args)
}

#[test]
fn vanilla() {
    let (args, rst_args) = init_test_env();
    test_grid_with_args(&args, rst_args.clone());
}

#[test]
fn sort_params() {
    let (mut args, mut rst_args) = init_test_env();
    args.sort_params_flag = true;
    rst_args.insert(0, "-isort".to_string());
    test_grid_with_args(&args, rst_args.clone());
}

#[test]
fn groundscatter_only() {
    let (mut args, mut rst_args) = init_test_env();
    args.ionosphere_only_flag = false;
    args.groundscatter_only_flag = true;
    rst_args.insert(0, "-gs".to_string());
    test_grid_with_args(&args, rst_args.clone());
}

#[test]
fn ionosphere_only() {
    let (mut args, mut rst_args) = init_test_env();
    args.ionosphere_only_flag = true;
    rst_args.insert(0, "-ion".to_string());
    test_grid_with_args(&args, rst_args.clone());
}

#[test]
fn ion_and_gs() {
    let (mut args, mut rst_args) = init_test_env();
    args.ionosphere_only_flag = false;
    args.all_data_flag = true;
    rst_args.insert(0, "-both".to_string());
    test_grid_with_args(&args, rst_args.clone());
}

#[test]
fn no_boxcar() {
    let (mut args, mut rst_args) = init_test_env();
    args.boxcar_filter_flag = false;
    rst_args.insert(0, "-nav".to_string());
    test_grid_with_args(&args, rst_args.clone());
}

#[test]
fn no_boxcar_short_duration() {
    let (mut args, mut rst_args) = init_test_env();
    args.boxcar_filter_flag = false;
    args.end_time = Some("01:00".to_string());
    args.start_time = Some("00:55".to_string());
    rst_args.insert(0, "-et".to_string());
    rst_args.insert(1, "01:00".to_string());
    rst_args.insert(0, "-st".to_string());
    rst_args.insert(1, "00:55".to_string());
    rst_args.insert(0, "-nav".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn scan_length() {
    let (mut args, mut rst_args) = init_test_env();
    args.scan_length = Some(60);
    rst_args.insert(0, "-tl".to_string());
    rst_args.insert(1, "60".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn start_time() {
    let (mut args, mut rst_args) = init_test_env();
    args.start_time = Some("01:30".to_string());
    rst_args.insert(0, "-st".to_string());
    rst_args.insert(1, "01:30".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn end_time() {
    let (mut args, mut rst_args) = init_test_env();
    args.end_time = Some("00:30".to_string());
    rst_args.insert(0, "-et".to_string());
    rst_args.insert(1, "00:30".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn start_and_end_time() {
    let (mut args, mut rst_args) = init_test_env();
    args.end_time = Some("01:00".to_string());
    args.start_time = Some("00:55".to_string());
    rst_args.insert(0, "-et".to_string());
    rst_args.insert(1, "01:00".to_string());
    rst_args.insert(0, "-st".to_string());
    rst_args.insert(1, "00:55".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn start_end_and_chisham() {
    let (mut args, mut rst_args) = init_test_env();
    args.end_time = Some("01:00".to_string());
    args.start_time = Some("00:55".to_string());
    args.chisham_flag = true;
    rst_args.insert(0, "-et".to_string());
    rst_args.insert(1, "01:00".to_string());
    rst_args.insert(0, "-st".to_string());
    rst_args.insert(1, "00:55".to_string());
    rst_args.insert(0, "-chisham".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn start_end_and_range() {
    let (mut args, mut rst_args) = init_test_env();
    args.end_time = Some("01:00".to_string());
    args.start_time = Some("00:50".to_string());
    args.min_range_gate = Some(10);
    args.max_range_gate = Some(20);
    rst_args.insert(0, "-et".to_string());
    rst_args.insert(1, "01:00".to_string());
    rst_args.insert(0, "-st".to_string());
    rst_args.insert(1, "00:50".to_string());
    rst_args.insert(0, "-minrng".to_string());
    rst_args.insert(1, "10".to_string());
    rst_args.insert(0, "-maxrng".to_string());
    rst_args.insert(1, "20".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn exclude_beams() {
    let (mut args, mut rst_args) = init_test_env();
    args.exclude_beams = Some(vec![0, 7, 14]);
    rst_args.insert(0, "-ebm".to_string());
    rst_args.insert(1, "0,7,14".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn interval() {
    let (mut args, mut rst_args) = init_test_env();
    args.interval = Some("01:30".to_string());
    rst_args.insert(0, "-ex".to_string());
    rst_args.insert(1, "01:30".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn record_interval() {
    let (mut args, mut rst_args) = init_test_env();
    args.record_interval = 30;
    rst_args.insert(0, "-i".to_string());
    rst_args.insert(1, "30".to_string());
    rst_args.insert(2, "-vb".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn channel() {
    let (mut args, mut rst_args) = init_test_env();
    args.channel = Some('a');
    rst_args.insert(0, "-cn".to_string());
    rst_args.insert(1, "a".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn channel_fix() {
    let (mut args, mut rst_args) = init_test_env();
    args.channel_fix = Some('a');
    rst_args.insert(0, "-cn_fix".to_string());
    rst_args.insert(1, "a".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn min_range_gate() {
    let (mut args, mut rst_args) = init_test_env();
    args.min_range_gate = Some(10);
    rst_args.insert(0, "-minrng".to_string());
    rst_args.insert(1, "10".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn max_range_gate() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_range_gate = Some(20);
    rst_args.insert(0, "-maxrng".to_string());
    rst_args.insert(1, "20".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn min_slant_range() {
    let (mut args, mut rst_args) = init_test_env();
    args.min_slant_range = Some(500.0);
    rst_args.insert(0, "-minsrng".to_string());
    rst_args.insert(1, "500".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn max_slant_range() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_slant_range = Some(500.0);
    rst_args.insert(0, "-maxsrng".to_string());
    rst_args.insert(1, "500".to_string());
    test_grid_with_args(&args, rst_args);
}

#[test]
fn filter_weighting() {
    let (mut args, mut rst_args) = init_test_env();
    args.filter_weighting = 4;
    rst_args.insert(0, "-fwgt".to_string());
    rst_args.insert(1, "4".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn max_power() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_power = 30.0;
    rst_args.insert(0, "-pmax".to_string());
    rst_args.insert(1, "30".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn max_velocity() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_velocity = 1200.0;
    rst_args.insert(0, "-vmax".to_string());
    rst_args.insert(1, "1200".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn max_spectral_width() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_spectral_width = 500.0;
    rst_args.insert(0, "-wmax".to_string());
    rst_args.insert(1, "500".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn max_velocity_error() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_velocity_error = 250.0;
    rst_args.insert(0, "-vemax".to_string());
    rst_args.insert(1, "250".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn min_power() {
    let (mut args, mut rst_args) = init_test_env();
    args.min_power = 4.5;
    rst_args.insert(0, "-pmin".to_string());
    rst_args.insert(1, "4.5".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn min_velocity() {
    let (mut args, mut rst_args) = init_test_env();
    args.min_velocity = 50.0;
    rst_args.insert(0, "-vmin".to_string());
    rst_args.insert(1, "50".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn min_spectral_width() {
    let (mut args, mut rst_args) = init_test_env();
    args.min_spectral_width = 30.0;
    rst_args.insert(0, "-wmin".to_string());
    rst_args.insert(1, "30".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn min_velocity_error() {
    let (mut args, mut rst_args) = init_test_env();
    args.min_velocity_error = 15.0;
    rst_args.insert(0, "-vemin".to_string());
    rst_args.insert(1, "15".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn altitude() {
    let (mut args, mut rst_args) = init_test_env();
    args.altitude = 250.0;
    rst_args.insert(0, "-alt".to_string());
    rst_args.insert(1, "250".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn max_frequency_var() {
    let (mut args, mut rst_args) = init_test_env();
    args.max_frequency_var = 100000;
    rst_args.insert(0, "-fmax".to_string());
    rst_args.insert(1, "100000".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn no_limits_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.no_limits_flag = true;
    rst_args.insert(0, "-nlm".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn op_param_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.op_param_flag = true;
    rst_args.insert(0, "-nb".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn exclude_neg_scan_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.exclude_neg_scan_flag = true;
    rst_args.insert(0, "-ns".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn extended_mode_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.extended_mode_flag = true;
    rst_args.insert(0, "-xtd".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn sort_params_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.sort_params_flag = true;
    rst_args.insert(0, "-isort".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn all_data_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.all_data_flag = true;
    args.ionosphere_only_flag = false;
    rst_args.insert(0, "-both".to_string());
    test_grid_with_args(&args, rst_args);
}
#[test]
fn inertial_frame_flag() {
    let (mut args, mut rst_args) = init_test_env();
    args.inertial_frame_flag = true;
    rst_args.insert(0, "-inertial".to_string());
    test_grid_with_args(&args, rst_args);
}
