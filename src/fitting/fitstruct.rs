// struct RangeNode {
//     range: i32,
//     cross_range_interference: &f64,
//     refractive_idx: f32,
//     llist alpha_2,
//     llist phases,
//     llist powers,
//     llist elev,
//     lin_pwr_fit: &FitData,
//     quad_pwr_fit: &FitData,
//     lin_pwr_fit_err: &FitData,
//     quad_pwr_fit_err: &FitData,
//     phase_fit: &FitData,
//     elev_fit: &FitData,
// }

struct PhaseNode {
    phi: f64,
    t: f64,
    std_dev: f64,
    lag_idx: i32, // TODO: Is this redundant with Alpha?
    alpha_2: f64, // TODO: Is this redundant with Alpha?
}

struct PowerNode {
    ln_power: f64,
    t: f64,
    std_dev: f64,
    lag_idx: i32, // TODO: Is this redundant with Alpha?
    alpha_2: f64, // TODO: Is this redundant with Alpha?
}

struct LagNode {
    lag_num: i32,
    pulses: [i32; 2],
    lag_idx: i32,
    sample_base_1: i32,
    sample_base_2: i32,
}

struct Alpha {
    lag_idx: i32,
    alpha_2: f64,
}

enum IntfPosition {
    Forward,
    Behind,
}

struct FitData<'a> {
    channel: i32,
    offset: i32, /* used for stereo badlags */
    cp: i32,
    xcf_flag: i32,
    transmit_freq: i32,
    noise: f32,
    num_ranges: i32,
    sample_sep: i32,
    num_avg: i32,
    num_lags: i32,
    multi_pulse_increment: i32,
    tx_pulse_len: i32,
    lag_to_first_range: i32,
    num_pulses: i32,
    beam_num: i32,
    old: i32,
    lag: &'a [i32; 2],
    pulse: &'a i32,
    pwr0: &'a f64,
    acfd: Vec<f64>,
    xcfd: Vec<f64>,
    maxbeam: i32,
    beam_offset: f64,
    beam_sep: f64,
    interferometer_offset: [f64; 3],
    phi_diff: f64,
    time_diff: f64,
    vel_dir: f64,
    time: Time,
}

struct Time {
    yr: i16,
    mo: i16,
    dy: i16,
    hr: i16,
    mt: i16,
    sc: i16,
    us: i32,
}

struct FittedData<'a> {
    delta: f64,
    a: f64,
    b: f64,
    variance_a: f64,
    variance_b: f64,
    delta_a: f64,
    delta_b: f64,
    covariance_ab: f64,
    residual_ab: f64,
    quality: f64,
    chi_squared: f64,
    sums: &'a Sums,
}

struct Sums {
    sum: f64,
    sum_x: f64,
    sum_y: f64,
    sum_x_squared: f64,
    sum_y_squared: f64,
}

enum FitType {
    Linear,
    Quadratic,
}
