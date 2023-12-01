use crate::fitting::fitacf3::fitacf_v3::Fitacf3Error;
use crate::fitting::fitacf3::fitstruct::RangeNode;
use crate::utils::dmap::convert_to_dmapvec;
use crate::utils::hdw::HdwInfo;
use dmap::formats::{FitacfRecord, RawacfRecord};
use dmap::DmapVec;
use std::f32::consts::PI as PI_f32;
use std::iter::zip;

pub const FITACF_REVISION_MAJOR: i32 = 3;
pub const FITACF_REVISION_MINOR: i32 = 0;
pub const V_MAX: f32 = 30.0;
pub const W_MAX: f32 = 90.0;

pub fn determinations(
    rec: &RawacfRecord,
    ranges: Vec<RangeNode>,
    noise_power: f32,
    hdw: &HdwInfo,
) -> Result<FitacfRecord, Fitacf3Error> {
    let range_list: Vec<i16> = ranges.iter().map(|r| r.range_num as i16).collect();
    let lag_0_power_db: Vec<f32> = rec
        .lag_zero_power
        .data
        .iter()
        .map(|p| {
            if p - noise_power > 0.0 {
                10.0 * ((p - noise_power) / noise_power).log10()
            } else {
                -50.0
            }
        })
        .collect();
    if range_list.is_empty() {
        Ok(FitacfRecord {
            radar_revision_major: rec.radar_revision_major,
            radar_revision_minor: rec.radar_revision_minor,
            origin_code: rec.origin_code,
            origin_time: "".to_string(),    // TODO: Get current time
            origin_command: "".to_string(), // TODO: Get this
            control_program: rec.control_program,
            station_id: rec.station_id,
            year: rec.year,
            month: rec.month,
            day: rec.day,
            hour: rec.hour,
            minute: rec.minute,
            second: rec.second,
            microsecond: rec.microsecond,
            tx_power: rec.tx_power,
            num_averages: rec.num_averages,
            attenuation: rec.attenuation,
            lag_to_first_range: rec.lag_to_first_range,
            sample_separation: rec.sample_separation,
            error_code: rec.error_code,
            agc_status: rec.agc_status,
            low_power_status: rec.low_power_status,
            search_noise: rec.search_noise,
            mean_noise: rec.mean_noise,
            channel: rec.channel,
            beam_num: rec.beam_num,
            beam_azimuth: rec.beam_azimuth,
            scan_flag: rec.scan_flag,
            offset: rec.offset,
            rx_rise_time: rec.rx_rise_time,
            intt_second: rec.intt_second,
            intt_microsecond: rec.intt_microsecond,
            tx_pulse_length: rec.tx_pulse_length,
            multi_pulse_increment: rec.multi_pulse_increment,
            num_pulses: rec.num_pulses,
            num_lags: rec.num_lags,
            num_lags_extras: rec.num_lags_extras,
            if_mode: rec.if_mode,
            num_ranges: rec.num_ranges,
            first_range: rec.first_range,
            range_sep: rec.range_sep,
            xcf_flag: rec.xcf_flag,
            tx_freq: rec.tx_freq,
            max_power: rec.max_power,
            max_noise_level: rec.max_noise_level,
            comment: rec.comment.clone(),
            algorithm: None,
            fitacf_revision_major: FITACF_REVISION_MAJOR,
            fitacf_revision_minor: FITACF_REVISION_MINOR,
            sky_noise: noise_power,
            lag_zero_noise: 0.0,
            velocity_noise: 0.0,
            tdiff: None,
            pulse_table: rec.pulse_table.clone(),
            lag_table: rec.lag_table.clone(),
            lag_zero_power: convert_to_dmapvec(lag_0_power_db),
            range_list: convert_to_dmapvec(vec![]),
            fitted_points: convert_to_dmapvec(vec![]),
            quality_flag: convert_to_dmapvec(vec![]),
            ground_flag: convert_to_dmapvec(vec![]),
            lambda_power: convert_to_dmapvec(vec![]),
            lambda_power_error: convert_to_dmapvec(vec![]),
            sigma_power: convert_to_dmapvec(vec![]),
            sigma_power_error: convert_to_dmapvec(vec![]),
            velocity: convert_to_dmapvec(vec![]),
            velocity_error: convert_to_dmapvec(vec![]),
            lambda_spectral_width: convert_to_dmapvec(vec![]),
            lambda_spectral_width_error: convert_to_dmapvec(vec![]),
            sigma_spectral_width: convert_to_dmapvec(vec![]),
            sigma_spectral_width_error: convert_to_dmapvec(vec![]),
            lambda_std_dev: convert_to_dmapvec(vec![]),
            sigma_std_dev: convert_to_dmapvec(vec![]),
            phi_std_dev: convert_to_dmapvec(vec![]),
            xcf_quality_flag: None,
            xcf_ground_flag: None,
            lambda_xcf_power: None,
            lambda_xcf_power_error: None,
            sigma_xcf_power: None,
            sigma_xcf_power_error: None,
            xcf_velocity: None,
            xcf_velocity_error: None,
            lambda_xcf_spectral_width: None,
            lambda_xcf_spectral_width_error: None,
            sigma_xcf_spectral_width: None,
            sigma_xcf_spectral_width_error: None,
            lag_zero_phi: None,
            lag_zero_phi_error: None,
            elevation: None,
            elevation_fitted: None,
            elevation_error: None,
            elevation_low: None,
            elevation_high: None,
            lambda_xcf_std_dev: None,
            sigma_xcf_std_dev: None,
            phi_xcf_std_dev: None,
        })
    } else {
        let num_lags: Vec<i16> = ranges
            .iter()
            .map(|r| r.powers.ln_power.len() as i16)
            .collect();
        let quality_flag: Vec<i8> = range_list.iter().map(|_| 1).collect();
        let noise_db: f32 = 10.0 * noise_power.log10();
        let power_linear: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * r
                    .lin_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf without linear fitted power")
                    .intercept as f32
                    / (10.0_f32).ln()
                    - noise_db
            })
            .collect();
        let power_linear_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * (r
                    .lin_pwr_fit_err
                    .as_ref()
                    .expect("Unable to make fitacf without linear fitted power error")
                    .variance_intercept as f32)
                    .sqrt()
                    / (10.0_f32).ln()
            })
            .collect();
        let power_quadratic: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * (r
                    .quad_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf without quadratic fitted power")
                    .intercept as f32)
                    / (10.0_f32).ln()
                    - noise_db
            })
            .collect();
        let power_quadratic_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                10.0 * (r
                    .quad_pwr_fit_err
                    .as_ref()
                    .expect("Unable to make fitacf without quadratic fitted power error")
                    .variance_intercept as f32)
                    .sqrt()
                    / (10.0_f32).ln()
            })
            .collect();
        let velocity_conversion: f32 =
            299792458.0 * hdw.velocity_sign / (4.0 * PI_f32 * rec.tx_freq as f32 * 1000.0);
        let velocity: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.phase_fit
                    .as_ref()
                    .expect("Unable to make fitacf without fitted velocity")
                    .slope as f32)
                    * velocity_conversion
            })
            .collect();
        let velocity_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.phase_fit
                    .as_ref()
                    .expect("Unable to make fitacf without fitted velocity")
                    .variance_slope as f32)
                    .sqrt()
                    * velocity_conversion
            })
            .collect();
        let width_conversion: f32 =
            299792458.0 * 2.0 / (4.0 * PI_f32 * rec.tx_freq as f32 * 1000.0);
        let spectral_width_linear: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.lin_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf spectral width without fitted power")
                    .slope as f32)
                    .abs()
                    * width_conversion
            })
            .collect();
        let spectral_width_linear_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.lin_pwr_fit_err
                    .as_ref()
                    .expect("Unable to make fitacf spectral width error without fitted power error")
                    .variance_slope as f32)
                    .sqrt()
                    * width_conversion
            })
            .collect();
        let quadratic_width_conversion: f32 =
            299792458.0 * (2.0_f32).ln().sqrt() / (PI_f32 * rec.tx_freq as f32 * 1000.0);
        let spectral_width_quadratic: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.quad_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf quadratic spectral width without fitted power")
                    .slope as f32)
                    .abs()
                    .sqrt()
                    * quadratic_width_conversion
            })
            .collect();
        let spectral_width_quadratic_error: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.quad_pwr_fit_err.as_ref().expect("Unable to make fitacf quadratic spectral width error without fitted power error")
                    .variance_slope as f32).sqrt() * quadratic_width_conversion /
                    ((r.quad_pwr_fit.as_ref().expect("Unable to make fitacf quadratic spectral width error without fitted power error")
                        .slope as f32).abs().sqrt() * 2.0)
            })
            .collect();
        let std_dev_linear: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.lin_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf linear std deviation without fitted power")
                    .chi_squared as f32
            })
            .collect();
        let std_dev_quadratic: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.quad_pwr_fit
                    .as_ref()
                    .expect("Unable to make fitacf quadratic std deviation without fitted power")
                    .chi_squared as f32
            })
            .collect();
        let std_dev_phi: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.phase_fit
                    .as_ref()
                    .expect("Unable to make fitacf phi std deviation")
                    .chi_squared as f32
            })
            .collect();
        let groundscatter_flag: Vec<i8> = zip(velocity.iter(), spectral_width_linear.iter())
            .map(|(v, w)| (v.abs() - (V_MAX - w * (V_MAX / W_MAX)) < 1.0) as i8)
            .collect();
        let xcfs = &rec
            .xcfs
            .as_ref()
            .expect("Unable to make fitacf xcf_phi0")
            .data;
        let xcf_phi0: Vec<f32> = ranges
            .iter()
            .map(|r| {
                xcfs[r.range_idx * rec.num_lags as usize * 2 + 1]
                    .atan2(xcfs[r.range_idx * rec.num_lags as usize * 2])
                    * hdw.phase_sign
            })
            .collect();
        let xcf_phi0_err: Vec<f32> = ranges
            .iter()
            .map(|r| {
                (r.elev_fit
                    .as_ref()
                    .expect("Unable to make fitacf xcf_phi0_err")
                    .variance_intercept as f32)
                    .sqrt()
            })
            .collect();
        let xcf_phi_std_dev: Vec<f32> = ranges
            .iter()
            .map(|r| {
                r.elev_fit
                    .as_ref()
                    .expect("Unable to make fitacf xcf_phi_std_dev")
                    .chi_squared as f32
            })
            .collect();
        let (elevation_low, elevation_normal, elevation_high) =
            calculate_elevation(&ranges, rec, &xcf_phi0, hdw);

        let float_zeros = DmapVec {
            data: quality_flag.iter().map(|_| 0.0_f32).collect(),
            dimensions: vec![quality_flag.len() as i32],
        };
        let i8_zeros = DmapVec {
            data: quality_flag.iter().map(|_| 0_i8).collect(),
            dimensions: vec![quality_flag.len() as i32],
        };

        Ok(FitacfRecord {
            radar_revision_major: rec.radar_revision_major,
            radar_revision_minor: rec.radar_revision_minor,
            origin_code: rec.origin_code,
            origin_time: "".to_string(),    // TODO: Get current time
            origin_command: "".to_string(), // TODO: Get this
            control_program: rec.control_program,
            station_id: rec.station_id,
            year: rec.year,
            month: rec.month,
            day: rec.day,
            hour: rec.hour,
            minute: rec.minute,
            second: rec.second,
            microsecond: rec.microsecond,
            tx_power: rec.tx_power,
            num_averages: rec.num_averages,
            attenuation: rec.attenuation,
            lag_to_first_range: rec.lag_to_first_range,
            sample_separation: rec.sample_separation,
            error_code: rec.error_code,
            agc_status: rec.agc_status,
            low_power_status: rec.low_power_status,
            search_noise: rec.search_noise,
            mean_noise: rec.mean_noise,
            channel: rec.channel,
            beam_num: rec.beam_num,
            beam_azimuth: rec.beam_azimuth,
            scan_flag: rec.scan_flag,
            offset: rec.offset,
            rx_rise_time: rec.rx_rise_time,
            intt_second: rec.intt_second,
            intt_microsecond: rec.intt_microsecond,
            tx_pulse_length: rec.tx_pulse_length,
            multi_pulse_increment: rec.multi_pulse_increment,
            num_pulses: rec.num_pulses,
            num_lags: rec.num_lags,
            num_lags_extras: rec.num_lags_extras,
            if_mode: rec.if_mode,
            num_ranges: rec.num_ranges,
            first_range: rec.first_range,
            range_sep: rec.range_sep,
            xcf_flag: rec.xcf_flag,
            tx_freq: rec.tx_freq,
            max_power: rec.max_power,
            max_noise_level: rec.max_noise_level,
            comment: rec.comment.clone(),
            algorithm: None,
            fitacf_revision_major: FITACF_REVISION_MAJOR,
            fitacf_revision_minor: FITACF_REVISION_MINOR,
            sky_noise: noise_power,
            lag_zero_noise: 0.0,
            velocity_noise: 0.0,
            tdiff: None,
            pulse_table: rec.pulse_table.clone(),
            lag_table: rec.lag_table.clone(),
            lag_zero_power: convert_to_dmapvec(lag_0_power_db),
            range_list: convert_to_dmapvec(range_list),
            fitted_points: convert_to_dmapvec(num_lags),
            quality_flag: convert_to_dmapvec(quality_flag),
            ground_flag: convert_to_dmapvec(groundscatter_flag),
            lambda_power: convert_to_dmapvec(power_linear),
            lambda_power_error: convert_to_dmapvec(power_linear_error),
            sigma_power: convert_to_dmapvec(power_quadratic),
            sigma_power_error: convert_to_dmapvec(power_quadratic_error),
            velocity: convert_to_dmapvec(velocity),
            velocity_error: convert_to_dmapvec(velocity_error),
            lambda_spectral_width: convert_to_dmapvec(spectral_width_linear),
            lambda_spectral_width_error: convert_to_dmapvec(spectral_width_linear_error),
            sigma_spectral_width: convert_to_dmapvec(spectral_width_quadratic),
            sigma_spectral_width_error: convert_to_dmapvec(spectral_width_quadratic_error),
            lambda_std_dev: convert_to_dmapvec(std_dev_linear),
            sigma_std_dev: convert_to_dmapvec(std_dev_quadratic),
            phi_std_dev: convert_to_dmapvec(std_dev_phi),
            xcf_quality_flag: Some(i8_zeros.clone()),
            xcf_ground_flag: Some(i8_zeros),
            lambda_xcf_power: Some(float_zeros.clone()),
            lambda_xcf_power_error: Some(float_zeros.clone()),
            sigma_xcf_power: Some(float_zeros.clone()),
            sigma_xcf_power_error: Some(float_zeros.clone()),
            xcf_velocity: Some(float_zeros.clone()),
            xcf_velocity_error: Some(float_zeros.clone()),
            lambda_xcf_spectral_width: Some(float_zeros.clone()),
            lambda_xcf_spectral_width_error: Some(float_zeros.clone()),
            sigma_xcf_spectral_width: Some(float_zeros.clone()),
            sigma_xcf_spectral_width_error: Some(float_zeros.clone()),
            lag_zero_phi: Some(convert_to_dmapvec(xcf_phi0)),
            lag_zero_phi_error: Some(convert_to_dmapvec(xcf_phi0_err)),
            elevation: Some(convert_to_dmapvec(elevation_normal)),
            elevation_fitted: None,
            elevation_error: None,
            elevation_low: Some(convert_to_dmapvec(elevation_low)),
            elevation_high: Some(convert_to_dmapvec(elevation_high)),
            lambda_xcf_std_dev: Some(float_zeros.clone()),
            sigma_xcf_std_dev: Some(float_zeros),
            phi_xcf_std_dev: Some(convert_to_dmapvec(xcf_phi_std_dev)),
        })
    }
}

fn calculate_elevation(
    ranges: &[RangeNode],
    rec: &RawacfRecord,
    xcf_phi0: &[f32],
    hdw: &HdwInfo,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let x = hdw.intf_offset_x;
    let y = hdw.intf_offset_y;
    let z = hdw.intf_offset_z;

    let array_separation: f32 = (x * x + y * y + z * z).sqrt();
    let mut elevation_corr = (z / array_separation).asin();
    let phi_sign: f32;
    if y > 0.0 {
        phi_sign = 1.0;
    } else {
        phi_sign = -1.0;
        elevation_corr *= -1.0;
    }
    let azimuth_offset = hdw.max_num_beams as f32 / 2.0 - 0.5;
    let phi_0 =
        (hdw.beam_separation * (rec.beam_num as f32 - azimuth_offset) * PI_f32 / 180.0).cos();
    let wave_num = 2.0 * PI_f32 * rec.tx_freq as f32 * 1000.0 / 299792458.0;
    let cable_offset = -2.0 * PI_f32 * rec.tx_freq as f32 * 1000.0 * hdw.tdiff_a * 1.0e-6;
    let phase_diff_max = phi_sign * wave_num * array_separation * phi_0 + cable_offset;
    let mut psi: Vec<f32> = ranges
        .iter()
        .map(|r| {
            let x = r
                .elev_fit
                .as_ref()
                .expect("Unable to find elevation without fitted elevation")
                .intercept as f32;
            let mut y =
                x + 2.0 * PI_f32 * ((phase_diff_max - x) / (2.0 * PI_f32)).floor() - cable_offset;
            if phi_sign < 0.0 {
                y += 2.0 * PI_f32;
            }
            y
        })
        .collect();
    let mut psi_kd: Vec<f32> = psi
        .iter()
        .map(|p| p / (wave_num * array_separation))
        .collect();
    let mut theta: Vec<f32> = psi_kd.iter().map(|p| phi_0 * phi_0 - p * p).collect();
    let elevation: Vec<f32> = theta
        .iter()
        .map(|&t| {
            if t < 0.0 || t.abs() > 1.0 {
                -elevation_corr
            } else {
                t.sqrt().asin()
            }
        })
        .collect();
    let elevation_high: Vec<f32> = elevation
        .iter()
        .map(|e| (e + elevation_corr) * 180.0 / PI_f32)
        .collect();
    let psi_k2d2: Vec<f32> = psi
        .iter()
        .map(|p| p / (wave_num * wave_num * array_separation * array_separation))
        .collect();
    let df_by_dy: Vec<f32> = zip(psi_k2d2.iter(), theta.iter())
        .map(|(p, t)| p / (t * (1.0 - t)).sqrt())
        .collect();
    let errors: Vec<f32> = ranges
        .iter()
        .map(|r| {
            r.elev_fit
                .as_ref()
                .expect("Unable to calculate elevation errors")
                .variance_intercept as f32
        })
        .collect();
    let elevations_low: Vec<f32> = zip(errors.iter(), df_by_dy.iter())
        .map(|(e, d)| e.sqrt() * d.abs() * 180.0 / PI_f32)
        .collect();

    // This time, use the xcf lag0 phase
    psi = xcf_phi0
        .iter()
        .map(|&x| {
            let mut y = x as f32
                + 2.0 * PI_f32 * ((phase_diff_max - x as f32) / (2.0 * PI_f32)).floor()
                - cable_offset;
            if phi_sign < 0.0 {
                y += 2.0 * PI_f32;
            }
            y
        })
        .collect();
    psi_kd = psi
        .iter()
        .map(|p| p / (wave_num * array_separation))
        .collect();
    theta = psi_kd.iter().map(|p| phi_0 * phi_0 - p * p).collect();
    let elevation_normal: Vec<f32> = theta
        .iter()
        .map(|&t| {
            if t < 0.0 || t.abs() > 1.0 {
                -180.0 / PI_f32 * elevation_corr
            } else {
                (t + elevation_corr).sqrt().asin() * 180.0 / PI_f32
            }
        })
        .collect();
    (elevations_low, elevation_normal, elevation_high)
}
