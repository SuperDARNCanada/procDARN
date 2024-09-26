use crate::fitting::fitacf3::fitacf_v3::Fitacf3Error;
use crate::fitting::fitacf3::fitstruct::RangeNode;
use crate::utils::hdw::HdwInfo;
use crate::utils::rawacf::Rawacf;
use dmap::formats::{dmap::Record, fitacf::FitacfRecord};
use dmap::types::DmapField;
use indexmap::IndexMap;
use numpy::ndarray::{Array, Array1, ArrayD};
use numpy::IxDyn;
use std::f32::consts::PI as PI_f32;
use std::iter::zip;
use chrono::Utc;

pub const FITACF_REVISION_MAJOR: i32 = 3;
pub const FITACF_REVISION_MINOR: i32 = 0;
pub const ORIGIN_CODE: i8 = 1;
pub const V_MAX: f32 = 30.0;
pub const W_MAX: f32 = 90.0;

pub(crate) fn determinations(
    rec: &Rawacf,
    ranges: Vec<RangeNode>,
    noise_power: f32,
    hdw: &HdwInfo,
) -> Result<FitacfRecord, Fitacf3Error> {
    let range_list: Vec<i16> = ranges.iter().map(|r| r.range_num as i16).collect();
    let lag_0_power_db: Array1<f32> = rec
        .pwr0
        .iter()
        .map(|p| {
            if p - noise_power > 0.0 {
                10.0 * ((p - noise_power) / noise_power).log10()
            } else {
                -50.0
            }
        })
        .collect();

    let mut fit_rec: IndexMap<String, DmapField> = IndexMap::new();

    fit_rec.insert(
        "radar.revision.major".to_string(),
        rec.radar_revision_major.into(),
    );
    fit_rec.insert(
        "radar.revision.minor".to_string(),
        rec.radar_revision_minor.into(),
    );
    fit_rec.insert("origin.code".to_string(), ORIGIN_CODE.into());
    let now: chrono::DateTime<Utc> = std::time::SystemTime::now().into();
    fit_rec.insert("origin.time".to_string(), format!("{}", now.format("%a %b %e %T %Y")).into());
    fit_rec.insert(
        "origin.command".to_string(),
        rec.origin_command.clone().into(),
    );
    fit_rec.insert("cp".to_string(), rec.cp.into());
    fit_rec.insert("stid".to_string(), rec.stid.into());
    fit_rec.insert("time.yr".to_string(), rec.time_yr.into());
    fit_rec.insert("time.mo".to_string(), rec.time_mo.into());
    fit_rec.insert("time.dy".to_string(), rec.time_dy.into());
    fit_rec.insert("time.hr".to_string(), rec.time_hr.into());
    fit_rec.insert("time.mt".to_string(), rec.time_mt.into());
    fit_rec.insert("time.sc".to_string(), rec.time_sc.into());
    fit_rec.insert("time.us".to_string(), rec.time_us.into());
    fit_rec.insert("txpow".to_string(), rec.txpow.into());
    fit_rec.insert("nave".to_string(), rec.nave.into());
    fit_rec.insert("atten".to_string(), rec.atten.into());
    fit_rec.insert("lagfr".to_string(), rec.lagfr.into());
    fit_rec.insert("smsep".to_string(), rec.smsep.into());
    fit_rec.insert("ercod".to_string(), rec.ercod.into());
    fit_rec.insert("stat.agc".to_string(), rec.stat_agc.into());
    fit_rec.insert("stat.lopwr".to_string(), rec.stat_lopwr.into());
    fit_rec.insert("noise.search".to_string(), rec.noise_search.into());
    fit_rec.insert("noise.mean".to_string(), rec.noise_mean.into());
    fit_rec.insert("channel".to_string(), rec.channel.into());
    fit_rec.insert("bmnum".to_string(), rec.bmnum.into());
    fit_rec.insert("bmazm".to_string(), rec.bmazm.into());
    fit_rec.insert("scan".to_string(), rec.scan.into());
    fit_rec.insert("offset".to_string(), rec.offset.into());
    fit_rec.insert("rxrise".to_string(), rec.rxrise.into());
    fit_rec.insert("intt.sc".to_string(), rec.intt_sc.into());
    fit_rec.insert("intt.us".to_string(), rec.intt_us.into());
    fit_rec.insert("txpl".to_string(), rec.txpl.into());
    fit_rec.insert("mpinc".to_string(), rec.mpinc.into());
    fit_rec.insert("mppul".to_string(), rec.mppul.into());
    fit_rec.insert("mplgs".to_string(), rec.mplgs.into());
    fit_rec.insert("nrang".to_string(), rec.nrang.into());
    fit_rec.insert("frang".to_string(), rec.frang.into());
    fit_rec.insert("rsep".to_string(), rec.rsep.into());
    fit_rec.insert("xcf".to_string(), rec.xcf.into());
    fit_rec.insert("tfreq".to_string(), rec.tfreq.into());
    fit_rec.insert("mxpwr".to_string(), rec.mxpwr.into());
    fit_rec.insert("lvmax".to_string(), rec.lvmax.into());
    fit_rec.insert("combf".to_string(), rec.combf.clone().into());
    fit_rec.insert("ptab".to_string(), rec.ptab.clone().into_dyn().into());
    fit_rec.insert("ltab".to_string(), rec.ltab.clone().into_dyn().into());
    fit_rec.insert(
        "fitacf.revision.major".to_string(),
        FITACF_REVISION_MAJOR.into(),
    );
    fit_rec.insert(
        "fitacf.revision.minor".to_string(),
        FITACF_REVISION_MINOR.into(),
    );
    fit_rec.insert("noise.lag0".to_string(), 0.0_f32.into());
    fit_rec.insert("noise.vel".to_string(), 0.0_f32.into());
    if let Some(x) = rec.ifmode {
        fit_rec.insert("ifmode".to_string(), x.into());
    }
    if let Some(x) = rec.mplgexs {
        fit_rec.insert("mplgexs".to_string(), x.into());
    }
    fit_rec.insert("noise.sky".to_string(), noise_power.into());
    fit_rec.insert("pwr0".to_string(), lag_0_power_db.into_dyn().into());

    if !range_list.is_empty() {
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
            299792458.0 * hdw.velocity_sign / (4.0 * PI_f32 * rec.tfreq as f32 * 1000.0);
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
        let width_conversion: f32 = 299792458.0 * 2.0 / (4.0 * PI_f32 * rec.tfreq as f32 * 1000.0);
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
            299792458.0 * 2.0_f32.ln().sqrt() / (PI_f32 * rec.tfreq as f32 * 1000.0);
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
        let xcfs = &rec.xcfd.as_ref().expect("Unable to make fitacf xcf_phi0");
        let xcf_phi0: Vec<f32> = ranges
            .iter()
            .map(|r| {
                xcfs[[r.range_idx, 0, 1]]
                    .atan2(xcfs[[r.range_idx, 0, 0]])
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

        let float_zeros: ArrayD<f32> = Array::zeros(IxDyn(&[quality_flag.len()]));
        let i8_zeros: ArrayD<i8> = Array::zeros(IxDyn(&[quality_flag.len()]));

        fit_rec.insert(
            "slist".to_string(),
            Array::from_vec(range_list).into_dyn().into(),
        );
        fit_rec.insert(
            "nlag".to_string(),
            Array::from_vec(num_lags).into_dyn().into(),
        );
        fit_rec.insert(
            "qflg".to_string(),
            Array::from_vec(quality_flag).into_dyn().into(),
        );
        fit_rec.insert(
            "gflg".to_string(),
            Array::from_vec(groundscatter_flag).into_dyn().into(),
        );
        fit_rec.insert(
            "p_l".to_string(),
            Array::from_vec(power_linear).into_dyn().into(),
        );
        fit_rec.insert(
            "p_l_e".to_string(),
            Array::from_vec(power_linear_error).into_dyn().into(),
        );
        fit_rec.insert(
            "p_s".to_string(),
            Array::from_vec(power_quadratic).into_dyn().into(),
        );
        fit_rec.insert(
            "p_s_e".to_string(),
            Array::from_vec(power_quadratic_error).into_dyn().into(),
        );
        fit_rec.insert("v".to_string(), Array::from_vec(velocity).into_dyn().into());
        fit_rec.insert(
            "v_e".to_string(),
            Array::from_vec(velocity_error).into_dyn().into(),
        );
        fit_rec.insert(
            "w_l".to_string(),
            Array::from_vec(spectral_width_linear).into_dyn().into(),
        );
        fit_rec.insert(
            "w_l_e".to_string(),
            Array::from_vec(spectral_width_linear_error)
                .into_dyn()
                .into(),
        );
        fit_rec.insert(
            "w_s".to_string(),
            Array::from_vec(spectral_width_quadratic).into_dyn().into(),
        );
        fit_rec.insert(
            "w_s_e".to_string(),
            Array::from_vec(spectral_width_quadratic_error)
                .into_dyn()
                .into(),
        );
        fit_rec.insert(
            "sd_l".to_string(),
            Array::from_vec(std_dev_linear).into_dyn().into(),
        );
        fit_rec.insert(
            "sd_s".to_string(),
            Array::from_vec(std_dev_quadratic).into_dyn().into(),
        );
        fit_rec.insert(
            "sd_phi".to_string(),
            Array::from_vec(std_dev_phi).into_dyn().into(),
        );
        fit_rec.insert("x_qflg".to_string(), i8_zeros.clone().into());
        fit_rec.insert("x_gflg".to_string(), i8_zeros.into());
        fit_rec.insert("x_p_l".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_p_l_e".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_p_s".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_p_s_e".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_v".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_v_e".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_w_l".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_w_l_e".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_w_s".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_w_s_e".to_string(), float_zeros.clone().into());
        fit_rec.insert(
            "phi0".to_string(),
            Array::from_vec(xcf_phi0).into_dyn().into(),
        );
        fit_rec.insert(
            "phi0_e".to_string(),
            Array::from_vec(xcf_phi0_err).into_dyn().into(),
        );
        fit_rec.insert(
            "elv".to_string(),
            Array::from_vec(elevation_normal).into_dyn().into(),
        );
        fit_rec.insert(
            "elv_low".to_string(),
            Array::from_vec(elevation_low).into_dyn().into(),
        );
        fit_rec.insert(
            "elv_high".to_string(),
            Array::from_vec(elevation_high).into_dyn().into(),
        );
        fit_rec.insert("x_sd_l".to_string(), float_zeros.clone().into());
        fit_rec.insert("x_sd_s".to_string(), float_zeros.into());
        fit_rec.insert(
            "x_sd_phi".to_string(),
            Array::from_vec(xcf_phi_std_dev).into_dyn().into(),
        );
    }
    let new_rec = FitacfRecord::new(&mut fit_rec).map_err(|e| {
        Fitacf3Error::Message(format!(
            "Could not create valid Fitacf record from results: {e}"
        ))
    })?;
    Ok(new_rec)
}

fn calculate_elevation(
    ranges: &[RangeNode],
    rec: &Rawacf,
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
    let phi_0 = (hdw.beam_separation * (rec.bmnum as f32 - azimuth_offset) * PI_f32 / 180.0).cos();
    let wave_num = 2.0 * PI_f32 * rec.tfreq as f32 * 1000.0 / 299792458.0;
    let cable_offset = -2.0 * PI_f32 * rec.tfreq as f32 * 1000.0 * hdw.tdiff_a * 1.0e-6;
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
            let mut y = x
                + 2.0 * PI_f32 * ((phase_diff_max - x) / (2.0 * PI_f32)).floor()
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
