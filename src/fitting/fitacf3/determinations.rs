use crate::fitting::fitacf3::fitacf_v3::Fitacf3Error;
use crate::fitting::fitacf3::fitstruct::RangeNode;
use crate::utils::hdw::HdwInfo;
use crate::utils::rawacf::Rawacf;
use chrono::Utc;
use dmap::formats::fitacf::FitacfRecord;
use dmap::record::Record;
use dmap::types::DmapField;
use indexmap::IndexMap;
use numpy::ndarray::{Array, Array1};
use std::f32::consts::PI as PI_f32;
use std::f64::consts::PI as PI_f64;
use std::iter::zip;

pub const FITACF_REVISION_MAJOR: i32 = 3;
pub const FITACF_REVISION_MINOR: i32 = 0;
const LIGHTSPEED: f32 = 299_792_458.0;
const KHZ_TO_HZ: f32 = 1000.0;
const US_TO_S: f32 = 1e-6;
pub const ORIGIN_CODE: i8 = 1;
pub const V_MAX: f32 = 30.0;
pub const W_MAX: f32 = 90.0;

pub(crate) fn determinations(
    rec: &Rawacf,
    ranges: &[RangeNode],
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
    fit_rec.insert(
        "origin.time".to_string(),
        format!("{}", now.format("%a %b %e %T %Y")).into(),
    );
    fit_rec.insert(
        "origin.command".to_string(),
        rec.origin_command.clone().into(),
    ); // todo: Get the invocation
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
    fit_rec.insert("tfreq".to_string(), (rec.tfreq.round() as i16).into());
    fit_rec.insert("mxpwr".to_string(), rec.mxpwr.into());
    fit_rec.insert("lvmax".to_string(), rec.lvmax.into());
    fit_rec.insert("combf".to_string(), rec.combf.clone().into());
    fit_rec.insert("ptab".to_string(), rec.ptab.clone().into_dyn().into());
    fit_rec.insert("ltab".to_string(), rec.ltab.clone().into_dyn().into());
    fit_rec.insert("algorithm".to_string(), "fitacf3".to_string().into());
    fit_rec.insert("tdiff".to_string(), hdw.tdiff_a.into());
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
    } else {
        fit_rec.insert("ifmode".to_string(), <DmapField as From<i16>>::from(0));
    }
    if let Some(x) = rec.mplgexs {
        fit_rec.insert("mplgexs".to_string(), x.into());
    } else {
        fit_rec.insert("mplgexs".to_string(), 0_i16.into());
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
                    / 10.0_f32.ln()
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
                    / 10.0_f32.ln()
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
                    / 10.0_f32.ln()
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
                    / 10.0_f32.ln()
            })
            .collect();
        let velocity_conversion: f32 =
            LIGHTSPEED * hdw.velocity_sign / (4.0 * PI_f32 * rec.tfreq * KHZ_TO_HZ);
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
        let width_conversion: f32 = LIGHTSPEED * 2.0 / (4.0 * PI_f32 * rec.tfreq * KHZ_TO_HZ);
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
            LIGHTSPEED * 2.0_f32.ln().sqrt() / (PI_f32 * rec.tfreq * KHZ_TO_HZ);
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
            .map(|(v, w)| i8::from(v.abs() - (V_MAX - w * (V_MAX / W_MAX)) < 0.0))
            .collect();
        let xcfs = &rec.xcfd.as_ref().expect("Unable to make fitacf xcf_phi0");
        let xcf_phi0: Vec<f32> = ranges
            .iter()
            .map(|r| xcfs[[r.range_idx, 0, 1]].atan2(xcfs[[r.range_idx, 0, 0]]) * hdw.phase_sign)
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
        let (elevation_phi0, elevation_intercept) =
            calculate_elevation_v2(ranges, rec, &xcf_phi0, hdw);
        let (_, elevation_intercept_error, _) = calculate_elevation(ranges, rec, &xcf_phi0, hdw);

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
            Array::from_vec(elevation_phi0).into_dyn().into(),
        );
        fit_rec.insert(
            "elv_error".to_string(),
            Array::from_vec(elevation_intercept_error).into_dyn().into(),
        );
        fit_rec.insert(
            "elv_fitted".to_string(),
            Array::from_vec(elevation_intercept).into_dyn().into(),
        );
        fit_rec.insert(
            "x_sd_phi".to_string(),
            Array::from_vec(xcf_phi_std_dev).into_dyn().into(),
        );
    }
    let new_rec = FitacfRecord::new(&mut fit_rec).map_err(|e| {
        Fitacf3Error::BadFit(format!(
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
    let x = hdw.intf_offset_x as f64;
    let y = hdw.intf_offset_y as f64;
    let z = hdw.intf_offset_z as f64;

    let array_separation: f64 = (x * x + y * y + z * z).sqrt();
    let mut elevation_corr = (z / array_separation).asin();
    let phi_sign: f64;
    if y > 0.0 {
        phi_sign = 1.0;
    } else {
        phi_sign = -1.0;
        elevation_corr *= -1.0;
    }
    let azimuth_offset = f32::from(hdw.max_num_beams) / 2.0 - 0.5;
    let cos_phi_0 = (hdw.boresight_shift
        + hdw.beam_separation * (f32::from(rec.bmnum) - azimuth_offset))
        .to_radians()
        .cos() as f64;
    let wave_num = 2.0 * PI_f64 * (rec.tfreq * KHZ_TO_HZ) as f64 / LIGHTSPEED as f64;
    let cable_offset =
        -2.0 * PI_f64 * (rec.tfreq * KHZ_TO_HZ) as f64 * (hdw.tdiff_a * US_TO_S) as f64;
    let phase_diff_max = phi_sign * wave_num * array_separation * cos_phi_0 + cable_offset;

    let mut psi: Vec<f64> = ranges
        .iter()
        .map(|r| {
            let a = r
                .elev_fit
                .as_ref()
                .expect("Unable to find elevation without fitted elevation")
                .intercept;
            let mut psi_raw = a + 2.0 * PI_f64 * ((phase_diff_max - a) / (2.0 * PI_f64)).floor();
            if phi_sign < 0.0 {
                psi_raw += 2.0 * PI_f64;
            }
            psi_raw - cable_offset
        })
        .collect();
    let mut psi_kd: Vec<f64> = psi
        .iter()
        .map(|p| p / (wave_num * array_separation))
        .collect();
    let mut theta: Vec<f64> = psi_kd
        .iter()
        .map(|p| cos_phi_0 * cos_phi_0 - p * p)
        .collect();
    let elevation_intercept: Vec<f32> = theta
        .iter()
        .map(|&t| {
            if t < 0.0 || t.abs() > 1.0 {
                -elevation_corr.to_degrees() as f32
            } else {
                t.sqrt().asin().to_degrees() as f32
            }
        })
        .collect();
    let psi_k2d2: Vec<f64> = psi
        .iter()
        .map(|p| p / (wave_num * wave_num * array_separation * array_separation))
        .collect();
    let df_by_dy: Vec<f64> = zip(psi_k2d2.iter(), theta.iter())
        .map(|(p, t)| p / (t - t * t).sqrt())
        .collect();
    let errors: Vec<f64> = ranges
        .iter()
        .map(|r| {
            r.elev_fit
                .as_ref()
                .expect("Unable to calculate elevation errors")
                .variance_intercept
        })
        .collect();
    let elevation_intercept_error: Vec<f32> = zip(errors.iter(), df_by_dy.iter())
        .map(|(e, d)| (e.sqrt() * d.abs()).to_degrees() as f32)
        .collect();

    // This time, use the xcf lag0 phase
    psi = xcf_phi0
        .iter()
        .map(|&p| {
            let mut ps = p as f64
                + 2.0 * PI_f64 * ((phase_diff_max - p as f64) / (2.0 * PI_f64)).floor()
                - cable_offset;
            if phi_sign < 0.0 {
                ps += 2.0 * PI_f64;
            }
            ps
        })
        .collect();
    psi_kd = psi
        .iter()
        .map(|p| p / (wave_num * array_separation))
        .collect();
    theta = psi_kd
        .iter()
        .map(|p| cos_phi_0 * cos_phi_0 - p * p)
        .collect();
    let elevation_phi0: Vec<f32> = theta
        .iter()
        .map(|&t| {
            if t < 0.0 || t.abs() > 1.0 {
                -elevation_corr.to_degrees() as f32
            } else {
                (t + elevation_corr).sqrt().asin().to_degrees() as f32
            }
        })
        .collect();
    (
        elevation_phi0,
        elevation_intercept_error,
        elevation_intercept,
    )
}

fn calculate_elevation_v2(
    ranges: &[RangeNode],
    rec: &Rawacf,
    xcf_phi0: &[f32],
    hdw: &HdwInfo,
) -> (Vec<f32>, Vec<f32>) {
    let x = hdw.intf_offset_x;
    let y = hdw.intf_offset_y;
    let z = hdw.intf_offset_z;

    let psi_sign: f32 = if y > 0.0 { 1.0 } else { -1.0 };

    let azimuth_offset = f32::from(hdw.max_num_beams) / 2.0 - 0.5;
    let phi_0 = (hdw.boresight_shift
        + hdw.beam_separation * (f32::from(rec.bmnum) - azimuth_offset))
        .to_radians();
    let cos_phi_0 = phi_0.cos(); // cp0
    let sin_phi_0 = phi_0.sin(); // sp0

    let wave_num = 2.0 * PI_f32 * rec.tfreq * KHZ_TO_HZ / LIGHTSPEED;
    let cable_offset_rad = -2.0 * PI_f32 * rec.tfreq * KHZ_TO_HZ * hdw.tdiff_a * US_TO_S; // psi_ele

    let mut elv_of_max_psi = (psi_sign * z * cos_phi_0 / (y * y + z * z).sqrt()).asin(); // a0
    if elv_of_max_psi < 0. {
        elv_of_max_psi = 0.0;
    }

    let cos_elv_of_max_psi = elv_of_max_psi.cos(); // ca0
    let sin_elv_of_max_psi = elv_of_max_psi.sin(); // sa0

    let psi_max = cable_offset_rad
        + wave_num
            * (x * sin_phi_0
                + y * (cos_elv_of_max_psi * cos_elv_of_max_psi - sin_phi_0 * sin_phi_0).sqrt()
                + z * sin_elv_of_max_psi);

    let num_phase_jump_func = {
        if y > 0.0 {
            f64::floor
        } else {
            f64::ceil
        }
    };
    let psi_calc = |p: &f32| -> f64 {
        let delta_psi = (psi_max - p) as f64;
        (p + 0.) as f64 + 2.0 * PI_f64 * num_phase_jump_func(delta_psi / (2.0 * PI_f64))
    };
    let e_calc = |p: &f64| -> f64 {
        (p / (2.0 * PI_f64 * (rec.tfreq * KHZ_TO_HZ) as f64) + (hdw.tdiff_a * US_TO_S) as f64)
            * LIGHTSPEED as f64
            - (x * sin_phi_0) as f64
    };
    let elv_calc = |e: &f64| -> f32 {
        (e * z as f64
            + (e * e * (z * z) as f64
                - ((y * y) as f64 + (z * z) as f64)
                    * (e * e - (y * y) as f64 * (cos_phi_0 * cos_phi_0) as f64))
                .sqrt()
                / ((y * y) as f64 + (z * z) as f64))
            .asin()
            .to_degrees() as f32
    };

    let psi_normal: Vec<f64> = xcf_phi0.iter().map(psi_calc).collect();
    let e_normal: Vec<f64> = psi_normal.iter().map(e_calc).collect();
    let elv_normal = e_normal // called alpha in RST
        .iter()
        .map(elv_calc)
        .collect();

    let psi_fitted: Vec<f64> = ranges
        .iter()
        .map(|r| {
            let p = r
                .elev_fit
                .as_ref()
                .expect("Unable to find elevation without fitted elevation")
                .intercept as f32
                * hdw.phase_sign;
            psi_calc(&p)
        })
        .collect();
    let e_fitted: Vec<f64> = psi_fitted.iter().map(e_calc).collect();
    let elv_fitted = e_fitted.iter().map(elv_calc).collect();

    (elv_normal, elv_fitted)
}
