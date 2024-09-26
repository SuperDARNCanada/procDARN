use dmap::error::DmapError;
use dmap::formats::rawacf::RawacfRecord;
use dmap::types::DmapField;
use numpy::{Ix1, Ix2, Ix3};
use numpy::ndarray::{Array1, Array2, Array3, ArrayD};

pub(crate) struct Rawacf {
    // Scalar fields
    pub radar_revision_major: i8,
    pub radar_revision_minor: i8,
    pub origin_code: i8,
    pub origin_time: String,
    pub origin_command: String,
    pub cp: i16,
    pub stid: i16,
    pub time_yr: i16,
    pub time_mo: i16,
    pub time_dy: i16,
    pub time_hr: i16,
    pub time_mt: i16,
    pub time_sc: i16,
    pub time_us: i32,
    pub txpow: i16,
    pub nave: i16,
    pub atten: i16,
    pub lagfr: i16,
    pub smsep: i16,
    pub ercod: i16,
    pub stat_agc: i16,
    pub stat_lopwr: i16,
    pub noise_search: f32,
    pub noise_mean: f32,
    pub channel: i16,
    pub bmnum: i16,
    pub bmazm: f32,
    pub scan: i16,
    pub offset: i16,
    pub rxrise: i16,
    pub intt_sc: i16,
    pub intt_us: i32,
    pub txpl: i16,
    pub mpinc: i16,
    pub mppul: i16,
    pub mplgs: i16,
    pub nrang: i16,
    pub frang: i16,
    pub rsep: i16,
    pub xcf: i16,
    pub tfreq: i16,
    pub mxpwr: i32,
    pub lvmax: i32,
    pub combf: String,
    // pub rawacf_revision_major: i32,
    // pub rawacf_revision_minor: i32,
    // pub thr: f32,

    // Optional scalar fields
    pub mplgexs: Option<i16>,
    pub ifmode: Option<i16>,

    // Vector fields
    pub ptab: Array1<i16>,
    pub ltab: Array2<i16>,
    pub pwr0: Array1<f32>,
    pub slist: Array1<i16>,
    pub acfd: Array3<f32>,

    // Optional vector fields
    pub xcfd: Option<Array3<f32>>,
}
impl TryFrom<&RawacfRecord> for Rawacf {
    type Error = DmapError;
    fn try_from(value: &RawacfRecord) -> Result<Self, Self::Error> {
        let scalar_getter = |key: &str| -> Result<&DmapField, DmapError> {
            value
                .get(&key.to_string())
                .ok_or_else(|| DmapError::InvalidScalar(key.to_string()))
        };
        let opt_scalar_getter = |key: &str| -> Option<&DmapField> { value.get(&key.to_string()) };
        let vector_getter = |key: &str| -> Result<&DmapField, DmapError> {
            value
                .get(&key.to_string())
                .ok_or_else(|| DmapError::InvalidVector(key.to_string()))
        };
        let opt_vector_getter = |key: &str| -> Option<&DmapField> { value.get(&key.to_string()) };
        Ok(Rawacf {
            radar_revision_major: scalar_getter("radar.revision.major")?.clone().try_into()?,
            radar_revision_minor: scalar_getter("radar.revision.minor")?.clone().try_into()?,
            origin_code: scalar_getter("origin.code")?.clone().try_into()?,
            origin_time: scalar_getter("origin.time")?.clone().try_into()?,
            origin_command: scalar_getter("origin.command")?.clone().try_into()?,
            cp: scalar_getter("cp")?.clone().try_into()?,
            stid: scalar_getter("stid")?.clone().try_into()?,
            time_yr: scalar_getter("time.yr")?.clone().try_into()?,
            time_mo: scalar_getter("time.mo")?.clone().try_into()?,
            time_dy: scalar_getter("time.dy")?.clone().try_into()?,
            time_hr: scalar_getter("time.hr")?.clone().try_into()?,
            time_mt: scalar_getter("time.mt")?.clone().try_into()?,
            time_sc: scalar_getter("time.sc")?.clone().try_into()?,
            time_us: scalar_getter("time.us")?.clone().try_into()?,
            txpow: scalar_getter("txpow")?.clone().try_into()?,
            nave: scalar_getter("nave")?.clone().try_into()?,
            atten: scalar_getter("atten")?.clone().try_into()?,
            lagfr: scalar_getter("lagfr")?.clone().try_into()?,
            smsep: scalar_getter("smsep")?.clone().try_into()?,
            ercod: scalar_getter("ercod")?.clone().try_into()?,
            stat_agc: scalar_getter("stat.agc")?.clone().try_into()?,
            stat_lopwr: scalar_getter("stat.lopwr")?.clone().try_into()?,
            noise_search: scalar_getter("noise.search")?.clone().try_into()?,
            noise_mean: scalar_getter("noise.mean")?.clone().try_into()?,
            channel: scalar_getter("channel")?.clone().try_into()?,
            bmnum: scalar_getter("bmnum")?.clone().try_into()?,
            bmazm: scalar_getter("bmazm")?.clone().try_into()?,
            scan: scalar_getter("scan")?.clone().try_into()?,
            offset: scalar_getter("offset")?.clone().try_into()?,
            rxrise: scalar_getter("rxrise")?.clone().try_into()?,
            intt_sc: scalar_getter("intt.sc")?.clone().try_into()?,
            intt_us: scalar_getter("intt.us")?.clone().try_into()?,
            txpl: scalar_getter("txpl")?.clone().try_into()?,
            mpinc: scalar_getter("mpinc")?.clone().try_into()?,
            mppul: scalar_getter("mppul")?.clone().try_into()?,
            mplgs: scalar_getter("mplgs")?.clone().try_into()?,
            nrang: scalar_getter("nrang")?.clone().try_into()?,
            frang: scalar_getter("frang")?.clone().try_into()?,
            rsep: scalar_getter("rsep")?.clone().try_into()?,
            xcf: scalar_getter("xcf")?.clone().try_into()?,
            tfreq: scalar_getter("tfreq")?.clone().try_into()?,
            mxpwr: scalar_getter("mxpwr")?.clone().try_into()?,
            lvmax: scalar_getter("lvmax")?.clone().try_into()?,
            combf: scalar_getter("combf")?.clone().try_into()?,
            // rawacf_revision_major: scalar_getter("rawacf.revision.major")?.clone().try_into()?,
            // rawacf_revision_minor: scalar_getter("rawacf.revision.minor")?.clone().try_into()?,
            // thr: scalar_getter("thr")?.clone().try_into()?,
            mplgexs: match opt_scalar_getter("mplgexs") {
                Some(x) => Some(x.clone().try_into()?),
                None => None,
            },
            ifmode: match opt_scalar_getter("ifmode") {
                Some(x) => Some(x.clone().try_into()?),
                None => None,
            },
            ptab: <DmapField as TryInto<ArrayD<i16>>>::try_into(vector_getter("ptab")?.clone())?.into_dimensionality::<Ix1>().map_err(|e| DmapError::InvalidVector(format!("Unable to map ptab to 1D vector: {e}")))?,
            ltab: <DmapField as TryInto<ArrayD<i16>>>::try_into(vector_getter("ltab")?.clone())?.into_dimensionality::<Ix2>().map_err(|e| DmapError::InvalidVector(format!("Unable to map ltab to 2D vector: {e}")))?,
            pwr0: <DmapField as TryInto<ArrayD<f32>>>::try_into(vector_getter("pwr0")?.clone())?.into_dimensionality::<Ix1>().map_err(|e| DmapError::InvalidVector(format!("Unable to map pwr0 to 1D vector: {e}")))?,
            slist: <DmapField as TryInto<ArrayD<i16>>>::try_into(vector_getter("slist")?.clone())?.into_dimensionality::<Ix1>().map_err(|e| DmapError::InvalidVector(format!("Unable to map slist to 1D vector: {e}")))?,
            acfd: <DmapField as TryInto<ArrayD<f32>>>::try_into(vector_getter("acfd")?.clone())?.into_dimensionality::<Ix3>().map_err(|e| DmapError::InvalidVector(format!("Unable to map acfd to 3D vector: {e}")))?,
            xcfd: match opt_vector_getter("xcfd") {
                Some(x) => Some(<DmapField as TryInto<ArrayD<f32>>>::try_into(x.clone())?.into_dimensionality::<Ix3>().map_err(|e| DmapError::InvalidVector(format!("Unable to map xcfd to 3D vector: {e}")))?),
                None => None,
            },
        })
    }
}
