
pub const VELOCITY_ERROR_MIN: f64 = 100.0;
pub const POWER_LIN_ERROR_MIN: f64 = 1.0;
pub const WIDTH_LIN_ERROR_MIN: f64 = 1.0;

#[derive(Debug)]
pub struct GridBeam {
    pub beam: i32,          // bm in RST
    pub first_range: i32,   // frang in RST, km
    pub range_sep: i32,     // rsep in RST, km
    pub rx_rise: i32,       // rxrise in RST, microseconds?
    pub num_ranges: i32,    // nrang in RST
    pub azimuth: Vec<f64>,  // azm in RST, degrees?
    pub ival: Vec<f64>,     // ival in RST
    pub index: Vec<i32>,    // inx in RST
}

#[derive(Debug)]
pub struct GridPoint {
    pub max: i32,                   // max in RST
    pub count: i32,                 // cnt in RST
    pub reference: i32,             // ref in RST
    pub magnetic_lat: f64,          // mlat in RST
    pub magnetic_lon: f64,          // mlon in RST
    pub azimuth: f64,               // azm in RST, degrees?
    pub velocity_median: f64,       // vel.median in RST, m/s
    pub velocity_median_n: f64,     // vel.median_n in RST, m/s
    pub velocity_median_e: f64,     // vel.median_e in RST, m/s
    pub velocity_stddev: f64,       // vel.sd in RST, m/s
    pub power_median: f64,          // pwr.median in RST, a.u. in linear scale
    pub power_stddev: f64,          // pwr.sd in RST, a.u. in linear scale
    pub spectral_width_median: f64, // wdt.median in RST, m/s
    pub spectral_width_stddev: f64, // wdt.sd in RST, m/s
}
impl GridPoint {
    pub fn clear(&mut self) {
        self.azimuth = 0.0;
        self.velocity_median_n = 0.0;
        self.velocity_median_e = 0.0;
        self.velocity_stddev = 0.0;
        self.power_median = 0.0;
        self.power_stddev = 0.0;
        self.spectral_width_median = 0.0;
        self.spectral_width_stddev = 0.0;
        self.count = 0;
    }
}

#[derive(Debug)]
pub struct GridTable {
    pub start_time: f64,            // st_time in RST
    pub end_time: f64,              // ed_time in RST
    pub channel: i32,               // chn in RST
    pub status: i32,                // status in RST
    pub station_id: i32,            // st_id in RST
    pub program_id: i32,            // prog_id in RST
    pub num_scans: i32,             // nscan in RST
    pub num_points_npnt: i32,       // npnt in RST
    pub freq: f64,                  // freq in RST
    pub noise_mean: f64,            // noise.mean in RST
    pub noise_stddev: f64,          // noise.sd in RST
    pub groundscatter: i32,         // gsct in RST
    pub min_power: f64,             // min[0] in RST, a.u. in linear scale
    pub min_velocity: f64,          // min[1] in RST, m/s
    pub min_spectral_width: f64,    // min[2] in RST, m/s
    pub min_velocity_error: f64,    // min[3] in RST, m/s
    pub max_power: f64,             // max[0] in RST, a.u. in linear scale
    pub max_velocity: f64,          // max[1] in RST, m/s
    pub max_spectral_width: f64,    // max[2] in RST, m/s
    pub max_velocity_error: f64,    // max[3] in RST, m/s
    pub num_beams: i32,             // bnum in RST
    pub beams: Vec<GridBeam>,       // bm in RST
    pub num_points_pnum: i32,       // pnum in RST
    pub points: Vec<GridPoint>,     // pnt in RST
}
impl GridTable {
    pub fn new(
        index: usize,
        range_num: usize,
        record: &RawacfRecord,
        lags: &[LagNode],
    ) -> Result<RangeNode, Fitacf3Error> {
        let cross_range_interference =
            RangeNode::calculate_cross_range_interference(range_num, record);
        let alpha_2 =
            RangeNode::calculate_alphas(range_num, &cross_range_interference, record, lags);
        let phases = PhaseNode::new(record, "acfd", lags, index)?;
        let elevations = PhaseNode::new(record, "xcfd", lags, index)?;
        let powers = PowerNode::new(record, lags, index, range_num, &alpha_2);
        Ok(RangeNode {
            range_idx: index,
            range_num,
            cross_range_interference,
            refractive_idx: 1.0,
            power_alpha_2: alpha_2.clone(),
            phase_alpha_2: alpha_2,
            phases,
            powers,
            elev: elevations,
            lin_pwr_fit: None,
            quad_pwr_fit: None,
            lin_pwr_fit_err: None,
            quad_pwr_fit_err: None,
            phase_fit: None,
            elev_fit: None,
        })
    }
    pub fn clear(&mut self) {
        for mut p in self.points {
            p.clear()
        }
    }
    pub fn test(&self, &scan: RadarScan) -> bool {
        
        true
    }
}