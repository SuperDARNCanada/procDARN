use std::collections::HashMap;
use std::env;
use std::fs::read_dir;

pub fn parse_hdw_files() -> HashMap<i16, HdwInfo> {
    let mut hdw_params = HashMap::new();
    let hdw_dir = env::var_os("OUT_DIR").unwrap();
    for entry in read_dir(hdw_dir)? {
        let entry = entry?;
        let path = entry.path();
        let station_id,

    }
}



pub struct HdwInfo {
    pub station_id: i16,
    pub year: i16,
    pub second: i16,
    pub latitude: f32,
    pub longitude: f32,
    pub altitude: f32,
    pub boresight: f32,
    pub boresight_shift: f32,
    pub beam_separation: f32,
    pub velocity_sign: f32,
    pub phase_sign: f32,
    pub tdiff_a: f32,
    pub tdiff_b: f32,
    pub intf_offset_x: f32,
    pub intf_offset_y: f32,
    pub intf_offset_z: f32,
    pub rx_rise_time: f32,
    pub rx_atten_step: f32,
    pub attenuation_stages: f32,
    pub max_num_ranges: i16,
    pub max_num_beams: i16
}
impl HdwInfo {
    pub fn new(station_id: i16) {
        // let site_name = match station_id {
        //     209 => "ade",
        //     208 => "adw",
        //     33 => "bks",
        //     24 => "bpk",
        //     66 => "cly",
        //     207 => "cve",
        //     206 => "cvw",
        //     96 => "dce",
        //     97 => "dcn",
        //     512 => "ekb",
        //     205 => "fhe",
        //     204 => "fhw",
        //     21 => "fir",
        //     1 => "gbr",
        //     4 => "hal",
        //     10 => "han",
        //     41 => "hkw",
        //     40 => "hok",
        //     211 => "ice",
        //
        // }
    }

    pub fn from_file(path: PathBuf) -> (i16, HdwInfo) {

    }
}