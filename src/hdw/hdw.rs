use crate::error::BackscatterError;
use chrono::NaiveDateTime;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug)]
pub struct HdwInfo {
    pub station_id: i16,
    pub valid_from: NaiveDateTime,
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
    pub max_num_beams: i16,
}

impl HdwInfo {
    pub fn new(station_id: i16, datetime: NaiveDateTime) -> Result<HdwInfo, BackscatterError> {
        let site_name = match station_id {
            209 => "ade",
            208 => "adw",
            33 => "bks",
            24 => "bpk",
            66 => "cly",
            207 => "cve",
            206 => "cvw",
            96 => "dce",
            97 => "dcn",
            512 => "ekb",
            205 => "fhe",
            204 => "fhw",
            21 => "fir",
            1 => "gbr",
            4 => "hal",
            10 => "han",
            41 => "hkw",
            40 => "hok",
            211 => "ice",
            210 => "icw",
            64 => "inv",
            50 => "jme",
            3 => "kap",
            15 => "ker",
            7 => "kod",
            16 => "ksr",
            90 => "lyr",
            20 => "mcm",
            6 => "pgr",
            9 => "pyk",
            65 => "rkn",
            11 => "san",
            5 => "sas",
            2 => "sch",
            22 => "sps",
            8 => "sto",
            13 => "sye",
            12 => "sys",
            14 => "tig",
            0 => "tst",
            18 => "unw",
            32 => "wal",
            19 => "zho",
            _ => Err(BackscatterError::new("Invalid station id"))?,
        };
        let raw_hdw_dir = env::var_os("HDW_DIR").unwrap();
        let hdw_dir = raw_hdw_dir.to_str().unwrap();
        let hdw_file = format!("{}hdw.dat.{}", hdw_dir, site_name);
        let mut hdw_params: Vec<HdwInfo> = vec![];
        let file =
            File::open(hdw_file).map_err(|_| BackscatterError::new("Unable to open hdw file"))?;
        let reader = BufReader::new(file).lines();
        for line in reader {
            let line =
                line.map_err(|_| BackscatterError::new("Unable to read line from hdw file"))?;
            if !line.starts_with('#') {
                let elements: Vec<&str> = line.split_whitespace().collect();
                let date = elements[2];
                let time = elements[3];
                let validity_date = NaiveDateTime::parse_from_str(
                    format!("{} {}", date, time).as_str(),
                    "%Y%m%d %H:%M:%S",
                )
                .map_err(|_| BackscatterError::new("Unable to read station id from hdw file"))?;

                if datetime < validity_date {
                    break;
                } //
                hdw_params.push(HdwInfo {
                    station_id: elements[0].parse::<i16>().map_err(|_| {
                        BackscatterError::new("Unable to read station id from hdw file")
                    })?,
                    valid_from: validity_date,
                    latitude: elements[4].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read latitude from hdw file")
                    })?,
                    longitude: elements[5].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read longitude from hdw file")
                    })?,
                    altitude: elements[6].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read altitude from hdw file")
                    })?,
                    boresight: elements[7].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read boresight from hdw file")
                    })?,
                    boresight_shift: elements[8].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read boresightshift from hdw file")
                    })?,
                    beam_separation: elements[9].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read beam separation from hdw file")
                    })?,
                    velocity_sign: elements[10].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read velocity sign from hdw file")
                    })?,
                    phase_sign: elements[11].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read phase sign from hdw file")
                    })?,
                    tdiff_a: elements[12].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read tdiff A from hdw file")
                    })?,
                    tdiff_b: elements[13].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read tdiff B from hdw file")
                    })?,
                    intf_offset_x: elements[14].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read intf offset X from hdw file")
                    })?,
                    intf_offset_y: elements[15].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read intf offset Y from hdw file")
                    })?,
                    intf_offset_z: elements[16].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read intf offset Z from hdw file")
                    })?,
                    rx_rise_time: elements[17].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read rx rise time from hdw file")
                    })?,
                    rx_atten_step: elements[18].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to read rx attenuation from hdw file")
                    })?,
                    attenuation_stages: elements[19].parse::<f32>().map_err(|_| {
                        BackscatterError::new("Unable to attenuation stages from hdw file")
                    })?,
                    max_num_ranges: elements[20].parse::<i16>().map_err(|_| {
                        BackscatterError::new("Unable to read max number of ranges from hdw file")
                    })?,
                    max_num_beams: elements[21].parse::<i16>().map_err(|_| {
                        BackscatterError::new("Unable to read max number of beams from hdw file")
                    })?,
                })
            }
        }
        hdw_params
            .pop()
            .ok_or_else(|| BackscatterError::new("No valid lines found in hdw file"))
    }
}
