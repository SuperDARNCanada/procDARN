
pub struct HdwInfo {
    pub station_id: i16,
    pub year: i16,
    pub second: i16,
    pub latitude: f32,
    pub longitude: f32,
    pub altitude: f32,
    pub boresight: f32,
    pub beam_separation: f32,
    pub velocity_sign: f32,
    pub rx_atten_step: f32,
    pub tdiff: f32,
    pub phase_sign: f32,
    pub intf_offset_x: f32,
    pub intf_offset_y: f32,
    pub intf_offset_z: f32,
    pub rx_rise_time: f32,
    pub attenuation_stages: f32,
    pub max_num_ranges: i16,
    pub max_num_beams: i16
}