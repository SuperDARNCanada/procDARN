use std::collections::HashMap;
use crate::dmap::DmapData;

/// Creates the lag table based on the data.
fn create_lag_list(record_hash: HashMap<String, DmapData>) -> Vec<[i32; 2]> {
    let num_lags = record_hash.get("mplgs").unwrap();
    let pulse_table = record_hash.get("ptab").unwrap();
    let num_pulses = record_hash.get("mppul").unwrap();
    let lag_table = record_hash.get("ltab").unwrap();
    let multi_pulse_increment = record_hash.get("mpinc").unwrap();
    let sample_separation = record_hash.get("smsep").unwrap();

    let mut lags = vec![];
    for i in 0..num_lags {
        let
    }
    vec![]
}