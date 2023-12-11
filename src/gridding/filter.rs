use crate::error::BackscatterError;
use crate::utils::scan::RadarScan;

pub const MAX_BEAM: i32 = 256;
pub const FILTER_HEIGHT: i32 = 3;
pub const FILTER_WIDTH: i32 = 3;
pub const FILTER_DEPTH: i32 = 3;

/// Performs median filtering on a sequence of RadarScans. The filter operates on each range/beam
/// cell, with a 3x3x3 weighted kernel of range/beam/time. If the weighted sum of valid cells in
/// the kernel exceeds a threshold, the median value of each parameter (velocity, power, and
/// spectral width) is determined from the kernel. Otherwise, the output cell is considered empty.
/// The associated parameter errors are calculated from the standard deviations of the input
/// parameters.
/// Called FilterRadarScan in filter.c of RST.
pub fn median_filter(
    mode: i32,
    depth: i32,
    index: i32,
    param: i32,
    isort: i32,
    scans: &[RadarScan],
) -> Result<RadarScan, BackscatterError> {
    let mut out_scan = RadarScan{ ..Default::default() };
    let cnum: i32 = 0;
    let count: i32 = 0;
    let mut max_beam: i32 = -1;
    let mut max_range: i32 = 1000;

    let filter_depth: usize;
    if depth > FILTER_DEPTH {
        filter_depth = FILTER_DEPTH as usize;
    } else {
        filter_depth = depth as usize;
    }

    // Find the largest beam number and range number in all the scans
    for i in 0..filter_depth {
        for beam in scans[i].beams.iter() {
            if beam.beam > max_beam {
                max_beam = beam.beam + 1;   // Add one since beam number is indexed from 0
            }
            if beam.num_ranges > max_range {
                max_range = beam.num_ranges;
            }
        }
    }

    // Calculate weight of each cell in the kernel.
    //   <---> beam
    //   1 1 1    2 2 2    1 1 1  ^
    //   1 2 1    2 4 2    1 2 1  | range
    //   1 1 1    2 2 2    1 1 1  ⌄
    //   <---------time--------> (previous scan, current scan, next scan)
    let mut weights: [[[i32; FILTER_DEPTH as usize]; FILTER_HEIGHT as usize]; FILTER_WIDTH as usize] = [];
    let mut f: i32;
    let mut w: i32 = 1;
    for z in 0..FILTER_DEPTH {
        if z == 1 {
            f = 2;
        } else {
            f = 1;
        }
        for y in 0..FILTER_HEIGHT {
            for x in 0..FILTER_WIDTH {
                if x == 1 && y == 1 {
                    w = 2;
                } else {
                    w = 1;
                }
                weights[x][y][z] = w * f;
            }
        }
    }

    // Add enough beams and ranges to the output RadarScan
    for beam in 0..max_beam {
        out_scan.add_beam(max_range);
        out_scan.beams[beam].beam = -1;
    }
    Ok(out_scan)
}
