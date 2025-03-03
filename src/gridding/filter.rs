use crate::error::BackscatterError;
use crate::utils::scan::{RadarBeam, RadarCell, RadarScan};

pub const MAX_BEAM: i32 = 256;
pub const FILTER_HEIGHT: i32 = 3;
pub const FILTER_WIDTH: i32 = 3;
pub const FILTER_DEPTH: i32 = 3;

/// Calculates the mean and standard deviation of a parameter from the vector `v`.
/// `f` is used to extract the parameter from an entry of `v`.
fn calculate_mean_sigma(v: &Vec<&RadarCell>, f: fn(&RadarCell) -> f64) -> (f64, f64) {
    let mut mean = 0.0;
    let mut variance = 0.0;

    // Calculate the mean value
    for &cell in v.iter() {
        mean += f(cell);
    }
    mean /= v.len();

    // Calculate the variance of the velocity values
    for &cell in v.iter() {
        variance += (f(cell) - mean) * (f(cell) - mean);
    }
    variance /= v.len();
    let sigma = variance.sqrt();

    (mean, sigma)
}

/// Calculates the median value of RadarCells in kernel, and the standard deviation
/// of those cells which are within two standard deviations of the mean of cells in kernel.
/// The parameter `f` is a function which extracts the parameter from an entry of kernel.
/// The parameter `g` is used to extract the parameter for sorting (which may be different than the
/// parameter having its median value calculated)
/// Returns the median and standard deviation.
fn calculate_median_sigma(
    kernel: &mut Vec<&RadarCell>,
    f: fn(&RadarCell) -> f64,
    g: fn(&RadarCell) -> f64,
) -> (f64, f64) {
    // Calculate mean and std deviation of kernel with respect to lambda power
    let (mean, sigma) = calculate_mean_sigma(&kernel, f);

    // Only keep values which fall within 2 std deviations of mean
    let mut valid_cells: Vec<&RadarCell> = vec![];
    for &cell in kernel.iter() {
        // If the cell deviates by more than 2 standard deviations from the mean, skip it
        if (f(cell) - mean).abs() > 2.0 * sigma {
            continue;
        }
        // Add the cell to the median structure
        valid_cells.push(cell);
    }
    // Sort cells in median by their value according to `g`
    valid_cells.sort_by(|&a, &b| g(a).partial_cmp(&g(b)).unwrap());

    // Set the beam/range cell in out_scan to the median value
    let median = f(valid_cells[valid_cells.len() / 2]);

    // Recalculate the standard deviation, this time including only valid cells
    let (_, sigma) = calculate_mean_sigma(&valid_cells, f);

    (median, sigma)
}

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
    isort: bool,
    scans: &[&RadarScan],
) -> Result<RadarScan, BackscatterError> {
    let mut out_scan = RadarScan {
        ..Default::default()
    };
    let mut max_beam: i32 = -1;
    let mut max_range: i32 = 1000;
    let threshold = &[12, 24];
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
                max_beam = beam.beam + 1; // Add one since beam number is indexed from 0
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
    let mut weights: [[[i32; FILTER_DEPTH as usize]; FILTER_HEIGHT as usize];
        FILTER_WIDTH as usize] = [];
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

    // [max_beams, depth, num_points] to store all observations grouped by beam number
    let mut beam_pointers: Vec<Vec<Vec<Option<RadarBeam>>>> = Vec::with_capacity(max_beam as usize);

    // The largest amount of observations for a given beam, for any scan
    let mut max_observations_for_a_beam: i32 = 0;

    // Add enough beams and ranges to the output RadarScan
    for beam in 0..max_beam {
        out_scan.add_beam(max_range);
        out_scan.beams[beam].beam = -1;

        // Initialize some vectors for storing observations along a beam direction
        for _ in 0..depth {
            // Adding an empty vector for this depth, containing an empty vector for the points
            beam_pointers.push(vec![vec![]])
        }
        // beam_pointers should now be [max_beams, depth, 0], where the last dimension is an empty Vec
    }

    for z in 0..depth {
        // Figure out if this scan is the current, previous, or next scan
        let mut i = index - (depth - 1) + z;
        if i < 0 {
            i += depth;
        }

        // Loop through the beams in this scan
        for beam in scans[i].beams.iter() {
            let beam_num = beam.beam;
            beam_pointers[beam_num][depth].push(beam);

            // Update the largest amount of observations seen
            if beam_pointers[beam_num][depth].len() > max_observations_for_a_beam {
                max_observations_for_a_beam = beam_pointers[beam_num][depth].len();
            }
        }
    }

    // Get index of center scan in temporal dimension
    let mut i = index - 1;
    if i < 0 {
        i += depth;
    }

    // Copy over parameters from center scan
    out_scan.station_id = scans[i].station_id;
    out_scan.version_major = scans[i].version_major;
    out_scan.version_minor = scans[i].version_minor;
    out_scan.start_time = scans[i].start_time;
    out_scan.end_time = scans[i].end_time;

    // If mode is a multiple of
    if mode % 4 == 0 {
        for beam_num in 0..max_beam {
            // If center scan doesn't have beams then skip this beam
            if beam_pointers[beam_num][depth / 2].len() == 0 {
                continue;
            }
            let beam = &beam_pointers[beam_num][depth / 2][0]; // First beam
            let b = &out_scan.beams[beam_num];

            // Copy radar operating parameters from first beam into corresponding beam of out_scan
            b.beam = beam_num;
            b.program_id = beam.program_id;
            b.time = beam.time;
            b.integration_time_s = beam.integration_time_s;
            b.integration_time_us = beam.integration_time_us;
            b.num_averages = beam.num_averages;
            b.first_range = beam.first_range;
            b.range_sep = beam.range_sep;
            b.rx_rise = beam.rx_rise;
            b.freq = beam.freq;
            b.noise = beam.noise;
            b.attenuation = beam.attenuation;
            b.channel = beam.channel;
            b.num_ranges = beam.num_ranges;
        }
    } else {
        for beam_num in 0..max_beam {
            let b = &mut out_scan.beams[beam_num];

            // Initialize radar operating parameters
            b.program_id = -1;
            b.time = 0;
            b.integration_time_s = 0;
            b.integration_time_us = 0;
            b.first_range = 0;
            b.range_sep = 0;
            b.rx_rise = 0;
            b.freq = 0;
            b.noise = 0;
            b.attenuation = 0;
            b.channel = -1;
            b.num_ranges = -1;
        }

        for z in 0..depth {
            for beam_num in 0..max_beam {
                // If no beams previously found, continue
                if beam_pointers[beam_num][z].len() == 0 {
                    continue;
                }

                // Corresponding beam in out_scan
                let mut out_beam = &mut out_scan.beams[beam_num];

                // Setting beam number in out_scan for this beam
                out_beam.beam = beam_num;

                // Go through all beams for this beam/time combo
                for in_beam in beam_pointers[beam_num][z].iter() {
                    // If this is the first beam then use it to set program_id for out_scan beam
                    if out_beam.program_id == -1 {
                        out_beam.program_id = in_beam.program_id
                    }

                    // Sum all the operating parameters, which will be averaged later once all beams
                    // have been added
                    out_beam.time += in_beam.time;
                    out_beam.integration_time_s += in_beam.integration_time_s;
                    out_beam.integration_time_us += in_beam.integration_time_us;
                    if out_beam.integration_time_us > 1_000_000 {
                        out_beam.integration_time_us += 1;
                        out_beam.integration_time_us -= 1_000_000;
                    }
                    out_beam.num_averages += in_beam.num_averages;
                    out_beam.first_range += in_beam.first_range;
                    out_beam.range_sep += in_beam.range_sep;
                    out_beam.rx_rise += in_beam.rx_rise;
                    out_beam.freq += in_beam.freq;
                    out_beam.noise += in_beam.noise;
                    out_beam.attenuation += in_beam.attenuation;
                    if out_beam.channel == 0 {
                        out_beam.channel = in_beam.channel;
                    }

                    // If this is the first beam in the time/beam combo then use max_range
                    // to set the number of range gates for the beam
                    if out_beam.num_ranges == -1 {
                        out_beam.num_ranges = max_range;
                    }
                }
            }
        }

        for beam_num in 0..max_beam {
            let mut count = 0;

            // Count all the observations for this beam, summing over all scans being averaged
            for z in 0..depth {
                count += beam_pointers[beam_num][z].len();
            }

            // Corresponding beam in out_scan
            let out_beam = &mut out_scan.beams[beam_num];

            out_beam.time /= count;
            out_beam.num_averages /= count;
            out_beam.first_range /= count;
            out_beam.range_sep /= count;
            out_beam.rx_rise /= count;
            out_beam.freq /= count;
            out_beam.noise /= count;
            out_beam.attenuation /= count;
            out_beam.integration_time_us /= count;
            let mut microseconds = (out_beam.integration_time_s * 1_000_000) / count;
            out_beam.integration_time_s /= count;
            microseconds -= out_beam.integration_time_s * 1_000_000;
            out_beam.integration_time_us += microseconds;
        }
    }

    // 3 x 3 x 3 kernel for storing all values of data for median filtering
    let mut kernel = vec![];

    for beam_num in 0..max_beam {
        for range in 0..max_range {
            // Set up the spatial 3x3 (beam by range) filtering boundaries
            let mut bmin = beam_num - FILTER_WIDTH / 2;
            let bbox = beam_num - FILTER_WIDTH / 2;
            let mut bmax = beam_num + FILTER_WIDTH / 2;
            let mut rmin = range - FILTER_HEIGHT / 2;
            let rbox = range - FILTER_HEIGHT / 2;
            let mut rmax = range + FILTER_HEIGHT / 2;

            // Set lower beam boundary to 0 when at edge of FOV
            if bmin < 0 {
                bmin = 0;
            }
            // Set upper beam boundary to highest beam when at other edge of FOV
            if bmax >= max_beam {
                bmax = max_beam - 1;
            }
            // Set lower range boundary to 0 when at nearest edge of FOV
            if rmin < 0 {
                rmin = 0;
            }
            // Set upper range boundary to furthest range gate when at other edge of FOV
            if rmax >= max_range {
                rmax = max_range - 1;
            }

            // Initialize center cell weight to zero
            let mut weight = 0;

            // Loop over beams
            for x in bmin..bmax {
                // Loop over ranges
                for y in rmin..rmax {
                    // Loop over time
                    for z in 0..depth {
                        // Loop over beams in time/beam combo
                        for beam in beam_pointers[x][z].iter() {
                            // Skip if this range gate is not in the beam
                            if y >= beam.scatter.len() {
                                continue;
                            }

                            // Check that there is scatter present in the beam/range/time cell
                            if beam.scatter[y] != 0 {
                                // Increment weight
                                weight += weights[x - bbox][y - rbox][z];
                                // Add this observation to the kernel
                                kernel.push(&beam.cells[y]);
                            }
                        }
                    }
                }
            }
            // If no cells with scatter found, continue
            if kernel.len() == 0 {
                continue;
            }

            // If the current beam is at the edge of the FOV then increase its weight by 50%
            // TODO: What about near/far range edges?
            // TODO: weight is an integer, this is kinda hacky
            if beam_num == 0 || beam_num == max_beam - 1 {
                weight = weight * 1.5;
            }

            // If the sum of weights of cells with scatter in the kernel is less than the threshold
            // then continue
            if weight <= threshold[mode % 2] {
                continue;
            }

            // Threshold was exceeded, so the output scan should have scatter in this beam/range cell
            let out_beam = &mut out_scan.beams[beam_num];
            out_beam.scatter[range] = 1;

            // Initialize observation parameters to zero
            let out_cell = &mut out_beam.cells[range];
            out_cell.groundscatter = 0;
            out_cell.power_lin = 0.0;
            out_cell.spectral_width_lin = 0.0;
            out_cell.velocity = 0.0;

            // TODO: Figure out how to properly check param (RST does bitwise checks)
            // Perform velocity median filtering if specified
            let mut compare_fn: fn(&RadarCell) -> f64 = |x| x.velocity;
            if param % 2 == 1 {
                (out_cell.velocity, out_cell.velocity_error) =
                    calculate_median_sigma(&kernel, |x| x.velocity, compare_fn);
            }

            // Perform lambda power median filtering if specified
            if param % 2 == 0 {
                if isort == true {
                    compare_fn = |x| x.power_lin;
                }
                (out_cell.power_lin, out_cell.power_lin_error) =
                    calculate_median_sigma(&kernel, |x| x.power_lin, compare_fn);
            }

            // Perform spectral width median filtering if specified
            if param % 4 == 0 {
                if isort == true {
                    compare_fn = |x| x.spectral_width_lin;
                }
                (
                    out_cell.spectral_width_lin,
                    out_cell.spectral_width_lin_error,
                ) = calculate_median_sigma(&kernel, |x| x.spectral_width_lin, compare_fn);
            }

            // Perform lag0 power median filtering if specified
            if param % 8 == 0 {
                if isort == true {
                    compare_fn = |x| x.power_lag_zero;
                }
                (out_cell.power_lag_zero, out_cell.power_error_lag_zero) =
                    calculate_median_sigma(&kernel, |x| x.power_lag_zero, compare_fn);
            }
        }
    }

    Ok(out_scan)
}
