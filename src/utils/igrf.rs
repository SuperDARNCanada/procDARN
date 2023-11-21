use crate::error::BackscatterError;
use geodesy::Coor4D;

pub fn igrf_magnetic_components(date: f64, loc: &Coor4D) -> Result<Coor4D, BackscatterError> {
    let b_field = igrf_call(date, loc)?;

    // Convert to local south/vertical (rather than north/down)
    b_field[0] *= -1;
    b_field[2] *= -1;

    b_field
}
