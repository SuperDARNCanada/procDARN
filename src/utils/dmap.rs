use dmap::{DmapVec, InDmap};

pub fn convert_to_dmapvec<T: InDmap>(vals: Vec<T>) -> DmapVec<T> {
    DmapVec {
        dimensions: vec![vals.len() as i32],
        data: vals,
    }
}
