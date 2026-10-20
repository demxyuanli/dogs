//! Round-trip count helper.

use super::super::prelude::*;
use super::*;

/// Round-trip helper: write `shape` to STEP, read it back, and count the
/// distinct vertices / edges / faces in the reconstructed model.
pub fn step_roundtrip_counts(shape: &TopoShape) -> Result<(usize, usize, usize), String> {
    let step = write_shape_step(shape);
    let model = read_step(&step)?;
    let mut nv = 0usize;
    let mut ne = 0usize;
    let mut nf = 0usize;
    for ms in &model.shapes {
        let c = crate::topo_tools::shape_counts(&ms.shape);
        nv += c.get(&ShapeType::Vertex).copied().unwrap_or(0);
        ne += c.get(&ShapeType::Edge).copied().unwrap_or(0);
        nf += c.get(&ShapeType::Face).copied().unwrap_or(0);
    }
    Ok((nv, ne, nf))
}
