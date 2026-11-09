//! `LocOpe::SampleEdges` (`LocOpe.cxx:195-236`).
//!
//! Collects interior points of every non-degenerated edge, then every vertex,
//! in `TopExp_Explorer` order, skipping repeats through a
//! `TopTools_ShapeMapHasher` set (edges and vertices share one set).

use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::locope::shape_maps::ShapeSet;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::explorer_shapes;

/// `NECHANT` (`LocOpe.cxx:39`).
const NECHANT: i32 = 10;

/// `LocOpe::SampleEdges(theShape, theSeq)` (`cxx:195-236`).
///
/// OCCT clears `theSeq` first and appends in place; the port returns the
/// collected sequence instead (`cxx:197`).
pub fn sample_edges(shape: &TopoShape) -> Vec<GpPnt> {
    let mut the_seq: Vec<GpPnt> = Vec::new();
    let mut the_map = ShapeSet::new();

    // `cxx:208-226`: interior points, extremities excluded.
    for sub in explorer_shapes(shape, ShapeType::Edge) {
        // `cxx:209-212`: `if (!theMap.Add(edg)) continue;`.
        if the_map.contains(&sub) {
            continue;
        }
        the_map.add(&sub);
        let Some(edg) = Edge::wrap(sub) else {
            continue;
        };
        // `cxx:213-224`.
        if !BRepTool::is_degenerated(&edg) {
            let Some(curve) = BRepTool::edge_curve_world(&edg) else {
                continue;
            };
            let (f, l) = BRepTool::edge_parameters(&edg);
            let delta = (l - f) / NECHANT as f64 * 0.123456;
            for i in 1..NECHANT {
                let prm = delta
                    + ((NECHANT - i) as f64 * f + i as f64 * l) / NECHANT as f64;
                the_seq.push(curve.d0(prm));
            }
        }
    }

    // `cxx:228-234`: every vertex not already counted as an edge.
    for sub in explorer_shapes(shape, ShapeType::Vertex) {
        // `cxx:230-233`.
        if the_map.contains(&sub) {
            continue;
        }
        the_map.add(&sub);
        if let Some(v) = Vertex::wrap(sub) {
            the_seq.push(BRepTool::vertex_point_world(&v));
        }
    }
    the_seq
}
