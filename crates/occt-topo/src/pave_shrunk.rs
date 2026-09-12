//! `BOPAlgo_PaveFiller::FillShrunkData(type1, type2)`
//! (`BOPAlgo_PaveFiller_9.cxx:65`).
//!
//! Collects pave blocks of the EDGE side of every interfering pair of the
//! given types, skips flagged (degenerated) edges and already-valid shrunk
//! data, then runs `IntTools_ShrunkRange` + `AnalyzeShrunkData`.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::pave_common::fill_shrunk_data_for_block;
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::collect_pairs;
use crate::shape::{Edge, Vertex};

/// `BOPDS_DS::IsValidShrunkData` reduced to the presence of shrunk data plus
/// a distance check of the shrunk ends against the bound vertices.
pub fn is_valid_shrunk_data(f: &PaveFiller, n_e: usize, pb: &crate::bopds::BopdsPaveBlock) -> bool {
    if !pb.has_shrunk_data() {
        return false;
    }
    let Some(es) = f.ds().shape(n_e) else {
        return false;
    };
    let edge = Edge(es.clone());
    let Some(curve) = BRepTool::edge_curve(&edge) else {
        return false;
    };
    let (ts1, ts2, _) = pb.shrunk_data();
    let (n_v1, n_v2) = pb.indices();
    let eps = BRepTool::edge_tolerance(&edge) * 0.01;
    for (t, n_v) in [(ts1, n_v1), (ts2, n_v2)] {
        let Some(vs) = f.ds().shape(n_v) else {
            return false;
        };
        let v = Vertex(vs.clone());
        let tol = BRepTool::vertex_tolerance(&v) + occt_core::precision::CONFUSION;
        let dist = BRepTool::vertex_point(&v).distance(&curve.d0(t));
        if tol - dist > eps {
            return false;
        }
    }
    true
}

/// `BOPAlgo_PaveFiller::FillShrunkData(aType1, aType2)`.
pub fn fill_shrunk_data_typed(
    f: &mut PaveFiller,
    a_type1: ShapeType,
    a_type2: ShapeType,
) -> Result<(), String> {
    let pairs = collect_pairs(f.ds(), a_type1, a_type2);
    if pairs.is_empty() {
        return Ok(());
    }
    let a_type = [a_type1, a_type2];
    let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut jobs: Vec<(usize, usize)> = Vec::new();
    for (n0, n1) in pairs {
        let n_s = [n0, n1];
        for i in 0..2 {
            let n_e = n_s[i];
            if a_type[i] != ShapeType::Edge || !a_mi.insert(n_e) {
                continue;
            }
            let Some(si_e) = f.ds().shape_info(n_e).cloned() else {
                continue;
            };
            if si_e.has_flag() {
                continue;
            }
            let n_pb = f.ds().pave_blocks(n_e).len();
            for slot in 0..n_pb {
                let pb = f.ds().pave_blocks(n_e)[slot].clone();
                if pb.has_shrunk_data() && is_valid_shrunk_data(f, n_e, &pb) {
                    continue;
                }
                jobs.push((n_e, slot));
            }
        }
    }
    let tol = f.fuzzy_value();
    for (n_e, slot) in jobs {
        let Some(es) = f.ds().shape(n_e).cloned() else {
            continue;
        };
        let a_e = Edge(es);
        let mut pb = f.ds().pave_blocks(n_e)[slot].clone();
        fill_shrunk_data_for_block(f, &a_e, tol, &mut pb);
        if let Some(slot_mut) = f.ds_mut().change_pave_blocks_mut(n_e).get_mut(slot) {
            *slot_mut = pb;
        }
    }
    Ok(())
}

/// `FillShrunkData(VERTEX, EDGE)` used by `PerformVE`.
pub fn fill_shrunk_data_ve(f: &mut PaveFiller) -> Result<(), String> {
    fill_shrunk_data_typed(f, ShapeType::Vertex, ShapeType::Edge)
}

/// `FillShrunkData(EDGE, EDGE)` used by `PerformEE`.
pub fn fill_shrunk_data_ee(f: &mut PaveFiller) -> Result<(), String> {
    fill_shrunk_data_typed(f, ShapeType::Edge, ShapeType::Edge)
}

/// `FillShrunkData(EDGE, FACE)` used by `PerformEF`.
pub fn fill_shrunk_data_ef(f: &mut PaveFiller) -> Result<(), String> {
    fill_shrunk_data_typed(f, ShapeType::Edge, ShapeType::Face)
}

/// `FillShrunkData(VERTEX, FACE)` used by `PerformVF` when shrinking
/// incident edges of the interfering vertices.
pub fn fill_shrunk_data_vf(f: &mut PaveFiller) -> Result<(), String> {
    fill_shrunk_data_typed(f, ShapeType::Vertex, ShapeType::Face)
}


/// `BOPAlgo_PaveFiller::AnalyzeShrunkData` (`BOPAlgo_PaveFiller_3.cxx:766`).
///
/// When the shrunk range failed or is not splittable, the warning shape is
/// the whole edge (if the pave block covers it) or a compound of the edge
/// plus its two bound vertices (`AlertTooSmallEdge` /
/// `AlertNotSplittableEdge` / `AlertBadPositioning`).
pub fn analyze_shrunk_data_full(
    f: &mut PaveFiller,
    pb: &mut crate::bopds::BopdsPaveBlock,
    sr: &crate::inttools_range::ShrunkRange,
    edge: &Edge,
) {
    let (a_e_first, a_e_last) = BRepTool::edge_parameters(edge);
    let (a_pb_first, a_pb_last) = pb.range();
    let b_whole_edge = a_pb_first <= a_e_first && a_pb_last >= a_e_last;
    let done = sr.is_done();
    let splittable = sr.is_splittable();
    if !done || !splittable {
        let warn = if b_whole_edge && pb.original_edge() != usize::MAX {
            edge.0.clone()
        } else {
            let (n_v1, n_v2) = pb.indices();
            let mut parts = vec![edge.0.clone()];
            if let Some(v1) = f.ds().shape(n_v1) {
                parts.push(v1.clone());
            }
            if let Some(v2) = f.ds().shape(n_v2) {
                parts.push(v2.clone());
            }
            crate::builder::TopoBuilder::new().make_compound_of(&parts).0
        };
        if !done {
            if b_whole_edge {
                f.add_warning(format!(
                    "too small edge: shrunk range cannot be computed ({})",
                    crate::tgeometry::GeometryRegistry::shape_key(&warn)
                ));
            } else {
                f.add_warning(format!(
                    "bad positioning: shrunk range cannot be computed ({})",
                    crate::tgeometry::GeometryRegistry::shape_key(&warn)
                ));
            }
            pb.clear_shrunk_data();
            return;
        }
        if b_whole_edge {
            f.add_warning(format!(
                "not splittable edge: shrunk range is too short ({})",
                crate::tgeometry::GeometryRegistry::shape_key(&warn)
            ));
        } else {
            f.add_warning(format!(
                "bad positioning: shrunk range is too short ({})",
                crate::tgeometry::GeometryRegistry::shape_key(&warn)
            ));
        }
    }
    match sr.shrunk_range() {
        Some((ts1, ts2)) => pb.set_shrunk_data(ts1, ts2, sr.is_splittable()),
        None => pb.set_shrunk_data(0.0, 0.0, false),
    }
}

/// `BOPAlgo_PaveFiller::AddIntersectionFailedWarning`.
pub fn add_intersection_failed_warning(f: &mut PaveFiller, s1: &crate::shape::TopoShape, s2: &crate::shape::TopoShape) {
    let wc = crate::builder::TopoBuilder::new().make_compound_of(&[s1.clone(), s2.clone()]);
    f.add_warning(format!(
        "intersection failed for sub-shapes ({})",
        crate::tgeometry::GeometryRegistry::shape_key(&wc.0)
    ));
}

