//! First `ProcessExistingPaveBlocks` overload — section edge vs ON/IN
//! pave blocks (`_6.cxx:3072`).
//!
//! AABB stand-in for `BOPTools_BoxTree`. A block already On/In both faces is
//! queued for post-treat. Otherwise `myDistances` is consulted for the
//! original-edge / opposite-face pair; an overlapping range whose distance
//! is within `tol(ES) + tol(EF)` queues the block and records the opposite
//! face on `thePBFacesMap`.

use std::collections::{HashMap, HashSet};

use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::BopdsPaveBlock;
use crate::brep_tool::BRepTool;
use crate::pave_ff::CoupleOfPaveBlocks;
use crate::pave_ff_exist::{pb_in_face, pb_key, prepare_post_treat_ff, PbKey};
use crate::pave_ff_is_exist::ranges_overlap;
use crate::pave_ff_update::PbFacesMap;
use crate::pave_filler::{EdgeRangeDistance, PaveFiller};
use crate::shape::{Edge, TopoShape};

fn edge_box_out(ds: &crate::bopds::BopdsDS, n_e: usize, other: &occt_core::bnd::BndBox) -> bool {
    ds.box_of(n_e).map(|b| b.is_out_box(other)).unwrap_or(false)
}

fn distance_on_overlap(
    distances: &HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    n_e: usize,
    n_f: usize,
    t1: f64,
    t2: f64,
) -> Option<f64> {
    let list = distances
        .get(&(n_e, n_f))
        .or_else(|| distances.get(&(n_f, n_e)))?;
    for r in list {
        if ranges_overlap(t1, t2, r.first, r.last) {
            return Some(r.distance);
        }
    }
    None
}

/// First `ProcessExistingPaveBlocks` overload (`_6.cxx:3072`).
pub fn process_existing_pave_blocks_es(
    f: &PaveFiller,
    the_int: usize,
    the_cur: usize,
    n_f1: usize,
    n_f2: usize,
    the_es: &TopoShape,
    on_in: &[BopdsPaveBlock],
    on1: &[(usize, f64, f64)],
    in1: &[(usize, f64, f64)],
    on2: &[(usize, f64, f64)],
    in2: &[(usize, f64, f64)],
    distances: &HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    lpbc: &mut Vec<BopdsPaveBlock>,
    mpb_add: &mut HashSet<PbKey>,
    pb_faces: &mut PbFacesMap,
) {
    let box_es = shape_bbox(the_es);
    let a_tol_es = BRepTool::edge_tolerance(&Edge(the_es.clone()));
    for a_pbf in on_in {
        if a_pbf.edge() == 0 {
            continue;
        }
        if edge_box_out(f.ds(), a_pbf.edge(), &box_es) {
            continue;
        }
        if mpb_add.contains(&pb_key(a_pbf)) {
            continue;
        }
        let b_in_f1 = pb_in_face(on1, in1, a_pbf);
        let b_in_f2 = pb_in_face(on2, in2, a_pbf);
        if b_in_f1 && b_in_f2 {
            mpb_add.insert(pb_key(a_pbf));
            prepare_post_treat_ff(f.ds(), the_int, the_cur, a_pbf, mscpb, lpbc);
            continue;
        }
        let n_f = if b_in_f1 { n_f2 } else { n_f1 };
        let (t1, t2) = a_pbf.range();
        let Some(a_dist) = distance_on_overlap(distances, a_pbf.original_edge(), n_f, t1, t2) else {
            continue;
        };
        let Some(ef_shape) = f.ds().shape(a_pbf.edge()) else {
            continue;
        };
        let a_tol_sum = a_tol_es + BRepTool::edge_tolerance(&Edge(ef_shape.clone()));
        if a_dist <= a_tol_sum {
            mpb_add.insert(pb_key(a_pbf));
            prepare_post_treat_ff(f.ds(), the_int, the_cur, a_pbf, mscpb, lpbc);
            let faces = pb_faces.entry(pb_key(a_pbf)).or_default();
            if !faces.contains(&n_f) {
                faces.push(n_f);
            }
        }
    }
}

/// Wire the first overload from `MakeBlocks` after a new section edge is born.
pub fn process_existing_after_make_edge(
    f: &PaveFiller,
    the_int: usize,
    the_cur: usize,
    n_f1: usize,
    n_f2: usize,
    the_es: &TopoShape,
    on_in: &[BopdsPaveBlock],
    on1: &[(usize, f64, f64)],
    in1: &[(usize, f64, f64)],
    on2: &[(usize, f64, f64)],
    in2: &[(usize, f64, f64)],
    distances: &HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    lpbc: &mut Vec<BopdsPaveBlock>,
    mpb_add: &mut HashSet<PbKey>,
    pb_faces: &mut PbFacesMap,
) {
    process_existing_pave_blocks_es(
        f, the_int, the_cur, n_f1, n_f2, the_es, on_in, on1, in1, on2, in2, distances, mscpb,
        lpbc, mpb_add, pb_faces,
    );
}
