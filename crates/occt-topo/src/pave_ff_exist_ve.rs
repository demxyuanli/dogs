//! Second `ProcessExistingPaveBlocks` overload — bound vertices vs ON/IN
//! pave blocks (`ComputeVE`).
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx:3171-3276`.
//!
//! For each bound vertex produced by `PutBoundPaveOnCurve`, scan ON/IN pave
//! blocks whose boxes overlap the vertex box (AABB stand-in for BVH). Skip
//! blocks that already use the vertex as an end, skip blocks already queued
//! for post-treat, then `ComputePE` (OCCT `ComputeVE`). A hit (`flag == 0`)
//! queues the block via `PreparePostTreatFF` and records the opposite face
//! on `thePBFacesMap` when the block is not already On/In both faces.

use std::collections::{HashMap, HashSet};

use crate::bopds::BopdsPaveBlock;

use crate::int_tools_full::IntToolsContext;
use crate::pave_ff::CoupleOfPaveBlocks;
use crate::pave_ff_exist::{pb_in_face, pb_key, prepare_post_treat_ff, PbKey};
use crate::pave_ff_update::PbFacesMap;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, TopoShape, Vertex};

fn edge_box_out_vertex(ds: &crate::bopds::BopdsDS, n_e: usize, n_v: usize) -> bool {
    let (Some(be), Some(bv)) = (ds.box_of(n_e), ds.box_of(n_v)) else {
        return false;
    };
    be.is_out_box(bv)
}

/// Bound-vertex overload of `ProcessExistingPaveBlocks` (`_6.cxx:3171`).
pub fn process_existing_pave_blocks_bound_ve(
    f: &mut PaveFiller,
    the_int: usize,
    n_f1: usize,
    n_f2: usize,
    on_in: &[BopdsPaveBlock],
    dmbv: &HashMap<usize, Vec<usize>>,
    on1: &[(usize, f64, f64)],
    in1: &[(usize, f64, f64)],
    on2: &[(usize, f64, f64)],
    in2: &[(usize, f64, f64)],
    mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    mpb_add: &mut HashSet<PbKey>,
    pb_faces: &mut PbFacesMap,
) {
    if dmbv.is_empty() {
        return;
    }
    let fuzzy = f.fuzzy_value();
    let ctx = IntToolsContext::new();
    let keys: Vec<(usize, Vec<usize>)> = dmbv.iter().map(|(&k, v)| (k, v.clone())).collect();
    for (i_c, a_lbv) in keys {
        let mut kept_extra: Vec<BopdsPaveBlock> = Vec::new();
        for n_v in a_lbv {
            let Some(v_shape) = f.ds().shape(n_v).cloned() else {
                continue;
            };
            let vertex = Vertex(v_shape);
            for a_pb in on_in {
                if a_pb.edge() == 0 {
                    continue;
                }
                let (n1, n2) = a_pb.indices();
                if n1 == n_v || n2 == n_v {
                    continue;
                }
                if mpb_add.contains(&pb_key(a_pb)) {
                    continue;
                }
                if edge_box_out_vertex(f.ds(), a_pb.edge(), n_v) {
                    continue;
                }
                let Some(e_shape) = f.ds().shape(a_pb.edge()).cloned() else {
                    continue;
                };
                let edge = Edge(e_shape);
                let flag = ctx.compute_pe(&vertex, &edge, fuzzy);
                if flag != 0 {
                    continue;
                }
                mpb_add.insert(pb_key(a_pb));
                prepare_post_treat_ff(f.ds(), the_int, i_c, a_pb, mscpb, &mut kept_extra);
                let b_in_f1 = pb_in_face(on1, in1, a_pb);
                let b_in_f2 = pb_in_face(on2, in2, a_pb);
                if !b_in_f1 || !b_in_f2 {
                    let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                    let faces = pb_faces.entry(pb_key(a_pb)).or_default();
                    if !faces.contains(&n_f) {
                        faces.push(n_f);
                    }
                }
            }
        }
        if !kept_extra.is_empty() {
            let lpbc = f.ds_mut().interf_ff_mut()[the_int]
                .change_curves()
                .get_mut(i_c);
            if let Some(list) = lpbc {
                list.change_pave_blocks().extend(kept_extra);
            }
        }
    }
}

/// Distance lookup used by the first ProcessExisting overload (`_6.cxx:3121`).
pub fn range_distance_on_overlap(
    distances: &HashMap<(usize, usize), Vec<crate::pave_filler::EdgeRangeDistance>>,
    n_e: usize,
    n_f: usize,
    t1: f64,
    t2: f64,
) -> Option<f64> {
    let list = distances.get(&(n_e, n_f)).or_else(|| distances.get(&(n_f, n_e)))?;
    for r in list {
        if (t1 <= r.first && r.first <= t2)
            || (t1 <= r.last && r.last <= t2)
            || (r.first <= t1 && t1 <= r.last)
            || (r.first <= t2 && t2 <= r.last)
        {
            return Some(r.distance);
        }
    }
    None
}
