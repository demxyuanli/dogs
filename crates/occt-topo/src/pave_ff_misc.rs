//! Remaining `BOPAlgo_PaveFiller_6.cxx` helpers used by MakeBlocks.
//!
//! Source: `GetStickVertices` at 2847, `GetFullShapeMap` at 2909,
//! `RemoveUsedVertices` at 2928, `EstimatePaveOnCurve` at 4056,
//! `CheckPlanes` at 3639, `UpdateBlocksWithSharedVertices` at 3946.

use std::collections::HashSet;

use crate::bopds::BopdsDS;
use crate::bopds_ff::BopdsCurve;


use crate::pave_filler::PaveFiller;



/// `BOPAlgo_PaveFiller::GetFullShapeMap`.
pub fn get_full_shape_map(ds: &BopdsDS, n_f: usize, mi: &mut HashSet<usize>) {
    mi.insert(n_f);
    if let Some(si) = ds.shape_info(n_f) {
        for &n in si.sub_shapes() {
            mi.insert(n);
        }
    }
}

/// `BOPAlgo_PaveFiller::GetStickVertices`.
pub fn get_stick_vertices(
    ds: &BopdsDS,
    n_f1: usize,
    n_f2: usize,
) -> (HashSet<usize>, HashSet<usize>, HashSet<usize>) {
    let mut mi = HashSet::new();
    get_full_shape_map(ds, n_f1, &mut mi);
    get_full_shape_map(ds, n_f2, &mut mi);
    let mut mv_stick = HashSet::new();
    let mut mv_ef = HashSet::new();
    let push_new = |arr: &[crate::bopds::BopdsInterf], mi: &HashSet<usize>, out: &mut HashSet<usize>| {
        for it in arr {
            let Some(n_new) = it.get_index_new() else {
                continue;
            };
            let (s1, s2) = it.indices();
            if mi.contains(&s1) && mi.contains(&s2) {
                out.insert(ds.get_same_domain_index(n_new));
            }
        }
    };
    push_new(ds.interf_vv(), &mi, &mut mv_stick);
    push_new(ds.interf_ve(), &mi, &mut mv_stick);
    push_new(ds.interf_ee(), &mi, &mut mv_stick);
    push_new(ds.interf_vf(), &mi, &mut mv_stick);
    for it in ds.interf_ef() {
        let Some(n_new) = it.get_index_new() else {
            continue;
        };
        let (s1, s2) = it.indices();
        if mi.contains(&s1) && mi.contains(&s2) {
            let n = ds.get_same_domain_index(n_new);
            mv_stick.insert(n);
            mv_ef.insert(n);
        }
    }
    (mv_stick, mv_ef, mi)
}

/// `BOPAlgo_PaveFiller::RemoveUsedVertices`.
pub fn remove_used_vertices(a_vc: &[BopdsCurve], mv: &mut HashSet<usize>) {
    for nc in a_vc {
        for pb in nc.pave_blocks() {
            mv.remove(&pb.index1);
            mv.remove(&pb.index2);
            for p in pb.ext_paves() {
                mv.remove(&p.index);
            }
        }
        if let Some(pb0) = nc.pave_blocks().first() {
            for p in pb0.ext_paves() {
                mv.remove(&p.index);
            }
        }
    }
}

/// `BOPAlgo_PaveFiller::EstimatePaveOnCurve`.
pub fn estimate_pave_on_curve(ds: &BopdsDS, n_v: usize, nc: &BopdsCurve, tol: f64) -> bool {
    crate::pave_ff_pave_stick::estimate_pave_on_curve(ds, n_v, nc, tol)
}

/// `BOPAlgo_PaveFiller::CheckPlanes`.
pub fn check_planes(ds: &BopdsDS, n_f1: usize, n_f2: usize) -> bool {
    let verts = |n_f: usize| -> HashSet<usize> {
        let mut s = HashSet::new();
        if let Some(fi) = ds.face_info(n_f) {
            for &v in fi.verts_on() {
                s.insert(ds.get_same_domain_index(v));
            }
            for &v in fi.verts_in() {
                s.insert(ds.get_same_domain_index(v));
            }
        }
        s
    };
    let a = verts(n_f1);
    let b = verts(n_f2);
    a.intersection(&b).count() > 1
}

/// `BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices`.
pub fn update_blocks_with_shared_vertices(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ff_shared::update_blocks_with_shared_vertices(f)
}
