//! `UpdateFaceInfo` remainder: replace On/In vertices and On/In/Sc pave
//! blocks after PostTreatFF produced SD vertices and replacement edges.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx:1673-1946` (step 2, `bVerts` / `bEdges`).
//! Step 1 (section edges/vertices into FaceInfo Sc) stays in
//! `pave_ff_update::update_face_info_ff`. This module walks every face that
//! participated in an F/F pair and substitutes `theDMV` vertices plus
//! `theDME` pave-block replacements, using `RealPaveBlock` and a fence so
//! each real block is added once.

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::bopds::BopdsPaveBlock;
use crate::pave_ff_exist::{pb_key, PbKey};
use crate::pave_ff_update::{DmExEdges, PbFacesMap};
use crate::pave_filler::PaveFiller;

fn replace_vert(list: &mut Vec<usize>, old: usize, new: usize) {
    let mut hit = false;
    list.retain(|&v| {
        if v == old {
            hit = true;
            false
        } else {
            true
        }
    });
    if hit && !list.contains(&new) {
        list.push(new);
    }
}

fn real_pb(ds: &crate::bopds::BopdsDS, pb: &BopdsPaveBlock) -> BopdsPaveBlock {
    ds.real_pave_block(pb)
}

/// Substitute SD vertices inside FaceInfo On/In maps (`_6.cxx:1880-1903`).
pub fn update_face_info_vertices(f: &mut PaveFiller, dm_v: &HashMap<usize, usize>) {
    if dm_v.is_empty() {
        return;
    }
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let is_face = f
            .ds()
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if !is_face || f.ds().face_info(i).is_none() {
            continue;
        }
        if let Some(fi) = f.ds_mut().face_info_mut(i) {
            for (&n_v1, &n_v2) in dm_v {
                replace_vert(&mut fi.verts_on, n_v1, n_v2);
                replace_vert(&mut fi.verts_in, n_v1, n_v2);
                replace_vert(&mut fi.verts_sc, n_v1, n_v2);
            }
        }
    }
}

/// Substitute replacement pave blocks inside FaceInfo On/In/Sc (`_6.cxx:1905-1943`).
pub fn update_face_info_pave_blocks(
    f: &mut PaveFiller,
    dm_e: &DmExEdges,
) {
    if dm_e.is_empty() {
        return;
    }
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let is_face = f
            .ds()
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if !is_face || f.ds().face_info(i).is_none() {
            continue;
        }
        let (on, inn, sc) = match f.ds().face_info(i) {
            Some(fi) => (
                fi.paves_on().to_vec(),
                fi.paves_in().to_vec(),
                fi.paves().to_vec(),
            ),
            None => continue,
        };
        let mut fence: HashSet<PbKey> = HashSet::new();
        let mut new_on = Vec::new();
        let mut new_in = Vec::new();
        let mut new_sc = Vec::new();
        let lists = [&on, &inn, &sc];
        let outs = [&mut new_on, &mut new_in, &mut new_sc];
        // Rewrite one list at a time without holding FaceInfo.
        drop(outs);
        new_on = rewrite_pave_list(f, &on, dm_e, &mut fence);
        fence.clear();
        new_in = rewrite_pave_list(f, &inn, dm_e, &mut fence);
        fence.clear();
        new_sc = rewrite_pave_list(f, &sc, dm_e, &mut fence);
        if let Some(fi) = f.ds_mut().face_info_mut(i) {
            fi.paves_on = new_on;
            fi.paves_in = new_in;
            fi.paves = new_sc;
        }
    }
}

fn rewrite_pave_list(
    f: &PaveFiller,
    src: &[(usize, f64, f64)],
    dm_e: &DmExEdges,
    fence: &mut HashSet<PbKey>,
) -> Vec<(usize, f64, f64)> {
    let mut out = Vec::new();
    for &(e, first, last) in src {
        let mut pb = BopdsPaveBlock::new();
        pb.set_edge(e);
        pb.set_original_edge(e);
        pb.set_range(first, last);
        let key = pb_key(&pb);
        if let Some(repl) = dm_e.get(&key) {
            for r in repl {
                let real = real_pb(f.ds(), r);
                let rk = pb_key(&real);
                if fence.insert(rk) {
                    out.push((real.edge(), real.first, real.last));
                }
            }
        } else {
            let real = real_pb(f.ds(), &pb);
            let rk = pb_key(&real);
            if fence.insert(rk) {
                out.push((real.edge(), real.first, real.last));
            }
        }
    }
    out
}

/// Apply vertex + pave substitutions, then attach extra IN faces from
/// `thePBFacesMap` (existing `update_face_info_ff` already does Sc + IN
/// attachment; this is the On/In rewrite half).
pub fn update_face_info_on_in(
    f: &mut PaveFiller,
    dm_e: &DmExEdges,
    dm_v: &HashMap<usize, usize>,
    _pb_faces: &PbFacesMap,
) {
    update_face_info_vertices(f, dm_v);
    update_face_info_pave_blocks(f, dm_e);
}
