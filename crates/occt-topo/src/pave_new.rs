//! `PerformNewVertices` and `TreatNewVertices`.
//!
//! Source: `BOPAlgo_PaveFiller_3.cxx` (`PerformNewVertices` at 594,
//! `TreatNewVertices` at 692).

use crate::algo_tools_range::make_vertex_from_list;
use crate::bopalgo_tools::{intersect_vertices, VertexTol};
use crate::bopds::BopdsPaveBlock;
use crate::pave_common::update_vertex_sd;
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::add_ext_pave;
use crate::shape::TopoShape;

/// Couple of pave blocks that produced a new vertex (`BOPDS_CoupleOfPaveBlocks`).
#[derive(Clone)]
pub struct NewVertexCpb {
    pub shape: TopoShape,
    pub tol: f64,
    pub pb1: BopdsPaveBlock,
    pub pb2: BopdsPaveBlock,
    pub n_e1: usize,
    pub n_e2: usize,
    pub t1: f64,
    pub t2: f64,
    pub range1_first: f64,
    pub range1_last: f64,
    pub is_ee: bool,
    /// `BOPDS_CoupleOfPaveBlocks::IndexInterf` — typed `InterfEE`/`InterfEF`
    /// slot that `PerformNewVertices` binds `IndexNew` onto.
    pub index_interf: Option<usize>,
}

/// `BOPAlgo_PaveFiller::TreatNewVertices`.
pub fn treat_new_vertices_cpb(
    the_mvcpb: &[NewVertexCpb],
    fuzzy: f64,
) -> Result<Vec<(TopoShape, Vec<usize>)>, String> {
    if the_mvcpb.is_empty() {
        return Ok(Vec::new());
    }
    let verts: Vec<VertexTol> = the_mvcpb
        .iter()
        .map(|c| VertexTol {
            shape: c.shape.clone(),
            extra_tol: c.tol,
        })
        .collect();
    let chains = intersect_vertices(&verts, fuzzy);
    let mut images: Vec<(TopoShape, Vec<usize>)> = Vec::new();
    for chain in chains {
        if chain.is_empty() {
            continue;
        }
        let fused = make_vertex_from_list(&chain)?;
        let mut members: Vec<usize> = Vec::new();
        for s in &chain {
            if let Some(i) = the_mvcpb.iter().position(|c| c.shape.same_tshape(s)) {
                members.push(i);
            }
        }
        images.push((fused, members));
    }
    Ok(images)
}

/// `BOPAlgo_PaveFiller::PerformNewVertices`.
///
/// `b_is_ee_intersection` selects EE vs EF `SetIndexNew` (`_3.cxx:650`).
pub fn perform_new_vertices(
    f: &mut PaveFiller,
    the_mvcpb: &[NewVertexCpb],
    b_is_ee_intersection: bool,
) -> Result<(), String> {
    if the_mvcpb.is_empty() {
        return Ok(());
    }
    let fuzzy = f.fuzzy_value();
    let images = treat_new_vertices_cpb(the_mvcpb, fuzzy)?;
    let mut created: Vec<usize> = Vec::new();
    let mut member_to_nv: Vec<Option<usize>> = vec![None; the_mvcpb.len()];
    for (fused, members) in images {
        let n_v = f.ds_mut().append(fused)?;
        created.push(n_v);
        for i in members {
            member_to_nv[i] = Some(n_v);
        }
    }
    let mut modified: Vec<usize> = Vec::new();
    for (i, cpb) in the_mvcpb.iter().enumerate() {
        let Some(n_v) = member_to_nv[i] else {
            continue;
        };
        if cpb.is_ee {
            let _ = f.ds_mut().bind_ee_new_vertex(cpb.n_e1, cpb.n_e2, n_v);
            let _ = f
                .ds_mut()
                .set_ee_common_range(cpb.n_e1, cpb.n_e2, cpb.range1_first, cpb.range1_last);
        } else if let Some(ix) = cpb.index_interf {
            // `_3.cxx:648-651`: `aEFs(iX).SetIndexNew(iV)` on the VERTEX
            // record only; a later EDGE record of the same pair stays without
            // `IndexNew` so `UpdateFaceInfoIn` can still add `PaveBlocksIn`.
            let _ = f.ds_mut().bind_ef_new_vertex_at(ix, n_v);
            let _ = f
                .ds_mut()
                .set_ef_common_range_at(ix, cpb.range1_first, cpb.range1_last);
        } else {
            let _ = f.ds_mut().bind_ef_new_vertex(cpb.n_e1, cpb.n_e2, n_v);
            let _ = f
                .ds_mut()
                .set_ef_common_range(cpb.n_e1, cpb.n_e2, cpb.range1_first, cpb.range1_last);
        }
        add_ext_pave(f.ds_mut(), cpb.n_e1, cpb.t1, n_v);
        modified.push(cpb.n_e1);
        // `_3.cxx:671-683`: extra paves go on `PaveBlock1`/`PaveBlock2`. For
        // EF both handles are the same pave block (`aPB, aPB`) so the loop
        // breaks after the first; `n_e2` is the face and is not a pave list.
        if cpb.is_ee && cpb.n_e2 != cpb.n_e1 {
            add_ext_pave(f.ds_mut(), cpb.n_e2, cpb.t2, n_v);
            modified.push(cpb.n_e2);
        }
        let _ = update_vertex_sd(f, n_v, cpb.tol);
        let _ = b_is_ee_intersection;
        let _ = (cpb.pb1.original_edge(), cpb.pb2.original_edge());
    }
    modified.sort_unstable();
    modified.dedup();
    if !modified.is_empty() {
        crate::pave_split_blocks::split_pave_blocks(f, &modified, false)?;
    }
    let _ = created;
    Ok(())
}
