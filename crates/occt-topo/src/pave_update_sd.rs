//! Remainder of `BOPAlgo_PaveFiller_10.cxx`.
//!
//! `UpdateVertex`, `UpdateEdgeTolerance`, `UpdatePaveBlocksWithSDVertices`,
//! `UpdateCommonBlocksWithSDVertices`, `UpdateInterfsWithSDVertices`.
//! `SetNonDestructive` auto-detection from `TopoDS_Shape::Locked` has no
//! equivalent on the Rust `TopoShape` and is left as a translation boundary
//! (the explicit `set_non_destructive` setter is used instead).

use std::collections::HashSet;

use occt_core::precision::CONFUSION;

use crate::algo_tools::AlgoTools;
use crate::bopds::BopdsPaveBlock;
use crate::brep_tool::BRepTool;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Vertex};

/// `BOPAlgo_PaveFiller::UpdateVertex` (`_10.cxx:105`).
///
/// New vertices, vertices that already have an SD representative, and the
/// destructive mode update the (representative) vertex in place and grow its
/// DS box. In non-destructive mode an original argument vertex is replaced by
/// a new vertex linked through the SD map, recorded in
/// `myVertsToAvoidExtension`. Vertices whose tolerance grew are recorded in
/// `myIncreasedSS`.
pub fn update_vertex(f: &mut PaveFiller, n_v: usize, a_tol_new: f64) -> Result<usize, String> {
    let n_v_new = f.ds().has_shape_sd(n_v).unwrap_or(n_v);
    let has_sd = f.ds().has_shape_sd(n_v).is_some();
    let is_new = f.ds().is_new_shape(n_v_new) || f.ds().is_new_shape(n_v);

    if is_new || has_sd || !f.non_destructive() {
        let shape = f
            .ds()
            .shape(n_v_new)
            .cloned()
            .ok_or("update_vertex: bad vertex index")?;
        let vtx = Vertex(shape);
        let a_tol_v = BRepTool::vertex_tolerance(&vtx);
        if a_tol_v < a_tol_new {
            vtx.set_tolerance(a_tol_new);
            f.ds_mut().refresh_vertex_box(n_v_new, a_tol_new + CONFUSION);
            f.ds_mut().increased_ss_mut().insert(n_v);
        }
        return Ok(n_v_new);
    }

    let (point, a_tol_v) = {
        let shape = f
            .ds()
            .shape(n_v)
            .cloned()
            .ok_or("update_vertex: bad vertex index")?;
        let vtx = Vertex(shape);
        (BRepTool::vertex_point(&vtx), BRepTool::vertex_tolerance(&vtx))
    };
    let v_new = AlgoTools::make_new_vertex(&point, a_tol_v.max(a_tol_new))?;
    let n_new = f.ds_mut().append_info(crate::bopds::BopdsShapeInfo::new(v_new));
    f.ds_mut().refresh_vertex_box(n_new, a_tol_v.max(a_tol_new) + CONFUSION);
    f.ds_mut().add_shape_sd(n_v, n_new);
    f.verts_to_avoid_extension_mut().insert(n_new);
    if a_tol_v < a_tol_new {
        f.ds_mut().increased_ss_mut().insert(n_v);
    }
    Ok(n_new)
}

/// `BOPAlgo_PaveFiller::UpdateEdgeTolerance` (`_10.cxx:63`).
///
/// In non-destructive mode an original argument edge is left untouched, as is
/// an edge whose bound vertices are still original (no SD replacement).
/// Otherwise the edge tolerance is written, the DS box is rebuilt, and each
/// bound vertex is passed to [`update_vertex`].
pub fn update_edge_tolerance(f: &mut PaveFiller, n_e: usize, the_tol: f64) -> Result<(), String> {
    let subs = f
        .ds()
        .shape_info(n_e)
        .map(|s| s.sub_shapes().to_vec())
        .ok_or("update_edge_tolerance: bad edge index")?;
    if f.non_destructive() && !f.ds().is_new_shape(n_e) {
        return Ok(());
    }
    if f.non_destructive() {
        for &n_v in &subs {
            if !f.ds().is_new_shape(n_v) && f.ds().has_shape_sd(n_v).is_none() {
                return Ok(());
            }
        }
    }
    {
        let shape = f
            .ds()
            .shape(n_e)
            .cloned()
            .ok_or("update_edge_tolerance: bad edge index")?;
        let edge = Edge(shape);
        if let Some(mut g) = crate::tgeometry::GeometryRegistry::global().edge_geom(&edge.0) {
            g.tolerance = the_tol;
            crate::tgeometry::GeometryRegistry::global().set_edge(&edge.0, g);
        }
    }
    f.ds_mut().refresh_shape_box(n_e, CONFUSION);
    for &n_v in &subs {
        update_vertex(f, n_v, the_tol)?;
    }
    Ok(())
}

/// `BOPDS_DS::UpdatePaveBlockWithSDVertices`.
pub fn update_pave_block_with_sd_vertices(f: &PaveFiller, pb: &mut BopdsPaveBlock) {
    let (n1, n2) = pb.indices();
    pb.set_indices(
        f.ds().get_same_domain_index(n1),
        f.ds().get_same_domain_index(n2),
    );
}

/// `BOPAlgo_PaveFiller::UpdatePaveBlocksWithSDVertices` (`_10.cxx:166`).
pub fn update_pave_blocks_with_sd_vertices(f: &mut PaveFiller) {
    f.ds_mut().update_pave_blocks_with_sd_vertices();
}

/// `BOPDS_DS::UpdateCommonBlockWithSDVertices`.
///
/// Redirects every member pave block of the common block to its SD vertices
/// and writes the updated members back onto the stored common block.
pub fn update_common_block_with_sd_vertices(f: &mut PaveFiller, cb_idx: usize) {
    let Some(cb) = f.ds().common_blocks().get(cb_idx).cloned() else {
        return;
    };
    let mut pbs = cb.pave_blocks().to_vec();
    for pb in &mut pbs {
        update_pave_block_with_sd_vertices(f, pb);
    }
    if let Some(stored) = f.ds_mut().common_block_mut(cb_idx) {
        stored.set_pave_blocks(pbs);
    }
}

/// `BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices` (`_10.cxx:173`).
///
/// In destructive mode this reduces to [`update_pave_blocks_with_sd_vertices`].
/// In non-destructive mode each unique common block first has its bound
/// vertices passed through [`update_vertex`] with `Precision::Confusion()`,
/// then the common-block members themselves are redirected to SD vertices.
pub fn update_common_blocks_with_sd_vertices(f: &mut PaveFiller) -> Result<(), String> {
    if !f.non_destructive() {
        update_pave_blocks_with_sd_vertices(f);
        return Ok(());
    }
    let a_nb_pbp = f.ds().pave_blocks_pool().len();
    if a_nb_pbp == 0 {
        return Ok(());
    }
    let a_tol_v = CONFUSION;
    let mut seen: HashSet<(Vec<usize>, Vec<usize>)> = HashSet::new();
    let mut jobs: Vec<(usize, usize)> = Vec::new();
    for list in f.ds().pave_blocks_pool() {
        for pb in list {
            let Some(cb) = f.ds().common_block(pb) else {
                continue;
            };
            let mut edges = cb.indices().to_vec();
            edges.sort_unstable();
            let mut faces = cb.faces().to_vec();
            faces.sort_unstable();
            if !seen.insert((edges, faces)) {
                continue;
            }
            let (n_v1, n_v2) = pb.indices();
            jobs.push((n_v1, n_v2));
        }
    }
    for (n_v1, n_v2) in jobs {
        update_vertex(f, n_v1, a_tol_v)?;
        update_vertex(f, n_v2, a_tol_v)?;
    }
    let n_cb = f.ds().common_blocks().len();
    for i in 0..n_cb {
        update_common_block_with_sd_vertices(f, i);
    }
    update_pave_blocks_with_sd_vertices(f);
    Ok(())
}

/// `BOPAlgo_PaveFiller::UpdateInterfsWithSDVertices` (`_10.cxx:248`).
pub fn update_interfs_with_sd_vertices(f: &mut PaveFiller) {
    f.ds_mut().update_interfs_with_sd_vertices();
}

/// `BOPAlgo_PaveFiller::SetNonDestructive` (`_10.cxx:41`).
///
/// OCCT walks the arguments and turns the flag on when any argument is
/// `Locked`. The Rust `TopoShape` has no lock bit, so the method is a no-op
/// unless the flag is already set (matching `if (!myIsPrimary || myNonDestructive) return`).
pub fn set_non_destructive(f: &mut PaveFiller) {
    if !f.is_primary() || f.non_destructive() {
        return;
    }
}
