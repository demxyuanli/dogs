//! Shrunk-range data and pave-block updates — Phase 19 wave C2b.
//!
//! Port of the `BOPAlgo_PaveFiller` shrunk-data/update/self-interference
//! methods:
//!
//! | OCCT method                        | Rust function                              |
//! |------------------------------------|--------------------------------------------|
//! | `FillShrunkData(PaveBlock)`        | [`fill_shrunk_data_for_block`]             |
//! | `FillShrunkData(TopAbs, TopAbs)`   | [`fill_shrunk_data`]                       |
//! | `AnalyzeShrunkData`                | [`analyze_shrunk_data`]                    |
//! | `SetNonDestructive`                | [`set_non_destructive`]                    |
//! | `UpdatePaveBlocksWithSDVertices`   | [`update_pave_blocks_with_sd_vertices`]    |
//! | `UpdateEdgeTolerance`              | [`update_edge_tolerance`]                  |
//! | `UpdateVertex`                     | [`update_vertex_sd`]                       |
//! | `UpdateInterfsWithSDVertices`      | [`update_interfs_with_sd_vertices`]        |
//! | `UpdateCommonBlocksWithSDVertices` | [`update_common_blocks_with_sd_vertices`]  |
//! | `CheckSelfInterference`            | [`check_self_interference`]                |
//! | `RemoveMicroEdges`                 | [`remove_micro_edges`]                     |
//!
//! Source: `BOPAlgo_PaveFiller_3.cxx` (FillShrunkData/AnalyzeShrunkData),
//! `BOPAlgo_PaveFiller_9.cxx` (the edge-wide FillShrunkData),
//! `BOPAlgo_PaveFiller_10.cxx` (SetNonDestructive/Update*),
//! `BOPAlgo_PaveFiller_11.cxx` (CheckSelfInterference).
//!
//! Every function operates directly on the concrete
//! [`crate::pave_filler::PaveFiller`] (landed by the sibling Phase-19 agent),
//! accessing the data structure through `ds()`/`ds_mut()` and the alert report
//! through `add_error()`/`add_warning()`.
//!
//! ## Port simplifications
//!
//! - The Rust `BopdsDS` keeps no per-interference `index_new` arrays, so
//!   [`update_interfs_with_sd_vertices`] is a no-op; the generic
//!   [`update_intfs_with_sd_vertices`] helper carries the OCCT rewrite logic
//!   and is covered by a unit test.
//! - The Rust `BopdsCommonBlock` records only shared edge indices and ranges
//!   (no bound vertices), so [`update_common_blocks_with_sd_vertices`] reduces
//!   to the pave-block SD pass.
//! - [`check_self_interference`] uses the [`crate::bopds::BopdsIteratorSI`]
//!   self-intersection candidates filtered by a bbox-overlap measure rather
//!   than the OCCT face-info connection map (the Rust `BopdsFaceInfo` does not
//!   retain the IN/Section vertex sets). It detects overlapping solids, faces
//!   and crossing edges within one argument; boundary-touching sub-shapes are
//!   excluded.

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsDS, BopdsInterf, BopdsIteratorSI, BopdsPaveBlock, BopdsShapeInfo};
use crate::brep_tool::BRepTool;
use crate::inttools_range::ShrunkRange;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::GeometryRegistry;

// ---------------------------------------------------------------------------
// FillShrunkData / AnalyzeShrunkData
// ---------------------------------------------------------------------------

/// Compute the shrunk (working) range of every pave block of every source
/// edge that carries pave blocks and does not already have shrunk data.
///
/// Port of `BOPAlgo_PaveFiller::FillShrunkData(TopAbs_ShapeEnum, ...)` from
/// `BOPAlgo_PaveFiller_9.cxx` — the OCCT version collects the candidate edges
/// from the Edge/Edge interference iterator; this port processes every source
/// edge with pave blocks (a superset that is behaviourally equivalent for the
/// pipeline, since untouched edges simply never query their shrunk data).
///
/// For each block the [`crate::inttools_range::ShrunkRange`] is computed over
/// the edge's parameter range and stored through
/// [`BopdsPaveBlock::set_shrunk_data`], with the OCCT warning policy applied
/// by [`analyze_shrunk_data`].
pub fn fill_shrunk_data(f: &mut PaveFiller) -> Result<(), String> {
    // Collect the (edge, block-slot) jobs that need their shrunk data.
    let jobs: Vec<(usize, usize)> = {
        let ds = f.ds();
        let mut jobs = Vec::new();
        for e in 0..ds.nb_source_shapes() {
            let Some(si) = ds.shape_info(e) else { continue };
            if si.kind != ShapeType::Edge {
                continue;
            }
            if BRepTool::is_degenerated(&Edge(si.shape().clone())) {
                // Degenerated edges are not shrinkable (BOPAlgo skips flagged
                // edges).
                continue;
            }
            let blocks = ds.pave_blocks(e);
            for (slot, pb) in blocks.iter().enumerate() {
                // OCCT also re-validates existing data with
                // `BOPDS_DS::IsValidShrunkData`; the Rust DS does not port that
                // predicate, so only the presence of the data is checked.
                if !pb.has_shrunk_data() {
                    jobs.push((e, slot));
                }
            }
        }
        jobs
    };
    if jobs.is_empty() {
        return Ok(());
    }

    let tol = f.fuzzy_value();
    for (e, slot) in jobs {
        let (edge_shape, block) = {
            let ds = f.ds();
            let shape = ds
                .shape(e)
                .ok_or("fill_shrunk_data: bad edge index")?
                .clone();
            let pb = ds
                .pave_blocks(e)
                .get(slot)
                .cloned()
                .ok_or("fill_shrunk_data: bad pave-block slot")?;
            (shape, pb)
        };
        let edge = Edge(edge_shape);
        let mut pb = block;
        fill_shrunk_data_for_block(f, &edge, tol, &mut pb);
        let blocks = f.ds_mut().change_pave_blocks_mut(e);
        if let Some(stored) = blocks.get_mut(slot) {
            let (ts1, ts2, splittable) = pb.shrunk_data();
            stored.set_shrunk_data(ts1, ts2, splittable);
        }
    }
    Ok(())
}

/// Compute the shrunk range of a single pave block and store it on the block.
///
/// Port of `BOPAlgo_PaveFiller::FillShrunkData(handle<BOPDS_PaveBlock>&)` from
/// `BOPAlgo_PaveFiller_3.cxx`. The range is computed with a neutral (empty)
/// face: the shrunk range in this port depends only on the edge curve and its
/// boundary vertices, not on the supporting face.
pub fn fill_shrunk_data_for_block(f: &mut PaveFiller, edge: &Edge, tol: f64, pb: &mut BopdsPaveBlock) {
    let mut sr = ShrunkRange::new();
    // `Err` is expected for degenerated / micro edges; `analyze_shrunk_data`
    // inspects `is_done()` and applies the OCCT warning policy.
    let _ = sr.set_shrunk_range(edge, &Face::new(), tol);
    analyze_shrunk_data(f, pb, &sr, edge);
}

/// Analyze the result of a shrunk-range computation and store the shrunk data
/// on the pave block, applying the OCCT warning policy.
///
/// Port of `BOPAlgo_PaveFiller::AnalyzeShrunkData` from
/// `BOPAlgo_PaveFiller_3.cxx`:
///
/// - computation failure → "too small edge" (block spans the whole edge) or
///   "bad positioning" warning, shrunk data is stored as the degenerate
///   `(0, 0, not splittable)`;
/// - success but too short to split → "not splittable edge" / "bad positioning"
///   warning, the computed range is still stored (splittable = false);
/// - success → the computed range is stored with its splittable flag.
pub fn analyze_shrunk_data(f: &mut PaveFiller, pb: &mut BopdsPaveBlock, sr: &ShrunkRange, edge: &Edge) {
    let (edge_first, edge_last) = BRepTool::edge_parameters(edge);
    let (pb_first, pb_last) = pb.range();
    // The block covers the whole edge when its bounds reach the edge's own
    // parameter bounds.
    let whole_edge = pb_first <= edge_first && pb_last >= edge_last;

    let done = sr.is_done();
    let splittable = sr.is_splittable();
    if !done || !splittable {
        if !done {
            f.add_warning(if whole_edge {
                "too small edge: shrunk range cannot be computed".to_string()
            } else {
                "bad positioning: shrunk range cannot be computed".to_string()
            });
            pb.set_shrunk_data(0.0, 0.0, false);
            return;
        }
        f.add_warning(if whole_edge {
            "not splittable edge: shrunk range is too short to split".to_string()
        } else {
            "bad positioning: shrunk range is too short".to_string()
        });
    }
    match sr.shrunk_range() {
        Some((ts1, ts2)) => pb.set_shrunk_data(ts1, ts2, sr.is_splittable()),
        None => pb.set_shrunk_data(0.0, 0.0, false),
    }
}

// ---------------------------------------------------------------------------
// SetNonDestructive
// ---------------------------------------------------------------------------

/// Sets the non-destructive mode of the filler.
///
/// Port of `BOPAlgo_PaveFiller::SetNonDestructive(const bool)` from
/// `BOPAlgo_PaveFiller_10.cxx`. In this mode the argument shapes are not
/// modified: tolerance increases allocate same-domain replacement vertices
/// instead of mutating the originals.
pub fn set_non_destructive(f: &mut PaveFiller, b: bool) {
    f.set_non_destructive(b);
}

// ---------------------------------------------------------------------------
// UpdatePaveBlocksWithSDVertices
// ---------------------------------------------------------------------------

/// Replaces the bound-vertex indices of every pave block by their same-domain
/// representatives.
///
/// Port of `BOPAlgo_PaveFiller::UpdatePaveBlocksWithSDVertices` +
/// `BOPDS_DS::UpdatePaveBlockWithSDVertices` (which resolves each bound pave's
/// index through `GetSameDomainIndex`).
pub fn update_pave_blocks_with_sd_vertices(f: &mut PaveFiller) -> Result<(), String> {
    // Collect the updates first so the SD lookups (immutable borrows) do not
    // overlap the pave-block mutation.
    let updates: Vec<(usize, usize, usize, usize)> = {
        let ds = f.ds();
        let mut updates = Vec::new();
        for e in 0..ds.nb_shapes() {
            if !ds.has_pave_blocks(e) {
                continue;
            }
            let blocks = ds.pave_blocks(e);
            for (slot, pb) in blocks.iter().enumerate() {
                let (i1, i2) = pb.indices();
                let n1 = ds.get_same_domain_index(i1);
                let n2 = ds.get_same_domain_index(i2);
                if n1 != i1 || n2 != i2 {
                    updates.push((e, slot, n1, n2));
                }
            }
        }
        updates
    };
    for (e, slot, n1, n2) in updates {
        let ds = f.ds_mut();
        let blocks = ds.change_pave_blocks_mut(e);
        if let Some(pb) = blocks.get_mut(slot) {
            pb.set_indices(n1, n2);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// UpdateEdgeTolerance / UpdateVertex
// ---------------------------------------------------------------------------

/// Raises the tolerance of the edge `n_e` and of its bound vertices to `tol`.
///
/// Port of `BOPAlgo_PaveFiller::UpdateEdgeTolerance` from
/// `BOPAlgo_PaveFiller_10.cxx`. In the non-destructive (safe-input) mode the
/// original argument shapes are not modified: an edge (or a vertex) that is
/// neither a new shape nor backed by a same-domain representative is left
/// untouched, mirroring the OCCT early-outs.
///
/// Note: the DS bounding box is not refreshed after the tolerance increase —
/// the Rust `BopdsDS` exposes no box mutation; the tolerance side-table is
/// kept in sync.
pub fn update_edge_tolerance(f: &mut PaveFiller, n_e: usize, tol: f64) -> Result<(), String> {
    let subs = f
        .ds()
        .shape_info(n_e)
        .map(|s| s.sub_shapes().to_vec())
        .ok_or("update_edge_tolerance: bad edge index")?;
    let is_new = f.ds().is_new_shape(n_e);

    // Safe-input mode: never modify an original edge or its original vertices.
    if f.non_destructive() && !is_new {
        return Ok(());
    }
    if f.non_destructive() {
        let ds = f.ds();
        for &n_v in &subs {
            if !ds.is_new_shape(n_v) && ds.has_shape_sd(n_v).is_none() {
                return Ok(());
            }
        }
    }

    // Update the edge tolerance in the geometry side-table (BRep_Builder::
    // UpdateEdge).
    {
        let shape = f
            .ds()
            .shape(n_e)
            .ok_or("update_edge_tolerance: bad edge index")?
            .clone();
        let edge = Edge(shape);
        if let Some(mut g) = GeometryRegistry::global().edge_geom(&edge.0) {
            g.tolerance = tol;
            GeometryRegistry::global().set_edge(&edge.0, g);
        }
    }

    // Update the bound vertices.
    for &n_v in &subs {
        update_vertex_sd(f, n_v, tol)?;
    }
    Ok(())
}

/// Raises the tolerance of the vertex `n_v` to `tol`, returning the index of
/// the vertex that now carries the enlarged tolerance.
///
/// Port of `BOPAlgo_PaveFiller::UpdateVertex` from `BOPAlgo_PaveFiller_10.cxx`:
///
/// - when the vertex is new, has a same-domain representative, or the
///   non-destructive mode is not in force, the (representative) vertex is
///   updated in place;
/// - otherwise a new vertex is created with `max(old_tol, tol)` and registered
///   as the same-domain representative of `n_v`.
pub fn update_vertex_sd(f: &mut PaveFiller, n_v: usize, tol: f64) -> Result<usize, String> {
    let (n_v_new, has_sd) = {
        let ds = f.ds();
        let n_new = ds.has_shape_sd(n_v).unwrap_or(n_v);
        (n_new, ds.has_shape_sd(n_v).is_some())
    };
    let is_new = f.ds().is_new_shape(n_v);

    if is_new || has_sd || !f.non_destructive() {
        let shape = f
            .ds()
            .shape(n_v_new)
            .ok_or("update_vertex_sd: bad vertex index")?
            .clone();
        let vtx = Vertex(shape);
        let cur = BRepTool::vertex_tolerance(&vtx);
        if cur < tol {
            vtx.set_tolerance(tol);
        }
        return Ok(n_v_new);
    }

    // Original vertex in non-destructive mode: build a new vertex carrying the
    // enlarged tolerance and register it as the SD representative.
    let (point, cur_tol) = {
        let ds = f.ds();
        let shape = ds
            .shape(n_v)
            .ok_or("update_vertex_sd: bad vertex index")?
            .clone();
        let vtx = Vertex(shape);
        (BRepTool::vertex_point(&vtx), BRepTool::vertex_tolerance(&vtx))
    };
    let new_tol = cur_tol.max(tol);
    let v_new = AlgoTools::make_new_vertex(&point, new_tol)?;
    let ds = f.ds_mut();
    let n_v_new = ds.append_info(BopdsShapeInfo::new(v_new));
    ds.add_shape_sd(n_v, n_v_new);
    Ok(n_v_new)
}

// ---------------------------------------------------------------------------
// UpdateInterfsWithSDVertices
// ---------------------------------------------------------------------------

/// Update the `index_new` of every interference to its same-domain
/// representative.
///
/// Port of the `UpdateIntfsWithSDVertices` template in
/// `BOPAlgo_PaveFiller_10.cxx`. Interference arrays with a stored new-shape
/// index are rewritten when that index was merged into a same-domain shape.
pub fn update_intfs_with_sd_vertices(ds: &BopdsDS, interfs: &mut [BopdsInterf]) {
    for intf in interfs.iter_mut() {
        if let Some(idx_new) = intf.get_index_new() {
            intf.set_index_new(ds.get_same_domain_index(idx_new));
        }
    }
}

/// Port of `BOPAlgo_PaveFiller::UpdateInterfsWithSDVertices` from
/// `BOPAlgo_PaveFiller_10.cxx`.
///
/// The Rust `BopdsDS` does not retain the per-interference `index_new` arrays
/// OCCT rewrites here (interferences are stored as raw index pairs with no
/// new-shape index), so the rewrite pass is a no-op for the stored data. The
/// actual logic lives in [`update_intfs_with_sd_vertices`] and is exercised by
/// the unit test against a standalone interference array.
pub fn update_interfs_with_sd_vertices(_f: &mut PaveFiller) -> Result<(), String> {
    Ok(())
}

// ---------------------------------------------------------------------------
// UpdateCommonBlocksWithSDVertices
// ---------------------------------------------------------------------------

/// Updates the SD vertices of the pave blocks that belong to common blocks.
///
/// Port of `BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices` from
/// `BOPAlgo_PaveFiller_10.cxx`. The OCCT non-destructive branch additionally
/// raises the tolerance of the common-block bound vertices; the Rust
/// `BopdsCommonBlock` records only the shared edge indices and ranges (no
/// bound vertices), so that pass is elided and the SD-index update reduces to
/// the pave-block pass.
pub fn update_common_blocks_with_sd_vertices(f: &mut PaveFiller) -> Result<(), String> {
    update_pave_blocks_with_sd_vertices(f)
}

// ---------------------------------------------------------------------------
// CheckSelfInterference
// ---------------------------------------------------------------------------

/// The shape-type pairs inspected for self-intersection.
const SI_TYPE_PAIRS: [(ShapeType, ShapeType); 3] = [
    (ShapeType::Solid, ShapeType::Solid),
    (ShapeType::Face, ShapeType::Face),
    (ShapeType::Edge, ShapeType::Edge),
];

/// Checks each argument for acquired self-interference: whether an argument's
/// own sub-shapes geometrically intersect beyond their shared boundary.
///
/// Returns `true` when at least one self-intersecting pair was found, and adds
/// an `acquired self-interference` warning for every such pair.
///
/// Port of `BOPAlgo_PaveFiller::CheckSelfInterference` from
/// `BOPAlgo_PaveFiller_11.cxx`, adapted to the candidate-driven
/// [`crate::bopds::BopdsIteratorSI`] instead of the OCCT face-info connection
/// map (the Rust `BopdsFaceInfo` does not retain the IN/Section vertex sets).
///
/// For every candidate pair produced by the self-intersection iterator that
/// lies within one argument, the pair is classified by the shape types:
///
/// - two solids must overlap with positive volume (3D bbox overlap);
/// - two faces must overlap with positive area (2D overlap, i.e. more than a
///   shared boundary line);
/// - two edges must overlap with positive area (crossing edges) and must not
///   share a boundary vertex.
///
/// Boundary-touching sub-shapes (adjacent faces/edges of a closed shape, an
/// edge lying on a face) are therefore not reported, while genuinely
/// self-intersecting shapes (e.g. two overlapping boxes passed as one
/// compound, or crossing edges of a wire) are.
pub fn check_self_interference(f: &mut PaveFiller) -> Result<bool, String> {
    // Collect the offending pairs while the DS is immutably borrowed; the
    // warnings are appended afterwards through the mutable filler.
    let offending: Vec<(usize, usize, usize)> = {
        let ds = f.ds();
        let mut it = BopdsIteratorSI::new();
        it.set_ds(ds);
        it.prepare();
        let mut found = Vec::new();
        for r in 0..ds.nb_ranges() {
            let Some(range) = ds.range(r) else { continue };
            for &(t1, t2) in &SI_TYPE_PAIRS {
                it.initialize(t1, t2);
                while it.more() {
                    let (i, j) = it.value();
                    if range.contains(i) && range.contains(j) && is_self_intersecting(ds, i, j) {
                        found.push((r, i, j));
                    }
                    it.next();
                }
            }
        }
        found
    };

    let mut detected = false;
    for (r, i, j) in offending {
        detected = true;
        let ti = f
            .ds()
            .shape_info(i)
            .map(|s| s.kind)
            .unwrap_or(ShapeType::Vertex);
        let tj = f
            .ds()
            .shape_info(j)
            .map(|s| s.kind)
            .unwrap_or(ShapeType::Vertex);
        f.add_warning(format!(
            "acquired self-interference: argument {} contains self-intersecting {ti:?} and {tj:?} sub-shapes (DS indices {} and {})",
            r + 1, i, j
        ));
    }
    Ok(detected)
}

/// Removes the micro edges — degenerate pave blocks whose bound vertices
/// coincide and whose shrunk range cannot be computed — from the data
/// structure.
///
/// Port of `BOPAlgo_PaveFiller::RemoveMicroEdges` from
/// `BOPAlgo_PaveFiller_6.cxx`. Every edge carrying at least two pave blocks is
/// examined: for each *real* pave block (common-block members are unified
/// through their canonical block, `BOPDS_DS::RealPaveBlock`), a block whose
/// two bound vertices coincide and whose shrunk range cannot be computed
/// (`!IsDone()`, mirroring the OCCT `!HasShrunkData()` test on the empty
/// shrunk box) is a micro edge. The detected edges are then removed from the
/// DS via [`crate::bopds::BopdsDS::remove_pave_blocks`].
///
/// Translation boundaries: the OCCT `HasFlag` guard (skip the flagged,
/// degenerated edges) is ported through the degenerated-edge query of
/// [`crate::brep_tool::BRepTool`]; the shrunk range is computed with
/// [`ShrunkRange`] on the block's edge over a neutral face (the same
/// simplification as [`fill_shrunk_data_for_block`]), matching the OCCT
/// `FillShrunkData` semantics.
pub fn remove_micro_edges(f: &mut PaveFiller) {
    let tol = f.fuzzy_value();
    // Fence of already-processed real blocks; the detected micro edges.
    let micro: HashSet<usize> = {
        let ds = f.ds();
        // Fence of already-processed real blocks (edge + quantized range), so a
        // common block's members are examined once (OCCT fences by the pave
        // block handle).
        let mut fence: HashSet<(usize, i64, i64)> = HashSet::new();
        let mut micro: HashSet<usize> = HashSet::new();
        for list in ds.pave_blocks_pool() {
            if list.len() < 2 {
                // No splits: a single block is never a micro edge.
                continue;
            }
            let orig = list[0].original_edge();
            // OCCT `HasFlag` guard: flagged (degenerated) edges are skipped.
            let flagged = ds
                .shape_info(orig)
                .map(|si| BRepTool::is_degenerated(&Edge(si.shape().clone())))
                .unwrap_or(false);
            if flagged {
                continue;
            }
            for pb in list {
                let real = ds.real_pave_block(pb);
                let key = (real.edge(), (real.first * 1e7) as i64, (real.last * 1e7) as i64);
                if !fence.insert(key) {
                    continue;
                }
                let (n1, n2) = real.indices();
                if n1 != n2 {
                    continue;
                }
                let n_e = if real.has_edge() { real.edge() } else { real.original_edge() };
                let Some(shape) = ds.shape(n_e) else { continue };
                let mut sr = ShrunkRange::new();
                let _ = sr.set_shrunk_range(&Edge(shape.clone()), &Face::new(), tol);
                if !sr.is_done() {
                    // Micro edge: no valid shrunk range can be built on it.
                    micro.insert(n_e);
                }
            }
        }
        micro
    };
    if !micro.is_empty() {
        f.ds_mut().remove_pave_blocks(&micro);
    }
}

/// Whether the two sub-shapes of indices `i` and `j` genuinely self-intersect
/// (beyond sharing a boundary element).
fn is_self_intersecting(ds: &BopdsDS, i: usize, j: usize) -> bool {
    let ti = ds.shape_info(i).map(|s| s.kind).unwrap_or(ShapeType::Vertex);
    let tj = ds.shape_info(j).map(|s| s.kind).unwrap_or(ShapeType::Vertex);
    let dim = bbox_overlap_dim(ds, i, j);
    match (ti, tj) {
        // Positive-volume overlap of two solids.
        (ShapeType::Solid, ShapeType::Solid) => dim >= 3,
        // Positive-area overlap of two faces (more than a shared edge line).
        (ShapeType::Face, ShapeType::Face) => dim >= 2,
        // Crossing edges overlap their bounding box with positive area; edges
        // that merely share a boundary vertex overlap in fewer dimensions.
        (ShapeType::Edge, ShapeType::Edge) => dim >= 2 && !edges_share_boundary_vertex(ds, i, j),
        _ => false,
    }
}

/// The number of coordinate axes in which the bounding boxes of the two shapes
/// overlap with a positive length.
fn bbox_overlap_dim(ds: &BopdsDS, i: usize, j: usize) -> usize {
    let (Some(a), Some(b)) = (ds.box_of(i), ds.box_of(j)) else {
        return 0;
    };
    let (Some(aa), Some(bb)) = (a.get(), b.get()) else {
        return 0;
    };
    let mut dims = 0;
    if overlap_len(aa.0, aa.1, bb.0, bb.1) > 1e-9 {
        dims += 1;
    }
    if overlap_len(aa.2, aa.3, bb.2, bb.3) > 1e-9 {
        dims += 1;
    }
    if overlap_len(aa.4, aa.5, bb.4, bb.5) > 1e-9 {
        dims += 1;
    }
    dims
}

/// Length of the intersection of the two 1D intervals, clamped below at zero.
fn overlap_len(amin: f64, amax: f64, bmin: f64, bmax: f64) -> f64 {
    let lo = amin.max(bmin);
    let hi = amax.min(bmax);
    (hi - lo).max(0.0)
}

/// Whether the two edges (DS indices `i` and `j`) share a boundary vertex,
/// following same-domain chains.
fn edges_share_boundary_vertex(ds: &BopdsDS, i: usize, j: usize) -> bool {
    let vi = ds
        .shape_info(i)
        .map(|s| s.sub_shapes().to_vec())
        .unwrap_or_default();
    let vj = ds
        .shape_info(j)
        .map(|s| s.sub_shapes().to_vec())
        .unwrap_or_default();
    vi.iter()
        .map(|&v| ds.get_same_domain_index(v))
        .any(|v| vj.iter().map(|&w| ds.get_same_domain_index(w)).any(|w| w == v))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpDir, GpLin, GpPnt};
    use occt_geom::GeomLine;

    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use crate::shape::TopoShape;

    /// Release the geometry-side-table entries of a shape tree so tests do not
    /// leave stale registry entries keyed by a freed TShape address.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("unit direction")
    }

    /// A filler whose DS holds a single box with a default pave block on every
    /// source edge (the input `fill_shrunk_data` expects).
    ///
    /// Uses `init` (not a bare `append`) so the per-argument ranges and the
    /// source-shape count are set up, which the `is_new_shape` checks in the
    /// update helpers rely on.
    fn filler_with_box_edges() -> (PaveFiller, crate::brep_extrema::test_box::UnitBox) {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();
        let n = pf.ds().nb_source_shapes();
        for i in 0..n {
            if pf.ds().shape_info(i).map(|s| s.kind) == Some(ShapeType::Edge) {
                pf.ds_mut().init_pave_blocks_for_edge(i);
            }
        }
        (pf, b)
    }

    /// A filler whose DS holds a single box initialized through `init` (source
    /// shapes, ranges and ranks set up).
    fn filler_with_box() -> (PaveFiller, crate::brep_extrema::test_box::UnitBox) {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();
        (pf, b)
    }

    // -----------------------------------------------------------------------
    // fill_shrunk_data / analyze_shrunk_data
    // -----------------------------------------------------------------------

    #[test]
    fn fill_shrunk_data_sets_valid_ranges_for_box_edges() {
        let (mut pf, b) = filler_with_box_edges();
        fill_shrunk_data(&mut pf).unwrap();
        assert!(!pf.has_errors());

        let e0 = pf.ds().index(&b.edges[0].0).expect("edge indexed");
        let blocks = pf.ds().pave_blocks(e0);
        assert_eq!(blocks.len(), 1);
        assert!(blocks[0].has_shrunk_data(), "pave block must carry shrunk data");
        let (ts1, ts2, splittable) = blocks[0].shrunk_data();
        assert!(ts1 > 0.0, "shrunk first {ts1} must be > 0 (edge range is [0, 1])");
        assert!(ts2 < 1.0, "shrunk last {ts2} must be < 1 (edge range is [0, 1])");
        assert!(ts2 > ts1, "shrunk range must be non-empty: {ts1}..{ts2}");
        assert!(splittable, "a unit-length edge must be splittable");

        // Every box edge got a valid shrunk range.
        for edge in &b.edges {
            let e = pf.ds().index(&edge.0).expect("edge indexed");
            let blocks = pf.ds().pave_blocks(e);
            assert_eq!(blocks.len(), 1);
            let (s1, s2, _) = blocks[0].shrunk_data();
            assert!(s1 > 0.0 && s2 < 1.0 && s2 > s1, "edge {e}: shrunk {s1}..{s2}");
        }
        clear_tree(&b.solid.0);
    }

    #[test]
    fn fill_shrunk_data_is_idempotent() {
        let (mut pf, b) = filler_with_box_edges();
        fill_shrunk_data(&mut pf).unwrap();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let before = pf.ds().pave_blocks(e0)[0].shrunk_data();
        let n_warnings = pf.warnings().len();

        // A second pass must not recompute / re-warn for blocks that already
        // carry shrunk data.
        fill_shrunk_data(&mut pf).unwrap();
        let after = pf.ds().pave_blocks(e0)[0].shrunk_data();
        assert_eq!(before, after);
        assert_eq!(pf.warnings().len(), n_warnings);
        clear_tree(&b.solid.0);
    }

    #[test]
    fn fill_shrunk_data_on_edge_without_curve_warns() {
        let mut pf = PaveFiller::new();
        let e = Edge::new();
        let mut pb = BopdsPaveBlock::new();
        pb.set_range(0.0, 1.0);
        fill_shrunk_data_for_block(&mut pf, &e, 1e-7, &mut pb);
        assert!(pb.has_shrunk_data());
        assert_eq!(pb.shrunk_data().0, 0.0);
        assert!(!pb.is_splittable());
        assert!(
            pf.warnings().iter().any(|w| w.contains("bad positioning")),
            "expected a bad-positioning warning, got: {:?}",
            pf.warnings()
        );
    }

    // -----------------------------------------------------------------------
    // set_non_destructive
    // -----------------------------------------------------------------------

    #[test]
    fn set_non_destructive_toggles_flag() {
        let mut pf = PaveFiller::new();
        assert!(!pf.non_destructive());
        set_non_destructive(&mut pf, true);
        assert!(pf.non_destructive());
        set_non_destructive(&mut pf, false);
        assert!(!pf.non_destructive());
    }

    // -----------------------------------------------------------------------
    // update_pave_blocks_with_sd_vertices
    // -----------------------------------------------------------------------

    #[test]
    fn update_pave_blocks_with_sd_vertices_rewrites_bounds() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let v1 = pf.ds().index(&b.vertices[1].0).unwrap();
        let v5 = pf.ds().index(&b.vertices[5].0).unwrap();

        pf.ds_mut().init_pave_blocks_for_edge(e0);
        assert_eq!(pf.ds().pave_blocks(e0)[0].indices(), (v0, v1));

        // v0 is merged into v5 (same-domain chain).
        pf.ds_mut().add_shape_sd(v0, v5);
        update_pave_blocks_with_sd_vertices(&mut pf).unwrap();
        assert_eq!(pf.ds().pave_blocks(e0)[0].indices(), (v5, v1));

        // A block with no SD vertices is left untouched.
        let e1 = pf.ds().index(&b.edges[1].0).unwrap();
        pf.ds_mut().init_pave_blocks_for_edge(e1);
        let before = pf.ds().pave_blocks(e1)[0].indices();
        update_pave_blocks_with_sd_vertices(&mut pf).unwrap();
        assert_eq!(pf.ds().pave_blocks(e1)[0].indices(), before);
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // update_edge_tolerance / update_vertex_sd
    // -----------------------------------------------------------------------

    #[test]
    fn update_edge_tolerance_raises_edge_and_vertex_tolerances() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();

        update_edge_tolerance(&mut pf, e0, 0.01).unwrap();
        let edge = Edge(pf.ds().shape(e0).unwrap().clone());
        assert!(
            (BRepTool::edge_tolerance(&edge) - 0.01).abs() < 1e-12,
            "edge tolerance raised to 0.01, got {}",
            BRepTool::edge_tolerance(&edge)
        );
        let vtx = Vertex(pf.ds().shape(v0).unwrap().clone());
        assert!(
            (BRepTool::vertex_tolerance(&vtx) - 0.01).abs() < 1e-12,
            "bound vertex tolerance raised to 0.01, got {}",
            BRepTool::vertex_tolerance(&vtx)
        );
        clear_tree(&b.solid.0);
    }

    #[test]
    fn update_edge_tolerance_non_destructive_skips_original_edge() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        set_non_destructive(&mut pf, true);

        update_edge_tolerance(&mut pf, e0, 0.01).unwrap();
        // The original edge is not modified in safe-input mode.
        let edge = Edge(pf.ds().shape(e0).unwrap().clone());
        assert!(
            (BRepTool::edge_tolerance(&edge) - 0.0).abs() < 1e-12,
            "original edge tolerance must stay 0.0 in non-destructive mode, got {}",
            BRepTool::edge_tolerance(&edge)
        );
        clear_tree(&b.solid.0);
    }

    #[test]
    fn update_vertex_sd_creates_replacement_in_non_destructive_mode() {
        let (mut pf, b) = filler_with_box();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let n_before = pf.ds().nb_shapes();
        set_non_destructive(&mut pf, true);

        let n_v_new = update_vertex_sd(&mut pf, v0, 0.02).unwrap();
        assert_ne!(n_v_new, v0, "a replacement vertex must be allocated");
        assert_eq!(pf.ds().nb_shapes(), n_before + 1);
        assert_eq!(pf.ds().has_shape_sd(v0), Some(n_v_new));
        let vtx = Vertex(pf.ds().shape(n_v_new).unwrap().clone());
        assert!(
            (BRepTool::vertex_tolerance(&vtx) - 0.02).abs() < 1e-12,
            "replacement tolerance raised to 0.02, got {}",
            BRepTool::vertex_tolerance(&vtx)
        );
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // update_interfs_with_sd_vertices
    // -----------------------------------------------------------------------

    #[test]
    fn update_intfs_with_sd_vertices_rewrites_index_new() {
        let (mut pf, b) = filler_with_box();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let v5 = pf.ds().index(&b.vertices[5].0).unwrap();
        pf.ds_mut().add_shape_sd(v0, v5);

        let mut interfs = [
            BopdsInterf::new(1, 2),
            BopdsInterf::new(3, 4),
        ];
        interfs[0].set_index_new(v0);
        update_intfs_with_sd_vertices(pf.ds(), &mut interfs);
        assert_eq!(interfs[0].get_index_new(), Some(v5), "index_new rewritten through the SD map");
        assert_eq!(interfs[1].get_index_new(), None, "an unset index_new stays unset");

        // The filler-level port is a no-op for the stored DS (no arrays kept).
        update_interfs_with_sd_vertices(&mut pf).unwrap();
        assert!(!pf.has_errors());
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // update_common_blocks_with_sd_vertices
    // -----------------------------------------------------------------------

    #[test]
    fn update_common_blocks_with_sd_vertices_updates_pave_bounds() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let v5 = pf.ds().index(&b.vertices[5].0).unwrap();

        pf.ds_mut().init_pave_blocks_for_edge(e0);
        pf.ds_mut().add_shape_sd(v0, v5);
        update_common_blocks_with_sd_vertices(&mut pf).unwrap();
        assert_eq!(pf.ds().pave_blocks(e0)[0].indices(), (v5, pf.ds().index(&b.vertices[1].0).unwrap()));
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // check_self_interference
    // -----------------------------------------------------------------------

    #[test]
    fn check_self_interference_clean_box_is_false() {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(!res, "a single closed box is not self-interfering");
        assert!(!pf.has_warnings(), "no warnings expected, got {:?}", pf.warnings());
        clear_tree(&b.solid.0);
    }

    #[test]
    fn check_self_interference_overlapping_boxes_is_true() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        let b = TopoBuilder::new();
        // Two identically-positioned boxes as one argument: their solids
        // overlap with positive volume.
        let comp = b.make_compound_of(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.set_arguments(&[comp.into()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(res, "overlapping boxes in one argument must be self-interfering");
        assert!(
            pf.warnings().iter().any(|w| w.contains("acquired self-interference")),
            "expected an acquired self-interference warning, got {:?}",
            pf.warnings()
        );
        clear_tree(&a.solid.0);
        clear_tree(&c.solid.0);
    }

    #[test]
    fn check_self_interference_crossing_edges_is_true() {
        let mut pf = PaveFiller::new();
        let b = TopoBuilder::new();
        // Edge (0,0,0)->(1,1,0) and edge (0,1,0)->(1,0,0) cross at (0.5, 0.5).
        let e1 = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.0), dir(1.0, 1.0, 0.0)))),
            0.0,
            2.0f64.sqrt(),
        );
        let e2 = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.0, 1.0, 0.0), dir(1.0, -1.0, 0.0)))),
            0.0,
            2.0f64.sqrt(),
        );
        let comp = b.make_compound_of(&[e1.0.clone(), e2.0.clone()]);
        pf.set_arguments(&[comp.into()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(res, "crossing edges of one argument must be self-interfering");
        assert!(pf.has_warnings());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn check_self_interference_separate_arguments_is_false() {
        // Two overlapping boxes passed as TWO arguments are not acquired
        // self-interference — the overlap is the intended operation.
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        pf.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(!res, "inter-argument overlap is not self-interference");
        assert!(!pf.has_warnings());
        clear_tree(&a.solid.0);
        clear_tree(&c.solid.0);
    }

    // -----------------------------------------------------------------------
    // remove_micro_edges
    // -----------------------------------------------------------------------

    #[test]
    fn remove_micro_edges_removes_degenerate_blocks() {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();
        let e = pf.ds().index(&b.edges[0].0).expect("edge indexed");
        // A curve-less edge appended to the DS: any block on it fails the
        // shrunk-range computation, so it is a micro edge.
        let n_empty = pf.ds_mut().append(Edge::new().0).unwrap();
        let subs = pf.ds().shape_info(e).unwrap().sub_shapes().to_vec();
        let (v1, v2) = (subs[0], subs[1]);
        {
            let blocks = pf.ds_mut().change_pave_blocks_mut(e);
            // Block 1: normal (distinct bounds on a real edge).
            let mut pb1 = BopdsPaveBlock::new();
            pb1.set_edge(e);
            pb1.set_original_edge(e);
            pb1.set_range(0.0, 0.4);
            pb1.set_indices(v1, v2);
            // Block 2: degenerate (coincident bounds) on a curve-less edge.
            let mut pb2 = BopdsPaveBlock::new();
            pb2.set_edge(n_empty);
            pb2.set_original_edge(e);
            pb2.set_range(0.4, 0.6);
            pb2.set_indices(n_empty, n_empty);
            blocks.clear();
            blocks.push(pb1);
            blocks.push(pb2);
        }

        remove_micro_edges(&mut pf);

        let blocks = pf.ds().pave_blocks(e);
        assert!(
            !blocks.iter().any(|pb| pb.edge() == n_empty),
            "the degenerate block must be removed, got: {:?}",
            blocks
        );
        assert!(blocks.iter().any(|pb| pb.edge() == e), "the normal block is kept");
        clear_tree(&b.solid.0);
    }
}
