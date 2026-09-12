use super::prelude::*;
use super::*;

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
    let (t1, t2) = pb.range();
    let (n1, n2) = pb.indices();
    let v1 = f.ds().shape(n1).cloned().map(Vertex);
    let v2 = f.ds().shape(n2).cloned().map(Vertex);
    // `Err` is expected for degenerated / micro edges; `analyze_shrunk_data`
    // inspects `is_done()` and applies the OCCT warning policy.
    let _ = sr.set_data(edge, t1, t2, v1.as_ref(), v2.as_ref(), tol);
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
            pb.clear_shrunk_data();
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
    crate::pave_update_sd::update_edge_tolerance(f, n_e, tol)
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
    crate::pave_update_sd::update_vertex(f, n_v, tol)
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
/// Walks the typed `InterfVV/VE/VF/EE/EF` arrays on the filler DS and
/// re-points each stored `index_new` at its same-domain representative.
pub fn update_interfs_with_sd_vertices(f: &mut PaveFiller) -> Result<(), String> {
    f.ds_mut().update_interfs_with_sd_vertices();
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
    crate::pave_update_sd::update_common_blocks_with_sd_vertices(f)
}

// ---------------------------------------------------------------------------
// CheckSelfInterference
// ---------------------------------------------------------------------------

/// The shape-type pairs inspected for self-intersection.
#[allow(dead_code)]
pub(super) const SI_TYPE_PAIRS: [(ShapeType, ShapeType); 3] = [
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
    crate::pave_self::check_self_interference(f)
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
    crate::pave_ff_se::remove_micro_edges(f);
}

/// Whether the two sub-shapes of indices `i` and `j` genuinely self-intersect
/// (beyond sharing a boundary element).
#[allow(dead_code)]
pub(super) fn is_self_intersecting(ds: &BopdsDS, i: usize, j: usize) -> bool {
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
pub(super) fn bbox_overlap_dim(ds: &BopdsDS, i: usize, j: usize) -> usize {
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
pub(super) fn overlap_len(amin: f64, amax: f64, bmin: f64, bmax: f64) -> f64 {
    let lo = amin.max(bmin);
    let hi = amax.min(bmax);
    (hi - lo).max(0.0)
}

/// Whether the two edges (DS indices `i` and `j`) share a boundary vertex,
/// following same-domain chains.
pub(super) fn edges_share_boundary_vertex(ds: &BopdsDS, i: usize, j: usize) -> bool {
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
