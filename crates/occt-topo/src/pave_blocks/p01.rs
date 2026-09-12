use super::prelude::*;
use super::*;

/// The context a pave-filler exposes to the block-building logic.
///
/// Implemented by the concrete `crate::pave_filler::PaveFiller` once that
/// module is filled; the local test stub provides the same surface so this
/// module is verifiable in isolation.

pub trait PaveFillerLike {
    /// Read access to the data structure.
    fn ds(&self) -> &BopdsDS;
    /// Mutable access to the data structure.
    fn ds_mut(&mut self) -> &mut BopdsDS;
    /// Record a fatal error (the parallel filler accumulates these).
    fn add_error(&mut self, msg: String);
    /// Record a non-fatal warning.
    fn add_warning(&mut self, msg: String);
    /// True when p-curve building is disabled (`myAvoidBuildPCurve`).
    fn avoid_build_pcurve(&self) -> bool {
        false
    }
    /// Section attributes: whether p-curves are requested on the first/second
    /// face of each face/face pair (`PCurveOnS1` / `PCurveOnS2`).
    fn section_pcurve_on(&self) -> (bool, bool) {
        (true, true)
    }
}

// ---------------------------------------------------------------------------
// MakeBlocks
// ---------------------------------------------------------------------------

/// Ensure every source edge carries a default pave block, order each edge's
/// blocks by their parametric range, and split every block that accumulated
/// extra paves during intersection.
///
/// This is the reduced form of the *EE* `MakeBlocks` organiser: the face/face
/// section-edge construction lives in [`crate::pave_ff::make_blocks_ff`]
/// (`BOPAlgo_PaveFiller_6.cxx`). This step only *organizes* the pave blocks the
/// non-FF intersections already produced.
pub fn make_blocks<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    // 1. Make sure every source edge has a default pave block. This is a no-op
    //    for edges already carrying blocks (idempotent).
    for i in 0..n {
        let is_edge = f.ds().shape_info(i).map(|s| s.shape_type() == ShapeType::Edge).unwrap_or(false);
        if is_edge && !f.ds().has_pave_blocks(i) {
            f.ds_mut().init_pave_blocks_for_edge(i);
        }
    }
    // 2. Order each edge's blocks by their interval endpoints so adjacent
    //    intervals line up (needed before the splits below are examined).
    let edge_indices: Vec<usize> = (0..f.ds().nb_shapes())
        .filter(|&i| {
            f.ds().shape_info(i).map(|s| s.shape_type() == ShapeType::Edge).unwrap_or(false)
                && f.ds().has_pave_blocks(i)
        })
        .collect();
    for e in edge_indices {
        let blocks = f.ds_mut().change_pave_blocks_mut(e);
        blocks.sort_by(|a, b| {
            let (a1, a2) = a.range();
            let (b1, b2) = b.range();
            a1.partial_cmp(&b1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a2.partial_cmp(&b2).unwrap_or(std::cmp::Ordering::Equal))
        });
    }
    // 3. Split every block that carries extra paves into elementary blocks.
    f.ds_mut().update_pave_blocks();
    Ok(())
}

// ---------------------------------------------------------------------------
// MakeBlocks (full, with CommonBlock detection)
// ---------------------------------------------------------------------------

/// Full `MakeBlocks`-equivalent for the non-FF pipeline.
///
/// [`make_blocks`] is the reduced form: it only organizes the pave blocks the
/// intersections already produced. `make_blocks_full` additionally detects the
/// pave blocks of different source edges that cover the same geometric
/// interval and records them as [`crate::bopds::BopdsCommonBlock`]s, so later
/// stages can treat the coincident edges as one entity.
///
/// The detection is interval-based (the edge geometry is not re-sampled): two
/// blocks share an interval when they are bounded by the same pair of vertices
/// ([`BopdsPaveBlock::has_same_bounds`]), or when their endpoint parameters
/// coincide within `Precision::PConfusion()` (in either orientation). For each
/// such group a common block is recorded through
/// [`crate::bopds::BopdsDS::update_common_block`]; the same edge may appear in
/// several groups when it carries several distinct intervals.
///
/// Steps:
/// 1. per source edge, sort the pave blocks by range and merge adjacent blocks
///    whose ranges touch at a shared *source* bound vertex (blocks carrying
///    extra paves are never merged — their split points are kept);
/// 2. group the blocks of different edges sharing a geometric interval into
///    common blocks;
/// 3. split every block that accumulated extra paves
///    ([`crate::bopds::BopdsDS::update_pave_blocks`]).
pub fn make_blocks_full<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    let tol = PCONFUSION;
    let n = f.ds().nb_source_shapes();
    // Source edges that already carry pave blocks. Only these can take part in
    // a common block: an edge without blocks has no interval to share.
    let edge_indices: Vec<usize> = (0..n)
        .filter(|&i| {
            f.ds().shape_info(i).map(|s| s.shape_type() == ShapeType::Edge).unwrap_or(false)
                && f.ds().has_pave_blocks(i)
        })
        .collect();

    // 1. Organize each edge: sort, then merge adjacent touching blocks.
    for &e in &edge_indices {
        let blocks = f.ds_mut().change_pave_blocks_mut(e);
        let old = std::mem::take(blocks);
        *blocks = merge_touching_blocks(old, n, tol);
    }

    // 2. Group blocks of different edges sharing the same interval. The
    //    grouping runs on a read-only snapshot; the common blocks are applied
    //    to the DS only after all groups are collected, so the `used` fence
    //    stays aligned with the block lists it indexes.
    let mut used: Vec<Vec<bool>> = edge_indices
        .iter()
        .map(|&e| vec![false; f.ds().pave_blocks(e).len()])
        .collect();
    let mut groups: Vec<crate::bopds::BopdsCommonBlock> = Vec::new();
    for (i, &e1) in edge_indices.iter().enumerate() {
        let blocks1 = f.ds().pave_blocks(e1).to_vec();
        for (k1, pb1) in blocks1.iter().enumerate() {
            if used[i][k1] {
                continue;
            }
            let mut cb = crate::bopds::BopdsCommonBlock::new();
            cb.add_range(pb1.first, pb1.last);
            cb.add_index(e1);
            used[i][k1] = true;
            for (j, &e2) in edge_indices.iter().enumerate() {
                if j <= i {
                    continue;
                }
                let blocks2 = f.ds().pave_blocks(e2).to_vec();
                for (k2, pb2) in blocks2.iter().enumerate() {
                    if used[j][k2] {
                        continue;
                    }
                    if blocks_match(pb1, pb2, tol) {
                        cb.add_index(e2);
                        if !cb.contains_range(pb2.first, pb2.last, tol) {
                            cb.add_range(pb2.first, pb2.last);
                        }
                        used[j][k2] = true;
                    }
                }
            }
            if cb.indices().len() >= 2 {
                cb.set_tolerance(tol);
                groups.push(cb);
            }
        }
    }
    for cb in groups {
        f.ds_mut().update_common_block(&cb);
    }

    // 3. Split every block that accumulated extra paves.
    f.ds_mut().update_pave_blocks();
    Ok(())
}

/// True when two pave blocks of different edges cover the same geometric
/// interval: they are bounded by the same pair of vertices, or their endpoint
/// parameters coincide within `tol` (in either orientation).
pub(super) fn blocks_match(a: &BopdsPaveBlock, b: &BopdsPaveBlock, tol: f64) -> bool {
    if a.edge_index == b.edge_index {
        return false;
    }
    if a.has_same_bounds(b) {
        return true;
    }
    let (a1, a2) = a.range();
    let (b1, b2) = b.range();
    ((a1 - b1).abs() <= tol && (a2 - b2).abs() <= tol)
        || ((a1 - b2).abs() <= tol && (a2 - b1).abs() <= tol)
}

/// Merge adjacent blocks of one edge whose ranges touch at a shared *source*
/// bound vertex.
///
/// A block is merged into its predecessor when the predecessor's last bound
/// touches the block's first bound (within `tol`) at the same vertex index,
/// that vertex is a source shape (so the merge does not erase an intersection
/// split), and neither block carries extra paves (which would be lost by the
/// merge). The predecessor is extended to cover the union interval; the shared
/// bound vertex disappears from the representation.
pub(super) fn merge_touching_blocks(
    blocks: Vec<BopdsPaveBlock>,
    nb_source: usize,
    tol: f64,
) -> Vec<BopdsPaveBlock> {
    let mut blocks = blocks;
    blocks.sort_by(|a, b| {
        a.first
            .partial_cmp(&b.first)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.last.partial_cmp(&b.last).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut merged: Vec<BopdsPaveBlock> = Vec::new();
    for pb in blocks {
        if let Some(last) = merged.last_mut() {
            let touching = (last.last - pb.first).abs() <= tol;
            let shared_bound = last.index2 == pb.index1;
            let safe = last.index2 < nb_source && !last.is_to_update() && !pb.is_to_update();
            if touching && shared_bound && safe {
                last.set_pave2(pb.pave2());
                continue;
            }
        }
        merged.push(pb);
    }
    merged
}

// ---------------------------------------------------------------------------
// FilterPavesOnCurves
// ---------------------------------------------------------------------------

/// Filter the paves on the curves of the data structure: deduplicate extra
/// paves and merge adjacent elementary blocks that touch at the same vertex.
///
/// A reduced port of `BOPAlgo_PaveFiller::FilterPavesOnCurves` (which walks
/// the section curves produced by face/face intersection). Here the same
/// cleanup is applied to every edge's pave blocks:
///
/// 1. an extra pave whose parameter collides (within `Precision::PConfusion()`)
///    with a bound pave or with a sibling extra pave is dropped;
/// 2. two adjacent blocks that share the same *source* boundary vertex *and*
///    the same boundary parameter, and carry no extra paves of their own, are
///    merged back into a single block spanning the union interval.
///
/// The merge is restricted to *source* vertices so a genuine intersection
/// split — always delimited by a newly-created vertex — is never undone.
pub fn filter_paves_on_curves<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    let tol = PCONFUSION;
    let nb_source = f.ds().nb_source_shapes();
    let edge_indices: Vec<usize> = (0..f.ds().nb_shapes())
        .filter(|&i| {
            f.ds().shape_info(i).map(|s| s.shape_type() == ShapeType::Edge).unwrap_or(false)
                && f.ds().has_pave_blocks(i)
        })
        .collect();
    for e in edge_indices {
        let blocks = f.ds_mut().change_pave_blocks_mut(e);
        // 1. Deduplicate extra paves by parameter.
        for pb in blocks.iter_mut() {
            let (t1, t2) = pb.range();
            let mut kept: Vec<f64> = vec![t1, t2];
            let mut paves: Vec<BopdsPave> = Vec::new();
            for p in pb.ext_paves().iter().cloned() {
                if kept.iter().any(|&s| (s - p.param).abs() <= tol) {
                    continue;
                }
                kept.push(p.param);
                paves.push(p);
            }
            set_ext_paves(pb, paves);
        }
        // 2. Merge adjacent blocks touching at the same source vertex/parameter.
        let mut merged: Vec<BopdsPaveBlock> = Vec::new();
        for pb in blocks.drain(..) {
            if let Some(last) = merged.last_mut() {
                if last.ext_paves().is_empty() && pb.ext_paves().is_empty() {
                    let p2 = last.pave2();
                    let p1 = pb.pave1();
                    if p2.index < nb_source && p2.index == p1.index && (p2.param - p1.param).abs() <= tol {
                        last.set_pave2(pb.pave2());
                        continue;
                    }
                }
            }
            merged.push(pb);
        }
        *blocks = merged;
    }
    Ok(())
}

/// Replace the extra paves of a block, rebuilding its vertex fence.
pub(super) fn set_ext_paves(pb: &mut BopdsPaveBlock, paves: Vec<BopdsPave>) {
    pb.change_ext_paves().clear();
    pb.fence.clear();
    for p in paves {
        pb.append_ext_pave1(p);
        pb.fence.insert(p.index);
    }
}

// ---------------------------------------------------------------------------
// MakeSplitEdges
// ---------------------------------------------------------------------------

/// Split every edge whose pave blocks are delimited by intersection vertices.
///
/// Port of `BOPAlgo_PaveFiller::MakeSplitEdges`. For each pave block of each
/// edge, when at least one bound vertex is a *new* shape (created by an
/// intersection) or the edge carries more than one block, a new sub-edge is
/// built on the original edge's curve bounded by the block's vertices at the
/// block's parameters ([`AlgoTools::make_split_edge`]), appended to the data
/// structure, and the block is re-pointed at it. Blocks whose bounds are both
/// source shapes on a single-block edge are left on the original edge.
pub fn make_split_edges<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    pub(super) struct Task {
        /// DS index of the shape whose pave-block list holds the block.
        pub(super) holder: usize,
        /// Position of the block inside the holder's list.
        pub(super) pos: usize,
        /// Original edge the split is taken from.
        pub(super) original: usize,
        pub(super) v1: usize,
        pub(super) t1: f64,
        pub(super) v2: usize,
        pub(super) t2: f64,
    }

    let mut tasks: Vec<Task> = Vec::new();
    let n = f.ds().nb_shapes();
    for i in 0..n {
        let Some(info) = f.ds().shape_info(i) else { continue };
        if info.shape_type() != ShapeType::Edge {
            continue;
        }
        if !f.ds().has_pave_blocks(i) {
            continue;
        }
        let blocks = f.ds().pave_blocks(i).to_vec();
        if blocks.is_empty() {
            continue;
        }
        // Skip degenerated edges.
        let original = blocks[0].original_edge();
        let orig_is_degenerated = f
            .ds()
            .shape(original)
            .map(|s| BRepTool::is_degenerated(&Edge(s.clone())))
            .unwrap_or(false);
        if orig_is_degenerated {
            continue;
        }
        let multi = blocks.len() > 1;
        for (k, pb) in blocks.iter().enumerate() {
            let (n_v1, n_v2) = pb.indices();
            let (a_t1, a_t2) = pb.range();
            // A zero-length block (two coincident bound vertices) cannot bound
            // a split piece — it would create a degenerate sub-edge. This can
            // happen when a new intersection vertex coincides with an existing
            // source vertex (e.g. a box corner on the other box's edge).
            if (a_t2 - a_t1).abs() <= PCONFUSION {
                continue;
            }
            // Already a split of this block (`aPB->Edge()` after MakeSplitEdges).
            if pb.edge() != 0 && pb.edge() != pb.original_edge() {
                continue;
            }
            let need_split = multi || f.ds().is_new_shape(n_v1) || f.ds().is_new_shape(n_v2);
            if need_split {
                tasks.push(Task {
                    holder: i,
                    pos: k,
                    original,
                    v1: n_v1,
                    t1: a_t1,
                    v2: n_v2,
                    t2: a_t2,
                });
            }
        }
    }

    for t in tasks {
        let e_shape = match f.ds().shape(t.original) {
            Some(s) => s.clone(),
            None => {
                let msg = format!("make_split_edges: edge {} not found", t.original);
                f.add_error(msg.clone());
                return Err(msg);
            }
        };
        let v1 = match f.ds().shape(t.v1) {
            Some(s) => s.clone(),
            None => return Err(format!("make_split_edges: vertex {} not found", t.v1)),
        };
        let v2 = match f.ds().shape(t.v2) {
            Some(s) => s.clone(),
            None => return Err(format!("make_split_edges: vertex {} not found", t.v2)),
        };
        let sp = AlgoTools::make_split_edge(&Edge(e_shape), Some(&v1), t.t1, Some(&v2), t.t2)
            .map_err(|e| {
                let msg = format!("make_split_edges: {e}");
                f.add_error(msg.clone());
                msg
            })?;
        let mut si = BopdsShapeInfo::new(sp.0);
        si.change_sub_shapes().extend_from_slice(&[t.v1, t.v2]);
        let n_sp = f.ds_mut().append_info(si);
        let blocks = f.ds_mut().change_pave_blocks_mut(t.holder);
        if let Some(pb) = blocks.get_mut(t.pos) {
            pb.set_edge(n_sp);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// MakePCurves
// ---------------------------------------------------------------------------

/// Build and attach the p-curve of every edge stored in the face-info pool
/// onto its face.
///
/// Port of `BOPAlgo_PaveFiller::MakePCurves` (`BOPAlgo_PaveFiller_7.cxx:589`),
/// both branches:
///
/// 1. **IN / ON pave blocks** — the faces' boundary (`On`) and coincident
///    (`In`) edges. Each edge without a p-curve on the face yet is projected
///    with [`pcurve_full::make_pcurve_full`] and stored in the
///    [`GeometryRegistry`] side-table via `set_edge_pcurve`. An `On` edge that
///    already carries a p-curve is skipped, matching the OCCT `bHasPC` test
///    (OCCT step 1: the `PaveBlocksIn` / `PaveBlocksOn` MPC batches).
/// 2. **Section pave blocks** — the face/face section edges. Their p-curves
///    are already attached by the F/F intersection ([`crate::pave_intersect`]),
///    so only the endpoint vertex tolerances are grown to cover the 3D/2D
///    deviation. This is the OCCT step 2, where the section MPCs carry the
///    section flag and `MPC::Perform` calls `UpdateVertices` even though no
///    new p-curve is built. When a section edge has no p-curve yet (a
///    degenerate pipeline state), it is built first, as OCCT does.
///
/// Translation boundaries: the common-block p-curve copy of the OCCT `On`
/// branch (`BOPTools_AlgoTools2D::AttachExistingPCurve`, reusing a coincident
/// edge's p-curve) is not ported — recomputation via [`pcurve_full`] is
/// geometrically equivalent for coincident edges and the Rust common block
/// stores edge indices + ranges rather than pave-block handles. The periodic
/// p-curve adjustment (`AdjustPCurveOnSurf`) of the `MPC::Perform` fast path
/// is likewise left to the F/F solver / trim step.
pub fn make_pcurves<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    if f.avoid_build_pcurve() {
        return Ok(());
    }
    let (on_s1, on_s2) = f.section_pcurve_on();
    if !on_s1 && !on_s2 {
        return Ok(());
    }
    let fi_pool = f.ds().face_info_pool().to_vec();
    for fi in fi_pool {
        let Some(face_shape) = f.ds().shape(fi.face_index).cloned() else { continue };
        let face = Face(face_shape);
        let face_key = GeometryRegistry::shape_key(&face.0);

        // 1. IN + ON pave blocks: build a p-curve for every edge that does not
        //    carry one on this face yet.
        let mut edges: Vec<usize> = fi.paves_in().iter().map(|&(e, _, _)| e).collect();
        for &(e, _, _) in fi.paves_on() {
            // OCCT skips the On blocks whose edge already has a curve on the face.
            let has_pc = f
                .ds()
                .shape(e)
                .map(|s| boptools_2d::curve_on_surface(&Edge(s.clone()), &face).is_some())
                .unwrap_or(false);
            if !has_pc {
                edges.push(e);
            }
        }
        edges.sort_unstable();
        edges.dedup();
        for edge_idx in edges {
            let Some(edge_shape) = f.ds().shape(edge_idx).cloned() else { continue };
            let edge = Edge(edge_shape);
            if boptools_2d::curve_on_surface(&edge, &face).is_some() {
                continue;
            }
            match pcurve_full::make_pcurve_full(&edge, &face) {
                Ok(pc) => {
                    GeometryRegistry::global().set_edge_pcurve(&edge.0, face_key, pc);
                }
                Err(e) => {
                    let msg = format!(
                        "make_pcurves: pcurve failed for edge {edge_idx} on face {}: {e}",
                        fi.face_index
                    );
                    f.add_warning(msg);
                }
            }
        }

        // 2. Section pave blocks: the p-curve is already attached by the F/F
        //    intersection — only the endpoint vertex tolerances are grown
        //    (`UpdateVertices`). An edge without a p-curve (degenerate pipeline
        //    state) gets one built first, as `MPC::Perform` does.
        if on_s1 || on_s2 {
            let mut section_edges: Vec<usize> = Vec::new();
            for &(edge_idx, _t1, _t2) in fi.paves() {
                if !section_edges.contains(&edge_idx) {
                    section_edges.push(edge_idx);
                }
            }
            for edge_idx in section_edges {
                let Some(edge_shape) = f.ds().shape(edge_idx).cloned() else { continue };
                let edge = Edge(edge_shape);
                let pc = match boptools_2d::curve_on_surface(&edge, &face) {
                    Some(pc) => pc,
                    None => match pcurve_full::make_pcurve_full(&edge, &face) {
                        Ok(pc) => {
                            GeometryRegistry::global().set_edge_pcurve(&edge.0, face_key, pc.clone());
                            pc
                        }
                        Err(e) => {
                            let msg = format!(
                                "make_pcurves: pcurve failed for edge {edge_idx} on face {}: {e}",
                                fi.face_index
                            );
                            f.add_warning(msg);
                            continue;
                        }
                    },
                };
                update_vertices(&edge, &face, &pc);
            }
        }
    }
    Ok(())
}

/// Grow the endpoint-vertex tolerances of `edge` so each vertex covers the
/// deviation between the 3D curve point and the surface point mapped by the
/// p-curve at the edge's boundary parameter (OCCT `UpdateVertices`).
pub(super) fn update_vertices(edge: &Edge, face: &Face, pc: &Arc<dyn Curve2d>) {
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else { return };
    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else { return };
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() {
        return;
    }
    let verts = edge_vertex_shapes(edge);
    for (j, t) in [(0usize, a), (1, b)] {
        if let Some(v) = verts.get(j) {
            let vtx = Vertex(v.clone());
            let tol = BRepTool::vertex_tolerance(&vtx);
            let p3d = c3d.d0(t);
            let q = pc.d0(t);
            let p3dx = surf.d0(q.x(), q.y());
            let d = p3d.distance(&p3dx);
            if d > tol {
                vtx.set_tolerance(d + crate::algo_tools::D_TOLERANCE);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// MakePCurves (full: trim + endpoint alignment + 2D tolerance)
// ---------------------------------------------------------------------------

/// Deepened [`make_pcurves`] — build the pcurve of every edge that needs one
/// on each adjacent face, then run the three OCCT post-processing steps:
///
/// * **Trim** — [`pcurve_full::trim_pcurve_to_face`] brings the pcurve inside
///   the face UV rectangle (whole-period shifts on periodic dimensions, a
///   clipped B-spline otherwise).
/// * **Align** — the pcurve endpoints are compared against the surface
///   projection of the edge's 3D endpoints; when the 2D deviation exceeds the
///   tolerance the pcurve is re-fitted so the endpoints land on the projected
///   values (`BRepLib::SameParameter`-style endpoint correction).
/// * **2D tolerance** — each boundary vertex's tolerance is grown to cover the
///   face tolerance so 2D comparisons on the face are consistent (OCCT
///   `UpdateVertices` plus a face-tolerance floor).
///
/// The edge/face pairs come from the [`crate::bopds::BopdsFaceInfo`] pool:
/// every source edge listed there plus the split edges its pave blocks now
/// reference, so each BOPDS edge carrying a p-curve requirement is covered.
pub fn make_pcurves_full<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    if f.avoid_build_pcurve() {
        return Ok(());
    }
    let (on_s1, on_s2) = f.section_pcurve_on();
    if !on_s1 && !on_s2 {
        return Ok(());
    }
    // 1. Collect every (edge, face) pair needing a pcurve: the source edges of
    //    the face-info pool plus the split edges reachable from their pave
    //    blocks (the BFS re-points to sub-edges when a split edge was itself
    //    split again).
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for fi in f.ds().face_info_pool() {
        let face = fi.face_index;
        let mut stack: Vec<usize> = fi.paves.iter().map(|&(e, _, _)| e).collect();
        let mut local: std::collections::HashSet<usize> = std::collections::HashSet::new();
        while let Some(e) = stack.pop() {
            if !local.insert(e) {
                continue;
            }
            if seen.insert((e, face)) {
                pairs.push((e, face));
            }
            for pb in f.ds().pave_blocks(e) {
                let n = pb.edge();
                if n != e
                    && n < f.ds().nb_shapes()
                    && f.ds()
                        .shape_info(n)
                        .map(|s| s.shape_type() == ShapeType::Edge)
                        .unwrap_or(false)
                {
                    stack.push(n);
                }
            }
        }
    }
    // 2. Build, trim, align, store and update tolerances for each pair.
    for (edge_idx, face_idx) in pairs {
        let Some(edge_shape) = f.ds().shape(edge_idx).cloned() else { continue };
        let Some(face_shape) = f.ds().shape(face_idx).cloned() else { continue };
        let edge = Edge(edge_shape);
        let face = Face(face_shape);
        if boptools_2d::curve_on_surface(&edge, &face).is_some() {
            continue;
        }
        let face_tol = BRepTool::face_tolerance(&face);
        let pc = match pcurve_full::make_pcurve_full(&edge, &face) {
            Ok(pc) => pc,
            Err(e) => {
                let msg = format!(
                    "make_pcurves_full: pcurve failed for edge {edge_idx} on face {face_idx}: {e}"
                );
                f.add_warning(msg);
                continue;
            }
        };
        let pc = match pcurve_full::trim_pcurve_to_face(&pc, &face, face_tol) {
            Ok(pc) => pc,
            Err(e) => {
                let msg = format!(
                    "make_pcurves_full: trim failed for edge {edge_idx} on face {face_idx}: {e}"
                );
                f.add_warning(msg);
                continue;
            }
        };
        let pc = match align_pcurve_endpoints(&pc, &edge, &face, face_tol) {
            Ok(pc) => pc,
            Err(e) => {
                let msg = format!(
                    "make_pcurves_full: alignment failed for edge {edge_idx} on face {face_idx}: {e}"
                );
                f.add_warning(msg);
                continue;
            }
        };
        let face_key = GeometryRegistry::shape_key(&face.0);
        GeometryRegistry::global().set_edge_pcurve(&edge.0, face_key, pc.clone());
        update_vertices_full(&edge, &face, &pc, face_tol);
    }
    Ok(())
}
