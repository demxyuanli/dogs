//! Pave-block building and edge splitting — Phase 19 wave C2b.
//!
//! Port of the `BOPAlgo_PaveFiller` block-building methods:
//!
//! | OCCT method                   | Rust function                          |
//! |-------------------------------|----------------------------------------|
//! | `MakeBlocks`                  | [`make_blocks`]                        |
//! | `MakeBlocks` (full)          | [`make_blocks_full`]                   |
//! | `FilterPavesOnCurves`         | [`filter_paves_on_curves`]             |
//! | `MakeSplitEdges`              | [`make_split_edges`]                   |
//! | `MakePCurves`                 | [`make_pcurves`]                       |
//! | `FindPaveBlocks`              | [`find_pave_blocks`]                   |
//! | `FillPaves`                   | [`fill_paves`]                         |
//! | `SplitEdge`                   | [`make_split_edge`]                    |
//! | `MakeSplitEdge` (ProcessDE)   | [`make_split_edge_de`]                 |
//! | `ProcessDE`                   | [`process_de`]                         |
//!
//! `MakeBlocks` in OCCT builds the face/face section edges (the `PerformFF`
//! step). That step is **not** part of this phase, so [`make_blocks`] is the
//! reduced form: it makes sure every source edge has a default pave block,
//! orders each edge's blocks by their parametric range, then splits any block
//! carrying extra paves (the OCCT `BOPDS_DS::UpdatePaveBlocks`). The FF-only
//! helpers (`CorrectToleranceOfSE`, `PutPavesOnCurve`, …) are omitted.
//!
//! Every function operates on a [`PaveFillerLike`] — the minimal contract the
//! concrete `crate::pave_filler::PaveFiller` (filled by a parallel agent)
//! satisfies. The port keeps the logic decoupled from that type so the two
//! can be developed and verified independently.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api::{intersect_curves, project_point_on_curve};
use occt_geom2d::Geom2dBSplineCurve;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsDS, BopdsPave, BopdsPaveBlock, BopdsShapeInfo};
use crate::boptools_2d;
use crate::brep_projection::project_point_on_face;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::pcurve_full;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

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
/// This is the reduced form of `BOPAlgo_PaveFiller::MakeBlocks` for the
/// non-FF pipeline: the face/face section-edge construction is not part of
/// this phase, so the step reduces to *organizing* the pave blocks the
/// intersections already produced.
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
fn blocks_match(a: &BopdsPaveBlock, b: &BopdsPaveBlock, tol: f64) -> bool {
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
fn merge_touching_blocks(
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
fn set_ext_paves(pb: &mut BopdsPaveBlock, paves: Vec<BopdsPave>) {
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
    struct Task {
        /// DS index of the shape whose pave-block list holds the block.
        holder: usize,
        /// Position of the block inside the holder's list.
        pos: usize,
        /// Original edge the split is taken from.
        original: usize,
        v1: usize,
        t1: f64,
        v2: usize,
        t2: f64,
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
fn update_vertices(edge: &Edge, face: &Face, pc: &Arc<dyn Curve2d>) {
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

/// Align the pcurve endpoints with the surface projection of the edge's 3D
/// endpoints (`BRepLib::SameParameter`-style endpoint correction).
///
/// When the pcurve spans the full edge range `[a, b]`, each endpoint is
/// compared against the projection of the edge's 3D curve point at `a` / `b`
/// onto the face surface. If the 2D deviation exceeds `tol`, the pcurve is
/// re-fitted as a degree-1 B-spline over `[a, b]` whose endpoints are the
/// projected values. A pcurve trimmed to a proper sub-range already ends on
/// the face boundary and is returned unchanged; periodic (circle) pcurves are
/// also left alone.
fn align_pcurve_endpoints(
    pc: &Arc<dyn Curve2d>,
    edge: &Edge,
    face: &Face,
    tol: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    if pc.is_periodic() {
        return Ok(pc.clone());
    }
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Ok(pc.clone());
    }
    let (p0, p1) = (pc.first_parameter(), pc.last_parameter());
    // A bounded pcurve that no longer spans the edge range was clipped by the
    // trim step — its endpoints are the intended face-boundary points.
    if p0.is_finite() && p1.is_finite() && ((a - p0).abs() > 1e-9 || (b - p1).abs() > 1e-9) {
        return Ok(pc.clone());
    }
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else { return Ok(pc.clone()) };
    let q0 = pc.d0(a);
    let q1 = pc.d0(b);
    let Some((u0, v0, _)) = project_point_on_face(face, &c3d.d0(a)) else { return Ok(pc.clone()) };
    let Some((u1, v1, _)) = project_point_on_face(face, &c3d.d0(b)) else { return Ok(pc.clone()) };
    let t0t = GpPnt2d::new(u0, v0);
    let t1t = GpPnt2d::new(u1, v1);
    if q0.distance(&t0t) <= tol && q1.distance(&t1t) <= tol {
        return Ok(pc.clone());
    }
    // Re-fit a degree-1 B-spline over [a, b] with the corrected endpoints.
    let n = 33;
    let mut pts: Vec<GpPnt2d> = (0..n)
        .map(|i| pc.d0(a + (b - a) * i as f64 / (n - 1) as f64))
        .collect();
    pts[0] = t0t;
    pts[n - 1] = t1t;
    let bs = bspline_from_pts_2d(&pts, a, b, 1)?;
    Ok(Arc::new(bs))
}

/// Fit a degree-`degree` clamped B-spline through `pts` over `[a, b]` — a
/// polyline for degree 1, mirroring the fitter used by
/// [`crate::pcurve_full`].
fn bspline_from_pts_2d(
    pts: &[GpPnt2d],
    a: f64,
    b: f64,
    degree: usize,
) -> Result<Geom2dBSplineCurve, String> {
    if pts.len() < 2 {
        return Err("bspline_from_pts_2d: need at least 2 points".into());
    }
    if degree >= pts.len() {
        return Err("bspline_from_pts_2d: degree below pole count".into());
    }
    let mut knots = occt_core::bspl::knots::build_uniform_knots(pts.len(), degree);
    for k in knots.iter_mut() {
        *k = a + (b - a) * *k;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots, degree).map_err(|e| e.to_string())
}

/// Update the endpoint-vertex tolerances of `edge` after a pcurve was built on
/// `face` (OCCT `UpdateVertices`): grow each vertex to cover the deviation
/// between the 3D curve point and the surface point mapped by the pcurve at
/// the boundary parameter, and floor it at the face 2D tolerance.
fn update_vertices_full(edge: &Edge, face: &Face, pc: &Arc<dyn Curve2d>, face_tol: f64) {
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
            let mut want = tol;
            if d > tol {
                want = d + crate::algo_tools::D_TOLERANCE;
            }
            if face_tol > want {
                want = face_tol;
            }
            if want > tol {
                vtx.set_tolerance(want);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ProcessDE helpers
// ---------------------------------------------------------------------------

/// Find the pave blocks of the edges of the face `face` that pass through the
/// vertex with DS index `v`.
///
/// Port of `BOPAlgo_PaveFiller::FindPaveBlocks`. The face's pave list (edge
/// index + parameter range) is scanned and every block of those edges whose
/// first or last bound references `v` is appended to `out`.
pub fn find_pave_blocks<F: PaveFillerLike>(f: &F, v: usize, face: usize, out: &mut Vec<BopdsPaveBlock>) {
    let Some(fi) = f.ds().face_info_pool().iter().find(|fi| fi.face_index == face) else {
        return;
    };
    for (edge, _t1, _t2) in &fi.paves {
        let blocks = f.ds().pave_blocks(*edge);
        for pb in blocks {
            let (n1, n2) = pb.indices();
            if n1 == v || n2 == v {
                out.push(pb.clone());
            }
        }
    }
}

/// Add split points to the pave block `pbd` of a degenerated edge: intersect
/// the 2D curve of the degenerated edge with the 2D curves of every passing
/// pave block and record the intersection parameters as extra paves.
///
/// Port of `BOPAlgo_PaveFiller::FillPaves`. When the 2D intersection produces
/// no point, the endpoint of the passing curve that corresponds to `v` is
/// projected onto the degenerated curve instead (the OCCT fallback).
pub fn fill_paves<F: PaveFillerLike>(
    f: &mut F,
    v: usize,
    e: usize,
    face: usize,
    lpb: &[BopdsPaveBlock],
    pbd: &mut BopdsPaveBlock,
) {
    let Some(de_shape) = f.ds().shape(e).cloned() else { return };
    let Some(face_shape) = f.ds().shape(face).cloned() else { return };
    let de = Edge(de_shape);
    let fa = Face(face_shape);
    let Ok(c2d_de) = boptools_2d::make_2d(&de, &fa) else { return };
    let tol = PCONFUSION;
    for pb in lpb {
        let n_e = pb.edge();
        if n_e >= f.ds().nb_shapes() {
            continue;
        }
        let Some(pb_shape) = f.ds().shape(n_e).cloned() else { continue };
        let pbe = Edge(pb_shape);
        let Ok(c2d) = boptools_2d::make_2d(&pbe, &fa) else { continue };
        let hits = intersect_curves(c2d_de.as_ref(), c2d.as_ref(), tol);
        if hits.is_empty() {
            let t = if v == pb.pave1().index() {
                pb.pave1().parameter()
            } else {
                pb.pave2().parameter()
            };
            let p2d = c2d.d0(t);
            if let Some(proj) = project_point_on_curve(c2d_de.as_ref(), &p2d, tol) {
                add_split_point(pbd, BopdsPave::new(v, proj.parameter), tol);
            }
        } else {
            for hit in hits {
                add_split_point(pbd, BopdsPave::new(v, hit.u1), tol);
            }
        }
    }
}

/// Validate and add `pave` as an extra pave of `pbd`: the parameter must be
/// strictly inside the block range and not collide with an existing pave.
/// Returns true when the pave was added (OCCT `AddSplitPoint`).
fn add_split_point(pbd: &mut BopdsPaveBlock, pave: BopdsPave, tol: f64) -> bool {
    let (td1, td2) = pbd.range();
    let t = pave.parameter();
    if t - td1 < tol || td2 - t < tol {
        return false;
    }
    if pbd.contains_parameter(t, tol).is_some() {
        return false;
    }
    pbd.append_ext_pave1(pave);
    true
}

/// Split the degenerated edge `de` on the face `df`, creating a new
/// (degenerated) sub-edge for each of its pave blocks delimited by a new
/// vertex. Port of `BOPAlgo_PaveFiller::MakeSplitEdge` (ProcessDE branch).
pub fn make_split_edge_de<F: PaveFillerLike>(f: &mut F, de: usize, df: usize) -> Result<(), String> {
    let Some(de_shape) = f.ds().shape(de).cloned() else { return Ok(()) };
    let blocks = f.ds().pave_blocks(de).to_vec();
    if blocks.is_empty() {
        return Ok(());
    }
    let multi = blocks.len() > 1;
    let mut split_indices: Vec<Option<usize>> = Vec::with_capacity(blocks.len());
    for pb in &blocks {
        let (n_v1, n_v2) = pb.indices();
        let (a_t1, a_t2) = pb.range();
        if f.ds().is_new_shape(n_v1) || multi {
            let v1 = f.ds().shape(n_v1).cloned().ok_or_else(|| format!("make_split_edge_de: vertex {n_v1} not found"))?;
            let v2 = f.ds().shape(n_v2).cloned().ok_or_else(|| format!("make_split_edge_de: vertex {n_v2} not found"))?;
            let sp = make_degenerate_split_edge(&Edge(de_shape.clone()), &v1, a_t1, &v2, a_t2)?;
            let mut si = BopdsShapeInfo::new(sp.0);
            si.change_sub_shapes().extend_from_slice(&[n_v1, n_v2]);
            let n_sp = f.ds_mut().append_info(si);
            split_indices.push(Some(n_sp));
        } else {
            // No split needed: drop the block (the whole list is cleared below).
            split_indices.push(None);
        }
    }
    // Clear the block list of the degenerated edge and re-point surviving blocks.
    if split_indices.iter().any(|s| s.is_some()) {
        let blocks = f.ds_mut().change_pave_blocks_mut(de);
        for (pb, si) in blocks.iter_mut().zip(split_indices.iter()) {
            if let Some(n_sp) = si {
                pb.set_edge(*n_sp);
            }
        }
    } else {
        {
            let si = f.ds_mut().change_shape_info(de).expect("de shape info");
            si.pb_reference = -1;
        }
        f.ds_mut().change_pave_blocks_mut(de).clear();
    }
    let _ = df;
    Ok(())
}

/// Build a degenerated sub-edge from the base edge `orig` bounded by `v1`/`v2`
/// (OCCT `MakeSplitEdge1`): a collapsed copy carrying the two vertices and the
/// degenerated flag.
fn make_degenerate_split_edge(
    orig: &Edge,
    v1: &TopoShape,
    t1: f64,
    v2: &TopoShape,
    t2: f64,
) -> Result<Edge, String> {
    let b = TopoBuilder::new();
    let p = if let Some(c) = BRepTool::edge_curve(orig) {
        c.d0(0.5 * (t1 + t2))
    } else {
        GpPnt::zero()
    };
    let mut e = b.make_edge_segment(&p, &p);
    if t1 < t2 {
        b.add(&mut e.0, &v1.oriented(Orientation::Forward));
        b.add(&mut e.0, &v2.oriented(Orientation::Reversed));
    } else {
        b.add(&mut e.0, &v1.oriented(Orientation::Reversed));
        b.add(&mut e.0, &v2.oriented(Orientation::Forward));
    }
    if let Some(mut g) = GeometryRegistry::global().edge_geom(&e.0) {
        g.degenerated = true;
        g.tolerance = CONFUSION;
        GeometryRegistry::global().set_edge(&e.0, g);
    }
    Ok(e)
}

/// Handle degenerated (collapsed) edges: collect split points on their pave
/// blocks from the edges crossing the degenerated edge's pole vertex, then
/// split them into degenerated sub-edges.
///
/// Port of `BOPAlgo_PaveFiller::ProcessDE`, adapted to the port's data model
/// (degenerated edges are detected through the [`GeometryRegistry`] flag rather
/// than the OCCT `HasFlag` sub-shape tag).
pub fn process_de<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    let de_edges: Vec<usize> = (0..n)
        .filter(|&i| {
            let is_edge = f.ds().shape_info(i).map(|s| s.shape_type() == ShapeType::Edge).unwrap_or(false);
            is_edge && f.ds().shape(i).map(|s| BRepTool::is_degenerated(&Edge(s.clone()))).unwrap_or(false)
        })
        .collect();
    for e in de_edges {
        let Some(e_shape) = f.ds().shape(e).cloned() else { continue };
        let verts = edge_vertex_shapes(&Edge(e_shape));
        let Some(n_v) = verts.first().and_then(|v| f.ds().index(v)) else { continue };
        let n_v = f.ds().get_same_domain_index(n_v);
        // Faces adjacent to the degenerated edge (their sub-shape list contains it).
        let faces: Vec<usize> = (0..n)
            .filter(|&fi| {
                let si = f.ds().shape_info(fi);
                si.map_or(false, |s| s.shape_type() == ShapeType::Face && s.sub_shapes().contains(&e))
            })
            .collect();
        for n_f in faces {
            let mut lpb: Vec<BopdsPaveBlock> = Vec::new();
            find_pave_blocks(f, n_v, n_f, &mut lpb);
            if !lpb.is_empty() {
                let pbd = {
                    let blocks = f.ds_mut().change_pave_blocks_mut(e);
                    let Some(first) = blocks.first().cloned() else { continue };
                    first
                };
                let mut pbd = pbd;
                fill_paves(f, n_v, e, n_f, &lpb, &mut pbd);
                let mut out = Vec::new();
                pbd.update(&mut out, true);
                let blocks = f.ds_mut().change_pave_blocks_mut(e);
                blocks.clear();
                blocks.extend(out);
            }
            make_split_edge_de(f, e, n_f)?;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// SplitEdge (single-edge split helper)
// ---------------------------------------------------------------------------

/// Create a new sub-edge of `edge` bounded by the vertices `v1`/`v2` at the
/// parameters `t1`/`t2` and append it to the data structure.
///
/// Port of `BOPAlgo_PaveFiller::SplitEdge`. Returns the DS index of the new
/// edge.
pub fn make_split_edge<F: PaveFillerLike>(
    f: &mut F,
    edge: usize,
    v1: usize,
    t1: f64,
    v2: usize,
    t2: f64,
) -> Result<usize, String> {
    let e_shape = f
        .ds()
        .shape(edge)
        .ok_or_else(|| format!("make_split_edge: edge {edge} not found"))?
        .clone();
    let v1_shape = f
        .ds()
        .shape(v1)
        .ok_or_else(|| format!("make_split_edge: vertex {v1} not found"))?
        .clone();
    let v2_shape = f
        .ds()
        .shape(v2)
        .ok_or_else(|| format!("make_split_edge: vertex {v2} not found"))?
        .clone();
    let sp = AlgoTools::make_split_edge(&Edge(e_shape), Some(&v1_shape), t1, Some(&v2_shape), t2)?;
    let mut si = BopdsShapeInfo::new(sp.0);
    si.change_sub_shapes().extend_from_slice(&[v1, v2]);
    Ok(f.ds_mut().append_info(si))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// The vertex child shapes of an edge (its two boundary vertices).
fn edge_vertex_shapes(edge: &Edge) -> Vec<TopoShape> {
    edge.0
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == ShapeType::Vertex)
        .cloned()
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;

    /// Minimal filler satisfying [`PaveFillerLike`] for isolated tests.
    #[derive(Debug, Default)]
    struct StubFiller {
        ds: BopdsDS,
        errors: Vec<String>,
        warnings: Vec<String>,
    }

    impl PaveFillerLike for StubFiller {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn ds_mut(&mut self) -> &mut BopdsDS {
            &mut self.ds
        }
        fn add_error(&mut self, msg: String) {
            self.errors.push(msg);
        }
        fn add_warning(&mut self, msg: String) {
            self.warnings.push(msg);
        }
    }

    /// A stub DS holding a unit box whose edge 0 runs (0,0,0) → (1,0,0).
    /// Uses `init` so the source/new-shape split (`nb_source_shapes`) is set.
    fn box_filler() -> (StubFiller, usize) {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        (f, e0)
    }

    /// Append a fresh vertex to the DS and return its index.
    fn push_vertex(f: &mut StubFiller, p: GpPnt) -> usize {
        let v = AlgoTools::make_new_vertex(&p, 1e-7).unwrap();
        f.ds.append(v).unwrap()
    }

    /// Build a pave block on `edge` bounded by `v1`/`v2` over `(t1, t2)`.
    fn mk_block(edge: usize, v1: usize, t1: f64, v2: usize, t2: f64) -> BopdsPaveBlock {
        let mut pb = BopdsPaveBlock::new();
        pb.set_edge(edge);
        pb.set_original_edge(edge);
        pb.set_pave1(BopdsPave::new(v1, t1));
        pb.set_pave2(BopdsPave::new(v2, t2));
        pb
    }

    /// A filler whose DS holds a unit box, with default pave blocks created on
    /// the box edges at the given indices. Returns the filler and the DS
    /// indices of the requested edges.
    fn box_edges(indices: &[usize]) -> (StubFiller, Vec<usize>) {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let inds: Vec<usize> = indices.iter().map(|&i| f.ds.index(&b.edges[i].0).unwrap()).collect();
        for &e in &inds {
            f.ds.init_pave_blocks_for_edge(e);
        }
        (f, inds)
    }

    #[test]
    fn make_blocks_splits_edge_with_two_extra_paves() {
        let (mut f, e0) = box_filler();
        // Edge 0 has one default block [0, 1]. Add two intersection paves.
        let i1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let i2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            assert_eq!(blocks.len(), 1);
            blocks[0].append_ext_pave(BopdsPave::new(i1, 0.25));
            blocks[0].append_ext_pave(BopdsPave::new(i2, 0.75));
            assert!(blocks[0].is_to_update());
        }

        make_blocks(&mut f).unwrap();

        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3, "expected 3 elementary blocks");
        assert_eq!(blocks[0].range(), (0.0, 0.25));
        assert_eq!(blocks[1].range(), (0.25, 0.75));
        assert_eq!(blocks[2].range(), (0.75, 1.0));
        // The middle block is bounded by the two new vertices.
        assert_eq!(blocks[1].indices(), (i1, i2));
    }

    #[test]
    fn make_split_edge_creates_sub_edge_with_correct_range() {
        let (mut f, e0) = box_filler();
        let v1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let v2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        let n_sp = make_split_edge(&mut f, e0, v1, 0.25, v2, 0.75).expect("split edge");
        assert!(n_sp >= f.ds.nb_source_shapes(), "new edge is a new shape");
        let new_shape = f.ds.shape(n_sp).cloned().expect("new edge shape");
        assert_eq!(new_shape.shape_type(), ShapeType::Edge);
        assert_eq!(BRepTool::edge_parameters(&Edge(new_shape)), (0.25, 0.75));
        // Sub-shapes of the new edge are the two bound vertices.
        let si = f.ds.shape_info(n_sp).unwrap();
        assert_eq!(si.sub_shapes(), &[v1, v2]);
        // 3D endpoints land on the original curve.
        let p0 = AlgoTools::point_on_edge(&Edge(f.ds.shape(n_sp).cloned().unwrap()), 0.25).unwrap();
        let p1 = AlgoTools::point_on_edge(&Edge(f.ds.shape(n_sp).cloned().unwrap()), 0.75).unwrap();
        assert!(p0.distance(&GpPnt::new(0.25, 0.0, 0.0)) < 1e-9, "start {p0:?}");
        assert!(p1.distance(&GpPnt::new(0.75, 0.0, 0.0)) < 1e-9, "end {p1:?}");
    }

    #[test]
    fn make_split_edges_splits_every_block_onto_new_edge() {
        let (mut f, e0) = box_filler();
        let i1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let i2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(i1, 0.25));
            blocks[0].append_ext_pave(BopdsPave::new(i2, 0.75));
        }
        make_blocks(&mut f).unwrap();
        let before = f.ds.nb_shapes();
        make_split_edges(&mut f).unwrap();
        assert_eq!(f.ds.nb_shapes(), before + 3, "3 sub-edges appended");

        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3);
        let expected = [(0.0, 0.25), (0.25, 0.75), (0.75, 1.0)];
        for (pb, (lo, hi)) in blocks.iter().zip(expected.iter()) {
            assert_eq!(pb.range(), (*lo, *hi));
            let n_sp = pb.edge();
            assert!(n_sp != e0, "block re-pointed to a new edge");
            let shape = f.ds.shape(n_sp).cloned().expect("split edge shape");
            assert_eq!(BRepTool::edge_parameters(&Edge(shape)), (*lo, *hi));
            assert_eq!(pb.original_edge(), e0, "original edge preserved");
        }
    }

    #[test]
    fn make_pcurves_stores_pcurve_in_registry() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.append(b.solid.0.clone()).unwrap();
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Register face 0 with edge 0 lying on it.
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        assert!(!boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some());

        make_pcurves(&mut f).unwrap();

        let pc = boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).expect("pcurve stored");
        // Face 0 (bottom, z = 0) parameterizes (u, v) → (v, u, 0): its u axis
        // runs +Y, so edge 0 (0,0,0)→(1,0,0) is a v-line at u = 0.
        let q0 = pc.d0(0.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {q0:?}");
        let q1 = pc.d0(1.0);
        assert!((q1.x() - 0.0).abs() < 1e-6 && (q1.y() - 1.0).abs() < 1e-6, "end {q1:?}");
    }

    #[test]
    fn make_pcurves_builds_pcurve_for_in_and_on_edges() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Edge 0 has no pcurve on face 0 yet.
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_none());
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave_on(e0, 0.0, 1.0);
            fi.add_pave_in(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }

        make_pcurves(&mut f).unwrap();

        // The IN/ON paves got a pcurve (deduplicated: the same edge appears in
        // both lists, the pcurve is built once).
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some(), "IN/ON pcurve built");
        assert!(f.warnings.is_empty(), "warnings: {:?}", f.warnings);
    }

    #[test]
    fn make_pcurves_builds_on_pcurve_once_when_already_present() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Edge 0 carries a pcurve on face 0 already.
        let pc = pcurve_full::make_pcurve_full(&b.edges[0], &b.faces[0]).unwrap();
        let key = GeometryRegistry::shape_key(&b.faces[0].0);
        GeometryRegistry::global().set_edge_pcurve(&b.edges[0].0, key, pc);
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave_on(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }

        make_pcurves(&mut f).unwrap();

        // The ON block whose edge already has a pcurve is skipped (OCCT bHasPC).
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some());
        assert!(f.warnings.is_empty(), "warnings: {:?}", f.warnings);
    }

    #[test]
    fn make_pcurves_updates_vertices_of_section_edge_with_existing_pcurve() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Store a deliberately wrong pcurve for edge 0 on face 0 (shifted by
        // (0.1, 0.1) in UV) — simulating a section pcurve that deviates from
        // the true surface mapping.
        let true_pc = pcurve_full::make_pcurve_full(&b.edges[0], &b.faces[0]).unwrap();
        let mut tr = occt_core::gp::GpTrsf2d::identity();
        tr.set_translation_vec(&occt_core::gp::GpVec2d::new(0.1, 0.1));
        let shifted: Arc<dyn Curve2d> = Arc::from(true_pc.transformed(&tr));
        let face_key = GeometryRegistry::shape_key(&b.faces[0].0);
        GeometryRegistry::global().set_edge_pcurve(&b.edges[0].0, face_key, shifted);
        // Register the edge as a section pave on the face.
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        let v0 = edge_vertex_shapes(&b.edges[0])[0].clone();
        let v1 = edge_vertex_shapes(&b.edges[0])[1].clone();
        let tol_before = BRepTool::vertex_tolerance(&Vertex(v0.clone()));

        make_pcurves(&mut f).unwrap();

        // The section branch calls UpdateVertices even though the pcurve already
        // exists: the boundary vertex tolerances grow to cover the 3D/2D gap.
        let tol_after = BRepTool::vertex_tolerance(&Vertex(v0.clone()));
        assert!(tol_after > tol_before, "vertex tolerance grew from {tol_before} to {tol_after}");
        let tol_after1 = BRepTool::vertex_tolerance(&Vertex(v1.clone()));
        assert!(tol_after1 > tol_before, "second vertex tolerance grew: {tol_after1}");
    }

    #[test]
    fn filter_paves_dedupes_bound_coincident_extra_pave() {
        let (mut f, e0) = box_filler();
        let i1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let i2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        // Split into three blocks.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(i1, 0.25));
            blocks[0].append_ext_pave(BopdsPave::new(i2, 0.75));
        }
        make_blocks(&mut f).unwrap();
        assert_eq!(f.ds.pave_blocks(e0).len(), 3);

        // Add a spurious extra pave at the exact upper bound of the first
        // block (a different vertex, so only the parameter collides).
        let i3 = push_vertex(&mut f, GpPnt::new(0.5, 0.5, 0.5));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(i3, 0.25));
            assert!(blocks[0].is_to_update());
        }
        filter_paves_on_curves(&mut f).unwrap();

        // The bound-coincident pave is dropped; the 3 elementary blocks remain
        // (their shared vertices are new, so no spurious merge happens).
        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3);
        for pb in blocks {
            assert!(pb.ext_paves().is_empty());
        }
    }

    #[test]
    fn filter_paves_merges_adjacent_touching_blocks() {
        let (mut f, e0) = box_filler();
        let v0 = f.ds.pave_blocks(e0)[0].indices().0;
        let v1 = f.ds.pave_blocks(e0)[0].indices().1;
        // Two adjacent blocks sharing the middle vertex/parameter.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            let mut a = BopdsPaveBlock::new();
            a.set_edge(e0);
            a.set_original_edge(e0);
            a.set_pave1(BopdsPave::new(v0, 0.0));
            a.set_pave2(BopdsPave::new(v1, 0.5));
            let mut b = BopdsPaveBlock::new();
            b.set_edge(e0);
            b.set_original_edge(e0);
            b.set_pave1(BopdsPave::new(v1, 0.5));
            b.set_pave2(BopdsPave::new(v0, 1.0));
            blocks.push(a);
            blocks.push(b);
        }
        filter_paves_on_curves(&mut f).unwrap();
        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 1, "adjacent blocks merged");
        assert_eq!(blocks[0].range(), (0.0, 1.0));
        assert_eq!(blocks[0].indices(), (v0, v0));
    }

    #[test]
    fn make_blocks_full_groups_parallel_edges_with_same_range() {
        let (mut f, inds) = box_edges(&[0, 4]);
        let (e0, e4) = (inds[0], inds[1]);
        assert_ne!(e0, e4);
        // Edge 0 (bottom front, along X) and edge 4 (top front, along X) are
        // parallel and equal-length. Artificially force the same pave-block
        // range on both — the CommonBlock detection must group them.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].set_range(0.25, 0.75);
            let blocks = f.ds_mut().change_pave_blocks_mut(e4);
            blocks[0].set_range(0.25, 0.75);
        }

        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 1, "expected one common block, got {}", cbs.len());
        let cb = &cbs[0];
        assert!(cb.contains_index(e0) && cb.contains_index(e4), "indices {:?}", cb.indices());
        assert!(cb.contains_range(0.25, 0.75, 1e-9), "ranges {:?}", cb.ranges());
    }

    #[test]
    fn make_blocks_full_does_not_group_different_ranges() {
        let (mut f, inds) = box_edges(&[0, 4]);
        let (e0, e4) = (inds[0], inds[1]);
        {
            // e4's block covers only the second half of its curve.
            let blocks = f.ds_mut().change_pave_blocks_mut(e4);
            blocks[0].set_range(0.5, 1.0);
        }

        make_blocks_full(&mut f).unwrap();

        assert!(
            f.ds.common_blocks().is_empty(),
            "no group expected, got {}",
            f.ds.common_blocks().len()
        );
        // The two edges still have their own (unmerged) block lists.
        assert_eq!(f.ds.pave_blocks(e0).len(), 1);
        assert_eq!(f.ds.pave_blocks(e4).len(), 1);
    }

    #[test]
    fn make_blocks_full_groups_by_shared_bound_vertices() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let e4 = f.ds.index(&b.edges[4].0).unwrap();
        let v4 = f.ds.index(&b.vertices[4].0).unwrap();
        let v5 = f.ds.index(&b.vertices[5].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        f.ds.init_pave_blocks_for_edge(e4);
        // Give e0's block the SAME bound vertices as e4's block, but a
        // different parametric range — the vertex match alone must group them.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].set_indices(v4, v5);
            blocks[0].set_range(0.2, 0.8);
        }

        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 1, "expected one common block, got {}", cbs.len());
        assert!(cbs[0].contains_index(e0) && cbs[0].contains_index(e4));
    }

    #[test]
    fn make_blocks_full_groups_three_edges_sharing_interval() {
        let (mut f, inds) = box_edges(&[0, 4, 6]);
        let (e0, e4, e6) = (inds[0], inds[1], inds[2]);
        // Edges 0, 4 and 6 are three parallel, equal-length box edges; their
        // default blocks all span [0, 1], so all three join one common block.
        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 1, "one common block for three coincident intervals");
        let cb = &cbs[0];
        assert!(
            cb.contains_index(e0) && cb.contains_index(e4) && cb.contains_index(e6),
            "indices {:?}",
            cb.indices()
        );
        assert!(cb.contains_range(0.0, 1.0, 1e-9));
    }

    #[test]
    fn make_blocks_full_records_two_intervals_on_same_edge() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let e4 = f.ds.index(&b.edges[4].0).unwrap();
        let v0 = f.ds.index(&b.vertices[0].0).unwrap();
        let v1 = f.ds.index(&b.vertices[1].0).unwrap();
        let v4 = f.ds.index(&b.vertices[4].0).unwrap();
        let v5 = f.ds.index(&b.vertices[5].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        f.ds.init_pave_blocks_for_edge(e4);
        // A new (intersection) vertex splits both edges at t = 0.5, so the two
        // halves are not mergeable and form two separate common intervals.
        let mid = push_vertex(&mut f, GpPnt::new(0.5, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            blocks.push(mk_block(e0, v0, 0.0, mid, 0.5));
            blocks.push(mk_block(e0, mid, 0.5, v1, 1.0));
            let blocks = f.ds_mut().change_pave_blocks_mut(e4);
            blocks.clear();
            blocks.push(mk_block(e4, v4, 0.0, mid, 0.5));
            blocks.push(mk_block(e4, mid, 0.5, v5, 1.0));
        }

        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 2, "two intervals on the same edges, got {}", cbs.len());
        for cb in cbs {
            assert!(cb.contains_index(e0) && cb.contains_index(e4), "indices {:?}", cb.indices());
        }
        let has_first = cbs.iter().any(|cb| cb.contains_range(0.0, 0.5, 1e-9));
        let has_second = cbs.iter().any(|cb| cb.contains_range(0.5, 1.0, 1e-9));
        assert!(has_first && has_second, "ranges: {:?}", cbs.iter().map(|c| c.ranges().to_vec()).collect::<Vec<_>>());
    }

    #[test]
    fn make_blocks_full_merges_touching_blocks() {
        let (mut f, e0) = box_filler();
        let (v0, v1) = {
            let blocks = f.ds.pave_blocks(e0);
            (blocks[0].indices().0, blocks[0].indices().1)
        };
        // Two clean adjacent blocks sharing source vertex v1 at t = 0.5.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            blocks.push(mk_block(e0, v0, 0.0, v1, 0.5));
            blocks.push(mk_block(e0, v1, 0.5, v0, 1.0));
        }

        make_blocks_full(&mut f).unwrap();

        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 1, "touching blocks merged into one");
        assert_eq!(blocks[0].range(), (0.0, 1.0));
    }

    #[test]
    fn make_blocks_full_does_not_merge_to_update_blocks() {
        let (mut f, e0) = box_filler();
        let (v0, v1) = {
            let blocks = f.ds.pave_blocks(e0);
            (blocks[0].indices().0, blocks[0].indices().1)
        };
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            blocks.push(mk_block(e0, v0, 0.0, v1, 0.5));
            // The second block carries an intersection split point, so it is
            // not merged away even though it touches the first one.
            let mut b = mk_block(e0, v1, 0.5, v0, 1.0);
            b.append_ext_pave(BopdsPave::new(v0, 0.75));
            blocks.push(b);
        }

        make_blocks_full(&mut f).unwrap();

        // No merge; the second block is split by update_pave_blocks.
        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[0].range(), (0.0, 0.5));
        assert_eq!(blocks[1].range(), (0.5, 0.75));
        assert_eq!(blocks[2].range(), (0.75, 1.0));
    }

    #[test]
    fn blocks_match_handles_edges_intervals_and_orientation() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let e4 = f.ds.index(&b.edges[4].0).unwrap();
        let v0 = f.ds.index(&b.vertices[0].0).unwrap();
        let v1 = f.ds.index(&b.vertices[1].0).unwrap();
        let v4 = f.ds.index(&b.vertices[4].0).unwrap();
        let v5 = f.ds.index(&b.vertices[5].0).unwrap();
        let a = mk_block(e0, v0, 0.0, v1, 1.0);
        // Blocks on the same edge never match (a self-overlap is not a common
        // block between distinct edges).
        assert!(!blocks_match(&a, &mk_block(e0, v0, 0.0, v1, 1.0), PCONFUSION), "same edge");
        // Same interval on a different edge matches.
        assert!(blocks_match(&a, &mk_block(e4, v4, 0.0, v5, 1.0), PCONFUSION), "same interval");
        // Disjoint intervals do not match.
        assert!(!blocks_match(&a, &mk_block(e4, v4, 0.5, v5, 1.0), PCONFUSION), "different interval");
        // A reversed orientation still covers the same interval.
        assert!(blocks_match(&a, &mk_block(e4, v5, 1.0, v4, 0.0), PCONFUSION), "reversed interval");
        // The same pair of bound vertices groups regardless of range.
        assert!(blocks_match(&a, &mk_block(e4, v0, 0.2, v1, 0.8), PCONFUSION), "shared bounds");
    }

    #[test]
    fn make_pcurves_full_stores_trimmed_aligned_pcurve_on_box_edge() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.append(b.solid.0.clone()).unwrap();
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_none());

        make_pcurves_full(&mut f).unwrap();

        let pc = boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).expect("pcurve stored");
        // Edge 0 runs (0,0,0)→(1,0,0) on the bottom face, whose (u,v)→(v,u,0)
        // parameterization maps the edge to the v-line u = 0, v ∈ [0, 1]. The
        // endpoints hit the face UV quad boundary.
        let q0 = pc.d0(0.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {q0:?}");
        let q1 = pc.d0(1.0);
        assert!((q1.x() - 0.0).abs() < 1e-6 && (q1.y() - 1.0).abs() < 1e-6, "end {q1:?}");
        // The endpoints equal the projection of the edge's 3D endpoints.
        let (u0, v0, _) = project_point_on_face(&b.faces[0], &GpPnt::new(0.0, 0.0, 0.0)).unwrap();
        assert!((q0.x() - u0).abs() < 1e-6 && (q0.y() - v0).abs() < 1e-6, "proj start");
        let (u1, v1, _) = project_point_on_face(&b.faces[0], &GpPnt::new(1.0, 0.0, 0.0)).unwrap();
        assert!((q1.x() - u1).abs() < 1e-6 && (q1.y() - v1).abs() < 1e-6, "proj end");
        // The trimmed pcurve lies within the face UV domain.
        let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(&b.faces[0]);
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let q = pc.d0(t);
            assert!(q.x() >= umin - 1e-6 && q.x() <= umax + 1e-6, "u {q:?} at {t}");
            assert!(q.y() >= vmin - 1e-6 && q.y() <= vmax + 1e-6, "v {q:?} at {t}");
        }
        assert!(f.warnings.is_empty(), "warnings: {:?}", f.warnings);
    }

    #[test]
    fn make_pcurves_full_covers_split_edges() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        // Split edge 0 at t = 0.5 with an intersection vertex.
        let mid = push_vertex(&mut f, GpPnt::new(0.5, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(mid, 0.5));
        }
        make_blocks(&mut f).unwrap();
        make_split_edges(&mut f).unwrap();
        let split: Vec<usize> = f.ds.pave_blocks(e0).iter().map(|pb| pb.edge()).collect();
        assert!(split.iter().any(|&e| e != e0), "split edges created: {split:?}");
        // The face info lists only the source edge.
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        make_pcurves_full(&mut f).unwrap();

        // The source edge got a pcurve.
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some());
        // Every split edge got a pcurve on the same face, with endpoints on the
        // face UV quad boundary (u = 0, v within [0, 1]).
        for &sp in &split {
            let Some(shape) = f.ds.shape(sp).cloned() else { continue };
            let e = Edge(shape);
            let pc = boptools_2d::curve_on_surface(&e, &b.faces[0]).expect("split edge pcurve");
            let (a, b) = BRepTool::edge_parameters(&e);
            let q0 = pc.d0(a);
            let q1 = pc.d0(b);
            assert!(q0.x().abs() < 1e-6, "split start {q0:?}");
            assert!(q1.x().abs() < 1e-6, "split end {q1:?}");
            assert!(q0.y() >= -1e-6 && q0.y() <= 1.0 + 1e-6, "split start y {q0:?}");
            assert!(q1.y() >= -1e-6 && q1.y() <= 1.0 + 1e-6, "split end y {q1:?}");
        }
    }

    #[test]
    fn make_pcurves_full_updates_vertex_2d_tolerance() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.append(b.solid.0.clone()).unwrap();
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Give the face a non-trivial tolerance.
        let reg = GeometryRegistry::global();
        let mut fg = reg.face_geom(&b.faces[0].0).expect("face geom");
        fg.tolerance = 0.05;
        reg.set_face(&b.faces[0].0, fg);
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        make_pcurves_full(&mut f).unwrap();

        // The boundary vertices of edge 0 carry a 2D tolerance ≥ face tolerance.
        let verts = edge_vertex_shapes(&b.edges[0]);
        for v in &verts {
            let tol = BRepTool::vertex_tolerance(&Vertex(v.clone()));
            assert!(tol >= 0.05 - 1e-12, "vertex tolerance {tol}");
        }
    }

    #[test]
    fn align_pcurve_endpoints_repairs_off_endpoints() {
        let b = unit_box();
        let edge = b.edges[0].clone();
        let face = b.faces[0].clone();
        let pc = pcurve_full::make_pcurve_full(&edge, &face).unwrap();
        // Shift the pcurve in UV so both endpoints leave the projected
        // positions by more than the alignment tolerance.
        let mut tr = occt_core::gp::GpTrsf2d::identity();
        tr.set_translation_vec(&occt_core::gp::GpVec2d::new(0.25, 0.0));
        let shifted: Arc<dyn Curve2d> = Arc::from(pc.transformed(&tr));
        let aligned = align_pcurve_endpoints(&shifted, &edge, &face, 1e-9).expect("aligned");
        // The endpoints snap back onto the projections of the edge's 3D ends.
        let q0 = aligned.d0(0.0);
        let q1 = aligned.d0(1.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {q0:?}");
        assert!((q1.x() - 0.0).abs() < 1e-6 && (q1.y() - 1.0).abs() < 1e-6, "end {q1:?}");
        // The interior is resampled from the shifted curve (still on the face).
        let qm = aligned.d0(0.5);
        assert!((qm.x() - 0.25).abs() < 1e-6, "mid {qm:?}");
    }

    #[test]
    fn align_pcurve_endpoints_leaves_trimmed_subrange_alone() {
        let b = unit_box();
        let edge = b.edges[0].clone();
        let face = b.faces[0].clone();
        let pc = pcurve_full::make_pcurve_full(&edge, &face).unwrap();
        // Simulate a trim that clipped the pcurve to a sub-range: the endpoints
        // are no longer at the edge boundary parameters, so alignment is skipped.
        let n = 17;
        let pts: Vec<GpPnt2d> = (0..n)
            .map(|i| pc.d0(0.25 + 0.5 * i as f64 / (n - 1) as f64))
            .collect();
        let bs = bspline_from_pts_2d(&pts, 0.25, 0.75, 1).unwrap();
        let trimmed: Arc<dyn Curve2d> = Arc::new(bs);
        let out = align_pcurve_endpoints(&trimmed, &edge, &face, 1e-9).expect("aligned");
        assert!(Arc::ptr_eq(&trimmed, &out), "trimmed pcurve returned unchanged");
    }
}
