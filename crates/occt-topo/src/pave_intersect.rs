//! Pairwise intersection of the PaveFiller — port of `BOPAlgo_PaveFiller_{1..4}.cxx`.
//!
//! This module contains the *intersection execution* of the Boolean component's
//! PaveFiller (Phase 19, wave C2b). The [`crate::pave_filler::PaveFiller`] owns
//! the data structure ([`BopdsDS`]) and the pipeline; the heavy per-pair
//! geometric work lives here:
//!
//! | Stage          | Source                           | Port                                                        |
//! |----------------|----------------------------------|-------------------------------------------------------------|
//! | Vertex/Vertex  | `PerformVV`                      | [`perform_vv`] — fuse coincident vertices into SD vertices  |
//! | Vertex/Edge    | `PerformVE` / `IntersectVE`      | [`perform_ve`] — project vertices onto edges, insert paves  |
//! | Edge/Edge      | `PerformEE` + `TreatNewVertices` | [`perform_ee`] — edge/edge intersection via [`EdgeEdge`]     |
//! | Vertex/Face    | `PerformVF`                      | [`perform_vf`] — classify vertices against faces            |
//! | Split          | `SplitPaveBlocks`                | [`split_pave_blocks`] — split pave blocks with extra paves  |
//!
//! The module works entirely on the DS through two handles:
//!
//! * [`BopdsDS`] — the shape registry + pave-block pool (read via
//!   [`crate::pave_filler::PaveFiller::ds`], write via `ds_mut`);
//! * [`FillCtx`] — the per-run options the intersection routines need
//!   (fuzzy value, glue mode, non-destructive flag) plus the error report.
//!
//! Every public function returns `Result<(), String>` and, on a hard failure,
//! records a fatal alert on the filler via
//! [`crate::pave_filler::PaveFiller::add_error`].

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::algo_tools::{AlgoTools, D_TOLERANCE};
use crate::bopds::{
    BopdsCommonBlock, BopdsDS, BopdsFaceInfo, BopdsIterator, BopdsPave, BopdsPaveBlock,
    BopdsShapeInfo,
};
use crate::brep_surface::surface_closest_params;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::edge_edge::EdgeEdge;
use crate::edge_face::EdgeFace;
use crate::fclass2d::FaceState;
use crate::int_face_face::FaceFace;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_data::{CommonPartType, IntRange};
use crate::pave_filler::{GlueEnum, PaveFiller};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::edges_of;

// ---------------------------------------------------------------------------
// FillCtx
// ---------------------------------------------------------------------------

/// The per-run options the intersection routines need, decoupled from the
/// [`PaveFiller`] so the core logic can be tested on a bare [`BopdsDS`].
#[derive(Debug, Clone)]
pub(crate) struct FillCtx {
    /// Additional tolerance for touching/coincidence detection, floored at
    /// `Precision::Confusion()` (the filler's `myFuzzyValue`).
    pub(crate) fuzzy: f64,
    /// Glue-full mode: coincident faces are treated as a single one without
    /// splitting them, so vertex/face intersection is skipped.
    pub(crate) glue_full: bool,
    /// Non-destructive mode: argument shapes are not modified; updated vertices
    /// are created as new shapes instead.
    pub(crate) non_destructive: bool,
    /// Fatal alerts collected during the run (the filler's `myError` list).
    pub(crate) errors: Vec<String>,
}

impl FillCtx {
    /// Build a context from the filler's options.
    pub(crate) fn from_filler(f: &PaveFiller) -> Self {
        Self {
            fuzzy: f.fuzzy_value(),
            glue_full: f.glue() == GlueEnum::Full,
            non_destructive: f.non_destructive(),
            errors: Vec::new(),
        }
    }

    /// Record a fatal alert.
    pub(crate) fn add_error(&mut self, msg: impl Into<String>) {
        self.errors.push(msg.into());
    }
}

// ---------------------------------------------------------------------------
// Small geometry / DS helpers
// ---------------------------------------------------------------------------

/// Merge a list of coincident vertex shapes into a single vertex.
///
/// Port of `BOPTools_AlgoTools::MakeVertex` + `BRepLib::BoundingVertex`: the
/// point is the average of the input points and the tolerance is the smallest
/// radius covering every input point sphere (max vertex tolerance + the
/// farthest point distance from the centre).
///
/// When `shapes` holds a single vertex it is returned unchanged.
fn merge_vertices(shapes: &[TopoShape]) -> Result<TopoShape, String> {
    if shapes.is_empty() {
        return Err("merge_vertices: empty vertex list".to_string());
    }
    if shapes.len() == 1 {
        return Ok(shapes[0].clone());
    }
    let mut cx: f64 = 0.0;
    let mut cy: f64 = 0.0;
    let mut cz: f64 = 0.0;
    let mut tol_max: f64 = 0.0;
    for s in shapes {
        let v = Vertex(s.clone());
        let p = BRepTool::vertex_point(&v);
        cx += p.x();
        cy += p.y();
        cz += p.z();
        tol_max = tol_max.max(BRepTool::vertex_tolerance(&v));
    }
    let n = shapes.len() as f64;
    let c = GpPnt::new(cx / n, cy / n, cz / n);
    let mut r: f64 = 0.0;
    for s in shapes {
        let v = Vertex(s.clone());
        r = r.max(c.distance(&BRepTool::vertex_point(&v)));
    }
    AlgoTools::make_new_vertex(&c, r + tol_max)
}

/// Make the SD (same-domain) vertex for a connected component of coincident
/// vertices.
///
/// Port of `BOPAlgo_PaveFiller::MakeSDVertices`. The given indices all point at
/// vertices that are coincident within the tolerance; a single vertex (new, or
/// the existing SD representative when the component already contains one) is
/// created and every input vertex is linked to it through
/// [`BopdsDS::add_shape_sd`]. When `add_interfs` is set, pairwise
/// [`BopdsDS::add_interf`] entries are recorded as well.
///
/// Returns the DS index of the SD vertex.
fn make_sd_vertices(
    ds: &mut BopdsDS,
    vert_indices: &[usize],
    add_interfs: bool,
) -> Result<usize, String> {
    if vert_indices.is_empty() {
        return Err("make_sd_vertices: empty vertex list".to_string());
    }
    // Collect every vertex shape, following existing SD chains so the merged
    // vertex absorbs the whole cluster.
    let mut n_sd: Option<usize> = None;
    let mut shapes: Vec<TopoShape> = Vec::new();
    for &n_x in vert_indices {
        let shape = ds
            .shape(n_x)
            .cloned()
            .ok_or_else(|| format!("make_sd_vertices: no shape at index {n_x}"))?;
        if let Some(n_sd1) = ds.has_shape_sd(n_x) {
            let sd_shape = ds
                .shape(n_sd1)
                .cloned()
                .ok_or_else(|| format!("make_sd_vertices: no SD shape at index {n_sd1}"))?;
            if n_sd.is_none() {
                n_sd = Some(n_sd1);
            } else {
                shapes.push(sd_shape);
            }
        }
        shapes.push(shape);
    }
    let v_new = merge_vertices(&shapes)?;

    // Reuse the existing SD vertex when the component already had one, updating
    // its point and tolerance in place; otherwise append a new vertex.
    let n_v = match n_sd {
        Some(n) => {
            if let Some(shape) = ds.shape(n).cloned() {
                let vtx = Vertex(shape);
                let p = BRepTool::vertex_point(&Vertex(v_new.clone()));
                let tol = BRepTool::vertex_tolerance(&Vertex(v_new.clone()));
                Vertex(vtx.0.clone()).set_point(p);
                Vertex(vtx.0.clone()).set_tolerance(tol);
            }
            n
        }
        None => ds.append(v_new)?,
    };

    // Fill the SD map and, on request, the interference table.
    for (i, &n1) in vert_indices.iter().enumerate() {
        ds.add_shape_sd(n1, n_v);
        if !add_interfs {
            continue;
        }
        for &n2 in &vert_indices[i + 1..] {
            ds.add_interf(n1, n2);
        }
    }
    Ok(n_v)
}

/// Update the tolerance of the vertex with index `n_v` so it reaches at least
/// `tol_new`, creating a new DS vertex when `n_v` is an argument sub-shape and
/// non-destructive mode is in force.
///
/// Port of `BOPAlgo_PaveFiller::UpdateVertex`. Returns the DS index of the
/// (possibly new) vertex.
fn update_vertex_ds(ds: &mut BopdsDS, n_v: usize, tol_new: f64, non_destructive: bool) -> usize {
    // New vertices and SD vertices are updated in place; old vertices are copied
    // in non-destructive mode only.
    if ds.is_new_shape(n_v) || ds.has_shape_sd(n_v).is_some() || !non_destructive {
        if let Some(shape) = ds.shape(n_v).cloned() {
            let vtx = Vertex(shape);
            let tol = BRepTool::vertex_tolerance(&vtx);
            if tol < tol_new {
                Vertex(vtx.0.clone()).set_tolerance(tol_new);
            }
        }
        return n_v;
    }
    // Old vertex in non-destructive mode: create a new vertex with the enlarged
    // tolerance and link the old one to it as an SD vertex.
    let Some(shape) = ds.shape(n_v).cloned() else { return n_v };
    let vtx = Vertex(shape);
    let p = BRepTool::vertex_point(&vtx);
    let tol = BRepTool::vertex_tolerance(&vtx);
    let v_new = TopoBuilder::new().make_vertex(p, tol.max(tol_new));
    match ds.append(v_new.into()) {
        Ok(n_new) => {
            ds.add_shape_sd(n_v, n_new);
            n_new
        }
        Err(_) => n_v,
    }
}

/// Whether the vertex `v` lies on the edge `e` within the summed tolerance.
///
/// Port of `IntTools_Context::ComputeVE`. Returns `Ok(Some((t, tol_v_new)))`
/// when the vertex projects onto the edge's 3D curve at parameter `t` and the
/// projection distance is not greater than
/// `fuzzy + edge tolerance + Precision::Confusion()`. `tol_v_new` is the vertex
/// tolerance raised to cover the projection distance. `Ok(None)` means the
/// vertex is separated from the edge (or the edge has no usable geometry).
fn vertex_on_edge(
    ctx: &IntToolsContext,
    v: &Vertex,
    e: &Edge,
    fuzzy: f64,
) -> Result<Option<(f64, f64)>, String> {
    let p = BRepTool::vertex_point(v);
    let Some(curve) = BRepTool::edge_curve(e) else {
        return Ok(None);
    };
    let Some(t) = ctx.project_point_on_edge(e, &p) else {
        return Ok(None);
    };
    let dist = p.distance(&curve.d0(t));
    let tol_sum = fuzzy.max(0.0) + BRepTool::edge_tolerance(e) + CONFUSION;
    if dist > tol_sum {
        return Ok(None);
    }
    let tol_new = BRepTool::vertex_tolerance(v).max(dist + D_TOLERANCE);
    Ok(Some((t, tol_new)))
}

/// Whether the vertex `v` lies on the face `f` within the summed tolerance.
///
/// Port of `IntTools_Context::ComputeVF`, extended to also accept vertices that
/// project onto the face *boundary* (state `On`). The vertex point is projected
/// onto the face surface; when the 3-D distance is within
/// `vertex tol + face tol + max(fuzzy, Confusion)` the UV point is classified by
/// [`FClass2d`]. `Ok(Some((u, v, tol_v_new)))` is returned for `In` and `On`
/// states, `Ok(None)` otherwise.
fn vertex_on_face(
    ctx: &mut IntToolsContext,
    v: &Vertex,
    f: &Face,
    fuzzy: f64,
) -> Result<Option<(f64, f64, f64)>, String> {
    let p = BRepTool::vertex_point(v);
    let Some(surf) = BRepTool::face_surface(f) else {
        return Ok(None);
    };
    let (u, vv) = surface_closest_params(surf.as_ref(), &p, 32, 32);
    let dist = surf.d0(u, vv).distance(&p);
    let tol_sum = BRepTool::vertex_tolerance(v) + BRepTool::face_tolerance(f) + fuzzy.max(CONFUSION);
    if dist > tol_sum {
        return Ok(None);
    }
    let classifier = ctx.fclass2d(f)?;
    let state = classifier.perform(GpPnt2d::new(u, vv));
    match state {
        FaceState::In | FaceState::On => {
            let tol_new = BRepTool::vertex_tolerance(v).max(dist + D_TOLERANCE);
            Ok(Some((u, vv, tol_new)))
        }
        FaceState::Out | FaceState::Unknown => Ok(None),
    }
}

/// Fill the shrunk-range data of a single pave block.
///
/// Port of `BOPAlgo_PaveFiller::FillShrunkData` (simplified): the bound
/// vertices are projected onto the edge's curve and their tolerance spheres
/// carve the shrunk interval `[ts1, ts2]` out of the block range. When a
/// projection fails the bound keeps the block bound. A block is *splittable*
/// when its shrunk interval is longer than `Precision::PConfusion()`.
fn fill_shrunk_data(ds: &BopdsDS, pb: &mut BopdsPaveBlock) {
    if pb.has_shrunk_data() {
        return;
    }
    let (n_v1, n_v2) = pb.indices();
    let (t1, t2) = pb.range();
    let n_e = if pb.has_edge() { pb.edge() } else { pb.original_edge() };
    let Some(edge_shape) = ds.shape(n_e) else { return };
    let e = Edge(edge_shape.clone());
    let ctx = IntToolsContext::new();
    let mut ts1 = t1;
    let mut ts2 = t2;
    if let Some(v1_shape) = ds.shape(n_v1) {
        let v1 = Vertex(v1_shape.clone());
        if let Some(tv) = ctx.project_point_on_edge(&e, &BRepTool::vertex_point(&v1)) {
            ts1 = (tv + BRepTool::vertex_tolerance(&v1)).min(t2);
        }
    }
    if let Some(v2_shape) = ds.shape(n_v2) {
        let v2 = Vertex(v2_shape.clone());
        if let Some(tv) = ctx.project_point_on_edge(&e, &BRepTool::vertex_point(&v2)) {
            ts2 = (tv - BRepTool::vertex_tolerance(&v2)).max(t1);
        }
    }
    if ts2 < ts1 {
        ts2 = ts1;
    }
    let splittable = (ts2 - ts1) > PCONFUSION;
    pb.set_shrunk_data(ts1, ts2, splittable);
}

/// Fill the shrunk-range data of every pave block of every source edge.
fn fill_shrunk_data_for_all_edges(ds: &mut BopdsDS) {
    let n = ds.nb_source_shapes();
    for i in 0..n {
        if ds.shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Edge) {
            continue;
        }
        let mut pbs = ds.pave_blocks(i).to_vec();
        if pbs.is_empty() {
            continue;
        }
        for pb in pbs.iter_mut() {
            fill_shrunk_data(ds, pb);
        }
        *ds.change_pave_blocks_mut(i) = pbs;
    }
}

/// Replace the bound indices of a pave block by their SD (same-domain)
/// representatives.
///
/// Port of `BOPDS_DS::UpdatePaveBlockWithSDVertices`.
fn update_pb_with_sd_vertices(ds: &BopdsDS, pb: &mut BopdsPaveBlock) {
    let (n1, n2) = pb.indices();
    let n1sd = ds.has_shape_sd(n1).unwrap_or(n1);
    let n2sd = ds.has_shape_sd(n2).unwrap_or(n2);
    pb.set_indices(n1sd, n2sd);
}

/// Add the vertex with index `n_v` as an extra pave at parameter `t` to the
/// pave block of the edge `n_e` that strictly contains `t`.
///
/// A no-op when no block strictly contains `t` (the point lies on a block
/// boundary or outside the edge).
fn add_ext_pave(ds: &mut BopdsDS, n_e: usize, t: f64, n_v: usize) {
    let pbs = ds.change_pave_blocks_mut(n_e);
    for pb in pbs.iter_mut() {
        let (f, l) = pb.range();
        if t > f && t < l {
            pb.append_ext_pave(BopdsPave::new(n_v, t));
            return;
        }
    }
}

/// Collect every interfering pair of shapes of the two given types from the DS,
/// using the inter-argument [`BopdsIterator`].
fn collect_pairs(ds: &BopdsDS, t1: ShapeType, t2: ShapeType) -> Vec<(usize, usize)> {
    let mut iter = BopdsIterator::new();
    iter.set_ds(ds);
    iter.prepare();
    iter.initialize(t1, t2);
    let mut out = Vec::new();
    while iter.more() {
        out.push(iter.value());
        iter.next();
    }
    out
}

/// Connected components of the undirected graph whose edges are the pairs
/// inserted by the caller (port of `BOPAlgo_Tools::FillMap` + `MakeBlocks`).
fn connected_components(adjacency: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
    let mut seen: HashSet<usize> = HashSet::new();
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for &start in adjacency.keys() {
        if !seen.insert(start) {
            continue;
        }
        let mut stack = vec![start];
        let mut comp = Vec::new();
        while let Some(u) = stack.pop() {
            comp.push(u);
            if let Some(neighbors) = adjacency.get(&u) {
                for &w in neighbors {
                    if seen.insert(w) {
                        stack.push(w);
                    }
                }
            }
        }
        blocks.push(comp);
    }
    blocks
}

// ---------------------------------------------------------------------------
// Vertex / Vertex
// ---------------------------------------------------------------------------

/// Core of [`perform_vv`]: fuse coincident vertices into same-domain vertices.
///
/// Every candidate V/V pair whose vertices coincide within the summed tolerance
/// is linked; the connected components of the link graph are then merged by
/// [`make_sd_vertices`].
fn perform_vv_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    let pairs = collect_pairs(ds, ShapeType::Vertex, ShapeType::Vertex);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for (n1, n2) in pairs {
        if ds.has_interf_pair(n1, n2) {
            adjacency.entry(n1).or_default().push(n2);
            adjacency.entry(n2).or_default().push(n1);
            continue;
        }
        let n1sd = ds.has_shape_sd(n1).unwrap_or(n1);
        let n2sd = ds.has_shape_sd(n2).unwrap_or(n2);
        let v1 = ds.shape(n1sd).cloned().ok_or("perform_vv: missing vertex 1")?;
        let v2 = ds.shape(n2sd).cloned().ok_or("perform_vv: missing vertex 2")?;
        let p2 = BRepTool::vertex_point(&Vertex(v2));
        if AlgoTools::compute_vv(&v1, &p2, ctx.fuzzy) == 1 {
            adjacency.entry(n1).or_default().push(n2);
            adjacency.entry(n2).or_default().push(n1);
        }
    }
    for block in connected_components(&adjacency) {
        make_sd_vertices(ds, &block, true)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Vertex / Edge
// ---------------------------------------------------------------------------

/// A vertex/edge intersection candidate collected during the read-only scan.
struct VeCandidate {
    /// Original vertex index.
    n_v: usize,
    /// Edge index.
    n_e: usize,
    /// Parameter of the vertex on the edge.
    t: f64,
    /// New tolerance the vertex must reach to cover the edge point.
    tol_v_new: f64,
}

/// Core of [`perform_ve`]: project every interfering vertex onto its edge and
/// insert paves at the projections.
fn perform_ve_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    fill_shrunk_data_for_all_edges(ds);
    let pairs = collect_pairs(ds, ShapeType::Vertex, ShapeType::Edge);
    if pairs.is_empty() {
        return Ok(());
    }
    let tools = IntToolsContext::new();
    let mut candidates: Vec<VeCandidate> = Vec::new();
    for (n_v, n_e) in pairs {
        let Some(si_e) = ds.shape_info(n_e) else { continue };
        if si_e.has_subshape(n_v) {
            // The vertex is a bound of the edge — no interior interference.
            continue;
        }
        if ds.has_interf_pair(n_v, n_e) {
            continue;
        }
        // The vertex already interferes with a sub-shape of the edge.
        let edge_subs: Vec<usize> = si_e.sub_shapes().to_vec();
        if edge_subs.iter().any(|&s| ds.has_interf_pair(n_v, s)) {
            continue;
        }
        let pbs = ds.pave_blocks(n_e);
        if pbs.is_empty() {
            continue;
        }
        // Micro edges (no splittable pave block) are ignored.
        if !pbs.iter().any(|pb| pb.is_splittable()) {
            continue;
        }
        let n_vsd = ds.has_shape_sd(n_v).unwrap_or(n_v);
        // The vertex is already a bound of some pave block of the edge.
        if pbs.iter().any(|pb| pb.pave1().index() == n_vsd || pb.pave2().index() == n_vsd) {
            continue;
        }
        let v = Vertex(ds.shape(n_vsd).cloned().ok_or("perform_ve: missing vertex")?);
        let e = Edge(ds.shape(n_e).cloned().ok_or("perform_ve: missing edge")?);
        if let Some((t, tol_v_new)) = vertex_on_edge(&tools, &v, &e, ctx.fuzzy)? {
            candidates.push(VeCandidate { n_v, n_e, t, tol_v_new });
        }
    }

    let mut modified: Vec<usize> = Vec::new();
    for c in candidates {
        let n_vx = update_vertex_ds(ds, c.n_v, c.tol_v_new, ctx.non_destructive);
        add_ext_pave(ds, c.n_e, c.t, n_vx);
        ds.add_interf(c.n_v, c.n_e);
        modified.push(c.n_e);
    }
    split_pave_blocks_impl(ds, ctx, &modified)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Edge / Edge
// ---------------------------------------------------------------------------

/// A ready-to-run edge/edge intersection of two pave blocks.
struct EeSolver {
    n_e1: usize,
    n_e2: usize,
    pb1: BopdsPaveBlock,
    pb2: BopdsPaveBlock,
    edge1: Edge,
    edge2: Edge,
    tol1: f64,
    tol2: f64,
}

/// The outcome of one edge/edge solver, ready for the write-back phase.
#[derive(Default)]
struct EeOutcome {
    n_e1: usize,
    n_e2: usize,
    /// Seeds for the new intersection vertices.
    seeds: Vec<NewVertexSeed>,
    /// Coincident overlap ranges `(e1_first, e1_last, e2_first, e2_last)` —
    /// each edge may parameterise the shared interval differently, so the
    /// common block must store the range per edge.
    common_ranges: Vec<(f64, f64, f64, f64)>,
}

/// Core of [`perform_ee`]: run edge/edge intersection on every interfering pair
/// of pave blocks and record the intersection vertices and common ranges.
fn perform_ee_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    fill_shrunk_data_for_all_edges(ds);
    let pairs = collect_pairs(ds, ShapeType::Edge, ShapeType::Edge);
    if pairs.is_empty() {
        return Ok(());
    }
    // --- read-only scan: build the solver list ---------------------------------
    let mut solvers: Vec<EeSolver> = Vec::new();
    for (n_e1, n_e2) in pairs {
        let edge1 = Edge(ds.shape(n_e1).cloned().ok_or("perform_ee: missing edge 1")?);
        let edge2 = Edge(ds.shape(n_e2).cloned().ok_or("perform_ee: missing edge 2")?);
        let pbs1 = ds.pave_blocks(n_e1).to_vec();
        let pbs2 = ds.pave_blocks(n_e2).to_vec();
        if pbs1.is_empty() || pbs2.is_empty() {
            continue;
        }
        let tol1 = BRepTool::edge_tolerance(&edge1);
        let tol2 = BRepTool::edge_tolerance(&edge2);
        for pb1 in &pbs1 {
            let (f1, l1) = pb1.range();
            if l1 - f1 <= PCONFUSION {
                continue;
            }
            for pb2 in &pbs2 {
                let (f2, l2) = pb2.range();
                if l2 - f2 <= PCONFUSION {
                    continue;
                }
                solvers.push(EeSolver {
                    n_e1,
                    n_e2,
                    pb1: pb1.clone(),
                    pb2: pb2.clone(),
                    edge1: edge1.clone(),
                    edge2: edge2.clone(),
                    tol1,
                    tol2,
                });
            }
        }
    }

    // --- run the solvers (no DS borrow) ---------------------------------------
    let mut outcomes: Vec<EeOutcome> = Vec::new();
    for s in &solvers {
        let mut ee = EdgeEdge::with_edges(s.edge1.clone(), s.edge2.clone());
        let (f1, l1) = s.pb1.range();
        let (f2, l2) = s.pb2.range();
        ee.set_range1(IntRange::new_unchecked(f1, l1));
        ee.set_range2(IntRange::new_unchecked(f2, l2));
        ee.set_fuzzy_value(ctx.fuzzy);
        if let Err(msg) = ee.perform() {
            ctx.add_error(format!("perform_ee: edge/edge intersection failed: {msg}"));
            continue;
        }
        let mut out = EeOutcome { n_e1: s.n_e1, n_e2: s.n_e2, ..Default::default() };
        // Discrete vertex hits: (parameter on edge 1, parameter on edge 2, point).
        for pt in ee.points() {
            out.seeds.push(NewVertexSeed {
                point: *pt.pnt1(),
                tol: ctx.fuzzy + s.tol1 + s.tol2,
                edge_a: s.n_e1,
                t_a: pt.uv1().0,
                edge_b: s.n_e2,
                t_b: pt.uv2().0,
            });
        }
        // Coincident edge overlaps: the shared interval in each edge's own
        // parameter space (the two edges may parameterise the same geometry
        // differently, e.g. two collinear segments of overlapping boxes).
        for cp in ee.common_parts() {
            if cp.part_type() == CommonPartType::Edge {
                if let Some((a, b, c, d)) = ee.coincident_ranges() {
                    if b - a > PCONFUSION {
                        out.common_ranges.push((a, b, c, d));
                    }
                }
            }
        }
        outcomes.push(out);
    }

    // --- write-back phase -----------------------------------------------------
    let mut seeds: Vec<NewVertexSeed> = Vec::new();
    for o in outcomes {
        if !o.seeds.is_empty() {
            ds.add_interf(o.n_e1, o.n_e2);
        }
        seeds.extend(o.seeds);
        for &(a, b, c, d) in &o.common_ranges {
            ds.add_interf(o.n_e1, o.n_e2);
            let mut cb = BopdsCommonBlock::new();
            cb.add_range(a, b);
            cb.add_index(o.n_e1);
            cb.add_range(c, d);
            cb.add_index(o.n_e2);
            ds.update_common_block(&cb);
        }
    }
    if !seeds.is_empty() {
        treat_new_vertices(ds, ctx.fuzzy, &seeds)?;
        // Split the pave blocks of the touched edges — OCCT's `SplitPaveBlocks`
        // at the end of `IntersectEE` (`BOPAlgo_PaveFiller_3.cxx`).
        let mut modified: Vec<usize> = Vec::new();
        for s in &seeds {
            modified.push(s.edge_a);
            if s.edge_b != s.edge_a {
                modified.push(s.edge_b);
            }
        }
        modified.sort_unstable();
        modified.dedup();
        split_pave_blocks_impl(ds, ctx, &modified)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// New vertices
// ---------------------------------------------------------------------------

/// A new vertex to create at an edge/edge (or edge/face) intersection point.
#[derive(Debug, Clone)]
pub(crate) struct NewVertexSeed {
    /// 3-D point of the intersection.
    pub point: GpPnt,
    /// Tolerance of the new vertex.
    pub tol: f64,
    /// First hit edge index and the parameter on it.
    pub edge_a: usize,
    pub t_a: f64,
    /// Second hit edge index and the parameter on it.
    pub edge_b: usize,
    pub t_b: f64,
}

/// Cluster coincident seeds by point proximity.
fn cluster_seeds(seeds: &[NewVertexSeed], fuzzy: f64) -> Vec<Vec<usize>> {
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut cluster_pts: Vec<GpPnt> = Vec::new();
    for (i, s) in seeds.iter().enumerate() {
        let tol = s.tol + fuzzy + CONFUSION;
        let tol2 = tol * tol;
        let mut placed = false;
        for (k, cpt) in cluster_pts.iter().enumerate() {
            if cpt.square_distance(&s.point) <= tol2 {
                clusters[k].push(i);
                placed = true;
                break;
            }
        }
        if !placed {
            clusters.push(vec![i]);
            cluster_pts.push(s.point);
        }
    }
    clusters
}

/// Fuse coincident new-vertex seeds, create one DS vertex per cluster and add
/// the vertex as an extra pave to every hit edge pave block.
///
/// Port of `BOPAlgo_PaveFiller::TreatNewVertices` + the append part of
/// `PerformNewVertices`. Returns the DS indices of the created vertices.
pub(crate) fn treat_new_vertices(
    ds: &mut BopdsDS,
    fuzzy: f64,
    seeds: &[NewVertexSeed],
) -> Result<Vec<usize>, String> {
    if seeds.is_empty() {
        return Ok(Vec::new());
    }
    let clusters = cluster_seeds(seeds, fuzzy);
    let mut out = Vec::new();
    for cluster in clusters {
        let mut cx: f64 = 0.0;
        let mut cy: f64 = 0.0;
        let mut cz: f64 = 0.0;
        let mut tol: f64 = 0.0;
        for &i in &cluster {
            let s = &seeds[i];
            cx += s.point.x();
            cy += s.point.y();
            cz += s.point.z();
            tol = tol.max(s.tol);
        }
        let n = cluster.len() as f64;
        let p = GpPnt::new(cx / n, cy / n, cz / n);
        // Reuse an existing source vertex at the same point (e.g. a box corner
        // of the other operand) instead of creating a coincident new one —
        // otherwise the two vertices split the same edge at the same parameter
        // into a zero-length block.
        let n_v = match find_source_vertex_at(ds, &p, tol + CONFUSION) {
            Some(existing) => existing,
            None => {
                let v = AlgoTools::make_new_vertex(&p, tol)?;
                ds.append(v)?
            }
        };
        for &i in &cluster {
            let s = &seeds[i];
            add_ext_pave(ds, s.edge_a, s.t_a, n_v);
            if s.edge_b != s.edge_a {
                add_ext_pave(ds, s.edge_b, s.t_b, n_v);
            }
        }
        out.push(n_v);
    }
    Ok(out)
}

/// The DS index of an existing source vertex within `tol` of `p`, if any.
fn find_source_vertex_at(ds: &BopdsDS, p: &GpPnt, tol: f64) -> Option<usize> {
    let n = ds.nb_source_shapes();
    let tol2 = tol * tol;
    (0..n).find(|&i| {
        ds.shape_info(i).map(|si| si.shape_type() == ShapeType::Vertex).unwrap_or(false)
            && ds.shape(i).and_then(|s| {
                if s.shape_type() == ShapeType::Vertex {
                    Some(crate::brep_tool::BRepTool::vertex_point(&Vertex(s.clone())))
                } else {
                    None
                }
            })
            .map(|q| q.square_distance(p) <= tol2)
            .unwrap_or(false)
    })
}

// ---------------------------------------------------------------------------
// Vertex / Face
// ---------------------------------------------------------------------------

/// Ensure the face with index `n_f` has a [`BopdsFaceInfo`] entry in the pool.
fn ensure_face_info(ds: &mut BopdsDS, n_f: usize) {
    let pool = ds.change_face_info_pool();
    if pool.iter().any(|fi| fi.face_index == n_f) {
        return;
    }
    pool.push(BopdsFaceInfo::new(n_f));
}

/// Record the vertex with index `n_v` (projecting to `(u, v)`) on the face `n_f`.
fn record_vertex_on_face(ds: &mut BopdsDS, n_f: usize, n_v: usize, u: f64, v: f64) {
    let pool = ds.change_face_info_pool();
    let info = match pool.iter_mut().find(|fi| fi.face_index == n_f) {
        Some(info) => info,
        None => {
            pool.push(BopdsFaceInfo::new(n_f));
            pool.last_mut().expect("just pushed")
        }
    };
    info.add_pave(n_v, u, v);
}

/// Core of [`perform_vf`]: classify every interfering vertex against its face
/// and record the vertices that lie on the face.
fn perform_vf_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    let pairs = collect_pairs(ds, ShapeType::Vertex, ShapeType::Face);
    if pairs.is_empty() {
        return Ok(());
    }
    if ctx.glue_full {
        // Glue-full mode: coincident faces are fused without splitting, so the
        // vertex/face intersection is skipped and the FaceInfo is initialized.
        for &(_, n_f) in &pairs {
            ensure_face_info(ds, n_f);
        }
        return Ok(());
    }
    let mut tools = IntToolsContext::new();
    let mut hits: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
    for (n_v, n_f) in pairs {
        let Some(si_f) = ds.shape_info(n_f) else { continue };
        if si_f.has_subshape(n_v) {
            continue;
        }
        if ds.has_interf_pair(n_v, n_f) {
            continue;
        }
        let n_vsd = ds.has_shape_sd(n_v).unwrap_or(n_v);
        let v = Vertex(ds.shape(n_vsd).cloned().ok_or("perform_vf: missing vertex")?);
        let f = Face(ds.shape(n_f).cloned().ok_or("perform_vf: missing face")?);
        if let Some((u, vv, tol_new)) = vertex_on_face(&mut tools, &v, &f, ctx.fuzzy)? {
            hits.push((n_v, n_f, u, vv, tol_new));
        }
    }
    for (n_v, n_f, u, vv, tol_new) in hits {
        let n_vx = update_vertex_ds(ds, n_v, tol_new, ctx.non_destructive);
        ds.add_interf(n_v, n_f);
        record_vertex_on_face(ds, n_f, n_vx, u, vv);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// SplitPaveBlocks
// ---------------------------------------------------------------------------

/// Core of [`split_pave_blocks`]: split the pave blocks of the given edges that
/// carry extra paves, unifying block vertices when a block has no valid range.
///
/// When `edges` is empty every source edge with update-able pave blocks is
/// processed.
fn split_pave_blocks_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    edges: &[usize],
) -> Result<(), String> {
    let to_split: Vec<usize> = if edges.is_empty() {
        (0..ds.nb_source_shapes())
            .filter(|&i| {
                ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge)
                    && ds.pave_blocks(i).iter().any(|pb| pb.is_to_update())
            })
            .collect()
    } else {
        edges.to_vec()
    };
    let mut unify_fence: HashSet<(usize, usize)> = HashSet::new();

    for n_e in to_split {
        if !ds.has_pave_blocks(n_e) {
            continue;
        }
        let old = ds.pave_blocks(n_e).to_vec();
        let mut new_list: Vec<BopdsPaveBlock> = Vec::new();
        for mut pb in old {
            if !pb.is_to_update() {
                new_list.push(pb);
                continue;
            }
            let mut out: Vec<BopdsPaveBlock> = Vec::new();
            pb.update(&mut out, true);
            for mut new_pb in out {
                // Resolve bound indices to their SD representatives and recompute
                // the shrunk range for the elementary block.
                update_pb_with_sd_vertices(ds, &mut new_pb);
                fill_shrunk_data(ds, &mut new_pb);
                let (ts1, ts2, _) = new_pb.shrunk_data();
                let has_valid = (ts2 - ts1) > PCONFUSION;
                let splittable = new_pb.is_splittable();
                let b_check_dist = has_valid && !splittable;
                if !has_valid || b_check_dist {
                    let (n_v1, n_v2) = new_pb.indices();
                    if n_v1 == n_v2 {
                        // Same vertex on both bounds — nothing to unify.
                        continue;
                    }
                    // Decide whether the vertices interfere; when they do, the
                    // block has no valid range and they must be unified.
                    let mut unify = !has_valid;
                    if b_check_dist {
                        if let (Some(a), Some(b)) = (ds.shape(n_v1).cloned(), ds.shape(n_v2).cloned()) {
                            let pb = BRepTool::vertex_point(&Vertex(b));
                            if AlgoTools::compute_vv(&a, &pb, ctx.fuzzy) == 1 {
                                unify = true;
                            }
                        }
                    }
                    if unify {
                        let key = (n_v1.min(n_v2), n_v1.max(n_v2));
                        if unify_fence.insert(key) {
                            make_sd_vertices(ds, &[n_v1, n_v2], true)?;
                        }
                        continue;
                    }
                }
                new_list.push(new_pb);
            }
        }
        *ds.change_pave_blocks_mut(n_e) = new_list;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Vertex/Vertex intersection: fuse coincident vertices of the arguments into
/// same-domain vertices.
///
/// Source: `BOPAlgo_PaveFiller::PerformVV` (`BOPAlgo_PaveFiller_1.cxx`).
pub fn perform_vv(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = perform_vv_impl(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Vertex/Edge intersection: project every interfering vertex onto its edge and
/// insert paves at the projections, splitting the affected pave blocks.
///
/// Source: `BOPAlgo_PaveFiller::PerformVE` / `IntersectVE`
/// (`BOPAlgo_PaveFiller_2.cxx`).
pub fn perform_ve(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = perform_ve_impl(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Edge/Edge intersection: intersect every interfering pair of pave blocks,
/// creating new vertices at the crossing points and recording coincident edge
/// common blocks.
///
/// Source: `BOPAlgo_PaveFiller::PerformEE` + `PerformNewVertices` /
/// `TreatNewVertices` (`BOPAlgo_PaveFiller_3.cxx`).
pub fn perform_ee(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = perform_ee_impl(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Vertex/Face intersection: classify every interfering vertex against its face
/// and record the vertices that project onto the face (inside or on the
/// boundary).
///
/// Source: `BOPAlgo_PaveFiller::PerformVF` + `TreatVerticesEE`
/// (`BOPAlgo_PaveFiller_4.cxx`).
pub fn perform_vf(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = perform_vf_impl(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Split the pave blocks of the DS that carry extra paves into elementary
/// blocks, unifying block vertices whose range collapsed.
///
/// Source: `BOPAlgo_PaveFiller::SplitPaveBlocks` (`BOPAlgo_PaveFiller_2.cxx`).
pub fn split_pave_blocks(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = split_pave_blocks_impl(f.ds_mut(), &mut ctx, &[]);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

// ---------------------------------------------------------------------------
// Edge / Face
// ---------------------------------------------------------------------------

/// An edge/face intersection hit collected during the read-only scan.
struct EfHit {
    /// Edge index.
    n_e: usize,
    /// Face index.
    n_f: usize,
    /// Seed for the new intersection vertex.
    seed: NewVertexSeed,
}

/// Mutable access to the face-info entry of the face `n_f`, creating it when
/// missing.
fn face_info_mut(ds: &mut BopdsDS, n_f: usize) -> &mut BopdsFaceInfo {
    let pool = ds.change_face_info_pool();
    if !pool.iter().any(|fi| fi.face_index == n_f) {
        pool.push(BopdsFaceInfo::new(n_f));
    }
    pool.iter_mut()
        .find(|fi| fi.face_index == n_f)
        .expect("just ensured")
}

/// Record the vertex `n_v` (whose 3-D point is `p`) on the face `n_f`, storing
/// its UV projection in the face-info pool.
fn record_vertex_point_on_face(ds: &mut BopdsDS, n_f: usize, n_v: usize, p: &GpPnt) {
    let Some(face_shape) = ds.shape(n_f).cloned() else { return };
    let face = Face(face_shape);
    let Some(surf) = BRepTool::face_surface(&face) else { return };
    let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
    record_vertex_on_face(ds, n_f, n_v, u, v);
}

/// Core of [`perform_ef`]: intersect every interfering edge with its face,
/// creating a new vertex at each piercing/touching point and recording
/// coincident sub-ranges as face paves.
fn perform_ef_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    fill_shrunk_data_for_all_edges(ds);
    let pairs = collect_pairs(ds, ShapeType::Edge, ShapeType::Face);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut hits: Vec<EfHit> = Vec::new();
    let mut face_paves: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut interfered: Vec<(usize, usize)> = Vec::new();

    for (n_e, n_f) in pairs {
        let Some(si_f) = ds.shape_info(n_f) else { continue };
        if si_f.has_subshape(n_e) {
            // The edge is a boundary of the face — no interior intersection.
            continue;
        }
        if ds.has_interf_pair(n_e, n_f) {
            continue;
        }
        let Some(edge_shape) = ds.shape(n_e).cloned() else { continue };
        let Some(face_shape) = ds.shape(n_f).cloned() else { continue };
        let edge = Edge(edge_shape);
        let face = Face(face_shape);
        if BRepTool::is_degenerated(&edge) || BRepTool::edge_curve(&edge).is_none() {
            continue;
        }
        let blocks = ds.pave_blocks(n_e);
        if blocks.is_empty() || !blocks.iter().any(|pb| (pb.range().1 - pb.range().0) > PCONFUSION) {
            continue;
        }
        let tol_e = BRepTool::edge_tolerance(&edge);
        let tol_f = BRepTool::face_tolerance(&face);

        let mut ef = EdgeFace::new();
        ef.set_edge(edge.clone());
        ef.set_face(face.clone());
        ef.set_fuzzy_value(ctx.fuzzy);
        ef.set_quick_coincidence_check(true);
        if let Err(msg) = ef.perform() {
            ctx.add_error(format!("perform_ef: edge/face intersection failed: {msg}"));
            continue;
        }
        if !ef.is_done() || ef.error_status() != 0 {
            continue;
        }
        let params = ef.point_parameters();
        if params.is_empty() {
            continue;
        }
        let curve = match BRepTool::edge_curve(&edge) {
            Some(c) => c,
            None => continue,
        };
        let tol_v = ctx.fuzzy + tol_e + tol_f;
        for (i, cp) in ef.common_parts().iter().enumerate() {
            match params.get(i).copied().flatten() {
                Some(t) => {
                    // The edge pierces/touches the face at a single parameter.
                    let p = curve.d0(t);
                    hits.push(EfHit {
                        n_e,
                        n_f,
                        seed: NewVertexSeed {
                            point: p,
                            tol: tol_v,
                            edge_a: n_e,
                            t_a: t,
                            edge_b: n_e,
                            t_b: t,
                        },
                    });
                    interfered.push((n_e, n_f));
                }
                None => {
                    // Coincident sub-range: the edge lies on the face over [f, l].
                    let r = cp.range();
                    if r.length() > PCONFUSION {
                        face_paves.push((n_e, n_f, r.first, r.last));
                        interfered.push((n_e, n_f));
                    }
                }
            }
        }
    }

    // Write back.
    for (n_e, n_f) in interfered {
        ds.add_interf(n_e, n_f);
    }
    for (n_e, n_f, t1, t2) in face_paves {
        let fi = face_info_mut(ds, n_f);
        // A coincident edge lying on the face is an IN pave block: it does not
        // cut the face (OCCT `FaceInfoIn`), so it is stored separately from the
        // section paves.
        if !fi
            .paves_in
            .iter()
            .any(|&(e, f, l)| e == n_e && (f - t1).abs() <= PCONFUSION && (l - t2).abs() <= PCONFUSION)
        {
            fi.add_pave_in(n_e, t1, t2);
        }
    }
    if !hits.is_empty() {
        let seeds: Vec<NewVertexSeed> = hits.iter().map(|h| h.seed.clone()).collect();
        let clusters = cluster_seeds(&seeds, ctx.fuzzy);
        let created = treat_new_vertices(ds, ctx.fuzzy, &seeds)?;
        for (k, cluster) in clusters.iter().enumerate() {
            let n_v = created[k];
            for &si in cluster {
                let h = &hits[si];
                record_vertex_point_on_face(ds, h.n_f, n_v, &h.seed.point);
            }
        }
        // Split the pave blocks of the edges pierced by the face — OCCT's
        // `SplitPaveBlocks` invoked through `PerformNewVertices` in
        // `IntersectEF` (`BOPAlgo_PaveFiller_5.cxx` → `_3.cxx`).
        let mut modified: Vec<usize> = Vec::new();
        for h in &hits {
            modified.push(h.n_e);
        }
        modified.sort_unstable();
        modified.dedup();
        split_pave_blocks_impl(ds, ctx, &modified)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Face / Face
// ---------------------------------------------------------------------------

/// A face/face intersection section-edge task, ready for the write-back phase.
struct FfTask {
    /// DS index of the sorted first face.
    n_fa: usize,
    /// DS index of the sorted second face.
    n_fb: usize,
    /// Sorted first face.
    face_a: Face,
    /// Sorted second face.
    face_b: Face,
    /// The 3-D intersection curve.
    curve: Arc<dyn Curve>,
    /// Parameter range of the curve valid on both faces.
    range: IntRange,
    /// Pcurve of the curve on face `a`.
    pcurve_a: Option<Arc<dyn Curve2d>>,
    /// Pcurve of the curve on face `b`.
    pcurve_b: Option<Arc<dyn Curve2d>>,
}

/// Core of [`perform_ff`]: intersect every interfering pair of faces and build
/// a section edge for each intersection curve.
fn perform_ff_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    let pairs = collect_pairs(ds, ShapeType::Face, ShapeType::Face);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut tasks: Vec<FfTask> = Vec::new();
    for (n_f1, n_f2) in pairs {
        if ds.has_interf_pair(n_f1, n_f2) {
            continue;
        }
        let Some(f1_shape) = ds.shape(n_f1).cloned() else { continue };
        let Some(f2_shape) = ds.shape(n_f2).cloned() else { continue };
        let face1 = Face(f1_shape);
        let face2 = Face(f2_shape);
        let tol =
            ctx.fuzzy + BRepTool::face_tolerance(&face1) + BRepTool::face_tolerance(&face2);
        let mut ff = FaceFace::new();
        ff.set_face1(face1);
        ff.set_face2(face2);
        ff.set_tolerance(tol);
        if let Err(msg) = ff.perform() {
            ctx.add_error(format!("perform_ff: face/face intersection failed: {msg}"));
            continue;
        }
        if !ff.is_done() {
            continue;
        }
        // The faces are internally sorted (higher analytic type first).
        let Some(fa) = ff.face1().cloned() else { continue };
        let Some(fb) = ff.face2().cloned() else { continue };
        let Some(n_fa) = ds.index(&fa.0) else { continue };
        let Some(n_fb) = ds.index(&fb.0) else { continue };
        let res = ff.result();
        for c in res.curves() {
            let r = c.range;
            if !r.is_valid() || r.length() <= PCONFUSION {
                continue;
            }
            tasks.push(FfTask {
                n_fa,
                n_fb,
                face_a: fa.clone(),
                face_b: fb.clone(),
                curve: c.curve.clone(),
                range: r,
                pcurve_a: c.pcurve1.clone(),
                pcurve_b: c.pcurve2.clone(),
            });
        }
    }
    if tasks.is_empty() {
        return Ok(());
    }

    let tol_edge = ctx.fuzzy + D_TOLERANCE;
    for t in tasks {
        let tol =
            tol_edge + BRepTool::face_tolerance(&t.face_a) + BRepTool::face_tolerance(&t.face_b);
        let p1 = t.curve.d0(t.range.first);
        let p2 = t.curve.d0(t.range.last);
        let v1 = match AlgoTools::make_new_vertex(&p1, tol) {
            Ok(v) => v,
            Err(msg) => {
                ctx.add_error(format!("perform_ff: {msg}"));
                continue;
            }
        };
        let v2 = match AlgoTools::make_new_vertex(&p2, tol) {
            Ok(v) => v,
            Err(msg) => {
                ctx.add_error(format!("perform_ff: {msg}"));
                continue;
            }
        };
        let n_v1 = match ds.append(v1) {
            Ok(n) => n,
            Err(msg) => {
                ctx.add_error(format!("perform_ff: {msg}"));
                continue;
            }
        };
        let n_v2 = match ds.append(v2) {
            Ok(n) => n,
            Err(msg) => {
                ctx.add_error(format!("perform_ff: {msg}"));
                continue;
            }
        };
        let e_shape1 = ds.shape(n_v1).cloned().ok_or_else(|| {
            format!("perform_ff: vertex {n_v1} not found after append")
        });
        let e_shape2 = ds.shape(n_v2).cloned().ok_or_else(|| {
            format!("perform_ff: vertex {n_v2} not found after append")
        });
        let (v1s, v2s) = match (e_shape1, e_shape2) {
            (Ok(a), Ok(b)) => (a, b),
            (Err(msg), _) | (_, Err(msg)) => {
                ctx.add_error(msg);
                continue;
            }
        };
        let edge = match AlgoTools::make_edge(
            t.curve.clone(),
            Some(&v1s),
            t.range.first,
            Some(&v2s),
            t.range.last,
            tol,
        ) {
            Ok(e) => e,
            Err(msg) => {
                ctx.add_error(format!("perform_ff: {msg}"));
                continue;
            }
        };
        // Attach the p-curves computed by the FaceFace solver.
        let reg = GeometryRegistry::global();
        if let Some(pc) = t.pcurve_a {
            let key = GeometryRegistry::shape_key(&t.face_a.0);
            reg.set_edge_pcurve(&edge, key, pc);
        }
        if let Some(pc) = t.pcurve_b {
            let key = GeometryRegistry::shape_key(&t.face_b.0);
            reg.set_edge_pcurve(&edge, key, pc);
        }
        let mut si = BopdsShapeInfo::new(edge);
        si.change_sub_shapes().extend_from_slice(&[n_v1, n_v2]);
        let n_edge = ds.append_info(si);
        // Record the section edge on both faces and the F/F interference.
        let fi_a = face_info_mut(ds, t.n_fa);
        fi_a.add_pave(n_edge, t.range.first, t.range.last);
        let fi_b = face_info_mut(ds, t.n_fb);
        fi_b.add_pave(n_edge, t.range.first, t.range.last);
        ds.add_interf(t.n_fa, t.n_fb);

        // Split the boundary edges of both faces at the section vertices so the
        // on-face edge endpoints connect to the face boundary — the
        // `UpdatePaveBlocks` step of `BOPAlgo_PaveFiller::PerformFF`. Without
        // it the boundary edges stay whole, the face-image WireSplitter sees an
        // open chain, and the split faces never close.
        let section_verts = [n_v1, n_v2];
        let mut modified: Vec<usize> = Vec::new();
        let tools = IntToolsContext::new();
        for n_f in [t.n_fa, t.n_fb] {
            let Some(f_shape) = ds.shape(n_f).cloned() else { continue };
            let boundary: Vec<usize> =
                edges_of(&Face(f_shape).0).iter().filter_map(|e| ds.index(&e.0)).collect();
            for n_e in boundary {
                for &n_v in &section_verts {
                    let Some(si_e) = ds.shape_info(n_e) else { continue };
                    if si_e.has_subshape(n_v) {
                        continue; // the vertex is a bound of the edge
                    }
                    let Some(v_shape) = ds.shape(n_v).cloned() else { continue };
                    let Some(e_shape) = ds.shape(n_e).cloned() else { continue };
                    if let Ok(Some((tt, _))) =
                        vertex_on_edge(&tools, &Vertex(v_shape), &Edge(e_shape), tol)
                    {
                        add_ext_pave(ds, n_e, tt, n_v);
                        modified.push(n_e);
                    }
                }
            }
        }
        if !modified.is_empty() {
            split_pave_blocks_impl(ds, ctx, &modified)?;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Edge/Face intersection: intersect every interfering edge with its face,
/// creating new vertices at the piercing points and recording coincident
/// sub-ranges as face paves.
///
/// Source: `BOPAlgo_PaveFiller::PerformEF` + `IntersectEF`
/// (`BOPAlgo_PaveFiller_5.cxx`).
pub fn perform_ef(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = perform_ef_impl(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Face/Face intersection: intersect every interfering pair of faces and build
/// a section edge for each intersection curve, recording the edge on both
/// faces and the F/F interference.
///
/// Source: `BOPAlgo_PaveFiller::PerformFF` (`BOPAlgo_PaveFiller_6.cxx`).
pub fn perform_ff(f: &mut PaveFiller) -> Result<(), String> {
    let mut ctx = FillCtx::from_filler(f);
    let r = perform_ff_impl(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use occt_core::gp::{GpAx3, GpDir, GpLin, GpVec};

    /// A unit box with its corner at `(ox, oy, oz)` (fresh shapes, so the DS
    /// sees two distinct arguments even when their corners coincide).
    struct BoxAt {
        solid: crate::shape::Solid,
        vertices: Vec<Vertex>,
        #[allow(dead_code)] // edge/face lists kept for future fixture use
        edges: Vec<Edge>,
        #[allow(dead_code)]
        faces: Vec<Face>,
    }

    /// Build a unit box with the corner at the given offset.
    fn box_at(ox: f64, oy: f64, oz: f64) -> BoxAt {
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(ox, oy, oz),
            GpPnt::new(1.0 + ox, oy, oz),
            GpPnt::new(1.0 + ox, 1.0 + oy, oz),
            GpPnt::new(ox, 1.0 + oy, oz),
            GpPnt::new(ox, oy, 1.0 + oz),
            GpPnt::new(1.0 + ox, oy, 1.0 + oz),
            GpPnt::new(1.0 + ox, 1.0 + oy, 1.0 + oz),
            GpPnt::new(ox, 1.0 + oy, 1.0 + oz),
        ];
        let vertices: Vec<Vertex> = corners.iter().map(|p| b.make_vertex(*p, 1e-7)).collect();
        let edge_idx: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges = Vec::new();
        for &(i, j) in &edge_idx {
            let p1 = corners[i];
            let p2 = corners[j];
            let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).expect("distinct");
            let lin = GpLin::from_pnt_dir(p1, dir);
            let mut e = b.make_edge(
                std::sync::Arc::new(occt_geom::GeomLine::new(lin)),
                0.0,
                p1.distance(&p2),
            );
            b.add(&mut e.0, &vertices[i].0);
            b.add(&mut e.0, &vertices[j].0);
            edges.push(e);
        }
        // Face planes: the same normal/u axes as the unit box, origins shifted.
        let face_planes: [(GpPnt, GpDir, GpDir); 6] = [
            (GpPnt::new(ox, oy, oz), GpDir::new(0.0, 0.0, -1.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
            (GpPnt::new(ox, oy, 1.0 + oz), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(ox, oy, oz), GpDir::new(0.0, -1.0, 0.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(ox, 1.0 + oy, oz), GpDir::new(0.0, 1.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(ox, oy, oz), GpDir::new(-1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(1.0 + ox, oy, oz), GpDir::new(1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
        ];
        let face_edge_sets: [[usize; 4]; 6] = [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 9, 4, 8],
            [2, 10, 6, 11],
            [3, 11, 7, 8],
            [1, 10, 5, 9],
        ];
        let mut faces = Vec::new();
        for fi in 0..6 {
            let (origin, normal, u_dir) = face_planes[fi];
            let ax3 = GpAx3::new(origin, normal, &u_dir).expect("perpendicular axes");
            let mut face = b.make_face_plane(&occt_core::gp::GpPln::new(ax3));
            let wire = b.make_wire(&face_edge_sets[fi].map(|ei| edges[ei].clone()));
            b.add_wire(&mut face, &wire);
            faces.push(face);
        }
        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        BoxAt { solid, vertices, edges, faces }
    }

    /// Initialize the pave blocks of every source edge of `ds`.
    fn init_all_edge_blocks(ds: &mut BopdsDS) {
        let n = ds.nb_source_shapes();
        for i in 0..n {
            if ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge) {
                ds.init_pave_blocks_for_edge(i);
            }
        }
    }

    #[test]
    fn merge_vertices_averages_points_and_takes_max_tolerance() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 1e-7);
        let v2 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 2e-7);
        let merged = merge_vertices(&[v1.0.clone(), v2.0.clone()]).unwrap();
        let mv = Vertex(merged);
        let p = BRepTool::vertex_point(&mv);
        assert!((p.x().abs() < 1e-9) && (p.y().abs() < 1e-9) && (p.z().abs() < 1e-9));
        assert!(BRepTool::vertex_tolerance(&mv) >= 2e-7);
    }

    #[test]
    fn make_sd_vertices_links_component_and_interferes() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 1e-7);
        let v2 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 1e-7);
        let mut ds = BopdsDS::new();
        let n1 = ds.append(v1.into()).unwrap();
        let n2 = ds.append(v2.into()).unwrap();
        let nv = make_sd_vertices(&mut ds, &[n1, n2], true).unwrap();
        assert_eq!(ds.get_same_domain_index(n1), nv);
        assert_eq!(ds.get_same_domain_index(n2), nv);
        assert!(ds.has_interf_pair(n1, n2));
    }

    #[test]
    fn perform_vv_merges_shared_box_corners() {
        let a = unit_box();
        let bbox = box_at(1.0, 0.0, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        // The 4 corners of the shared face (x = 1) coincide.
        let idx_a = [
            a.vertices[1].0.clone(), // (1,0,0)
            a.vertices[2].0.clone(), // (1,1,0)
            a.vertices[5].0.clone(), // (1,0,1)
            a.vertices[6].0.clone(), // (1,1,1)
        ];
        let idx_b = [
            bbox.vertices[0].0.clone(), // (1,0,0)
            bbox.vertices[3].0.clone(), // (1,1,0)
            bbox.vertices[4].0.clone(), // (1,0,1)
            bbox.vertices[7].0.clone(), // (1,1,1)
        ];
        perform_vv(&mut f).unwrap();
        for i in 0..4 {
            let na = f.ds().index(&idx_a[i]).unwrap();
            let nb = f.ds().index(&idx_b[i]).unwrap();
            assert!(
                f.ds().has_shape_sd(na).is_some() || f.ds().has_shape_sd(nb).is_some(),
                "corner pair {i} merged into an SD vertex"
            );
            assert_eq!(
                f.ds().get_same_domain_index(na),
                f.ds().get_same_domain_index(nb),
                "corner pair {i} share the same SD vertex"
            );
        }
        // A non-coincident corner pair is not merged.
        let n0a = f.ds().index(&a.vertices[0].0).unwrap(); // (0,0,0)
        let n0b = f.ds().index(&bbox.vertices[1].0).unwrap(); // (2,0,0)
        assert_ne!(f.ds().get_same_domain_index(n0a), f.ds().get_same_domain_index(n0b));
    }

    #[test]
    fn vertex_on_edge_finds_interior_parameter() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let v = b.make_vertex(GpPnt::new(0.5, 0.0, 0.0), 1e-7);
        let ctx = IntToolsContext::new();
        let (t, tol) = vertex_on_edge(&ctx, &v, &e, 1e-7).unwrap().unwrap();
        assert!((t - 0.5).abs() < 1e-6, "t = {t}");
        assert!(tol >= 1e-7);
        // A vertex separated from the edge is not found.
        let far = b.make_vertex(GpPnt::new(0.5, 5.0, 0.0), 1e-7);
        assert!(vertex_on_edge(&ctx, &far, &e, 1e-7).unwrap().is_none());
    }

    #[test]
    fn perform_ve_paves_vertex_on_edge() {
        // A = [0,1]^3, B shifted by (0.5, 0, 0). B's corner (0.5, 0, 0) lies on
        // A's bottom-front edge (0,0,0)-(1,0,0), so perform_ve must insert a
        // pave at t = 0.5 and split A's edge into two blocks.
        let a = unit_box();
        let bbox = box_at(0.5, 0.0, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        init_all_edge_blocks(f.ds_mut());
        perform_ve(&mut f).unwrap();

        let na_e = f.ds().index(&a.edges[0].0).unwrap();
        let nb_v = f.ds().index(&bbox.vertices[0].0).unwrap(); // (0.5, 0, 0)
        // The edge was split into two blocks by the interior pave.
        let blocks = f.ds().pave_blocks(na_e);
        assert_eq!(blocks.len(), 2, "edge split at the vertex projection");
        // The interior block references B's corner vertex.
        let mid = blocks
            .iter()
            .find(|pb| pb.pave1().index() == nb_v || pb.pave2().index() == nb_v)
            .expect("a block is bounded by B's corner vertex");
        let (t1, t2) = mid.range();
        assert!(t1.abs() < 1e-6 || (t2 - 1.0).abs() < 1e-6, "range {t1}..{t2} touches the split");
        // The V/E interference was recorded.
        assert!(f.ds().has_interf_pair(nb_v, na_e));
    }

    #[test]
    fn perform_ee_creates_vertex_at_crossing() {
        // A = [0,1]^3, B shifted by (0.5, 0.5, 0). A's right-bottom edge
        // (1,0,0)-(1,1,0) crosses B's bottom-front edge (0.5,0.5,0)-(1.5,0.5,0)
        // at (1, 0.5, 0), interior to both.
        let a = unit_box();
        let bbox = box_at(0.5, 0.5, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        init_all_edge_blocks(f.ds_mut());
        perform_ee(&mut f).unwrap();

        // A new vertex close to (1, 0.5, 0) must be present among the new
        // shapes appended after the source shapes.
        let n_src = f.ds().nb_source_shapes();
        let mut found = false;
        for i in n_src..f.ds().nb_shapes() {
            if f.ds().shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Vertex) {
                continue;
            }
            let p = BRepTool::vertex_point(&Vertex(f.ds().shape(i).unwrap().clone()));
            if p.distance(&GpPnt::new(1.0, 0.5, 0.0)) < 1e-6 {
                found = true;
                break;
            }
        }
        assert!(found, "crossing vertex created in the DS");

        // The crossing edge of A is split at the crossing: the block carrying
        // the extra pave was replaced by elementary blocks (OCCT
        // `SplitPaveBlocks` at the end of `IntersectEE`).
        let na_e = f.ds().index(&a.edges[1].0).unwrap();
        let blocks = f.ds().pave_blocks(na_e);
        assert!(blocks.len() >= 2, "A's right-bottom edge is split at the crossing");
        assert!(
            blocks.iter().all(|pb| pb.ext_paves().is_empty()),
            "all extra paves consumed into elementary blocks"
        );
    }

    #[test]
    fn split_pave_blocks_splits_block_with_two_extra_paves() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let v1 = b.make_vertex(GpPnt::new(1.0, 0.0, 0.0), 1e-7);
        let v2 = b.make_vertex(GpPnt::new(2.0, 0.0, 0.0), 1e-7);
        let mut ds = BopdsDS::new();
        let n_e = ds.append(e.0.clone()).unwrap();
        let n1 = ds.append(v1.into()).unwrap();
        let n2 = ds.append(v2.into()).unwrap();
        ds.init_pave_blocks_for_edge(n_e);
        // Add two extra paves inside the block.
        {
            let pbs = ds.change_pave_blocks_mut(n_e);
            pbs[0].append_ext_pave(BopdsPave::new(n1, 1.0));
            pbs[0].append_ext_pave(BopdsPave::new(n2, 2.0));
        }
        assert_eq!(ds.pave_blocks(n_e).len(), 1);
        assert!(ds.pave_blocks(n_e)[0].is_to_update());

        let mut ctx = FillCtx { fuzzy: 1e-7, glue_full: false, non_destructive: false, errors: Vec::new() };
        split_pave_blocks_impl(&mut ds, &mut ctx, &[n_e]).unwrap();
        let blocks = ds.pave_blocks(n_e);
        assert_eq!(blocks.len(), 3, "2 extra paves split 1 block into 3");
        let ranges: Vec<(f64, f64)> = blocks.iter().map(|pb| pb.range()).collect();
        assert_eq!(ranges[0], (0.0, 1.0));
        assert_eq!(ranges[1], (1.0, 2.0));
        assert_eq!(ranges[2], (2.0, 3.0));
    }

    #[test]
    fn perform_vf_classifies_inside_on_and_out() {
        let box_solid = unit_box();
        let b = TopoBuilder::new();
        // Vertices exactly on the bottom-face plane (z = 0) so the bounding-box
        // candidate filter (which has no tolerance gap) still pairs them.
        let inside = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 1e-7);
        let on_edge = b.make_vertex(GpPnt::new(0.5, 0.0, 0.0), 1e-7);
        let mid_box = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 1e-7);
        // Package the three vertices as a second argument (compound).
        let comp = b.make_compound_of(&[inside.0.clone(), on_edge.0.clone(), mid_box.0.clone()]);

        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone(), comp.0.clone()]);
        f.init().unwrap();
        perform_vf(&mut f).unwrap();

        // The bottom face of the box is face 0 (z = 0).
        let n_face = f.ds().index(&box_solid.faces[0].0).unwrap();
        let pool = f.ds().face_info_pool();
        let info = pool.iter().find(|fi| fi.face_index == n_face);
        assert!(info.is_some(), "face info initialized");
        let info = info.unwrap();
        let n_inside = f.ds().index(&inside.0).unwrap();
        let n_on = f.ds().index(&on_edge.0).unwrap();
        let n_mid = f.ds().index(&mid_box.0).unwrap();
        let recorded: Vec<usize> = info.paves().iter().map(|p| p.0).collect();
        assert!(recorded.contains(&n_inside), "inside vertex recorded");
        assert!(recorded.contains(&n_on), "on-boundary vertex recorded");
        assert!(!recorded.contains(&n_mid), "floating vertex not recorded");
        // The interference table records the inside/on pairs.
        assert!(f.ds().has_interf_pair(n_inside, n_face));
        assert!(f.ds().has_interf_pair(n_on, n_face));
    }

    #[test]
    fn perform_vf_glue_full_only_initializes_face_info() {
        let box_solid = unit_box();
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 1e-7);
        let comp = b.make_compound_of(&[v.0.clone()]);
        let mut f = PaveFiller::new();
        f.set_glue(GlueEnum::Full);
        f.set_arguments(&[box_solid.solid.0.clone(), comp.0.clone()]);
        f.init().unwrap();
        perform_vf(&mut f).unwrap();
        let n_face = f.ds().index(&box_solid.faces[0].0).unwrap();
        let pool = f.ds().face_info_pool();
        assert!(
            pool.iter().any(|fi| fi.face_index == n_face),
            "face info initialized in glue-full mode"
        );
        let n_v = f.ds().index(&v.0).unwrap();
        assert!(!f.ds().has_interf_pair(n_v, n_face));
    }

    #[test]
    fn treat_new_vertices_fuses_coincident_seeds() {
        let mut ds = BopdsDS::new();
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let ne1 = ds.append(e1.0.clone()).unwrap();
        let ne2 = ds.append(e2.0.clone()).unwrap();
        ds.init_pave_blocks_for_edge(ne1);
        ds.init_pave_blocks_for_edge(ne2);
        let seeds = vec![
            NewVertexSeed {
                point: GpPnt::new(0.5, 0.0, 0.0),
                tol: 1e-7,
                edge_a: ne1,
                t_a: 0.5,
                edge_b: ne2,
                t_b: 0.5,
            },
            NewVertexSeed {
                point: GpPnt::new(0.5, 0.0, 0.0),
                tol: 1e-7,
                edge_a: ne1,
                t_a: 0.5,
                edge_b: ne2,
                t_b: 0.5,
            },
        ];
        let idx = treat_new_vertices(&mut ds, 1e-7, &seeds).unwrap();
        assert_eq!(idx.len(), 1, "coincident seeds fuse into one vertex");
        assert_eq!(ds.pave_blocks(ne1)[0].ext_paves().len(), 1, "single extra pave on edge 1");
        assert_eq!(ds.pave_blocks(ne2)[0].ext_paves().len(), 1, "single extra pave on edge 2");
    }

    // -----------------------------------------------------------------------
    // perform_ef
    // -----------------------------------------------------------------------

    #[test]
    fn perform_ef_creates_vertex_at_piercing_point() {
        // A unit box + an edge piercing its bottom face (z = 0) at (0.5, 0.5, 0).
        let box_solid = unit_box();
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        let comp = b.make_compound_of(&[e.0.clone()]);

        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone(), comp.0.clone()]);
        f.init().unwrap();
        init_all_edge_blocks(f.ds_mut());
        perform_ef(&mut f).unwrap();
        assert!(!f.has_errors(), "errors: {:?}", f.errors());

        // A new vertex close to (0.5, 0.5, 0) must exist among the new shapes.
        let n_src = f.ds().nb_source_shapes();
        let mut found = false;
        for i in n_src..f.ds().nb_shapes() {
            if f.ds().shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Vertex) {
                continue;
            }
            let p = BRepTool::vertex_point(&Vertex(f.ds().shape(i).unwrap().clone()));
            if p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-4 {
                found = true;
                break;
            }
        }
        assert!(found, "piercing vertex created in the DS");

        // The E/F interference is recorded and the edge is split at the
        // piercing point (the extra pave was consumed into elementary blocks).
        let n_e = f.ds().index(&e.0).unwrap();
        let n_f = f.ds().index(&box_solid.faces[0].0).unwrap();
        assert!(f.ds().has_interf_pair(n_e, n_f));
        let blocks = f.ds().pave_blocks(n_e);
        assert!(
            blocks.iter().any(|pb| pb.pave1().parameter() == 1.0 || pb.pave2().parameter() == 1.0),
            "piercing point is a bound of an elementary block"
        );
        assert!(
            blocks.iter().all(|pb| pb.ext_paves().is_empty()),
            "all extra paves consumed into elementary blocks"
        );
    }

    #[test]
    fn perform_ef_ignores_face_boundary_edge() {
        // Edge 0 of the box IS a boundary of face 0 — no new vertex must be
        // created for the pair (the EF pass skips boundary sub-shapes).
        let box_solid = unit_box();
        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone()]);
        f.init().unwrap();
        let n_before = f.ds().nb_shapes();
        init_all_edge_blocks(f.ds_mut());
        perform_ef(&mut f).unwrap();
        assert!(!f.has_errors());
        // One argument → no inter-argument pairs; the DS is untouched.
        assert_eq!(f.ds().nb_shapes(), n_before);
    }

    // -----------------------------------------------------------------------
    // perform_ff
    // -----------------------------------------------------------------------

    #[test]
    fn perform_ff_creates_section_edge_for_intersecting_faces() {
        // Box A = [0,1]^3, box B shifted by (0.5, 0.5, 0). A's bottom face
        // (z = 0) and B's front face (y = 0.5) meet in the segment
        // (0.5, 0.5, 0)-(1, 0.5, 0).
        let a = unit_box();
        let bbox = box_at(0.5, 0.5, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        perform_ff(&mut f).unwrap();
        assert!(!f.has_errors(), "errors: {:?}", f.errors());

        // A new edge whose mid-point is (0.75, 0.5, 0) must exist.
        let n_src = f.ds().nb_source_shapes();
        let mut found = false;
        for i in n_src..f.ds().nb_shapes() {
            if f.ds().shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Edge) {
                continue;
            }
            let e = Edge(f.ds().shape(i).unwrap().clone());
            let Some(curve) = BRepTool::edge_curve(&e) else { continue };
            let (t1, t2) = BRepTool::edge_parameters(&e);
            if !t1.is_finite() || !t2.is_finite() || t2 <= t1 {
                continue;
            }
            let mid = curve.d0(0.5 * (t1 + t2));
            if mid.distance(&GpPnt::new(0.75, 0.5, 0.0)) < 1e-4 {
                found = true;
                break;
            }
        }
        assert!(found, "section edge created in the DS");

        // The F/F interference is recorded for the intersecting faces.
        let n_fa = f.ds().index(&a.faces[0].0).unwrap(); // A bottom (z = 0)
        let n_fb = f.ds().index(&bbox.faces[2].0).unwrap(); // B front (y = 0.5)
        assert!(f.ds().has_interf_pair(n_fa, n_fb));
    }

    #[test]
    fn perform_ff_disjoint_faces_produce_no_edges() {
        // A single box: no inter-argument face pairs, so no section edges.
        let box_solid = unit_box();
        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone()]);
        f.init().unwrap();
        let n_before = f.ds().nb_shapes();
        perform_ff(&mut f).unwrap();
        assert!(!f.has_errors());
        assert_eq!(f.ds().nb_shapes(), n_before);
    }
}
