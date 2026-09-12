use super::prelude::*;
use super::*;

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
pub(super) fn merge_vertices(shapes: &[TopoShape]) -> Result<TopoShape, String> {
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
pub(crate) fn make_sd_vertices(
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
            ds.add_interf_vv(n1, n2, Some(n_v));
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
pub(super) fn update_vertex_ds(ds: &mut BopdsDS, n_v: usize, tol_new: f64, non_destructive: bool) -> usize {
    // New vertices and SD vertices are updated in place; old vertices are copied
    // in non-destructive mode only.
    if ds.is_new_shape(n_v) || ds.has_shape_sd(n_v).is_some() || !non_destructive {
        if let Some(shape) = ds.shape(n_v).cloned() {
            let vtx = Vertex(shape);
            let tol = BRepTool::vertex_tolerance(&vtx);
            if tol < tol_new {
                Vertex(vtx.0.clone()).set_tolerance(tol_new);
                // Record the vertex for the repeat-intersection stage and grow
                // its DS box (OCCT `UpdateVertex` → `myIncreasedSS` +
                // `BRepBndLib::Add`/`SetGap`).
                ds.increased_ss_mut().insert(n_v);
                ds.refresh_vertex_box(n_v, tol_new);
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
            if tol < tol_new {
                ds.increased_ss_mut().insert(n_v);
                ds.refresh_vertex_box(n_new, tol_new);
            }
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
pub(crate) fn vertex_on_edge(
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
pub(super) fn vertex_on_face(
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
pub(super) fn fill_shrunk_data(ds: &BopdsDS, pb: &mut BopdsPaveBlock) {
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

/// `BOPDS_DS::IsValidShrunkData`: shrunk ends must sit outside the vertex
/// tolerance spheres (within 1% of the edge tolerance).
pub(super) fn is_valid_shrunk_data(ds: &BopdsDS, pb: &BopdsPaveBlock) -> bool {
    if !pb.has_shrunk_data() {
        return false;
    }
    let n_e = if pb.original_edge() != 0 {
        pb.original_edge()
    } else {
        pb.edge()
    };
    if n_e == 0 {
        return false;
    }
    let Some(edge_shape) = ds.shape(n_e) else { return false };
    let edge = Edge(edge_shape.clone());
    let Some(curve) = BRepTool::edge_curve(&edge) else { return false };
    let (ts1, ts2, _) = pb.shrunk_data();
    let (n_v1, n_v2) = pb.indices();
    let eps = BRepTool::edge_tolerance(&edge) * 0.01;
    for (t, n_v) in [(ts1, n_v1), (ts2, n_v2)] {
        let Some(v_shape) = ds.shape(n_v) else { return false };
        let v = Vertex(v_shape.clone());
        let tol = BRepTool::vertex_tolerance(&v) + CONFUSION;
        let dist = BRepTool::vertex_point(&v).distance(&curve.d0(t));
        if tol - dist > eps {
            return false;
        }
    }
    true
}

pub(super) fn ensure_shrunk_data(ds: &BopdsDS, pb: &mut BopdsPaveBlock) {
    if pb.has_shrunk_data() && is_valid_shrunk_data(ds, pb) {
        return;
    }
    pb.has_shrunk = false;
    fill_shrunk_data(ds, pb);
}

/// Fill the shrunk-range data of every pave block of every source edge.
pub(crate) fn fill_shrunk_data_for_all_edges(ds: &mut BopdsDS) {
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
pub(super) fn update_pb_with_sd_vertices(ds: &BopdsDS, pb: &mut BopdsPaveBlock) {
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
pub(crate) fn add_ext_pave(ds: &mut BopdsDS, n_e: usize, t: f64, n_v: usize) {
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
pub(crate) fn collect_pairs(ds: &BopdsDS, t1: ShapeType, t2: ShapeType) -> Vec<(usize, usize)> {
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
pub(super) fn connected_components(adjacency: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
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
///
/// When `pairs_override` is set (the repeat-intersection stage) the given
/// pairs are used *instead of* the iterator's regular inter-argument pairs —
/// OCCT re-initializes its iterator with the extended pair lists (the ones
/// involving the vertices with increased tolerance).
pub(super) fn perform_vv_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    pairs_override: Option<&[(usize, usize)]>,
) -> Result<(), String> {
    let pairs: Vec<(usize, usize)> = match pairs_override {
        Some(p) => p.to_vec(),
        None => collect_pairs(ds, ShapeType::Vertex, ShapeType::Vertex),
    };
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
pub(super) struct VeCandidate {
    /// Original vertex index.
    pub(super) n_v: usize,
    /// Edge index.
    pub(super) n_e: usize,
    /// Parameter of the vertex on the edge.
    pub(super) t: f64,
    /// New tolerance the vertex must reach to cover the edge point.
    pub(super) tol_v_new: f64,
}

/// Core of [`perform_ve`]: project every interfering vertex onto its edge and
/// insert paves at the projections.
///
/// `pairs_override` has the same semantics as in [`perform_vv_impl`].
pub(super) fn perform_ve_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    pairs_override: Option<&[(usize, usize)]>,
) -> Result<(), String> {
    fill_shrunk_data_for_all_edges(ds);
    let pairs: Vec<(usize, usize)> = match pairs_override {
        Some(p) => p.to_vec(),
        None => collect_pairs(ds, ShapeType::Vertex, ShapeType::Edge),
    };
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
        ds.add_interf_ve(c.n_v, c.n_e, Some(n_vx));
        modified.push(c.n_e);
    }
    split_pave_blocks_impl(ds, ctx, &modified)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Edge / Edge
// ---------------------------------------------------------------------------

/// A ready-to-run edge/edge intersection of two pave blocks.
pub(super) struct EeSolver {
    pub(super) n_e1: usize,
    pub(super) n_e2: usize,
    pub(super) pb1: BopdsPaveBlock,
    pub(super) pb2: BopdsPaveBlock,
    pub(super) edge1: Edge,
    pub(super) edge2: Edge,
    pub(super) tol1: f64,
    pub(super) tol2: f64,
}

/// The outcome of one edge/edge solver, ready for the write-back phase.
#[derive(Default)]
pub(super) struct EeOutcome {
    pub(super) n_e1: usize,
    pub(super) n_e2: usize,
    /// Seeds for the new intersection vertices.
    pub(super) seeds: Vec<NewVertexSeed>,
    /// Coincident overlap ranges `(e1_first, e1_last, e2_first, e2_last)` —
    /// each edge may parameterise the shared interval differently, so the
    /// common block must store the range per edge.
    pub(super) common_ranges: Vec<(f64, f64, f64, f64)>,
}
