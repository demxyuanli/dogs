use super::prelude::*;

// ---------------------------------------------------------------------------
// FillCtx (test-only)
// ---------------------------------------------------------------------------

/// The per-run options the DS-only cores need, decoupled from the
/// [`PaveFiller`] so the core logic can be tested on a bare [`BopdsDS`].
///
/// The live pipeline reads the fuzzy value straight off the filler, so this
/// bundle only carries what the isolated tests drive.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct FillCtx {
    /// Additional tolerance for touching/coincidence detection, floored at
    /// `Precision::Confusion()` (the filler's `myFuzzyValue`).
    pub(crate) fuzzy: f64,
}

// ---------------------------------------------------------------------------
// Small geometry / DS helpers
// ---------------------------------------------------------------------------

/// Merge a list of coincident vertex shapes into a single vertex.
///
/// Port of `BOPTools_AlgoTools::MakeVertex` + `BRepLib::BoundingVertex`: the
/// point is the average of the input points and the tolerance is the smallest
/// radius covering every input point sphere (max vertex tolerance + the
/// farthest point distance from the center).
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
pub(crate) fn intersect_make_sd_vertices(
    ds: &mut BopdsDS,
    vert_indices: &[usize],
    add_interfs: bool,
) -> Result<usize, String> {
    if vert_indices.is_empty() {
        return Err("intersect_make_sd_vertices: empty vertex list".to_string());
    }
    // Collect every vertex shape, following existing SD chains so the merged
    // vertex absorbs the whole cluster.
    let mut n_sd: Option<usize> = None;
    let mut shapes: Vec<TopoShape> = Vec::new();
    for &n_x in vert_indices {
        let shape = ds
            .shape(n_x)
            .cloned()
            .ok_or_else(|| format!("intersect_make_sd_vertices: no shape at index {n_x}"))?;
        if let Some(n_sd1) = ds.has_shape_sd(n_x) {
            let sd_shape = ds
                .shape(n_sd1)
                .cloned()
                .ok_or_else(|| format!("intersect_make_sd_vertices: no SD shape at index {n_sd1}"))?;
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
///
/// The live pipeline updates tolerances through
/// [`crate::pave_common::update_vertex_sd`] on the filler; this DS-only variant
/// stays as the isolated-test entry point.
#[cfg(test)]
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


/// Fill the shrunk-range data of a single pave block.
///
/// Port of `BOPAlgo_PaveFiller::FillShrunkData` (simplified): the bound
/// vertices are projected onto the edge's curve and their tolerance spheres
/// carve the shrunk interval `[ts1, ts2]` out of the block range. When a
/// projection fails the bound keeps the block bound. A block is *splittable*
/// when its shrunk interval is longer than `Precision::PConfusion()`.
#[cfg(test)]
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

/// Replace the bound indices of a pave block by their SD (same-domain)
/// representatives.
///
/// Port of `BOPDS_DS::UpdatePaveBlockWithSDVertices`.
#[cfg(test)]
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
