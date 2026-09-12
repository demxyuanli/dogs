use super::prelude::*;
use super::*;

/// Core of [`perform_ee`]: run edge/edge intersection on every interfering pair
/// of pave blocks and record the intersection vertices and common ranges.
pub(super) fn perform_ee_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
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
                range1_first: pt.uv1().0,
                range1_last: pt.uv1().0,
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
            ds.add_interf_ee(o.n_e1, o.n_e2, None);
        }
        seeds.extend(o.seeds);
        for &(a, b, c, d) in &o.common_ranges {
            ds.add_interf_ee(o.n_e1, o.n_e2, None);
            let mut cb = BopdsCommonBlock::new();
            cb.add_range(a, b);
            cb.add_index(o.n_e1);
            cb.add_range(c, d);
            cb.add_index(o.n_e2);
            ds.update_common_block(&cb);
        }
    }
    if !seeds.is_empty() {
        let created = treat_new_vertices(ds, ctx.fuzzy, &seeds)?;
        // `PerformNewVertices` writes `InterfEE::SetIndexNew` on the record
        // that produced each fused vertex (`BOPAlgo_PaveFiller_3.cxx`).
        let clusters = cluster_seeds(&seeds, ctx.fuzzy);
        for (k, cluster) in clusters.iter().enumerate() {
            let n_v = created[k];
            for &si in cluster {
                let s = &seeds[si];
                ds.bind_ee_new_vertex(s.edge_a, s.edge_b, n_v);
                ds.set_ee_common_range(s.edge_a, s.edge_b, s.range1_first, s.range1_last);
            }
        }
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
    /// `IntTools_CommonPrt::Range1` on the first edge.
    pub range1_first: f64,
    pub range1_last: f64,
}

/// Cluster coincident seeds by point proximity.
pub(super) fn cluster_seeds(seeds: &[NewVertexSeed], fuzzy: f64) -> Vec<Vec<usize>> {
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
pub(super) fn find_source_vertex_at(ds: &BopdsDS, p: &GpPnt, tol: f64) -> Option<usize> {
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
pub(super) fn ensure_face_info(ds: &mut BopdsDS, n_f: usize) {
    let pool = ds.change_face_info_pool();
    if pool.iter().any(|fi| fi.face_index == n_f) {
        return;
    }
    pool.push(BopdsFaceInfo::new(n_f));
}

/// Record the vertex with index `n_v` (projecting to `(u, v)`) on the face `n_f`.
///
/// The vertex is stored in the face-info `verts` list, not `paves`: the
/// section-edge list (`PaveBlocksSc`) must hold only edge indices (OCCT keeps
/// on-face vertices in the separate `VerticesSc`/`VerticesOn`/`VerticesIn`
/// maps of `BOPDS_FaceInfo`).
pub(super) fn record_vertex_on_face(ds: &mut BopdsDS, n_f: usize, n_v: usize, u: f64, v: f64) {
    let pool = ds.change_face_info_pool();
    let info = match pool.iter_mut().find(|fi| fi.face_index == n_f) {
        Some(info) => info,
        None => {
            pool.push(BopdsFaceInfo::new(n_f));
            pool.last_mut().expect("just pushed")
        }
    };
    info.add_vert(n_v, u, v);
}

/// Core of [`perform_vf`]: classify every interfering vertex against its face
/// and record the vertices that lie on the face.
///
/// `pairs_override` has the same semantics as in [`perform_vv_impl`].
pub(super) fn perform_vf_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    pairs_override: Option<&[(usize, usize)]>,
) -> Result<(), String> {
    let pairs: Vec<(usize, usize)> = match pairs_override {
        Some(p) => p.to_vec(),
        None => collect_pairs(ds, ShapeType::Vertex, ShapeType::Face),
    };
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
        ds.add_interf_vf(n_v, n_f, Some(n_vx));
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
pub(crate) fn split_pave_blocks_impl(
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

/// Run an intersection core against the filler, forwarding the recorded errors.
pub(super) fn run_with<F>(f: &mut PaveFiller, body: F) -> Result<(), String>
where
    F: FnOnce(&mut BopdsDS, &mut FillCtx) -> Result<(), String>,
{
    let mut ctx = FillCtx::from_filler(f);
    let r = body(f.ds_mut(), &mut ctx);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Vertex/Vertex intersection: fuse coincident vertices of the arguments into
/// same-domain vertices.
///
/// Source: `BOPAlgo_PaveFiller::PerformVV` (`BOPAlgo_PaveFiller_1.cxx`).
pub fn perform_vv(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_vv::perform_vv(f)
}

/// Vertex/Vertex intersection restricted to the given pairs (used by the
/// repeat-intersection stage — `RepeatIntersection`).
pub(crate) fn perform_vv_pairs(f: &mut PaveFiller, pairs: &[(usize, usize)]) -> Result<(), String> {
    crate::pave_vv::perform_vv_pairs(f, pairs)
}

/// Vertex/Edge intersection: project every interfering vertex onto its edge and
/// insert paves at the projections, splitting the affected pave blocks.
///
/// Source: `BOPAlgo_PaveFiller::PerformVE` / `IntersectVE`
/// (`BOPAlgo_PaveFiller_2.cxx`).
pub fn perform_ve(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ve::perform_ve(f)
}

/// Vertex/Edge intersection restricted to the given pairs (used by the
/// repeat-intersection stage — `RepeatIntersection`).
pub(crate) fn perform_ve_pairs(f: &mut PaveFiller, pairs: &[(usize, usize)]) -> Result<(), String> {
    crate::pave_ve::perform_ve_pairs(f, pairs)
}

/// Edge/Edge intersection: intersect every interfering pair of pave blocks,
/// creating new vertices at the crossing points and recording coincident edge
/// common blocks.
///
/// Source: `BOPAlgo_PaveFiller::PerformEE` + `PerformNewVertices` /
/// `TreatNewVertices` (`BOPAlgo_PaveFiller_3.cxx`).
pub fn perform_ee(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ee_perform::perform_ee(f)
}

/// Vertex/Face intersection: classify every interfering vertex against its face
/// and record the vertices that project onto the face (inside or on the
/// boundary).
///
/// Source: `BOPAlgo_PaveFiller::PerformVF` + `TreatVerticesEE`
/// (`BOPAlgo_PaveFiller_4.cxx`).
pub fn perform_vf(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_vf::perform_vf(f)
}

/// Vertex/Face intersection restricted to the given pairs (used by the
/// repeat-intersection stage — `RepeatIntersection`).
pub(crate) fn perform_vf_pairs(f: &mut PaveFiller, pairs: &[(usize, usize)]) -> Result<(), String> {
    crate::pave_vf::perform_vf_pairs(f, pairs)
}

/// Split the pave blocks of the DS that carry extra paves into elementary
/// blocks, unifying block vertices whose range collapsed.
///
/// Source: `BOPAlgo_PaveFiller::SplitPaveBlocks` (`BOPAlgo_PaveFiller_2.cxx`).
pub fn split_pave_blocks(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_split_blocks::split_pave_blocks(f, &[], true)
}

// ---------------------------------------------------------------------------
// Edge / Face
// ---------------------------------------------------------------------------

/// An edge/face intersection hit collected during the read-only scan.
pub(super) struct EfHit {
    /// Edge index.
    pub(super) n_e: usize,
    /// Face index.
    pub(super) n_f: usize,
    /// Seed for the new intersection vertex.
    pub(super) seed: NewVertexSeed,
}

/// Mutable access to the face-info entry of the face `n_f`, creating it when
/// missing.
pub(super) fn face_info_mut(ds: &mut BopdsDS, n_f: usize) -> &mut BopdsFaceInfo {
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
pub(super) fn record_vertex_point_on_face(ds: &mut BopdsDS, n_f: usize, n_v: usize, p: &GpPnt) {
    let Some(face_shape) = ds.shape(n_f).cloned() else { return };
    let face = Face(face_shape);
    let Some(surf) = BRepTool::face_surface(&face) else { return };
    let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
    record_vertex_on_face(ds, n_f, n_v, u, v);
}

/// Core of [`perform_ef`]: intersect every interfering edge with its face,
/// creating a new vertex at each piercing/touching point and recording
/// coincident sub-ranges as face paves.
pub(super) fn perform_ef_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    distances: &mut HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
) -> Result<(), String> {
    fill_shrunk_data_for_all_edges(ds);
    let pairs = collect_pairs(ds, ShapeType::Edge, ShapeType::Face);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut hits: Vec<EfHit> = Vec::new();
    let mut face_paves: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut interfered: Vec<(usize, usize)> = Vec::new();
    let mut mpbl = crate::bopalgo_tools::PbFaceListMap::new();

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
        let blocks = ds.pave_blocks(n_e).to_vec();
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
        if ef.common_parts().is_empty() {
            let dist = ef.minimal_distance();
            for pb in &blocks {
                crate::pave_ef::record_ef_distance(distances, n_e, n_f, pb, dist, tol_e, tol_f);
            }
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
                            range1_first: t,
                            range1_last: t,
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
        ds.add_interf_ef(n_e, n_f, None);
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
        // `FillMap(aPB, nF)` + `PerformCommonBlocks`: the coincident pave
        // range belongs to a common block that lists this face.
        ds.add_face_to_common_block(n_e, t1, t2, n_f);
        if let Some(pb) = ds.pave_blocks(n_e).iter().find(|pb| {
            let (f0, l0) = pb.range();
            (f0 - t1).abs() <= PCONFUSION && (l0 - t2).abs() <= PCONFUSION
        }) {
            crate::bopalgo_tools::fill_map_pb_face(pb, n_f, &mut mpbl);
        }
    }
    if !mpbl.is_empty() {
        let ictx = IntToolsContext::new();
        crate::bopalgo_tools::perform_common_blocks_faces(&mpbl, ds, &ictx);
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
                ds.bind_ef_new_vertex(h.n_e, h.n_f, n_v);
                ds.set_ef_common_range(h.n_e, h.n_f, h.seed.range1_first, h.seed.range1_last);
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
