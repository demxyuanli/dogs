use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Face / Face
// ---------------------------------------------------------------------------

/// Whether two planar faces share more than one On/In vertex, so a FaceFace
/// intersection is warranted.
///
/// Source: `BOPAlgo_PaveFiller::CheckPlanes` (`BOPAlgo_PaveFiller_6.cxx`).
pub(super) fn check_planes(ds: &BopdsDS, n_f1: usize, n_f2: usize) -> bool {
    let verts_of = |n_f: usize| -> HashSet<usize> {
        let mut s = HashSet::new();
        if let Some(shape) = ds.shape(n_f) {
            for v in vertices_of(shape) {
                if let Some(i) = ds.index(&v.0) {
                    s.insert(ds.get_same_domain_index(i));
                }
            }
        }
        if let Some(fi) = ds.face_info_pool().iter().find(|fi| fi.face_index == n_f) {
            for &(nv, _, _) in fi.verts() {
                s.insert(ds.get_same_domain_index(nv));
            }
            for &(e, _, _) in fi.paves_in() {
                if let Some(si) = ds.shape_info(e) {
                    for &sub in si.sub_shapes() {
                        if ds.shape_info(sub).map(|x| x.shape_type()) == Some(ShapeType::Vertex) {
                            s.insert(ds.get_same_domain_index(sub));
                        }
                    }
                }
            }
        }
        s
    };
    verts_of(n_f1).intersection(&verts_of(n_f2)).count() > 1
}

/// Core of [`perform_ff`]: intersect every interfering pair of faces and build
/// a section edge for each intersection curve.
pub(super) fn perform_ff_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    starts: &HashMap<(usize, usize), Vec<(f64, f64, f64, f64)>>,
) -> Result<(), String> {
    let pairs = collect_pairs(ds, ShapeType::Face, ShapeType::Face);
    // `BOPAlgo_PaveFiller::PerformFF` (`BOPAlgo_PaveFiller_6.cxx`): update
    // FaceInfo On/In for every F/F pair and every already-touched face
    // (FaceInfo already allocated), then return when there are no F/F pairs.
    let mut fence: HashSet<usize> = HashSet::new();
    for &(n_f1, n_f2) in &pairs {
        fence.insert(n_f1);
        fence.insert(n_f2);
    }
    for i in 0..ds.nb_source_shapes() {
        let is_face = ds
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if is_face && ds.face_info_pool().iter().any(|fi| fi.face_index == i) {
            fence.insert(i);
        }
    }
    ds.update_face_info_on_faces(&fence);
    ds.update_face_info_in_faces(&fence);
    if pairs.is_empty() {
        return Ok(());
    }
    for (n_f1, n_f2) in pairs {
        if ds.has_interf_pair(n_f1, n_f2) {
            continue;
        }
        let Some(f1_shape) = ds.shape(n_f1).cloned() else { continue };
        let Some(f2_shape) = ds.shape(n_f2).cloned() else { continue };
        let face1 = Face(f1_shape);
        let face2 = Face(f2_shape);
        // `BOPAlgo_PaveFiller::CheckPlanes`: two planes are only sent to
        // FaceFace when they already share more than one On/In vertex. A
        // covering pair (cylinder base on box top) is handled by EF IN edges
        // instead of a coplanar FaceFace curve.
        if face_is_planar(&face1) && face_is_planar(&face2) && !check_planes(ds, n_f1, n_f2) {
            // OCCT: Append InterfFF with Init(0,0) and skip FaceFace.
            ds.append_interf_ff(BopdsInterfFf::new(n_f1, n_f2));
            ds.add_interf(n_f1, n_f2);
            continue;
        }
        let tol =
            ctx.fuzzy + BRepTool::face_tolerance(&face1) + BRepTool::face_tolerance(&face2);
        // `IntTools_Curve::TangentialTolerance` (FaceFace ComputeTolReached3d /
        // plane-plane MakeCurve). Copied onto `BOPDS_Curve` for CorrectToleranceOfSE.
        let tang = curve_tangential_tolerance(&face1, &face2);
        let mut ff = FaceFace::new();
        ff.set_face1(face1);
        ff.set_face2(face2);
        ff.set_tolerance(tol);
        if let Some(list) = starts.get(&(n_f1, n_f2)) {
            ff.set_list(list.clone());
        }
        if let Err(msg) = ff.perform() {
            ctx.add_error(format!("perform_ff: face/face intersection failed: {msg}"));
            ds.append_interf_ff(BopdsInterfFf::new(n_f1, n_f2));
            continue;
        }
        if !ff.is_done() {
            ds.append_interf_ff(BopdsInterfFf::new(n_f1, n_f2));
            continue;
        }
        let res = ff.result();
        // `PerformFF` stores `BOPDS_Curve` on `InterfFF`; section edges are
        // born later in `pave_ff::make_blocks_ff` (`MakeBlocks` / `PostTreatFF`).
        let mut rec = BopdsInterfFf::new(n_f1, n_f2);
        rec.init(res.nb_curves(), 0);
        for c in res.curves() {
            let r = c.range;
            if !r.is_valid() || r.length() <= PCONFUSION {
                continue;
            }
            let mut nc = BopdsCurve::new();
            nc.set_curve(c.curve.clone());
            nc.set_pcurves(c.pcurve1.clone(), c.pcurve2.clone());
            nc.set_range(r.first, r.last);
            nc.set_tolerance(tol);
            nc.set_tangential_tolerance(tang);
            let mut bx = occt_core::bnd::BndBox::new();
            bx.add_point(&c.curve.d0(r.first));
            bx.add_point(&c.curve.d0(r.last));
            bx.enlarge(tol);
            nc.set_box(bx);
            nc.init_pave_block1();
            rec.change_curves().push(nc);
        }
        ds.append_interf_ff(rec);
    }
    Ok(())
}

/// `IntTools_Tools::ComputeIntRange` (`IntTools_Tools.cxx:783`).
pub(super) fn compute_int_range(tol1: f64, tol2: f64, angle: f64) -> f64 {
    use occt_core::precision::ANGULAR;
    use std::f64::consts::PI;
    if (PI * 0.5 - angle).abs() < ANGULAR {
        return tol2;
    }
    let an_angle = if angle > PI * 0.5 { PI - angle } else { angle };
    let a1 = tol1 * (PI * 0.5 - an_angle).tan();
    let a2 = tol2 / an_angle.sin();
    a1 + a2
}

/// `IntTools_FaceFace` tangential tolerance of an intersection curve.
/// Non-planar: `max(TolF1, TolF2)`. Plane/plane: `sqrt(Dt^2 + TolF1^2)`,
/// then floored at `max(TolF1, TolF2)` (`ComputeTolReached3d`).
pub(super) fn curve_tangential_tolerance(f1: &Face, f2: &Face) -> f64 {
    let t1 = BRepTool::face_tolerance(f1);
    let t2 = BRepTool::face_tolerance(f2);
    let t_max = t1.max(t2);
    if face_is_planar(f1) && face_is_planar(f2) {
        if let (Some(p1), Some(p2)) = (
            crate::brep_surface::face_plane(f1),
            crate::brep_surface::face_plane(f2),
        ) {
            let ax1 = p1.axis();
            let ax2 = p2.axis();
            let angle = ax1.direction().angle(ax2.direction());
            let dt = compute_int_range(t1, t2, angle);
            return (dt * dt + t1 * t1).sqrt().max(t_max);
        }
    }
    t_max
}

// ---------------------------------------------------------------------------
// ForceInterfEE / ForceInterfEF
// ---------------------------------------------------------------------------

/// Port of the tangent-alignment probe used by OCCT to decide whether the
/// increased vertex tolerance may be used as a fuzzy value when re-intersecting
/// two edges that share their bounding vertices. Two *straight lines* that
/// share both endpoints are automatically coincident, so no probe is needed;
/// for any other pair of curves the tangents at the shared middle point must
/// be nearly parallel (angle within 25°), or the extra tolerance would fuse
/// two merely *touching* edges.
pub(super) fn ee_use_add_tol(
    tools: &mut IntToolsContext,
    curve1: &dyn Curve,
    curve2: &dyn Curve,
    e2: &Edge,
    mid: f64,
) -> bool {
    if curve_is_line(curve1) && curve_is_line(curve2) {
        return true;
    }
    let (_, vt1) = curve1.d1(mid);
    if vt1.square_magnitude() < RESOLUTION {
        return false;
    }
    let vt1 = vt1.normalized();
    let pm = curve1.d0(mid);
    let Some(t2) = tools.project_point_on_edge(e2, &pm) else {
        return false;
    };
    let (_, vt2) = curve2.d1(t2);
    if vt2.square_magnitude() < RESOLUTION {
        return false;
    }
    let cos = vt1.dot(&vt2.normalized());
    cos.abs() >= 0.9063 // cos(25°)
}

/// True if the curve is geometrically a straight line (sampled collinearity).
pub(super) fn curve_is_line(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let p0 = c.d0(a);
    let p1 = c.d0(0.5 * (a + b));
    let p2 = c.d0(b);
    let v1 = GpVec::from_pnts(&p0, &p1);
    let v2 = GpVec::from_pnts(&p0, &p2);
    let m = v1.crossed(&v2).magnitude();
    let scale = v1.magnitude() * v2.magnitude();
    m < 1e-6 * scale.max(1e-12)
}

/// Tolerance of the vertex with index `n_v` (0 when missing).
pub(super) fn vertex_tolerance_of(ds: &BopdsDS, n_v: usize) -> f64 {
    ds.shape(n_v)
        .map(|s| BRepTool::vertex_tolerance(&Vertex(s.clone())))
        .unwrap_or(0.0)
}

/// Index of the common block containing the pave block `pb`, if any.
pub(super) fn pb_common_block_idx(ds: &BopdsDS, pb: &BopdsPaveBlock) -> Option<usize> {
    let (f, l) = pb.range();
    let e = pb.original_edge();
    ds.common_blocks().iter().position(|cb| {
        cb.contains_index(e)
            && cb.ranges()
                .iter()
                .any(|&(a, b)| (a - f).abs() <= PCONFUSION && (b - l).abs() <= PCONFUSION)
    })
}

/// Append one `(edge, range)` entry to a common block, skipping exact duplicates.
pub(super) fn union_pb_into(cb: &mut BopdsCommonBlock, e: usize, f: f64, l: f64) {
    for k in 0..cb.ranges().len() {
        if cb.indices().get(k) == Some(&e) {
            let (a, b) = cb.ranges()[k];
            if (a - f).abs() <= PCONFUSION && (b - l).abs() <= PCONFUSION {
                return;
            }
        }
    }
    cb.add_range(f, l);
    cb.add_index(e);
}

/// Merge every `(edge, range)` entry of `from` into `into`.
pub(super) fn union_common_block(into: &mut BopdsCommonBlock, from: &BopdsCommonBlock) {
    for k in 0..from.ranges().len() {
        if let Some(&e) = from.indices().get(k) {
            let (f, l) = from.ranges()[k];
            union_pb_into(into, e, f, l);
        }
    }
    for &n_f in from.faces() {
        into.add_face(n_f);
    }
}

/// Core of [`force_interf_ee`].
pub(super) fn force_interf_ee_impl(ds: &mut BopdsDS, ctx: &mut FillCtx) -> Result<(), String> {
    // Build the map from bounding-vertex pair to the real pave blocks having
    // those vertices (port of `ForceInterfEE`'s first pass). Only the blocks of
    // non-degenerated edges that actually carry pave blocks participate.
    let mut pb_map: HashMap<(usize, usize), Vec<BopdsPaveBlock>> = HashMap::new();
    let mut fence: HashSet<(usize, usize, usize)> = HashSet::new();
    let n = ds.nb_source_shapes();
    for i in 0..n {
        let Some(si) = ds.shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Edge {
            continue;
        }
        let Some(edge_shape) = ds.shape(i).cloned() else { continue };
        if BRepTool::is_degenerated(&Edge(edge_shape)) {
            continue;
        }
        let pbs = ds.pave_blocks(i).to_vec();
        if pbs.is_empty() {
            continue;
        }
        for pb in pbs {
            let (n_v1, n_v2) = pb.indices();
            if fence.insert((i, n_v1, n_v2)) {
                pb_map.entry((n_v1, n_v2)).or_default().push(pb);
            }
        }
    }

    let mut tools = IntToolsContext::new();
    for ((n_v1, n_v2), group) in pb_map {
        if group.len() < 2 {
            continue;
        }
        // Use the max tolerance of the bounding vertices as the fuzzy addition
        // (the `bSICheckMode` branch — myFuzzyValue — is skipped: the
        // self-interference check pipeline is not ported).
        let a_tol_add = 2.0 * vertex_tolerance_of(ds, n_v1).max(vertex_tolerance_of(ds, n_v2));

        // Check every pair combined from the group.
        for (k1, pb1) in group.iter().enumerate() {
            let n_e1 = pb1.original_edge();
            let (t11, t12) = pb1.range();
            let i_r1 = ds.rank(n_e1);
            let Some(e1_shape) = ds.shape(n_e1).cloned() else { continue };
            let e1 = Edge(e1_shape);
            let Some(curve1) = BRepTool::edge_curve(&e1) else { continue };
            for pb2 in group.iter().skip(k1 + 1) {
                let n_e2 = pb2.original_edge();
                if n_e1 == n_e2 {
                    continue;
                }
                let (t21, t22) = pb2.range();
                let i_r2 = ds.rank(n_e2);

                // Skip pairs of edges from the same argument unless the vertex
                // sharing was acquired during the operation.
                if i_r1 == i_r2 {
                    if (!ds.is_new_shape(n_v1) && ds.rank(n_v1) == i_r1)
                        || (!ds.is_new_shape(n_v2) && ds.rank(n_v2) == i_r2)
                    {
                        continue;
                    }
                }

                // Skip pairs already forming a common block together.
                if let (Some(c1), Some(c2)) =
                    (pb_common_block_idx(ds, pb1), pb_common_block_idx(ds, pb2))
                {
                    if c1 == c2 {
                        continue;
                    }
                }

                let Some(e2_shape) = ds.shape(n_e2).cloned() else { continue };
                let e2 = Edge(e2_shape);
                let Some(curve2) = BRepTool::edge_curve(&e2) else { continue };

                // Tangent gate: only use the increased tolerance when the edges
                // are nearly parallel at the shared middle point.
                let use_add_tol = ee_use_add_tol(
                    &mut tools,
                    curve1.as_ref(),
                    curve2.as_ref(),
                    &e2,
                    (t11 + t12) * 0.5,
                );

                let mut ee = EdgeEdge::with_edges(e1.clone(), e2.clone());
                ee.set_range1(IntRange::new_unchecked(t11, t12));
                ee.set_range2(IntRange::new_unchecked(t21, t22));
                ee.set_fuzzy_value(if use_add_tol { ctx.fuzzy + a_tol_add } else { ctx.fuzzy });
                if let Err(msg) = ee.perform() {
                    ctx.add_error(format!("force_interf_ee: edge/edge intersection failed: {msg}"));
                    continue;
                }
                if !ee.is_done() {
                    continue;
                }
                let cps = ee.common_parts();
                if cps.len() != 1 || cps[0].part_type() != CommonPartType::Edge {
                    continue;
                }

                ds.add_interf_ee(n_e1, n_e2, None);

                // Merge the pair (and the common blocks it already belongs to)
                // into one common block — port of `BOPAlgo_Tools::FillMap` +
                // `PerformCommonBlocks` restricted to this pair.
                let mut cb = BopdsCommonBlock::new();
                union_pb_into(&mut cb, n_e1, t11, t12);
                union_pb_into(&mut cb, n_e2, t21, t22);
                if let Some(c1) = pb_common_block_idx(ds, pb1) {
                    union_common_block(&mut cb, &ds.common_blocks()[c1]);
                }
                if let Some(c2) = pb_common_block_idx(ds, pb2) {
                    union_common_block(&mut cb, &ds.common_blocks()[c2]);
                }
                ds.update_common_block(&cb);
            }
        }
    }
    Ok(())
}

/// Force intersection of the edges after the increase of the tolerance values
/// of their vertices.
///
/// Source: `BOPAlgo_PaveFiller::ForceInterfEE` (`BOPAlgo_PaveFiller_3.cxx`).
/// All real intersections already happened; this stage only looks for
/// *additional common blocks* among the pairs of pave blocks bounded by the
/// same two vertices, using the increased vertex tolerance as a fuzzy value.
pub(crate) fn force_interf_ee(f: &mut PaveFiller) -> Result<(), String> {
    run_with(f, |ds, ctx| force_interf_ee_impl(ds, ctx))
}

/// Read-only access to the face-info entry of the face `n_f`, if any.
pub(super) fn face_info_of(ds: &BopdsDS, n_f: usize) -> Option<&BopdsFaceInfo> {
    ds.face_info_pool().iter().find(|fi| fi.face_index == n_f)
}

/// Add the bounding-vertex indices of the pave block of edge `e` with range
/// `(f, l)` to `out`, when that block exists.
pub(super) fn add_pb_bounds(ds: &BopdsDS, out: &mut HashSet<usize>, e: usize, f: f64, l: f64) {
    for pb in ds.pave_blocks(e) {
        let (f0, l0) = pb.range();
        if (f0 - f).abs() <= PCONFUSION && (l0 - l).abs() <= PCONFUSION {
            let (n_v1, n_v2) = pb.indices();
            out.insert(n_v1);
            out.insert(n_v2);
            return;
        }
    }
}

pub(super) fn pb_already_on_face(fi: &BopdsFaceInfo, pb: &BopdsPaveBlock) -> bool {
    let edges = [pb.edge(), pb.original_edge()];
    let (f0, l0) = pb.range();
    let hit = |e: usize, list: &[(usize, f64, f64)]| {
        e != 0
            && list.iter().any(|&(ee, f, l)| {
                ee == e && (f - f0).abs() <= PCONFUSION && (l - l0).abs() <= PCONFUSION
            })
    };
    for e in edges {
        if hit(e, fi.paves()) || hit(e, fi.paves_in()) || hit(e, fi.paves_on()) {
            return true;
        }
    }
    false
}

pub(super) fn collect_face_verts(ds: &BopdsDS, fi: &BopdsFaceInfo) -> HashSet<usize> {
    let mut verts: HashSet<usize> = HashSet::new();
    for &(v, _, _) in fi.verts() {
        verts.insert(v);
    }
    for &v in fi.verts_in() {
        verts.insert(v);
    }
    for &(e, f, l) in fi.paves().iter().chain(fi.paves_in()).chain(fi.paves_on()) {
        add_pb_bounds(ds, &mut verts, e, f, l);
    }
    verts
}

pub(super) fn collect_source_edge_pbs(ds: &BopdsDS) -> Vec<BopdsPaveBlock> {
    let mut seen: HashSet<(usize, u64, u64)> = HashSet::new();
    let mut out = Vec::new();
    let n = ds.nb_source_shapes();
    for i in 0..n {
        let Some(si) = ds.shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Edge {
            continue;
        }
        if !ds.has_pave_blocks(i) {
            continue;
        }
        let Some(edge_shape) = ds.shape(i) else { continue };
        if BRepTool::is_degenerated(&Edge(edge_shape.clone())) {
            continue;
        }
        for pb in ds.pave_blocks(i) {
            let real = ds.real_pave_block(pb);
            let key = (real.edge(), real.first.to_bits(), real.last.to_bits());
            if seen.insert(key) {
                out.push(real);
            }
        }
    }
    out
}

/// Core of `BOPAlgo_PaveFiller::ForceInterfEF(theMPB, theAddInterf)`.
pub(super) fn force_interf_ef_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    mut the_mpb: Vec<BopdsPaveBlock>,
    add_interf: bool,
    si_check: bool,
) -> Result<(), String> {
    if the_mpb.is_empty() {
        return Ok(());
    }
    for pb in &mut the_mpb {
        ensure_shrunk_data(ds, pb);
    }
    the_mpb.retain(|pb| pb.has_shrunk_data());
    if the_mpb.is_empty() {
        return Ok(());
    }

    let mut tools = IntToolsContext::new();
    let mut hits: Vec<(usize, usize, BopdsPaveBlock, f64)> = Vec::new();
    let n = ds.nb_source_shapes();
    for n_f in 0..n {
        let Some(si_f) = ds.shape_info(n_f) else { continue };
        if si_f.shape_type() != ShapeType::Face {
            continue;
        }
        let Some(fi) = face_info_of(ds, n_f).cloned() else { continue };
        let Some(face_shape) = ds.shape(n_f).cloned() else { continue };
        let face = Face(face_shape);
        let Some(surf) = BRepTool::face_surface(&face) else { continue };
        let face_verts = collect_face_verts(ds, &fi);
        let surf_plane = crate::brep_surface::classify_surface(surf.as_ref())
            == crate::brep_surface::SurfaceKind::Plane;

        for pb in &the_mpb {
            if pb_already_on_face(&fi, pb) {
                continue;
            }
            let (n_v1, n_v2) = pb.indices();
            if !face_verts.contains(&n_v1) || !face_verts.contains(&n_v2) {
                continue;
            }
            let n_e = if pb.edge() != 0 {
                pb.edge()
            } else {
                let orig = pb.original_edge();
                if orig == 0 {
                    continue;
                }
                if ds.rank(n_f) == ds.rank(orig) {
                    continue;
                }
                orig
            };
            if let (Some(box_f), Some(box_e)) = (ds.box_of(n_f), ds.box_of(n_e)) {
                if box_f.is_out_box(box_e) {
                    continue;
                }
            }
            let Some(edge_shape) = ds.shape(n_e).cloned() else { continue };
            let edge = Edge(edge_shape);
            if BRepTool::is_degenerated(&edge) {
                continue;
            }
            let Some(curve) = BRepTool::edge_curve(&edge) else { continue };
            let (ts1, ts2, _) = pb.shrunk_data();
            let mid = intermediate_point(ts1, ts2);
            let (pm, vtg) = curve.d1(mid);
            if vtg.square_magnitude() < RESOLUTION {
                continue;
            }
            let vtg = vtg.normalized();

            let a_tol_check = if si_check {
                ctx.fuzzy
            } else {
                2.0 * vertex_tolerance_of(ds, n_v1).max(vertex_tolerance_of(ds, n_v2))
            };
            let (u, v) = match tools.project_point_on_face(&face, &pm) {
                Ok(uv) => uv,
                Err(_) => continue,
            };
            if surf.d0(u, v).distance(&pm) > a_tol_check + ctx.fuzzy {
                continue;
            }
            if !tools.is_point_in_face(&face, &pm, Some((u, v)), ctx.fuzzy)? {
                continue;
            }

            let mut use_add_tol = true;
            if !surf_plane || !curve_is_line(curve.as_ref()) {
                let p_on_s = surf.d0(u, v);
                let norm = GpVec::from_pnts(&p_on_s, &pm);
                if norm.square_magnitude() > RESOLUTION {
                    let cos = norm.normalized().dot(&vtg);
                    if cos.abs() > 0.4226 {
                        use_add_tol = false;
                    }
                }
            }

            let mut a_tol_add = 0.0;
            if use_add_tol {
                for t in [ts1, ts2] {
                    let p = curve.d0(t);
                    if let Ok((u2, v2)) = tools.project_point_on_face(&face, &p) {
                        let d = surf.d0(u2, v2).distance(&p);
                        if d < a_tol_check && d > a_tol_add {
                            a_tol_add = d;
                        }
                    }
                }
                if a_tol_add > 0.0 {
                    a_tol_add -= BRepTool::edge_tolerance(&edge) + BRepTool::face_tolerance(&face);
                    if a_tol_add < 0.0 {
                        a_tol_add = 0.0;
                    }
                }
            }

            // `myFPBDone` is not stored; `has_interf_pair` stands in for the
            // OCCT Seek miss (`!pMPB || !Contains`).
            let b_intersect = a_tol_add > 0.0 || !ds.has_interf_pair(n_e, n_f);
            if b_intersect {
                hits.push((n_e, n_f, pb.clone(), a_tol_add));
            }
        }
    }

    let mut interfered: Vec<(usize, usize)> = Vec::new();
    let mut new_in: Vec<(usize, usize, f64, f64, bool, BopdsPaveBlock)> = Vec::new();
    for (n_e, n_f, pb, a_tol_add) in hits {
        let (t1, t2) = pb.range();
        let Some(edge_shape) = ds.shape(n_e).cloned() else { continue };
        let Some(face_shape) = ds.shape(n_f).cloned() else { continue };
        let mut ef = EdgeFace::new();
        ef.set_edge(Edge(edge_shape));
        ef.set_face(Face(face_shape));
        ef.set_range(t1, t2);
        ef.set_fuzzy_value(ctx.fuzzy + a_tol_add);
        ef.set_quick_coincidence_check(true);
        if let Err(msg) = ef.perform() {
            ctx.add_error(format!("force_interf_ef: edge/face intersection failed: {msg}"));
            continue;
        }
        if !ef.is_done() || ef.error_status() != 0 {
            continue;
        }
        let cps = ef.common_parts();
        if cps.len() != 1 || cps[0].part_type() != CommonPartType::Edge {
            continue;
        }
        if add_interf {
            interfered.push((n_e, n_f));
        }
        new_in.push((n_e, n_f, t1, t2, add_interf, pb));
    }
    for (n_e, n_f) in interfered {
        ds.add_interf_ef(n_e, n_f, None);
    }
    let mut mpbl = crate::bopalgo_tools::PbFaceListMap::new();
    for (n_e, n_f, f, l, do_cb, pb) in new_in {
        let fi = face_info_mut(ds, n_f);
        if !fi
            .paves_in()
            .iter()
            .any(|&(e, f0, l0)| e == n_e && (f0 - f).abs() <= PCONFUSION && (l0 - l).abs() <= PCONFUSION)
        {
            fi.add_pave_in(n_e, f, l);
        }
        if do_cb {
            ds.add_face_to_common_block(n_e, f, l, n_f);
            crate::bopalgo_tools::fill_map_pb_face(&pb, n_f, &mut mpbl);
        }
    }
    if !mpbl.is_empty() {
        let ictx = IntToolsContext::new();
        crate::bopalgo_tools::perform_common_blocks_faces(&mpbl, ds, &ictx);
    }
    Ok(())
}

/// Force edge/face intersection of the given pave blocks
/// (`BOPAlgo_PaveFiller::ForceInterfEF(theMPB, theAddInterf)`).
pub(crate) fn force_interf_ef_on(
    f: &mut PaveFiller,
    pbs: Vec<BopdsPaveBlock>,
    add_interf: bool,
) -> Result<(), String> {
    let si_check = f.arguments().len() == 1;
    let mut ctx = FillCtx::from_filler(f);
    let r = force_interf_ef_impl(f.ds_mut(), &mut ctx, pbs, add_interf, si_check);
    for e in ctx.errors {
        f.add_error(e);
    }
    r
}

/// Force edge/face intersection after the increase of the tolerance values of
/// their vertices.
///
/// Source: `BOPAlgo_PaveFiller::ForceInterfEF` (`BOPAlgo_PaveFiller_5.cxx`).
/// Looks for additional edge/face common blocks among the pairs of pave blocks
/// whose bounding vertices lie on the face, using the increased vertex
/// tolerance as a fuzzy value.
pub(crate) fn force_interf_ef(f: &mut PaveFiller) -> Result<(), String> {
    // `BOPAlgo_PaveFiller::ForceInterfEF`: nested PostTreatFF fillers are not
    // primary and must not hunt extra E/F common blocks.
    if !f.is_primary() {
        return Ok(());
    }
    let pbs = collect_source_edge_pbs(f.ds());
    force_interf_ef_on(f, pbs, true)
}
