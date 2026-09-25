use super::prelude::*;
use super::*;

/// Interior 2D sample points of a UV region: points along the segments from the
/// area-weighted centroid to each vertex (plus triangle centroids), filtered out
/// of holes. Used to classify the region against the other solid. The centroid
/// rays reach deep into thin strips whose vertices all lie on one boundary edge
/// (a sphere latitude cap), so the probes are not ambiguous points on the edge.
pub(super) fn region_interior_points_2d(region: &[Vec<GpPnt2d>]) -> Vec<GpPnt2d> {
    let outer = match region.first() {
        Some(o) if o.len() >= 3 => o,
        _ => return Vec::new(),
    };
    let c = polygon_centroid_2d(outer);
    let mut pts: Vec<GpPnt2d> = Vec::new();
    for v in outer {
        for t in [0.3, 0.6] {
            let p = GpPnt2d::new(c.x() + (v.x() - c.x()) * t, c.y() + (v.y() - c.y()) * t);
            if region_contains(region, &p) {
                pts.push(p);
            }
        }
    }
    let outer3: Vec<GpPnt> = outer.iter().map(|p| GpPnt::new(p.x(), p.y(), 0.0)).collect();
    let tris = triangulate_polygon(&outer3);
    for (a, b, c) in tris {
        let pa = outer[a];
        let pb = outer[b];
        let pc = outer[c];
        let centroid = GpPnt2d::new((pa.x() + pb.x() + pc.x()) / 3.0, (pa.y() + pb.y() + pc.y()) / 3.0);
        if region_contains(region, &centroid) {
            pts.push(centroid);
        }
    }
    if pts.is_empty() && region_contains(region, &outer[0]) {
        pts.push(outer[0]);
    }
    pts
}

/// Classify a UV region of `face` against `other` by sampling interior points.
pub(super) fn region_inside_other(
    _face: &Face,
    surf: &dyn Surface,
    region: &[Vec<GpPnt2d>],
    other: &TopoShape,
    tol: f64,
) -> FaceRegion {
    let pts2d = region_interior_points_2d(region);
    let mut inside = 0usize;
    let mut total = 0usize;
    for p2 in &pts2d {
        let p = surf.d0(p2.x(), p2.y());
        if point_in_solid_curved(p, other, tol) {
            inside += 1;
        }
        total += 1;
    }
    if total == 0 {
        FaceRegion::On
    } else if inside == 0 {
        FaceRegion::Outside
    } else if inside == total {
        FaceRegion::Inside
    } else if inside * 2 >= total {
        FaceRegion::Inside
    } else {
        FaceRegion::Outside
    }
}

/// Map a UV region loop to a 3D point loop on the surface. Planar faces keep
/// one point per vertex (straight edges); curved faces sub-sample UV edges so
/// curved boundaries (e.g. a cylinder's top rim) resolve correctly.
pub(super) fn uv_polygon_to_3d_loop(surf: &dyn Surface, loop_uv: &[GpPnt2d], rect: &[GpPnt2d], _tol: f64) -> Vec<GpPnt> {
    let planar = classify_surface(surf) == SurfaceKind::Plane;
    let step = poly_diag2d(rect) / 32.0;
    let mut out: Vec<GpPnt> = Vec::new();
    for i in 0..loop_uv.len() {
        let a = loop_uv[i];
        let b = loop_uv[(i + 1) % loop_uv.len()];
        let seg_len = a.distance(&b);
        let n = if planar {
            1
        } else {
            (seg_len / step.max(1e-9)).round().max(1.0) as usize
        };
        let n = n.max(1).min(48);
        for k in 0..n {
            let t = k as f64 / n as f64;
            let u = a.x() + (b.x() - a.x()) * t;
            let v = a.y() + (b.y() - a.y()) * t;
            out.push(surf.d0(u, v));
        }
    }
    out
}

/// Build a face from 3D boundary loops (first loop = outer wire, the rest are
/// holes).
pub(super) fn build_face_from_loops(surf: Arc<dyn Surface>, loops_3d: &[Vec<GpPnt>]) -> Result<Face, String> {
    let bld = TopoBuilder::new();
    let mut wires: Vec<Wire> = Vec::new();
    for lp in loops_3d {
        if lp.len() < 3 {
            continue;
        }
        let mut edges: Vec<Edge> = Vec::new();
        for i in 0..lp.len() {
            let j = (i + 1) % lp.len();
            if lp[i].distance(&lp[j]) > 1e-12 {
                edges.push(bld.make_edge_segment(&lp[i], &lp[j]));
            }
        }
        if edges.len() >= 2 {
            let w = bld.make_wire(&edges);
            w.set_closed(true);
            wires.push(w);
        }
    }
    Ok(bld.make_face(surf, &wires))
}

/// The finite UV window of a face as a rectangle polygon.
pub(super) fn face_uv_window_rect(face: &Face) -> Result<Vec<GpPnt2d>, String> {
    let (u0, u1, v0, v1) = face_uv_window_local(face);
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) || u1 <= u0 || v1 <= v0 {
        return Err("face_uv_window_rect: unbounded or degenerate UV window".into());
    }
    Ok(vec![GpPnt2d::new(u0, v0), GpPnt2d::new(u1, v0), GpPnt2d::new(u1, v1), GpPnt2d::new(u0, v1)])
}

/// Trim `face`'s UV window along an intersection polyline (2D points in the
/// face's UV space), keeping the side selected by `keep_inside`.
///
/// `keep_inside = true` keeps the region enclosed by a closed loop (or the
/// left-hand side of a directed open polyline). The returned face keeps the
/// original analytic surface; only the boundary wire is replaced.
pub fn trim_face_to_region(face: &Face, keep_inside: bool, intersection_uvs: &[GpPnt2d], tol: f64) -> Result<Option<Face>, String> {
    let rect = face_uv_window_rect(face)?;
    let (s1, s2) = split_uv_region_by_polyline(&rect, intersection_uvs, tol)?;
    let kept = if keep_inside { s1 } else { s2 };
    if kept.is_empty() {
        return Ok(None);
    }
    let surf = BRepTool::face_surface(face).ok_or("trim_face_to_region: no face surface")?;
    let loops_3d: Vec<Vec<GpPnt>> = kept.iter().map(|lp| uv_polygon_to_3d_loop(surf.as_ref(), lp, &rect, tol)).collect();
    build_face_from_loops(surf, &loops_3d).map(Some)
}

/// Whether a curve is (approximately) a straight line on `[t0, t1]`.
pub(super) fn curve_is_line(c: &dyn Curve, t0: f64, t1: f64) -> bool {
    let p0 = c.d0(t0);
    let p1 = c.d0(0.5 * (t0 + t1));
    let p2 = c.d0(t1);
    let d1 = GpVec::from_pnts(&p0, &p1);
    let d2 = GpVec::from_pnts(&p0, &p2);
    d1.cross_magnitude(&d2) < 1e-6 * (d1.magnitude() * d2.magnitude()).max(1e-12)
}


/// Sample a whole face's boundary into contiguous 3D loops.
///
/// A face whose wire is a single closed edge (a cap) is sampled along that
/// edge. Any other face (planar quads, the cylinder lateral, …) is rebuilt from
/// its finite UV-window rectangle mapped through the surface, which yields the
/// correct quad for planar faces and the correct seam+circle boundary for a
/// cylinder.
pub(super) fn face_boundary_loops_3d(face: &Face) -> Vec<Vec<GpPnt>> {
    let mut out: Vec<Vec<GpPnt>> = Vec::new();
    for w in wires_of_face(face) {
        let edges = edges_of_wire(&w);
        if edges.len() == 1 {
            let s = edge_samples_3d(&edges[0]);
            if !s.is_empty() {
                out.push(s);
            }
        } else if let (Ok(rect), Some(surf)) = (face_uv_window_rect(face), BRepTool::face_surface(face)) {
            let loop3d = uv_polygon_to_3d_loop(surf.as_ref(), &rect, &rect, 1e-6);
            if loop3d.len() >= 3 {
                out.push(loop3d);
            }
        } else {
            // Fallback: dump each edge's samples (used for unbounded faces with
            // no finite UV window).
            for e in &edges {
                let s = edge_samples_3d(e);
                if !s.is_empty() {
                    out.push(s);
                }
            }
        }
    }
    out
}

/// Sample an edge's curve into 3D points (straight lines → 2, otherwise 24).
pub(super) fn edge_samples_3d(e: &Edge) -> Vec<GpPnt> {
    if let Some(c) = BRepTool::edge_curve(e) {
        let (t0, t1) = BRepTool::edge_parameters(e);
        if t0.is_finite() && t1.is_finite() && t1 > t0 {
            let n = if curve_is_line(c.as_ref(), t0, t1) { 2 } else { 24 };
            let mut pts = Vec::with_capacity(n);
            for k in 0..n {
                let t = t0 + (t1 - t0) * k as f64 / (n - 1) as f64;
                pts.push(c.d0(t));
            }
            return pts;
        }
    }
    if let Some((a, b)) = BRepTool::edge_vertices(e) {
        return vec![a, b];
    }
    Vec::new()
}

/// Weld the boundary loops of all kept faces into shared edges, then rebuild
/// every face with the welded wires (so adjacent faces close the shell).
pub(super) fn weld_loops_into_faces(faces_data: Vec<(Arc<dyn Surface>, Vec<Vec<GpPnt>>)>, tol: f64) -> Vec<Face> {
    let bld = TopoBuilder::new();
    // Relative to the geometry scale so the two faces' projections of the same
    // shared trace points (each off by ~1e-3) merge, while distinct curve
    // samples (spaced by the trace grid, ~5e-2) stay separate.
    let mut diag = 0.0f64;
    for (_, loops) in &faces_data {
        for lp in loops {
            if lp.is_empty() {
                continue;
            }
            let mut lo = lp[0];
            let mut hi = lp[0];
            for p in lp {
                lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
                hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
            }
            diag = diag.max(lo.distance(&hi));
        }
    }
    let weld_tol = tol.max(diag * 2e-3).max(1e-4);
    let cell = |p: &GpPnt| -> (i64, i64, i64) {
        (
            f64::floor(p.x() / weld_tol) as i64,
            f64::floor(p.y() / weld_tol) as i64,
            f64::floor(p.z() / weld_tol) as i64,
        )
    };
    let mut points: Vec<GpPnt> = Vec::new();
    let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> = std::collections::HashMap::new();
    let mut idx_loops: Vec<Vec<Vec<usize>>> = Vec::new();
    for (_, loops) in &faces_data {
        let mut face_idx: Vec<Vec<usize>> = Vec::new();
        for loop_pts in loops {
            let mut idx: Vec<usize> = Vec::with_capacity(loop_pts.len());
            for p in loop_pts {
                let k = cell(p);
                let mut found = None;
                'search: for dx in -1i64..=1 {
                    for dy in -1i64..=1 {
                        for dz in -1i64..=1 {
                            if let Some(bucket) = grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                                for &j in bucket {
                                    if points[j].distance(p) <= weld_tol {
                                        found = Some(j);
                                        break 'search;
                                    }
                                }
                            }
                        }
                    }
                }
                match found {
                    Some(j) => idx.push(j),
                    None => {
                        let j = points.len();
                        points.push(*p);
                        grid.entry(k).or_default().push(j);
                        idx.push(j);
                    }
                }
            }
            face_idx.push(idx);
        }
        idx_loops.push(face_idx);
    }
    let mut edge_map: std::collections::HashMap<(usize, usize), Edge> = std::collections::HashMap::new();
    let mut face_edges: Vec<Vec<Vec<Edge>>> = Vec::new();
    for face_idx in &idx_loops {
        let mut loops_edges: Vec<Vec<Edge>> = Vec::new();
        for idx in face_idx {
            let mut edges: Vec<Edge> = Vec::new();
            let m = idx.len();
            for i in 0..m {
                let j = idx[i];
                let k = idx[(i + 1) % m];
                if j == k {
                    continue;
                }
                let key = (j.min(k), j.max(k));
                let e = edge_map.entry(key).or_insert_with(|| bld.make_edge_segment(&points[j], &points[k])).clone();
                edges.push(e);
            }
            loops_edges.push(edges);
        }
        face_edges.push(loops_edges);
    }
    let mut faces: Vec<Face> = Vec::new();
    for (i, (surf, _)) in faces_data.iter().enumerate() {
        let wires: Vec<Wire> = face_edges[i]
            .iter()
            .map(|edges| {
                let w = bld.make_wire(edges);
                w.set_closed(true);
                w
            })
            .collect();
        faces.push(bld.make_face(surf.clone(), &wires));
    }
    faces
}

/// Boolean producing trimmed B-Rep faces: each kept face keeps its original
/// NURBS/analytic surface, with the boundary wire trimmed along the shared
/// face–face intersection curves. Planar inputs delegate to
/// [`crate::bop_builder::boolean`].
pub fn general_boolean_trimmed(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<TopoShape, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol).map(|r| r.shape);
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op).shape);
    }
    let fa = faces_of(a);
    let fb = faces_of(b);
    if fa.is_empty() || fb.is_empty() {
        return Err("general_boolean_trimmed: input has no faces".into());
    }
    let pair_curves = general_pair_curves_windowed(&fa, &fb, tol);
    let mut faces_data: Vec<(Arc<dyn Surface>, Vec<Vec<GpPnt>>)> = Vec::new();

    // Faces of A, then faces of B.
    for (from_a, src, other, other_len, row) in [
        (true, a, b, fb.len(), &pair_curves),
        (false, b, a, fa.len(), &pair_curves),
    ] {
        let faces = faces_of(src);
        for (i, f) in faces.iter().enumerate() {
            let Some(surf) = BRepTool::face_surface(f) else { continue };
            let mut polylines_3d: Vec<Vec<GpPnt>> = Vec::new();
            for j in 0..other_len {
                for poly in if from_a { &row[i][j] } else { &row[j][i] } {
                    polylines_3d.push(poly.clone());
                }
            }
            let r = classify_face_general(f, other, tol);
            if r != FaceRegion::On || polylines_3d.is_empty() {
                if select_keep(r, from_a, op) {
                    faces_data.push((surf.clone(), face_boundary_loops_3d(f)));
                }
                continue;
            }
            let rect = face_uv_window_rect(f)?;
            let polylines_uv: Vec<Vec<GpPnt2d>> = polylines_3d
                .iter()
                .map(|poly| {
                    let raw: Vec<GpPnt2d> = poly
                        .iter()
                        .map(|p| {
                            let (u, v) = crate::intpatch::project_params(surf.as_ref(), p);
                            GpPnt2d::new(u, v)
                        })
                        .collect();
                    unwrap_uv_polyline(&raw, surf.as_ref())
                })
                .collect();
            let regions = split_regions_by_polylines(&rect, &polylines_uv, tol);
            for region in &regions {
                let reg_class = region_inside_other(f, surf.as_ref(), region, other, tol);
                if select_keep(reg_class, from_a, op) {
                    let loops_3d: Vec<Vec<GpPnt>> = region
                        .iter()
                        .map(|lp| uv_polygon_to_3d_loop(surf.as_ref(), lp, &rect, tol))
                        .collect();
                    faces_data.push((surf.clone(), loops_3d));
                }
            }
        }
    }

    let bld = TopoBuilder::new();
    if faces_data.is_empty() {
        let comp = bld.make_compound_of(&[]);
        return Ok(comp.0);
    }
    let faces = weld_loops_into_faces(faces_data, tol);
    let shell = bld.make_shell(&faces);
    let closed = crate::shell_check::shell_is_closed(&shell);
    let solid = if closed { Some(bld.make_solid(&[shell.clone()])) } else { None };
    Ok(solid.map(|s| s.0).unwrap_or_else(|| shell.0.clone()))
}

/// The full boolean dispatcher:
///
/// 1. all faces planar → `crate::bop_builder::boolean` (`BOPAlgo_BOP`);
/// 2. any general (non-planar, non-sphere) curved face → [`general_boolean_trimmed`]
///    (trimmed B-Rep faces that preserve the analytic surface);
/// 3. otherwise (spheres, and quadric-only solids) → [`curved_boolean`].
pub fn curved_boolean_full(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let result = if all_faces_planar(a) && all_faces_planar(b) {
        crate::bop_builder::boolean(a, b, op, tol)?
    } else if has_general_curved_face(a) || has_general_curved_face(b) {
        let shape = general_boolean_trimmed(a, b, op, tol)?;
        let faces = faces_of(&shape);
        boolean_result_from_shape(shape, faces, vec![])
    } else {
        curved_boolean(a, b, op, tol)?
    };
    provision_face_pcurves(&result.shape);
    Ok(result)
}

/// OCCT's boolean result carries a pcurve on every edge of every face: the
/// splitter copies each `BRep_GCurve` onto the trimmed faces it produces
/// (`BOPAlgo_Builder::BuildSplitFaces` / `BRepTools_Modifier`), and
/// `ShapeFix_Edge::FixAddPCurve` (`ShapeFix_Edge.cxx:517-534`) *projects* one
/// whenever an edge of the result has none. The port's trimmed faces are built
/// from scratch, so run that projection pass here — without it the analytic
/// passes (`BRepGProp`) fall back to the 3D→UV map on the curved faces and
/// mis-integrate the result.
fn provision_face_pcurves(shape: &TopoShape) {
    for f in faces_of(shape) {
        let mut wires = wires_of_face(&f);
        for w in wires.iter_mut() {
            crate::shhealing::check_pcurves_and_shift(w, &f, occt_core::precision::CONFUSION);
        }
    }
}
