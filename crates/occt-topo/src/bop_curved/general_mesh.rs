use super::prelude::*;
use super::*;

/// A more accurate distance-to-surface than `intpatch::distance_to_surface`:
/// a fine grid search followed by several Newton refinements from neighbouring
/// starts. The intpatch projection has ~2e-3 error near a small sphere, which
/// would swamp the level-set tolerance of the trace.
pub(super) fn distance_to_surface_fine(p: &GpPnt, s: &dyn Surface) -> f64 {
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(s);
    let (nu, nv) = (48usize, 48usize);
    let mut bu = u0;
    let mut bv = v0;
    let mut bd = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).square_distance(p);
            if d < bd {
                bd = d;
                bu = u;
                bv = v;
            }
        }
    }
    let refine = |u: f64, v: f64| {
        let (u2, v2, q) = crate::intpatch::refine_point_on_surface(s, *p, u, v, 30);
        q.distance(p).min(s.d0(u2, v2).distance(p))
    };
    let mut best = refine(bu, bv);
    for (du, dv) in [(0.03, 0.0), (-0.03, 0.0), (0.0, 0.03), (0.0, -0.03)] {
        best = best.min(refine(bu + du, bv + dv));
    }
    best
}

/// Crossing of a segment `a → b` (with signed field values `va, vb`) at the
/// zero level set (used by the windowed marching-squares tracer).
pub(super) fn edge_crossing_2(pa: &GpPnt, va: f64, pb: &GpPnt, vb: f64) -> Option<GpPnt> {
    if va * vb > 0.0 {
        return None;
    }
    if va.abs() < 1e-12 && vb.abs() < 1e-12 {
        return None;
    }
    if va.abs() < 1e-12 {
        return Some(*pa);
    }
    if vb.abs() < 1e-12 {
        return Some(*pb);
    }
    let t = va / (va - vb);
    Some(GpPnt::new(
        pa.x() + t * (pb.x() - pa.x()),
        pa.y() + t * (pb.y() - pa.y()),
        pa.z() + t * (pb.z() - pa.z()),
    ))
}

/// Trace the intersection of `sa` (over the explicit UV window) with `sb`,
/// marching on `sa`'s grid exactly like the general intpatch tracer but over
/// the face's real trimmed window (the clamped `sample_bounds` used by
/// `intersection_curve_points` misses offset faces whose intersection lies
/// outside [-1,1]²).
pub(super) fn face_intersection_polylines(
    sa: &dyn Surface,
    window: (f64, f64, f64, f64),
    sb: &dyn Surface,
    tol: f64,
) -> Vec<Vec<GpPnt>> {
    let (u0, u1, v0, v1) = window;
    if !(u1 > u0 && v1 > v0) {
        return Vec::new();
    }
    let trace_tol = tol.max(1e-4);
    let (nu, nv) = (40usize, 40usize);
    let mut field = vec![vec![0.0f64; nv + 1]; nu + 1];
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = sa.d0(u, v);
            field[i][j] = distance_to_surface_fine(&p, sb) - trace_tol;
        }
    }
    let mut pts: Vec<GpPnt> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let p = [
                sa.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                sa.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                sa.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
                sa.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
            ];
            let f = [field[i][j], field[i + 1][j], field[i + 1][j + 1], field[i][j + 1]];
            let e = [
                edge_crossing_2(&p[0], f[0], &p[1], f[1]),
                edge_crossing_2(&p[1], f[1], &p[2], f[2]),
                edge_crossing_2(&p[2], f[2], &p[3], f[3]),
                edge_crossing_2(&p[3], f[3], &p[0], f[0]),
            ];
            let mut seg = Vec::new();
            for q in e.into_iter().flatten() {
                seg.push(q);
            }
            if seg.len() == 2 {
                pts.push(seg[0]);
                pts.push(seg[1]);
            } else if seg.len() >= 4 {
                pts.push(seg[0]);
                pts.push(seg[2]);
                pts.push(seg[1]);
                pts.push(seg[3]);
            }
        }
    }
    // Link the level-set fragments with a tolerance that covers the marching
    // grid spacing.
    let grid_step = ((u1 - u0) / nu as f64).max((v1 - v0) / nv as f64);
    let link_tol = trace_tol.max(2.0 * grid_step);
    crate::intpatch::chain_intersection_points(&pts, link_tol)
}

/// Axis-aligned bounding box of a face over its real trimmed UV window (unlike
/// [`face_bbox`], which samples the clamped `sample_bounds` and misplaces
/// offset faces).
pub(super) fn face_bbox_local(face: &Face) -> BndBox {
    let mut b = BndBox::new();
    if let Some(s) = BRepTool::face_surface(face) {
        let (u0, u1, v0, v1) = face_uv_window_local(face);
        if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() && u1 > u0 && v1 > v0 {
            for i in 0..=8usize {
                for j in 0..=8usize {
                    let u = u0 + (u1 - u0) * i as f64 / 8.0;
                    let v = v0 + (v1 - v0) * j as f64 / 8.0;
                    b.add_point(&s.d0(u, v));
                }
            }
        }
    }
    b
}

/// Face–face intersection polylines for the trimmed-face path: like
/// [`general_pair_curves`] but only face pairs whose real UV-window bboxes
/// overlap are traced, and the analytic intersection dispatcher
/// (`surface_surface_intersection`) is used first so plane∩sphere and
/// plane∩cylinder give exact curves (the general grid tracer's
/// `distance_to_surface` projection noise swamps the level set near a small
/// sphere).
pub(super) fn general_pair_curves_windowed(fa: &[Face], fb: &[Face], tol: f64) -> Vec<Vec<Vec<Vec<GpPnt>>>> {
    let mut pair_curves: Vec<Vec<Vec<Vec<GpPnt>>>> = vec![vec![Vec::new(); fb.len()]; fa.len()];
    for (i, f_i) in fa.iter().enumerate() {
        let Some(sa) = BRepTool::face_surface(f_i) else { continue };
        let ba = face_bbox_local(f_i);
        let Ok(win_a) = face_uv_window_rect(f_i).map(|r| rect_bounds(&r)) else { continue };
        for (j, f_j) in fb.iter().enumerate() {
            let Some(sb) = BRepTool::face_surface(f_j) else { continue };
            let bb = face_bbox_local(f_j);
            if !ba.is_void() && !bb.is_void() && ba.is_out_box(&bb) {
                continue;
            }
            let mut curves: Vec<Vec<GpPnt>> = Vec::new();
            // The surface_surface_intersection dispatcher handles plane∩sphere,
            // sphere∩sphere and plane∩plane; add the analytic plane∩cylinder
            // (and cone/torus) closed forms so those traces are exact too.
            let cyl_a = crate::int_face_face::cylinder_from_surface(sa.as_ref());
            let cyl_b = crate::int_face_face::cylinder_from_surface(sb.as_ref());
            let plane_a = crate::intpatch::plane_from_surface(sa.as_ref());
            let plane_b = crate::intpatch::plane_from_surface(sb.as_ref());
            // T-28 step 4: the plane∩cylinder closed form lives in `intana`
            // (`IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cylinder)`).
            let cyl_pair = match (plane_a, cyl_b) {
                (Some(pln), Some(cyl)) => Some((pln, cyl)),
                _ => match (cyl_a, plane_b) {
                    (Some(cyl), Some(pln)) => Some((pln, cyl)),
                    _ => None,
                },
            };
            let analytic = cyl_pair.map(|(pln, cyl)| {
                crate::intpatch::ic_list_from_quadric(
                    occt_geom::intana::quadric_quadric_plane_cylinder(&pln, &cyl, 1e-12, 1e-7),
                    sa.as_ref(),
                    sb.as_ref(),
                )
            });
            if let Some(ics) = analytic.filter(|v| !v.is_empty()) {
                curves = ics.iter().map(|ic| ic.points.clone()).collect();
            } else {
                match crate::intpatch::surface_surface_intersection(sa.as_ref(), sb.as_ref(), tol) {
                    crate::intpatch::SurfaceIntersection::Curves(ics) => {
                        curves = ics.iter().map(|ic| ic.points.clone()).collect();
                    }
                    _ => {}
                }
            }
            if curves.is_empty() {
                // Fall back to the windowed marching-squares tracer.
                curves = face_intersection_polylines(sa.as_ref(), win_a, sb.as_ref(), tol);
            }
            if !curves.is_empty() {
                pair_curves[i][j] = curves;
            }
        }
    }
    pair_curves
}

/// The retained triangle mesh of the general curved boolean.
///
/// Returns `None` when no face crosses the boundary (the analytic whole-face
/// path applies); `Some(mesh)` otherwise.
pub(super) fn general_boolean_mesh(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<Option<TriMesh>, String> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    if fa.is_empty() || fb.is_empty() {
        return Err("general_curved_boolean: input has no faces".into());
    }
    let pair_curves = general_pair_curves(&fa, &fb, tol);
    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();
    let has_on = regions_a.contains(&FaceRegion::On) || regions_b.contains(&FaceRegion::On);
    if !has_on {
        return Ok(None);
    }
    boolean_mesh_curves(a, b, op, tol, &pair_curves).map(Some)
}

/// Boolean on solids with general (B-spline / non-analytic) curved faces.
///
/// Computes face–face intersection polylines (via
/// [`crate::intpatch::intersection_curve_points`]) between every overlapping
/// face pair, classifies each face, then rebuilds the result:
///
/// * faces that do not cross the boundary are kept whole (analytic faces stay
///   analytic);
/// * crossing faces are grid-cell classified (a cell is kept when its centroid
///   classifies per the operation) with grid vertices near an intersection
///   polyline snapped onto the curve, so both solids' retained meshes close
///   along the shared intersection;
/// * the welded retained mesh is rebuilt into a faceted BRep solid.
///
/// Planar inputs are delegated to `crate::bop_builder::boolean`.
pub fn general_curved_boolean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }
    match general_boolean_mesh(a, b, op, tol)? {
        Some(mesh) => Ok(result_from_mesh(mesh)),
        None => {
            let fa = faces_of(a);
            let fb = faces_of(b);
            let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
            let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();
            Ok(result_from_kept_faces(&fa, &fb, &regions_a, &regions_b, op))
        }
    }
}

/// The enclosed volume of the general curved boolean result.
///
/// Delegates to [`curved_boolean_volume`] when neither solid has a general
/// curved face; otherwise returns the retained welded mesh volume (exact for a
/// closed, consistently-oriented mesh).
pub fn general_curved_boolean_volume(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<f64, String> {
    if !has_general_curved_face(a) && !has_general_curved_face(b) {
        return curved_boolean_volume(a, b, op, tol);
    }
    if all_faces_planar(a) && all_faces_planar(b) {
        let r = crate::bop_builder::boolean(a, b, op, tol)?;
        return Ok(crate::shape_mesh::shape_volume(&r.shape, 0.02));
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return match op {
            BoolOp::Fuse => Ok(crate::shape_mesh::shape_volume(a, 0.02) + crate::shape_mesh::shape_volume(b, 0.02)),
            BoolOp::Cut => Ok(crate::shape_mesh::shape_volume(a, 0.02)),
            BoolOp::Common => Ok(0.0),
        };
    }
    match general_boolean_mesh(a, b, op, tol)? {
        Some(mesh) => Ok(mesh_volume(&mesh)),
        None => {
            let fa = faces_of(a);
            let fb = faces_of(b);
            let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
            let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();
            let r = result_from_kept_faces(&fa, &fb, &regions_a, &regions_b, op);
            Ok(crate::shape_mesh::shape_volume(&r.shape, 0.02))
        }
    }
}

/// The explicit dispatcher for the curved boolean:
///
/// 1. all faces planar → `crate::bop_builder::boolean` (`BOPAlgo_BOP`);
/// 2. any general (non-planar, non-sphere) curved face → `general_boolean_trimmed`
///    (trimmed B-Rep faces that preserve the analytic surface);
/// 3. otherwise (spheres, and quadric-only solids) → [`curved_boolean`].
pub fn curved_boolean_ext(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }
    if has_general_curved_face(a) || has_general_curved_face(b) {
        let shape = general_boolean_trimmed(a, b, op, tol)?;
        let faces = faces_of(&shape);
        return Ok(boolean_result_from_shape(shape, faces, vec![]));
    }
    curved_boolean(a, b, op, tol)
}

// ---------------------------------------------------------------------------
// Trimmed-face (real B-Rep) general boolean
// ---------------------------------------------------------------------------

/// A face that was trimmed along an intersection curve: the original analytic
/// surface is preserved and the boundary wire follows the kept UV region.
pub struct TrimmedFace {
    pub face: crate::shape::Face,
    pub surface: Arc<dyn Surface>,
    pub kept_region_polygon: Vec<GpPnt2d>,
    pub wire: crate::shape::Wire,
}

/// Assemble a [`BooleanResult`] from an already-built result shape.
pub fn boolean_result_from_shape(shape: TopoShape, faces: Vec<Face>, warnings: Vec<String>) -> BooleanResult {
    let shells: Vec<Shell> = shapes_of(&shape, ShapeType::Shell).into_iter().map(Shell).collect();
    let solid = Solid::wrap(shape.clone());
    BooleanResult { shape, solid, shells, faces, warnings }
}

/// The region of `face` relative to `solid` (reuses [`classify_face_general`]).
pub fn face_kept_region(face: &Face, solid: &TopoShape, op: BoolOp, tol: f64) -> FaceRegion {
    let _ = op;
    classify_face_general(face, solid, tol)
}

// -- 2D region primitives ------------------------------------------------

pub(super) fn unit2d(dx: f64, dy: f64) -> (f64, f64) {
    let len = (dx * dx + dy * dy).sqrt().max(1e-30);
    (dx / len, dy / len)
}

pub(super) fn poly_diag2d(poly: &[GpPnt2d]) -> f64 {
    if poly.is_empty() {
        return 0.0;
    }
    let mut lo = poly[0];
    let mut hi = poly[0];
    for p in poly {
        lo = GpPnt2d::new(lo.x().min(p.x()), lo.y().min(p.y()));
        hi = GpPnt2d::new(hi.x().max(p.x()), hi.y().max(p.y()));
    }
    lo.distance(&hi)
}

/// Unwrap the periodic first UV coordinate (e.g. sphere/cylinder `u` over
/// `[0, 2π]`) so a seam-crossing polyline becomes continuous: each point's `u`
/// is shifted by a multiple of the period to be closest to the previous point.
/// Without this a latitude circle's closing point projects to `u ≈ 0` instead
/// of `2π`, degenerating the polyline.
pub(super) fn unwrap_uv_polyline(pts: &[GpPnt2d], surf: &dyn Surface) -> Vec<GpPnt2d> {
    let (u0, u1) = surf.u_range();
    if !(u0.is_finite() && u1.is_finite()) || u1 - u0 < 3.0 || pts.is_empty() {
        return pts.to_vec();
    }
    let period = u1 - u0;
    let mut out = vec![pts[0]];
    let mut cur = pts[0].x();
    for i in 1..pts.len() {
        let raw = pts[i].x();
        let mut best = raw;
        let mut bd = (raw - cur).abs();
        for k in -2i32..=2 {
            let cand = raw + k as f64 * period;
            let d = (cand - cur).abs();
            if d < bd {
                bd = d;
                best = cand;
            }
        }
        cur = best;
        out.push(GpPnt2d::new(best, pts[i].y()));
    }
    // Place the unwrapped run inside the surface's periodic window. A loop that
    // crosses the seam can unwrap to `[-period, 0]` instead of `[u0, u1]`; it
    // then lies outside the face window (touching it in one point only), so the
    // rect-boundary extension in `cut_two_sides` retraces the polyline itself
    // and both sides degenerate. Aligning the run's centre with the window
    // centre modulo the period puts it back in the window while leaving a run
    // that already sits in the window untouched.
    let lo = out.iter().map(|p| p.x()).fold(f64::INFINITY, f64::min);
    let hi = out.iter().map(|p| p.x()).fold(f64::NEG_INFINITY, f64::max);
    let k = ((0.5 * (u0 + u1) - 0.5 * (lo + hi)) / period).round();
    if k != 0.0 {
        for p in &mut out {
            *p = GpPnt2d::new(p.x() + k * period, p.y());
        }
    }
    out
}

/// Whether a UV polyline is closed: its endpoint gap is at most the trace
/// spacing (the gap between the last and first point of a discretized closed
/// loop is one segment length).
pub(super) fn polyline_is_closed(poly: &[GpPnt2d], tol: f64) -> bool {
    if poly.len() < 2 {
        return false;
    }
    let diag = poly_diag2d(poly);
    let mut seg_sum = 0.0;
    for i in 0..poly.len() - 1 {
        seg_sum += poly[i].distance(&poly[i + 1]);
    }
    let avg_seg = seg_sum / (poly.len() - 1) as f64;
    let gap = poly.first().unwrap().distance(poly.last().unwrap());
    gap <= tol.max(diag * 0.05).max(1.5 * avg_seg)
}



/// Is `p` on the segment `a-b` (within 1e-9)?
pub(super) fn on_segment2d(a: &GpPnt2d, b: &GpPnt2d, p: &GpPnt2d) -> bool {
    let cross = (b.x() - a.x()) * (p.y() - a.y()) - (b.y() - a.y()) * (p.x() - a.x());
    if cross.abs() > 1e-9 {
        return false;
    }
    p.x() >= a.x().min(b.x()) - 1e-9
        && p.x() <= a.x().max(b.x()) + 1e-9
        && p.y() >= a.y().min(b.y()) - 1e-9
        && p.y() <= a.y().max(b.y()) + 1e-9
}

pub(super) fn rect_bounds(rect: &[GpPnt2d]) -> (f64, f64, f64, f64) {
    let u0 = rect.iter().map(|q| q.x()).fold(f64::INFINITY, f64::min);
    let u1 = rect.iter().map(|q| q.x()).fold(f64::NEG_INFINITY, f64::max);
    let v0 = rect.iter().map(|q| q.y()).fold(f64::INFINITY, f64::min);
    let v1 = rect.iter().map(|q| q.y()).fold(f64::NEG_INFINITY, f64::max);
    (u0, u1, v0, v1)
}

/// Extend the ray `p + t·d` (t ≥ 0, `d` a unit 2D direction) until it exits the
/// axis-aligned rectangle `rect`.
pub(super) fn extend_to_rect_boundary(p: GpPnt2d, d: (f64, f64), rect: &[GpPnt2d]) -> GpPnt2d {
    let (u0, u1, v0, v1) = rect_bounds(rect);
    let mut tmin = f64::INFINITY;
    if d.0.abs() > 1e-15 {
        for u in [u0, u1] {
            let t = (u - p.x()) / d.0;
            if t > 1e-12 {
                tmin = tmin.min(t);
            }
        }
    }
    if d.1.abs() > 1e-15 {
        for v in [v0, v1] {
            let t = (v - p.y()) / d.1;
            if t > 1e-12 {
                tmin = tmin.min(t);
            }
        }
    }
    let t = if tmin.is_finite() { tmin } else { 0.0 };
    GpPnt2d::new(p.x() + t * d.0, p.y() + t * d.1)
}

/// The rectangle edge (0..=3) that the boundary point `p` lies on.
pub(super) fn edge_containing(rect: &[GpPnt2d], p: GpPnt2d) -> usize {
    let n = 4;
    for i in 0..n {
        if p.distance(&rect[i]) < 1e-9 {
            return i;
        }
    }
    for i in 0..n {
        if on_segment2d(&rect[i], &rect[(i + 1) % n], &p) {
            return i;
        }
    }
    0
}

/// Walk the rectangle boundary from `from` to `to` (both on the boundary),
/// returning the corner points passed (plus the two endpoints).
pub(super) fn boundary_path(rect: &[GpPnt2d], from: GpPnt2d, to: GpPnt2d, cw: bool) -> Vec<GpPnt2d> {
    let n = 4;
    let ef = edge_containing(rect, from);
    let et = edge_containing(rect, to);
    let mut out = vec![from];
    let mut e = ef;
    loop {
        if e == et {
            break;
        }
        let next = if cw { (e + n - 1) % n } else { (e + 1) % n };
        let corner = if cw { rect[e] } else { rect[(e + 1) % n] };
        if corner.distance(out.last().unwrap()) > 1e-9 {
            out.push(corner);
        }
        e = next;
    }
    if to.distance(out.last().unwrap()) > 1e-9 {
        out.push(to);
    }
    out
}

/// Build the closed cut polygon for a (possibly seam-wrapping) UV polyline:
/// Build the two closed cut polygons for an open UV polyline: extend both
/// endpoints to the rectangle boundary and close along the boundary in each
/// direction. Returns `(left, right)` of the directed polyline, both carrying
/// the polyline's points on the shared edge (so a welded boundary can align the
/// two faces).
pub(super) fn cut_two_sides(rect: &[GpPnt2d], poly: &[GpPnt2d], _tol: f64) -> (Vec<GpPnt2d>, Vec<GpPnt2d>) {
    let p0 = poly[0];
    let p1 = poly[1];
    let pn = poly[poly.len() - 1];
    let pnm1 = poly[poly.len() - 2];
    let d_back = unit2d(p0.x() - p1.x(), p0.y() - p1.y());
    let d_fwd = unit2d(pn.x() - pnm1.x(), pn.y() - pnm1.y());
    let a = extend_to_rect_boundary(p0, d_back, rect);
    let b = extend_to_rect_boundary(pn, d_fwd, rect);
    let base: Vec<GpPnt2d> = std::iter::once(a).chain(poly.iter().copied()).chain(std::iter::once(b)).collect();
    let mut cut1 = base.clone();
    cut1.extend(boundary_path(rect, b, a, true));
    let mut cut2 = base;
    cut2.extend(boundary_path(rect, b, a, false));
    let (mx, my) = ((p0.x() + p1.x()) * 0.5, (p0.y() + p1.y()) * 0.5);
    let (dx, dy) = (p1.x() - p0.x(), p1.y() - p0.y());
    let len = (dx * dx + dy * dy).sqrt().max(1e-12);
    let off = poly_diag2d(rect) * 1e-3;
    let probe = GpPnt2d::new(mx + (-dy / len) * off, my + (dx / len) * off);
    let in1 = point_in_polygon2d(&cut1, &probe);
    let in2 = point_in_polygon2d(&cut2, &probe);
    if in1 && !in2 {
        (cut1, cut2)
    } else if in2 && !in1 {
        (cut2, cut1)
    } else {
        (cut1, cut2)
    }
}

/// Whether a single-loop region is the full rectangle (the initial split).
pub(super) fn region_is_rect(region: &[Vec<GpPnt2d>], rect: &[GpPnt2d]) -> bool {
    if region.len() != 1 || region[0].len() != 4 {
        return false;
    }
    for p in &region[0] {
        if !rect.iter().any(|q| q.distance(p) < 1e-9) {
            return false;
        }
    }
    true
}

/// Split the UV rectangle `rect` by the intersection polyline `poly` into the
/// two region loop-sets `(side_a, side_b)`.
///
/// * a proper closed loop → `side_a` = loop interior, `side_b` = loop exterior;
/// * an open polyline (or a degenerate flat loop, e.g. a latitude circle on a
///   sphere) → `side_a` = the left-hand side of the directed polyline.
pub(super) fn split_uv_region_by_polyline(
    rect: &[GpPnt2d],
    poly: &[GpPnt2d],
    tol: f64,
) -> Result<(Vec<Vec<GpPnt2d>>, Vec<Vec<GpPnt2d>>), String> {
    if poly.len() < 3 || rect.len() < 4 {
        return Err("split_uv_region_by_polyline: degenerate input".into());
    }
    let closed = polyline_is_closed(poly, tol);
    if closed && signed_area2d(poly).abs() > 1e-12 {
        let s1 = polygon_boolean(rect, poly, PolygonBoolOp::Intersect);
        let s2 = polygon_boolean(rect, poly, PolygonBoolOp::Difference);
        return Ok((s1, s2));
    }
    // Open polyline or degenerate flat loop: build the two half-regions directly
    // (they carry the polyline's points on the shared edge, which polygon_boolean's
    // convex fast path would collapse into a straight line).
    let (left, right) = cut_two_sides(rect, poly, tol);
    let mut s1 = Vec::new();
    let mut s2 = Vec::new();
    if left.len() >= 3 {
        s1.push(left);
    }
    if right.len() >= 3 {
        s2.push(right);
    }
    Ok((s1, s2))
}

/// Split a multi-loop region (outer loop + holes) by one UV polyline.
pub(super) fn split_multipolygon_by_polyline(
    region: &[Vec<GpPnt2d>],
    poly: &[GpPnt2d],
    rect: &[GpPnt2d],
    tol: f64,
) -> (Vec<Vec<GpPnt2d>>, Vec<Vec<GpPnt2d>>) {
    let mut a: Vec<Vec<GpPnt2d>> = Vec::new();
    let mut b: Vec<Vec<GpPnt2d>> = Vec::new();
    let closed = polyline_is_closed(poly, tol);
    if closed && signed_area2d(poly).abs() > 1e-12 {
        for loop_poly in region {
            a.extend(polygon_boolean(loop_poly, poly, PolygonBoolOp::Intersect));
            b.extend(polygon_boolean(loop_poly, poly, PolygonBoolOp::Difference));
        }
    } else if region_is_rect(region, rect) {
        let (left, right) = cut_two_sides(rect, poly, tol);
        a.push(left);
        b.push(right);
    } else {
        // A sub-region split by an open polyline: clip the cut to each loop.
        // (The convex fast path may collapse the polyline detail here; this only
        // affects the interior arrangement of a multiply-crossed face.)
        let (left, right) = cut_two_sides(rect, poly, tol);
        for loop_poly in region {
            a.extend(polygon_boolean(loop_poly, &left, PolygonBoolOp::Intersect));
            b.extend(polygon_boolean(loop_poly, &right, PolygonBoolOp::Intersect));
        }
    }
    (a, b)
}

/// Split a face's UV window into an arrangement by all its intersection
/// polylines. Returns the list of region loop-sets.
pub(super) fn split_regions_by_polylines(rect: &[GpPnt2d], polylines_uv: &[Vec<GpPnt2d>], tol: f64) -> Vec<Vec<Vec<GpPnt2d>>> {
    let mut regions: Vec<Vec<Vec<GpPnt2d>>> = vec![vec![rect.to_vec()]];
    for poly in polylines_uv {
        let mut next: Vec<Vec<Vec<GpPnt2d>>> = Vec::new();
        for region in &regions {
            let (a, b) = split_multipolygon_by_polyline(region, poly, rect, tol);
            if !a.is_empty() {
                next.push(a);
            }
            if !b.is_empty() {
                next.push(b);
            }
        }
        regions = next;
        if regions.is_empty() {
            break;
        }
    }
    regions
}

/// Whether a 2D point is inside the multi-loop region (outer + holes).
pub(super) fn region_contains(region: &[Vec<GpPnt2d>], p: &GpPnt2d) -> bool {
    if region.is_empty() || region[0].len() < 3 {
        return false;
    }
    if !point_in_polygon2d(&region[0], p) {
        return false;
    }
    for hole in &region[1..] {
        if hole.len() >= 3 && point_in_polygon2d(hole, p) {
            return false;
        }
    }
    true
}

/// Area-weighted centroid of a 2D polygon (shoelace).
pub(super) fn polygon_centroid_2d(poly: &[GpPnt2d]) -> GpPnt2d {
    let n = poly.len();
    if n == 0 {
        return GpPnt2d::zero();
    }
    let mut a = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let f = poly[i].x() * poly[j].y() - poly[j].x() * poly[i].y();
        a += f;
        cx += (poly[i].x() + poly[j].x()) * f;
        cy += (poly[i].y() + poly[j].y()) * f;
    }
    if a.abs() < 1e-12 {
        return poly[0];
    }
    GpPnt2d::new(cx / (3.0 * a), cy / (3.0 * a))
}
