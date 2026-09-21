use super::prelude::*;
use super::*;

/// A point where two edges meet: the parameter on each edge and the 3D point.
#[derive(Debug, Clone, Copy, PartialEq)]

pub struct EdgeEdgeHit {
    pub u1: f64,
    pub u2: f64,
    pub point: GpPnt,
}

// ---------------------------------------------------------------------------
// Edge geometry helpers (edge-parameter-space evaluation)
// ---------------------------------------------------------------------------

/// Evaluate the edge's geometry at a parameter in the *edge's* parameter
/// space (`BRepTool::edge_parameters`). `GeomTrimmedCurve` reparametrizes its
/// basis to `[0, 1]`, so those are mapped back into the edge's stored range.
pub(super) fn eval_edge(e: &Edge, u: f64) -> GpPnt {
    let Some(curve) = BRepTool::edge_curve(e) else {
        return GpPnt::zero();
    };
    let (ef, el) = BRepTool::edge_parameters(e);
    if curve.first_parameter() == 0.0 && curve.last_parameter() == 1.0 {
        let t = if (el - ef).abs() > 1e-15 {
            ((u - ef) / (el - ef)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        curve.d0(t)
    } else {
        curve.d0(u)
    }
}

/// Map a *curve*-space parameter onto the edge's parameter space. Only
/// `GeomTrimmedCurve` (range `[0, 1]`) differs from the edge range here.
pub(super) fn curve_to_edge_param(curve: &dyn Curve, e_first: f64, e_last: f64, u_curve: f64) -> f64 {
    if curve.first_parameter() == 0.0 && curve.last_parameter() == 1.0 {
        e_first + u_curve * (e_last - e_first)
    } else {
        u_curve
    }
}

/// Whether a curve is geometrically a straight line. Unbounded parameter
/// ranges belong only to `GeomLine` in this port; bounded curves are checked
/// for collinearity of eight samples against the first–last chord.
pub(super) fn is_line_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return true;
    }
    if (b - a).abs() <= 1e-15 {
        return false;
    }
    let n = 8;
    let p0 = c.d0(a);
    let pl = c.d0(b);
    // Curve extent: a closed curve (full circle) has a ~zero chord, so it can
    // never be line-like.
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let d0 = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-6 * size;
    for i in 1..n {
        let p = c.d0(a + (b - a) * i as f64 / n as f64);
        if GpVec::from_pnts(&p0, &p).cross_magnitude(&d0) > tol * size {
            return false;
        }
    }
    true
}

/// Whether a curve is geometrically a (planar) circle: every sample is
/// coplanar and equidistant from a common center.
pub(super) fn is_circle_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = match first_three_spanning(&pts) {
        Some(x) => x,
        None => return false,
    };
    let center = match circumcenter(&p0, &p1, &p2) {
        Some(c) => c,
        None => return false,
    };
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return false;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return false;
    }
    let n = nrm.divided(m);
    let scale = radius.max(1.0).max((b - a).abs());
    let tol = 1e-6 * scale;
    for p in pts {
        let v = GpVec::from_pnts(&center, &p);
        if v.dot(&n).abs() > tol {
            return false;
        }
        if (v.magnitude() - radius).abs() > tol {
            return false;
        }
    }
    true
}

/// `(center, radius, unit plane normal)` of a circle-like curve.
pub(super) fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = first_three_spanning(&pts)?;
    let center = circumcenter(&p0, &p1, &p2)?;
    let radius = p0.distance(&center);
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return None;
    }
    Some((center, radius, nrm.divided(m)))
}

/// First three non-collinear points of a sample set, if they exist.
pub(super) fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
    let p0 = pts[0];
    let mut i1 = None;
    for (i, p) in pts.iter().enumerate().skip(1) {
        if GpVec::from_pnts(&p0, p).magnitude() > 1e-9 {
            i1 = Some(i);
            break;
        }
    }
    let i1 = i1?;
    let p1 = pts[i1];
    let d0 = GpVec::from_pnts(&p0, &p1);
    for p in pts.iter().skip(i1 + 1) {
        if GpVec::from_pnts(&p0, p).cross_magnitude(&d0) > 1e-9 * d0.magnitude().max(1e-9) {
            return Some((p0, p1, *p));
        }
    }
    None
}

pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

pub(super) fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear points (perpendicular-bisector system).
pub(super) fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.crossed(&d2);
    if n.magnitude() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.xyz().x, n.xyz().y, n.xyz().z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n.xyz())];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

/// Parameter (in edge space) of the point nearest to `p` on the edge, when it
/// lies on the edge within tolerance.
pub(super) fn param_on_edge(e: &Edge, curve: &dyn Curve, p: &GpPnt, a: f64, b: f64, tol: f64) -> Option<f64> {
    let etol = tol.max(1e-7);
    let u_curve = if is_line_like(curve) {
        let pa = eval_edge(e, a);
        let pb = eval_edge(e, b);
        let d = GpVec::from_pnts(&pa, &pb);
        let ll = d.square_magnitude();
        if ll <= 1e-30 {
            return None;
        }
        let t = (GpVec::from_pnts(&pa, p).dot(&d) / ll).clamp(0.0, 1.0);
        let proj = pa.translated_vec(&d.multiplied_scalar(t));
        if proj.distance(p) > etol {
            return None;
        }
        a + t * (b - a)
    } else {
        let pc = geom_api::project_point_on_curve(curve, p, etol)?;
        if pc.distance > etol {
            return None;
        }
        curve_to_edge_param(curve, a, b, pc.parameter)
    };
    let (lo, hi) = (a.min(b), a.max(b));
    if u_curve < lo - etol || u_curve > hi + etol {
        return None;
    }
    let u = u_curve.clamp(lo, hi);
    if eval_edge(e, u).distance(p) > etol {
        return None;
    }
    Some(u)
}

// ---------------------------------------------------------------------------
// Edge–edge intersection
// ---------------------------------------------------------------------------

/// All intersection points of two edges, as `(param_on_1, param_on_2, point)`.
///
/// Line–line, line–circle and coplanar circle–circle are solved exactly; every
/// other combination runs the faithful `IntTools_EdgeEdge`
/// ([`crate::edge_edge::EdgeEdge`]) — OCCT's only route for those pairs
/// (`IntTools_EdgeEdge.cxx:185-243`).
pub fn edge_edge_intersections(e1: &Edge, e2: &Edge, tol: f64) -> Vec<EdgeEdgeHit> {
    let Some(c1) = BRepTool::edge_curve(e1) else {
        return Vec::new();
    };
    let Some(c2) = BRepTool::edge_curve(e2) else {
        return Vec::new();
    };
    let (a1, b1) = BRepTool::edge_parameters(e1);
    let (a2, b2) = BRepTool::edge_parameters(e2);
    if !a1.is_finite() || !b1.is_finite() || !a2.is_finite() || !b2.is_finite() {
        return Vec::new();
    }

    let l1 = is_line_like(&*c1);
    let l2 = is_line_like(&*c2);
    let ci1 = !l1 && is_circle_like(&*c1);
    let ci2 = !l2 && is_circle_like(&*c2);

    let mut points: Vec<GpPnt> = Vec::new();
    match (l1, l2) {
        (true, true) => {
            if let Some(p) = segment_segment_intersection_3d(
                &eval_edge(e1, a1),
                &eval_edge(e1, b1),
                &eval_edge(e2, a2),
                &eval_edge(e2, b2),
            ) {
                points.push(p);
            }
        }
        (true, false) => {
            if ci2 {
                line_circle_hits(e1, &*c2, tol, &mut points);
            } else {
                points = edge_edge_general_points(e1, e2);
            }
        }
        (false, true) => {
            if ci1 {
                line_circle_hits(e2, &*c1, tol, &mut points);
            } else {
                points = edge_edge_general_points(e1, e2);
            }
        }
        (false, false) => {
            if ci1 && ci2 {
                match circle_circle_exact(&*c1, &*c2, tol) {
                    Some(pts) => points = pts,
                    None => points = edge_edge_general_points(e1, e2),
                }
            } else {
                points = edge_edge_general_points(e1, e2);
            }
        }
    }

    let mut hits = Vec::new();
    for p in points {
        let u1 = param_on_edge(e1, &*c1, &p, a1, b1, tol);
        let u2 = param_on_edge(e2, &*c2, &p, a2, b2, tol);
        if let (Some(u1), Some(u2)) = (u1, u2) {
            hits.push(EdgeEdgeHit { u1, u2, point: p });
        }
    }
    dedupe_hits(hits, tol)
}

/// Line–circle intersection points (both coplanar and transverse cases).
pub(super) fn line_circle_hits(e_line: &Edge, c_circ: &dyn Curve, tol: f64, out: &mut Vec<GpPnt>) {
    let (al, bl) = BRepTool::edge_parameters(e_line);
    let pa = eval_edge(e_line, al);
    let pb = eval_edge(e_line, bl);
    let d = GpVec::from_pnts(&pa, &pb);
    let l = d.magnitude();
    if l <= 1e-30 {
        return;
    }
    let du = d.divided(l);
    let (center, radius, n) = match circle_geometry(c_circ) {
        Some(x) => x,
        None => return,
    };
    let etol = tol.max(1e-9);
    let denom = n.dot(&du);
    if denom.abs() <= 1e-12 {
        // Line parallel to the circle's plane: intersect only if it lies in it.
        if n.dot(&GpVec::from_pnts(&center, &pa)).abs() > etol {
            return;
        }
        let ac = GpVec::from_pnts(&center, &pa);
        let bq = du.dot(&ac); // du·(pa - center)
        let cq = ac.square_magnitude() - radius * radius;
        let disc = bq * bq - cq; // (b/2)² − c with leading coefficient 1
        if disc < 0.0 {
            return;
        }
        let sq = disc.sqrt();
        for s in [-bq - sq, -bq + sq] {
            if s >= -etol && s <= l + etol {
                out.push(pa.translated_vec(&du.multiplied_scalar(s)));
            }
        }
    } else {
        // Transverse: the line meets the circle's plane at one point; the
        // point is on the circle iff it lies on the rim.
        let s = n.dot(&GpVec::from_pnts(&pa, &center)) / denom;
        let p0 = pa.translated_vec(&du.multiplied_scalar(s));
        if (p0.distance(&center) - radius).abs() <= etol {
            out.push(p0);
        }
    }
}

/// Coplanar circle–circle intersection via the radical line. `None` when the
/// circles are not coplanar (caller falls back to sampling).
pub(super) fn circle_circle_exact(c1: &dyn Curve, c2: &dyn Curve, tol: f64) -> Option<Vec<GpPnt>> {
    let (c1c, r1, n1) = circle_geometry(c1)?;
    let (c2c, r2, n2) = circle_geometry(c2)?;
    let etol = tol.max(1e-9);
    if n1.crossed(&n2).magnitude() > 1e-9 {
        return None; // planes cross → non-coplanar, not handled exactly here
    }
    if n1.dot(&GpVec::from_pnts(&c1c, &c2c)).abs() > etol {
        return Some(Vec::new()); // parallel, distinct planes
    }
    // Coplanar: radical line in the shared plane.
    let cc = GpVec::from_pnts(&c1c, &c2c);
    let dist = cc.magnitude();
    if dist <= 1e-30 {
        return Some(Vec::new()); // concentric
    }
    let u = cc.divided(dist);
    let x = (r1 * r1 - r2 * r2 + dist * dist) / (2.0 * dist);
    let h2 = r1 * r1 - x * x;
    if h2 < -etol * r1.max(1.0) {
        return Some(Vec::new());
    }
    let h = if h2 > 0.0 { h2.sqrt() } else { 0.0 };
    let w = u.crossed(&n1); // unit, perpendicular to the center line in-plane
    let base = c1c.translated_vec(&u.multiplied_scalar(x));
    let mut out = vec![base.translated_vec(&w.multiplied_scalar(h))];
    if h > etol {
        out.push(base.translated_vec(&w.multiplied_scalar(-h)));
    }
    Some(out)
}

pub(super) fn dedupe_hits(mut hits: Vec<EdgeEdgeHit>, tol: f64) -> Vec<EdgeEdgeHit> {
    let etol = tol.max(1e-9);
    let mut out: Vec<EdgeEdgeHit> = Vec::new();
    for h in hits.drain(..) {
        let dup = out
            .iter()
            .any(|o| o.point.distance(&h.point) <= etol && (o.u1 - h.u1).abs() <= etol);
        if !dup {
            out.push(h);
        }
    }
    out
}

/// `IntTools_EdgeEdge` for the non-analytic curve pairs: the faithful
/// `FindSolutions` parameter-box recursion of [`crate::edge_edge::EdgeEdge`]
/// (`IntTools_EdgeEdge.cxx:185-243`).
///
/// Replaces `geom_api::curve_curve_intersections` (a 256×256 sampler with
/// alternating 1-D minimization, audit A16). The edge tolerances and the fuzzy
/// value drive the engine exactly as OCCT's `Prepare` does; the caller's `tol`
/// is not an `IntTools_EdgeEdge` input.
fn edge_edge_general_points(e1: &Edge, e2: &Edge) -> Vec<GpPnt> {
    let mut ee = crate::edge_edge::EdgeEdge::with_edges(e1.clone(), e2.clone());
    if ee.perform().is_err() {
        return Vec::new();
    }
    ee.points().iter().map(|p| p.pnt1).collect()
}

/// Non-planar edge–face points via the faithful curve/surface intersector.
///
/// `IntTools_EdgeFace.cxx:426-445`: `IntCurveSurface_HInter anExactIntersector;
/// anExactIntersector.Perform(aCurve, aSurface);` and every point whose `W()`
/// lies in `[aTF, aTL]` is kept. The port calls
/// [`crate::intcurvesurface::perform_curve_surface`] (that same algorithm) and
/// leaves the range filter to the caller's `param_on_edge`, which requires the
/// point to project back onto the edge inside `[a, b]`.
///
/// Replaces `geom_api::curve_surface_intersections` (curve sampling + distance
/// dip + golden section, audit A16).
fn curve_surface_points(
    curve: &dyn Curve,
    face: &Face,
    surf: &dyn Surface,
    a: f64,
    b: f64,
) -> Vec<GpPnt> {
    // OCCT hands the intersector a bare `GeomAdaptor_Surface` (the surface's own
    // bounds); the port's engine needs a finite window, so it gets the face's UV
    // box — the same window `crate::int_curves_face` uses.
    let (u0, u1, v0, v1) = crate::int_curves_face::finite_uv(face);
    // `Perform(aCurve, aSurface)` uses the curve adaptor's own range; an
    // unbounded curve (an edge on an infinite line) is restricted to the edge
    // range, which the caller's filter keeps anyway.
    let (c0, c1) = (curve.first_parameter(), curve.last_parameter());
    let cu_range = if c0.is_finite() && c1.is_finite() { (c0, c1) } else { (a, b) };
    match crate::intcurvesurface::perform_curve_surface(curve, surf, cu_range, (u0, u1, v0, v1)) {
        Ok(res) => res.points().iter().map(|p| p.pnt).collect(),
        Err(_) => Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Edge–face intersection
// ---------------------------------------------------------------------------

/// All points where an edge meets a face: `(param_on_edge, point)`.
///
/// Planar faces use the exact line/circle–plane solves and keep only the hits
/// inside the face's 2D boundary polygon. Non-planar faces use the faithful
/// `IntCurveSurface_HInter` ([`crate::intcurvesurface::perform_curve_surface`]),
/// exactly as `IntTools_EdgeFace::Perform` does
/// (`IntTools_EdgeFace.cxx:426-445`).
pub fn edge_face_intersections(e: &Edge, f: &Face, tol: f64) -> Vec<(f64, GpPnt)> {
    let Some(curve) = BRepTool::edge_curve(e) else {
        return Vec::new();
    };
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() {
        return Vec::new();
    }

    if let Some(pln) = face_plane(f) {
        let poly = face_polygon_2d_with_plane(f, &pln);
        let mut points: Vec<GpPnt> = Vec::new();
        if is_line_like(&*curve) {
            line_plane_hits(e, &pln, a, b, tol, poly.as_deref(), &mut points);
        } else if is_circle_like(&*curve) {
            circle_plane_hits(e, &*curve, &pln, a, b, tol, &mut points);
        } else {
            general_plane_hits(e, &pln, a, b, tol, &mut points);
        }
        let mut hits = Vec::new();
        for p in points {
            if let Some(u) = param_on_edge(e, &*curve, &p, a, b, tol) {
                let inside = poly
                    .as_ref()
                    .map_or(true, |poly| point_in_polygon2d(poly, &project_to_2d(&pln, &p)));
                if inside {
                    hits.push((u, p));
                }
            }
        }
        dedupe_edge_face(hits, tol)
    } else {
        // Non-planar face.
        let Some(surf) = BRepTool::face_surface(f) else {
            return Vec::new();
        };
        let pts = curve_surface_points(&*curve, f, &*surf, a, b);
        let mut hits = Vec::new();
        for p in pts {
            if let Some(u) = param_on_edge(e, &*curve, &p, a, b, tol) {
                hits.push((u, p));
            }
        }
        dedupe_edge_face(hits, tol)
    }
}

/// Intersect a line-like edge with a plane; `poly` (in the plane's frame)
/// restricts an in-plane edge to its polygon-overlapping segment.
pub(super) fn line_plane_hits(
    e: &Edge,
    pln: &GpPln,
    a: f64,
    b: f64,
    tol: f64,
    poly: Option<&[GpPnt2d]>,
    out: &mut Vec<GpPnt>,
) {
    let pa = eval_edge(e, a);
    let pb = eval_edge(e, b);
    let d = GpVec::from_pnts(&pa, &pb);
    let l = d.magnitude();
    if l <= 1e-30 {
        return;
    }
    let du = d.divided(l);
    let n = GpVec::from_xyz(pln.axis().direction().xyz());
    let q = pln.location();
    let etol = tol.max(1e-9);
    let denom = n.dot(&du);
    if denom.abs() <= 1e-12 {
        // Parallel: intersects only if it lies in the plane.
        if n.dot(&GpVec::from_pnts(&q, &pa)).abs() > etol {
            return;
        }
        match poly {
            Some(poly) => {
                for (t0, t1) in line_polygon_t_intervals(pln, &pa, &du, poly) {
                    let u0 = a + t0.max(0.0);
                    let u1 = a + t1.min(l);
                    if u1 - u0 > etol {
                        out.push(eval_edge(e, u0));
                        out.push(eval_edge(e, u1));
                    }
                }
            }
            None => {
                out.push(pa);
                out.push(pb);
            }
        }
        return;
    }
    let s = n.dot(&GpVec::from_pnts(&pa, &q)) / denom;
    if s >= -etol && s <= l + etol {
        out.push(pa.translated_vec(&du.multiplied_scalar(s)));
    }
}

/// Intersect a circle-like edge with a plane: the circle's plane and the face
/// plane meet in a line (or coincide), giving up to two rim points.
pub(super) fn circle_plane_hits(
    e: &Edge,
    curve: &dyn Curve,
    pln: &GpPln,
    _a: f64,
    _b: f64,
    tol: f64,
    out: &mut Vec<GpPnt>,
) {
    let _ = e;
    let (center, radius, n) = match circle_geometry(curve) {
        Some(x) => x,
        None => return,
    };
    let np = GpVec::from_xyz(pln.axis().direction().xyz());
    let q = pln.location();
    let etol = tol.max(1e-9);
    let cross = n.crossed(&np);
    if cross.magnitude() <= 1e-9 {
        // Circle's plane parallel to the face plane: coplanar → the whole
        // circle lies in the face (deferred), otherwise no intersection.
        return;
    }
    // Line L = intersection of the two planes, through a point P0 on both.
    let d = cross.normalized();
    let v = d.crossed(&n); // perpendicular to d, in the circle's plane
    let denom = np.dot(&v);
    if denom.abs() <= 1e-30 {
        return;
    }
    let s = np.dot(&GpVec::from_pnts(&center, &q)) / denom;
    let p0 = center.translated_vec(&v.multiplied_scalar(s));
    let w = GpVec::from_pnts(&center, &p0);
    let a_coef = d.dot(&d);
    let b_coef = 2.0 * d.dot(&w);
    let c_coef = w.square_magnitude() - radius * radius;
    let disc = b_coef * b_coef - 4.0 * a_coef * c_coef;
    if disc < -etol * radius.max(1.0) {
        return;
    }
    let sq = if disc > 0.0 { disc.sqrt() } else { 0.0 };
    for t in [(-b_coef - sq) / (2.0 * a_coef), (-b_coef + sq) / (2.0 * a_coef)] {
        out.push(p0.translated_vec(&d.multiplied_scalar(t)));
    }
}

/// Intersect a general curve with a plane by sampled sign-change bisection.
pub(super) fn general_plane_hits(e: &Edge, pln: &GpPln, a: f64, b: f64, tol: f64, out: &mut Vec<GpPnt>) {
    let n = GpVec::from_xyz(pln.axis().direction().xyz());
    let q = pln.location();
    let (lo, hi) = (a.min(b), a.max(b));
    if (hi - lo).abs() <= 1e-15 {
        return;
    }
    let samples = 64;
    let step = (hi - lo) / samples as f64;
    let sdist = |u: f64| n.dot(&GpVec::from_pnts(&q, &eval_edge(e, u)));
    let mut prev_u = lo;
    let mut prev_s = sdist(lo);
    for i in 1..=samples {
        let u = lo + step * i as f64;
        let s = sdist(u);
        if prev_s * s < 0.0 {
            let (mut ua, mut ub) = (prev_u, u);
            let mut sa = prev_s;
            for _ in 0..48 {
                let um = 0.5 * (ua + ub);
                let sm = sdist(um);
                if sa * sm <= 0.0 {
                    ub = um;
                } else {
                    ua = um;
                    sa = sm;
                }
                if ub - ua < 1e-10 {
                    break;
                }
            }
            out.push(eval_edge(e, 0.5 * (ua + ub)));
        } else if s.abs() <= tol.max(1e-9) && prev_s.abs() > s.abs() {
            let g = |u: f64| sdist(u).abs();
            let (um, dm) = golden_1d(&g, prev_u, u, 1e-10);
            if dm <= tol.max(1e-9) {
                out.push(eval_edge(e, um));
            }
        }
        prev_u = u;
        prev_s = s;
    }
}

pub(super) fn dedupe_edge_face(mut hits: Vec<(f64, GpPnt)>, tol: f64) -> Vec<(f64, GpPnt)> {
    let etol = tol.max(1e-9);
    let mut out: Vec<(f64, GpPnt)> = Vec::new();
    for (u, p) in hits.drain(..) {
        let dup = out
            .iter()
            .any(|(ou, op)| (*ou - u).abs() <= etol && op.distance(&p) <= etol);
        if !dup {
            out.push((u, p));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Face–face intersection
// ---------------------------------------------------------------------------

/// Intersection of two faces as one or more 3D segments.
///
/// Planar faces produce the exact plane–plane line clipped to both boundary
/// polygons. Non-planar faces fall back to a sampled curve approximation.
pub fn face_face_intersection_segments(f1: &Face, f2: &Face, tol: f64) -> Vec<(GpPnt, GpPnt)> {
    let (p1, p2) = match (face_plane_from_face(f1), face_plane_from_face(f2)) {
        (Some(p1), Some(p2)) => (p1, p2),
        _ => return non_planar_face_face(f1, f2, tol),
    };
    let (origin, dir) = match plane_plane_intersection(&p1, &p2) {
        Some(x) => x,
        None => return Vec::new(), // parallel (coplanar included)
    };
    let poly1 = match face_polygon_2d_with_plane(f1, &p1) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let poly2 = match face_polygon_2d_with_plane(f2, &p2) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let iv1 = line_polygon_t_intervals(&p1, &origin, &dir, &poly1);
    let iv2 = line_polygon_t_intervals(&p2, &origin, &dir, &poly2);
    let mut out = Vec::new();
    for (t0, t1) in intersect_intervals(&iv1, &iv2, tol) {
        out.push((
            origin.translated_vec(&dir.multiplied_scalar(t0)),
            origin.translated_vec(&dir.multiplied_scalar(t1)),
        ));
    }
    out
}

/// Approximate face–face intersection for non-planar faces: a polyline of
/// short segments chaining the sampled near-coincident points.
pub(super) fn non_planar_face_face(f1: &Face, f2: &Face, tol: f64) -> Vec<(GpPnt, GpPnt)> {
    match face_face_intersection(f1, f2, tol) {
        FaceIntersect::Curve { mut points } => {
            // ponytail: the grid sampler returns an unordered cloud; chaining
            // by the dominant axis is a rough polyline. A marching-cubes or
            // subdivision contour tracer is the upgrade path.
            if points.len() < 2 {
                return Vec::new();
            }
            let mut min = points[0];
            let mut max = points[0];
            for p in &points {
                min = GpPnt::new(min.x().min(p.x()), min.y().min(p.y()), min.z().min(p.z()));
                max = GpPnt::new(max.x().max(p.x()), max.y().max(p.y()), max.z().max(p.z()));
            }
            let (ex, ey, ez) = (max.x() - min.x(), max.y() - min.y(), max.z() - min.z());
            if ex >= ey && ex >= ez {
                points.sort_by(|a, b| a.x().partial_cmp(&b.x()).unwrap_or(std::cmp::Ordering::Equal));
            } else if ey >= ez {
                points.sort_by(|a, b| a.y().partial_cmp(&b.y()).unwrap_or(std::cmp::Ordering::Equal));
            } else {
                points.sort_by(|a, b| a.z().partial_cmp(&b.z()).unwrap_or(std::cmp::Ordering::Equal));
            }
            points.windows(2).map(|w| (w[0], w[1])).collect()
        }
        _ => Vec::new(),
    }
}
