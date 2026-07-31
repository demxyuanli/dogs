//! Exact boolean foundation — edge/face intersection queries.
//!
//! Ports the `IntTools_EdgeEdge`, `IntTools_EdgeFace` and `IntTools_FaceFace`
//! algorithms (TKBO) at the level the exact boolean needs:
//!
//! - **edge–edge**: exact for line–line (segment intersection), line–circle
//!   (quadratic in the circle's plane) and coplanar circle–circle (radical
//!   line); all other combinations fall back to a sampling+refine solver.
//! - **edge–face**: exact for line–plane and circle–plane; general curves
//!   cross a plane by sampled sign-change bisection; non-planar faces use a
//!   curve–surface sampler. Results are filtered to the face's 2D boundary
//!   polygon.
//! - **face–face**: exact for planar faces (plane–plane line clipped to both
//!   boundary polygons); non-planar faces yield an approximate polyline.
//!
//! `bop_builder` consumes these exact signatures to split edges at their
//! intersection parameters and build the boolean result.

use occt_core::geom::polygon_ops::point_in_polygon2d;
use occt_core::gp::{GpPln, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::int::curve_curve::segment_segment_intersection_3d;
use occt_geom::geom_api;
use occt_geom::Curve;

use crate::brep_surface::face_plane;
use crate::brep_tool::BRepTool;
use crate::edge_split;
use crate::face_face::{plane_plane_intersection, FaceIntersect, face_face_intersection};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, vertex_position, wires_of_face};

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
fn eval_edge(e: &Edge, u: f64) -> GpPnt {
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
fn curve_to_edge_param(curve: &dyn Curve, e_first: f64, e_last: f64, u_curve: f64) -> f64 {
    if curve.first_parameter() == 0.0 && curve.last_parameter() == 1.0 {
        e_first + u_curve * (e_last - e_first)
    } else {
        u_curve
    }
}

/// Whether a curve is geometrically a straight line. Unbounded parameter
/// ranges belong only to `GeomLine` in this port; bounded curves are checked
/// for collinearity of eight samples against the first–last chord.
fn is_line_like(c: &dyn Curve) -> bool {
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
fn is_circle_like(c: &dyn Curve) -> bool {
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
fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
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
fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
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

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
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
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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
fn param_on_edge(e: &Edge, curve: &dyn Curve, p: &GpPnt, a: f64, b: f64, tol: f64) -> Option<f64> {
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
/// Line–line, line–circle and coplanar circle–circle are solved exactly; other
/// combinations fall back to the sampling solver in `geom_api`.
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
                points = geom_api::curve_curve_intersections(&*c1, &*c2, tol);
            }
        }
        (false, true) => {
            if ci1 {
                line_circle_hits(e2, &*c1, tol, &mut points);
            } else {
                points = geom_api::curve_curve_intersections(&*c1, &*c2, tol);
            }
        }
        (false, false) => {
            if ci1 && ci2 {
                match circle_circle_exact(&*c1, &*c2, tol) {
                    Some(pts) => points = pts,
                    None => points = geom_api::curve_curve_intersections(&*c1, &*c2, tol),
                }
            } else {
                points = geom_api::curve_curve_intersections(&*c1, &*c2, tol);
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
fn line_circle_hits(e_line: &Edge, c_circ: &dyn Curve, tol: f64, out: &mut Vec<GpPnt>) {
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
fn circle_circle_exact(c1: &dyn Curve, c2: &dyn Curve, tol: f64) -> Option<Vec<GpPnt>> {
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

fn dedupe_hits(mut hits: Vec<EdgeEdgeHit>, tol: f64) -> Vec<EdgeEdgeHit> {
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

// ---------------------------------------------------------------------------
// Edge–face intersection
// ---------------------------------------------------------------------------

/// All points where an edge meets a face: `(param_on_edge, point)`.
///
/// Planar faces use the exact line/circle–plane solves and keep only the hits
/// inside the face's 2D boundary polygon. Non-planar faces are handled by a
/// sampled curve–surface solver.
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
        let pts = geom_api::curve_surface_intersections(&*curve, &*surf, tol, 200);
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
fn line_plane_hits(
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
fn circle_plane_hits(
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
fn general_plane_hits(e: &Edge, pln: &GpPln, a: f64, b: f64, tol: f64, out: &mut Vec<GpPnt>) {
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

fn dedupe_edge_face(mut hits: Vec<(f64, GpPnt)>, tol: f64) -> Vec<(f64, GpPnt)> {
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
fn non_planar_face_face(f1: &Face, f2: &Face, tol: f64) -> Vec<(GpPnt, GpPnt)> {
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

/// Parameter intervals of the line `p0 + t·dir` that lie inside a 2D polygon
/// (in the plane's frame). `t` matches the 3D line's arc-length parameter.
fn line_polygon_t_intervals(pln: &GpPln, p0: &GpPnt, dir: &GpVec, poly: &[GpPnt2d]) -> Vec<(f64, f64)> {
    let n = poly.len();
    if n < 3 {
        return Vec::new();
    }
    let ax = pln.position();
    let xd = GpVec::from_xyz(ax.x_direction().xyz());
    let yd = GpVec::from_xyz(ax.y_direction().xyz());
    let v0 = GpVec::from_pnts(&pln.location(), p0);
    let p0_2 = GpPnt2d::new(v0.dot(&xd), v0.dot(&yd));
    let d_2 = GpVec2d::new(dir.dot(&xd), dir.dot(&yd));
    let (dx, dy) = (d_2.x(), d_2.y());
    let len2 = dx * dx + dy * dy;
    if len2 <= 1e-30 {
        return Vec::new();
    }
    let s = |p: &GpPnt2d| dx * (p.y() - p0_2.y()) - dy * (p.x() - p0_2.x());
    let proj = |p: &GpPnt2d| ((p.x() - p0_2.x()) * dx + (p.y() - p0_2.y()) * dy) / len2;
    let eps = 1e-9 * len2.sqrt().max(1.0);

    let mut crossings: Vec<f64> = Vec::new();
    let mut collinear: Vec<(f64, f64)> = Vec::new();
    for i in 0..n {
        let a = poly[i];
        let c = poly[(i + 1) % n];
        let sa = s(&a);
        let sc = s(&c);
        if sa.abs() <= eps && sc.abs() <= eps {
            let (ta, tc) = (proj(&a), proj(&c));
            collinear.push((ta.min(tc), ta.max(tc)));
        } else if sa * sc < 0.0 {
            let w = sa / (sa - sc);
            crossings.push(proj(&a) + (proj(&c) - proj(&a)) * w);
        }
    }
    crossings.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|x, y| (*x - *y).abs() <= eps);

    let mut intervals: Vec<(f64, f64)> = collinear;
    for w in crossings.windows(2) {
        let (t0, t1) = (w[0], w[1]);
        if t1 - t0 <= eps {
            continue;
        }
        let mid = GpPnt2d::new(p0_2.x() + dx * 0.5 * (t0 + t1), p0_2.y() + dy * 0.5 * (t0 + t1));
        if point_in_polygon2d(poly, &mid) {
            intervals.push((t0, t1));
        }
    }
    merge_intervals(intervals)
}

fn merge_intervals(mut v: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    v.retain(|(a, b)| b - a > 1e-12);
    if v.is_empty() {
        return v;
    }
    v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut out = vec![v[0]];
    for &(a, b) in v.iter().skip(1) {
        let last = out.last_mut().unwrap();
        if a <= last.1 + 1e-9 {
            last.1 = last.1.max(b);
        } else {
            out.push((a, b));
        }
    }
    out
}

fn intersect_intervals(a: &[(f64, f64)], b: &[(f64, f64)], tol: f64) -> Vec<(f64, f64)> {
    let eps = tol.max(1e-9);
    let mut out = Vec::new();
    for &(a0, a1) in a {
        for &(b0, b1) in b {
            let lo = a0.max(b0);
            let hi = a1.min(b1);
            if hi - lo > eps {
                out.push((lo, hi));
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Face classification helpers
// ---------------------------------------------------------------------------

/// The plane underlying a planar face, if it is planar.
pub fn face_plane_from_face(f: &Face) -> Option<GpPln> {
    face_plane(f)
}

/// Whether `p` is on the face: within `tol` of the surface and inside its
/// 2D boundary polygon (planar faces), or within `tol` of the surface
/// (non-planar, boundary not checked).
pub fn point_on_face(f: &Face, p: &GpPnt, tol: f64) -> bool {
    if let Some(pln) = face_plane_from_face(f) {
        let n = GpVec::from_xyz(pln.axis().direction().xyz());
        let dist = n.dot(&GpVec::from_pnts(&pln.location(), p)).abs();
        if dist > tol {
            return false;
        }
        match face_polygon_2d(f) {
            Some(poly) => point_in_polygon2d(&poly, &project_to_2d(&pln, p)),
            None => true, // unbounded face
        }
    } else {
        let Some(surf) = BRepTool::face_surface(f) else {
            return false;
        };
        geom_api::project_point_on_surface(&*surf, p, tol).map_or(false, |ps| ps.distance <= tol)
    }
}

/// The outer-wire vertex positions of a planar face projected onto the face
/// plane as 2D points.
pub fn face_polygon_2d(f: &Face) -> Option<Vec<GpPnt2d>> {
    let pln = face_plane_from_face(f)?;
    face_polygon_2d_with_plane(f, &pln)
}

/// [`face_polygon_2d`] using an explicitly supplied plane (must match the one
/// used for projections so the frames agree).
///
/// The wire edges are chained geometrically end-to-end rather than trusting
/// each edge's stored vertex order: a shared edge may be stored in either
/// direction relative to the wire traversal (OCCT uses an orientation flag,
/// which this port does not model on the edge).
fn face_polygon_2d_with_plane(f: &Face, pln: &GpPln) -> Option<Vec<GpPnt2d>> {
    let wire = wires_of_face(f).into_iter().next()?;
    let edges = edges_of_wire(&wire);
    if edges.len() < 3 {
        return None;
    }
    let mut segs: Vec<(GpPnt, GpPnt)> = Vec::with_capacity(edges.len());
    for e in edges {
        let (a, b) = edge_vertices(&e);
        segs.push((vertex_position(&a?), vertex_position(&b?)));
    }
    let mut pts: Vec<GpPnt> = Vec::with_capacity(segs.len() + 1);
    pts.push(segs[0].0);
    pts.push(segs[0].1);
    let mut used = vec![false; segs.len()];
    used[0] = true;
    for _ in 1..segs.len() {
        let tail = *pts.last().unwrap();
        let mut extended = false;
        for (i, (a, b)) in segs.iter().enumerate() {
            if used[i] {
                continue;
            }
            if a.distance(&tail) <= 1e-9 {
                pts.push(*b);
                used[i] = true;
                extended = true;
                break;
            }
            if b.distance(&tail) <= 1e-9 {
                pts.push(*a);
                used[i] = true;
                extended = true;
                break;
            }
        }
        if !extended {
            return None; // wire does not chain end-to-end
        }
    }
    // Drop the closing duplicate of a closed wire.
    if pts.len() > 1 && pts[0].distance(&pts[pts.len() - 1]) <= 1e-9 {
        pts.pop();
    }
    if pts.len() < 3 {
        return None;
    }
    Some(pts.iter().map(|p| project_to_2d(pln, p)).collect())
}

/// Project a 3D point onto a plane's 2D frame.
fn project_to_2d(pln: &GpPln, p: &GpPnt) -> GpPnt2d {
    let ax = pln.position();
    let xd = GpVec::from_xyz(ax.x_direction().xyz());
    let yd = GpVec::from_xyz(ax.y_direction().xyz());
    let v = GpVec::from_pnts(&pln.location(), p);
    GpPnt2d::new(v.dot(&xd), v.dot(&yd))
}

// ---------------------------------------------------------------------------
// Edge splitting / batch collection
// ---------------------------------------------------------------------------

/// Split an edge at the given interior parameters (`split_edge` wrapper).
pub fn split_edges_at_params(e: &Edge, params: &[f64]) -> Vec<Edge> {
    edge_split::split_edge(e, params)
}

/// For every edge index in `edges`, every intersection point with an edge of
/// the solids in `shapes`: `(edge_index, param_on_edge, point)`.
pub fn collect_edge_intersection_points(
    edges: &[Edge],
    shapes: &[TopoShape],
    tol: f64,
) -> Vec<(usize, f64, GpPnt)> {
    let mut out = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        for s in shapes {
            for oe in edges_of(s) {
                if oe.0.same_tshape(&e.0) {
                    continue;
                }
                for hit in edge_edge_intersections(e, &oe, tol) {
                    out.push((i, hit.u1, hit.point));
                }
            }
        }
    }
    let etol = tol.max(1e-9);
    out.sort_by(|a, b| {
        a.0.cmp(&b.0).then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    out.dedup_by(|a, b| a.0 == b.0 && (a.1 - b.1).abs() <= etol && a.2.distance(&b.2) <= etol);
    out
}

// ---------------------------------------------------------------------------
// Small numeric helpers
// ---------------------------------------------------------------------------

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    const GOLD: f64 = 0.6180339887498949;
    let mut a = lo;
    let mut b = hi;
    let mut c = b - GOLD * (b - a);
    let mut d = a + GOLD * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    while (b - a) > eps {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - GOLD * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GOLD * (b - a);
            fd = f(d);
        }
    }
    let x = 0.5 * (a + b);
    (x, f(x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::geom::polygon_ops::polygon_area2d;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir};
    use occt_geom::{GeomCircle, GeomPlane};

    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::shape::{Compound, TopoShape};
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::faces_of;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    /// Planar square face `origin + [0,size]·u_dir + [0,size]·v_dir`.
    fn make_square_face(b: &TopoBuilder, origin: GpPnt, u_dir: GpDir, v_dir: GpDir, size: f64) -> Face {
        let uv = GpVec::from_xyz(u_dir.xyz()).multiplied_scalar(size);
        let vv = GpVec::from_xyz(v_dir.xyz()).multiplied_scalar(size);
        let c0 = origin;
        let c1 = c0.translated_vec(&uv);
        let c2 = c1.translated_vec(&vv);
        let c3 = c0.translated_vec(&vv);
        let e1 = b.make_edge_segment(&c0, &c1);
        let e2 = b.make_edge_segment(&c1, &c2);
        let e3 = b.make_edge_segment(&c2, &c3);
        let e4 = b.make_edge_segment(&c3, &c0);
        let wire = b.make_wire(&[e1, e2, e3, e4]);
        let normal = GpVec::from_xyz(u_dir.xyz()).crossed(&GpVec::from_xyz(v_dir.xyz()));
        let n = GpDir::from_vec(&normal).expect("normal");
        let ax3 = GpAx3::new(origin, n, &u_dir).expect("frame");
        b.make_face(Arc::new(GeomPlane::new(GpPln::new(ax3))), &[wire])
    }

    #[test]
    fn crossing_line_edges_single_hit() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let hits = edge_edge_intersections(&e1, &e2, 1e-9);
        assert_eq!(hits.len(), 1, "hits: {hits:?}");
        assert!(hits[0].point.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6);
        let d = 2.0f64.sqrt();
        assert!((hits[0].u1 - d).abs() < 1e-6, "u1={}", hits[0].u1);
        assert!((hits[0].u2 - d).abs() < 1e-6, "u2={}", hits[0].u2);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_edge_through_planar_face() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        let hits = edge_face_intersections(&e, &face, 1e-9);
        assert_eq!(hits.len(), 1, "hits: {hits:?}");
        assert!(hits[0].1.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6);
        assert!((hits[0].0 - 1.0).abs() < 1e-6, "u={}", hits[0].0);
        clear_tree(&face.0);
        clear_tree(&e.0);
    }

    #[test]
    fn line_edge_missing_face() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let e = b.make_edge_segment(&GpPnt::new(5.0, 5.0, -1.0), &GpPnt::new(5.0, 5.0, 1.0));
        assert!(edge_face_intersections(&e, &face, 1e-9).is_empty());
        clear_tree(&face.0);
        clear_tree(&e.0);
    }

    #[test]
    fn perpendicular_faces_intersection_segment() {
        let b = TopoBuilder::new();
        let f1 = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let f2 = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), 1.0);
        let segs = face_face_intersection_segments(&f1, &f2, 1e-9);
        assert_eq!(segs.len(), 1, "segs: {segs:?}");
        let (p, q) = segs[0];
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, 0.0, 0.0);
        assert!((p.distance(&a) < 1e-6 && q.distance(&c) < 1e-6) || (p.distance(&c) < 1e-6 && q.distance(&a) < 1e-6));
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn parallel_faces_no_segment() {
        let b = TopoBuilder::new();
        let f1 = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let f2 = make_square_face(&b, GpPnt::new(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        assert!(face_face_intersection_segments(&f1, &f2, 1e-9).is_empty());
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn point_on_face_inside_outside() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        assert!(point_on_face(&face, &GpPnt::new(0.5, 0.5, 0.0), 1e-9));
        assert!(!point_on_face(&face, &GpPnt::new(2.0, 0.5, 0.0), 1e-9));
        assert!(!point_on_face(&face, &GpPnt::new(0.5, 0.5, 1.0), 1e-9));
        clear_tree(&face.0);
    }

    #[test]
    fn face_polygon_2d_unit_square() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let poly = face_polygon_2d(&face).expect("polygon");
        assert_eq!(poly.len(), 4);
        assert!((polygon_area2d(&poly).abs() - 1.0).abs() < 1e-9, "area={}", polygon_area2d(&poly));
        clear_tree(&face.0);
    }

    #[test]
    fn split_edge_at_params_wrapper() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
        let parts = split_edges_at_params(&e, &[3.0, 7.0]);
        assert_eq!(parts.len(), 3);
        let ranges: Vec<(f64, f64)> = parts.iter().map(|p| BRepTool::edge_parameters(p)).collect();
        assert_eq!(ranges, vec![(0.0, 3.0), (3.0, 7.0), (7.0, 10.0)]);
        for p in &parts {
            clear_tree(&p.0);
        }
        clear_tree(&e.0);
    }

    #[test]
    fn line_circle_edges_coplanar_intersect() {
        let b = TopoBuilder::new();
        let ce = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let le = b.make_edge_segment(&GpPnt::new(-1.5, 0.0, 0.0), &GpPnt::new(1.5, 0.0, 0.0));
        let hits = edge_edge_intersections(&le, &ce, 1e-6);
        assert_eq!(hits.len(), 2, "hits: {hits:?}");
        for h in &hits {
            assert!((h.point.coord.modulus() - 1.0).abs() < 1e-6, "point {:?}", h.point);
        }
        clear_tree(&ce.0);
        clear_tree(&le.0);
    }

    #[test]
    fn circle_circle_edges_coplanar_intersect() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let hits = edge_edge_intersections(&c1, &c2, 1e-6);
        assert_eq!(hits.len(), 2, "hits: {hits:?}");
        for h in &hits {
            assert!((h.point.x() - 0.5).abs() < 1e-6, "x={}", h.point.x());
            assert!((h.point.y().abs() - 0.75f64.sqrt()).abs() < 1e-4, "y={}", h.point.y());
        }
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    fn find_face_by_normal(bx: &crate::primitives::BRepPrimBox, n: (f64, f64, f64)) -> Face {
        faces_of(&bx.solid.0)
            .into_iter()
            .find(|f| {
                face_plane_from_face(f).map_or(false, |pln| {
                    let d = GpVec::from_xyz(pln.axis().direction().xyz());
                    (d.x() - n.0).abs() < 1e-9 && (d.y() - n.1).abs() < 1e-9 && (d.z() - n.2).abs() < 1e-9
                })
            })
            .expect("box face with normal")
    }

    #[test]
    fn box_face_polygon_chains_shared_edges() {
        let bx = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&bx.solid.0);
        assert_eq!(faces.len(), 6);
        for f in faces {
            let poly = face_polygon_2d(&f).expect("box face polygon");
            assert_eq!(poly.len(), 4, "poly: {poly:?}");
            assert!((polygon_area2d(&poly).abs() - 1.0).abs() < 1e-9, "area={}", polygon_area2d(&poly));
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn adjacent_box_faces_intersect_along_edge() {
        let bx = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let fx = find_face_by_normal(&bx, (1.0, 0.0, 0.0));
        let fy = find_face_by_normal(&bx, (0.0, 1.0, 0.0));
        let segs = face_face_intersection_segments(&fx, &fy, 1e-9);
        assert_eq!(segs.len(), 1, "segs: {segs:?}");
        let (p, q) = segs[0];
        assert!((p.x() - 1.0).abs() < 1e-6 && (p.y() - 1.0).abs() < 1e-6);
        assert!((q.x() - 1.0).abs() < 1e-6 && (q.y() - 1.0).abs() < 1e-6);
        let zs = [p.z(), q.z()];
        assert!(zs.contains(&0.0) && zs.contains(&1.0), "z endpoints {zs:?}");
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn collect_edge_intersection_points_across_shapes() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let mut comp = Compound::new();
        b.add_compound(&mut comp, &e2.0);
        let res = collect_edge_intersection_points(&[e1.clone()], &[comp.0.clone()], 1e-9);
        assert_eq!(res.len(), 1, "res: {res:?}");
        assert_eq!(res[0].0, 0);
        assert!(res[0].2.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
        clear_tree(&comp.0);
    }
}
