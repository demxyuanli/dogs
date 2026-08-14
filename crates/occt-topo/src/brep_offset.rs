//! Parallel offsets of curves, faces, shells and solids.
//! Source: `BRepOffsetAPI_MakeOffsetShape`, `BRepOffsetAPI_MakeOffset`,
//! `BRepOffset_Offset`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpVec, GpXyz};
use occt_geom::surface_fit::fit_plane;
use occt_geom::{Curve, GeomCircle, GeomCylinder, GeomPlane, GeomSphere, Surface};

use crate::abs::ShapeType;
use crate::brep_builder_api::make_edge_arc;
use crate::brep_surface::{face_is_planar, face_plane, is_planar, sphere_center, surface_normal};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::shape_ops::translated_copy;
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// 3×3 determinant (Cramer's-rule helper).
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system via Cramer's rule.
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

/// Circumcenter of three non-collinear 3D points, if it exists.
fn circumcenter3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.modulus() < 1e-20 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.x, n.y, n.z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n)];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

/// Intersection of three planes: solve `n_i · x = d_i` (Cramer's rule).
/// Returns `None` when any two planes are parallel (singular system).
pub fn plane_plane_plane_intersection(p1: &GpPln, p2: &GpPln, p3: &GpPln) -> Option<GpPnt> {
    let n1 = *p1.axis().direction().xyz();
    let n2 = *p2.axis().direction().xyz();
    let n3 = *p3.axis().direction().xyz();
    let d1 = n1.dot(&p1.location().coord);
    let d2 = n2.dot(&p2.location().coord);
    let d3 = n3.dot(&p3.location().coord);
    let mat = [[n1.x, n1.y, n1.z], [n2.x, n2.y, n2.z], [n3.x, n3.y, n3.z]];
    let rhs = [d1, d2, d3];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

// ---------------------------------------------------------------------------
// 1. Curve offset (2D, in a plane)
// ---------------------------------------------------------------------------

/// Miter offset of a closed polygon in `plane`: every edge is shifted
/// perpendicular in-plane by `distance` and consecutive offset edges are
/// re-intersected, producing sharp (miter) corners. Positive `distance`
/// expands a CCW polygon (relative to the plane normal); negative shrinks it.
pub fn offset_polygon(points: &[GpPnt], plane: &GpPln, distance: f64) -> Result<Vec<GpPnt>, String> {
    let n = points.len();
    if n < 3 {
        return Err("offset_polygon: need at least 3 points".into());
    }
    let nxyz = *plane.axis().direction().xyz();
    let plane_n = GpVec::from_xyz(&nxyz);
    // Signed area relative to the plane normal decides which perpendicular
    // direction is "outward" (CCW polygon → interior on the left).
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += points[i].coord.crossed(&points[j].coord).dot(&nxyz);
    }
    let side = if area >= 0.0 { 1.0 } else { -1.0 };

    let edge_dir = |a: &GpPnt, b: &GpPnt| -> Option<GpVec> {
        let e = GpVec::from_pnts(a, b).normalized();
        if e.xyz().square_modulus() < 1e-24 {
            None
        } else {
            Some(e)
        }
    };
    let outward = |a: &GpPnt, b: &GpPnt| -> Option<GpVec> {
        let e = edge_dir(a, b)?;
        let mut o = e.xyz().crossed(&nxyz);
        if side < 0.0 {
            o = o.multiplied(-1.0);
        }
        let om = o.modulus();
        if om < 1e-30 {
            None
        } else {
            Some(GpVec::from_xyz(&o.divided(om)))
        }
    };

    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let a = &points[i];
        let b = &points[(i + 1) % n];
        let j = (i + n - 1) % n;
        let a_prev = &points[j];
        let b_prev = &points[i];
        let o_cur = outward(a, b).ok_or("offset_polygon: degenerate edge")?;
        let o_prev = outward(a_prev, b_prev).ok_or("offset_polygon: degenerate edge")?;
        let e_cur = edge_dir(a, b).ok_or("offset_polygon: degenerate edge")?;
        let e_prev = edge_dir(a_prev, b_prev).ok_or("offset_polygon: degenerate edge")?;
        let a1 = a.translated_vec(&o_cur.multiplied_scalar(distance));
        let a2 = a_prev.translated_vec(&o_prev.multiplied_scalar(distance));
        let pt = line_line_intersect_3d(&a1, &e_cur, &a2, &e_prev, &plane_n)
            .ok_or("offset_polygon: parallel offset edges (miter fails)")?;
        result.push(pt);
    }
    Ok(result)
}

/// Intersection of two coplanar lines `a1 + t·d1` and `a2 + s·d2`.
fn line_line_intersect_3d(a1: &GpPnt, d1: &GpVec, a2: &GpPnt, d2: &GpVec, n: &GpVec) -> Option<GpPnt> {
    let num = GpVec::from_pnts(a1, a2).xyz().crossed(d2.xyz()).dot(n.xyz());
    let den = d1.xyz().crossed(d2.xyz()).dot(n.xyz());
    if den.abs() < 1e-12 {
        return None;
    }
    let t = num / den;
    Some(a1.translated_vec(&d1.multiplied_scalar(t)))
}

/// Offset of a circle in `plane`: radius ± distance (same center).
pub fn offset_circle(center: GpPnt, radius: f64, plane: &GpPln, distance: f64) -> Result<GpCirc, String> {
    let new_r = radius + distance;
    if new_r <= 0.0 {
        return Err("offset_circle: radius would become non-positive".into());
    }
    let ax = plane.position();
    let ax2 = GpAx2::new(center, ax.direction(), *ax.x_direction()).map_err(|e| e.to_string())?;
    Ok(GpCirc::new(ax2, new_r))
}

/// Project a 3D point onto the plane's `(u, v)` frame.
fn to_2d(plane: &GpPln, p: &GpPnt) -> (f64, f64) {
    let ax = plane.position();
    let d = GpVec::from_pnts(&ax.location(), p);
    to_2d_vec(plane, &d)
}

/// Project a 3D vector onto the plane's in-plane axes.
fn to_2d_vec(plane: &GpPln, v: &GpVec) -> (f64, f64) {
    let ax = plane.position();
    (v.xyz().dot(ax.x_direction().xyz()), v.xyz().dot(ax.y_direction().xyz()))
}

/// Reconstruct a 3D point from plane-frame coordinates.
fn to_3d(plane: &GpPln, u: f64, v: f64) -> GpPnt {
    let ax = plane.position();
    let o = ax.location();
    let x = GpVec::from_xyz(ax.x_direction().xyz());
    let y = GpVec::from_xyz(ax.y_direction().xyz());
    o.translated_vec(&x.multiplied_scalar(u)).translated_vec(&y.multiplied_scalar(v))
}

/// Intersection of two infinite 2D lines.
fn line_line_2d(a1: (f64, f64), d1: (f64, f64), a2: (f64, f64), d2: (f64, f64)) -> Option<(f64, f64)> {
    let den = d1.0 * d2.1 - d1.1 * d2.0;
    if den.abs() < 1e-12 {
        return None;
    }
    let t = ((a2.0 - a1.0) * d2.1 - (a2.1 - a1.1) * d2.0) / den;
    Some((a1.0 + t * d1.0, a1.1 + t * d1.1))
}

/// Intersection of a 2D line and a 2D circle (0, 1 or 2 points).
fn line_circle_2d(a: (f64, f64), d: (f64, f64), c: (f64, f64), r: f64) -> Vec<(f64, f64)> {
    let dx = a.0 - c.0;
    let dy = a.1 - c.1;
    let qa = d.0 * d.0 + d.1 * d.1;
    let qb = 2.0 * (dx * d.0 + dy * d.1);
    let qc = dx * dx + dy * dy - r * r;
    let disc = qb * qb - 4.0 * qa * qc;
    if disc < -1e-9 {
        return Vec::new();
    }
    if disc < 0.0 {
        let t = -qb / (2.0 * qa);
        return vec![(a.0 + t * d.0, a.1 + t * d.1)];
    }
    let s = disc.sqrt();
    let two_a = 2.0 * qa;
    vec![
        (a.0 + (-qb + s) * d.0 / two_a, a.1 + (-qb + s) * d.1 / two_a),
        (a.0 + (-qb - s) * d.0 / two_a, a.1 + (-qb - s) * d.1 / two_a),
    ]
}

/// Intersection of two 2D circles (0, 1 or 2 points).
fn circle_circle_2d(c1: (f64, f64), r1: f64, c2: (f64, f64), r2: f64) -> Vec<(f64, f64)> {
    let dx = c2.0 - c1.0;
    let dy = c2.1 - c1.1;
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-12 {
        return Vec::new();
    }
    let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let h2 = r1 * r1 - a * a;
    if h2 < -1e-9 {
        return Vec::new();
    }
    let h = if h2 < 0.0 { 0.0 } else { h2.sqrt() };
    let mx = c1.0 + a * dx / d;
    let my = c1.1 + a * dy / d;
    let px = -dy / d;
    let py = dx / d;
    vec![(mx + h * px, my + h * py), (mx - h * px, my - h * py)]
}

/// Geometric kind of an edge's curve, detected by sampling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OffsetCurveKind {
    Line,
    Circle,
}

/// Per-edge classification for wire offsetting.
struct EdgeInfo {
    kind: OffsetCurveKind,
    curve: Arc<dyn Curve>,
    first: f64,
    last: f64,
    center: GpPnt,
    radius: f64,
}

/// Line parameters from samples (a point and a unit direction), if collinear.
fn line_params(pts: &[GpPnt]) -> Option<(GpPnt, GpVec)> {
    let mut d0: Option<GpVec> = None;
    for i in 1..pts.len() {
        let d = GpVec::from_pnts(&pts[0], &pts[i]);
        if d.xyz().square_modulus() < 1e-24 {
            continue;
        }
        d0 = Some(d.normalized());
        break;
    }
    let d0 = d0?;
    for p in pts {
        let d = GpVec::from_pnts(&pts[0], p);
        if d.xyz().crossed(d0.xyz()).modulus() > 1e-6 * d.xyz().modulus().max(1e-30) {
            return None;
        }
    }
    Some((pts[0], d0))
}

/// Circle parameters from samples (center and radius), if circular.
fn circle_params(pts: &[GpPnt]) -> Option<(GpPnt, f64)> {
    let n = pts.len();
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let v1 = GpVec::from_pnts(&pts[i], &pts[j]);
                let v2 = GpVec::from_pnts(&pts[i], &pts[k]);
                if v1.xyz().crossed(v2.xyz()).modulus() > 1e-12 {
                    let c = circumcenter3(&pts[i], &pts[j], &pts[k])?;
                    let r = c.distance(&pts[i]);
                    if pts.iter().all(|p| (p.distance(&c) - r).abs() < 1e-6 * r.max(1.0)) {
                        return Some((c, r));
                    }
                    return None;
                }
            }
        }
    }
    None
}

/// Classify an edge's curve as a line segment or a circular arc.
fn classify_edge(e: &Edge) -> Result<EdgeInfo, String> {
    let curve = BRepTool::edge_curve(e).ok_or("offset_wire_2d: edge has no curve")?;
    let (first, last) = BRepTool::edge_parameters(e);
    if !first.is_finite() || !last.is_finite() || last <= first {
        return Err("offset_wire_2d: edge has unbounded or empty parameter range".into());
    }
    let ns = 9;
    let mut samples = Vec::with_capacity(ns);
    for i in 0..ns {
        samples.push(curve.d0(first + (last - first) * i as f64 / (ns - 1) as f64));
    }
    if line_params(&samples).is_some() {
        return Ok(EdgeInfo { kind: OffsetCurveKind::Line, curve, first, last, center: GpPnt::zero(), radius: 0.0 });
    }
    if let Some((c, r)) = circle_params(&samples) {
        return Ok(EdgeInfo { kind: OffsetCurveKind::Circle, curve, first, last, center: c, radius: r });
    }
    Err("offset_wire_2d: unsupported curve kind (must be a line or a circular arc)".into())
}

/// An offset edge curve, in plane-frame coordinates.
enum OffsetCurve2 {
    Line { a: (f64, f64), d: (f64, f64) },
    Circle { c: (f64, f64), r: f64 },
}

/// Build the offset of a single edge. Returns the offset curve plus the offset
/// arc midpoint (for circle edges) used when rebuilding the arc.
fn offset_curve_for_edge(
    info: &EdgeInfo,
    a: &GpPnt,
    b: &GpPnt,
    plane: &GpPln,
    side: f64,
    distance: f64,
) -> Result<(OffsetCurve2, Option<GpPnt>), String> {
    let nrm = GpVec::from_xyz(plane.axis().direction().xyz());
    match info.kind {
        OffsetCurveKind::Line => {
            let p1 = info.curve.d0(info.first);
            let p2 = info.curve.d0(info.last);
            let e = GpVec::from_pnts(&p1, &p2).normalized();
            if e.xyz().square_modulus() < 1e-24 {
                return Err("offset_wire_2d: degenerate line edge".into());
            }
            let mut out = e.xyz().crossed(&nrm.xyz());
            if side < 0.0 {
                out = out.multiplied(-1.0);
            }
            let om = out.modulus();
            if om < 1e-30 {
                return Err("offset_wire_2d: degenerate offset direction".into());
            }
            let out = GpVec::from_xyz(&out.divided(om));
            let anchor = p1.translated_vec(&out.multiplied_scalar(distance));
            Ok((OffsetCurve2::Line { a: to_2d(plane, &anchor), d: to_2d_vec(plane, &e) }, None))
        }
        OffsetCurveKind::Circle => {
            let mid_u = 0.5 * (info.first + info.last);
            let m = info.curve.d0(mid_u);
            let rdir = GpVec::from_pnts(&info.center, &m).normalized();
            let t = GpVec::from_pnts(a, b).normalized();
            let mut outward = t.xyz().crossed(&nrm.xyz());
            if side < 0.0 {
                outward = outward.multiplied(-1.0);
            }
            let om = outward.modulus();
            if om < 1e-30 {
                return Err("offset_wire_2d: degenerate offset direction".into());
            }
            let outward = GpVec::from_xyz(&outward.divided(om));
            // A bulge pointing outward grows the radius; an inward notch shrinks it.
            let sign = if rdir.xyz().dot(&outward.xyz()) >= 0.0 { 1.0 } else { -1.0 };
            let new_r = info.radius + distance * sign;
            if new_r <= 0.0 {
                return Err("offset_wire_2d: arc radius would become non-positive".into());
            }
            let mid_off = m.translated_vec(&rdir.multiplied_scalar(distance * sign));
            Ok((OffsetCurve2::Circle { c: to_2d(plane, &info.center), r: new_r }, Some(mid_off)))
        }
    }
}

/// Intersection of two offset curves, choosing the branch nearest `orig`.
fn intersect_offset_curves(a: &OffsetCurve2, b: &OffsetCurve2, orig: &GpPnt, plane: &GpPln) -> Option<GpPnt> {
    let candidates: Vec<(f64, f64)> = match (a, b) {
        (OffsetCurve2::Line { a: a1, d: d1 }, OffsetCurve2::Line { a: a2, d: d2 }) => {
            line_line_2d(*a1, *d1, *a2, *d2).into_iter().collect()
        }
        (OffsetCurve2::Line { a: la, d: ld }, OffsetCurve2::Circle { c, r }) => line_circle_2d(*la, *ld, *c, *r),
        (OffsetCurve2::Circle { c, r }, OffsetCurve2::Line { a: la, d: ld }) => line_circle_2d(*la, *ld, *c, *r),
        (OffsetCurve2::Circle { c: c1, r: r1 }, OffsetCurve2::Circle { c: c2, r: r2 }) => circle_circle_2d(*c1, *r1, *c2, *r2),
    };
    let o = to_2d(plane, orig);
    let d2 = |p: &(f64, f64)| {
        let dx = p.0 - o.0;
        let dy = p.1 - o.1;
        dx * dx + dy * dy
    };
    candidates
        .into_iter()
        .min_by(|p, q| d2(p).partial_cmp(&d2(q)).unwrap_or(std::cmp::Ordering::Equal))
        .map(|p| to_3d(plane, p.0, p.1))
}

/// Rebuild a wire that is a single full circle at the new radius.
fn rebuild_circle_wire(circ: &GpCirc, first: f64, last: f64) -> Result<Wire, String> {
    let b = TopoBuilder::new();
    let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ.clone()));
    let mut e = b.make_edge(curve.clone(), first, last);
    let seam = curve.d0(first);
    let v = b.make_vertex(seam, 0.0);
    b.add(&mut e.0, &v.0);
    b.add(&mut e.0, &v.0);
    let w = b.make_wire(&[e]);
    w.set_closed(true);
    Ok(w)
}

/// Offset a planar wire in its plane: line edges shift perpendicular,
/// circular-arc edges grow/shrink their radius, and consecutive offset edges
/// are re-intersected to rebuild the wire with sharp miter corners.
pub fn offset_wire_2d(wire: &Wire, plane: &GpPln, distance: f64) -> Result<Wire, String> {
    let edges = edges_of_wire(wire);
    if edges.is_empty() {
        return Err("offset_wire_2d: wire has no edges".into());
    }
    let infos: Vec<EdgeInfo> = edges.iter().map(classify_edge).collect::<Result<_, _>>()?;

    // A single closed circle offsets to a single circle of the new radius.
    if edges.len() == 1 && infos[0].kind == OffsetCurveKind::Circle {
        let info = &infos[0];
        let circ = offset_circle(info.center, info.radius, plane, distance)?;
        return rebuild_circle_wire(&circ, info.first, info.last);
    }

    let n = edges.len();
    let pts: Vec<GpPnt> = infos.iter().map(|i| i.curve.d0(i.first)).collect();
    let nxyz = *plane.axis().direction().xyz();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += pts[i].coord.crossed(&pts[j].coord).dot(&nxyz);
    }
    let side = if area >= 0.0 { 1.0 } else { -1.0 };

    let mut curves = Vec::with_capacity(n);
    let mut arc_mids = Vec::with_capacity(n);
    for i in 0..n {
        let (c, mid) = offset_curve_for_edge(&infos[i], &pts[i], &pts[(i + 1) % n], plane, side, distance)?;
        curves.push(c);
        arc_mids.push(mid);
    }

    let mut new_pts = Vec::with_capacity(n);
    for i in 0..n {
        let prev = (i + n - 1) % n;
        let p = intersect_offset_curves(&curves[prev], &curves[i], &pts[i], plane)
            .ok_or("offset_wire_2d: consecutive offset edges do not intersect")?;
        new_pts.push(p);
    }

    let b = TopoBuilder::new();
    let mut new_edges: Vec<Edge> = Vec::with_capacity(n);
    for i in 0..n {
        let p1 = new_pts[i];
        let p2 = new_pts[(i + 1) % n];
        if p1.distance(&p2) < 1e-9 {
            continue;
        }
        match infos[i].kind {
            OffsetCurveKind::Line => new_edges.push(b.make_edge_segment(&p1, &p2)),
            OffsetCurveKind::Circle => {
                let mid = arc_mids[i].ok_or("offset_wire_2d: missing arc midpoint")?;
                new_edges.push(make_edge_arc(&p1, &mid, &p2)?);
            }
        }
    }
    if new_edges.len() < 2 {
        return Err("offset_wire_2d: degenerate offset wire".into());
    }
    let w = b.make_wire(&new_edges);
    w.set_closed(true);
    Ok(w)
}

/// Offset a planar 3D wire by sampling its points, fitting the plane via
/// least squares, then delegating to [`offset_wire_2d`].
pub fn offset_wire_3d(wire: &Wire, distance: f64) -> Result<Wire, String> {
    let mut pts = Vec::new();
    for e in edges_of_wire(wire) {
        if let Some(c) = BRepTool::edge_curve(&e) {
            let (a, b) = BRepTool::edge_parameters(&e);
            if a.is_finite() && b.is_finite() && b > a {
                for i in 0..4 {
                    pts.push(c.d0(a + (b - a) * i as f64 / 3.0));
                }
            }
        }
    }
    if pts.len() < 3 {
        return Err("offset_wire_3d: wire has too few sample points".into());
    }
    let pln = fit_plane(&pts).ok_or("offset_wire_3d: wire is not planar")?;
    offset_wire_2d(wire, &pln, distance)
}

// ---------------------------------------------------------------------------
// 2. Face offset (parallel surface)
// ---------------------------------------------------------------------------

/// Distance from a point to an axis line.
fn dist_point_axis(p: &GpPnt, loc: &GpPnt, dir: &GpDir) -> f64 {
    let v = GpVec::from_pnts(loc, p);
    let d = dir.xyz();
    let proj = d.multiplied(v.xyz().dot(d));
    GpVec::from_xyz(&v.xyz().subtracted(&proj)).xyz().modulus()
}

/// Extract cylinder parameters (a point on the axis, the axis direction, the
/// radius) by fitting two circles to rings at different heights.
fn extract_cylinder(s: &dyn Surface) -> Option<(GpPnt, GpDir, f64)> {
    let (u0, u1, v0, v1) = bounds(s);
    let va = v0 + 0.25 * (v1 - v0);
    let vb = v0 + 0.75 * (v1 - v0);
    let us = [0.0, 0.25, 0.5];
    let ring = |v: f64| -> Vec<GpPnt> {
        us.iter().map(|&f| s.d0(u0 + (u1 - u0) * f, v)).collect()
    };
    let pa = ring(va);
    let pb = ring(vb);
    let ca = circumcenter3(&pa[0], &pa[1], &pa[2])?;
    let cb = circumcenter3(&pb[0], &pb[1], &pb[2])?;
    let dir = GpDir::from_vec(&GpVec::from_pnts(&ca, &cb)).ok()?;
    let r = ca.distance(&pa[0]);
    Some((ca, dir, r))
}

/// Whether every sampled point of `s` is `r` from `center` (spherical test).
fn is_spherical(s: &dyn Surface, center: &GpPnt, r: f64) -> bool {
    let (u0, u1, v0, v1) = bounds(s);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            if (s.d0(u, v).distance(center) - r).abs() > 1e-4 * r.max(1.0) {
                return false;
            }
        }
    }
    true
}

/// Whether every sampled point of `s` is `r` from the axis (cylindrical test).
fn is_cylindrical(s: &dyn Surface, loc: &GpPnt, dir: &GpDir, r: f64) -> bool {
    let (u0, u1, v0, v1) = bounds(s);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let p = s.d0(u, v);
            if (dist_point_axis(&p, loc, dir) - r).abs() > 1e-4 * r.max(1.0) {
                return false;
            }
        }
    }
    true
}

/// Offset surface of a plane: translate by `distance` along its normal.
fn offset_plane_surface(s: &dyn Surface, distance: f64) -> Result<Box<dyn Surface>, String> {
    let (u0, _, v0, _) = bounds(s);
    let o = s.d0(u0, v0);
    let n = surface_normal(s, u0, v0);
    let nm = n.xyz().modulus();
    if nm < 1e-30 {
        return Err("offset_plane_surface: degenerate plane normal".into());
    }
    let n = GpVec::new(n.x() / nm, n.y() / nm, n.z() / nm);
    let nd = GpDir::from_vec(&n).map_err(|e| e.to_string())?;
    let z = GpDir::new(0.0, 0.0, 1.0).unwrap();
    let x_dir = if nd.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
    let pln = GpPln::new(GpAx3::new(o, nd, &x_dir).map_err(|e| e.to_string())?);
    Ok(Box::new(GeomPlane::new(pln.translated_vec(&n.multiplied_scalar(distance)))))
}

/// Offset surface of a general surface by sampling along its normal and
/// rebuilding a degree-1 bilinear B-spline surface (best-effort for
/// cone/torus/other surfaces).
fn offset_by_sampling(s: &dyn Surface, distance: f64) -> Result<Box<dyn Surface>, String> {
    let (u0, u1, v0, v1) = bounds(s);
    let (nu, nv) = (16usize, 16usize);
    let mut poles = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            let p = s.d0(u, v);
            let n = surface_normal(s, u, v);
            if n.xyz().square_modulus() < 1e-24 {
                return Err("offset_by_sampling: degenerate surface normal".into());
            }
            poles.push(p.translated_vec(&n.multiplied_scalar(distance)));
        }
    }
    let knots_u = occt_core::bspl::knots::build_uniform_knots(nu, 1);
    let knots_v = occt_core::bspl::knots::build_uniform_knots(nv, 1);
    let bs = crate::brep_faces::GeomBSplineSurface::new(poles, nu, nv, knots_u, knots_v, 1, 1)?;
    Ok(Box::new(bs))
}

/// The offset surface of `s` by `distance` along its own normal. Planes are
/// translated, spheres/cylinders change radius; cone/torus/other surfaces are
/// sampled and rebuilt as a bilinear B-spline surface.
pub fn offset_face_surface(s: &dyn Surface, distance: f64) -> Result<Box<dyn Surface>, String> {
    // Plane.
    if is_planar(s, 8, 8, 1e-6) {
        return offset_plane_surface(s, distance);
    }
    // Sphere.
    if let Some(center) = sphere_center(s) {
        let r = s.d0(0.5, 0.5).distance(&center);
        if r > 0.0 && is_spherical(s, &center, r) {
            let new_r = r + distance;
            if new_r <= 0.0 {
                return Err("offset_face_surface: sphere radius would become non-positive".into());
            }
            let pos = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
                .map_err(|e| e.to_string())?;
            let sph = GpSphere::new(pos, new_r).map_err(|e| e.to_string())?;
            return Ok(Box::new(GeomSphere::new(sph)));
        }
    }
    // Cylinder.
    if let Some((loc, dir, r)) = extract_cylinder(s) {
        if is_cylindrical(s, &loc, &dir, r) {
            let new_r = r + distance;
            if new_r <= 0.0 {
                return Err("offset_face_surface: cylinder radius would become non-positive".into());
            }
            let z = GpDir::new(0.0, 0.0, 1.0).unwrap();
            let x_dir = if dir.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
            let pos = GpAx3::new(loc, dir, &x_dir).map_err(|e| e.to_string())?;
            let cyl = GpCylinder::new(pos, new_r).map_err(|e| e.to_string())?;
            return Ok(Box::new(GeomCylinder::new(cyl)));
        }
    }
    offset_by_sampling(s, distance)
}

/// Offset a face along its surface normal by `distance`. The boundary wires of
/// a planar face are translated rigidly; other faces keep their wires as-is
/// (the edges are not re-projected onto the offset surface).
pub fn offset_face(face: &Face, distance: f64) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("offset_face: face has no surface")?;
    let new_surf = offset_face_surface(surf.as_ref(), distance)?;
    let b = TopoBuilder::new();
    let wires = wires_of_face(face);
    let mut new_wires = Vec::with_capacity(wires.len());
    if is_planar(new_surf.as_ref(), 8, 8, 1e-6) {
        // Rigid translation keeps the offset boundary on the offset plane.
        let n = surface_normal(new_surf.as_ref(), 0.0, 0.0);
        let v = n.multiplied_scalar(distance);
        for w in &wires {
            let w_edges = edges_of_wire(w);
            let mut new_edges = Vec::with_capacity(w_edges.len());
            for e in w_edges {
                let te = translated_copy(&e.0, &v).map_err(|e| e.to_string())?;
                new_edges.push(Edge(te));
            }
            new_wires.push(b.make_wire(&new_edges));
        }
    } else {
        for w in &wires {
            new_wires.push(w.clone());
        }
    }
    Ok(b.make_face(Arc::from(new_surf), &new_wires))
}

// ---------------------------------------------------------------------------
// 3. Shell / solid offset (BRepOffsetAPI_MakeOffsetShape)
// ---------------------------------------------------------------------------

/// Offset every vertex of a convex planar-face polyhedron by intersecting the
/// offset planes of the incident faces, then rebuild the boundary tree with
/// shared edges so the result is a closed shell.
fn offset_convex_polyhedron(shape: &TopoShape, faces: &[Face], distance: f64) -> Result<TopoShape, String> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return Err("offset_convex_polyhedron: no vertices".into());
    }
    let mut acc = GpXyz::zero();
    for v in &verts {
        acc = acc.added(&BRepTool::vertex_point(v).coord);
    }
    let centroid = GpPnt::from_xyz(&acc.divided(verts.len() as f64));

    struct FaceOff {
        pln: GpPln,
        offset_pln: GpPln,
        normal: GpVec,
    }
    let mut face_offs = Vec::with_capacity(faces.len());
    for f in faces {
        let pln = face_plane(f).ok_or("offset_convex_polyhedron: non-planar face")?;
        let n0 = *pln.axis().direction();
        let probe = GpVec::from_pnts(&centroid, &pln.location());
        let n = if n0.xyz().dot(&probe.xyz()) >= 0.0 { n0 } else { n0.reversed() };
        let nv = GpVec::from_xyz(n.xyz());
        for v in &verts {
            let p = BRepTool::vertex_point(v);
            let d = GpVec::from_pnts(&pln.location(), &p).xyz().dot(&n.xyz());
            if d > 1e-6 {
                return Err("offset_convex_polyhedron: shape is not convex".into());
            }
        }
        let offset_pln = pln.translated_vec(&nv.multiplied_scalar(distance));
        face_offs.push(FaceOff { pln, offset_pln, normal: nv });
    }

    let b = TopoBuilder::new();
    let mut new_verts: Vec<(GpPnt, Vertex, GpPnt)> = Vec::with_capacity(verts.len());
    for v in &verts {
        let p = BRepTool::vertex_point(v);
        let mut incident: Vec<&GpPln> = Vec::new();
        for fo in &face_offs {
            let d = GpVec::from_pnts(&fo.pln.location(), &p).xyz().dot(&fo.normal.xyz());
            if d.abs() < 1e-6 {
                incident.push(&fo.offset_pln);
            }
        }
        if incident.len() < 3 {
            return Err("offset_convex_polyhedron: vertex has fewer than 3 incident planes".into());
        }
        let np = plane_plane_plane_intersection(incident[0], incident[1], incident[2])
            .ok_or("offset_convex_polyhedron: parallel offset planes")?;
        new_verts.push((p, b.make_vertex(np, 0.0), np));
    }
    let find = |p: &GpPnt| new_verts.iter().find(|(orig, _, _)| orig.distance(p) < 1e-6);

    let edges = edges_of(shape);
    let mut emap: HashMap<usize, Edge> = HashMap::new();
    for e in &edges {
        let (a, z) = crate::topo_tools_full::edge_vertices(e);
        let (Some(va), Some(vb)) = (a, z) else {
            return Err("offset_convex_polyhedron: edge missing endpoint vertices".into());
        };
        let p1 = BRepTool::vertex_point(&va);
        let p2 = BRepTool::vertex_point(&vb);
        let (_, _, np1) = find(&p1).ok_or("offset_convex_polyhedron: endpoint not matched")?;
        let (_, _, np2) = find(&p2).ok_or("offset_convex_polyhedron: endpoint not matched")?;
        let seg = b.make_edge_segment(&np1, &np2);
        emap.insert(std::sync::Arc::as_ptr(&e.tshape) as usize, seg);
    }

    let mut new_faces = Vec::with_capacity(faces.len());
    for (f, fo) in faces.iter().zip(&face_offs) {
        let wires = wires_of_face(f);
        let mut new_wires = Vec::with_capacity(wires.len());
        for w in wires {
            let w_edges = edges_of_wire(&w);
            let mapped: Vec<Edge> = w_edges
                .iter()
                .map(|e| emap.get(&(std::sync::Arc::as_ptr(&e.tshape) as usize)).cloned())
                .collect::<Option<_>>()
                .ok_or("offset_convex_polyhedron: face edge not mapped")?;
            let nw = b.make_wire(&mapped);
            nw.set_closed(true);
            new_wires.push(nw);
        }
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(fo.offset_pln.clone()));
        new_faces.push(b.make_face(surface, &new_wires));
    }

    let shell = b.make_shell(&new_faces);
    if shape.shape_type() == ShapeType::Solid {
        Ok(b.make_solid(&[shell]).into())
    } else {
        Ok(shell.into())
    }
}

/// Offset a shell or solid. Convex planar-face shapes (boxes, prisms) are
/// rebuilt by intersecting the offset face planes; a single curved face (a
/// sphere) offsets in place. Other shapes return an error rather than a
/// broken result.
pub fn offset_shell(shape: &TopoShape, distance: f64) -> Result<TopoShape, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("offset_shell: no faces".into());
    }
    if faces.iter().all(face_is_planar) {
        if let Ok(s) = offset_convex_polyhedron(shape, &faces, distance) {
            return Ok(s);
        }
    }
    if faces.len() == 1 {
        let f = offset_face(&faces[0], distance)?;
        let b = TopoBuilder::new();
        let shell = b.make_shell(&[f]);
        if shape.shape_type() == ShapeType::Solid {
            return Ok(b.make_solid(&[shell]).into());
        }
        return Ok(shell.into());
    }
    Err("offset_shell: non-convex or curved offset unsupported".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    use crate::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};
    use crate::shape::Shell;
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6 * b.abs().max(1.0)
    }

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn z0_plane() -> GpPln {
        GpPln::new(GpAx3::standard())
    }

    #[test]
    fn offset_polygon_square_outward() {
        let plane = z0_plane();
        let pts = [
            GpPnt::new(-0.5, -0.5, 0.0),
            GpPnt::new(0.5, -0.5, 0.0),
            GpPnt::new(0.5, 0.5, 0.0),
            GpPnt::new(-0.5, 0.5, 0.0),
        ];
        let out = offset_polygon(&pts, &plane, 0.5).expect("offset");
        assert_eq!(out.len(), 4);
        for p in &out {
            assert!(approx(p.x().abs(), 1.0), "x {:?}", p);
            assert!(approx(p.y().abs(), 1.0), "y {:?}", p);
            assert!(approx(p.z(), 0.0));
        }
    }

    #[test]
    fn offset_polygon_square_inward() {
        let plane = z0_plane();
        let pts = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
        ];
        let out = offset_polygon(&pts, &plane, -0.5).expect("offset");
        assert_eq!(out.len(), 4);
        for p in &out {
            assert!(approx(p.x().abs(), 0.5), "x {:?}", p);
            assert!(approx(p.y().abs(), 0.5), "y {:?}", p);
        }
    }

    #[test]
    fn offset_circle_radius() {
        let plane = z0_plane();
        let c = offset_circle(GpPnt::zero(), 1.0, &plane, 0.5).expect("offset +");
        assert!(approx(c.radius(), 1.5));
        let c2 = offset_circle(GpPnt::zero(), 1.0, &plane, -0.5).expect("offset -");
        assert!(approx(c2.radius(), 0.5));
    }

    #[test]
    fn offset_wire_2d_closed_square() {
        let b = TopoBuilder::new();
        let pts = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
        ];
        let edges = [
            b.make_edge_segment(&pts[0], &pts[1]),
            b.make_edge_segment(&pts[1], &pts[2]),
            b.make_edge_segment(&pts[2], &pts[3]),
            b.make_edge_segment(&pts[3], &pts[0]),
        ];
        let wire = b.make_wire(&edges);
        let out = offset_wire_2d(&wire, &z0_plane(), 0.3).expect("offset");
        for v in vertices_of(&out.0) {
            let p = BRepTool::vertex_point(&v);
            assert!(approx(p.x().abs(), 1.3), "x {}", p.x());
            assert!(approx(p.y().abs(), 1.3), "y {}", p.y());
        }
        clear_tree(&out.0);
    }

    #[test]
    fn offset_face_plane_translates() {
        let pts = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
        ];
        let face = crate::brep_builder_api::make_face_from_polygon(&pts).expect("square face");
        let off = offset_face(&face, 1.0).expect("offset");
        let s = BRepTool::face_surface(&off).expect("surface");
        let p = s.d0(0.0, 0.0);
        assert!(approx(p.z(), 1.0), "z {}", p.z());
        clear_tree(&off.0);
    }

    #[test]
    fn offset_face_sphere_grows() {
        let s = BRepPrimSphere::make_sphere(1.0);
        let face = faces_of(&s.solid.0)[0].clone();
        let off = offset_face(&face, 0.5).expect("offset");
        let surf = BRepTool::face_surface(&off).expect("surface");
        let c = sphere_center(surf.as_ref()).expect("center");
        assert!(approx(c.distance(&GpPnt::zero()), 0.0), "center {:?}", c);
        let r = surf.d0(0.0, 0.0).distance(&c);
        assert!(approx(r, 1.5), "radius {}", r);
        clear_tree(&off.0);
    }

    #[test]
    fn offset_face_cylinder_radius() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let lateral = faces_of(&c.solid.0)
            .into_iter()
            .find(|f| !face_is_planar(f))
            .expect("lateral face");
        let off = offset_face(&lateral, 0.2).expect("offset");
        let surf = BRepTool::face_surface(&off).expect("surface");
        let (_, _, r) = extract_cylinder(surf.as_ref()).expect("cylinder params");
        assert!(approx(r, 1.2), "radius {}", r);
        clear_tree(&off.0);
    }

    #[test]
    fn offset_solid_box_grows() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let result = offset_shell(&b.solid.0, 0.5).expect("offset");
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut zs = Vec::new();
        for v in vertices_of(&result) {
            let p = BRepTool::vertex_point(&v);
            xs.push(p.x());
            ys.push(p.y());
            zs.push(p.z());
        }
        let (min, max) = (
            |v: &Vec<f64>| v.iter().cloned().fold(f64::INFINITY, f64::min),
            |v: &Vec<f64>| v.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        assert!(approx(min(&xs), -0.5) && approx(max(&xs), 2.5), "x [{},{}]", min(&xs), max(&xs));
        assert!(approx(min(&ys), -0.5) && approx(max(&ys), 2.5));
        assert!(approx(min(&zs), -0.5) && approx(max(&zs), 2.5));
        // Side length 3.0.
        assert!(approx(max(&xs) - min(&xs), 3.0));
        // Closed manifold shell.
        let shell = Shell(result.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "offset box must be a closed shell");
        clear_tree(&result);
    }

    #[test]
    fn offset_solid_box_shrinks() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let result = offset_shell(&b.solid.0, -0.4).expect("offset");
        let xs: Vec<f64> = vertices_of(&result)
            .iter()
            .map(|v| BRepTool::vertex_point(v).x())
            .collect();
        let ys: Vec<f64> = vertices_of(&result)
            .iter()
            .map(|v| BRepTool::vertex_point(v).y())
            .collect();
        let mn = |v: &Vec<f64>| v.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = |v: &Vec<f64>| v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(approx(mn(&xs), 0.4) && approx(mx(&xs), 1.6), "x [{},{}]", mn(&xs), mx(&xs));
        assert!(approx(mn(&ys), 0.4) && approx(mx(&ys), 1.6));
        assert!(approx(mx(&xs) - mn(&xs), 1.2), "side length");
        clear_tree(&result);
    }

    #[test]
    fn offset_negative_shrinks_sphere() {
        let s = BRepPrimSphere::make_sphere(1.0);
        let result = offset_shell(&s.solid.0, -0.3).expect("offset");
        let face = faces_of(&result)[0].clone();
        let surf = BRepTool::face_surface(&face).expect("surface");
        let c = sphere_center(surf.as_ref()).expect("center");
        let r = surf.d0(0.0, 0.0).distance(&c);
        assert!(approx(r, 0.7), "radius {}", r);
        clear_tree(&result);
    }

    #[test]
    fn plane_plane_plane_corner() {
        let px = GpPln::new(
            GpAx3::new(GpPnt::new(1.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap(), &GpDir::new(0.0, 1.0, 0.0).unwrap())
                .expect("x plane"),
        );
        let py = GpPln::new(
            GpAx3::new(GpPnt::new(0.0, 1.0, 0.0), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(0.0, 0.0, 1.0).unwrap())
                .expect("y plane"),
        );
        let pz = GpPln::new(
            GpAx3::new(GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
                .expect("z plane"),
        );
        let pt = plane_plane_plane_intersection(&px, &py, &pz).expect("corner");
        assert!(pt.distance(&GpPnt::new(1.0, 1.0, 1.0)) < 1e-9, "{:?}", pt);
        // Two parallel planes (x=1 and x=2) with z=1: no common point.
        let px2 = GpPln::new(
            GpAx3::new(GpPnt::new(2.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap(), &GpDir::new(0.0, 1.0, 0.0).unwrap())
                .expect("x=2 plane"),
        );
        assert!(plane_plane_plane_intersection(&px, &px2, &pz).is_none());
    }

    #[test]
    fn offset_wire_3d_circle() {
        let b = TopoBuilder::new();
        let ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("circle frame");
        let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(ax3.ax2(), 1.0)));
        let mut e = b.make_edge(curve.clone(), 0.0, 2.0 * PI);
        let seam = b.make_vertex(curve.d0(0.0), 0.0);
        b.add(&mut e.0, &seam.0);
        b.add(&mut e.0, &seam.0);
        let wire = b.make_wire(&[e]);
        let out = offset_wire_3d(&wire, 0.2).expect("offset");
        for v in vertices_of(&out.0) {
            let p = BRepTool::vertex_point(&v);
            assert!(approx(p.distance(&GpPnt::zero()), 1.2), "radius {}", p.distance(&GpPnt::zero()));
        }
        clear_tree(&out.0);
    }
}
