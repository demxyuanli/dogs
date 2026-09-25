use super::prelude::*;
use super::*;

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).

pub(super) fn bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// 3×3 determinant (Cramer's-rule helper).
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system via Cramer's rule.
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

/// Circumcenter of three non-collinear 3D points, if it exists.
pub(super) fn circumcenter3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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
pub(super) fn line_line_intersect_3d(a1: &GpPnt, d1: &GpVec, a2: &GpPnt, d2: &GpVec, n: &GpVec) -> Option<GpPnt> {
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
pub(super) fn to_2d(plane: &GpPln, p: &GpPnt) -> (f64, f64) {
    let ax = plane.position();
    let d = GpVec::from_pnts(&ax.location(), p);
    to_2d_vec(plane, &d)
}

/// Project a 3D vector onto the plane's in-plane axes.
pub(super) fn to_2d_vec(plane: &GpPln, v: &GpVec) -> (f64, f64) {
    let ax = plane.position();
    (v.xyz().dot(ax.x_direction().xyz()), v.xyz().dot(ax.y_direction().xyz()))
}

/// Reconstruct a 3D point from plane-frame coordinates.
pub(super) fn to_3d(plane: &GpPln, u: f64, v: f64) -> GpPnt {
    let ax = plane.position();
    let o = ax.location();
    let x = GpVec::from_xyz(ax.x_direction().xyz());
    let y = GpVec::from_xyz(ax.y_direction().xyz());
    o.translated_vec(&x.multiplied_scalar(u)).translated_vec(&y.multiplied_scalar(v))
}

/// Intersection of two infinite 2D lines.
pub(super) fn line_line_2d(a1: (f64, f64), d1: (f64, f64), a2: (f64, f64), d2: (f64, f64)) -> Option<(f64, f64)> {
    let den = d1.0 * d2.1 - d1.1 * d2.0;
    if den.abs() < 1e-12 {
        return None;
    }
    let t = ((a2.0 - a1.0) * d2.1 - (a2.1 - a1.1) * d2.0) / den;
    Some((a1.0 + t * d1.0, a1.1 + t * d1.1))
}

/// Intersection of a 2D line and a 2D circle (0, 1 or 2 points).
pub(super) fn line_circle_2d(a: (f64, f64), d: (f64, f64), c: (f64, f64), r: f64) -> Vec<(f64, f64)> {
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
pub(super) fn circle_circle_2d(c1: (f64, f64), r1: f64, c2: (f64, f64), r2: f64) -> Vec<(f64, f64)> {
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
pub(super) enum OffsetCurveKind {
    Line,
    Circle,
}

/// Per-edge classification for wire offsetting.
pub(super) struct EdgeInfo {
    pub(super) kind: OffsetCurveKind,
    pub(super) curve: Arc<dyn Curve>,
    pub(super) first: f64,
    pub(super) last: f64,
    pub(super) center: GpPnt,
    pub(super) radius: f64,
}

/// Line parameters from samples (a point and a unit direction), if collinear.
pub(super) fn line_params(pts: &[GpPnt]) -> Option<(GpPnt, GpVec)> {
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
pub(super) fn circle_params(pts: &[GpPnt]) -> Option<(GpPnt, f64)> {
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
pub(super) fn classify_edge(e: &Edge) -> Result<EdgeInfo, String> {
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
pub(super) enum OffsetCurve2 {
    Line { a: (f64, f64), d: (f64, f64) },
    Circle { c: (f64, f64), r: f64 },
}

/// Build the offset of a single edge. Returns the offset curve plus the offset
/// arc midpoint (for circle edges) used when rebuilding the arc.
pub(super) fn offset_curve_for_edge(
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
pub(super) fn intersect_offset_curves(a: &OffsetCurve2, b: &OffsetCurve2, orig: &GpPnt, plane: &GpPln) -> Option<GpPnt> {
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
pub(super) fn rebuild_circle_wire(circ: &GpCirc, first: f64, last: f64) -> Result<Wire, String> {
    let b = TopoBuilder::new();
    let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ.clone()));
    let mut e = b.make_edge(curve.clone(), first, last);
    let seam = curve.d0(first);
    let v = b.make_vertex(seam, 0.0);
    b.add_edge_vertices(&mut e, &v, &v);
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
pub(super) fn dist_point_axis(p: &GpPnt, loc: &GpPnt, dir: &GpDir) -> f64 {
    let v = GpVec::from_pnts(loc, p);
    let d = dir.xyz();
    let proj = d.multiplied(v.xyz().dot(d));
    GpVec::from_xyz(&v.xyz().subtracted(&proj)).xyz().modulus()
}

/// Extract cylinder parameters (a point on the axis, the axis direction, the
/// radius) by fitting two circles to rings at different heights.
pub(super) fn extract_cylinder(s: &dyn Surface) -> Option<(GpPnt, GpDir, f64)> {
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
pub(super) fn is_spherical(s: &dyn Surface, center: &GpPnt, r: f64) -> bool {
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
pub(super) fn is_cylindrical(s: &dyn Surface, loc: &GpPnt, dir: &GpDir, r: f64) -> bool {
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
pub(super) fn offset_plane_surface(s: &dyn Surface, distance: f64) -> Result<Box<dyn Surface>, String> {
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
pub(super) fn offset_by_sampling(s: &dyn Surface, distance: f64) -> Result<Box<dyn Surface>, String> {
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
    } else if let (Some(c_src), Some(c_new)) = (
        sphere_center(surf.as_ref()),
        sphere_center(new_surf.as_ref()),
    ) {
        // `BRepOffset`'s analytic branch for `GeomAbs_Sphere`: a spherical face
        // offsets to a *concentric* sphere and every boundary curve maps with it
        // (`BRepOffset_MakeOffset` → `BRepOffset::MakeOffset`; the map is
        // `gp_Trsf::SetScale(centre, (r + d) / r)`). The source wire must not be
        // reused: it lies on the *source* surface and carries no
        // (edge, offset-face) pcurve, so the analytic passes (`BRepGProp`) see a
        // pcurve-less face and report a zero volume.
        let r_src = surf.d0(0.0, 0.0).distance(&c_src);
        let r_new = new_surf.d0(0.0, 0.0).distance(&c_new);
        let mut rebuilt = false;
        if r_src > 0.0 && r_new > 0.0 {
            let mut t = occt_core::gp::GpTrsf::identity();
            t.set_scale(&c_src, r_new / r_src).map_err(|e| e.to_string())?;
            let mut copied: std::collections::HashMap<usize, Edge> =
                std::collections::HashMap::new();
            for w in &wires {
                let w_edges = edges_of_wire(w);
                let mut new_edges = Vec::with_capacity(w_edges.len());
                for e in w_edges {
                    // A seam edge appears twice in the wire; both occurrences must
                    // stay the *same* TShape (only their orientation differs) or
                    // `ShapeFix_Wire` no longer sees a seam and hands both sides
                    // the same pcurve, whose boundary terms then cancel
                    // (`data/Offset.step`'s periodic faces behave the same way).
                    let k = std::sync::Arc::as_ptr(&e.0.tshape) as usize;
                    let mut occ = match copied.get(&k) {
                        Some(c) => c.clone(),
                        None => {
                            let te = crate::shape_ops::transformed_copy(&e.0, &t)
                                .map_err(|e| e.to_string())?;
                            let c = Edge(te);
                            copied.insert(k, c.clone());
                            c
                        }
                    };
                    occ.0.set_orientation(e.orientation());
                    new_edges.push(occ);
                }
                new_wires.push(b.make_wire(&new_edges));
            }
            rebuilt = true;
        }
        if !rebuilt {
            for w in &wires {
                new_wires.push(w.clone());
            }
        }
    } else {
        for w in &wires {
            new_wires.push(w.clone());
        }
    }
    let face_out = b.make_face(Arc::from(new_surf), &new_wires);
    // `BRepOffset_MakeOffset` runs `BRepLib::SameParameter` / `ShapeFix` on its
    // result, so every edge of a rebuilt face carries a pcurve on the *offset*
    // surface (`ShapeFix_Edge::FixAddPCurve`, `ShapeFix_Edge.cxx:517-534`).
    // Only the rebuilt wires are touched: the reused-clone branch shares its
    // TShapes with the source shape, and repairing those in place would mutate
    // the input.
    {
        let mut ws = wires_of_face(&face_out);
        for w in ws.iter_mut() {
            crate::shhealing::check_pcurves_and_shift(w, &face_out, occt_core::precision::CONFUSION);
        }
    }
    Ok(face_out)
}
