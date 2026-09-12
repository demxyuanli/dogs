use super::prelude::*;
use super::*;

/// Parameter intervals of the line `p0 + t·dir` that lie inside a 2D polygon
/// (in the plane's frame). `t` matches the 3D line's arc-length parameter.
pub(super) fn line_polygon_t_intervals(pln: &GpPln, p0: &GpPnt, dir: &GpVec, poly: &[GpPnt2d]) -> Vec<(f64, f64)> {
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

pub(super) fn merge_intervals(mut v: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
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

pub(super) fn intersect_intervals(a: &[(f64, f64)], b: &[(f64, f64)], tol: f64) -> Vec<(f64, f64)> {
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
pub(super) fn face_polygon_2d_with_plane(f: &Face, pln: &GpPln) -> Option<Vec<GpPnt2d>> {
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
pub(super) fn project_to_2d(pln: &GpPln, p: &GpPnt) -> GpPnt2d {
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
pub(super) fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    pub(super) const GOLD: f64 = 0.6180339887498949;
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
