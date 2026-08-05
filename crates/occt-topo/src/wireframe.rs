//! Wireframe & face tessellation.
//!
//! Port of the *meshing* half of `BRepMesh_IncrementalMesh`: an edge is
//! turned into a deflection-bounded polyline (`GCPnts_UniformDeflection`) and a
//! face into a UV-grid triangle soup. Geometry is read from the side-table
//! registry (`tgeometry::GeometryRegistry`) — the same data `BRep_Tool` reads.

use occt_core::gp::{GpPnt, GpPnt2d, GpXyz};
use occt_core::gcpnts::{CurveSample, UniformDeflection, UniformPoints};
use occt_core::poly::triangulation::Triangle;
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::brep_surface::{face_is_planar, face_plane};
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, wires_of_face};

fn reg() -> &'static GeometryRegistry {
    GeometryRegistry::global()
}

/// Adapter making a `&dyn Curve` sampleable by the `gcpnts` samplers.
struct CurveAdapter<'a> {
    curve: &'a dyn Curve,
}
impl CurveSample for CurveAdapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.curve.d0(u)
    }
}

/// Discretize an edge into a polyline of 3D points.
///
/// The underlying curve is sampled with [`UniformDeflection`] over the edge's
/// parameter range; both endpoints are always included. If the edge has no
/// registered curve (or an unbounded range) an empty polyline is returned;
/// a fallback of 64 uniform samples covers the degenerate-adaptive case.
pub fn edge_to_polyline(e: &Edge, deflection: f64) -> Vec<GpPnt> {
    let Some(curve) = reg().edge_curve(e) else { return Vec::new() };
    let (f0, f1) = reg().edge_parameters(e);
    let (a, b) = if f0.is_finite() && f1.is_finite() && f1 >= f0 { (f0, f1) } else { return Vec::new() };
    if b - a < 1e-12 {
        return vec![curve.d0(a)];
    }
    let tol = deflection.max(1e-9);
    let sampled = UniformDeflection::from_curve_with_deflection(&CurveAdapter { curve: curve.as_ref() }, a, b, tol);
    if sampled.points.len() >= 2 {
        sampled.points
    } else {
        UniformPoints::from_curve(&CurveAdapter { curve: curve.as_ref() }, a, b, 64).points
    }
}

/// Maximum deviation of a sampled polyline from the true edge curve.
///
/// Each polyline segment is mapped to a parameter span by chord-length
/// parametrization and the curve is sampled at 64 interior points per span;
/// the furthest distance to the chord is returned.
pub fn edge_chord_error(e: &Edge, polyline: &[GpPnt]) -> f64 {
    let n = polyline.len();
    if n < 2 { return 0.0; }
    let Some(curve) = reg().edge_curve(e) else { return 0.0 };
    let (f0, f1) = reg().edge_parameters(e);
    if !(f0.is_finite() && f1.is_finite()) || f1 <= f0 { return 0.0; }

    // Cumulative chord lengths → curve parameter per polyline point.
    let mut s = Vec::with_capacity(n);
    s.push(0.0);
    for i in 1..n {
        s.push(s[i - 1] + polyline[i].distance(&polyline[i - 1]));
    }
    let total = *s.last().unwrap();
    if total <= 1e-12 { return 0.0; }

    let mut max_dev = 0.0f64;
    for i in 0..n - 1 {
        let ua = f0 + (f1 - f0) * s[i] / total;
        let ub = f0 + (f1 - f0) * s[i + 1] / total;
        let pa = &polyline[i];
        let pb = &polyline[i + 1];
        for k in 1..64 {
            let u = ua + (ub - ua) * k as f64 / 64.0;
            let d = point_segment_dist(&curve.d0(u), pa, pb);
            if d > max_dev { max_dev = d; }
        }
    }
    max_dev
}

/// Wind a triangle so its normal agrees with `n`.
pub(crate) fn orient3(pts: &[GpPnt], i0: usize, i1: usize, i2: usize, n: &GpXyz) -> Triangle {
    let pa = pts[i0].coord;
    let pb = pts[i1].coord;
    let pc = pts[i2].coord;
    let nn = pb.subtracted(&pa).crossed(&pc.subtracted(&pa));
    if n.dot(&nn) < 0.0 {
        Triangle::new(i0, i2, i1)
    } else {
        Triangle::new(i0, i1, i2)
    }
}

/// Whether a (planar, non-self-intersecting) polygon is convex when seen along
/// `n` (all consecutive-triple signed areas share a sign).
pub(crate) fn is_convex(pts: &[GpPnt], n: &GpXyz) -> bool {
    let mut sign: Option<f64> = None;
    let m = pts.len();
    for i in 0..m {
        let a = pts[i].coord;
        let b = pts[(i + 1) % m].coord;
        let c = pts[(i + 2) % m].coord;
        let d = b.subtracted(&a).crossed(&c.subtracted(&a)).dot(n);
        if d.abs() > 1e-12 {
            match sign {
                None => sign = Some(d.signum()),
                Some(s) if s != d.signum() => return false,
                _ => {}
            }
        }
    }
    true
}

/// Whether `p` is STRICTLY inside the 2D triangle `(a, b, c)` (boundary and
/// collinear points excluded — ear clipping must not let an edge point block an
/// ear).
fn point_in_triangle_2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let s = |u: &GpPnt2d, v: &GpPnt2d, w: &GpPnt2d| {
        (v.x() - u.x()) * (w.y() - u.y()) - (v.y() - u.y()) * (w.x() - u.x())
    };
    let d1 = s(p, a, b);
    let d2 = s(p, b, c);
    let d3 = s(p, c, a);
    if d1 == 0.0 || d2 == 0.0 || d3 == 0.0 {
        return false;
    }
    (d1 > 0.0 && d2 > 0.0 && d3 > 0.0) || (d1 < 0.0 && d2 < 0.0 && d3 < 0.0)
}

/// Ear-clipping triangulation of a simple (possibly non-convex) 2D polygon.
/// Returns `(i0, i1, i2)` vertex-index triples. `None` when degenerate.
fn ear_clip(poly: &[GpPnt2d]) -> Option<Vec<(usize, usize, usize)>> {
    if poly.len() < 3 {
        return None;
    }
    // Overall orientation sign (cross-product of the first convex vertex).
    let area = {
        let mut a = 0.0;
        for i in 0..poly.len() {
            let p = poly[i];
            let q = poly[(i + 1) % poly.len()];
            a += p.x() * q.y() - q.x() * p.y();
        }
        a
    };
    if area.abs() < 1e-12 {
        return None;
    }
    let sign = area.signum();
    let mut idx: Vec<usize> = (0..poly.len()).collect();
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();
    let mut guard = 0;
    while idx.len() > 3 && guard < poly.len() * poly.len() {
        guard += 1;
        let m = idx.len();
        let mut clipped = false;
        for k in 0..m {
            let (a, b, c) = (idx[k], idx[(k + 1) % m], idx[(k + 2) % m]);
            let cross = (poly[b].x() - poly[a].x()) * (poly[c].y() - poly[b].y())
                - (poly[b].y() - poly[a].y()) * (poly[c].x() - poly[b].x());
            if cross * sign <= 0.0 {
                continue; // reflex or degenerate vertex
            }
            let mut has_inside = false;
            for &i in &idx {
                if i != a && i != b && i != c
                    && point_in_triangle_2d(&poly[a], &poly[b], &poly[c], &poly[i])
                {
                    has_inside = true;
                    break;
                }
            }
            if has_inside {
                continue;
            }
            tris.push((a, b, c));
            idx.remove((k + 1) % m);
            clipped = true;
            break;
        }
        if !clipped {
            return None;
        }
    }
    if idx.len() == 3 {
        tris.push((idx[0], idx[1], idx[2]));
        Some(tris)
    } else {
        None
    }
}

/// Bridge every hole to the outer boundary by its closest vertex pair, yielding
/// a single simple polygon (the face's material) suitable for ear clipping.
fn bridge_holes(outer: &[GpPnt2d], holes: &[&[GpPnt2d]]) -> Vec<GpPnt2d> {
    let signed_area = |pts: &[GpPnt2d]| -> f64 {
        let mut a = 0.0;
        for i in 0..pts.len() {
            let p = pts[i];
            let q = pts[(i + 1) % pts.len()];
            a += p.x() * q.y() - q.x() * p.y();
        }
        a
    };
    let outer_sign = signed_area(outer).signum();
    let mut poly = outer.to_vec();
    for hole in holes {
        // Traverse the hole OPPOSITE to the outer so the bridged polygon is a
        // simple ring cut open along the bridge (not a self-crossing loop).
        let reverse = signed_area(hole).signum() == outer_sign;
        let (mut oi, mut hi) = (0usize, 0usize);
        let mut best = f64::INFINITY;
        for (i, op) in poly.iter().enumerate() {
            for (j, hp) in hole.iter().enumerate() {
                let d = op.distance(hp);
                if d < best {
                    best = d;
                    oi = i;
                    hi = j;
                }
            }
        }
        let hlen = hole.len();
        let mut new_poly = Vec::with_capacity(poly.len() + hlen + 2);
        new_poly.extend_from_slice(&poly[..=oi]);
        for k in 0..=hlen {
            let idx = if reverse { (hi + hlen - k) % hlen } else { (hi + k) % hlen };
            new_poly.push(hole[idx]);
        }
        new_poly.extend_from_slice(&poly[oi..]);
        poly = new_poly;
    }
    poly
}

/// Triangulate a convex planar face from its boundary polygon.
///
/// Boundary points are collected from each wire edge with the given deflection,
/// deduplicated, and ordered radially around the face centroid (recovering the
/// boundary order of a convex polygon), then fan-triangulated with outward
/// winding along the plane axis. A face with holes (outer + hole wires, e.g. a
/// box top ring cut by a boss cylinder's base) is triangulated by ear-clipping
/// the hole-bridged polygon. Returns `None` when the boundary is degenerate or
/// non-convex — the caller then falls back to a UV grid.
pub(crate) fn planar_polygon_triangulate(
    face: &Face,
    bound_def: f64,
) -> Option<(Vec<GpPnt>, Vec<Triangle>)> {
    // Per-wire boundary point loops.
    let mut loops: Vec<Vec<GpPnt>> = Vec::new();
    for wire in wires_of_face(face) {
        let mut loop_pts: Vec<GpPnt> = Vec::new();
        for e in edges_of_wire(&wire) {
            for p in edge_to_polyline(&e, bound_def) {
                if !loop_pts.iter().any(|q| q.distance(&p) < 1e-9) {
                    loop_pts.push(p);
                }
            }
        }
        if loop_pts.len() >= 3 {
            loops.push(loop_pts);
        }
    }
    if loops.is_empty() {
        return None;
    }

    let pln = face_plane(face).unwrap_or_else(occt_core::gp::GpPln::default);
    let n = *pln.axis().direction().xyz();
    let xdir = *pln.x_axis().direction().xyz();
    let ydir = xdir.crossed(&n);

    // Order each loop radially around the face centroid (convex recovery).
    let mut all_pts: Vec<GpPnt> = Vec::new();
    for l in &loops {
        all_pts.extend(l.iter().cloned());
    }
    let centroid = all_pts
        .iter()
        .fold(GpPnt::new(0.0, 0.0, 0.0), |acc, p| {
            GpPnt::new(acc.x() + p.x(), acc.y() + p.y(), acc.z() + p.z())
        });
    let inv = 1.0 / all_pts.len() as f64;
    let c = GpPnt::new(centroid.x() * inv, centroid.y() * inv, centroid.z() * inv);
    let order = |l: &mut Vec<GpPnt>| {
        l.sort_by(|p, q| {
            let ap = p.coord.subtracted(&c.coord);
            let aq = q.coord.subtracted(&c.coord);
            let ang_p = ap.dot(&ydir).atan2(ap.dot(&xdir));
            let ang_q = aq.dot(&ydir).atan2(aq.dot(&xdir));
            ang_p.partial_cmp(&ang_q).unwrap_or(std::cmp::Ordering::Equal)
        });
    };
    let mut pts = loops[0].clone();
    order(&mut pts);
    let to2d = |p: &GpPnt| -> GpPnt2d {
        let v = p.coord.subtracted(&pln.location().coord);
        GpPnt2d::new(v.dot(&xdir), v.dot(&ydir))
    };

    if loops.len() == 1 {
        // Single boundary: convex fan, else ear-clip.
        if !is_convex(&pts, &n) {
            let poly2d: Vec<GpPnt2d> = pts.iter().map(to2d).collect();
            let tris2 = ear_clip(&poly2d)?;
            let tris: Vec<Triangle> = tris2
                .iter()
                .map(|&(a, b, c2)| orient3(&pts, a, b, c2, &n))
                .collect();
            let pts3 = pts;
            return Some((pts3, tris));
        }
    } else {
        // Outer + holes: identify the outer (largest area), bridge holes, ear-clip.
        let mut area_of: Vec<f64> = Vec::new();
        for l in &loops {
            let poly2d: Vec<GpPnt2d> = l.iter().map(to2d).collect();
            let mut a = 0.0;
            for i in 0..poly2d.len() {
                let p = poly2d[i];
                let q = poly2d[(i + 1) % poly2d.len()];
                a += p.x() * q.y() - q.x() * p.y();
            }
            area_of.push(a.abs());
        }
        let outer_i = area_of
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)?;
        let outer_pts = loops[outer_i].clone();
        let outer2d: Vec<GpPnt2d> = outer_pts.iter().map(to2d).collect();
        let holes2d: Vec<Vec<GpPnt2d>> = loops
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != outer_i)
            .map(|(_, l)| l.iter().map(to2d).collect())
            .collect();
        let holes_refs: Vec<&[GpPnt2d]> = holes2d.iter().map(|v| v.as_slice()).collect();
        let bridged = bridge_holes(&outer2d, &holes_refs);
        let tris2 = ear_clip(&bridged)?;
        // Map triangulated 2D points back to 3D.
        let pts3: Vec<GpPnt> = bridged.iter().map(|p| {
            let o = pln.location().coord;
            let v = xdir.multiplied(p.x()).added(&ydir.multiplied(p.y()));
            GpPnt::new(o.x() + v.x(), o.y() + v.y(), o.z() + v.z())
        }).collect();
        let tris: Vec<Triangle> = tris2
            .iter()
            .map(|&(a, b, c2)| orient3(&pts3, a, b, c2, &n))
            .collect();
        return Some((pts3, tris));
    }

    let mut tris = Vec::with_capacity(pts.len() - 2);
    for i in 1..pts.len() - 1 {
        tris.push(orient3(&pts, 0, i, i + 1, &n));
    }
    // A REVERSED face bounds the volume with inward-pointing parametric normals
    // (face_plane's axis is not orientation-aware); flip so the mesh points
    // outward, matching `face_to_triangles`' UV-grid convention.
    if face.0.orientation() == crate::abs::Orientation::Reversed {
        for t in tris.iter_mut() {
            std::mem::swap(&mut t.n1, &mut t.n2);
        }
    }
    Some((pts, tris))
}

/// Tessellate a face into a UV-grid triangle soup (vertices + triangles).
///
/// Planar faces triangulate exactly from their boundary polygon (see
/// [`planar_polygon_triangulate`]) — a UV grid of the bounding box would
/// over-cover a disk cap as its circumscribing square. Curved faces fall back
/// to the grid below. Grid resolution is driven by `deflection` relative to
/// the UV domain size (`nu = clamp(ceil(du/deflection)+1, 3, 64)`). Triangle
/// winding follows the surface normal from `d1` so the resulting mesh points
/// outward. Degenerate cells are skipped. Missing surface or unbounded domain
/// yields an empty soup.
pub fn face_to_triangles(f: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
    if face_is_planar(f) {
        if let Some(m) = planar_polygon_triangulate(f, deflection.max(0.01)) {
            return m;
        }
    }
    let Some(surface) = reg().face_surface(f) else { return (Vec::new(), Vec::new()) };
    let (u0, u1, v0, v1) = face_uv_bounds(f, surface.as_ref());
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        return (Vec::new(), Vec::new());
    }
    let (du, dv) = (u1 - u0, v1 - v0);
    if du <= 0.0 || dv <= 0.0 { return (Vec::new(), Vec::new()); }

    let def = deflection.max(1e-9);
    let nu = ((du / def).ceil() as usize + 1).clamp(3, 64);
    let nv = ((dv / def).ceil() as usize + 1).clamp(3, 64);

    let mut pts = Vec::with_capacity(nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let u = u0 + du * i as f64 / (nu - 1) as f64;
            let v = v0 + dv * j as f64 / (nv - 1) as f64;
            pts.push(surface.d0(u, v));
        }
    }

    let idx = |i: usize, j: usize| j * nu + i;
    let mut tris = Vec::new();
    for j in 0..nv - 1 {
        for i in 0..nu - 1 {
            let a = idx(i, j);
            let b = idx(i + 1, j);
            let c = idx(i, j + 1);
            let d = idx(i + 1, j + 1);
            let uc = u0 + du * (i as f64 + 0.5) / (nu - 1) as f64;
            let vc = v0 + dv * (j as f64 + 0.5) / (nv - 1) as f64;
            let (_, su, sv) = surface.d1(uc, vc);
            let ns = su.coord.crossed(&sv.coord);
            let t1 = orient(&pts, a, b, c, &ns);
            let t2 = orient(&pts, b, d, c, &ns);
            if triangle_area(&pts, &t1) > 1e-14 && triangle_area(&pts, &t2) > 1e-14 {
                tris.push(t1);
                tris.push(t2);
            }
        }
    }
    // A REVERSED face bounds the volume with inward-pointing surface normals;
    // the triangle winding must be flipped so the signed volume is computed
    // with the correct orientation.
    if f.0.orientation() == crate::abs::Orientation::Reversed {
        for t in tris.iter_mut() {
            std::mem::swap(&mut t.n1, &mut t.n2);
        }
    }
    (pts, tris)
}

/// Sum of triangle areas of a face's tessellation.
pub fn face_area(f: &Face) -> f64 {
    let (pts, tris) = face_to_triangles(f, 0.01);
    tris.iter().map(|t| triangle_area(&pts, t)).sum()
}

/// UV domain of a face.
///
/// Prefers the surface's own finite natural range; otherwise derives a box by
/// inverting the face's boundary edge curves (wire children) onto the surface
/// with a Newton iterate on `d1`. Falls back to `[0,1]²` when nothing can be
/// determined.
pub(crate) fn face_uv_bounds(f: &Face, surface: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() && u1 > u0 && v1 > v0 {
        return (u0, u1, v0, v1);
    }

    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    // The face's UV domain is bounded by the UV image of its boundary edges
    // (`BRepAdaptor_Surface::UVBounds` / `BRep_Tool` semantics). Each edge has
    // a pcurve on this face; its endpoints' `(u, v)` delimit the edge's UV
    // span. Inverting 3D boundary points numerically (`invert_uv`) is fragile
    // on revolution surfaces (a cylinder's v-range is infinite on the surface),
    // so we take the pcurve's UV directly.
    //
    // Periodic u (seam): a full circle edge spans one period; its pcurve u
    // runs 0→−2π. The u domain must therefore cover a full period, not
    // collapse to ~0, so we unwrap each edge's u span to the surface's
    // [u0, u0+period] window.
    let (su0, su1) = surface.u_range();
    let u_period = if su0.is_finite() && su1.is_finite() && su1 > su0 { su1 - su0 } else { 0.0 };
    let mut full_period = false;
    for w in f.tshape.read().unwrap().children.iter() {
        let w = TopoShape::from_handle(w.clone());
        if w.shape_type() != ShapeType::Wire { continue; }
        for eh in w.tshape.read().unwrap().children.iter() {
            let e = TopoShape::from_handle(eh.clone());
            if e.shape_type() != ShapeType::Edge { continue; }
            let edge = Edge(e.clone());
            let Ok(pc) = crate::pcurve_full::make_pcurve_full(&edge, f) else { continue };
            let (a0, a1) = reg().edge_parameters(&e);
            if !(a0.is_finite() && a1.is_finite() && a1 > a0) { continue; }
            // Sample the pcurve along the edge. The pcurve of a full-circle
            // edge is unwrapped monotonically (u runs 0 → −2π), so its UV
            // endpoints differ by the period rather than coinciding; the u-span
            // must be measured in continuously-unwrapped coordinates.
            let mut cu_min = f64::INFINITY;
            let mut cu_max = f64::NEG_INFINITY;
            let mut prev: Option<f64> = None;
            for k in 0..=8 {
                let t = a0 + (a1 - a0) * k as f64 / 8.0;
                let uv = pc.d0(t);
                let mut u = uv.x();
                if u_period > 0.0 {
                    // Continuously unwrap relative to the previous sample so a
                    // seam-wrapping edge keeps its true span instead of
                    // collapsing modulo the period.
                    if let Some(p) = prev {
                        while u - p > u_period * 0.5 { u -= u_period; }
                        while u - p < -u_period * 0.5 { u += u_period; }
                    }
                    prev = Some(u);
                }
                cu_min = cu_min.min(u);
                cu_max = cu_max.max(u);
                vmin = vmin.min(uv.y());
                vmax = vmax.max(uv.y());
            }
            umin = umin.min(cu_min);
            umax = umax.max(cu_max);
            // A boundary edge spanning a full u-period (a full-circle / seam
            // loop) forces the face's u-domain to the whole period — otherwise
            // a full cylinder's u-extent collapses to [0, 7π/4] (the sampled
            // max never reaches the period boundary) and the mesh loses the
            // final π/4 wedge (~12.5% of the lateral area).
            if u_period > 0.0 && cu_max - cu_min >= u_period - 1e-6 {
                full_period = true;
            }
        }
    }
    if full_period {
        (su0, su0 + u_period, vmin, vmax)
    } else if umin.is_finite() && umax > umin && vmin.is_finite() && vmax > vmin {
        (umin, umax, vmin, vmax)
    } else {
        (0.0, 1.0, 0.0, 1.0)
    }
}

/// Invert a 3D point lying on a surface to its `(u, v)` parameters via Newton
/// iteration on the surface `d1` Gram matrix. Returns `None` if it fails to
/// converge (point off-surface, singular Jacobian).
fn invert_uv(surface: &dyn Surface, p: &GpPnt) -> Option<(f64, f64)> {
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let (mut u, mut v) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
    if !(u.is_finite() && v.is_finite()) {
        (u, v) = (0.0, 0.0);
    }
    for _ in 0..12 {
        let (sp, su, sv) = surface.d1(u, v);
        let r = p.coord.subtracted(&sp.coord);
        if r.modulus() < 1e-9 { return Some((u, v)); }
        let (gx, gy) = (su.coord, sv.coord);
        let g11 = gx.dot(&gx);
        let g12 = gx.dot(&gy);
        let g22 = gy.dot(&gy);
        let det = g11 * g22 - g12 * g12;
        if det.abs() < 1e-24 { return None; }
        let du = (g22 * gx.dot(&r) - g12 * gy.dot(&r)) / det;
        let dv = (-g12 * gx.dot(&r) + g11 * gy.dot(&r)) / det;
        if !du.is_finite() || !dv.is_finite() { return None; }
        u += du;
        v += dv;
        if du.abs() < 1e-12 && dv.abs() < 1e-12 { break; }
    }
    if p.distance(&surface.d0(u, v)) < 1e-6 { Some((u, v)) } else { None }
}

/// Triangle area by cross-product magnitude.
fn triangle_area(pts: &[GpPnt], t: &Triangle) -> f64 {
    let ab = pts[t.n1].coord.subtracted(&pts[t.n0].coord);
    let ac = pts[t.n2].coord.subtracted(&pts[t.n0].coord);
    0.5 * ab.crossed(&ac).modulus()
}

/// Wind the triangle so its normal agrees with the surface normal `ns`.
fn orient(pts: &[GpPnt], i0: usize, i1: usize, i2: usize, ns: &occt_core::gp::GpXyz) -> Triangle {
    let pa = pts[i0].coord;
    let pb = pts[i1].coord;
    let pc = pts[i2].coord;
    let n = pb.subtracted(&pa).crossed(&pc.subtracted(&pa));
    if ns.dot(&n) < 0.0 {
        Triangle::new(i0, i2, i1)
    } else {
        Triangle::new(i0, i1, i2)
    }
}

/// Distance from `p` to the segment `a..b`.
fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::{GpDir, GpLin, GpPln};
    use occt_geom::{GeomLine, GeomPlane};
    use crate::builder::TopoBuilder;
    use crate::shape::Wire;
    use crate::tgeometry::{EdgeGeom, FaceGeom};

    // ---- shared fixtures (also used by sibling modules' tests) ----

    /// Edge along `(x0,y0)` → `(x1,y1)` in the z=0 plane, parameterized [0,1].
    pub(crate) fn segment_edge(x0: f64, y0: f64, x1: f64, y1: f64) -> Edge {
        let e = Edge::new();
        let dir = GpDir::new(x1 - x0, y1 - y0, 0.0).unwrap();
        let lin = GpLin::from_pnt_dir(GpPnt::new(x0, y0, 0.0), dir);
        reg().set_edge(&e.0, EdgeGeom::new(Arc::new(GeomLine::new(lin)), 0.0, 1.0));
        e
    }

    /// Unit-square face in the z=0 plane (surface registered, 4-edge wire).
    pub(crate) fn square_face() -> Face {
        let b = TopoBuilder::new();
        let mut face = Face::new();
        reg().set_face(&face.0, FaceGeom::new(Arc::new(GeomPlane::new(GpPln::default()))));
        let mut wire = Wire::new();
        for (x0, y0, x1, y1) in [
            (0.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 1.0, 1.0),
            (1.0, 1.0, 0.0, 1.0),
            (0.0, 1.0, 0.0, 0.0),
        ] {
            let e = segment_edge(x0, y0, x1, y1);
            b.add_edge(&mut wire, &e);
        }
        b.add_wire(&mut face, &wire);
        face
    }

    /// A face with an outer 2×2 square and a 24-gon hole (a box-top ring cut by
    /// a boss cylinder) triangulates to the ring area, not the full square.
    #[test]
    fn face_with_hole_triangulates_ring_area() {
        use occt_core::gp::{GpAx3, GpPnt};
        let b = TopoBuilder::new();
        let mut face = Face::new();
        let pln = GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap());
        reg().set_face(&face.0, FaceGeom::new(Arc::new(GeomPlane::new(pln.clone()))));
        // Outer 2×2 square.
        let mut outer = Wire::new();
        for (x0, y0, x1, y1) in [(0.0, 0.0, 2.0, 0.0), (2.0, 0.0, 2.0, 2.0), (2.0, 2.0, 0.0, 2.0), (0.0, 2.0, 0.0, 0.0)] {
            let e = b.make_edge_segment(&GpPnt::new(x0, y0, 1.0), &GpPnt::new(x1, y1, 1.0));
            b.add_edge(&mut outer, &e);
        }
        b.add_wire(&mut face, &outer);
        // Hole: 24-gon at (1,1), r=0.25.
        let r = 0.25;
        let n = 24usize;
        let mut hole = Wire::new();
        for i in 0..n {
            let t0 = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            let t1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / n as f64;
            let e = b.make_edge_segment(
                &GpPnt::new(1.0 + r * t0.cos(), 1.0 + r * t0.sin(), 1.0),
                &GpPnt::new(1.0 + r * t1.cos(), 1.0 + r * t1.sin(), 1.0),
            );
            b.add_edge(&mut hole, &e);
        }
        b.add_wire(&mut face, &hole);
        let (pts, tris) = face_to_triangles(&face, 0.05);
        let mut area = 0.0;
        for t in &tris {
            let a = pts[t.n0].coord;
            let x = pts[t.n1].coord.subtracted(&a);
            let y = pts[t.n2].coord.subtracted(&a);
            area += 0.5 * x.crossed(&y).modulus();
        }
        let expect = 4.0 - std::f64::consts::PI * r * r;
        assert!(
            (area - expect).abs() < 0.05,
            "ring area {area} vs expected {expect} (hole not subtracted?)"
        );
    }

    #[test]
    fn edge_polyline_endpoints_match() {
        let e = segment_edge(0.0, 0.0, 1.0, 0.0);
        let poly = edge_to_polyline(&e, 0.01);
        assert!(poly.len() >= 2, "poly len {}", poly.len());
        assert!(poly.first().unwrap().distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9);
        assert!(poly.last().unwrap().distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn edge_chord_error_straight_line_is_zero() {
        let e = segment_edge(0.0, 0.0, 1.0, 0.0);
        let poly = edge_to_polyline(&e, 0.01);
        assert!(edge_chord_error(&e, &poly) < 1e-9, "err {}", edge_chord_error(&e, &poly));
    }

    #[test]
    fn face_triangulates_to_unit_area() {
        let f = square_face();
        let (pts, tris) = face_to_triangles(&f, 0.25);
        // Planar faces triangulate exactly from their boundary polygon: a
        // square is 2 triangles (not the dense UV grid).
        assert_eq!(tris.len(), 2, "tris {}", tris.len());
        let area: f64 = tris.iter().map(|t| triangle_area(&pts, t)).sum();
        assert!((area - 1.0).abs() < 1e-6, "area {area}");
    }

    #[test]
    fn face_area_of_unit_square() {
        let f = square_face();
        assert!((face_area(&f) - 1.0).abs() < 1e-6, "area {}", face_area(&f));
    }
}
