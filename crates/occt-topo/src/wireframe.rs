//! Wireframe & face tessellation.
//!
//! Port of the *meshing* half of `BRepMesh_IncrementalMesh`: an edge is
//! turned into a deflection-bounded polyline (`GCPnts_UniformDeflection`) and a
//! face into a UV-grid triangle soup. Geometry is read from the side-table
//! registry (`tgeometry::GeometryRegistry`) — the same data `BRep_Tool` reads.

use occt_core::gp::GpPnt;
use occt_core::gcpnts::{CurveSample, UniformDeflection, UniformPoints};
use occt_core::poly::triangulation::Triangle;
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

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

/// Tessellate a face into a UV-grid triangle soup (vertices + triangles).
///
/// Grid resolution is driven by `deflection` relative to the UV domain size
/// (`nu = clamp(ceil(du/deflection)+1, 3, 64)`). Triangle winding follows the
/// surface normal from `d1` so the resulting mesh points outward. Degenerate
/// cells are skipped. Missing surface or unbounded domain yields an empty soup.
pub fn face_to_triangles(f: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
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
            // Sample the pcurve along the edge (a full circle edge's endpoints
            // coincide in UV, so endpoints alone would collapse the domain).
            for k in 0..=8 {
                let t = a0 + (a1 - a0) * k as f64 / 8.0;
                let uv = pc.d0(t);
                let u0a = if u_period > 0.0 {
                    (uv.x() - su0).rem_euclid(u_period) + su0
                } else {
                    uv.x()
                };
                umin = umin.min(u0a);
                umax = umax.max(u0a);
                vmin = vmin.min(uv.y());
                vmax = vmax.max(uv.y());
            }
        }
    }
    if umin.is_finite() && umax > umin && vmin.is_finite() && vmax > vmin {
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
        assert!(tris.len() >= 6, "tris {}", tris.len());
        let area: f64 = tris.iter().map(|t| triangle_area(&pts, t)).sum();
        assert!((area - 1.0).abs() < 1e-6, "area {area}");
    }

    #[test]
    fn face_area_of_unit_square() {
        let f = square_face();
        assert!((face_area(&f) - 1.0).abs() < 1e-6, "area {}", face_area(&f));
    }
}
