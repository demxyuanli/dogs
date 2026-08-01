//! Analytic mass properties and analytic extrema — exact `BRepGProp` + `Extrema`.
//!
//! Two families of routines:
//!
//! * **Mass properties** — exact surface area, volume and centroid for shapes
//!   built from analytic surfaces (planes, spheres, cylinders, cones, toruses).
//!   A box reports `2(dx·dy + dy·dz + dx·dz)` and `dx·dy·dz` exactly, unlike the
//!   tessellation-based [`crate::brep_gprop`] which inflates curved faces. Faces
//!   whose surface is not one of the analytic types fall back to the mesh area.
//!
//! * **Extrema** — closed-form nearest-point solutions for point/line/plane/
//!   sphere/cylinder/cone/torus (no grid search). Unsupported surfaces fall back
//!   to [`occt_geom::extrema::point_surface_extrema`].
//!
//! `Surface` trait objects cannot be downcast (the port keeps no `Any`-typed
//! geometry), so analytic parameters (radius, axis, …) are reconstructed from
//! sampled geometry invariants rather than read from a concrete `Gp*` struct.

use std::f64::consts::PI;

use occt_core::geom::polygon_ops::polygon_area3d;
use occt_core::gp::{GpDir, GpPln, GpPnt, GpVec};
use occt_geom::extrema::{point_surface_extrema, ExtremaPair};
use occt_geom::{GeomCircle, GeomLine, GeomPlane, GeomSphere, Surface};
use occt_math::integrate;

use crate::brep_surface::{
    classify_surface, face_is_planar, face_plane, sphere_center, surface_normal, SurfaceKind,
};
use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::topo_tools_full::{
    edge_vertices, edges_of_wire, faces_of, vertex_position, vertices_of, wires_of_face,
};

/// Combined analytic mass properties of a closed solid.
#[derive(Debug, Clone)]
pub struct AnalyticProps {
    pub surface_area: f64,
    pub volume: f64,
    pub centroid: GpPnt,
    /// `true` when every face was handled analytically (no mesh fallback).
    pub exact: bool,
}

// ---------------------------------------------------------------------------
// Small vector helpers
// ---------------------------------------------------------------------------

fn pnt_add_vec(p: GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::from_xyz(&p.coord.added(&v.coord))
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new((a.x() + b.x()) / 2.0, (a.y() + b.y()) / 2.0, (a.z() + b.z()) / 2.0)
}

/// Distance from `p` to the infinite line through `c` along unit `ax`.
fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

/// A unit vector perpendicular to `v`.
fn perpendicular(v: &GpVec) -> GpVec {
    let a = GpVec::new(1.0, 0.0, 0.0);
    let b = GpVec::new(0.0, 1.0, 0.0);
    let cand = if v.cross_magnitude(&a) > 1e-9 { v.crossed(&a) } else { v.crossed(&b) };
    cand.normalized()
}

// ---------------------------------------------------------------------------
// Analytic surface classification (extends `brep_surface::classify_surface`)
// ---------------------------------------------------------------------------

/// Classify a surface, adding cylinder / torus / cone detection to the
/// plane / sphere detection already provided by `brep_surface`.
fn classify_surface_full(s: &dyn Surface) -> SurfaceKind {
    let base = classify_surface(s);
    if base == SurfaceKind::Plane || base == SurfaceKind::Sphere {
        return base;
    }
    if cylinder_params(s).is_some() {
        return SurfaceKind::Cylinder;
    }
    if torus_params(s).is_some() {
        return SurfaceKind::Torus;
    }
    if cone_params(s).is_some() {
        return SurfaceKind::Cone;
    }
    SurfaceKind::Other
}

/// Cylinder parameters from a `dyn Surface`:
/// `(center, unit axis, radius)`, reconstructed from sampled invariants.
fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    // The axis is along ∂S/∂v (constant for a cylinder): S(u,v+1) − S(u,v).
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axis.divided(m);
    // Opposite points on a ring give 2r; their midpoint is on the axis.
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let r = q0.distance(&q1) / 2.0;
    if r <= 1e-9 {
        return None;
    }
    let center = midpoint(&q0, &q1);
    // Verify every sampled point sits at distance r from the axis.
    let (nu, nv) = (8, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = (j as f64 - nv as f64 / 2.0) / 2.0;
            let p = s.d0(u, v);
            if (dist_to_axis(&p, &center, &ax) - r).abs() > 1e-3 * r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, r))
}

/// Torus parameters: `(center, axis, major radius, minor radius)`.
fn torus_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64, f64)> {
    // v = ±π/2 are the top/bottom of the tube; their difference is the axis.
    let top = s.d0(0.0, std::f64::consts::FRAC_PI_2);
    let bot = s.d0(0.0, -std::f64::consts::FRAC_PI_2);
    let axv = GpVec::from_pnts(&bot, &top);
    let m = axv.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axv.divided(m);
    // Outer equator (v=0): opposite points give the center and major radius.
    let e0 = s.d0(0.0, 0.0);
    let e1 = s.d0(PI, 0.0);
    let center = midpoint(&e0, &e1);
    let rel = GpVec::from_pnts(&center, &top);
    let rel_axis = rel.dot(&ax);
    let in_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let big_r = in_plane.magnitude();
    let small_r = rel_axis.abs();
    if big_r <= 1e-9 || small_r <= 1e-9 {
        return None;
    }
    // Verify sampled points lie on the tube: sqrt((ρ−R)² + z²) = r.
    let (nu, nv) = (6, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = 2.0 * PI * j as f64 / nv as f64;
            let p = s.d0(u, v);
            let rr = GpVec::from_pnts(&center, &p);
            let ra = rr.dot(&ax);
            let rho = rr.subtracted(&ax.multiplied_scalar(ra)).magnitude();
            let err = ((rho - big_r).powi(2) + ra.powi(2)).sqrt();
            if (err - small_r).abs() > 1e-3 * small_r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, big_r, small_r))
}

/// Closest point between two lines `p1 + t·d1` and `p2 + s·d2`.
/// `None` when the lines are parallel (denominator ~ 0) or do not meet within
/// tolerance.
fn closest_point_lines(p1: GpPnt, d1: GpVec, p2: GpPnt, d2: GpVec) -> Option<GpPnt> {
    // `from_pnts(p1, p2)` = p2 − p1, so the closest-point system is
    //   a·t − b·s = d1·(p2−p1),   b·t − c·s = d2·(p2−p1).
    let r = GpVec::from_pnts(&p1, &p2);
    let a = d1.dot(&d1);
    let b = d1.dot(&d2);
    let c = d2.dot(&d2);
    let e = d1.dot(&r);
    let f = d2.dot(&r);
    let denom = a * c - b * b;
    if denom.abs() < 1e-12 {
        return None; // parallel
    }
    let t = (e * c - b * f) / denom;
    let s = (b * e - a * f) / denom;
    let q1 = pnt_add_vec(p1, &d1.multiplied_scalar(t));
    let q2 = pnt_add_vec(p2, &d2.multiplied_scalar(s));
    if q1.distance(&q2) < 1e-6 * (1.0 + p1.distance(&p2).max(1.0)) {
        Some(q1)
    } else {
        None
    }
}

/// Cone parameters: `(apex, unit axis, semi-angle)`.
///
/// A cone's fixed-`u` isocurves are straight generatrices that all meet at the
/// apex, so the apex is the intersection of two of them; the axis is the line
/// through the apex and the midpoints of the circular cross-sections; the
/// semi-angle is the (constant) angle every surface point makes with the axis.
fn cone_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let g0a = s.d0(0.0, 0.0);
    let g0b = s.d0(0.0, 1.0);
    let d0 = GpVec::from_pnts(&g0a, &g0b);
    if d0.magnitude() < 1e-9 {
        return None;
    }
    let g1a = s.d0(PI / 2.0, 0.0);
    let g1b = s.d0(PI / 2.0, 1.0);
    let d1 = GpVec::from_pnts(&g1a, &g1b);
    if d1.magnitude() < 1e-9 {
        return None;
    }
    let apex = closest_point_lines(g0a, d0, g1a, d1)?;
    // Axis: midpoint of opposite points on a cross-section lies on the axis.
    let m0 = midpoint(&s.d0(0.0, 0.0), &s.d0(PI, 0.0));
    let axis = GpVec::from_pnts(&apex, &m0);
    let am = axis.magnitude();
    if am < 1e-9 {
        return None;
    }
    let ax = axis.divided(am);
    let alpha = GpVec::from_pnts(&apex, &g0a).angle(&ax);
    if alpha < 1e-6 || alpha > std::f64::consts::FRAC_PI_2 - 1e-6 {
        return None;
    }
    // Verify every sampled point makes the same angle with the axis.
    let (nu, nv) = (6, 5);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = -1.0 + 2.0 * j as f64 / nv as f64;
            let p = s.d0(u, v);
            let ang = GpVec::from_pnts(&apex, &p).angle(&ax);
            if (ang - alpha).abs() > 1e-3 * alpha.max(0.01) {
                return None;
            }
        }
    }
    Some((apex, ax, alpha))
}

// ---------------------------------------------------------------------------
// Face polygon + planar area
// ---------------------------------------------------------------------------

/// Ordered boundary points of a face's outer wire, following the wire's
/// traversal regardless of each edge's curve direction. The closing duplicate
/// is dropped.
fn polygon_chain(face: &Face) -> Vec<GpPnt> {
    let Some(w) = wires_of_face(face).first().cloned() else { return Vec::new() };
    let edges = edges_of_wire(&w);
    let mut segs: Vec<(GpPnt, GpPnt)> = Vec::new();
    for e in &edges {
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            segs.push((vertex_position(&va), vertex_position(&vb)));
        }
    }
    if segs.is_empty() {
        return Vec::new();
    }
    let mut pts = vec![segs[0].0, segs[0].1];
    let mut used = vec![false; segs.len()];
    used[0] = true;
    let mut current = segs[0].1;
    loop {
        let mut progressed = false;
        for i in 0..segs.len() {
            if used[i] {
                continue;
            }
            if segs[i].0.distance(&current) < 1e-9 {
                current = segs[i].1;
                pts.push(current);
                used[i] = true;
                progressed = true;
                break;
            } else if segs[i].1.distance(&current) < 1e-9 {
                current = segs[i].0;
                pts.push(current);
                used[i] = true;
                progressed = true;
                break;
            }
        }
        if !progressed {
            break;
        }
    }
    if pts.len() > 1 && pts.last().unwrap().distance(&pts[0]) < 1e-9 {
        pts.pop();
    }
    pts
}

/// Exact area of a planar face: the boundary integral
/// `area = ½ ∮ (u dv − v du)` over the wire's edge curves projected onto the
/// plane, with per-edge orientation taken from the wire traversal. Handles
/// polygonal faces (a box face → its quad area) and disk caps (a circle edge →
/// πr²) alike.
fn planar_face_area(face: &Face, pln: &GpPln) -> f64 {
    let Some(w) = wires_of_face(face).first().cloned() else { return 0.0 };
    let edges = edges_of_wire(&w);
    if edges.is_empty() {
        return 0.0;
    }
    let chain = polygon_chain(face);
    let n = chain.len();
    let xd = pln.pos.x_direction().xyz();
    let yd = pln.pos.y_direction().xyz();
    let loc = pln.location();

    let mut total = 0.0;
    let mut any = false;
    for e in &edges {
        let Some(curve) = BRepTool::edge_curve(e) else { continue };
        let (a, b) = BRepTool::edge_parameters(e);
        if !(a.is_finite() && b.is_finite() && b > a) {
            continue;
        }
        let pa = curve.d0(a);
        let pb = curve.d0(b);
        let sign = if n >= 3 {
            let ia = chain.iter().position(|q| q.distance(&pa) < 1e-6);
            let ib = chain.iter().position(|q| q.distance(&pb) < 1e-6);
            match (ia, ib) {
                (Some(i), Some(j)) => {
                    if (j + 1) % n == i {
                        -1.0
                    } else {
                        1.0
                    }
                }
                _ => 1.0,
            }
        } else {
            1.0
        };
        let f = |t: f64| {
            let p = curve.d0(t);
            let (_, tg) = curve.d1(t);
            let u = (p.x() - loc.x()) * xd.x + (p.y() - loc.y()) * xd.y + (p.z() - loc.z()) * xd.z;
            let v = (p.x() - loc.x()) * yd.x + (p.y() - loc.y()) * yd.y + (p.z() - loc.z()) * yd.z;
            let du = tg.x() * xd.x + tg.y() * xd.y + tg.z() * xd.z;
            let dv = tg.x() * yd.x + tg.y() * yd.y + tg.z() * yd.z;
            u * dv - v * du
        };
        total += sign * integrate(&f, a, b, 24);
        any = true;
    }
    if !any {
        return polygon_area3d(&chain);
    }
    0.5 * total.abs()
}

/// Extent of a face's boundary along a unit axis (used for the cylinder height).
fn face_axial_extent(face: &Face, ax: &GpVec) -> f64 {
    let mut tmin = f64::INFINITY;
    let mut tmax = f64::NEG_INFINITY;
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            let Some(curve) = BRepTool::edge_curve(&e) else { continue };
            let (a, b) = BRepTool::edge_parameters(&e);
            if !(a.is_finite() && b.is_finite() && b > a) {
                continue;
            }
            for i in 0..=8 {
                let p = curve.d0(a + (b - a) * i as f64 / 8.0);
                let t = p.x() * ax.x() + p.y() * ax.y() + p.z() * ax.z();
                tmin = tmin.min(t);
                tmax = tmax.max(t);
            }
        }
    }
    if tmin.is_finite() && tmax.is_finite() {
        tmax - tmin
    } else {
        0.0
    }
}

/// Extent of all a shape's vertices along a unit axis.
fn shape_axial_extent(shape: &TopoShape, ax: &GpVec) -> f64 {
    let mut tmin = f64::INFINITY;
    let mut tmax = f64::NEG_INFINITY;
    for v in vertices_of(shape) {
        let p = vertex_position(&v);
        let t = p.x() * ax.x() + p.y() * ax.y() + p.z() * ax.z();
        tmin = tmin.min(t);
        tmax = tmax.max(t);
    }
    if tmin.is_finite() && tmax.is_finite() {
        tmax - tmin
    } else {
        0.0
    }
}

/// For a cone-like solid, return `(base radius, height)` derived from the
/// planar base face's area and the vertex extent along the base normal.
fn cone_radius_height(shape: &TopoShape) -> Option<(f64, f64)> {
    let faces = faces_of(shape);
    let base = faces.iter().find(|f| face_is_planar(f))?;
    let pln = face_plane(base)?;
    let base_area = planar_face_area(base, &pln);
    if base_area <= 0.0 {
        return None;
    }
    let r = (base_area / PI).sqrt();
    let ax1 = pln.axis();
    let ax = ax1.direction().xyz();
    let h = shape_axial_extent(shape, &GpVec::from_xyz(&ax));
    if h <= 0.0 {
        return None;
    }
    Some((r, h))
}

// ---------------------------------------------------------------------------
// Surface area
// ---------------------------------------------------------------------------

/// Exact surface area of a shape built from analytic surfaces.
///
/// * Plane — area of the face's UV window (boundary integral);
/// * Sphere — `r² Δu (sin v₁ − sin v₀)`;
/// * Cylinder — `r Δu Δv` with the height read from the face boundary;
/// * Cone — `π r ℓ` (full revolution) with `ℓ` the slant length;
/// * Torus — `r (v₁−v₀)(R (u₁−u₀) + r (sin u₁ − sin u₀))`;
/// * anything else — mesh area ([`crate::brep_gprop::surface_area`]).
pub fn analytic_surface_area(shape: &TopoShape) -> Result<f64, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("analytic_surface_area: shape has no faces".into());
    }
    let cone_rh = cone_radius_height(shape);
    let mut total = 0.0;
    for f in &faces {
        let Some(surf) = BRepTool::face_surface(f) else { continue };
        let a = match classify_surface_full(surf.as_ref()) {
            SurfaceKind::Plane => match face_plane(f) {
                Some(pln) => planar_face_area(f, &pln),
                None => crate::brep_gprop::surface_area(&f.0, 0.05),
            },
            SurfaceKind::Sphere => {
                let center = sphere_center(surf.as_ref()).unwrap_or(GpPnt::zero());
                let r = surf.d0(0.0, 0.0).distance(&center);
                let (u0, u1) = surf.u_range();
                let (v0, v1) = surf.v_range();
                if u1.is_finite() && v1.is_finite() {
                    r * r * (u1 - u0) * (v1.sin() - v0.sin())
                } else {
                    crate::brep_gprop::surface_area(&f.0, 0.05)
                }
            }
            SurfaceKind::Cylinder => match cylinder_params(surf.as_ref()) {
                Some((_, ax, r)) => {
                    let (u0, u1) = surf.u_range();
                    let du = if u1.is_finite() { u1 - u0 } else { 2.0 * PI };
                    let dv = face_axial_extent(f, &ax);
                    r * du * dv
                }
                None => crate::brep_gprop::surface_area(&f.0, 0.05),
            },
            SurfaceKind::Cone => match &cone_rh {
                Some((r, h)) => {
                    let (u0, u1) = surf.u_range();
                    let du = if u1.is_finite() { u1 - u0 } else { 2.0 * PI };
                    let l = (r * r + h * h).sqrt();
                    0.5 * r * l * du
                }
                None => crate::brep_gprop::surface_area(&f.0, 0.05),
            },
            SurfaceKind::Torus => match torus_params(surf.as_ref()) {
                Some((_, _, big_r, small_r)) => {
                    let (u0, u1) = surf.u_range();
                    let (v0, v1) = surf.v_range();
                    if u1.is_finite() && v1.is_finite() {
                        small_r * (v1 - v0)
                            * (big_r * (u1 - u0) + small_r * (u1.sin() - u0.sin()))
                    } else {
                        crate::brep_gprop::surface_area(&f.0, 0.05)
                    }
                }
                None => crate::brep_gprop::surface_area(&f.0, 0.05),
            },
            SurfaceKind::Other => crate::brep_gprop::surface_area(&f.0, 0.05),
        };
        total += a;
    }
    Ok(total)
}

// ---------------------------------------------------------------------------
// Volume
// ---------------------------------------------------------------------------

/// Signed volume of a closed convex polyhedron by the divergence theorem:
/// each planar face is fan-triangulated and `Σ (1/6) a·(b×c)` accumulated.
fn polyhedron_volume(shape: &TopoShape) -> Option<f64> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return None;
    }
    let mut vol = 0.0;
    for f in &faces {
        let pln = face_plane(f)?;
        let pts = polygon_chain(f);
        if pts.len() < 3 {
            continue;
        }
        for i in 1..(pts.len() - 1) {
            let a = pts[0].coord;
            let b = pts[i].coord;
            let c = pts[i + 1].coord;
            vol += a.dot_cross(&b, &c);
        }
        let _ = pln;
    }
    if vol.abs() < 1e-12 {
        None
    } else {
        Some((vol / 6.0).abs())
    }
}

fn sphere_volume_from_shape(shape: &TopoShape) -> Option<f64> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            let center = sphere_center(surf.as_ref())?;
            let r = surf.d0(0.0, 0.0).distance(&center);
            if r > 0.0 {
                return Some(4.0 / 3.0 * PI * r.powi(3));
            }
        }
    }
    None
}

fn cylinder_volume_from_shape(shape: &TopoShape) -> Option<f64> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((_, ax, r)) = cylinder_params(surf.as_ref()) {
                let h = face_axial_extent(&f, &ax);
                if h > 0.0 {
                    return Some(PI * r * r * h);
                }
            }
        }
    }
    None
}

fn torus_volume_from_shape(shape: &TopoShape) -> Option<f64> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((_, _, big_r, small_r)) = torus_params(surf.as_ref()) {
                return Some(2.0 * PI * PI * big_r * small_r * small_r);
            }
        }
    }
    None
}

/// Exact volume of an analytic solid.
///
/// All-planar closed solids are integrated by the divergence theorem (a box →
/// `dx·dy·dz` exactly). Spheres, cylinders, cones and toruses use their closed
/// forms. Anything else returns `Err` rather than a wrong number.
pub fn analytic_volume(shape: &TopoShape) -> Result<f64, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("analytic_volume: shape has no faces".into());
    }
    if faces.iter().all(face_is_planar) {
        if let Some(v) = polyhedron_volume(shape) {
            return Ok(v);
        }
    }
    let kinds: Vec<SurfaceKind> = faces
        .iter()
        .filter_map(|f| BRepTool::face_surface(f))
        .map(|s| classify_surface_full(s.as_ref()))
        .collect();
    if kinds.is_empty() {
        return Err("analytic_volume: no face surfaces".into());
    }
    let all_in = |k: SurfaceKind| kinds.iter().all(|&x| x == SurfaceKind::Plane || x == k);
    if all_in(SurfaceKind::Sphere) && kinds.contains(&SurfaceKind::Sphere) {
        if let Some(v) = sphere_volume_from_shape(shape) {
            return Ok(v);
        }
    }
    if all_in(SurfaceKind::Cylinder) && kinds.contains(&SurfaceKind::Cylinder) {
        if let Some(v) = cylinder_volume_from_shape(shape) {
            return Ok(v);
        }
    }
    if all_in(SurfaceKind::Cone) && kinds.contains(&SurfaceKind::Cone) {
        if let Some((r, h)) = cone_radius_height(shape) {
            return Ok(PI * r * r * h / 3.0);
        }
    }
    if all_in(SurfaceKind::Torus) && kinds.contains(&SurfaceKind::Torus) {
        if let Some(v) = torus_volume_from_shape(shape) {
            return Ok(v);
        }
    }
    Err("analytic_volume: unsupported shape".into())
}

// ---------------------------------------------------------------------------
// Centroid
// ---------------------------------------------------------------------------

/// Exact centroid of an analytic solid: bbox center for boxes, surface-center
/// for spheres/toruses, mid-axis point for cylinders, ¾-height for cones.
/// Falls back to the mesh centroid otherwise.
pub fn analytic_centroid(shape: &TopoShape) -> Result<GpPnt, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("analytic_centroid: shape has no faces".into());
    }
    // All-planar solid (a box): center of the vertex bounding box.
    if faces.iter().all(face_is_planar) {
        let verts = vertices_of(shape);
        if verts.len() >= 4 {
            let mut mn = [f64::INFINITY; 3];
            let mut mx = [f64::NEG_INFINITY; 3];
            for v in &verts {
                let p = vertex_position(v);
                mn[0] = mn[0].min(p.x());
                mn[1] = mn[1].min(p.y());
                mn[2] = mn[2].min(p.z());
                mx[0] = mx[0].max(p.x());
                mx[1] = mx[1].max(p.y());
                mx[2] = mx[2].max(p.z());
            }
            if mx[0] > mn[0] && mx[1] > mn[1] && mx[2] > mn[2] {
                return Ok(GpPnt::new(
                    (mn[0] + mx[0]) / 2.0,
                    (mn[1] + mx[1]) / 2.0,
                    (mn[2] + mx[2]) / 2.0,
                ));
            }
        }
    }
    // Sphere → center.
    for f in &faces {
        if let Some(surf) = BRepTool::face_surface(f) {
            if classify_surface_full(surf.as_ref()) == SurfaceKind::Sphere {
                if let Some(c) = sphere_center(surf.as_ref()) {
                    return Ok(c);
                }
            }
        }
    }
    // Cylinder → mid-axis point.
    for f in &faces {
        if let Some(surf) = BRepTool::face_surface(f) {
            if let Some((center, ax, _)) = cylinder_params(surf.as_ref()) {
                let h = shape_axial_extent(shape, &ax);
                return Ok(pnt_add_vec(center, &ax.multiplied_scalar(h / 2.0)));
            }
        }
    }
    // Cone → ¾ height from the base toward the apex.
    if let Some((r, _h)) = cone_radius_height(shape) {
        if let Some(pln) = faces.iter().find(|f| face_is_planar(f)).and_then(face_plane) {
            let ax1 = pln.axis();
            let ax = GpVec::from_xyz(ax1.direction().xyz());
            let base_center = pln.location();
            // The apex is the vertex farthest from the base center along the axis.
            let mut t_extreme = 0.0f64;
            for v in vertices_of(shape) {
                let p = vertex_position(&v);
                let t = (p.x() - base_center.x()) * ax.x()
                    + (p.y() - base_center.y()) * ax.y()
                    + (p.z() - base_center.z()) * ax.z();
                if t.abs() > t_extreme.abs() {
                    t_extreme = t;
                }
            }
            if t_extreme.abs() > 1e-12 {
                let _ = r;
                return Ok(pnt_add_vec(base_center, &ax.multiplied_scalar(t_extreme * 0.75)));
            }
        }
    }
    // Torus → center.
    for f in &faces {
        if let Some(surf) = BRepTool::face_surface(f) {
            if let Some((center, _, _, _)) = torus_params(surf.as_ref()) {
                return Ok(center);
            }
        }
    }
    // General fallback: mesh centroid.
    Ok(crate::brep_gprop::centroid(shape, 0.05).unwrap_or(GpPnt::zero()))
}

/// Combined exact mass properties.
pub fn analytic_properties(shape: &TopoShape) -> Result<AnalyticProps, String> {
    let surface_area = analytic_surface_area(shape)?;
    let volume = analytic_volume(shape)?;
    let centroid = analytic_centroid(shape)?;
    let exact = is_analytic(shape) && volume > 0.0;
    Ok(AnalyticProps { surface_area, volume, centroid, exact })
}

/// Whether every face of `shape` classifies as an analytic surface
/// (plane / sphere / cylinder / cone / torus).
pub fn is_analytic(shape: &TopoShape) -> bool {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return false;
    }
    faces.iter().all(|f| match BRepTool::face_surface(f) {
        Some(surf) => classify_surface_full(surf.as_ref()) != SurfaceKind::Other,
        None => false,
    })
}

// ---------------------------------------------------------------------------
// Analytic extrema
// ---------------------------------------------------------------------------

/// Exact closest point from `p` to an infinite line.
pub fn extrema_point_line(p: GpPnt, line: &GeomLine) -> ExtremaPair {
    let lin = line.lin();
    let p0 = lin.location();
    let d = *lin.direction().xyz();
    let v = GpVec::from_pnts(&p0, &p);
    let t = v.dot(&GpVec::from_xyz(&d));
    let q = GpPnt::from_xyz(&p0.coord.added(&d.multiplied(t)));
    ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: t,
        v1: None,
        u2: t,
        v2: None,
    }
}

/// Exact closest point from `p` to a circle.
///
/// Projects `p` onto the circle's plane and solves the in-plane angle. When `p`
/// lies on the axis every circle point is equidistant; angle 0 is chosen.
pub fn extrema_point_circle(p: GpPnt, circle: &GeomCircle) -> ExtremaPair {
    let c = circle.circ();
    let center = c.location();
    let ax2 = c.position();
    let n = GpVec::from_xyz(ax2.direction().xyz());
    let xd = GpVec::from_xyz(ax2.x_direction().xyz());
    let yd = GpVec::from_xyz(ax2.y_direction().xyz());
    let r = c.radius();
    let rel = GpVec::from_pnts(&center, &p);
    let proj = rel.subtracted(&n.multiplied_scalar(rel.dot(&n)));
    let in_r = proj.magnitude();
    let angle = if in_r < 1e-12 {
        0.0
    } else {
        proj.dot(&yd).atan2(proj.dot(&xd))
    };
    let q = pnt_add_vec(
        center,
        &xd.multiplied_scalar(r * angle.cos()).add(&yd.multiplied_scalar(r * angle.sin())),
    );
    ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: angle,
        v1: None,
        u2: angle,
        v2: None,
    }
}

/// Exact signed distance from `p` to a plane (closest point along the normal).
pub fn extrema_point_plane(p: GpPnt, plane: &GeomPlane) -> ExtremaPair {
    let pl = plane.pln();
    let loc = pl.location();
    let n = GpVec::from_xyz(pl.axis().direction().xyz());
    let xd = GpVec::from_xyz(pl.pos.x_direction().xyz());
    let yd = GpVec::from_xyz(pl.pos.y_direction().xyz());
    let rel = GpVec::from_pnts(&loc, &p);
    let d = rel.dot(&n);
    let q = GpPnt::from_xyz(&p.coord.subtracted(&n.coord.multiplied(d)));
    let u = rel.dot(&xd);
    let v = rel.dot(&yd);
    ExtremaPair {
        p1: p,
        p2: q,
        distance: d.abs(),
        u1: u,
        v1: Some(v),
        u2: u,
        v2: Some(v),
    }
}

/// Exact closest point from `p` to a sphere surface.
pub fn extrema_point_sphere(p: GpPnt, sphere: &GeomSphere) -> ExtremaPair {
    let center = sphere_center(sphere).unwrap_or_else(|| sphere.d0(0.0, 0.0));
    let r = sphere.d0(0.0, 0.0).distance(&center);
    let v = GpVec::from_pnts(&center, &p);
    let q = if v.magnitude() > 1e-12 {
        pnt_add_vec(center, &v.normalized().multiplied_scalar(r))
    } else {
        GpPnt::new(center.x() + r, center.y(), center.z())
    };
    let dist = (v.magnitude() - r).abs();
    let qrel = GpVec::from_pnts(&center, &q);
    let u = qrel.y().atan2(qrel.x());
    let vv = (qrel.z() / r.max(1e-30)).clamp(-1.0, 1.0).asin();
    ExtremaPair {
        p1: p,
        p2: q,
        distance: dist,
        u1: u,
        v1: Some(vv),
        u2: u,
        v2: Some(vv),
    }
}

/// Exact closest points between two (possibly skew) lines. Parallel lines
/// return any point pair at the perpendicular distance.
pub fn extrema_line_line(a: &GeomLine, b: &GeomLine) -> ExtremaPair {
    let la = a.lin();
    let lb = b.lin();
    let p1 = la.location();
    let p2 = lb.location();
    let d1 = GpVec::from_xyz(la.direction().xyz());
    let d2 = GpVec::from_xyz(lb.direction().xyz());
    let w0 = GpVec::from_pnts(&p2, &p1); // p1 − p2
    let aa = d1.dot(&d1);
    let bb = d1.dot(&d2);
    let cc = d2.dot(&d2);
    let dd = d1.dot(&w0);
    let ee = d2.dot(&w0);
    let denom = aa * cc - bb * bb;
    let (t, s) = if denom.abs() < 1e-12 {
        // parallel: any t; s projects p1 onto line b.
        (0.0, if cc > 1e-12 { ee / cc } else { 0.0 })
    } else {
        ((bb * ee - cc * dd) / denom, (aa * ee - bb * dd) / denom)
    };
    let q1 = GpPnt::from_xyz(&p1.coord.added(&d1.coord.multiplied(t)));
    let q2 = GpPnt::from_xyz(&p2.coord.added(&d2.coord.multiplied(s)));
    ExtremaPair {
        p1: q1,
        p2: q2,
        distance: q1.distance(&q2),
        u1: t,
        v1: None,
        u2: s,
        v2: None,
    }
}

/// Exact minimum distance between a line and a sphere:
/// `max(0, dist(line, center) − r)`. `None` when the sphere geometry is not
/// resolvable.
pub fn extrema_line_sphere(l: &GeomLine, s: &GeomSphere) -> Option<ExtremaPair> {
    let center = sphere_center(s)?;
    let r = s.d0(0.0, 0.0).distance(&center);
    let lin = l.lin();
    let p0 = lin.location();
    let d = GpVec::from_xyz(lin.direction().xyz());
    let v = GpVec::from_pnts(&p0, &center);
    let t = v.dot(&d);
    let q_line = pnt_add_vec(p0, &d.multiplied_scalar(t));
    let w = GpVec::from_pnts(&q_line, &center); // center − q_line
    let wm = w.magnitude();
    let q_sphere = if wm > 1e-12 {
        pnt_add_vec(center, &w.divided(wm).multiplied_scalar(-r))
    } else {
        GpPnt::new(center.x() + r, center.y(), center.z())
    };
    let dist = (wm - r).max(0.0);
    Some(ExtremaPair {
        p1: q_line,
        p2: q_sphere,
        distance: dist,
        u1: t,
        v1: None,
        u2: 0.0,
        v2: None,
    })
}

/// Exact nearest point on a cylinder surface to `p`.
fn extrema_point_cylinder(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    let (center, ax, r) = cylinder_params(s).ok_or("extrema: cylinder params")?;
    let rel = GpVec::from_pnts(&center, &p);
    let rel_axis = rel.dot(&ax);
    let rel_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let in_r = rel_plane.magnitude();
    let s0 = s.d0(0.0, 0.0);
    let xd_v = GpVec::from_pnts(&center, &s0);
    let xd = if xd_v.magnitude() > 1e-12 {
        xd_v.divided(xd_v.magnitude())
    } else {
        perpendicular(&ax)
    };
    let yd = ax.crossed(&xd).normalized();
    let angle = if in_r > 1e-12 {
        rel_plane.dot(&yd).atan2(rel_plane.dot(&xd))
    } else {
        0.0
    };
    let q = pnt_add_vec(
        center,
        &xd.multiplied_scalar(r * angle.cos())
            .add(&yd.multiplied_scalar(r * angle.sin()))
            .add(&ax.multiplied_scalar(rel_axis)),
    );
    let u = angle.rem_euclid(2.0 * PI);
    Ok(ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: Some(rel_axis),
        u2: u,
        v2: Some(rel_axis),
    })
}

/// Exact nearest point on a cone surface to `p`, solved in the axial plane.
fn extrema_point_cone(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    let (apex, ax, alpha) = cone_params(s).ok_or("extrema: cone params")?;
    let rel = GpVec::from_pnts(&apex, &p);
    let z_p = rel.dot(&ax);
    let radial = rel.subtracted(&ax.multiplied_scalar(z_p));
    let r_p = radial.magnitude();
    let e_radial = if r_p > 1e-12 {
        radial.divided(r_p)
    } else {
        perpendicular(&ax)
    };
    let (sa, ca) = alpha.sin_cos();
    let t = r_p * sa + z_p * ca;
    let q = pnt_add_vec(
        apex,
        &ax.multiplied_scalar(t * ca).add(&e_radial.multiplied_scalar(t * sa)),
    );
    let u = if r_p > 1e-12 {
        let xd = e_radial;
        let yd = ax.crossed(&xd).normalized();
        let ux = e_radial.dot(&xd);
        let uy = e_radial.dot(&yd);
        uy.atan2(ux)
    } else {
        0.0
    };
    Ok(ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: Some(t),
        u2: u,
        v2: Some(t),
    })
}

/// Exact nearest point on a torus surface to `p`: project onto the center
/// circle, then onto the tube circle in the radial/axial plane.
fn extrema_point_torus(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    let (center, ax, big_r, small_r) = torus_params(s).ok_or("extrema: torus params")?;
    let rel = GpVec::from_pnts(&center, &p);
    let rel_axis = rel.dot(&ax);
    let rel_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let in_r = rel_plane.magnitude();
    let e_radial = if in_r > 1e-12 {
        rel_plane.divided(in_r)
    } else {
        perpendicular(&ax)
    };
    let q_center = pnt_add_vec(center, &e_radial.multiplied_scalar(big_r));
    let w = GpVec::from_pnts(&q_center, &p);
    let wm = w.magnitude();
    let q = if wm > 1e-12 {
        pnt_add_vec(q_center, &w.divided(wm).multiplied_scalar(small_r))
    } else {
        pnt_add_vec(q_center, &e_radial.multiplied_scalar(small_r))
    };
    let xd = {
        let s0 = s.d0(0.0, 0.0);
        let v = GpVec::from_pnts(&center, &s0);
        let vp = v.subtracted(&ax.multiplied_scalar(v.dot(&ax)));
        if vp.magnitude() > 1e-12 {
            vp.divided(vp.magnitude())
        } else {
            perpendicular(&ax)
        }
    };
    let yd = ax.crossed(&xd).normalized();
    let u = e_radial.dot(&yd).atan2(e_radial.dot(&xd));
    let w_dir = if wm > 1e-12 { w.divided(wm) } else { e_radial };
    let vv = w_dir.dot(&ax).atan2(w_dir.dot(&e_radial));
    Ok(ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: Some(vv),
        u2: u,
        v2: Some(vv),
    })
}

/// Exact nearest point from `p` to a surface, dispatching on the analytic
/// surface type. Unsupported surfaces fall back to the grid+descent
/// [`occt_geom::extrema::point_surface_extrema`].
pub fn extrema_point_surface_exact(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    match classify_surface_full(s) {
        SurfaceKind::Plane => {
            let origin = s.d0(0.0, 0.0);
            let n = surface_normal(s, 0.0, 0.0);
            if n.xyz().square_modulus() < 1e-30 {
                return Err("extrema: degenerate plane normal".into());
            }
            let z_axis = GpDir::new(0.0, 0.0, 1.0).map_err(|e| e.to_string())?;
            let xd_v = if n.is_parallel(&GpVec::from_xyz(z_axis.xyz())) {
                GpVec::new(1.0, 0.0, 0.0)
            } else {
                GpVec::from_xyz(z_axis.xyz())
            };
            let yd = n.crossed(&xd_v.normalized()).normalized();
            let rel = GpVec::from_pnts(&origin, &p);
            let d = rel.dot(&n);
            let q = GpPnt::from_xyz(&p.coord.subtracted(&n.coord.multiplied(d)));
            let u = rel.dot(&xd_v.normalized());
            let v = rel.dot(&yd);
            Ok(ExtremaPair {
                p1: p,
                p2: q,
                distance: d.abs(),
                u1: u,
                v1: Some(v),
                u2: u,
                v2: Some(v),
            })
        }
        SurfaceKind::Sphere => {
            let center = sphere_center(s).ok_or("extrema: sphere center")?;
            let r = s.d0(0.0, 0.0).distance(&center);
            let v = GpVec::from_pnts(&center, &p);
            let q = if v.magnitude() > 1e-12 {
                pnt_add_vec(center, &v.normalized().multiplied_scalar(r))
            } else {
                GpPnt::new(center.x() + r, center.y(), center.z())
            };
            let qrel = GpVec::from_pnts(&center, &q);
            let u = qrel.y().atan2(qrel.x());
            let vv = (qrel.z() / r.max(1e-30)).clamp(-1.0, 1.0).asin();
            Ok(ExtremaPair {
                p1: p,
                p2: q,
                distance: p.distance(&q),
                u1: u,
                v1: Some(vv),
                u2: u,
                v2: Some(vv),
            })
        }
        SurfaceKind::Cylinder => extrema_point_cylinder(s, p),
        SurfaceKind::Cone => extrema_point_cone(s, p),
        SurfaceKind::Torus => extrema_point_torus(s, p),
        SurfaceKind::Other => Ok(point_surface_extrema(s, &p)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpDir, GpLin, GpPln, GpSphere as GpSphereT};
    use occt_geom::GeomCylinder;
    use crate::primitives::{
        BRepPrimBox, BRepPrimCone, BRepPrimCylinder, BRepPrimSphere, BRepPrimTorus,
    };

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn line(p: GpPnt, d: GpDir) -> GeomLine {
        GeomLine::new(GpLin::from_pnt_dir(p, d))
    }

    #[test]
    fn box_area_exact() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let a = analytic_surface_area(&b.solid.0).unwrap();
        assert!(approx(a, 52.0, 1e-9), "area {a}");
    }

    #[test]
    fn box_volume_exact() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let v = analytic_volume(&b.solid.0).unwrap();
        assert!(approx(v, 24.0, 1e-9), "volume {v}");
    }

    #[test]
    fn sphere_area_volume() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let a = analytic_surface_area(&s.solid.0).unwrap();
        assert!(approx(a, 4.0 * PI * 4.0, 1e-6), "area {a}");
        let v = analytic_volume(&s.solid.0).unwrap();
        assert!(approx(v, 4.0 / 3.0 * PI * 8.0, 1e-6), "volume {v}");
    }

    #[test]
    fn cylinder_volume() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let v = analytic_volume(&c.solid.0).unwrap();
        assert!(approx(v, 3.0 * PI, 1e-6), "volume {v}");
    }

    #[test]
    fn cone_volume() {
        let c = BRepPrimCone::make_cone(2.0, 6.0);
        let v = analytic_volume(&c.solid.0).unwrap();
        assert!(approx(v, PI * 4.0 * 6.0 / 3.0, 1e-6), "volume {v}");
    }

    #[test]
    fn box_centroid() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let c = analytic_centroid(&b.solid.0).unwrap();
        assert!(c.distance(&GpPnt::new(1.0, 1.0, 1.0)) < 1e-9, "centroid {c:?}");
    }

    #[test]
    fn analytic_properties_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let p = analytic_properties(&b.solid.0).unwrap();
        assert!(p.exact);
        assert!(approx(p.volume, 24.0, 1e-9));
        assert!(approx(p.surface_area, 52.0, 1e-9));
    }

    #[test]
    fn extrema_point_line_exact() {
        let l = line(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e = extrema_point_line(GpPnt::new(3.0, 4.0, 0.0), &l);
        assert!(approx(e.distance, 4.0, 1e-12), "dist {}", e.distance);
        assert!(approx(e.p2.x(), 3.0, 1e-12));
    }

    #[test]
    fn extrema_point_circle_exact() {
        let c = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let e = extrema_point_circle(GpPnt::new(3.0, 0.0, 5.0), &c);
        let want = (29.0f64).sqrt();
        assert!(approx(e.distance, want, 1e-9), "dist {} want {want}", e.distance);
        assert!(approx(e.p2.x(), 1.0, 1e-9));
    }

    #[test]
    fn extrema_point_plane() {
        let pl = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let e = super::extrema_point_plane(GpPnt::new(1.0, 2.0, 5.0), &pl);
        assert!(approx(e.distance, 5.0, 1e-12), "dist {}", e.distance);
        assert!(approx(e.p2.z(), 0.0, 1e-12));
    }

    #[test]
    fn extrema_point_sphere() {
        let s = GeomSphere::new(GpSphereT::new(GpAx3::standard(), 1.0).unwrap());
        let e = super::extrema_point_sphere(GpPnt::new(3.0, 0.0, 0.0), &s);
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
        assert!(approx(e.p2.x(), 1.0, 1e-9));
    }

    #[test]
    fn extrema_line_line_skew() {
        let a = line(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let b = line(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let e = extrema_line_line(&a, &b);
        assert!(approx(e.distance, 3.0, 1e-9), "dist {}", e.distance);
    }

    #[test]
    fn extrema_line_sphere() {
        let s = GeomSphere::new(GpSphereT::new(GpAx3::standard(), 1.0).unwrap());
        // Line at y=3 → distance 2.
        let l = line(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e = super::extrema_line_sphere(&l, &s).expect("extrema");
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
        // Line through the center → distance 0.
        let l2 = line(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e2 = super::extrema_line_sphere(&l2, &s).expect("extrema");
        assert!(approx(e2.distance, 0.0, 1e-9), "dist {}", e2.distance);
    }

    #[test]
    fn extrema_point_cylinder() {
        let cyl = GeomCylinder::new(GpCylinder::new(GpAx3::standard(), 1.0).unwrap());
        let e = extrema_point_surface_exact(&cyl, GpPnt::new(3.0, 0.0, 4.0)).unwrap();
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
        assert!(approx(e.p2.z(), 4.0, 1e-9));
    }

    #[test]
    fn is_analytic_box_sphere() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(is_analytic(&b.solid.0));
        let s = BRepPrimSphere::make_sphere(1.0);
        assert!(is_analytic(&s.solid.0));
    }

    #[test]
    fn torus_area_known() {
        let t = BRepPrimTorus::make_torus(1.0, 1.0);
        let a = analytic_surface_area(&t.solid.0).unwrap();
        assert!(approx(a, 4.0 * PI * PI * 1.0 * 1.0, 1e-3), "area {a}");
    }

    #[test]
    fn extrema_point_surface_sphere() {
        let s = GeomSphere::new(GpSphereT::new(GpAx3::standard(), 1.0).unwrap());
        let e = extrema_point_surface_exact(&s, GpPnt::new(0.0, 3.0, 0.0)).unwrap();
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
    }

    #[test]
    fn torus_volume_known() {
        let t = BRepPrimTorus::make_torus(2.0, 1.0);
        let v = analytic_volume(&t.solid.0).unwrap();
        assert!(approx(v, 2.0 * PI * PI * 2.0 * 1.0, 1e-6), "volume {v}");
    }
}
