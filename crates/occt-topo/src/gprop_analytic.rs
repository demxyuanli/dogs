//! Analytic mass properties and analytic extrema — exact `BRepGProp` + `Extrema`.
//!
//! Three families of routines:
//!
//! * **Mass properties** — exact surface area, volume and centroid for shapes
//!   built from analytic surfaces (planes, spheres, cylinders, cones, toruses).
//!   A box reports `2(dx·dy + dy·dz + dx·dz)` and `dx·dy·dz` exactly, unlike the
//!   tessellation-based [`crate::brep_gprop`] which inflates curved faces. Faces
//!   whose surface is not one of the analytic types fall back to the mesh area.
//!
//! * **Higher-order properties** — the 3×3 inertia tensor about the centroid
//!   ([`inertia_tensor`]), the principal moments and principal axes
//!   ([`principal_inertia`]), the density-weighted mass ([`analytic_mass`]), and
//!   the exact analytic curve length ([`analytic_curve_length`]). Analytic
//!   solids (box, sphere, cylinder, cone, torus) are integrated in closed form;
//!   other solids fall back to a vertex-mass approximation. The full set is
//!   assembled by [`full_analytic_properties`].
//!
//! * **Extrema** — closed-form nearest-point solutions for point/line/plane/
//!   sphere/cylinder/cone/torus (no grid search). Unsupported surfaces fall back
//!   to [`occt_geom::extrema::point_surface_extrema`].
//!
//! `Surface` trait objects cannot be downcast (the port keeps no `Any`-typed
//! geometry), so analytic parameters (radius, axis, …) are reconstructed from
//! sampled geometry invariants rather than read from a concrete `Gp*` struct.
//!
//! The inertia tensor is reported about the shape's centroid. For non-analytic
//! solids the vertex-mass fallback distributes the total mass equally over the
//! vertices and shifts the origin tensor to the centroid by the parallel-axis
//! theorem — a coarse approximation, documented per call site.

use std::f64::consts::PI;

use occt_core::geom::polygon_ops::polygon_area3d;
use occt_core::gp::{GpDir, GpPln, GpPnt, GpVec};
use occt_geom::extrema::{point_surface_extrema, ExtremaPair};
use occt_geom::{GeomCircle, GeomLine, GeomPlane, GeomSphere, Surface};
use occt_math::eigen_ext::jacobi_eigen_symmetric;
use occt_math::{integrate, MathMatrix};

use crate::brep_surface::{
    classify_surface, face_is_planar, face_plane, sphere_center, surface_normal, SurfaceKind,
};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, vertex_position, vertices_of, wires_of_face,
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

/// Volume of a closed convex polyhedron by the divergence theorem: each planar
/// face is fan-triangulated and `Σ (1/6) a·(b×c)` accumulated.
///
/// The sign of each face's contribution is fixed so its polygon-chain normal
/// points away from the shape's interior (the average of the vertices). This
/// makes the result independent of how the wire chains were oriented when the
/// solid was built — a translated box gives the same volume as one at the
/// origin, whereas the raw signed sum is only correct when every face happens
/// to be wound consistently outward.
fn polyhedron_volume(shape: &TopoShape) -> Option<f64> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return None;
    }
    let verts = vertices_of(shape);
    if verts.len() < 4 {
        return None;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for v in &verts {
        let p = vertex_position(v);
        cx += p.x();
        cy += p.y();
        cz += p.z();
    }
    let n = verts.len() as f64;
    let center = GpVec::new(cx / n, cy / n, cz / n);

    let mut vol = 0.0;
    for f in &faces {
        let pts = polygon_chain(f);
        if pts.len() < 3 {
            continue;
        }
        let mut face_vol = 0.0;
        for i in 1..(pts.len() - 1) {
            let a = pts[0].coord;
            let b = pts[i].coord;
            let c = pts[i + 1].coord;
            face_vol += a.dot_cross(&b, &c);
        }
        // Polygon-chain area normal (right-hand rule) and the face centroid.
        let mut nrm = GpVec::zero();
        let mut fc = GpVec::zero();
        for i in 0..pts.len() {
            let j = (i + 1) % pts.len();
            nrm = nrm.added(&GpVec::from_xyz(&pts[i].coord.crossed(&pts[j].coord)));
            fc = fc.added(&GpVec::from_xyz(&pts[i].coord));
        }
        fc = fc.multiplied_scalar(1.0 / pts.len() as f64);
        // The face is outward when its chain normal points away from the
        // interior; flip the signed contribution otherwise.
        if nrm.dot(&fc.subtracted(&center)) < 0.0 {
            face_vol = -face_vol;
        }
        vol += face_vol;
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
// Inertia tensor + principal axes
// ---------------------------------------------------------------------------

/// Symmetric 3×3 inertia tensor `(Ixx, Iyy, Izz, Ixy, Ixz, Iyz)` about the
/// centroid (off-diagonal entries are symmetric: `Ixy = Iyx`, etc.).
#[derive(Debug, Clone, Copy)]
pub struct InertiaTensor {
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixy: f64,
    pub ixz: f64,
    pub iyz: f64,
}

impl InertiaTensor {
    /// The 3×3 symmetric matrix form `[[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]`.
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        [
            [self.ixx, self.ixy, self.ixz],
            [self.ixy, self.iyy, self.iyz],
            [self.ixz, self.iyz, self.izz],
        ]
    }

    /// The trace `Ixx + Iyy + Izz` (twice the sum of the squared radii of
    /// gyration about the coordinate axes).
    pub fn trace(&self) -> f64 {
        self.ixx + self.iyy + self.izz
    }

    /// The zero tensor.
    pub fn zero() -> InertiaTensor {
        InertiaTensor { ixx: 0.0, iyy: 0.0, izz: 0.0, ixy: 0.0, ixz: 0.0, iyz: 0.0 }
    }

    /// A pure diagonal tensor (a solid whose coordinate axes are already
    /// principal).
    pub fn diagonal(ixx: f64, iyy: f64, izz: f64) -> InertiaTensor {
        InertiaTensor { ixx, iyy, izz, ixy: 0.0, ixz: 0.0, iyz: 0.0 }
    }

    /// Build from a 3×3 matrix, symmetrizing the off-diagonal entries
    /// (`M[i][j]` and `M[j][i]` are averaged).
    pub fn from_matrix(m: &[[f64; 3]; 3]) -> InertiaTensor {
        InertiaTensor {
            ixx: m[0][0],
            iyy: m[1][1],
            izz: m[2][2],
            ixy: 0.5 * (m[0][1] + m[1][0]),
            ixz: 0.5 * (m[0][2] + m[2][0]),
            iyz: 0.5 * (m[1][2] + m[2][1]),
        }
    }

    /// Element-wise sum — combine the inertia tensors of disjoint pieces
    /// (both must be about the same point).
    pub fn sum(&self, other: &InertiaTensor) -> InertiaTensor {
        InertiaTensor {
            ixx: self.ixx + other.ixx,
            iyy: self.iyy + other.iyy,
            izz: self.izz + other.izz,
            ixy: self.ixy + other.ixy,
            ixz: self.ixz + other.ixz,
            iyz: self.iyz + other.iyz,
        }
    }

    /// Scale all moments by `factor` (e.g. a density ratio).
    pub fn scaled(&self, factor: f64) -> InertiaTensor {
        InertiaTensor {
            ixx: self.ixx * factor,
            iyy: self.iyy * factor,
            izz: self.izz * factor,
            ixy: self.ixy * factor,
            ixz: self.ixz * factor,
            iyz: self.iyz * factor,
        }
    }

    /// Shift the tensor to a parallel point displaced by `offset` via the
    /// parallel-axis theorem: `I' = I + m·(d²δ − d·dᵀ)`.
    pub fn translated(&self, mass: f64, offset: &GpVec) -> InertiaTensor {
        let (x, y, z) = (offset.x(), offset.y(), offset.z());
        let d2 = x * x + y * y + z * z;
        InertiaTensor {
            ixx: self.ixx + mass * (d2 - x * x),
            iyy: self.iyy + mass * (d2 - y * y),
            izz: self.izz + mass * (d2 - z * z),
            ixy: self.ixy + mass * (-x * y),
            ixz: self.ixz + mass * (-x * z),
            iyz: self.iyz + mass * (-y * z),
        }
    }

    /// Principal moments (eigenvalues, descending) and principal axes
    /// (eigenvectors) of this tensor, via the symmetric Jacobi solver.
    pub fn principal(&self) -> Result<(Vec<f64>, Vec<GpVec>), String> {
        let mat = self.matrix();
        let mut m = MathMatrix::new(1, 3, 1, 3);
        for i in 0..3 {
            for j in 0..3 {
                m.set_value(i + 1, j + 1, mat[i][j]);
            }
        }
        let (vals, vecs) = jacobi_eigen_symmetric(&m, 1e-9)?;
        let mut pairs: Vec<(f64, GpVec)> = (0..3)
            .map(|i| (vals[i], GpVec::new(vecs[i].value(1), vecs[i].value(2), vecs[i].value(3))))
            .collect();
        pairs.sort_by(|a, b| b.0.total_cmp(&a.0));
        Ok((
            pairs.iter().map(|(v, _)| *v).collect(),
            pairs.iter().map(|(_, v)| *v).collect(),
        ))
    }

    /// Radii of gyration about the three coordinate axes: `sqrt(Ixx/m)`,
    /// `sqrt(Iyy/m)`, `sqrt(Izz/m)`.
    pub fn radius_of_gyration(&self, mass: f64) -> (f64, f64, f64) {
        let m = mass.max(1e-30);
        ((self.ixx / m).sqrt(), (self.iyy / m).sqrt(), (self.izz / m).sqrt())
    }

    /// The moment of inertia about a unit direction `u`: `I = u·(I·u)`.
    pub fn moment_about_axis(&self, u: &GpVec) -> f64 {
        let (x, y, z) = (u.x(), u.y(), u.z());
        self.ixx * x * x + self.iyy * y * y + self.izz * z * z
            + 2.0 * (self.ixy * x * y + self.ixz * x * z + self.iyz * y * z)
    }

    /// The inertia matrix applied to a vector (matrix-vector product `I·v`).
    pub fn apply(&self, v: &GpVec) -> GpVec {
        let (x, y, z) = (v.x(), v.y(), v.z());
        GpVec::new(
            self.ixx * x + self.ixy * y + self.ixz * z,
            self.ixy * x + self.iyy * y + self.iyz * z,
            self.ixz * x + self.iyz * y + self.izz * z,
        )
    }

    /// Angular momentum of a rigid body with this inertia tensor and angular
    /// velocity `omega`: `L = I·ω`.
    pub fn angular_momentum(&self, omega: &GpVec) -> GpVec {
        self.apply(omega)
    }
}

/// Axis-aligned box extents `(dx, dy, dz)` from the vertex bounding box, for
/// a solid with the box's structural signature (8 vertices / 12 edges / 6 faces).
fn box_dimensions(shape: &TopoShape) -> Option<(f64, f64, f64)> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return None;
    }
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
        Some((mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]))
    } else {
        None
    }
}

/// Radius of a sphere primitive, from any spherical face.
fn sphere_radius_from_shape(shape: &TopoShape) -> Option<f64> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if classify_surface_full(surf.as_ref()) == SurfaceKind::Sphere {
                let center = sphere_center(surf.as_ref())?;
                let r = surf.d0(0.0, 0.0).distance(&center);
                if r > 0.0 {
                    return Some(r);
                }
            }
        }
    }
    None
}

/// Cylinder `(axis, radius, height)` from the lateral surface and the axial
/// extent of the whole shape.
fn cylinder_axis_radius_height(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((_, ax, r)) = cylinder_params(surf.as_ref()) {
                let h = shape_axial_extent(shape, &ax);
                if h > 0.0 {
                    return Some((ax, r, h));
                }
            }
        }
    }
    None
}

/// Cone `(axis, base radius, height)`; the axis is the planar base face normal.
fn cone_axis_radius_height(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    let (r, h) = cone_radius_height(shape)?;
    let faces = faces_of(shape);
    let base = faces.iter().find(|f| face_is_planar(f))?;
    let pln = face_plane(base)?;
    let ax = GpVec::from_xyz(pln.axis().direction().xyz());
    Some((ax, r, h))
}

/// Torus `(axis, major radius, minor radius)`.
fn torus_axis_radii(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((_, ax, big_r, small_r)) = torus_params(surf.as_ref()) {
                return Some((ax, big_r, small_r));
            }
        }
    }
    None
}

/// Axis-aligned box extents `(dx, dy, dz)` of an analytic box solid.
pub fn analytic_box_dimensions(shape: &TopoShape) -> Option<(f64, f64, f64)> {
    box_dimensions(shape)
}

/// Radius of an analytic sphere solid.
pub fn analytic_sphere_radius(shape: &TopoShape) -> Option<f64> {
    sphere_radius_from_shape(shape)
}

/// `(center, axis, radius, height)` of an analytic cylinder solid.
pub fn analytic_cylinder_parameters(shape: &TopoShape) -> Option<(GpPnt, GpVec, f64, f64)> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((center, ax, r)) = cylinder_params(surf.as_ref()) {
                let h = shape_axial_extent(shape, &ax);
                if h > 0.0 {
                    return Some((center, ax, r, h));
                }
            }
        }
    }
    None
}

/// `(axis, base radius, height)` of an analytic cone solid.
pub fn analytic_cone_parameters(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    cone_axis_radius_height(shape)
}

/// `(axis, major radius, minor radius)` of an analytic torus solid.
pub fn analytic_torus_parameters(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    torus_axis_radii(shape)
}

/// Build the global-frame inertia tensor from the moments about the symmetry
/// axis (`i_axial`) and about any axis perpendicular to it (`i_radial`):
/// `I = I_rad·δ + (I_ax − I_rad)·(a⊗a)` for a unit axis `a`.
fn tensor_from_axis(ax: &GpVec, i_radial: f64, i_axial: f64) -> InertiaTensor {
    let (x, y, z) = (ax.x(), ax.y(), ax.z());
    let d = i_axial - i_radial;
    InertiaTensor {
        ixx: i_radial + d * x * x,
        iyy: i_radial + d * y * y,
        izz: i_radial + d * z * z,
        ixy: d * x * y,
        ixz: d * x * z,
        iyz: d * y * z,
    }
}

/// Fallback inertia tensor: treat each vertex as carrying an equal share of the
/// total mass and accumulate `Σ mᵢ·(rᵢ²δ − rᵢ rᵢᵀ)` about the origin, then shift
/// to the centroid by the parallel-axis theorem. This is an approximation (the
/// exact analytic forms above are used whenever the shape classifies) and is
/// only a coarse estimate for non-analytic solids.
fn vertex_mass_tensor(shape: &TopoShape, density: f64, volume: f64) -> Option<InertiaTensor> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return None;
    }
    let m = volume * density;
    let m_i = m / verts.len() as f64;
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    let mut ixy = 0.0;
    let mut ixz = 0.0;
    let mut iyz = 0.0;
    for v in &verts {
        let p = vertex_position(v);
        let (x, y, z) = (p.x(), p.y(), p.z());
        let r2 = x * x + y * y + z * z;
        ixx += m_i * (r2 - x * x);
        iyy += m_i * (r2 - y * y);
        izz += m_i * (r2 - z * z);
        ixy += m_i * (-x * y);
        ixz += m_i * (-x * z);
        iyz += m_i * (-y * z);
    }
    // Parallel-axis theorem: I_cm = I_origin − m·(d²δ − d·dᵀ) with d = centroid.
    let c = analytic_centroid(shape).unwrap_or(GpPnt::zero());
    let (cx, cy, cz) = (c.x(), c.y(), c.z());
    let d2 = cx * cx + cy * cy + cz * cz;
    ixx -= m * (d2 - cx * cx);
    iyy -= m * (d2 - cy * cy);
    izz -= m * (d2 - cz * cz);
    ixy -= m * (-cx * cy);
    ixz -= m * (-cx * cz);
    iyz -= m * (-cy * cz);
    Some(InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz })
}

/// Whether `shape` is an axis-aligned box solid: all-planar with the box's
/// structural signature (8 vertices / 12 edges / 6 faces).
pub fn shape_is_box(shape: &TopoShape) -> bool {
    let faces = faces_of(shape);
    faces.len() == 6
        && vertices_of(shape).len() == 8
        && edges_of(shape).len() == 12
        && faces.iter().all(face_is_planar)
}

/// Analytic kind of a closed solid, matching the inertia integration paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidKind {
    Box,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

/// Classify a closed analytic solid by the analytic surface it is built from.
pub fn classify_solid(shape: &TopoShape) -> SolidKind {
    if shape_is_box(shape) {
        return SolidKind::Box;
    }
    if sphere_radius_from_shape(shape).is_some() {
        return SolidKind::Sphere;
    }
    if cylinder_axis_radius_height(shape).is_some() {
        return SolidKind::Cylinder;
    }
    if cone_axis_radius_height(shape).is_some() {
        return SolidKind::Cone;
    }
    if torus_axis_radii(shape).is_some() {
        return SolidKind::Torus;
    }
    SolidKind::Other
}

/// 3×3 inertia tensor of `shape` about its centroid, with uniform `density`.
///
/// Analytic solids are integrated in closed form:
///
/// * box `dx×dy×dz` → `Ixx = m/12·(dy²+dz²)`, …;
/// * sphere radius `r` → `Ixx=Iyy=Izz = (2/5)·m·r²`;
/// * cylinder `r×h` → `Izz = (1/2)·m·r²`, `Ixx=Iyy = m·(3r²+h²)/12`;
/// * cone `r×h` → `Izz = (3/10)·m·r²`, `Ixx=Iyy = (3/20)·m·r² + (3/80)·m·h²`;
/// * torus `R×r` → `Izz = m·(R² + 3r²/4)`, `Ixx=Iyy = m·(R²/2 + 5r²/8)`.
///
/// Everything else falls back to the vertex-mass approximation. All results are
/// about the shape's centroid (no further offset is applied).
pub fn inertia_tensor(shape: &TopoShape, density: f64) -> Result<InertiaTensor, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("inertia_tensor: shape has no faces".into());
    }
    if !density.is_finite() {
        return Err("inertia_tensor: density must be finite".into());
    }
    let volume = analytic_volume(shape)?;
    if volume <= 0.0 {
        return Err("inertia_tensor: shape has no volume".into());
    }
    let m = volume * density;

    // Box: axis-aligned rectangular solid.
    if shape_is_box(shape) {
        if let Some((dx, dy, dz)) = box_dimensions(shape) {
            return Ok(InertiaTensor {
                ixx: m / 12.0 * (dy * dy + dz * dz),
                iyy: m / 12.0 * (dx * dx + dz * dz),
                izz: m / 12.0 * (dx * dx + dy * dy),
                ixy: 0.0,
                ixz: 0.0,
                iyz: 0.0,
            });
        }
    }
    // Sphere.
    if let Some(r) = sphere_radius_from_shape(shape) {
        let s = (2.0 / 5.0) * m * r * r;
        return Ok(InertiaTensor { ixx: s, iyy: s, izz: s, ixy: 0.0, ixz: 0.0, iyz: 0.0 });
    }
    // Cylinder.
    if let Some((ax, r, h)) = cylinder_axis_radius_height(shape) {
        let i_axial = 0.5 * m * r * r;
        let i_radial = m * (3.0 * r * r + h * h) / 12.0;
        return Ok(tensor_from_axis(&ax, i_radial, i_axial));
    }
    // Cone.
    if let Some((ax, r, h)) = cone_axis_radius_height(shape) {
        let i_axial = (3.0 / 10.0) * m * r * r;
        let i_radial = (3.0 / 20.0) * m * r * r + (3.0 / 80.0) * m * h * h;
        return Ok(tensor_from_axis(&ax, i_radial, i_axial));
    }
    // Torus.
    if let Some((ax, big_r, small_r)) = torus_axis_radii(shape) {
        let i_axial = m * (big_r * big_r + 0.75 * small_r * small_r);
        let i_radial = m * (0.5 * big_r * big_r + (5.0 / 8.0) * small_r * small_r);
        return Ok(tensor_from_axis(&ax, i_radial, i_axial));
    }

    // Fallback: vertex-mass approximation.
    if let Some(t) = vertex_mass_tensor(shape, density, volume) {
        return Ok(t);
    }
    Err("inertia_tensor: unsupported shape".into())
}

/// The 3×3 symmetric inertia matrix of `shape` about its centroid:
/// `[[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]`.
pub fn inertia_matrix(shape: &TopoShape, density: f64) -> Result<[[f64; 3]; 3], String> {
    Ok(inertia_tensor(shape, density)?.matrix())
}

/// Principal moments of inertia (eigenvalues, sorted descending) and the
/// corresponding principal axes (eigenvectors) of the inertia matrix.
///
/// The eigensolver is the symmetric Jacobi rotation method.
pub fn principal_inertia(shape: &TopoShape, density: f64) -> Result<(Vec<f64>, Vec<GpVec>), String> {
    let mat = inertia_matrix(shape, density)?;
    let mut m = MathMatrix::new(1, 3, 1, 3);
    for i in 0..3 {
        for j in 0..3 {
            m.set_value(i + 1, j + 1, mat[i][j]);
        }
    }
    let (vals, vecs) = jacobi_eigen_symmetric(&m, 1e-9)?;
    let mut pairs: Vec<(f64, GpVec)> = (0..3)
        .map(|i| {
            let v = GpVec::new(vecs[i].value(1), vecs[i].value(2), vecs[i].value(3));
            (vals[i], v)
        })
        .collect();
    pairs.sort_by(|a, b| b.0.total_cmp(&a.0));
    let moments: Vec<f64> = pairs.iter().map(|(v, _)| *v).collect();
    let axes: Vec<GpVec> = pairs.iter().map(|(_, v)| *v).collect();
    Ok((moments, axes))
}

/// The 3×3 rotation matrix whose columns are the shape's principal axes
/// (an orthonormal frame for the principal coordinate system).
pub fn principal_axes_matrix(shape: &TopoShape, density: f64) -> Result<[[f64; 3]; 3], String> {
    let (_, axes) = principal_inertia(shape, density)?;
    Ok([
        [axes[0].x(), axes[1].x(), axes[2].x()],
        [axes[0].y(), axes[1].y(), axes[2].y()],
        [axes[0].z(), axes[1].z(), axes[2].z()],
    ])
}

/// Inertia tensor about an arbitrary point `p`, obtained from the centroid
/// tensor via the parallel-axis theorem: `I_p = I_cm + m·(d²δ − d·dᵀ)` with
/// `d = p − centroid`.
pub fn inertia_tensor_at(shape: &TopoShape, density: f64, p: &GpPnt) -> Result<InertiaTensor, String> {
    let t = inertia_tensor(shape, density)?;
    let volume = analytic_volume(shape)?;
    let mass = volume * density;
    let c = analytic_centroid(shape)?;
    let offset = GpVec::from_pnts(&c, p);
    Ok(t.translated(mass, &offset))
}

/// Inertia tensor of a set of disjoint solids about their combined centre of
/// mass.
///
/// Each solid's tensor is computed about its own centroid (via
/// [`inertia_tensor`]), shifted to the combined centre of mass by the
/// parallel-axis theorem, and summed element-wise. Returns the combined tensor
/// and the combined centroid.
pub fn inertia_tensor_composite(
    shapes: &[&TopoShape],
    density: f64,
) -> Result<(InertiaTensor, GpPnt), String> {
    if shapes.is_empty() {
        return Err("inertia_tensor_composite: no shapes".into());
    }
    let mut total_mass = 0.0;
    let mut sum_mc = occt_core::gp::GpXyz::zero();
    let mut entries: Vec<(InertiaTensor, f64, GpPnt)> = Vec::with_capacity(shapes.len());
    for s in shapes {
        let volume = analytic_volume(s)?;
        let mass = volume * density;
        let c = analytic_centroid(s)?;
        sum_mc = sum_mc.added(&c.coord.multiplied(mass));
        total_mass += mass;
        entries.push((inertia_tensor(s, density)?, mass, c));
    }
    if total_mass <= 0.0 {
        return Err("inertia_tensor_composite: no mass".into());
    }
    let center = GpPnt::from_xyz(&sum_mc.divided(total_mass));
    let mut acc = InertiaTensor::zero();
    for (t, mass, c) in &entries {
        let offset = GpVec::from_pnts(c, &center);
        acc = acc.sum(&t.translated(*mass, &offset));
    }
    Ok((acc, center))
}

/// Density-weighted mass of an analytic solid: `density × volume`.
pub fn analytic_mass(shape: &TopoShape, density: f64) -> Result<f64, String> {
    if !density.is_finite() || density < 0.0 {
        return Err("analytic_mass: density must be finite and non-negative".into());
    }
    Ok(analytic_volume(shape)? * density)
}

// ---------------------------------------------------------------------------
// Analytic curve length
// ---------------------------------------------------------------------------

/// Analytic kind of a 3D curve, reconstructed from sampled geometry invariants
/// (the port keeps no `Any`-typed geometry, so the concrete `Geom*` type cannot
/// be downcast).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    /// Straight segment: constant tangent direction, zero second derivative.
    Line,
    /// Circular arc: constant tangent speed and curvature, planar binormal.
    Circle,
    /// Bounded planar curve that is not a circle (ellipse, etc.).
    Ellipse,
    /// Anything else (B-spline, helix, …), handled by chord-length sampling.
    Other,
}

/// Classify a curve over the parameter window `[a, b]`:
///
/// * `Line` — the second derivative vanishes at every sample;
/// * `Circle` — constant tangent speed `r` and constant nonzero curvature with
///   a planar binormal direction, so the arc length is exactly `r·|b−a|`;
/// * `Ellipse` — the sampled points are coplanar but the tangent speed is not
///   constant (an ellipse's arc length needs an elliptic integral);
/// * `Other` — non-planar or unbounded curves.
pub fn classify_curve(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> CurveKind {
    if is_line_curve(curve, a, b) {
        CurveKind::Line
    } else if circle_radius(curve, a, b).is_some() {
        CurveKind::Circle
    } else if curve_is_planar(curve, a, b) {
        CurveKind::Ellipse
    } else {
        CurveKind::Other
    }
}

/// Whether the curve's tangent direction is constant over `[a, b]` (a straight
/// line): every sampled second derivative is ~0.
fn is_line_curve(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> bool {
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = curve.d2(t);
        if d2.square_magnitude() > 1e-18 {
            return false;
        }
    }
    true
}

/// Radius of a circular arc over `[a, b]`, `None` when the sampled points are
/// not on a circle (constant tangent speed and curvature, planar binormal).
/// The tangent speed of a `GeomCircle` equals its radius, so the arc length is
/// `radius · |b − a|`.
fn circle_radius(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> Option<f64> {
    let mut speeds = Vec::with_capacity(6);
    let mut accs = Vec::with_capacity(6);
    let mut binorms: Vec<GpVec> = Vec::with_capacity(6);
    for i in 0..=5 {
        let t = a + (b - a) * i as f64 / 5.0;
        let (_, d1, d2) = curve.d2(t);
        let s = d1.magnitude();
        if s < 1e-12 {
            return None;
        }
        speeds.push(s);
        accs.push(d2.magnitude());
        binorms.push(d1.crossed(&d2).normalized());
    }
    let s0 = speeds[0];
    if speeds.iter().any(|&s| (s - s0).abs() > 1e-6 * s0.max(1.0)) {
        return None;
    }
    let a0 = accs[0];
    if a0 < 1e-12 || accs.iter().any(|&a| (a - a0).abs() > 1e-6 * a0.max(1.0)) {
        return None;
    }
    for bv in binorms.iter().skip(1) {
        if bv.cross_magnitude(&binorms[0]) > 1e-6 {
            return None;
        }
    }
    Some(s0)
}

/// Whether sampled points of the curve over `[a, b]` lie in a common plane.
fn curve_is_planar(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> bool {
    let n = 6;
    let pts: Vec<GpPnt> = (0..=n).map(|i| curve.d0(a + (b - a) * i as f64 / n as f64)).collect();
    if pts.len() < 3 {
        return false;
    }
    let v1 = GpVec::from_pnts(&pts[0], &pts[1]);
    let v2 = GpVec::from_pnts(&pts[0], &pts[pts.len() / 2]);
    let nrm = v1.crossed(&v2);
    let nm = nrm.magnitude();
    if nm < 1e-12 {
        return false;
    }
    pts.iter().all(|p| {
        GpVec::from_pnts(&pts[0], p).dot(&nrm).abs() <= 1e-6 * nm.max(1.0)
    })
}

/// Perimeter of a full ellipse with semi-axes `a` and `b` by the Ramanujan
/// approximation (the documented closed-form alternative to the arc-length
/// quadrature used for partial arcs):
/// `π·(3(a+b) − √((3a+b)·(a+3b)))`.
///
/// Exact for circles and accurate to roughly `1e-7` of the true perimeter for
/// any aspect ratio.
pub fn ellipse_length_ramanujan(a: f64, b: f64) -> f64 {
    let s = 3.0 * (a + b) - ((3.0 * a + b) * (a + 3.0 * b)).sqrt();
    PI * s
}

/// Exact length of one edge's curve.
///
/// * line — distance between the endpoints;
/// * circle — `radius · |b − a|` (the `GeomCircle` parameter is in radians);
/// * ellipse / other bounded planar curve — Gauss-Legendre quadrature of the
///   tangent speed (the elliptic integral of the second kind);
/// * anything else (B-spline, helix, …) — chord-length sampling via
///   [`crate::brep_measure::edge_length`].
fn exact_edge_length(e: &Edge) -> f64 {
    let Some(curve) = BRepTool::edge_curve(e) else { return 0.0 };
    let (a, b) = BRepTool::edge_parameters(e);
    if !(a.is_finite() && b.is_finite() && b > a) {
        return 0.0;
    }
    match classify_curve(curve.as_ref(), a, b) {
        CurveKind::Line => curve.d0(a).distance(&curve.d0(b)),
        CurveKind::Circle => match circle_radius(curve.as_ref(), a, b) {
            Some(r) => r * (b - a),
            None => crate::brep_measure::edge_length(e, 32),
        },
        CurveKind::Ellipse => {
            let speed = |t: f64| curve.d1(t).1.magnitude();
            integrate(&speed, a, b, 24)
        }
        CurveKind::Other => crate::brep_measure::edge_length(e, 32),
    }
}

/// Exact arc length of a single edge: line/circle are computed in closed form,
/// an ellipse by Gauss-Legendre quadrature, and anything else by chord-length
/// sampling. See [`classify_curve`] for the classification rules.
pub fn edge_length_exact(e: &Edge) -> f64 {
    exact_edge_length(e)
}

/// Total length of all `shape`'s edges, exact for analytic curves (a box's 12
/// edges sum to `4·(dx+dy+dz)`, a cylinder's cap circles to `2πr` each, …).
pub fn analytic_curve_length(shape: &TopoShape) -> Result<f64, String> {
    let edges = edges_of(shape);
    if edges.is_empty() {
        return Err("analytic_curve_length: shape has no edges".into());
    }
    Ok(edges.iter().map(edge_length_exact).sum())
}

/// Total length of all `shape`'s edges; `0.0` when the shape has no edges
/// (unlike [`analytic_curve_length`], which errors).
pub fn analytic_perimeter(shape: &TopoShape) -> f64 {
    edges_of(shape).iter().map(edge_length_exact).sum()
}

/// Exact perimeter of a face's outer wire (sum of its edges' exact lengths).
pub fn face_perimeter_exact(face: &Face) -> f64 {
    let Some(w) = wires_of_face(face).first().cloned() else { return 0.0 };
    edges_of_wire(&w).iter().map(edge_length_exact).sum()
}

// ---------------------------------------------------------------------------
// Combined analytic mass properties
// ---------------------------------------------------------------------------

/// Full analytic mass properties of a solid: mass, volume, surface area,
/// centroid, inertia tensor, principal moments and principal axes.
#[derive(Debug, Clone)]
pub struct FullAnalyticProps {
    pub mass: f64,
    pub volume: f64,
    pub surface_area: f64,
    pub centroid: GpPnt,
    pub inertia: InertiaTensor,
    pub principal_moments: Vec<f64>,
    pub principal_axes: Vec<GpVec>,
    /// `true` when every face was handled analytically (no mesh fallback).
    pub exact: bool,
}

impl FullAnalyticProps {
    /// The inertia tensor about an arbitrary point `p`, shifted from the
    /// centroid by the parallel-axis theorem.
    pub fn inertia_about(&self, p: &GpPnt) -> InertiaTensor {
        let offset = GpVec::from_pnts(&self.centroid, p);
        self.inertia.translated(self.mass, &offset)
    }

    /// Radii of gyration about the principal axes: `sqrt(I/m)` per principal
    /// moment (descending order matches [`Self::principal_moments`]).
    pub fn radius_of_gyration(&self) -> Vec<f64> {
        let m = self.mass.max(1e-30);
        self.principal_moments.iter().map(|i| (i / m).sqrt()).collect()
    }
}

/// One-stop analytic mass properties: exact area/volume/centroid plus the
/// inertia tensor and principal moments/axes about the centroid.
pub fn full_analytic_properties(shape: &TopoShape, density: f64) -> Result<FullAnalyticProps, String> {
    let props = analytic_properties(shape)?;
    let inertia = inertia_tensor(shape, density)?;
    let (principal_moments, principal_axes) = principal_inertia(shape, density)?;
    Ok(FullAnalyticProps {
        mass: props.volume * density,
        volume: props.volume,
        surface_area: props.surface_area,
        centroid: props.centroid,
        inertia,
        principal_moments,
        principal_axes,
        exact: props.exact,
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

    #[test]
    fn box_inertia_diagonal() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0); // dx=2, dy=3, dz=4
        let t = inertia_tensor(&b.solid.0, 1.0).unwrap();
        let m = 24.0;
        // Ixx about the x-axis = m/12·(dy²+dz²) = 2·(9+16) = 50.
        assert!(approx(t.ixx, m / 12.0 * (3.0 * 3.0 + 4.0 * 4.0), 1e-9), "ixx {}", t.ixx);
        assert!(approx(t.iyy, m / 12.0 * (2.0 * 2.0 + 4.0 * 4.0), 1e-9), "iyy {}", t.iyy);
        assert!(approx(t.izz, m / 12.0 * (2.0 * 2.0 + 3.0 * 3.0), 1e-9), "izz {}", t.izz);
        assert!(t.ixy.abs() < 1e-12 && t.ixz.abs() < 1e-12 && t.iyz.abs() < 1e-12);
        // inertia_matrix mirrors the tensor.
        let mat = inertia_matrix(&b.solid.0, 1.0).unwrap();
        assert!(approx(mat[0][0], t.ixx, 1e-12) && approx(mat[1][1], t.iyy, 1e-12) && approx(mat[2][2], t.izz, 1e-12));
        assert!(approx(mat[0][1], t.ixy, 1e-12) && approx(mat[0][2], t.ixz, 1e-12) && approx(mat[1][2], t.iyz, 1e-12));
    }

    #[test]
    fn sphere_inertia() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let t = inertia_tensor(&s.solid.0, 1.0).unwrap();
        let m = 4.0 / 3.0 * PI * 8.0;
        let want = (2.0 / 5.0) * m * 4.0;
        assert!(approx(t.ixx, want, 1e-6), "ixx {} want {}", t.ixx, want);
        assert!(approx(t.iyy, want, 1e-6));
        assert!(approx(t.izz, want, 1e-6));
        assert!(t.ixy.abs() < 1e-9);
    }

    #[test]
    fn cylinder_inertia() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let t = inertia_tensor(&c.solid.0, 1.0).unwrap();
        let m = PI * 3.0;
        assert!(approx(t.izz, 0.5 * m * 1.0, 1e-6), "izz {}", t.izz);
        let ix_want = m * (3.0 * 1.0 + 9.0) / 12.0;
        assert!(approx(t.ixx, ix_want, 1e-6), "ixx {} want {}", t.ixx, ix_want);
        assert!(approx(t.iyy, ix_want, 1e-6), "iyy {}", t.iyy);
        assert!(t.ixy.abs() < 1e-9 && t.ixz.abs() < 1e-9 && t.iyz.abs() < 1e-9);
    }

    #[test]
    fn principal_inertia_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (moments, _) = principal_inertia(&b.solid.0, 1.0).unwrap();
        assert_eq!(moments.len(), 3);
        // Sorted descending: largest moment about the smallest cross-section
        // axis (the 2 side → Ixx = 50), smallest about the largest (4 → Izz = 26).
        assert!(moments[0] >= moments[1] && moments[1] >= moments[2]);
        assert!(approx(moments[0], 50.0, 1e-6), "largest {}", moments[0]);
        assert!(approx(moments[1], 40.0, 1e-6), "middle {}", moments[1]);
        assert!(approx(moments[2], 26.0, 1e-6), "smallest {}", moments[2]);
    }

    #[test]
    fn principal_axes_orthogonal() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (_, axes) = principal_inertia(&b.solid.0, 1.0).unwrap();
        assert_eq!(axes.len(), 3);
        for i in 0..3 {
            assert!((axes[i].magnitude() - 1.0).abs() < 1e-6, "axis {i} not unit");
            for j in (i + 1)..3 {
                assert!(axes[i].dot(&axes[j]).abs() < 1e-6, "dot {i},{j} = {}", axes[i].dot(&axes[j]));
            }
        }
    }

    /// Build a properly-wired axis-aligned box solid with one corner at `p0`
    /// and extents `dx × dy × dz` (mirrors `BRepPrimBox` but at an arbitrary
    /// origin), so its centroid is not the global origin.
    fn box_solid_at(p0: &GpPnt, dx: f64, dy: f64, dz: f64) -> crate::shape::Solid {
        use std::sync::Arc;
        let b = crate::builder::TopoBuilder::new();
        let corners = [
            GpPnt::new(p0.x(), p0.y(), p0.z()),
            GpPnt::new(p0.x() + dx, p0.y(), p0.z()),
            GpPnt::new(p0.x() + dx, p0.y() + dy, p0.z()),
            GpPnt::new(p0.x(), p0.y() + dy, p0.z()),
            GpPnt::new(p0.x(), p0.y(), p0.z() + dz),
            GpPnt::new(p0.x() + dx, p0.y(), p0.z() + dz),
            GpPnt::new(p0.x() + dx, p0.y() + dy, p0.z() + dz),
            GpPnt::new(p0.x(), p0.y() + dy, p0.z() + dz),
        ];
        let verts: Vec<crate::shape::Vertex> =
            corners.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let box_edges: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges: Vec<crate::shape::Edge> = Vec::new();
        for &(i, j) in &box_edges {
            let p1 = corners[i];
            let p2 = corners[j];
            let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).expect("box edge direction");
            let mut e = b.make_edge(
                Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir))),
                0.0,
                p1.distance(&p2),
            );
            b.add(&mut e.0, &verts[i].0);
            b.add(&mut e.0, &verts[j].0);
            edges.push(e);
        }
        let edge_index = |a: usize, b: usize| {
            box_edges
                .iter()
                .position(|&(i, j)| (i == a && j == b) || (i == b && j == a))
                .expect("box edge")
        };
        let plane = |origin: GpPnt, normal: GpDir| {
            let z_axis = GpDir::new(0.0, 0.0, 1.0).unwrap();
            let xd = if normal.is_normal(&z_axis) { z_axis } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
            GpPln::new(GpAx3::new(origin, normal, &xd).unwrap())
        };
        let faces_def: [(GpPnt, GpDir, [usize; 4]); 6] = [
            (corners[0], GpDir::new(0.0, 0.0, -1.0).unwrap(), [0, 3, 2, 1]),
            (corners[4], GpDir::new(0.0, 0.0, 1.0).unwrap(), [4, 5, 6, 7]),
            (corners[0], GpDir::new(0.0, -1.0, 0.0).unwrap(), [0, 1, 5, 4]),
            (corners[3], GpDir::new(0.0, 1.0, 0.0).unwrap(), [3, 2, 6, 7]),
            (corners[0], GpDir::new(-1.0, 0.0, 0.0).unwrap(), [0, 4, 7, 3]),
            (corners[1], GpDir::new(1.0, 0.0, 0.0).unwrap(), [1, 2, 6, 5]),
        ];
        let mut faces = Vec::new();
        for (origin, normal, cycle) in faces_def {
            let quad: Vec<crate::shape::Edge> = [0, 1, 2, 3]
                .iter()
                .map(|&k| edges[edge_index(cycle[k], cycle[(k + 1) % 4])].clone())
                .collect();
            let wire = b.make_wire(&quad);
            let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane(origin, normal)));
            faces.push(b.make_face(surface, &[wire]));
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    #[test]
    fn inertia_about_centroid() {
        // A box spanning [0,dx]×[0,dy]×[0,dz] and the same box translated: the
        // tensor about each box's own centroid is identical (the analytic
        // formulas are inherently about the centroid, so no origin offset leaks
        // in and no parallel-axis term needs to be subtracted).
        let a = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let b = box_solid_at(&GpPnt::new(1.0, -2.0, 5.0), 2.0, 3.0, 4.0);
        let ta = inertia_tensor(&a.solid.0, 1.0).unwrap();
        let tb = inertia_tensor(&b.0, 1.0).unwrap();
        assert!(approx(ta.ixx, tb.ixx, 1e-9), "ixx {} vs {}", ta.ixx, tb.ixx);
        assert!(approx(ta.iyy, tb.iyy, 1e-9));
        assert!(approx(ta.izz, tb.izz, 1e-9));
        assert!(ta.ixy.abs() < 1e-9 && tb.ixy.abs() < 1e-9);
        // The value is the centroid tensor (Ixx = m/12(dy²+dz²) = 50). Had the
        // tensor been about the global origin the Ixx would be 200.
        assert!(approx(ta.ixx, 50.0, 1e-9), "centroid Ixx {}", ta.ixx);
    }

    #[test]
    fn analytic_curve_length_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let l = analytic_curve_length(&b.solid.0).unwrap();
        assert!(approx(l, 4.0 * (2.0 + 3.0 + 4.0), 1e-9), "length {l}");
    }

    #[test]
    fn circle_length() {
        // A cylinder: two full cap circles (r=1) + one seam line (h=1).
        let c = BRepPrimCylinder::make_cylinder(1.0, 1.0);
        let l = analytic_curve_length(&c.solid.0).unwrap();
        let want = 2.0 * 2.0 * PI * 1.0 + 1.0;
        assert!(approx(l, want, 1e-6), "length {l} want {want}");

        // Half circle edge: params [0, π] on a unit circle → π.
        let b = crate::builder::TopoBuilder::new();
        let half = b.make_edge_circle(&GpAx2::standard(), 1.0, 0.0, PI);
        let len = exact_edge_length(&half);
        assert!(approx(len, PI, 1e-9), "half circle {len}");
    }

    #[test]
    fn line_length() {
        let b = crate::builder::TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(3.0, 4.0, 0.0));
        let len = exact_edge_length(&e);
        assert!(approx(len, 5.0, 1e-9), "line length {len}");
    }

    #[test]
    fn full_analytic_props_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let p = full_analytic_properties(&b.solid.0, 1.0).unwrap();
        assert!(approx(p.mass, 24.0, 1e-9), "mass {}", p.mass);
        assert!(approx(p.volume, 24.0, 1e-9));
        assert!(approx(p.surface_area, 52.0, 1e-9));
        assert!(p.exact);
        assert_eq!(p.principal_moments.len(), 3);
        assert_eq!(p.principal_axes.len(), 3);
        assert!(approx(p.inertia.ixx, 50.0, 1e-9));
        assert!(approx(p.centroid.x(), 1.0, 1e-9));
        assert!(approx(p.centroid.y(), 1.5, 1e-9));
        assert!(approx(p.centroid.z(), 2.0, 1e-9));
    }

    #[test]
    fn mass_density() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let m = analytic_mass(&b.solid.0, 2.0).unwrap();
        assert!(approx(m, 2.0, 1e-9), "mass {m}");
    }

    #[test]
    fn moment_about_axis_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let t = inertia_tensor(&b.solid.0, 1.0).unwrap();
        // Moment about each coordinate axis is the corresponding diagonal entry.
        assert!(approx(t.moment_about_axis(&GpVec::new(1.0, 0.0, 0.0)), 50.0, 1e-9));
        assert!(approx(t.moment_about_axis(&GpVec::new(0.0, 1.0, 0.0)), 40.0, 1e-9));
        assert!(approx(t.moment_about_axis(&GpVec::new(0.0, 0.0, 1.0)), 26.0, 1e-9));
        // The axis (1,1,1)/√3 is equally weighted: I = (Ixx+Iyy+Izz)/3.
        let u = GpVec::new(1.0, 1.0, 1.0).normalized();
        assert!(approx(t.moment_about_axis(&u), (50.0 + 40.0 + 26.0) / 3.0, 1e-9));
        // Inertia matrix applied to a unit vector recovers the same value.
        let v = t.apply(&u);
        assert!(approx(u.dot(&v), (50.0 + 40.0 + 26.0) / 3.0, 1e-9));
    }

    #[test]
    fn inertia_tensor_at_origin() {
        // The tensor about the centroid, shifted to the origin by the
        // parallel-axis theorem, matches the closed-form origin tensor.
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let origin = GpPnt::zero();
        let t_origin = inertia_tensor_at(&b.solid.0, 1.0, &origin).unwrap();
        // Box 2×3×4 spanning [0,2]×[0,3]×[0,4], centroid (1,1.5,2), mass 24.
        // Ixx_origin = Ixx_cm + m(d²−dx²) = 50 + 24·(7.25−1) = 200.
        assert!(approx(t_origin.ixx, 200.0, 1e-9), "ixx {}", t_origin.ixx);
        assert!(approx(t_origin.iyy, 40.0 + 24.0 * (7.25 - 2.25), 1e-9));
        assert!(approx(t_origin.izz, 26.0 + 24.0 * (7.25 - 4.0), 1e-9));
        // Subtracting the parallel-axis term (negative mass) returns to the
        // centroid tensor.
        let back = t_origin.translated(-24.0, &GpVec::new(-1.0, -1.5, -2.0));
        assert!(approx(back.ixx, 50.0, 1e-9), "round-trip ixx {}", back.ixx);
    }

    #[test]
    fn inertia_composite_boxes() {
        // Two unit boxes side by side along x. Combined centroid at (1,0.5,0.5).
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = box_solid_at(&GpPnt::new(1.0, 0.0, 0.0), 1.0, 1.0, 1.0);
        let (t, c) = inertia_tensor_composite(&[&a.solid.0, &b.0], 1.0).unwrap();
        assert!(approx(c.x(), 1.0, 1e-9) && approx(c.y(), 0.5, 1e-9) && approx(c.z(), 0.5, 1e-9));
        // Each 1×1×1 box about its own centroid: Ixx=Iyy=Izz=1/6, m=1.
        // Shifted to the combined centroid: Ixx = 1/3, Iyy = Izz = 5/6.
        assert!(approx(t.ixx, 1.0 / 3.0, 1e-9), "ixx {}", t.ixx);
        assert!(approx(t.iyy, 5.0 / 6.0, 1e-9), "iyy {}", t.iyy);
        assert!(approx(t.izz, 5.0 / 6.0, 1e-9), "izz {}", t.izz);
        assert!(t.ixy.abs() < 1e-9);
    }

    #[test]
    fn curve_classify_line_circle() {
        let b = crate::builder::TopoBuilder::new();
        let line_edge = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(5.0, 0.0, 0.0));
        let curve = BRepTool::edge_curve(&line_edge).unwrap();
        let (a, bb) = BRepTool::edge_parameters(&line_edge);
        assert_eq!(classify_curve(curve.as_ref(), a, bb), CurveKind::Line);

        let circle_edge = b.make_edge_circle(&GpAx2::standard(), 2.0, 0.0, PI);
        let curve = BRepTool::edge_curve(&circle_edge).unwrap();
        let (a, bb) = BRepTool::edge_parameters(&circle_edge);
        assert_eq!(classify_curve(curve.as_ref(), a, bb), CurveKind::Circle);
        assert!(approx(edge_length_exact(&circle_edge), 2.0 * PI, 1e-9));
    }

    #[test]
    fn ellipse_ramanujan_circle_limit() {
        // A circle is an ellipse with a == b: the Ramanujan perimeter gives 2πa.
        assert!(approx(ellipse_length_ramanujan(2.0, 2.0), 4.0 * PI, 1e-9));
        // Degenerate flat ellipse a=b? Rather, an extreme a=10, b=1 gives a
        // value between 4a and 2π·(a+b)/... just check it is sane.
        let p = ellipse_length_ramanujan(10.0, 1.0);
        assert!(p > 40.0 && p < 44.0, "ellipse perimeter {p}");
    }

    #[test]
    fn classify_solid_analytics() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(shape_is_box(&b.solid.0));
        assert_eq!(classify_solid(&b.solid.0), SolidKind::Box);
        let s = BRepPrimSphere::make_sphere(1.0);
        assert_eq!(classify_solid(&s.solid.0), SolidKind::Sphere);
        let c = BRepPrimCylinder::make_cylinder(1.0, 1.0);
        assert_eq!(classify_solid(&c.solid.0), SolidKind::Cylinder);
        let k = BRepPrimCone::make_cone(1.0, 1.0);
        assert_eq!(classify_solid(&k.solid.0), SolidKind::Cone);
        let t = BRepPrimTorus::make_torus(1.0, 1.0);
        assert_eq!(classify_solid(&t.solid.0), SolidKind::Torus);
    }

    #[test]
    fn full_props_inertia_about() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let p = full_analytic_properties(&b.solid.0, 1.0).unwrap();
        // inertia_about(centroid) is the stored centroid tensor.
        let at_c = p.inertia_about(&p.centroid);
        assert!(approx(at_c.ixx, p.inertia.ixx, 1e-12));
        // Radius of gyration about the x principal axis: sqrt(50/24).
        let rg = p.radius_of_gyration();
        assert_eq!(rg.len(), 3);
        assert!(approx(rg[0], (50.0f64 / 24.0f64).sqrt(), 1e-9), "rg {}", rg[0]);
    }
}
