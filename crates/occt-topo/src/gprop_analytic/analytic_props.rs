use super::prelude::*;
use super::*;

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

pub(super) fn pnt_add_vec(p: GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::from_xyz(&p.coord.added(&v.coord))
}

pub(super) fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new((a.x() + b.x()) / 2.0, (a.y() + b.y()) / 2.0, (a.z() + b.z()) / 2.0)
}

/// Distance from `p` to the infinite line through `c` along unit `ax`.
pub(super) fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

/// A unit vector perpendicular to `v`.
pub(super) fn perpendicular(v: &GpVec) -> GpVec {
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
pub(super) fn classify_surface_full(s: &dyn Surface) -> SurfaceKind {
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
pub(super) fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
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
pub(super) fn torus_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64, f64)> {
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
pub(super) fn closest_point_lines(p1: GpPnt, d1: GpVec, p2: GpPnt, d2: GpVec) -> Option<GpPnt> {
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
pub(super) fn cone_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
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
pub(super) fn polygon_chain(face: &Face) -> Vec<GpPnt> {
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
pub(super) fn planar_face_area(face: &Face, pln: &GpPln) -> f64 {
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
pub(super) fn face_axial_extent(face: &Face, ax: &GpVec) -> f64 {
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
pub(super) fn shape_axial_extent(shape: &TopoShape, ax: &GpVec) -> f64 {
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
pub(super) fn cone_radius_height(shape: &TopoShape) -> Option<(f64, f64)> {
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
pub(super) fn polyhedron_volume(shape: &TopoShape) -> Option<f64> {
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

pub(super) fn sphere_volume_from_shape(shape: &TopoShape) -> Option<f64> {
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

pub(super) fn cylinder_volume_from_shape(shape: &TopoShape) -> Option<f64> {
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

pub(super) fn torus_volume_from_shape(shape: &TopoShape) -> Option<f64> {
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
