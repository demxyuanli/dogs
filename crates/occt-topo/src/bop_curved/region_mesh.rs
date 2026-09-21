use super::prelude::*;
use super::*;

/// Classification of a face relative to the other solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub enum FaceRegion {
    /// Every sampled face point is inside the other solid.
    Inside,
    /// Every sampled face point is outside the other solid.
    Outside,
    /// The face crosses the other solid's boundary.
    On,
}

/// A triangle mesh fragment (vertices + index triples).
#[derive(Debug, Clone, Default)]
pub struct TriMesh {
    pub verts: Vec<GpPnt>,
    pub tris: Vec<(usize, usize, usize)>,
}

// ---------------------------------------------------------------------------
// Point / ray classification
// ---------------------------------------------------------------------------

/// Möller–Trumbore ray/triangle intersection for a ray from `origin` in
/// direction `dir`. Returns the hit parameter `t` when the ray crosses the
/// triangle (excluding hits behind the origin).
pub(super) fn ray_triangle_hit(a: &GpPnt, b: &GpPnt, c: &GpPnt, origin: &GpPnt, dir: &GpVec) -> Option<f64> {
    pub(super) const EPS: f64 = 1e-12;
    let d = dir.xyz();
    let e1 = b.coord.subtracted(&a.coord);
    let e2 = c.coord.subtracted(&a.coord);
    let h = d.crossed(&e2);
    let det = e1.dot(&h);
    if det.abs() < EPS {
        return None;
    }
    let inv = 1.0 / det;
    let s = origin.coord.subtracted(&a.coord);
    let u = s.dot(&h) * inv;
    if u < -EPS || u > 1.0 + EPS {
        return None;
    }
    let q = s.crossed(&e1);
    let v = d.dot(&q) * inv;
    if v < -EPS || u + v > 1.0 + EPS {
        return None;
    }
    let t = e2.dot(&q) * inv;
    if t > EPS { Some(t) } else { None }
}

/// Number of crossings of the ray `p + t·dir` (unit `dir`) against the faces of
/// `solid`. Each face is grid-triangulated on its UV window; the count is odd
/// for points inside a closed solid.
pub fn ray_hits_solid(p: GpPnt, solid: &TopoShape, dir: GpVec, tol: f64) -> usize {
    let dir = if dir.magnitude() > 1e-30 { dir.normalized() } else { GpVec::new(1.0, 0.0, 0.0) };
    let _ = tol;
    let jitter = GpVec::new(1e-7, 1e-7, -2e-7);
    let origin = p.translated_vec(&jitter);
    let mut hits = 0usize;
    let meshes = crate::brep_extrema::mesh_faces(solid, 8, 8);
    for (verts, tris) in &meshes {
        for (i, j, k) in tris {
            if ray_triangle_hit(&verts[*i], &verts[*j], &verts[*k], &origin, &dir).is_some() {
                hits += 1;
            }
        }
    }
    hits
}

/// Even-odd point-in-solid test.
///
/// For solids made entirely of spheres the test is analytic (distance to each
/// sphere center), which is exact and avoids the coarse ray-cast mesh that
/// would misclassify points near the boundary. Other solids use
/// `brep_extrema::is_inside` (ray-cast against a face grid mesh).
pub fn point_in_solid(p: GpPnt, solid: &TopoShape, tol: f64) -> bool {
    let fa = faces_of(solid);
    if !fa.is_empty() {
        let spheres: Vec<(GpPnt, f64)> = fa.iter().filter_map(face_sphere).collect();
        if spheres.len() == fa.len() {
            return spheres.iter().any(|(c, r)| p.distance(c) <= r + tol);
        }
    }
    let _ = tol;
    crate::brep_extrema::is_inside(solid, &p)
}

/// Point-in-solid for the general curved boolean path: prefers analytic
/// sphere *and* cylinder tests (the coarse ray-cast mesh misclassifies points
/// near a cylindrical boundary) before falling back to the ray-cast mesh.
pub(super) fn point_in_solid_curved(p: GpPnt, solid: &TopoShape, tol: f64) -> bool {
    let fa = faces_of(solid);
    if !fa.is_empty() {
        let spheres: Vec<(GpPnt, f64)> = fa.iter().filter_map(face_sphere).collect();
        if spheres.len() == fa.len() {
            return spheres.iter().any(|(c, r)| p.distance(c) <= r + tol);
        }
        if let Some((c, dir, r, t0, t1)) = solid_cylinder(&fa) {
            let w = GpVec::from_pnts(&c, &p);
            let t = w.dot(&dir);
            let radial = w.subtracted(&dir.multiplied_scalar(t)).magnitude();
            return radial <= r + tol && t >= t0 - tol && t <= t1 + tol;
        }
    }
    crate::brep_extrema::is_inside(solid, &p)
}

/// Circumcenter of three points (used to locate a cylinder axis).
pub(super) fn circumcenter3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.crossed(&d2);
    if n.magnitude() < 1e-20 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    // O·d1 = (|b|²−|a|²)/2 ; O·d2 = (|c|²−|a|²)/2 ; O·n = a·n.
    let rhs1 = (n2(b) - n2(a)) * 0.5;
    let rhs2 = (n2(c) - n2(a)) * 0.5;
    let rhs3 = a.coord.dot(&n.coord);
    let m = GpMat3x3::new(d1, d2, n);
    m.solve(rhs1, rhs2, rhs3).map(|xyz| GpPnt::from_xyz(&xyz))
}

/// Small 3×3 matrix solve used by [`circumcenter3`].
pub(super) struct GpMat3x3 {
    pub(super) r: [[f64; 3]; 3],
}

impl GpMat3x3 {
    pub(super) fn new(c0: GpVec, c1: GpVec, c2: GpVec) -> Self {
        // Columns c0, c1, c2.
        Self {
            r: [
                [c0.x(), c1.x(), c2.x()],
                [c0.y(), c1.y(), c2.y()],
                [c0.z(), c1.z(), c2.z()],
            ],
        }
    }
    pub(super) fn det(&self) -> f64 {
        let r = &self.r;
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }
    /// Solve M·x = b by Cramer's rule.
    pub(super) fn solve(&self, b0: f64, b1: f64, b2: f64) -> Option<occt_core::gp::GpXyz> {
        let d = self.det();
        if d.abs() < 1e-20 {
            return None;
        }
        let r = &self.r;
        let det0 = b0 * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (b1 * r[2][2] - r[1][2] * b2)
            + r[0][2] * (b1 * r[2][1] - r[1][1] * b2);
        let det1 = r[0][0] * (b1 * r[2][2] - r[1][2] * b2)
            - b0 * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * b2 - b1 * r[2][0]);
        let det2 = r[0][0] * (r[1][1] * b2 - b1 * r[2][1])
            - r[0][1] * (r[1][0] * b2 - b1 * r[2][0])
            + b0 * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
        Some(occt_core::gp::GpXyz::new(det0 / d, det1 / d, det2 / d))
    }
}

/// Center, axis direction, radius and along-axis extent `(t0, t1)` of a solid
/// cylinder, detected as one cylindrical lateral face plus two planar caps.
pub(super) fn solid_cylinder(fa: &[Face]) -> Option<(GpPnt, GpVec, f64, f64, f64)> {
    let mut lateral: Option<&Face> = None;
    let mut caps: Vec<&Face> = Vec::new();
    for f in fa {
        let s = BRepTool::face_surface(f)?;
        match classify_surface(s.as_ref()) {
            SurfaceKind::Plane => caps.push(f),
            SurfaceKind::Sphere => return None,
            _ => {
                if lateral.is_some() {
                    return None; // more than one curved face → not a simple cylinder
                }
                lateral = Some(f);
            }
        }
    }
    let lat = lateral?;
    if caps.len() != 2 {
        return None;
    }
    let lat_surf = BRepTool::face_surface(lat)?;
    let (c, dir, r) = surface_cylinder_params(lat_surf.as_ref())?;
    let mut ts = Vec::new();
    for cap in &caps {
        let s = BRepTool::face_surface(cap)?;
        let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(s.as_ref());
        let (uc, vc) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let p = s.d0(uc, vc);
        ts.push(GpVec::from_pnts(&c, &p).dot(&dir));
    }
    let (t0, t1) = if ts[0] <= ts[1] { (ts[0], ts[1]) } else { (ts[1], ts[0]) };
    Some((c, dir, r, t0, t1))
}

/// Axis (point + direction) and radius of a cylindrical surface.
pub(super) fn surface_cylinder_params(s: &dyn occt_geom::Surface) -> Option<(GpPnt, GpVec, f64)> {
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(s);
    let v_mid = 0.5 * (v0 + v1);
    let (ua, ub, uc) = (u0, u0 + (u1 - u0) / 3.0, u0 + 2.0 * (u1 - u0) / 3.0);
    let p0 = s.d0(ua, v_mid);
    let p1 = s.d0(ub, v_mid);
    let p2 = s.d0(uc, v_mid);
    let c = circumcenter3(&p0, &p1, &p2)?;
    let dir = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2)).normalized();
    let r = p0.distance(&c);
    if r < 1e-9 {
        return None;
    }
    // Verify cylindrical: every sample at distance r from the axis.
    for i in 0..6usize {
        for j in 0..6usize {
            let u = u0 + (u1 - u0) * i as f64 / 5.0;
            let v = v0 + (v1 - v0) * j as f64 / 5.0;
            let p = s.d0(u, v);
            let w = GpVec::from_pnts(&c, &p);
            let radial = w.subtracted(&dir.multiplied_scalar(w.dot(&dir))).magnitude();
            if (radial - r).abs() > 1e-4 * r.max(1.0) {
                return None;
            }
        }
    }
    Some((c, dir, r))
}

/// Classify a face against the other solid by sampling its surface on a 5×5
/// grid and testing each sample with point-in-solid.
pub fn classify_face(f: &Face, solid: &TopoShape, tol: f64) -> FaceRegion {
    let surf = match BRepTool::face_surface(f) {
        Some(s) => s,
        None => return FaceRegion::Outside,
    };
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(surf.as_ref());
    let (nu, nv) = (5usize, 5usize);
    let mut inside = 0usize;
    let mut total = 0usize;
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let v = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            let p = surf.d0(u, v);
            if point_in_solid(p, solid, tol) {
                inside += 1;
            }
            total += 1;
        }
    }
    if inside == 0 {
        FaceRegion::Outside
    } else if inside == total {
        FaceRegion::Inside
    } else {
        FaceRegion::On
    }
}

/// Classify a general (non-planar, non-sphere) curved face against the other
/// solid by sampling its surface on a denser 8×8 grid. Falls back to the same
/// all-inside / all-outside / mixed rule as [`classify_face`] but with more
/// samples, which matters for strongly-curved `Other` faces whose sign can vary
/// within one coarse cell.
pub fn classify_face_general(f: &Face, solid: &TopoShape, tol: f64) -> FaceRegion {
    let surf = match BRepTool::face_surface(f) {
        Some(s) => s,
        None => return FaceRegion::Outside,
    };
    let (u0, u1, v0, v1) = face_uv_window_local(f);
    let (nu, nv) = (8usize, 8usize);
    let mut inside = 0usize;
    let mut total = 0usize;
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let v = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            let p = surf.d0(u, v);
            if point_in_solid_curved(p, solid, tol) {
                inside += 1;
            }
            total += 1;
        }
    }
    if inside == 0 {
        FaceRegion::Outside
    } else if inside == total {
        FaceRegion::Inside
    } else {
        FaceRegion::On
    }
}

/// Whether any face of `shape` carries a general (non-planar, non-sphere)
/// curved surface — i.e. the general boolean path must handle it.
pub(super) fn has_general_curved_face(shape: &TopoShape) -> bool {
    faces_of(shape).iter().any(|f| match BRepTool::face_surface(f) {
        Some(s) => {
            let k = classify_surface(s.as_ref());
            k != SurfaceKind::Plane && k != SurfaceKind::Sphere
        }
        None => false,
    })
}

// ---------------------------------------------------------------------------
// Mesh helpers
// ---------------------------------------------------------------------------

/// Signed divergence-theorem volume of a triangle mesh. For a closed,
/// consistently-oriented mesh this is the enclosed volume.
pub fn mesh_volume(m: &TriMesh) -> f64 {
    let mut vol = 0.0;
    for &(i, j, k) in &m.tris {
        let a = m.verts[i].coord;
        let b = m.verts[j].coord;
        let c = m.verts[k].coord;
        vol += a.dot_cross(&b, &c);
    }
    (vol / 6.0).abs()
}

/// Re-orient every triangle so its geometric normal points away from `center`
/// (used for spherical caps, where the grid winding near the poles is
/// ambiguous).
pub(super) fn orient_outward(m: &mut TriMesh, center: &GpPnt) {
    for t in &mut m.tris {
        let (i, j, k) = *t;
        let a = m.verts[i].coord;
        let b = m.verts[j].coord;
        let c = m.verts[k].coord;
        let centroid = a.added(&b).added(&c).multiplied(1.0 / 3.0);
        let n = b.subtracted(&a).crossed(&c.subtracted(&a));
        let out = centroid.subtracted(&center.coord);
        if n.dot(&out) < 0.0 {
            *t = (i, k, j);
        }
    }
}

/// Average 3D position of a solid's face-surface samples — the reference point
/// used to orient each face mesh outward (mirroring `orient_outward`'s use for
/// spherical caps).
pub(super) fn solid_centroid(shape: &TopoShape) -> GpPnt {
    let mut sum = occt_core::gp::GpXyz::zero();
    let mut n = 0usize;
    for f in faces_of(shape) {
        if let Some(s) = BRepTool::face_surface(&f) {
            let (u0, u1, v0, v1) = face_uv_window_local(&f);
            for i in 0..=4usize {
                for j in 0..=4usize {
                    let p = s.d0(u0 + (u1 - u0) * i as f64 / 4.0, v0 + (v1 - v0) * j as f64 / 4.0);
                    sum = sum.added(&p.coord);
                    n += 1;
                }
            }
        }
    }
    if n > 0 {
        GpPnt::from_xyz(&sum.multiplied(1.0 / n as f64))
    } else {
        GpPnt::zero()
    }
}

/// Grid-triangulate a face's surface on an `(nu+1)×(nv+1)` UV grid.
pub(super) fn surface_grid_mesh(surf: &dyn occt_geom::Surface, nu: usize, nv: usize) -> TriMesh {
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(surf);
    surface_grid_mesh_uv(surf, u0, u1, v0, v1, nu, nv)
}

/// Grid-triangulate a surface over an explicit `(u0, u1, v0, v1)` window.
pub(super) fn surface_grid_mesh_uv(surf: &dyn occt_geom::Surface, u0: f64, u1: f64, v0: f64, v1: f64, nu: usize, nv: usize) -> TriMesh {
    let (nu, nv) = (nu.max(2), nv.max(2));
    let mut verts = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            verts.push(surf.d0(u, v));
        }
    }
    let stride = nv + 1;
    let mut tris = Vec::with_capacity(2 * nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            // Wind so the triangle normal follows the surface normal (from d1).
            let (_, su, sv) = surf.d1(
                u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64,
                v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64,
            );
            let ns = su.coord.crossed(&sv.coord);
            if ns.square_modulus() > 1e-30 {
                let n1 = verts[b].coord.subtracted(&verts[a].coord).crossed(&verts[c].coord.subtracted(&verts[a].coord));
                if ns.dot(&n1) < 0.0 {
                    tris.push((a, c, b));
                    tris.push((c, d, b));
                } else {
                    tris.push((a, b, c));
                    tris.push((b, d, c));
                }
            } else {
                // Surfaces with zero d1 (spheres/cylinders) are wound by the
                // natural grid orientation, which is outward for the standard
                // parametric orientation.
                tris.push((a, b, c));
                tris.push((b, d, c));
            }
        }
    }
    TriMesh { verts, tris }
}

/// Mesh of a whole face (its full surface).
pub(super) fn face_full_mesh(face: &Face, nu: usize, nv: usize) -> TriMesh {
    match BRepTool::face_surface(face) {
        Some(s) => surface_grid_mesh(s.as_ref(), nu, nv),
        None => TriMesh::default(),
    }
}

/// Center + radius of a face whose surface is spherical.
///
/// `classify_surface` must confirm the surface is a sphere first: the
/// circumcenter heuristic alone returns a point for planar faces too (any three
/// coplanar points have a circumcenter), which would corrupt the point-in-solid
/// analytic fast path.
pub(super) fn face_sphere(f: &Face) -> Option<(GpPnt, f64)> {
    let s = BRepTool::face_surface(f)?;
    if classify_surface(s.as_ref()) != SurfaceKind::Sphere {
        return None;
    }
    crate::intpatch::sphere_params(s.as_ref())
}

/// The spheres present among a solid's faces.
pub(super) fn solid_spheres(solid: &TopoShape) -> Vec<(GpPnt, f64)> {
    let mut out = Vec::new();
    for f in faces_of(solid) {
        if let Some(sp) = face_sphere(&f) {
            out.push(sp);
        }
    }
    out
}

/// Intersection circle of two spheres, as `(center, radius, unit normal)`.
pub(super) fn sphere_sphere_circle(c1: GpPnt, r1: f64, c2: GpPnt, r2: f64) -> Option<(GpPnt, f64, GpVec)> {
    let d = c1.distance(&c2);
    if d <= 1e-12 {
        return None;
    }
    let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    if a > r1 || a < -r1 {
        return None;
    }
    let dir = GpVec::from_pnts(&c1, &c2).divided(d);
    let radius = (r1 * r1 - a * a).max(0.0).sqrt();
    let center = c1.translated_vec(&dir.multiplied_scalar(a));
    Some((center, radius, dir))
}

/// An in-plane orthonormal basis `(e1, e2)` perpendicular to `n`.
pub(super) fn perpendicular_basis(n: &GpVec) -> (GpVec, GpVec) {
    let z = GpVec::new(0.0, 0.0, 1.0);
    let e1 = if n.cross_magnitude(&z) > 1e-9 {
        n.crossed(&z).normalized()
    } else {
        GpVec::new(1.0, 0.0, 0.0)
    };
    let e2 = n.crossed(&e1).normalized();
    (e1, e2)
}

/// A discretized intersection circle shared by the two caps of an intersecting
/// pair, plus the in-plane basis used to build it (both caps must use the SAME
/// basis so their interior rows stay angularly aligned with the shared rim).
#[derive(Clone)]
pub struct IntersectionRim {
    /// The rim polygon vertices (in boundary order).
    pub points: Vec<GpPnt>,
    /// In-plane X basis of the rim circle.
    pub e1: GpVec,
    /// In-plane Y basis of the rim circle.
    pub e2: GpVec,
}

/// Discretize the intersection circle into `nu` points `(center + R·(cos·e1 +
/// sin·e2))`, shared by both spherical caps.
pub(super) fn rim_points(center: GpPnt, radius: f64, normal: &GpVec, nu: usize) -> IntersectionRim {
    let (e1, e2) = perpendicular_basis(normal);
    let points = (0..nu)
        .map(|i| {
            let a = 2.0 * PI * i as f64 / nu as f64;
            center
                .translated_vec(&e1.multiplied_scalar(radius * a.cos()))
                .translated_vec(&e2.multiplied_scalar(radius * a.sin()))
        })
        .collect();
    IntersectionRim { points, e1, e2 }
}

/// Mesh a spherical cap: the portion of the sphere `(c, r)` with polar angle
/// from `pole` in `[0, π/2 − v_lo]` (i.e. `v ∈ [v_lo, π/2]`). The rim row uses
/// the shared `rim` polygon (and its in-plane basis) so two caps of an
/// intersecting pair close along it with angularly aligned interior rows.
pub(super) fn sphere_cap_mesh(
    c: GpPnt,
    r: f64,
    pole: GpVec,
    v_lo: f64,
    rim: &IntersectionRim,
    nu: usize,
    nv: usize,
) -> TriMesh {
    let e1 = rim.e1;
    let e2 = rim.e2;
    let pole_u = pole.normalized();
    let (nu, nv) = (nu.max(4), nv.max(2));
    let v_lo = v_lo.clamp(-FRAC_PI_2, FRAC_PI_2);

    let mut verts: Vec<GpPnt> = Vec::new();
    // Rim row (row 0): the shared polygon.
    for p in rim.points.iter().take(nu) {
        verts.push(*p);
    }
    // Interior rows 1..=nv (row nv is the pole).
    for j in 1..=nv {
        let v = v_lo + (FRAC_PI_2 - v_lo) * j as f64 / nv as f64;
        let (cv, sv) = (v.cos(), v.sin());
        if j == nv {
            verts.push(c.translated_vec(&pole_u.multiplied_scalar(r)));
        } else {
            for i in 0..nu {
                let a = 2.0 * PI * i as f64 / nu as f64;
                let p = c
                    .translated_vec(&e1.multiplied_scalar(r * cv * a.cos()))
                    .translated_vec(&e2.multiplied_scalar(r * cv * a.sin()))
                    .translated_vec(&pole_u.multiplied_scalar(r * sv));
                verts.push(p);
            }
        }
    }
    let pole_idx = verts.len() - 1;
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();
    for j in 0..nv {
        let r0 = j * nu;
        let r1 = (j + 1) * nu;
        let r1_is_pole = j + 1 == nv;
        if r1_is_pole {
            for i in 0..nu {
                let a = r0 + i;
                let b = r0 + (i + 1) % nu;
                tris.push((a, b, pole_idx));
            }
        } else {
            for i in 0..nu {
                let a = r0 + i;
                let b = r0 + (i + 1) % nu;
                let cc = r1 + i;
                let dd = r1 + (i + 1) % nu;
                tris.push((a, b, cc));
                tris.push((b, dd, cc));
            }
        }
    }
    let mut mesh = TriMesh { verts, tris };
    orient_outward(&mut mesh, &c);
    mesh
}

/// Weld near-duplicate vertices (spatial hash) so shared rim edges close the
/// mesh. Returns remapped triangle indices.
pub(super) fn weld_mesh(m: &TriMesh, tol: f64) -> TriMesh {
    let tol = tol.max(1e-9);
    let cell = |p: &GpPnt| -> (i64, i64, i64) {
        (
            f64::floor(p.x() / tol) as i64,
            f64::floor(p.y() / tol) as i64,
            f64::floor(p.z() / tol) as i64,
        )
    };
    let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> = std::collections::HashMap::new();
    let mut unique: Vec<GpPnt> = Vec::new();
    let mut remap = vec![0usize; m.verts.len()];
    for (i, p) in m.verts.iter().enumerate() {
        let k = cell(p);
        let mut found = None;
        'search: for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                for dz in -1i64..=1 {
                    if let Some(bucket) = grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        for &j in bucket {
                            if p.distance(&unique[j]) <= tol {
                                found = Some(j);
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        match found {
            Some(j) => remap[i] = j,
            None => {
                let j = unique.len();
                unique.push(*p);
                grid.entry(k).or_default().push(j);
                remap[i] = j;
            }
        }
    }
    let tris = m
        .tris
        .iter()
        .map(|&(a, b, c)| (remap[a], remap[b], remap[c]))
        .filter(|&(a, b, c)| a != b && b != c && a != c)
        .collect();
    TriMesh { verts: unique, tris }
}

// ---------------------------------------------------------------------------
// Selection & assembly
// ---------------------------------------------------------------------------

/// Whether a face with the given classification is kept for the operation.
pub(super) fn select_keep(region: FaceRegion, from_a: bool, op: BoolOp) -> bool {
    match op {
        BoolOp::Fuse => region != FaceRegion::Inside,
        BoolOp::Cut => {
            if from_a {
                region != FaceRegion::Inside
            } else {
                region == FaceRegion::Inside
            }
        }
        BoolOp::Common => region != FaceRegion::Outside,
    }
}

/// Whether the retained side of a crossing face is toward the other solid.
pub(super) fn keep_inside(from_a: bool, op: BoolOp) -> bool {
    match op {
        BoolOp::Fuse => false,
        BoolOp::Cut => !from_a,
        BoolOp::Common => true,
    }
}

/// Are all faces of `shape` planar (eligible for the exact planar boolean)?
pub(super) fn all_faces_planar(shape: &TopoShape) -> bool {
    let fa = faces_of(shape);
    if fa.is_empty() {
        return false;
    }
    fa.iter().all(|f| match BRepTool::face_surface(f) {
        Some(s) => classify_surface(s.as_ref()) == SurfaceKind::Plane,
        None => false,
    })
}
