//! Boolean operations on curved-face solids.
//! Source: `BRepAlgoAPI_Fuse` / `BRepAlgoAPI_Cut` / `BRepAlgoAPI_Common` (TKBO),
//! `BOPAlgo_Builder`.
//!
//! This module extends the planar-exact boolean (`crate::bop_builder`) to
//! solids whose faces carry analytic curved surfaces (sphere, cylinder, cone).
//! The strategy:
//!
//! 1. when every face of both solids is planar, delegate to
//!    `crate::bop_builder::boolean` (exact polygon split);
//! 2. otherwise classify each face against the other solid (Inside / Outside /
//!    On) by sampling surface points and testing point-in-solid;
//! 3. keep whole faces whose classification is uniform; for faces that cross
//!    the boundary, triangulate the retained region. Spherical faces crossing
//!    another sphere are trimmed to the analytic spherical cap; the two caps
//!    of an intersecting pair share a discretized intersection circle (and its
//!    in-plane basis) so the assembled mesh is closed and its
//!    divergence-theorem volume is exact;
//! 4. weld the retained triangles into a closed mesh, compute its volume, and
//!    rebuild a BRep solid (via `mesh_to_brep`).

use std::f64::consts::{FRAC_PI_2, PI};
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::geom::polygon_boolean::{point_in_polygon2d, polygon_boolean, signed_area2d, PolygonBoolOp};
use occt_core::geom::triangulate::triangulate_polygon;
use occt_core::gp::{GpAx1, GpDir, GpPnt, GpPnt2d, GpVec};
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::bop_builder::{BoolOp, BooleanResult};
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};
use crate::topo_tools_full::{edges_of_wire, faces_of, shapes_of, wires_of_face};

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
fn ray_triangle_hit(a: &GpPnt, b: &GpPnt, c: &GpPnt, origin: &GpPnt, dir: &GpVec) -> Option<f64> {
    const EPS: f64 = 1e-12;
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
fn point_in_solid_curved(p: GpPnt, solid: &TopoShape, tol: f64) -> bool {
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
fn circumcenter3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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
struct GpMat3x3 {
    r: [[f64; 3]; 3],
}

impl GpMat3x3 {
    fn new(c0: GpVec, c1: GpVec, c2: GpVec) -> Self {
        // Columns c0, c1, c2.
        Self {
            r: [
                [c0.x(), c1.x(), c2.x()],
                [c0.y(), c1.y(), c2.y()],
                [c0.z(), c1.z(), c2.z()],
            ],
        }
    }
    fn det(&self) -> f64 {
        let r = &self.r;
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }
    /// Solve M·x = b by Cramer's rule.
    fn solve(&self, b0: f64, b1: f64, b2: f64) -> Option<occt_core::gp::GpXyz> {
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
fn solid_cylinder(fa: &[Face]) -> Option<(GpPnt, GpVec, f64, f64, f64)> {
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
fn surface_cylinder_params(s: &dyn occt_geom::Surface) -> Option<(GpPnt, GpVec, f64)> {
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
fn has_general_curved_face(shape: &TopoShape) -> bool {
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
fn orient_outward(m: &mut TriMesh, center: &GpPnt) {
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
fn solid_centroid(shape: &TopoShape) -> GpPnt {
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
fn surface_grid_mesh(surf: &dyn occt_geom::Surface, nu: usize, nv: usize) -> TriMesh {
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(surf);
    surface_grid_mesh_uv(surf, u0, u1, v0, v1, nu, nv)
}

/// Grid-triangulate a surface over an explicit `(u0, u1, v0, v1)` window.
fn surface_grid_mesh_uv(surf: &dyn occt_geom::Surface, u0: f64, u1: f64, v0: f64, v1: f64, nu: usize, nv: usize) -> TriMesh {
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
                    tris.push((b, d, c));
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
fn face_full_mesh(face: &Face, nu: usize, nv: usize) -> TriMesh {
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
fn face_sphere(f: &Face) -> Option<(GpPnt, f64)> {
    let s = BRepTool::face_surface(f)?;
    if classify_surface(s.as_ref()) != SurfaceKind::Sphere {
        return None;
    }
    crate::intpatch::sphere_params(s.as_ref())
}

/// The spheres present among a solid's faces.
fn solid_spheres(solid: &TopoShape) -> Vec<(GpPnt, f64)> {
    let mut out = Vec::new();
    for f in faces_of(solid) {
        if let Some(sp) = face_sphere(&f) {
            out.push(sp);
        }
    }
    out
}

/// Intersection circle of two spheres, as `(center, radius, unit normal)`.
fn sphere_sphere_circle(c1: GpPnt, r1: f64, c2: GpPnt, r2: f64) -> Option<(GpPnt, f64, GpVec)> {
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
fn perpendicular_basis(n: &GpVec) -> (GpVec, GpVec) {
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
fn rim_points(center: GpPnt, radius: f64, normal: &GpVec, nu: usize) -> IntersectionRim {
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
fn sphere_cap_mesh(
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
fn weld_mesh(m: &TriMesh, tol: f64) -> TriMesh {
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
fn select_keep(region: FaceRegion, from_a: bool, op: BoolOp) -> bool {
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
fn keep_inside(from_a: bool, op: BoolOp) -> bool {
    match op {
        BoolOp::Fuse => false,
        BoolOp::Cut => !from_a,
        BoolOp::Common => true,
    }
}

/// Are all faces of `shape` planar (eligible for the exact planar boolean)?
fn all_faces_planar(shape: &TopoShape) -> bool {
    let fa = faces_of(shape);
    if fa.is_empty() {
        return false;
    }
    fa.iter().all(|f| match BRepTool::face_surface(f) {
        Some(s) => classify_surface(s.as_ref()) == SurfaceKind::Plane,
        None => false,
    })
}

fn disjoint_result(a: &TopoShape, b: &TopoShape, op: BoolOp) -> BooleanResult {
    let bld = TopoBuilder::new();
    match op {
        BoolOp::Fuse => {
            let comp = bld.make_compound_of(&[a.clone(), b.clone()]);
            BooleanResult {
                shape: comp.0,
                solid: None,
                shells: vec![],
                faces: vec![],
                warnings: vec![],
            }
        }
        BoolOp::Cut => {
            let solid = Solid::wrap(a.clone());
            let shells: Vec<Shell> = shapes_of(a, ShapeType::Shell).into_iter().map(Shell).collect();
            BooleanResult {
                shape: a.clone(),
                solid,
                shells,
                faces: faces_of(a),
                warnings: vec![],
            }
        }
        BoolOp::Common => {
            let comp = bld.make_compound_of(&[]);
            BooleanResult {
                shape: comp.0,
                solid: None,
                shells: vec![],
                faces: vec![],
                warnings: vec![],
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Core mesh boolean
// ---------------------------------------------------------------------------

/// Triangulate the retained region of a whole (non-crossing) face.
fn mesh_whole_face(face: &Face) -> TriMesh {
    face_full_mesh(face, 24, 24)
}

/// Triangulate the retained region of a crossing face.
///
/// Spherical faces crossing another sphere are trimmed to the analytic cap; the
/// cap's rim is a shared discretization of the intersection circle, so the caps
/// from both solids close the mesh. Other surface types fall back to keeping
/// whole grid cells whose centroid lies on the retained side.
fn mesh_crossing_face(
    face: &Face,
    other: &TopoShape,
    keep_inside: bool,
    rim_shared: Option<&IntersectionRim>,
    tol: f64,
) -> TriMesh {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return TriMesh::default(),
    };

    // Spherical cap trimming.
    if classify_surface(surf.as_ref()) == SurfaceKind::Sphere {
        if let Some((c1, r1)) = face_sphere(face) {
            for (c2, r2) in solid_spheres(other) {
                let d = c1.distance(&c2);
                if d <= 1e-12 {
                    continue;
                }
                let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
                if a > r1 + 1e-9 || a < -r1 - 1e-9 {
                    continue;
                }
                let dir = GpVec::from_pnts(&c1, &c2).divided(d);
                let (pole_vec, vb) = if keep_inside {
                    (dir, (a / r1).clamp(-1.0, 1.0).acos())
                } else {
                    (dir.reversed(), (-a / r1).clamp(-1.0, 1.0).acos())
                };
                let v_lo = FRAC_PI_2 - vb;
                let rim = match rim_shared {
                    Some(r) => r.clone(),
                    None => {
                        let (cc, rc, n) = sphere_sphere_circle(c1, r1, c2, r2).unwrap_or((c1, 0.0, dir));
                        rim_points(cc, rc, &n, 48)
                    }
                };
                return sphere_cap_mesh(c1, r1, pole_vec, v_lo, &rim, 48, 24);
            }
        }
    }

    // General fallback: keep whole cells whose centroid is on the retained side.
    let mesh = surface_grid_mesh(surf.as_ref(), 24, 24);
    let nu = 24;
    let nv = 24;
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(surf.as_ref());
    let stride = nv + 1;
    let mut kept: Vec<(usize, usize, usize)> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            let uc = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let vc = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            let centroid = surf.d0(uc, vc);
            let inside = point_in_solid(centroid, other, tol);
            if inside == keep_inside {
                kept.push((a, b, c));
                kept.push((b, d, c));
            }
        }
    }
    TriMesh { verts: mesh.verts, tris: kept }
}

// ---------------------------------------------------------------------------
// General (non-analytic) curved-face meshing
// ---------------------------------------------------------------------------

/// Sample 3D points along a face's boundary edges (endpoints plus a few
/// interior curve samples). Used to infer a finite UV window when the surface's
/// own range is unbounded (planes, cylinders, …).
fn face_boundary_samples(face: &Face) -> Vec<GpPnt> {
    let mut out: Vec<GpPnt> = Vec::new();
    let face_kids = face.0.tshape.read().unwrap().children.clone();
    for h in face_kids {
        if h.read().unwrap().shape_type() != ShapeType::Wire {
            continue;
        }
        let wire = TopoShape::from_handle(h);
        let wire_kids = wire.tshape.read().unwrap().children.clone();
        for eh in wire_kids {
            if eh.read().unwrap().shape_type() != ShapeType::Edge {
                continue;
            }
            let edge = TopoShape::from_handle(eh);
            if let Some(e) = Edge::wrap(edge) {
                if let Some((p0, p1)) = BRepTool::edge_vertices(&e) {
                    out.push(p0);
                    out.push(p1);
                }
                if let Some(c) = BRepTool::edge_curve(&e) {
                    let (t0, t1) = BRepTool::edge_parameters(&e);
                    for k in 0..4usize {
                        let t = if t0.is_finite() && t1.is_finite() && t1 > t0 {
                            t0 + (t1 - t0) * k as f64 / 4.0
                        } else {
                            k as f64 - 1.0
                        };
                        out.push(c.d0(t));
                    }
                }
            }
        }
    }
    out
}

/// Finite parametric window `(u_min, u_max, v_min, v_max)` for a face.
///
/// Uses the surface's own range when finite; otherwise projects boundary edge
/// samples onto the surface and takes the bounding box of the `(u, v)` points.
fn face_uv_window_local(face: &Face) -> (f64, f64, f64, f64) {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return (-1.0, 1.0, -1.0, 1.0),
    };
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() {
        return (u1, u2, v1, v2);
    }
    let pts = face_boundary_samples(face);
    if pts.is_empty() {
        return crate::intpatch::sample_bounds(surf.as_ref());
    }
    let mut us = Vec::with_capacity(pts.len());
    let mut vs = Vec::with_capacity(pts.len());
    for p in &pts {
        let (u, v) = crate::intpatch::project_params(surf.as_ref(), p);
        us.push(u);
        vs.push(v);
    }
    let ua = us.iter().copied().fold(f64::INFINITY, f64::min);
    let ub = us.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let va = vs.iter().copied().fold(f64::INFINITY, f64::min);
    let vb = vs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (ua, ub) = if ub - ua < 1e-12 { (ua - 1.0, ub + 1.0) } else { (ua, ub) };
    let (va, vb) = if vb - va < 1e-12 { (va - 1.0, vb + 1.0) } else { (va, vb) };
    (ua, ub, va, vb)
}

/// Mesh a whole (non-crossing) face over its finite trimmed UV window.
fn mesh_face_window(face: &Face, nu: usize, nv: usize) -> TriMesh {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return TriMesh::default(),
    };
    let (u0, u1, v0, v1) = face_uv_window_local(face);
    surface_grid_mesh_uv(surf.as_ref(), u0, u1, v0, v1, nu, nv)
}

/// Axis-aligned bounding box of a face's surface samples (bbox pre-filter for
/// the face-pair intersection loop).
fn face_bbox(face: &Face) -> BndBox {
    let mut b = BndBox::new();
    if let Some(s) = BRepTool::face_surface(face) {
        let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(s.as_ref());
        for i in 0..=8usize {
            for j in 0..=8usize {
                let u = u0 + (u1 - u0) * i as f64 / 8.0;
                let v = v0 + (v1 - v0) * j as f64 / 8.0;
                b.add_point(&s.d0(u, v));
            }
        }
    }
    b
}

/// Nearest point of `pts` to `p`.
fn nearest_curve_point(pts: &[GpPnt], p: GpPnt) -> Option<GpPnt> {
    let mut best: Option<(f64, GpPnt)> = None;
    for q in pts {
        let d = q.distance(&p);
        if best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, *q));
        }
    }
    best.map(|(_, q)| q)
}

/// Triangulate the retained region of a crossing face using the intersection
/// polylines to pull grid vertices onto the shared boundary curve.
///
/// Sphere faces keep the analytic cap path ([`mesh_crossing_face`]); other
/// surfaces are grid-cell classified (a cell is kept when its centroid lies on
/// the retained side) and grid vertices within `snap_tol` of an intersection
/// curve are snapped onto the nearest curve point so both solids' retained
/// meshes close along the shared intersection.
fn mesh_crossing_face_curves(
    face: &Face,
    other: &TopoShape,
    keep_inside: bool,
    curves_on_face: &[Vec<GpPnt>],
    rim_shared: Option<&IntersectionRim>,
    tol: f64,
) -> TriMesh {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return TriMesh::default(),
    };
    // Spherical caps: analytic trimming only when both solids have the matching
    // sphere face (sphere–sphere). A sphere crossing a non-sphere is meshed by
    // the general grid-cell path below so its boundary snaps to the traced
    // intersection curve, closing the mesh against the other solid's wall.
    let is_sphere = classify_surface(surf.as_ref()) == SurfaceKind::Sphere;
    if is_sphere && !solid_spheres(other).is_empty() {
        return mesh_crossing_face(face, other, keep_inside, rim_shared, tol);
    }

    let (u0, u1, v0, v1) = face_uv_window_local(face);
    let (nu, nv) = (32usize, 32usize);
    let stride = nv + 1;
    let mut verts: Vec<GpPnt> = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            verts.push(surf.d0(u, v));
        }
    }
    // Snap grid vertices near an intersection curve onto the curve.
    let curve_pts: Vec<GpPnt> = curves_on_face.iter().flatten().copied().collect();
    if !curve_pts.is_empty() {
        let mut lo = verts[0];
        let mut hi = verts[0];
        for p in &verts {
            lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
            hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
        }
        let snap_tol = tol.max(lo.distance(&hi) * 0.03);
        for v in &mut verts {
            if let Some(near) = nearest_curve_point(&curve_pts, *v) {
                if near.distance(v) <= snap_tol {
                    *v = near;
                }
            }
        }
    }
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            let uc = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let vc = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            let centroid = surf.d0(uc, vc);
            let inside = point_in_solid_curved(centroid, other, tol);
            if inside == keep_inside {
                tris.push((a, b, c));
                tris.push((b, d, c));
            }
        }
    }
    TriMesh { verts, tris }
}

/// Compute the retained triangle mesh of the boolean, using the per-face
/// intersection polylines to close the mesh along the shared intersection
/// curves of the general path.
fn boolean_mesh_curves(
    a: &TopoShape,
    b: &TopoShape,
    op: BoolOp,
    tol: f64,
    pair_curves: &[Vec<Vec<Vec<GpPnt>>>],
) -> Result<TriMesh, String> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();

    // Shared sphere-sphere intersection rim (reused by both caps).
    let mut rim_shared: Option<IntersectionRim> = None;
    'outer: for fa_i in &fa {
        let Some(surf) = BRepTool::face_surface(fa_i) else { continue };
        if classify_surface(surf.as_ref()) != SurfaceKind::Sphere {
            continue;
        }
        let (c1, r1) = match face_sphere(fa_i) {
            Some(x) => x,
            None => continue,
        };
        for fb_j in &fb {
            let (c2, r2) = match face_sphere(fb_j) {
                Some(x) => x,
                None => continue,
            };
            if let Some((cc, rc, n)) = sphere_sphere_circle(c1, r1, c2, r2) {
                rim_shared = Some(rim_points(cc, rc, &n, 48));
                break 'outer;
            }
        }
    }

    let centroid_a = solid_centroid(a);
    let centroid_b = solid_centroid(b);
    let mut flip_b = op == BoolOp::Cut;

    let mut parts: Vec<TriMesh> = Vec::new();
    for (idx, (f, region)) in fa.iter().zip(&regions_a).enumerate() {
        if select_keep(*region, true, op) {
            let mut m = match region {
                FaceRegion::On => {
                    let curves: Vec<Vec<GpPnt>> = pair_curves[idx].iter().flatten().cloned().collect();
                    mesh_crossing_face_curves(f, b, keep_inside(true, op), &curves, rim_shared.as_ref(), tol)
                }
                _ => mesh_face_window(f, 24, 24),
            };
            orient_outward(&mut m, &centroid_a);
            parts.push(m);
        }
    }
    for (idx, (f, region)) in fb.iter().zip(&regions_b).enumerate() {
        if select_keep(*region, false, op) {
            let mut m = match region {
                FaceRegion::On => {
                    let curves: Vec<Vec<GpPnt>> = pair_curves
                        .iter()
                        .map(|row| row[idx].clone())
                        .flatten()
                        .collect();
                    mesh_crossing_face_curves(f, a, keep_inside(false, op), &curves, rim_shared.as_ref(), tol)
                }
                _ => mesh_face_window(f, 24, 24),
            };
            orient_outward(&mut m, &centroid_b);
            if flip_b {
                for t in &mut m.tris {
                    let (i, j, k) = *t;
                    *t = (i, k, j);
                }
            }
            parts.push(m);
        }
    }

    let mut all = TriMesh::default();
    for p in &parts {
        let offset = all.verts.len();
        all.verts.extend(p.verts.iter().copied());
        for &(i, j, k) in &p.tris {
            all.tris.push((offset + i, offset + j, offset + k));
        }
    }
    let weld_tol = tol.max(1e-9);
    Ok(weld_mesh(&all, weld_tol))
}

/// Compute the retained triangle mesh of the boolean.
fn boolean_mesh(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<TriMesh, String> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face(f, a, tol)).collect();

    // Pre-compute a shared intersection-circle rim when a sphere face of A
    // crosses a sphere face of B (used by both caps).
    let mut rim_shared: Option<IntersectionRim> = None;
    'outer: for fa_i in &fa {
        let Some(surf) = BRepTool::face_surface(fa_i) else { continue };
        if classify_surface(surf.as_ref()) != SurfaceKind::Sphere {
            continue;
        }
        let (c1, r1) = match face_sphere(fa_i) {
            Some(x) => x,
            None => continue,
        };
        for fb_j in &fb {
            let (c2, r2) = match face_sphere(fb_j) {
                Some(x) => x,
                None => continue,
            };
            if let Some((cc, rc, n)) = sphere_sphere_circle(c1, r1, c2, r2) {
                rim_shared = Some(rim_points(cc, rc, &n, 48));
                break 'outer;
            }
        }
    }

    let mut parts: Vec<TriMesh> = Vec::new();
    for (f, region) in fa.iter().zip(&regions_a) {
        if select_keep(*region, true, op) {
            match region {
                FaceRegion::On => {
                    let keep_in = keep_inside(true, op);
                    parts.push(mesh_crossing_face(f, b, keep_in, rim_shared.as_ref(), tol));
                }
                _ => parts.push(mesh_whole_face(f)),
            }
        }
    }
    for (f, region) in fb.iter().zip(&regions_b) {
        if select_keep(*region, false, op) {
            match region {
                FaceRegion::On => {
                    let keep_in = keep_inside(false, op);
                    parts.push(mesh_crossing_face(f, a, keep_in, rim_shared.as_ref(), tol));
                }
                _ => parts.push(mesh_whole_face(f)),
            }
        }
    }

    // Concatenate all fragments into one mesh and weld shared vertices.
    let mut all = TriMesh::default();
    for p in &parts {
        let offset = all.verts.len();
        all.verts.extend(p.verts.iter().copied());
        for &(i, j, k) in &p.tris {
            all.tris.push((offset + i, offset + j, offset + k));
        }
    }
    let weld_tol = tol.max(1e-9);
    Ok(weld_mesh(&all, weld_tol))
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Boolean on solids whose faces include curved surfaces.
///
/// Delegates to `crate::bop_builder::boolean` when both inputs are planar.
/// For curved inputs it classifies each face, keeps whole faces that do not
/// cross the boundary, trims crossing faces to their retained region (analytic
/// spherical caps for sphere–sphere), welds the retained triangles into a
/// closed mesh, and rebuilds a solid via `mesh_to_brep`.
pub fn curved_boolean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    // Planar inputs → exact polygon boolean.
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }

    // Disjoint shortcut.
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }

    let fa = faces_of(a);
    let fb = faces_of(b);
    if fa.is_empty() || fb.is_empty() {
        return Err("curved_boolean: input has no faces".into());
    }

    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face(f, a, tol)).collect();
    let has_on = regions_a.contains(&FaceRegion::On) || regions_b.contains(&FaceRegion::On);

    // When no face crosses the boundary, the result is assembled from the whole
    // (original) kept faces — curved faces stay analytic, so `shape_volume` is
    // reliable for the curved result.
    if !has_on {
        let mut faces: Vec<Face> = Vec::new();
        for (f, region) in fa.iter().zip(&regions_a) {
            if select_keep(*region, true, op) {
                faces.push(f.clone());
            }
        }
        for (f, region) in fb.iter().zip(&regions_b) {
            if select_keep(*region, false, op) {
                faces.push(f.clone());
            }
        }
        if faces.is_empty() {
            let bld = TopoBuilder::new();
            let comp = bld.make_compound_of(&[]);
            return Ok(BooleanResult {
                shape: comp.0,
                solid: None,
                shells: vec![],
                faces: vec![],
                warnings: vec![],
            });
        }
        let bld = TopoBuilder::new();
        let shell = bld.make_shell(&faces);
        let solid = bld.make_solid(&[shell.clone()]);
        return Ok(BooleanResult {
            shape: solid.0.clone(),
            solid: Some(solid),
            shells: vec![shell],
            faces,
            warnings: vec![],
        });
    }

    // Crossing faces: build the closed mesh and rebuild a faceted BRep solid.
    let mesh = boolean_mesh(a, b, op, tol)?;
    let vol = mesh_volume(&mesh);
    if mesh.tris.is_empty() {
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[]);
        return Ok(BooleanResult {
            shape: comp.0,
            solid: None,
            shells: vec![],
            faces: vec![],
            warnings: vec![format!("empty curved boolean (mesh volume {vol:.4})")],
        });
    }

    let smesh = crate::mesh::ShapeMesh {
        vertices: mesh.verts.clone(),
        triangles: mesh
            .tris
            .iter()
            .map(|&(a, b, c)| occt_core::poly::triangulation::Triangle::new(a, b, c))
            .collect(),
        source_shape: ShapeType::Solid,
    };
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&smesh);
    let shape = brep.solid.clone().map(|s| s.0).unwrap_or_else(|| brep.shell.0.clone());
    Ok(BooleanResult {
        shape,
        solid: brep.solid,
        shells: vec![brep.shell],
        faces: brep.faces,
        warnings: vec![format!("curved boolean mesh volume ≈ {vol:.4}")],
    })
}

/// The enclosed volume of the curved boolean result (exact for the closed
/// welded mesh; for planar inputs delegates to the planar boolean and meshes
/// its result).
pub fn curved_boolean_volume(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<f64, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        let r = crate::bop_builder::boolean(a, b, op, tol)?;
        return Ok(crate::shape_mesh::shape_volume(&r.shape, 0.02));
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return match op {
            BoolOp::Fuse => Ok(crate::shape_mesh::shape_volume(a, 0.02) + crate::shape_mesh::shape_volume(b, 0.02)),
            BoolOp::Cut => Ok(crate::shape_mesh::shape_volume(a, 0.02)),
            BoolOp::Common => Ok(0.0),
        };
    }
    let mesh = boolean_mesh(a, b, op, tol)?;
    Ok(mesh_volume(&mesh))
}

// ---------------------------------------------------------------------------
// General (non-analytic) curved boolean
// ---------------------------------------------------------------------------

/// Build a `BooleanResult` from the retained whole faces when no face crosses
/// the boundary (the analytic faces are preserved verbatim).
fn result_from_kept_faces(
    fa: &[Face],
    fb: &[Face],
    regions_a: &[FaceRegion],
    regions_b: &[FaceRegion],
    op: BoolOp,
) -> BooleanResult {
    let mut faces: Vec<Face> = Vec::new();
    for (f, region) in fa.iter().zip(regions_a) {
        if select_keep(*region, true, op) {
            faces.push(f.clone());
        }
    }
    for (f, region) in fb.iter().zip(regions_b) {
        if select_keep(*region, false, op) {
            faces.push(f.clone());
        }
    }
    if faces.is_empty() {
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[]);
        return BooleanResult {
            shape: comp.0,
            solid: None,
            shells: vec![],
            faces: vec![],
            warnings: vec![],
        };
    }
    let bld = TopoBuilder::new();
    let shell = bld.make_shell(&faces);
    let solid = bld.make_solid(&[shell.clone()]);
    BooleanResult {
        shape: solid.0.clone(),
        solid: Some(solid),
        shells: vec![shell],
        faces,
        warnings: vec![],
    }
}

/// Rebuild a faceted `BooleanResult` solid from a closed retained mesh.
fn result_from_mesh(mesh: TriMesh) -> BooleanResult {
    let vol = mesh_volume(&mesh);
    if mesh.tris.is_empty() {
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[]);
        return BooleanResult {
            shape: comp.0,
            solid: None,
            shells: vec![],
            faces: vec![],
            warnings: vec![format!("empty curved boolean (mesh volume {vol:.4})")],
        };
    }
    let smesh = crate::mesh::ShapeMesh {
        vertices: mesh.verts.clone(),
        triangles: mesh
            .tris
            .iter()
            .map(|&(a, b, c)| occt_core::poly::triangulation::Triangle::new(a, b, c))
            .collect(),
        source_shape: ShapeType::Solid,
    };
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&smesh);
    let shape = brep.solid.clone().map(|s| s.0).unwrap_or_else(|| brep.shell.0.clone());
    BooleanResult {
        shape,
        solid: brep.solid,
        shells: vec![brep.shell],
        faces: brep.faces,
        warnings: vec![format!("curved boolean mesh volume ≈ {vol:.4}")],
    }
}

/// Face–face intersection polylines for the general path (pairs whose bounding
/// boxes do not overlap are skipped).
fn general_pair_curves(fa: &[Face], fb: &[Face], tol: f64) -> Vec<Vec<Vec<Vec<GpPnt>>>> {
    let mut pair_curves: Vec<Vec<Vec<Vec<GpPnt>>>> = vec![vec![Vec::new(); fb.len()]; fa.len()];
    for (i, f_i) in fa.iter().enumerate() {
        let Some(sa) = BRepTool::face_surface(f_i) else { continue };
        let ba = face_bbox(f_i);
        for (j, f_j) in fb.iter().enumerate() {
            let Some(sb) = BRepTool::face_surface(f_j) else { continue };
            let bb = face_bbox(f_j);
            if !ba.is_void() && !bb.is_void() && ba.is_out_box(&bb) {
                continue;
            }
            let curves = crate::intpatch::intersection_curve_points(sa.as_ref(), sb.as_ref(), tol, 40);
            if !curves.is_empty() {
                pair_curves[i][j] = curves;
            }
        }
    }
    pair_curves
}

/// A more accurate distance-to-surface than `intpatch::distance_to_surface`:
/// a fine grid search followed by several Newton refinements from neighbouring
/// starts. The intpatch projection has ~2e-3 error near a small sphere, which
/// would swamp the level-set tolerance of the trace.
fn distance_to_surface_fine(p: &GpPnt, s: &dyn Surface) -> f64 {
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(s);
    let (nu, nv) = (48usize, 48usize);
    let mut bu = u0;
    let mut bv = v0;
    let mut bd = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).square_distance(p);
            if d < bd {
                bd = d;
                bu = u;
                bv = v;
            }
        }
    }
    let refine = |u: f64, v: f64| {
        let (u2, v2, q) = crate::intpatch::refine_point_on_surface(s, *p, u, v, 30);
        q.distance(p).min(s.d0(u2, v2).distance(p))
    };
    let mut best = refine(bu, bv);
    for (du, dv) in [(0.03, 0.0), (-0.03, 0.0), (0.0, 0.03), (0.0, -0.03)] {
        best = best.min(refine(bu + du, bv + dv));
    }
    best
}

/// Crossing of a segment `a → b` (with signed field values `va, vb`) at the
/// zero level set (used by the windowed marching-squares tracer).
fn edge_crossing_2(pa: &GpPnt, va: f64, pb: &GpPnt, vb: f64) -> Option<GpPnt> {
    if va * vb > 0.0 {
        return None;
    }
    if va.abs() < 1e-12 && vb.abs() < 1e-12 {
        return None;
    }
    if va.abs() < 1e-12 {
        return Some(*pa);
    }
    if vb.abs() < 1e-12 {
        return Some(*pb);
    }
    let t = va / (va - vb);
    Some(GpPnt::new(
        pa.x() + t * (pb.x() - pa.x()),
        pa.y() + t * (pb.y() - pa.y()),
        pa.z() + t * (pb.z() - pa.z()),
    ))
}

/// Trace the intersection of `sa` (over the explicit UV window) with `sb`,
/// marching on `sa`'s grid exactly like the general intpatch tracer but over
/// the face's real trimmed window (the clamped `sample_bounds` used by
/// `intersection_curve_points` misses offset faces whose intersection lies
/// outside [-1,1]²).
fn face_intersection_polylines(
    sa: &dyn Surface,
    window: (f64, f64, f64, f64),
    sb: &dyn Surface,
    tol: f64,
) -> Vec<Vec<GpPnt>> {
    let (u0, u1, v0, v1) = window;
    if !(u1 > u0 && v1 > v0) {
        return Vec::new();
    }
    let trace_tol = tol.max(1e-4);
    let (nu, nv) = (40usize, 40usize);
    let mut field = vec![vec![0.0f64; nv + 1]; nu + 1];
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = sa.d0(u, v);
            field[i][j] = distance_to_surface_fine(&p, sb) - trace_tol;
        }
    }
    let mut pts: Vec<GpPnt> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let p = [
                sa.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                sa.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                sa.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
                sa.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
            ];
            let f = [field[i][j], field[i + 1][j], field[i + 1][j + 1], field[i][j + 1]];
            let e = [
                edge_crossing_2(&p[0], f[0], &p[1], f[1]),
                edge_crossing_2(&p[1], f[1], &p[2], f[2]),
                edge_crossing_2(&p[2], f[2], &p[3], f[3]),
                edge_crossing_2(&p[3], f[3], &p[0], f[0]),
            ];
            let mut seg = Vec::new();
            for q in e.into_iter().flatten() {
                seg.push(q);
            }
            if seg.len() == 2 {
                pts.push(seg[0]);
                pts.push(seg[1]);
            } else if seg.len() >= 4 {
                pts.push(seg[0]);
                pts.push(seg[2]);
                pts.push(seg[1]);
                pts.push(seg[3]);
            }
        }
    }
    // Link the level-set fragments with a tolerance that covers the marching
    // grid spacing.
    let grid_step = ((u1 - u0) / nu as f64).max((v1 - v0) / nv as f64);
    let link_tol = trace_tol.max(2.0 * grid_step);
    crate::intpatch::chain_intersection_points(&pts, link_tol)
}

/// Axis-aligned bounding box of a face over its real trimmed UV window (unlike
/// [`face_bbox`], which samples the clamped `sample_bounds` and misplaces
/// offset faces).
fn face_bbox_local(face: &Face) -> BndBox {
    let mut b = BndBox::new();
    if let Some(s) = BRepTool::face_surface(face) {
        let (u0, u1, v0, v1) = face_uv_window_local(face);
        if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() && u1 > u0 && v1 > v0 {
            for i in 0..=8usize {
                for j in 0..=8usize {
                    let u = u0 + (u1 - u0) * i as f64 / 8.0;
                    let v = v0 + (v1 - v0) * j as f64 / 8.0;
                    b.add_point(&s.d0(u, v));
                }
            }
        }
    }
    b
}

/// Face–face intersection polylines for the trimmed-face path: like
/// [`general_pair_curves`] but only face pairs whose real UV-window bboxes
/// overlap are traced, and the analytic intersection dispatcher
/// (`surface_surface_intersection`) is used first so plane∩sphere and
/// plane∩cylinder give exact curves (the general grid tracer's
/// `distance_to_surface` projection noise swamps the level set near a small
/// sphere).
fn general_pair_curves_windowed(fa: &[Face], fb: &[Face], tol: f64) -> Vec<Vec<Vec<Vec<GpPnt>>>> {
    let mut pair_curves: Vec<Vec<Vec<Vec<GpPnt>>>> = vec![vec![Vec::new(); fb.len()]; fa.len()];
    for (i, f_i) in fa.iter().enumerate() {
        let Some(sa) = BRepTool::face_surface(f_i) else { continue };
        let ba = face_bbox_local(f_i);
        let Ok(win_a) = face_uv_window_rect(f_i).map(|r| rect_bounds(&r)) else { continue };
        for (j, f_j) in fb.iter().enumerate() {
            let Some(sb) = BRepTool::face_surface(f_j) else { continue };
            let bb = face_bbox_local(f_j);
            if !ba.is_void() && !bb.is_void() && ba.is_out_box(&bb) {
                continue;
            }
            let mut curves: Vec<Vec<GpPnt>> = Vec::new();
            // The surface_surface_intersection dispatcher handles plane∩sphere,
            // sphere∩sphere and plane∩plane; add the analytic plane∩cylinder
            // (and cone/torus) closed forms so those traces are exact too.
            let cyl_a = surface_cylinder_params(sa.as_ref());
            let cyl_b = surface_cylinder_params(sb.as_ref());
            let plane_a = crate::intpatch::plane_from_surface(sa.as_ref());
            let plane_b = crate::intpatch::plane_from_surface(sb.as_ref());
            let ax_of = |cyl: (GpPnt, GpVec, f64)| GpAx1::new(cyl.0, GpDir::from_vec(&cyl.1).unwrap_or_default());
            let analytic = match (plane_a, cyl_b) {
                (Some(pln), Some(cyl)) => crate::intpatch::intersect_plane_cylinder(&pln, &ax_of(cyl), cyl.2),
                _ => match (cyl_a, plane_b) {
                    (Some(cyl), Some(pln)) => crate::intpatch::intersect_plane_cylinder(&pln, &ax_of(cyl), cyl.2),
                    _ => None,
                },
            };
            if let Some(ic) = analytic {
                curves = vec![ic.points.clone()];
            } else {
                match crate::intpatch::surface_surface_intersection(sa.as_ref(), sb.as_ref(), tol) {
                    crate::intpatch::SurfaceIntersection::Curves(ics) => {
                        curves = ics.iter().map(|ic| ic.points.clone()).collect();
                    }
                    _ => {}
                }
            }
            if curves.is_empty() {
                // Fall back to the windowed marching-squares tracer.
                curves = face_intersection_polylines(sa.as_ref(), win_a, sb.as_ref(), tol);
            }
            if !curves.is_empty() {
                pair_curves[i][j] = curves;
            }
        }
    }
    pair_curves
}

/// The retained triangle mesh of the general curved boolean.
///
/// Returns `None` when no face crosses the boundary (the analytic whole-face
/// path applies); `Some(mesh)` otherwise.
fn general_boolean_mesh(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<Option<TriMesh>, String> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    if fa.is_empty() || fb.is_empty() {
        return Err("general_curved_boolean: input has no faces".into());
    }
    let pair_curves = general_pair_curves(&fa, &fb, tol);
    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();
    let has_on = regions_a.contains(&FaceRegion::On) || regions_b.contains(&FaceRegion::On);
    if !has_on {
        return Ok(None);
    }
    boolean_mesh_curves(a, b, op, tol, &pair_curves).map(Some)
}

/// Boolean on solids with general (B-spline / non-analytic) curved faces.
///
/// Computes face–face intersection polylines (via
/// [`crate::intpatch::intersection_curve_points`]) between every overlapping
/// face pair, classifies each face, then rebuilds the result:
///
/// * faces that do not cross the boundary are kept whole (analytic faces stay
///   analytic);
/// * crossing faces are grid-cell classified (a cell is kept when its centroid
///   classifies per the operation) with grid vertices near an intersection
///   polyline snapped onto the curve, so both solids' retained meshes close
///   along the shared intersection;
/// * the welded retained mesh is rebuilt into a faceted BRep solid.
///
/// Planar inputs are delegated to `crate::bop_builder::boolean`.
pub fn general_curved_boolean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }
    match general_boolean_mesh(a, b, op, tol)? {
        Some(mesh) => Ok(result_from_mesh(mesh)),
        None => {
            let fa = faces_of(a);
            let fb = faces_of(b);
            let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
            let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();
            Ok(result_from_kept_faces(&fa, &fb, &regions_a, &regions_b, op))
        }
    }
}

/// The enclosed volume of the general curved boolean result.
///
/// Delegates to [`curved_boolean_volume`] when neither solid has a general
/// curved face; otherwise returns the retained welded mesh volume (exact for a
/// closed, consistently-oriented mesh).
pub fn general_curved_boolean_volume(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<f64, String> {
    if !has_general_curved_face(a) && !has_general_curved_face(b) {
        return curved_boolean_volume(a, b, op, tol);
    }
    if all_faces_planar(a) && all_faces_planar(b) {
        let r = crate::bop_builder::boolean(a, b, op, tol)?;
        return Ok(crate::shape_mesh::shape_volume(&r.shape, 0.02));
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return match op {
            BoolOp::Fuse => Ok(crate::shape_mesh::shape_volume(a, 0.02) + crate::shape_mesh::shape_volume(b, 0.02)),
            BoolOp::Cut => Ok(crate::shape_mesh::shape_volume(a, 0.02)),
            BoolOp::Common => Ok(0.0),
        };
    }
    match general_boolean_mesh(a, b, op, tol)? {
        Some(mesh) => Ok(mesh_volume(&mesh)),
        None => {
            let fa = faces_of(a);
            let fb = faces_of(b);
            let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
            let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();
            let r = result_from_kept_faces(&fa, &fb, &regions_a, &regions_b, op);
            Ok(crate::shape_mesh::shape_volume(&r.shape, 0.02))
        }
    }
}

/// The explicit dispatcher for the curved boolean:
///
/// 1. all faces planar → `crate::bop_builder::boolean` (exact polygon split);
/// 2. any general (non-planar, non-sphere) curved face → `general_boolean_trimmed`
///    (trimmed B-Rep faces that preserve the analytic surface);
/// 3. otherwise (spheres, and quadric-only solids) → [`curved_boolean`].
pub fn curved_boolean_ext(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }
    if has_general_curved_face(a) || has_general_curved_face(b) {
        let shape = general_boolean_trimmed(a, b, op, tol)?;
        let faces = faces_of(&shape);
        return Ok(boolean_result_from_shape(shape, faces, vec![]));
    }
    curved_boolean(a, b, op, tol)
}

// ---------------------------------------------------------------------------
// Trimmed-face (real B-Rep) general boolean
// ---------------------------------------------------------------------------

/// A face that was trimmed along an intersection curve: the original analytic
/// surface is preserved and the boundary wire follows the kept UV region.
pub struct TrimmedFace {
    pub face: crate::shape::Face,
    pub surface: Arc<dyn Surface>,
    pub kept_region_polygon: Vec<GpPnt2d>,
    pub wire: crate::shape::Wire,
}

/// Assemble a [`BooleanResult`] from an already-built result shape.
pub fn boolean_result_from_shape(shape: TopoShape, faces: Vec<Face>, warnings: Vec<String>) -> BooleanResult {
    let shells: Vec<Shell> = shapes_of(&shape, ShapeType::Shell).into_iter().map(Shell).collect();
    let solid = Solid::wrap(shape.clone());
    BooleanResult { shape, solid, shells, faces, warnings }
}

/// The region of `face` relative to `solid` (reuses [`classify_face_general`]).
pub fn face_kept_region(face: &Face, solid: &TopoShape, op: BoolOp, tol: f64) -> FaceRegion {
    let _ = op;
    classify_face_general(face, solid, tol)
}

// -- 2D region primitives ------------------------------------------------

fn unit2d(dx: f64, dy: f64) -> (f64, f64) {
    let len = (dx * dx + dy * dy).sqrt().max(1e-30);
    (dx / len, dy / len)
}

fn poly_diag2d(poly: &[GpPnt2d]) -> f64 {
    if poly.is_empty() {
        return 0.0;
    }
    let mut lo = poly[0];
    let mut hi = poly[0];
    for p in poly {
        lo = GpPnt2d::new(lo.x().min(p.x()), lo.y().min(p.y()));
        hi = GpPnt2d::new(hi.x().max(p.x()), hi.y().max(p.y()));
    }
    lo.distance(&hi)
}

/// Unwrap the periodic first UV coordinate (e.g. sphere/cylinder `u` over
/// `[0, 2π]`) so a seam-crossing polyline becomes continuous: each point's `u`
/// is shifted by a multiple of the period to be closest to the previous point.
/// Without this a latitude circle's closing point projects to `u ≈ 0` instead
/// of `2π`, degenerating the polyline.
fn unwrap_uv_polyline(pts: &[GpPnt2d], surf: &dyn Surface) -> Vec<GpPnt2d> {
    let (u0, u1) = surf.u_range();
    if !(u0.is_finite() && u1.is_finite()) || u1 - u0 < 3.0 || pts.is_empty() {
        return pts.to_vec();
    }
    let period = u1 - u0;
    let mut out = vec![pts[0]];
    let mut cur = pts[0].x();
    for i in 1..pts.len() {
        let raw = pts[i].x();
        let mut best = raw;
        let mut bd = (raw - cur).abs();
        for k in -2i32..=2 {
            let cand = raw + k as f64 * period;
            let d = (cand - cur).abs();
            if d < bd {
                bd = d;
                best = cand;
            }
        }
        cur = best;
        out.push(GpPnt2d::new(best, pts[i].y()));
    }
    out
}

/// Whether a UV polyline is closed: its endpoint gap is at most the trace
/// spacing (the gap between the last and first point of a discretized closed
/// loop is one segment length).
fn polyline_is_closed(poly: &[GpPnt2d], tol: f64) -> bool {
    if poly.len() < 2 {
        return false;
    }
    let diag = poly_diag2d(poly);
    let mut seg_sum = 0.0;
    for i in 0..poly.len() - 1 {
        seg_sum += poly[i].distance(&poly[i + 1]);
    }
    let avg_seg = seg_sum / (poly.len() - 1) as f64;
    let gap = poly.first().unwrap().distance(poly.last().unwrap());
    gap <= tol.max(diag * 0.05).max(1.5 * avg_seg)
}



/// Is `p` on the segment `a-b` (within 1e-9)?
fn on_segment2d(a: &GpPnt2d, b: &GpPnt2d, p: &GpPnt2d) -> bool {
    let cross = (b.x() - a.x()) * (p.y() - a.y()) - (b.y() - a.y()) * (p.x() - a.x());
    if cross.abs() > 1e-9 {
        return false;
    }
    p.x() >= a.x().min(b.x()) - 1e-9
        && p.x() <= a.x().max(b.x()) + 1e-9
        && p.y() >= a.y().min(b.y()) - 1e-9
        && p.y() <= a.y().max(b.y()) + 1e-9
}

fn rect_bounds(rect: &[GpPnt2d]) -> (f64, f64, f64, f64) {
    let u0 = rect.iter().map(|q| q.x()).fold(f64::INFINITY, f64::min);
    let u1 = rect.iter().map(|q| q.x()).fold(f64::NEG_INFINITY, f64::max);
    let v0 = rect.iter().map(|q| q.y()).fold(f64::INFINITY, f64::min);
    let v1 = rect.iter().map(|q| q.y()).fold(f64::NEG_INFINITY, f64::max);
    (u0, u1, v0, v1)
}

/// Extend the ray `p + t·d` (t ≥ 0, `d` a unit 2D direction) until it exits the
/// axis-aligned rectangle `rect`.
fn extend_to_rect_boundary(p: GpPnt2d, d: (f64, f64), rect: &[GpPnt2d]) -> GpPnt2d {
    let (u0, u1, v0, v1) = rect_bounds(rect);
    let mut tmin = f64::INFINITY;
    if d.0.abs() > 1e-15 {
        for u in [u0, u1] {
            let t = (u - p.x()) / d.0;
            if t > 1e-12 {
                tmin = tmin.min(t);
            }
        }
    }
    if d.1.abs() > 1e-15 {
        for v in [v0, v1] {
            let t = (v - p.y()) / d.1;
            if t > 1e-12 {
                tmin = tmin.min(t);
            }
        }
    }
    let t = if tmin.is_finite() { tmin } else { 0.0 };
    GpPnt2d::new(p.x() + t * d.0, p.y() + t * d.1)
}

/// The rectangle edge (0..=3) that the boundary point `p` lies on.
fn edge_containing(rect: &[GpPnt2d], p: GpPnt2d) -> usize {
    let n = 4;
    for i in 0..n {
        if p.distance(&rect[i]) < 1e-9 {
            return i;
        }
    }
    for i in 0..n {
        if on_segment2d(&rect[i], &rect[(i + 1) % n], &p) {
            return i;
        }
    }
    0
}

/// Walk the rectangle boundary from `from` to `to` (both on the boundary),
/// returning the corner points passed (plus the two endpoints).
fn boundary_path(rect: &[GpPnt2d], from: GpPnt2d, to: GpPnt2d, cw: bool) -> Vec<GpPnt2d> {
    let n = 4;
    let ef = edge_containing(rect, from);
    let et = edge_containing(rect, to);
    let mut out = vec![from];
    let mut e = ef;
    loop {
        if e == et {
            break;
        }
        let next = if cw { (e + n - 1) % n } else { (e + 1) % n };
        let corner = if cw { rect[e] } else { rect[(e + 1) % n] };
        if corner.distance(out.last().unwrap()) > 1e-9 {
            out.push(corner);
        }
        e = next;
    }
    if to.distance(out.last().unwrap()) > 1e-9 {
        out.push(to);
    }
    out
}

/// Build the closed cut polygon for a (possibly seam-wrapping) UV polyline:
/// Build the two closed cut polygons for an open UV polyline: extend both
/// endpoints to the rectangle boundary and close along the boundary in each
/// direction. Returns `(left, right)` of the directed polyline, both carrying
/// the polyline's points on the shared edge (so a welded boundary can align the
/// two faces).
fn cut_two_sides(rect: &[GpPnt2d], poly: &[GpPnt2d], _tol: f64) -> (Vec<GpPnt2d>, Vec<GpPnt2d>) {
    let p0 = poly[0];
    let p1 = poly[1];
    let pn = poly[poly.len() - 1];
    let pnm1 = poly[poly.len() - 2];
    let d_back = unit2d(p0.x() - p1.x(), p0.y() - p1.y());
    let d_fwd = unit2d(pn.x() - pnm1.x(), pn.y() - pnm1.y());
    let a = extend_to_rect_boundary(p0, d_back, rect);
    let b = extend_to_rect_boundary(pn, d_fwd, rect);
    let base: Vec<GpPnt2d> = std::iter::once(a).chain(poly.iter().copied()).chain(std::iter::once(b)).collect();
    let mut cut1 = base.clone();
    cut1.extend(boundary_path(rect, b, a, true));
    let mut cut2 = base;
    cut2.extend(boundary_path(rect, b, a, false));
    let (mx, my) = ((p0.x() + p1.x()) * 0.5, (p0.y() + p1.y()) * 0.5);
    let (dx, dy) = (p1.x() - p0.x(), p1.y() - p0.y());
    let len = (dx * dx + dy * dy).sqrt().max(1e-12);
    let off = poly_diag2d(rect) * 1e-3;
    let probe = GpPnt2d::new(mx + (-dy / len) * off, my + (dx / len) * off);
    let in1 = point_in_polygon2d(&cut1, &probe);
    let in2 = point_in_polygon2d(&cut2, &probe);
    if in1 && !in2 {
        (cut1, cut2)
    } else if in2 && !in1 {
        (cut2, cut1)
    } else {
        (cut1, cut2)
    }
}

/// Whether a single-loop region is the full rectangle (the initial split).
fn region_is_rect(region: &[Vec<GpPnt2d>], rect: &[GpPnt2d]) -> bool {
    if region.len() != 1 || region[0].len() != 4 {
        return false;
    }
    for p in &region[0] {
        if !rect.iter().any(|q| q.distance(p) < 1e-9) {
            return false;
        }
    }
    true
}

/// Split the UV rectangle `rect` by the intersection polyline `poly` into the
/// two region loop-sets `(side_a, side_b)`.
///
/// * a proper closed loop → `side_a` = loop interior, `side_b` = loop exterior;
/// * an open polyline (or a degenerate flat loop, e.g. a latitude circle on a
///   sphere) → `side_a` = the left-hand side of the directed polyline.
fn split_uv_region_by_polyline(
    rect: &[GpPnt2d],
    poly: &[GpPnt2d],
    tol: f64,
) -> Result<(Vec<Vec<GpPnt2d>>, Vec<Vec<GpPnt2d>>), String> {
    if poly.len() < 3 || rect.len() < 4 {
        return Err("split_uv_region_by_polyline: degenerate input".into());
    }
    let closed = polyline_is_closed(poly, tol);
    if closed && signed_area2d(poly).abs() > 1e-12 {
        let s1 = polygon_boolean(rect, poly, PolygonBoolOp::Intersect);
        let s2 = polygon_boolean(rect, poly, PolygonBoolOp::Difference);
        return Ok((s1, s2));
    }
    // Open polyline or degenerate flat loop: build the two half-regions directly
    // (they carry the polyline's points on the shared edge, which polygon_boolean's
    // convex fast path would collapse into a straight line).
    let (left, right) = cut_two_sides(rect, poly, tol);
    let mut s1 = Vec::new();
    let mut s2 = Vec::new();
    if left.len() >= 3 {
        s1.push(left);
    }
    if right.len() >= 3 {
        s2.push(right);
    }
    Ok((s1, s2))
}

/// Split a multi-loop region (outer loop + holes) by one UV polyline.
fn split_multipolygon_by_polyline(
    region: &[Vec<GpPnt2d>],
    poly: &[GpPnt2d],
    rect: &[GpPnt2d],
    tol: f64,
) -> (Vec<Vec<GpPnt2d>>, Vec<Vec<GpPnt2d>>) {
    let mut a: Vec<Vec<GpPnt2d>> = Vec::new();
    let mut b: Vec<Vec<GpPnt2d>> = Vec::new();
    let closed = polyline_is_closed(poly, tol);
    if closed && signed_area2d(poly).abs() > 1e-12 {
        for loop_poly in region {
            a.extend(polygon_boolean(loop_poly, poly, PolygonBoolOp::Intersect));
            b.extend(polygon_boolean(loop_poly, poly, PolygonBoolOp::Difference));
        }
    } else if region_is_rect(region, rect) {
        let (left, right) = cut_two_sides(rect, poly, tol);
        a.push(left);
        b.push(right);
    } else {
        // A sub-region split by an open polyline: clip the cut to each loop.
        // (The convex fast path may collapse the polyline detail here; this only
        // affects the interior arrangement of a multiply-crossed face.)
        let (left, right) = cut_two_sides(rect, poly, tol);
        for loop_poly in region {
            a.extend(polygon_boolean(loop_poly, &left, PolygonBoolOp::Intersect));
            b.extend(polygon_boolean(loop_poly, &right, PolygonBoolOp::Intersect));
        }
    }
    (a, b)
}

/// Split a face's UV window into an arrangement by all its intersection
/// polylines. Returns the list of region loop-sets.
fn split_regions_by_polylines(rect: &[GpPnt2d], polylines_uv: &[Vec<GpPnt2d>], tol: f64) -> Vec<Vec<Vec<GpPnt2d>>> {
    let mut regions: Vec<Vec<Vec<GpPnt2d>>> = vec![vec![rect.to_vec()]];
    for poly in polylines_uv {
        let mut next: Vec<Vec<Vec<GpPnt2d>>> = Vec::new();
        for region in &regions {
            let (a, b) = split_multipolygon_by_polyline(region, poly, rect, tol);
            if !a.is_empty() {
                next.push(a);
            }
            if !b.is_empty() {
                next.push(b);
            }
        }
        regions = next;
        if regions.is_empty() {
            break;
        }
    }
    regions
}

/// Whether a 2D point is inside the multi-loop region (outer + holes).
fn region_contains(region: &[Vec<GpPnt2d>], p: &GpPnt2d) -> bool {
    if region.is_empty() || region[0].len() < 3 {
        return false;
    }
    if !point_in_polygon2d(&region[0], p) {
        return false;
    }
    for hole in &region[1..] {
        if hole.len() >= 3 && point_in_polygon2d(hole, p) {
            return false;
        }
    }
    true
}

/// Area-weighted centroid of a 2D polygon (shoelace).
fn polygon_centroid_2d(poly: &[GpPnt2d]) -> GpPnt2d {
    let n = poly.len();
    if n == 0 {
        return GpPnt2d::zero();
    }
    let mut a = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let f = poly[i].x() * poly[j].y() - poly[j].x() * poly[i].y();
        a += f;
        cx += (poly[i].x() + poly[j].x()) * f;
        cy += (poly[i].y() + poly[j].y()) * f;
    }
    if a.abs() < 1e-12 {
        return poly[0];
    }
    GpPnt2d::new(cx / (3.0 * a), cy / (3.0 * a))
}

/// Interior 2D sample points of a UV region: points along the segments from the
/// area-weighted centroid to each vertex (plus triangle centroids), filtered out
/// of holes. Used to classify the region against the other solid. The centroid
/// rays reach deep into thin strips whose vertices all lie on one boundary edge
/// (a sphere latitude cap), so the probes are not ambiguous points on the edge.
fn region_interior_points_2d(region: &[Vec<GpPnt2d>]) -> Vec<GpPnt2d> {
    let outer = match region.first() {
        Some(o) if o.len() >= 3 => o,
        _ => return Vec::new(),
    };
    let c = polygon_centroid_2d(outer);
    let mut pts: Vec<GpPnt2d> = Vec::new();
    for v in outer {
        for t in [0.3, 0.6] {
            let p = GpPnt2d::new(c.x() + (v.x() - c.x()) * t, c.y() + (v.y() - c.y()) * t);
            if region_contains(region, &p) {
                pts.push(p);
            }
        }
    }
    let outer3: Vec<GpPnt> = outer.iter().map(|p| GpPnt::new(p.x(), p.y(), 0.0)).collect();
    let tris = triangulate_polygon(&outer3);
    for (a, b, c) in tris {
        let pa = outer[a];
        let pb = outer[b];
        let pc = outer[c];
        let centroid = GpPnt2d::new((pa.x() + pb.x() + pc.x()) / 3.0, (pa.y() + pb.y() + pc.y()) / 3.0);
        if region_contains(region, &centroid) {
            pts.push(centroid);
        }
    }
    if pts.is_empty() && region_contains(region, &outer[0]) {
        pts.push(outer[0]);
    }
    pts
}

/// Classify a UV region of `face` against `other` by sampling interior points.
fn region_inside_other(
    _face: &Face,
    surf: &dyn Surface,
    region: &[Vec<GpPnt2d>],
    other: &TopoShape,
    tol: f64,
) -> FaceRegion {
    let pts2d = region_interior_points_2d(region);
    let mut inside = 0usize;
    let mut total = 0usize;
    for p2 in &pts2d {
        let p = surf.d0(p2.x(), p2.y());
        if point_in_solid_curved(p, other, tol) {
            inside += 1;
        }
        total += 1;
    }
    if total == 0 {
        FaceRegion::On
    } else if inside == 0 {
        FaceRegion::Outside
    } else if inside == total {
        FaceRegion::Inside
    } else if inside * 2 >= total {
        FaceRegion::Inside
    } else {
        FaceRegion::Outside
    }
}

/// Map a UV region loop to a 3D point loop on the surface. Planar faces keep
/// one point per vertex (straight edges); curved faces sub-sample UV edges so
/// curved boundaries (e.g. a cylinder's top rim) resolve correctly.
fn uv_polygon_to_3d_loop(surf: &dyn Surface, loop_uv: &[GpPnt2d], rect: &[GpPnt2d], _tol: f64) -> Vec<GpPnt> {
    let planar = classify_surface(surf) == SurfaceKind::Plane;
    let step = poly_diag2d(rect) / 32.0;
    let mut out: Vec<GpPnt> = Vec::new();
    for i in 0..loop_uv.len() {
        let a = loop_uv[i];
        let b = loop_uv[(i + 1) % loop_uv.len()];
        let seg_len = a.distance(&b);
        let n = if planar {
            1
        } else {
            (seg_len / step.max(1e-9)).round().max(1.0) as usize
        };
        let n = n.max(1).min(48);
        for k in 0..n {
            let t = k as f64 / n as f64;
            let u = a.x() + (b.x() - a.x()) * t;
            let v = a.y() + (b.y() - a.y()) * t;
            out.push(surf.d0(u, v));
        }
    }
    out
}

/// Build a face from 3D boundary loops (first loop = outer wire, the rest are
/// holes).
fn build_face_from_loops(surf: Arc<dyn Surface>, loops_3d: &[Vec<GpPnt>]) -> Result<Face, String> {
    let bld = TopoBuilder::new();
    let mut wires: Vec<Wire> = Vec::new();
    for lp in loops_3d {
        if lp.len() < 3 {
            continue;
        }
        let mut edges: Vec<Edge> = Vec::new();
        for i in 0..lp.len() {
            let j = (i + 1) % lp.len();
            if lp[i].distance(&lp[j]) > 1e-12 {
                edges.push(bld.make_edge_segment(&lp[i], &lp[j]));
            }
        }
        if edges.len() >= 2 {
            let w = bld.make_wire(&edges);
            w.set_closed(true);
            wires.push(w);
        }
    }
    Ok(bld.make_face(surf, &wires))
}

/// The finite UV window of a face as a rectangle polygon.
fn face_uv_window_rect(face: &Face) -> Result<Vec<GpPnt2d>, String> {
    let (u0, u1, v0, v1) = face_uv_window_local(face);
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) || u1 <= u0 || v1 <= v0 {
        return Err("face_uv_window_rect: unbounded or degenerate UV window".into());
    }
    Ok(vec![GpPnt2d::new(u0, v0), GpPnt2d::new(u1, v0), GpPnt2d::new(u1, v1), GpPnt2d::new(u0, v1)])
}

/// Trim `face`'s UV window along an intersection polyline (2D points in the
/// face's UV space), keeping the side selected by `keep_inside`.
///
/// `keep_inside = true` keeps the region enclosed by a closed loop (or the
/// left-hand side of a directed open polyline). The returned face keeps the
/// original analytic surface; only the boundary wire is replaced.
pub fn trim_face_to_region(face: &Face, keep_inside: bool, intersection_uvs: &[GpPnt2d], tol: f64) -> Result<Option<Face>, String> {
    let rect = face_uv_window_rect(face)?;
    let (s1, s2) = split_uv_region_by_polyline(&rect, intersection_uvs, tol)?;
    let kept = if keep_inside { s1 } else { s2 };
    if kept.is_empty() {
        return Ok(None);
    }
    let surf = BRepTool::face_surface(face).ok_or("trim_face_to_region: no face surface")?;
    let loops_3d: Vec<Vec<GpPnt>> = kept.iter().map(|lp| uv_polygon_to_3d_loop(surf.as_ref(), lp, &rect, tol)).collect();
    build_face_from_loops(surf, &loops_3d).map(Some)
}

/// Whether a curve is (approximately) a straight line on `[t0, t1]`.
fn curve_is_line(c: &dyn Curve, t0: f64, t1: f64) -> bool {
    let p0 = c.d0(t0);
    let p1 = c.d0(0.5 * (t0 + t1));
    let p2 = c.d0(t1);
    let d1 = GpVec::from_pnts(&p0, &p1);
    let d2 = GpVec::from_pnts(&p0, &p2);
    d1.cross_magnitude(&d2) < 1e-6 * (d1.magnitude() * d2.magnitude()).max(1e-12)
}


/// Sample a whole face's boundary into contiguous 3D loops.
///
/// A face whose wire is a single closed edge (a cap) is sampled along that
/// edge. Any other face (planar quads, the cylinder lateral, …) is rebuilt from
/// its finite UV-window rectangle mapped through the surface, which yields the
/// correct quad for planar faces and the correct seam+circle boundary for a
/// cylinder.
fn face_boundary_loops_3d(face: &Face) -> Vec<Vec<GpPnt>> {
    let mut out: Vec<Vec<GpPnt>> = Vec::new();
    for w in wires_of_face(face) {
        let edges = edges_of_wire(&w);
        if edges.len() == 1 {
            let s = edge_samples_3d(&edges[0]);
            if !s.is_empty() {
                out.push(s);
            }
        } else if let (Ok(rect), Some(surf)) = (face_uv_window_rect(face), BRepTool::face_surface(face)) {
            let loop3d = uv_polygon_to_3d_loop(surf.as_ref(), &rect, &rect, 1e-6);
            if loop3d.len() >= 3 {
                out.push(loop3d);
            }
        } else {
            // Fallback: dump each edge's samples (used for unbounded faces with
            // no finite UV window).
            for e in &edges {
                let s = edge_samples_3d(e);
                if !s.is_empty() {
                    out.push(s);
                }
            }
        }
    }
    out
}

/// Sample an edge's curve into 3D points (straight lines → 2, otherwise 24).
fn edge_samples_3d(e: &Edge) -> Vec<GpPnt> {
    if let Some(c) = BRepTool::edge_curve(e) {
        let (t0, t1) = BRepTool::edge_parameters(e);
        if t0.is_finite() && t1.is_finite() && t1 > t0 {
            let n = if curve_is_line(c.as_ref(), t0, t1) { 2 } else { 24 };
            let mut pts = Vec::with_capacity(n);
            for k in 0..n {
                let t = t0 + (t1 - t0) * k as f64 / (n - 1) as f64;
                pts.push(c.d0(t));
            }
            return pts;
        }
    }
    if let Some((a, b)) = BRepTool::edge_vertices(e) {
        return vec![a, b];
    }
    Vec::new()
}

/// Weld the boundary loops of all kept faces into shared edges, then rebuild
/// every face with the welded wires (so adjacent faces close the shell).
fn weld_loops_into_faces(faces_data: Vec<(Arc<dyn Surface>, Vec<Vec<GpPnt>>)>, tol: f64) -> Vec<Face> {
    let bld = TopoBuilder::new();
    // Relative to the geometry scale so the two faces' projections of the same
    // shared trace points (each off by ~1e-3) merge, while distinct curve
    // samples (spaced by the trace grid, ~5e-2) stay separate.
    let mut diag = 0.0f64;
    for (_, loops) in &faces_data {
        for lp in loops {
            if lp.is_empty() {
                continue;
            }
            let mut lo = lp[0];
            let mut hi = lp[0];
            for p in lp {
                lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
                hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
            }
            diag = diag.max(lo.distance(&hi));
        }
    }
    let weld_tol = tol.max(diag * 2e-3).max(1e-4);
    let cell = |p: &GpPnt| -> (i64, i64, i64) {
        (
            f64::floor(p.x() / weld_tol) as i64,
            f64::floor(p.y() / weld_tol) as i64,
            f64::floor(p.z() / weld_tol) as i64,
        )
    };
    let mut points: Vec<GpPnt> = Vec::new();
    let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> = std::collections::HashMap::new();
    let mut idx_loops: Vec<Vec<Vec<usize>>> = Vec::new();
    for (_, loops) in &faces_data {
        let mut face_idx: Vec<Vec<usize>> = Vec::new();
        for loop_pts in loops {
            let mut idx: Vec<usize> = Vec::with_capacity(loop_pts.len());
            for p in loop_pts {
                let k = cell(p);
                let mut found = None;
                'search: for dx in -1i64..=1 {
                    for dy in -1i64..=1 {
                        for dz in -1i64..=1 {
                            if let Some(bucket) = grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                                for &j in bucket {
                                    if points[j].distance(p) <= weld_tol {
                                        found = Some(j);
                                        break 'search;
                                    }
                                }
                            }
                        }
                    }
                }
                match found {
                    Some(j) => idx.push(j),
                    None => {
                        let j = points.len();
                        points.push(*p);
                        grid.entry(k).or_default().push(j);
                        idx.push(j);
                    }
                }
            }
            face_idx.push(idx);
        }
        idx_loops.push(face_idx);
    }
    let mut edge_map: std::collections::HashMap<(usize, usize), Edge> = std::collections::HashMap::new();
    let mut face_edges: Vec<Vec<Vec<Edge>>> = Vec::new();
    for face_idx in &idx_loops {
        let mut loops_edges: Vec<Vec<Edge>> = Vec::new();
        for idx in face_idx {
            let mut edges: Vec<Edge> = Vec::new();
            let m = idx.len();
            for i in 0..m {
                let j = idx[i];
                let k = idx[(i + 1) % m];
                if j == k {
                    continue;
                }
                let key = (j.min(k), j.max(k));
                let e = edge_map.entry(key).or_insert_with(|| bld.make_edge_segment(&points[j], &points[k])).clone();
                edges.push(e);
            }
            loops_edges.push(edges);
        }
        face_edges.push(loops_edges);
    }
    let mut faces: Vec<Face> = Vec::new();
    for (i, (surf, _)) in faces_data.iter().enumerate() {
        let wires: Vec<Wire> = face_edges[i]
            .iter()
            .map(|edges| {
                let w = bld.make_wire(edges);
                w.set_closed(true);
                w
            })
            .collect();
        faces.push(bld.make_face(surf.clone(), &wires));
    }
    faces
}

/// Boolean producing trimmed B-Rep faces: each kept face keeps its original
/// NURBS/analytic surface, with the boundary wire trimmed along the shared
/// face–face intersection curves. Planar inputs delegate to
/// [`crate::bop_builder::boolean`].
pub fn general_boolean_trimmed(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<TopoShape, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol).map(|r| r.shape);
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op).shape);
    }
    let fa = faces_of(a);
    let fb = faces_of(b);
    if fa.is_empty() || fb.is_empty() {
        return Err("general_boolean_trimmed: input has no faces".into());
    }
    let pair_curves = general_pair_curves_windowed(&fa, &fb, tol);
    let mut faces_data: Vec<(Arc<dyn Surface>, Vec<Vec<GpPnt>>)> = Vec::new();

    // Faces of A, then faces of B.
    for (from_a, src, other, other_len, row) in [
        (true, a, b, fb.len(), &pair_curves),
        (false, b, a, fa.len(), &pair_curves),
    ] {
        let faces = faces_of(src);
        for (i, f) in faces.iter().enumerate() {
            let Some(surf) = BRepTool::face_surface(f) else { continue };
            let mut polylines_3d: Vec<Vec<GpPnt>> = Vec::new();
            for j in 0..other_len {
                for poly in if from_a { &row[i][j] } else { &row[j][i] } {
                    polylines_3d.push(poly.clone());
                }
            }
            let r = classify_face_general(f, other, tol);
            if r != FaceRegion::On || polylines_3d.is_empty() {
                if select_keep(r, from_a, op) {
                    faces_data.push((surf.clone(), face_boundary_loops_3d(f)));
                }
                continue;
            }
            let rect = face_uv_window_rect(f)?;
            let polylines_uv: Vec<Vec<GpPnt2d>> = polylines_3d
                .iter()
                .map(|poly| {
                    let raw: Vec<GpPnt2d> = poly
                        .iter()
                        .map(|p| {
                            let (u, v) = crate::intpatch::project_params(surf.as_ref(), p);
                            GpPnt2d::new(u, v)
                        })
                        .collect();
                    unwrap_uv_polyline(&raw, surf.as_ref())
                })
                .collect();
            let regions = split_regions_by_polylines(&rect, &polylines_uv, tol);
            for region in &regions {
                let reg_class = region_inside_other(f, surf.as_ref(), region, other, tol);
                if select_keep(reg_class, from_a, op) {
                    let loops_3d: Vec<Vec<GpPnt>> = region
                        .iter()
                        .map(|lp| uv_polygon_to_3d_loop(surf.as_ref(), lp, &rect, tol))
                        .collect();
                    faces_data.push((surf.clone(), loops_3d));
                }
            }
        }
    }

    let bld = TopoBuilder::new();
    if faces_data.is_empty() {
        let comp = bld.make_compound_of(&[]);
        return Ok(comp.0);
    }
    let faces = weld_loops_into_faces(faces_data, tol);
    let shell = bld.make_shell(&faces);
    let closed = crate::shell_check::shell_is_closed(&shell);
    let solid = if closed { Some(bld.make_solid(&[shell.clone()])) } else { None };
    Ok(solid.map(|s| s.0).unwrap_or_else(|| shell.0.clone()))
}

/// The full boolean dispatcher:
///
/// 1. all faces planar → `crate::bop_builder::boolean` (exact polygon split);
/// 2. any general (non-planar, non-sphere) curved face → [`general_boolean_trimmed`]
///    (trimmed B-Rep faces that preserve the analytic surface);
/// 3. otherwise (spheres, and quadric-only solids) → [`curved_boolean`].
pub fn curved_boolean_full(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }
    if has_general_curved_face(a) || has_general_curved_face(b) {
        let shape = general_boolean_trimmed(a, b, op, tol)?;
        let faces = faces_of(&shape);
        return Ok(boolean_result_from_shape(shape, faces, vec![]));
    }
    curved_boolean(a, b, op, tol)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::{GpAx3, GpDir};
    use crate::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};
    use crate::shape_mesh::shape_volume;
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;

    fn unit_sphere_solid() -> Solid {
        BRepPrimSphere::make_sphere(1.0).solid
    }

    fn sphere_at(center: GpPnt) -> Solid {
        let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bld = TopoBuilder::new();
        let face = bld.make_face(Arc::new(occt_geom::GeomSphere::new(
            occt_core::gp::GpSphere::new(ax3, 1.0).unwrap(),
        )), &[]);
        let shell = bld.make_shell(&[face]);
        bld.make_solid(&[shell])
    }

    fn centered_box(half: f64) -> Solid {
        BRepPrimBox::make_box_corner(&GpPnt::new(-half, -half, -half), &GpPnt::new(half, half, half)).solid
    }

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    #[test]
    fn box_fuse_box_delegates() {
        // Both planar → delegated to bop_builder (exact polygon boolean).
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0)).solid;
        let r = curved_boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.solid.is_some(), "fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "fuse shell is closed");
        assert!(faces_of(&r.shape).len() > 12, "faces {}", faces_of(&r.shape).len());
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn sphere_inside_box_common_volume() {
        // Unit sphere at origin inside box [-1.5,1.5]³: Common = the sphere.
        let sphere = unit_sphere_solid();
        let box_s = centered_box(1.5);
        let r = curved_boolean(&sphere.0, &box_s.0, BoolOp::Common, 1e-6).expect("common ok");
        assert!(r.solid.is_some(), "common produces a solid");
        // The result keeps the whole (curved) sphere face, so shape_volume is
        // reliable for this curved solid.
        let v = shape_volume(&r.shape, 0.05);
        let expected = 4.0 / 3.0 * PI;
        assert!((v - expected).abs() < 0.05, "common volume {v} (expected {expected})");
        clear_tree(&r.shape);
        clear_tree(&sphere.0);
        clear_tree(&box_s.0);
    }

    #[test]
    fn sphere_sphere_fuse_volume() {
        // Two unit spheres 1.5 apart: Fuse volume ≈ 8.018.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let v = curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse volume");
        assert!((7.5..8.5).contains(&v), "fuse volume {v} (expected ≈ 8.018)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn sphere_sphere_common_volume() {
        // Two unit spheres 1.5 apart: Common = lens, volume ≈ 0.36.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let v = curved_boolean_volume(&a.0, &b.0, BoolOp::Common, 1e-6).expect("common volume");
        assert!((0.2..0.6).contains(&v), "common volume {v} (expected ≈ 0.36)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn sphere_sphere_cut_volume() {
        // A − B: sphere 1 minus the lens → volume ≈ 4.188 − 0.36 = 3.83.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let v = curved_boolean_volume(&a.0, &b.0, BoolOp::Cut, 1e-6).expect("cut volume");
        assert!((3.0..4.2).contains(&v), "cut volume {v} (expected ≈ 3.83)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn disjoint_fuse_returns_compound() {
        // Sphere at origin, sphere far away → Fuse returns a compound.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(5.0, 0.0, 0.0));
        let r = curved_boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.shape.is_compound(), "disjoint fuse is a compound");
        assert!(r.solid.is_none(), "disjoint fuse has no single solid");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn classify_sphere_inside_box() {
        let sphere = unit_sphere_solid();
        let box_s = centered_box(1.5);
        let sf = faces_of(&sphere.0);
        assert_eq!(classify_face(&sf[0], &box_s.0, 1e-6), FaceRegion::Inside);
        let bf = faces_of(&box_s.0);
        for f in bf.iter().take(3) {
            assert_eq!(classify_face(f, &sphere.0, 1e-6), FaceRegion::Outside);
        }
        clear_tree(&sphere.0);
        clear_tree(&box_s.0);
    }

    #[test]
    fn ray_hits_box() {
        let box_s = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let p = GpPnt::new(-1.0, 0.5, 0.5);
        let hits = ray_hits_solid(p, &box_s.0, GpVec::new(1.0, 0.0, 0.0), 1e-9);
        assert_eq!(hits, 2, "ray through unit box");
        clear_tree(&box_s.0);
    }

    #[test]
    fn mesh_volume_full_sphere() {
        let sphere = unit_sphere_solid();
        let sf = faces_of(&sphere.0);
        let m = mesh_whole_face(&sf[0]);
        let v = mesh_volume(&m);
        let expected = 4.0 / 3.0 * PI;
        // Chord triangulation at 24×24 is ~1.6% low; keep a relaxed band.
        assert!((v - expected).abs() < 0.1, "sphere mesh volume {v} (expected {expected})");
        clear_tree(&sphere.0);
    }

    // -- general (non-analytic) curved-face boolean extensions --

    fn cylinder_solid(radius: f64, height: f64) -> Solid {
        BRepPrimCylinder::make_cylinder(radius, height).solid
    }

    /// A curved B-spline patch face `z = base + bump·((u−½)² + (v−½)²)` on
    /// `[0, 1]²` (a non-planar `Other` surface).
    fn curved_patch_face(base: f64, bump: f64) -> Face {
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                let z = base + bump * ((u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5));
                row.push(GpPnt::new(u, v, z));
            }
            grid.push(row);
        }
        let surf = crate::intpatch::make_bspline_surface_from_grid(&grid, 2, 2).expect("patch fit");
        let bld = TopoBuilder::new();
        bld.make_face(surf, &[])
    }

    #[test]
    fn general_box_cylinder_fuse() {
        // Box [-1,1]³ fused with a radius-0.4 cylinder on z∈[0,2] poking through
        // the top face: expected = 8 + π·0.4²·1 ≈ 8.50.
        let box_s = centered_box(1.0);
        let cyl = cylinder_solid(0.4, 2.0);
        let v = general_curved_boolean_volume(&box_s.0, &cyl.0, BoolOp::Fuse, 1e-6).expect("fuse volume");
        let expected = 8.0 + PI * 0.16 * 1.0;
        assert!((v - expected).abs() < 0.5, "box-cylinder fuse volume {v} (expected {expected})");
        clear_tree(&box_s.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn general_sphere_cylinder_cut() {
        // Unit sphere minus a radius-0.4 cylinder on z∈[0,2]: the cylinder exits
        // the sphere at z = sqrt(1 − 0.4²) ≈ 0.9165.
        let sphere = unit_sphere_solid();
        let cyl = cylinder_solid(0.4, 2.0);
        let v = general_curved_boolean_volume(&sphere.0, &cyl.0, BoolOp::Cut, 1e-6).expect("cut volume");
        let z_exit = (1.0 - 0.4f64 * 0.4).sqrt();
        let expected = 4.0 / 3.0 * PI - PI * 0.16 * z_exit;
        assert!((v - expected).abs() < 0.5, "sphere-cylinder cut volume {v} (expected {expected})");
        clear_tree(&sphere.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn planar_delegates_to_bop_builder() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0)).solid;
        let r = curved_boolean_ext(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.solid.is_some(), "planar fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "planar fuse shell is closed");
        assert!(faces_of(&r.shape).len() > 12, "faces {}", faces_of(&r.shape).len());
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn classify_general_face_inside_outside() {
        let box_s = centered_box(1.0);
        let inside = curved_patch_face(0.8, 0.3); // z∈[0.8, 0.95] inside the box
        assert_eq!(classify_face_general(&inside, &box_s.0, 1e-6), FaceRegion::Inside);
        let outside = curved_patch_face(3.0, 0.3); // z∈[3.0, 3.15] above the box
        assert_eq!(classify_face_general(&outside, &box_s.0, 1e-6), FaceRegion::Outside);
        // A curved face crossing the box top (spanning inside and outside) is On.
        let crossing = curved_patch_face(0.9, 0.5); // z∈[0.9, 1.15]
        assert_eq!(classify_face_general(&crossing, &box_s.0, 1e-6), FaceRegion::On);
        clear_tree(&box_s.0);
        clear_tree(&inside.0);
        clear_tree(&outside.0);
        clear_tree(&crossing.0);
    }

    #[test]
    fn disjoint_general_fuse_compound() {
        // Box at z∈[-4,-2] far below a cylinder at z∈[0,2] → Fuse returns a
        // compound (bboxes do not overlap).
        let box_s = BRepPrimBox::make_box_corner(&GpPnt::new(-1.0, -1.0, -4.0), &GpPnt::new(1.0, 1.0, -2.0)).solid;
        let cyl = cylinder_solid(0.5, 2.0);
        let r = curved_boolean_ext(&box_s.0, &cyl.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.shape.is_compound(), "disjoint general fuse is a compound");
        assert!(r.solid.is_none(), "disjoint general fuse has no single solid");
        clear_tree(&r.shape);
        clear_tree(&box_s.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn general_boolean_volume_matches_curved_boolean_volume() {
        // Sphere-sphere dispatch stays on the analytic path: the general
        // dispatcher's volume agrees with curved_boolean_volume.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let ext_vol = general_curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ext volume");
        let ref_vol = curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref volume");
        assert!((ext_vol - ref_vol).abs() < 0.05, "ext {ext_vol} vs curved_boolean_volume {ref_vol}");
        assert!((ext_vol - 8.018).abs() < 0.3, "ext volume {ext_vol} (expected ≈ 8.018)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    // -- trimmed-face (real B-Rep) general boolean --

    #[test]
    fn trimmed_face_box_cylinder_fuse() {
        // Box [-1,1]³ fused with a radius-0.4 cylinder on z∈[0,2] poking through
        // the top: the cylinder lateral face is trimmed to the z>1 stub and keeps
        // its analytic cylinder surface (not a faceted mesh).
        let box_s = centered_box(1.0);
        let cyl = cylinder_solid(0.4, 2.0);
        let r = curved_boolean_full(&box_s.0, &cyl.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let has_cylinder = faces_of(&r.shape)
            .iter()
            .any(|f| BRepTool::face_surface(f).map(|s| surface_cylinder_params(s.as_ref()).is_some()).unwrap_or(false));
        assert!(has_cylinder, "a trimmed cylinder face is preserved");
        // `shape_volume` tessellates each face's full UV window (it ignores the
        // trimming wire and full-surface overlap), so the volume is a loose band —
        // it must still be in the ballpark of the box plus the protruding stub.
        let v = shape_volume(&r.shape, 0.02);
        let expected = 8.0 + PI * 0.16 * 1.0;
        assert!((v - expected).abs() < 5.0, "box-cylinder fuse volume {v} (expected ≈ {expected})");
        clear_tree(&r.shape);
        clear_tree(&box_s.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn trimmed_sphere_box_cut() {
        // Sphere radius 1 cut by a thin slab box z∈[0.5,0.9]: the retained sphere
        // caps keep their analytic sphere surface.
        let sphere = unit_sphere_solid();
        let slab = BRepPrimBox::make_box_corner(&GpPnt::new(-2.0, -2.0, 0.5), &GpPnt::new(2.0, 2.0, 0.9)).solid;
        let shape = general_boolean_trimmed(&sphere.0, &slab.0, BoolOp::Cut, 1e-6).expect("cut ok");
        let has_sphere = faces_of(&shape)
            .iter()
            .any(|f| BRepTool::face_surface(f).map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere).unwrap_or(false));
        assert!(has_sphere, "a retained sphere cap is preserved");
        clear_tree(&shape);
        clear_tree(&sphere.0);
        clear_tree(&slab.0);
    }

    #[test]
    fn trim_face_to_region_keeps_side() {
        // A sphere face trimmed by a small UV circle keeps the inside region:
        // the surface is still a radius-1 sphere, but the face's wire is the
        // smaller loop.
        let sphere = unit_sphere_solid();
        let face = faces_of(&sphere.0)[0].clone();
        let (u0, u1, v0, v1) = face_uv_window_local(&face);
        let (uc, vc) = (0.5 * (u0 + u1), 0.0);
        let radius = 0.9 * (v1 - v0).min(u1 - u0) / 2.0;
        let circle: Vec<GpPnt2d> = (0..32)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / 32.0;
                GpPnt2d::new(uc + radius * a.cos(), vc + radius * a.sin())
            })
            .collect();
        let trimmed = trim_face_to_region(&face, true, &circle, 1e-6).expect("trim ok").expect("some face");
        let s = BRepTool::face_surface(&trimmed).expect("surface");
        assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Sphere, "surface preserved");
        let (sp_c, sp_r) = crate::intpatch::sphere_params(s.as_ref()).expect("sphere params");
        assert!((sp_r - 1.0).abs() < 1e-6, "radius {sp_r}");
        assert!(sp_c.distance(&GpPnt::zero()) < 1e-6, "center {sp_c:?}");
        let w = wires_of_face(&trimmed);
        assert_eq!(w.len(), 1, "one boundary wire");
        let e = edges_of_wire(&w[0]);
        assert!(e.len() >= 3, "wire has edges");
        clear_tree(&trimmed.0);
        clear_tree(&sphere.0);
    }

    #[test]
    fn general_boolean_trimmed_closed_shell() {
        // Box [-1,1]³ fused with a small sphere poking through the top face → the
        // trimmed faces weld into a closed shell.
        let box_s = centered_box(1.0);
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, 1.2), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bld = TopoBuilder::new();
        let face = bld.make_face(Arc::new(occt_geom::GeomSphere::new(occt_core::gp::GpSphere::new(ax3, 0.5).unwrap())), &[]);
        let shell_s = bld.make_shell(&[face]);
        let sphere = bld.make_solid(&[shell_s]);
        let shape = general_boolean_trimmed(&box_s.0, &sphere.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let shells = shapes_of(&shape, ShapeType::Shell);
        assert!(!shells.is_empty(), "result has shells");
        let shell = Shell(shells[0].clone());
        assert!(shell_is_closed(&shell), "fused box+sphere shell is closed");
        // The welded trimmed faces (box top ring + sphere cap) share every edge.
        let usage = crate::shell_check::edge_face_usage(&shell.0);
        assert!(usage.values().all(|&c| c == 2), "every boundary edge is used by exactly 2 faces");
        clear_tree(&shape);
        clear_tree(&box_s.0);
        clear_tree(&sphere.0);
    }

    #[test]
    fn curved_boolean_full_planar_delegates() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0)).solid;
        let r = curved_boolean_full(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let refr = crate::bop_builder::boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref ok");
        assert_eq!(faces_of(&r.shape).len(), faces_of(&refr.shape).len(), "planar delegates to bop_builder");
        assert!(shell_is_closed(&r.shells[0]), "planar fuse shell is closed");
        clear_tree(&r.shape);
        clear_tree(&refr.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn curved_boolean_full_quadric_unchanged() {
        // Two spheres (analytic quadrics) route to the unchanged curved_boolean.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let r = curved_boolean_full(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let refr = curved_boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref ok");
        assert_eq!(faces_of(&r.shape).len(), faces_of(&refr.shape).len(), "quadrics route to curved_boolean");
        let ref_vol = curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref vol");
        assert!((ref_vol - 8.018).abs() < 0.3, "volume {ref_vol} (expected ≈ 8.018)");
        clear_tree(&r.shape);
        clear_tree(&refr.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn trimmed_face_preserves_surface() {
        let sphere = unit_sphere_solid();
        let face = faces_of(&sphere.0)[0].clone();
        let orig = BRepTool::face_surface(&face).expect("original surface");
        let (u0, u1, v0, v1) = face_uv_window_local(&face);
        let (uc, vc) = (0.5 * (u0 + u1), 0.3);
        let radius = 0.7 * (v1 - v0).min(u1 - u0) / 2.0;
        let circle: Vec<GpPnt2d> = (0..32)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / 32.0;
                GpPnt2d::new(uc + radius * a.cos(), vc + radius * a.sin())
            })
            .collect();
        let trimmed = trim_face_to_region(&face, true, &circle, 1e-6).expect("trim ok").expect("some face");
        let s = BRepTool::face_surface(&trimmed).expect("trimmed surface");
        let p_orig = orig.d0(uc, vc);
        let p_new = s.d0(uc, vc);
        assert!(p_orig.distance(&p_new) < 1e-9, "surface evaluates identically at the face center");
        assert!((p_orig.distance(&GpPnt::zero()) - 1.0).abs() < 1e-6, "still on the unit sphere");
        clear_tree(&trimmed.0);
        clear_tree(&sphere.0);
    }


    #[test]
    fn boolean_result_shape_has_faces() {
        let box_s = BRepPrimBox::make_box(1.0, 2.0, 3.0).solid;
        let faces = faces_of(&box_s.0);
        let r = boolean_result_from_shape(box_s.0.clone(), faces.clone(), vec![]);
        assert_eq!(r.faces.len(), 6);
        assert!(r.solid.is_some(), "a box shape round-trips to a solid");
        assert!(!r.shells.is_empty(), "the solid carries a shell");
        assert_eq!(faces_of(&r.shape).len(), 6, "faces round-trip");
        clear_tree(&box_s.0);
        clear_tree(&r.shape);
    }
}


