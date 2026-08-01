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

use occt_core::gp::{GpPnt, GpVec};

use crate::abs::ShapeType;
use crate::bop_builder::{BoolOp, BooleanResult};
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Face, Shell, Solid, TopoShape};
use crate::topo_tools_full::{faces_of, shapes_of};

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

/// Grid-triangulate a face's surface on an `(nu+1)×(nv+1)` UV grid.
fn surface_grid_mesh(surf: &dyn occt_geom::Surface, nu: usize, nv: usize) -> TriMesh {
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(surf);
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
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::{GpAx3, GpDir};
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
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
}
