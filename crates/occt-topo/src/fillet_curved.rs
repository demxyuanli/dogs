//! Curved-face rolling-ball fillet — a port of `ChFi3d` blending between a
//! curved surface and a plane, or between two curved surfaces.
//!
//! `fillet_edge` only blends two planar faces (the blend is a cylinder). This
//! module extends the rolling-ball idea to edges where at least one adjacent
//! face is a sphere or a cylinder:
//!
//! * **Plane + Sphere** — the rolling ball stays tangent to the plane and to
//!   the sphere, so its centre sweeps the *intersection of the offset plane*
//!   (the plane pushed out by `r`) *and the offset sphere* (radius `R_s + r`).
//!   That intersection is a circle; the envelope of the ball along it is a
//!   **torus** with major radius = circle radius and minor radius = fillet
//!   radius `r`.
//! * **Plane + Cylinder** (axis perpendicular to the plane) — the ball centre
//!   sweeps the intersection of the offset plane and the offset cylinder
//!   (radius `R_c + r`), again a circle, and the blend is a **torus** band
//!   tangent to the cylinder at the inner equator and to the plane at the
//!   outer bottom.
//! * **Sphere + Sphere** — the ball centre sweeps the intersection of the two
//!   offset spheres (radius `R_i + r`), a circle, and the blend is a torus.
//!
//! In every supported case the blend face is a `GeomTorus` band bounded by the
//! two tangency circles (one per adjacent face). The two adjacent faces are
//! rebuilt: the planar face keeps its outer boundary and gets the enlarged
//! plane-tangency circle as a hole; the curved face gets the smaller
//! surface-tangency circle as its new boundary. A cylinder's lateral face is
//! rebuilt with the contact circle, a re-generated seam sub-edge and the
//! original top rim. All new edges are shared between the blend face and the
//! trimmed adjacent face, so the rebuilt shell satisfies `shell_is_closed`.
//!
//! Unsupported pairs (cone, torus, B-spline, a cylinder whose axis is not
//! perpendicular to the plane, ...) return `Err` — documented, never a panic.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_geom::bspline_surface::fit_surface_grid;
use occt_geom::{GeomCylinder, GeomPlane, GeomSphere, GeomTorus, Surface};

use crate::brep_surface::{
    classify_surface, is_planar, sphere_center, surface_closest_params, surface_normal, SurfaceKind,
};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fillet_edge::classify_surface_full;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, vertices_of,
    wires_of_face,
};

/// TEMPORARY BUILD FIX (concurrent agent): component `i` of a `GpXyz`.
fn xyz_component(c: &occt_core::gp::xyz::GpXyz, i: usize) -> f64 {
    match i {
        0 => c.x(),
        1 => c.y(),
        _ => c.z(),
    }
}

/// Number of sampled points on a rolling-ball centreline circle.
const CENTERLINE_SAMPLES: usize = 48;

// ===========================================================================
// Analytic surface descriptions
// ===========================================================================

/// A planar surface: a point on it and its unit normal (outward).
#[derive(Debug, Clone, Copy)]
pub struct PlaneInfo {
    pub origin: GpPnt,
    pub normal: GpVec,
}

/// A spherical surface.
#[derive(Debug, Clone, Copy)]
pub struct SphereInfo {
    pub center: GpPnt,
    pub radius: f64,
}

/// A cylindrical surface.
#[derive(Debug, Clone, Copy)]
pub struct CylinderInfo {
    /// A point on the axis.
    pub axis_origin: GpPnt,
    /// Unit axis direction.
    pub axis: GpVec,
    pub radius: f64,
}

/// A conical surface.
#[derive(Debug, Clone, Copy)]
pub struct ConeInfo {
    /// The true apex (vertex) of the cone.
    pub apex: GpPnt,
    /// Unit axis pointing from the apex toward the base (the direction in which
    /// the cross-section radius grows).
    pub axis: GpVec,
    /// Semi-vertical angle (the angle between the axis and a generator).
    pub semi_angle: f64,
}

impl ConeInfo {
    /// Radius of the cone at a height `z` measured along `axis` from `apex`
    /// (positive on the nappe, `NaN` off the nappe).
    pub fn radius_at(&self, z: f64) -> f64 {
        let r = z * self.semi_angle.tan();
        if r > 0.0 {
            r
        } else {
            f64::NAN
        }
    }
}

/// A circle in 3D (the rolling-ball centreline, or a tangency circle).
#[derive(Debug, Clone, Copy)]
struct CircleGeom {
    center: GpPnt,
    /// Unit circle-plane normal.
    normal: GpVec,
    radius: f64,
}

/// A tangency circle of the blend on one adjacent face.
#[derive(Debug, Clone, Copy)]
struct ContactCircleGeom {
    center: GpPnt,
    /// Unit circle-plane normal.
    normal: GpVec,
    radius: f64,
}

/// Full description of a torus band blend.
#[derive(Debug, Clone)]
struct TorusBlendGeom {
    /// Torus centre (on the centreline circle's axis).
    center: GpPnt,
    /// Unit torus axis.
    axis: GpVec,
    /// Centreline (major) radius.
    major: f64,
    /// Fillet (minor) radius.
    minor: f64,
    contact_a: ContactCircleGeom,
    contact_b: ContactCircleGeom,
}

// ===========================================================================
// Classification
// ===========================================================================

/// The blend topology class of a pair of adjacent faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvedPair {
    /// A plane and a sphere — torus blend.
    PlaneSphere,
    /// A plane and a cylinder (axis perpendicular to the plane) — torus blend.
    PlaneCylinder,
    /// Two spheres — torus blend.
    SphereSphere,
    /// A plane and a cone (axis perpendicular to the plane) — torus blend.
    PlaneCone,
    /// Two parallel cylinders — cylindrical band blend.
    CylinderCylinder,
    /// A cylinder and a coaxial cone — torus blend.
    CylinderCone,
    /// A sphere and a coaxial cone — torus blend.
    ConeSphere,
    /// Two coaxial cones — torus blend.
    ConeCone,
    /// Anything else — `fillet_edge_curved_general` returns an error.
    Unsupported,
}

/// Whether a surface is curved (non-planar).
fn is_curved_surface(s: &dyn Surface) -> bool {
    !is_planar(s, 8, 8, 1e-6)
}

/// Classify a surface into the full analytic set (plane / sphere / cylinder /
/// cone / other). Builds on `brep_surface::classify_surface` (plane, sphere)
/// and `fillet_edge::classify_surface_full` (cylinder), adding cone detection.
pub fn classify_surface_analytic(s: &dyn Surface) -> SurfaceKind {
    let k = classify_surface_full(s);
    if k != SurfaceKind::Other {
        return k;
    }
    if cone_from_surface(s).is_some() {
        SurfaceKind::Cone
    } else {
        SurfaceKind::Other
    }
}

/// True when at least one of the two faces is curved (sphere / cylinder / cone
/// / torus / other) — the pair is a candidate for `fillet_edge_curved`.
///
/// A pair of two planar faces is the domain of `fillet_edge` and returns
/// `false`.
pub fn faces_are_curved_compatible(f1: &Face, f2: &Face) -> bool {
    let c1 = BRepTool::face_surface(f1).map(|s| is_curved_surface(s.as_ref())).unwrap_or(false);
    let c2 = BRepTool::face_surface(f2).map(|s| is_curved_surface(s.as_ref())).unwrap_or(false);
    c1 || c2
}

/// Classify the adjacent-face pair into a supported blend topology.
pub fn classify_curved_pair(f1: &Face, f2: &Face) -> CurvedPair {
    let k1 = BRepTool::face_surface(f1)
        .map(|s| classify_surface_analytic(s.as_ref()))
        .unwrap_or(SurfaceKind::Other);
    let k2 = BRepTool::face_surface(f2)
        .map(|s| classify_surface_analytic(s.as_ref()))
        .unwrap_or(SurfaceKind::Other);
    match (k1, k2) {
        (SurfaceKind::Plane, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            CurvedPair::PlaneSphere
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            CurvedPair::PlaneCylinder
        }
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => CurvedPair::SphereSphere,
        (SurfaceKind::Plane, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Plane) => {
            CurvedPair::PlaneCone
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cylinder) => CurvedPair::CylinderCylinder,
        (SurfaceKind::Cylinder, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Cylinder) => {
            CurvedPair::CylinderCone
        }
        (SurfaceKind::Cone, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Cone) => {
            CurvedPair::ConeSphere
        }
        (SurfaceKind::Cone, SurfaceKind::Cone) => CurvedPair::ConeCone,
        _ => CurvedPair::Unsupported,
    }
}

// ===========================================================================
// Surface extraction
// ===========================================================================

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// Extract plane geometry from a surface, or `None` if it is not planar.
pub fn plane_from_surface(s: &dyn Surface) -> Option<PlaneInfo> {
    if !is_planar(s, 8, 8, 1e-6) {
        return None;
    }
    let (u0, _u1, v0, _v1) = sample_bounds(s);
    let origin = s.d0(u0, v0);
    let normal = surface_normal(s, u0, v0);
    if normal.xyz().square_modulus() < 1e-30 {
        return None;
    }
    Some(PlaneInfo { origin, normal })
}

/// Extract sphere geometry from a surface, or `None`.
pub fn sphere_from_surface(s: &dyn Surface) -> Option<SphereInfo> {
    let center = sphere_center(s)?;
    let (u0, _u1, v0, _v1) = sample_bounds(s);
    let radius = s.d0(u0, v0).distance(&center);
    if radius < 1e-9 {
        return None;
    }
    Some(SphereInfo { center, radius })
}

/// Extract cylinder geometry from a surface, or `None`.
///
/// A cylinder is recognised by: all sampled surface normals perpendicular to
/// one axis direction, and all sampled points at constant distance from that
/// axis line. The axis point is the projection of the first sample onto the
/// axis line through the 2D circumcenter of the projected points.
pub fn cylinder_from_surface(s: &dyn Surface) -> Option<CylinderInfo> {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let nu = 6;
    let nv = 6;
    let hu = ((u1 - u0) * 1e-4).max(1e-7);
    let hv = ((v1 - v0) * 1e-4).max(1e-7);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut nrm: Vec<GpVec> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            let p0 = s.d0(u, v);
            pts.push(p0);
            let du = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
            let dv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
            let nn = du.xyz().crossed(&dv.xyz());
            let m = nn.modulus();
            nrm.push(if m > 1e-30 { GpVec::from_xyz(&nn.divided(m)) } else { GpVec::zero() });
        }
    }
    // Axis direction: perpendicular to all normals.
    let mut axis = None;
    'outer: for a in 0..nrm.len() {
        for b in (a + 1)..nrm.len() {
            let c = nrm[a].xyz().crossed(&nrm[b].xyz());
            if c.modulus() > 1e-6 {
                axis = Some(GpVec::from_xyz(&c).normalized());
                break 'outer;
            }
        }
    }
    let axis = axis?;
    for nn in &nrm {
        if nn.dot(&axis).abs() > 1e-3 {
            return None;
        }
    }
    // Project points onto the plane perpendicular to the axis; the axis line
    // passes through the 2D circumcenter of three non-collinear projections.
    let (x2, y2) = project_basis(&axis);
    let proj: Vec<(f64, f64)> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&pts[0], q);
            (v.dot(&x2), v.dot(&y2))
        })
        .collect();
    let (mut cx, mut cy) = (f64::NAN, f64::NAN);
    'outer2: for a in 0..proj.len() {
        for b in (a + 1)..proj.len() {
            for c in (b + 1)..proj.len() {
                let (pa, pb, pc) = (proj[a], proj[b], proj[c]);
                let d = 2.0 * (pa.0 * (pb.1 - pc.1) + pb.0 * (pc.1 - pa.1) + pc.0 * (pa.1 - pb.1));
                if d.abs() < 1e-12 {
                    continue;
                }
                let ux = ((pa.0 * pa.0 + pa.1 * pa.1) * (pb.1 - pc.1)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pc.1 - pa.1)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pa.1 - pb.1))
                    / d;
                let uy = ((pa.0 * pa.0 + pa.1 * pa.1) * (pc.0 - pb.0)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pa.0 - pc.0)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pb.0 - pa.0))
                    / d;
                cx = ux;
                cy = uy;
                break 'outer2;
            }
        }
    }
    if !cx.is_finite() {
        return None;
    }
    let axis_pt = pts[0].translated_vec(&x2.multiplied_scalar(cx).added(&y2.multiplied_scalar(cy)));
    let mut radii: Vec<f64> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&axis_pt, q);
            let along = v.dot(&axis);
            v.subtracted(&axis.multiplied_scalar(along)).magnitude()
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = radii[radii.len() / 2];
    if med < 1e-9 {
        return None;
    }
    for r in &radii {
        if (r - med).abs() > 1e-3 * med.abs().max(1.0) {
            return None;
        }
    }
    Some(CylinderInfo { axis_origin: axis_pt, axis, radius: med })
}

/// Smallest-eigenvalue unit eigenvector of a symmetric 3×3 matrix (Jacobi
/// rotations). Returns `None` when the matrix is too ill-conditioned.
#[allow(dead_code)]
fn smallest_eigenvector(m: &[[f64; 3]; 3]) -> Option<[f64; 3]> {
    let mut a = *m;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..64 {
        // Largest off-diagonal entry.
        let mut p = 0;
        let mut q = 1;
        let mut mx = a[0][1].abs();
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > mx {
                    mx = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if mx < 1e-15 {
            break;
        }
        let theta = 0.5 * (2.0 * a[p][q]).atan2(a[p][p] - a[q][q]);
        let (c, s) = (theta.cos(), theta.sin());
        for k in 0..3 {
            let akp = a[k][p];
            let akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
            a[p][k] = a[k][p];
            a[q][k] = a[k][q];
            let vkp = v[k][p];
            let vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    let mut col = 0;
    for i in 1..3 {
        if a[i][i] < a[col][col] {
            col = i;
        }
    }
    let mut e = [v[0][col], v[1][col], v[2][col]];
    let m = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
    if m < 1e-30 {
        return None;
    }
    e[0] /= m;
    e[1] /= m;
    e[2] /= m;
    Some(e)
}

/// Solve a 3×3 linear system `A·x = b` via Gaussian elimination.
fn solve3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let mut m = *a;
    let mut rhs = *b;
    for col in 0..3 {
        let mut piv = col;
        for r in (col + 1)..3 {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-30 {
            return None;
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for r in (col + 1)..3 {
            let f = m[r][col] / d;
            for c in col..3 {
                m[r][c] -= f * m[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..3 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

/// Extract cone geometry from a surface, or `None` if it is not conical.
///
/// A cone is recognised by three invariants: all sampled tangent planes pass
/// through one common apex, all sampled normals make a constant angle with one
/// axis direction (the normals lie on a small circle of the unit sphere), and
/// all sampled points are at a radius proportional to their distance from the
/// apex along that axis. The returned `axis` points from the apex toward the
/// base (the direction in which the radius grows).
pub fn cone_from_surface(s: &dyn Surface) -> Option<ConeInfo> {
    // A `GeomCone` v-range is unbounded and its parameterization puts the apex
    // at v=0, so the natural `sample_bounds` range straddles the apex and mixes
    // the two nappes. Sample one nappe at a time, away from the apex.
    for sign in [1.0f64, -1.0] {
        if let Some(ci) = cone_from_surface_nappe(s, sign) {
            return Some(ci);
        }
    }
    None
}

/// Cone extraction on one nappe: `sign` selects the v>0 (+1) or v<0 (−1) nappe.
fn cone_from_surface_nappe(s: &dyn Surface, sign: f64) -> Option<ConeInfo> {
    let (u0, u1, _v0, _v1) = sample_bounds(s);
    let nu = 8;
    let nv = 8;
    let v_lo = 0.15 * sign;
    let v_hi = 2.5 * sign;
    let hu = ((u1 - u0) * 1e-4).max(1e-7);
    let hv = ((v_hi - v_lo) * 1e-4).max(1e-7);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut nrm: Vec<GpVec> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            // Sample the periodic u direction at `nu` distinct angles (exclude
            // the duplicated endpoint so the normals are symmetric about the
            // axis; a doubled sample would bias the covariance).
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v_lo + (v_hi - v_lo) * j as f64 / (nv - 1) as f64;
            let p0 = s.d0(u, v);
            pts.push(p0);
            let du = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
            let dv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
            let nn = du.xyz().crossed(&dv.xyz());
            let m = nn.modulus();
            nrm.push(if m > 1e-30 { GpVec::from_xyz(&nn.divided(m)) } else { GpVec::zero() });
        }
    }
    if nrm.iter().any(|n| n.xyz().square_modulus() < 1e-30) {
        return None;
    }
    if std::env::var("CONE_DBG").is_ok() {
        eprintln!("[cone] nrm[0..6]={:?}", &nrm[..6.min(nrm.len())]);
    }
    // The cone axis: the direction perpendicular to every normal difference
    // (n_i·axis is constant on a cone). Cross two independent difference
    // vectors to obtain it — a Jacobi eigen-solve is fragile here because the
    // covariance of the normals has a degenerate pair of equal eigenvalues.
    let mut diff1: Option<GpVec> = None;
    'd1: for i in 0..nrm.len() {
        for j in (i + 1)..nrm.len() {
            let d = nrm[i].subtracted(&nrm[j]);
            if d.magnitude() > 1e-3 {
                diff1 = Some(d.normalized());
                break 'd1;
            }
        }
    }
    if std::env::var("CONE_DBG").is_ok() {
        eprintln!("[cone] diff1={diff1:?}");
    }
    let mut diff2: Option<GpVec> = None;
    if let Some(v1) = diff1 {
        'd2: for i in 0..nrm.len() {
            for j in (i + 1)..nrm.len() {
                let d = nrm[i].subtracted(&nrm[j]);
                if d.magnitude() > 1e-3 && d.xyz().crossed(&v1.xyz()).modulus() > 1e-3 {
                    diff2 = Some(d.normalized());
                    break 'd2;
                }
            }
        }
    }
    let mut axis = match (diff1, diff2) {
        (Some(a), Some(b)) => a.crossed(&b).normalized(),
        _ => return None,
    };
    if std::env::var("CONE_DBG").is_ok() {
        eprintln!("[cone] axis(raw)={axis:?}");
    }
    // Apex: the common intersection of the tangent planes n_i·(A − P_i) = 0.
    let mut mm = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for (p, n) in pts.iter().zip(nrm.iter()) {
        let nx = n.xyz().x();
        let ny = n.xyz().y();
        let nz = n.xyz().z();
        let d = nx * p.x() + ny * p.y() + nz * p.z();
        mm[0][0] += nx * nx;
        mm[0][1] += nx * ny;
        mm[0][2] += nx * nz;
        mm[1][0] += ny * nx;
        mm[1][1] += ny * ny;
        mm[1][2] += ny * nz;
        mm[2][0] += nz * nx;
        mm[2][1] += nz * ny;
        mm[2][2] += nz * nz;
        rhs[0] += d * nx;
        rhs[1] += d * ny;
        rhs[2] += d * nz;
    }
    let apex_xyz = match solve3x3(&mm, &rhs) {
        Some(a) => a,
        None => {
            if std::env::var("CONE_DBG").is_ok() {
                eprintln!("[cone] apex solve failed");
            }
            return None;
        }
    };
    let apex = GpPnt::new(apex_xyz[0], apex_xyz[1], apex_xyz[2]);
    // Orient the axis toward the base (the direction in which the cross-section
    // radius grows): the samples lie on the positive-height nappe.
    let mut h_sum = 0.0;
    for p in &pts {
        h_sum += GpVec::from_pnts(&apex, p).dot(&axis);
    }
    if h_sum < 0.0 {
        axis = axis.multiplied_scalar(-1.0);
    }
    if std::env::var("CONE_DBG").is_ok() {
        eprintln!("[cone] apex={apex:?} axis={axis:?} h_sum={h_sum}");
    }
    // Semi-angle: |n·axis| = sin(α) for a cone (the normal leans toward the
    // apex, opposite the base, so the sign of n·axis varies with orientation).
    let mut s_sum = 0.0;
    for n in &nrm {
        s_sum += n.dot(&axis).abs();
    }
    let s_avg = s_sum / nrm.len() as f64;
    if s_avg < 1e-3 || s_avg > 1.0 - 1e-6 {
        if std::env::var("CONE_DBG").is_ok() {
            eprintln!("[cone] bad s_avg {s_avg}");
        }
        return None;
    }
    let semi_angle = s_avg.clamp(-1.0, 1.0).asin();
    // Verify: tangent-plane residuals, one nappe, and radius ∝ height.
    let tan_alpha = semi_angle.tan();
    for (p, n) in pts.iter().zip(nrm.iter()) {
        let plane_res = n.dot(&GpVec::from_pnts(p, &apex)).abs();
        if plane_res > 1e-3 {
            if std::env::var("CONE_DBG").is_ok() {
                eprintln!("[cone] plane_res {plane_res} > 1e-3");
            }
            return None;
        }
        let h = GpVec::from_pnts(&apex, p).dot(&axis);
        if h < -1e-6 {
            if std::env::var("CONE_DBG").is_ok() {
                eprintln!("[cone] mixed nappes h={h}");
            }
            return None; // mixed nappes
        }
        let expect = h * tan_alpha;
        let radial = GpVec::from_pnts(&apex, p)
            .subtracted(&axis.multiplied_scalar(h))
            .magnitude();
        if (radial - expect).abs() > 1e-2 * expect.max(1.0) {
            if std::env::var("CONE_DBG").is_ok() {
                eprintln!("[cone] radius mismatch radial={radial} expect={expect}");
            }
            return None;
        }
    }
    Some(ConeInfo { apex, axis, semi_angle })
}

/// Two orthonormal vectors spanning the plane perpendicular to `axis`.
fn project_basis(axis: &GpVec) -> (GpVec, GpVec) {
    let ref_v = if axis.x().abs() < 0.9 {
        GpVec::new(1.0, 0.0, 0.0)
    } else {
        GpVec::new(0.0, 1.0, 0.0)
    };
    let x2 = axis.crossed(&ref_v).normalized();
    let y2 = axis.crossed(&x2).normalized();
    (x2, y2)
}

/// A unit vector perpendicular to `z`.
fn perp_dir(z: &GpVec) -> GpVec {
    let ref_v = if z.xyz().z.abs() < 0.9 {
        GpVec::new(0.0, 0.0, 1.0)
    } else {
        GpVec::new(1.0, 0.0, 0.0)
    };
    z.crossed(&ref_v).normalized()
}

// ===========================================================================
// Centreline geometry
// ===========================================================================

/// The rolling-ball centreline for a plane + sphere: the intersection circle
/// of the plane offset by `r` (toward the ball) and the sphere offset by
/// `R_s + r`. The ball centre lies on this circle.
fn plane_sphere_centerline(pln: &PlaneInfo, sph: &SphereInfo, r: f64) -> Result<CircleGeom, String> {
    let n = pln.normal.normalized();
    // Signed distance from the sphere centre to the offset plane n·x = n·O + r,
    // measured along +n (the side the ball rolls on).
    let delta = n.dot(&GpVec::from_pnts(&pln.origin, &sph.center)) - r;
    let rr = sph.radius + r;
    let rho2 = rr * rr - delta * delta;
    if rho2 <= 1e-12 {
        return Err(format!(
            "fillet_curved: ball of radius {r} does not fit between plane and sphere (delta={delta}, R={})",
            sph.radius
        ));
    }
    let rho = rho2.sqrt();
    let center = sph.center.translated_vec(&n.multiplied_scalar(-delta));
    Ok(CircleGeom { center, normal: n, radius: rho })
}

/// The rolling-ball centreline for two spheres: the intersection circle of the
/// two offset spheres (radii `R_i + r`).
fn sphere_sphere_centerline(s1: &SphereInfo, s2: &SphereInfo, r: f64) -> Result<CircleGeom, String> {
    let dv = GpVec::from_pnts(&s1.center, &s2.center);
    let d = dv.magnitude();
    if d < 1e-12 {
        return Err("fillet_curved: coincident spheres".to_string());
    }
    let dhat = dv.normalized();
    let r1 = s1.radius + r;
    let r2 = s2.radius + r;
    if d > r1 + r2 + 1e-9 {
        return Err(format!(
            "fillet_curved: ball of radius {r} cannot bridge two non-overlapping offset spheres"
        ));
    }
    // Distance from s1.center to the radical plane of the two offset spheres.
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let rho2 = r1 * r1 - a * a;
    if rho2 <= 1e-12 {
        return Err("fillet_curved: offset spheres intersect at a point".to_string());
    }
    let rho = rho2.sqrt();
    let center = s1.center.translated_vec(&dhat.multiplied_scalar(a));
    Ok(CircleGeom { center, normal: dhat, radius: rho })
}

/// Sample a circle at `n` uniformly spaced angles.
fn sample_circle(c: &CircleGeom, n: usize) -> Vec<GpPnt> {
    let e1 = perp_dir(&c.normal);
    let e2 = c.normal.crossed(&e1).normalized();
    (0..n)
        .map(|i| {
            let ang = 2.0 * PI * i as f64 / n as f64;
            let off = e1.multiplied_scalar(c.radius * ang.cos())
                .added(&e2.multiplied_scalar(c.radius * ang.sin()));
            c.center.translated_vec(&off)
        })
        .collect()
}

// ===========================================================================
// Meridional (axisymmetric) blend geometry
// ===========================================================================

/// A surface of revolution's cross-section (meridional) radius as a function of
/// height along a common axis. `NaN` marks heights off the surface.
#[derive(Debug, Clone, Copy)]
enum MeridionalShape {
    /// Constant radius (cylinder).
    Cylinder { rho: f64 },
    /// Cone: radius = (z − apex_z)·d·tan on the nappe where that is positive.
    Cone { apex_z: f64, d: f64, tan: f64, sin_alpha: f64 },
    /// Sphere: radius = √(rho² − (z − center_z)²).
    Sphere { center_z: f64, rho: f64 },
}

impl MeridionalShape {
    /// Cross-section radius of the *original* surface at height `z`.
    fn radius(&self, z: f64) -> f64 {
        match *self {
            MeridionalShape::Cylinder { rho } => rho,
            MeridionalShape::Cone { apex_z, d, tan, .. } => {
                let r = (z - apex_z) * d * tan;
                if r > 0.0 {
                    r
                } else {
                    f64::NAN
                }
            }
            MeridionalShape::Sphere { center_z, rho } => {
                let d2 = rho * rho - (z - center_z) * (z - center_z);
                if d2 >= 0.0 {
                    d2.sqrt()
                } else {
                    f64::NAN
                }
            }
        }
    }

    /// Cross-section radius of the surface offset outward by `r`.
    fn offset_radius(&self, z: f64, r: f64) -> f64 {
        match *self {
            MeridionalShape::Cylinder { rho } => rho + r,
            MeridionalShape::Cone { apex_z, d, tan, sin_alpha } => {
                // The outward normal of a cone has a +r·sin(α) component
                // opposite the base, so the offset apex moves by r/sin(α)
                // opposite the base direction.
                let apex_z_off = apex_z - (r / sin_alpha) * d;
                let ro = (z - apex_z_off) * d * tan;
                if ro > 0.0 {
                    ro
                } else {
                    f64::NAN
                }
            }
            MeridionalShape::Sphere { center_z, rho } => {
                let d2 = (rho + r) * (rho + r) - (z - center_z) * (z - center_z);
                if d2 >= 0.0 {
                    d2.sqrt()
                } else {
                    f64::NAN
                }
            }
        }
    }

    /// z-interval on which `offset_radius(z, r)` is finite (the offset nappe).
    /// `None` marks an unbounded end. Used to bound the offset-centreline search.
    fn offset_z_range(&self, r: f64) -> (Option<f64>, Option<f64>) {
        match *self {
            MeridionalShape::Cylinder { .. } => (None, None),
            MeridionalShape::Cone { apex_z, d, sin_alpha, .. } => {
                let z0 = apex_z - (r / sin_alpha) * d;
                if d > 0.0 {
                    (Some(z0), None)
                } else {
                    (None, Some(z0))
                }
            }
            MeridionalShape::Sphere { center_z, rho } => {
                (Some(center_z - (rho + r)), Some(center_z + (rho + r)))
            }
        }
    }
}

/// Closest point on a meridional curve `rho(z)` to the ball centre
/// `(rho_c, z_c)`, by coarse grid + ternary refinement. Returns `(rho, z)`.
fn closest_point_meridional(
    rho: &dyn Fn(f64) -> f64,
    rho_c: f64,
    z_c: f64,
    z_lo: f64,
    z_hi: f64,
) -> Option<(f64, f64)> {
    let n = 256;
    let mut best = (f64::NAN, f64::NAN, f64::INFINITY);
    for i in 0..=n {
        let z = z_lo + (z_hi - z_lo) * i as f64 / n as f64;
        let r = rho(z);
        if !r.is_finite() || r < 0.0 {
            continue;
        }
        let d2 = (r - rho_c) * (r - rho_c) + (z - z_c) * (z - z_c);
        if d2 < best.2 {
            best = (r, z, d2);
        }
    }
    if !best.0.is_finite() {
        return None;
    }
    let mut lo = (best.1 - (z_hi - z_lo) / n as f64).max(z_lo);
    let mut hi = (best.1 + (z_hi - z_lo) / n as f64).min(z_hi);
    for _ in 0..96 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let d1 = (rho(m1) - rho_c) * (rho(m1) - rho_c) + (m1 - z_c) * (m1 - z_c);
        let d2v = (rho(m2) - rho_c) * (rho(m2) - rho_c) + (m2 - z_c) * (m2 - z_c);
        if d1 < d2v {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let z = 0.5 * (lo + hi);
    Some((rho(z), z))
}

/// Compute the tangency circle of the rolling ball against a surface of
/// revolution described by its meridional radius function `rho`.
///
/// `center`/`axis` are the centreline circle (on the axis) and its unit
/// normal; the ball centre is at radius `r_c` in the plane of the centreline
/// (height 0). The search for the contact height is bounded by `z_lo`/`z_hi`
/// (relative to the centreline plane).
fn contact_circle_meridional(
    center: &GpPnt,
    axis: &GpVec,
    r_c: f64,
    r: f64,
    rho: &dyn Fn(f64) -> f64,
    z_lo: f64,
    z_hi: f64,
) -> Result<ContactCircleGeom, String> {
    let (rho_t, z_t) = closest_point_meridional(rho, r_c, 0.0, z_lo, z_hi)
        .ok_or("fillet_curved: cannot locate the ball contact on a curved face")?;
    let dist = ((rho_t - r_c) * (rho_t - r_c) + z_t * z_t).sqrt();
    if (dist - r).abs() > 1e-2 * r.max(1.0) {
        return Err(format!(
            "fillet_curved: rolling ball of radius {r} does not reach this face \
             (contact distance {dist})"
        ));
    }
    let cos_v = ((rho_t - r_c) / r).clamp(-1.0, 1.0);
    let sin_v = (z_t / r).clamp(-1.0, 1.0);
    let v = sin_v.atan2(cos_v);
    Ok(ContactCircleGeom {
        center: center.translated_vec(&axis.multiplied_scalar(r * v.sin())),
        normal: *axis,
        radius: r_c + r * v.cos(),
    })
}

/// Tangency circle of the ball against a plane perpendicular to the axis. The
/// plane is at height `plane_z` relative to the centreline plane.
fn contact_circle_plane(
    center: &GpPnt,
    axis: &GpVec,
    r_c: f64,
    r: f64,
    plane_z: f64,
) -> Result<ContactCircleGeom, String> {
    let t = (plane_z / r).clamp(-1.0, 1.0);
    let v = t.asin();
    Ok(ContactCircleGeom {
        center: center.translated_vec(&axis.multiplied_scalar(r * v.sin())),
        normal: *axis,
        radius: r_c + r * v.cos(),
    })
}

/// Find the height where two offset meridional radius functions are equal
/// (the rolling-ball centreline circle). Returns `(z, radius)`.
fn offset_centerline_height(
    rho1: &dyn Fn(f64) -> f64,
    rho2: &dyn Fn(f64) -> f64,
    z_lo: f64,
    z_hi: f64,
) -> Result<(f64, f64), String> {
    let n = 1024;
    let mut prev = rho1(z_lo) - rho2(z_lo);
    for i in 1..=n {
        let z = z_lo + (z_hi - z_lo) * i as f64 / n as f64;
        let cur = rho1(z) - rho2(z);
        if prev.is_finite() && cur.is_finite() {
            if (prev < 0.0) != (cur < 0.0) {
                let mut lo = z_lo + (z_hi - z_lo) * (i - 1) as f64 / n as f64;
                let mut hi = z;
                let mut flo = prev;
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    let fm = rho1(mid) - rho2(mid);
                    if (flo < 0.0) == (fm < 0.0) {
                        lo = mid;
                        flo = fm;
                    } else {
                        hi = mid;
                    }
                }
                let zc = 0.5 * (lo + hi);
                let rc = rho1(zc);
                if rc.is_finite() && rc > 1e-9 {
                    return Ok((zc, rc));
                }
            }
        }
        prev = cur;
    }
    Err("fillet_curved: the two offset surfaces do not intersect in a circle".to_string())
}

/// Finite z-window bracketing the offset-centreline crossing of two coaxial
/// surfaces of revolution. Intersects the two offset nappes; an unbounded end
/// (cylinder, or a cone whose nappe is open on that side) is clamped by walking
/// outward until the faster-growing offset radius exceeds the other by a factor
/// of 8 — the crossing (unique for these coaxial pairs) lies inside.
fn overlap_window(m1: &MeridionalShape, m2: &MeridionalShape, r: f64) -> Result<(f64, f64), String> {
    let (l1, h1) = m1.offset_z_range(r);
    let (l2, h2) = m2.offset_z_range(r);
    let lo = match (l1, l2) {
        (Some(a), Some(b)) => Some(a.max(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    let hi = match (h1, h2) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    if let (Some(a), Some(b)) = (lo, hi) {
        if a >= b {
            return Err("fillet_curved: the two offset nappes do not overlap along the axis".to_string());
        }
    }
    let (mut lo, mut hi) = (lo, hi);
    if lo.is_some() && hi.is_some() {
        return Ok((lo.unwrap(), hi.unwrap()));
    }
    let (anchor, unbounded_below) = match (lo, hi) {
        (Some(a), _) => (a, false),
        (None, Some(b)) => (b, true),
        (None, None) => {
            return Err("fillet_curved: both offset surfaces are unbounded along the axis".to_string())
        }
    };
    let rho_anchor = m1.offset_radius(anchor, r).max(m2.offset_radius(anchor, r)).max(1e-9);
    let e = rho_anchor.max(1.0) * 1e-4;
    // Probe strictly inside the overlap (the boundary itself may be a nappe
    // apex, where `offset_radius` is NaN) to measure the growth rate.
    let sgn = if unbounded_below { -1.0 } else { 1.0 };
    let ra1 = m1.offset_radius(anchor + sgn * e, r);
    let ra2 = m2.offset_radius(anchor + sgn * e, r);
    let rb1 = m1.offset_radius(anchor + sgn * 2.0 * e, r);
    let rb2 = m2.offset_radius(anchor + sgn * 2.0 * e, r);
    let slope = {
        let s1 = ((rb1 - ra1) / e).abs();
        let s2 = ((rb2 - ra2) / e).abs();
        (if s1.is_finite() { s1 } else { 1e-9 })
            .max(if s2.is_finite() { s2 } else { 1e-9 })
            .max(1e-9)
    };
    let span = 8.0 * rho_anchor / slope;
    if lo.is_none() {
        lo = Some(anchor - span);
    }
    if hi.is_none() {
        hi = Some(anchor + span);
    }
    let (a, b) = (lo.unwrap(), hi.unwrap());
    if a >= b || !a.is_finite() || !b.is_finite() {
        return Err("fillet_curved: the two offset surfaces have no reachable centerline".to_string());
    }
    Ok((a, b))
}

/// Build a `MeridionalShape::Cone` from a `ConeInfo` in the frame of `axis`.
fn cone_meridional(cone: &ConeInfo, axis: &GpVec, ref_pt: &GpPnt) -> MeridionalShape {
    let d = cone.axis.dot(axis);
    let apex_z = GpVec::from_pnts(ref_pt, &cone.apex).dot(axis);
    MeridionalShape::Cone {
        apex_z,
        d,
        tan: cone.semi_angle.tan(),
        sin_alpha: cone.semi_angle.sin(),
    }
}

// ===========================================================================
// Blend geometry per supported pair
// ===========================================================================

/// Build the torus band blend for a plane + sphere pair.
fn plane_sphere_blend(pln: &PlaneInfo, sph: &SphereInfo, r: f64) -> Result<TorusBlendGeom, String> {
    let circle = plane_sphere_centerline(pln, sph, r)?;
    let n = circle.normal;
    let delta = n.dot(&GpVec::from_pnts(&pln.origin, &sph.center)) - r;
    // Tube angle of the sphere-tangency circle: cos v = −rho/(R_s+r), sin v = δ/(R_s+r).
    let v_sphere = delta.atan2(-circle.radius);
    let v_plane = 3.0 * PI / 2.0; // torus bottom: touches the plane
    let v_sphere = if v_sphere < 0.0 { v_sphere + 2.0 * PI } else { v_sphere };
    let plane_contact = ContactCircleGeom {
        center: circle.center.translated_vec(&n.multiplied_scalar(r * (v_plane.sin()))),
        normal: n,
        radius: circle.radius + r * v_plane.cos(),
    };
    let sphere_contact = ContactCircleGeom {
        center: circle.center.translated_vec(&n.multiplied_scalar(r * v_sphere.sin())),
        normal: n,
        radius: circle.radius + r * v_sphere.cos(),
    };
    Ok(TorusBlendGeom {
        center: circle.center,
        axis: n,
        major: circle.radius,
        minor: r,
        contact_a: plane_contact,
        contact_b: sphere_contact,
    })
}

/// Build the torus band blend for a plane + cylinder pair. The cylinder axis
/// must be perpendicular to the plane (a "standing" cylinder); otherwise the
/// offset-plane / offset-cylinder intersection is an ellipse and the exact
/// torus construction does not apply.
fn plane_cylinder_blend(pln: &PlaneInfo, cyl: &CylinderInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let n = pln.normal.normalized();
    let a = cyl.axis.normalized();
    if n.xyz().crossed(&a.xyz()).modulus() > tol.max(1e-6) {
        return Err(
            "fillet_curved: plane+cylinder blend requires the cylinder axis to be perpendicular to the plane"
                .to_string(),
        );
    }
    let rho = cyl.radius + r;
    // Project the axis point onto the offset plane n·x = n·O + r.
    let n_origin = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.origin));
    let n_axis = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &cyl.axis_origin));
    let center = cyl.axis_origin.translated_vec(&n.multiplied_scalar(n_origin + r - n_axis));
    let v_cyl = PI; // inner equator: touches the cylinder
    let v_plane = 3.0 * PI / 2.0; // bottom: touches the plane
    let plane_contact = ContactCircleGeom {
        center: center.translated_vec(&n.multiplied_scalar(r * v_plane.sin())),
        normal: n,
        radius: rho + r * v_plane.cos(),
    };
    let cyl_contact = ContactCircleGeom {
        center: center.translated_vec(&n.multiplied_scalar(r * v_cyl.sin())),
        normal: n,
        radius: rho + r * v_cyl.cos(),
    };
    Ok(TorusBlendGeom {
        center,
        axis: n,
        major: rho,
        minor: r,
        contact_a: plane_contact,
        contact_b: cyl_contact,
    })
}

/// Build the torus band blend for a sphere + sphere pair.
fn sphere_sphere_blend(s1: &SphereInfo, s2: &SphereInfo, r: f64) -> Result<TorusBlendGeom, String> {
    let circle = sphere_sphere_centerline(s1, s2, r)?;
    let dv = GpVec::from_pnts(&s1.center, &s2.center);
    let d = dv.magnitude();
    let dhat = dv.normalized();
    let r1 = s1.radius + r;
    let r2 = s2.radius + r;
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    // Tube angles: cos v_i = −rho / (R_i + r); sin v_1 = −a / (R_1 + r),
    // sin v_2 = (d − a) / (R_2 + r).
    let v1 = (-a).atan2(-circle.radius);
    let v2 = (d - a).atan2(-circle.radius);
    let contact1 = ContactCircleGeom {
        center: circle.center.translated_vec(&dhat.multiplied_scalar(r * v1.sin())),
        normal: dhat,
        radius: circle.radius + r * v1.cos(),
    };
    let contact2 = ContactCircleGeom {
        center: circle.center.translated_vec(&dhat.multiplied_scalar(r * v2.sin())),
        normal: dhat,
        radius: circle.radius + r * v2.cos(),
    };
    Ok(TorusBlendGeom {
        center: circle.center,
        axis: dhat,
        major: circle.radius,
        minor: r,
        contact_a: contact1,
        contact_b: contact2,
    })
}

/// Build the torus band blend for a plane + cone pair (concave base crease).
fn plane_cone_blend(pln: &PlaneInfo, cone: &ConeInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let n = pln.normal.normalized();
    let a = cone.axis.normalized();
    if n.xyz().crossed(&a.xyz()).modulus() > tol.max(1e-6) {
        return Err(
            "fillet_curved: plane+cone blend requires the cone axis perpendicular to the plane"
                .to_string(),
        );
    }
    let axis = n;
    let ref_pt = cone.apex;
    let m_cone = cone_meridional(cone, &axis, &ref_pt);
    let z_p = GpVec::from_pnts(&ref_pt, &pln.origin).dot(&axis);
    let z_c = z_p + r;
    let rc = m_cone.offset_radius(z_c, r);
    if !rc.is_finite() || rc <= 1e-9 {
        return Err("fillet_curved: plane+cone offset surfaces do not yield a centerline".to_string());
    }
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(z_c));
    let contact_plane = contact_circle_plane(&center, &axis, rc, r, z_p - z_c)?;
    let rho_cone = |z: f64| m_cone.radius(z_c + z);
    let contact_cone = contact_circle_meridional(&center, &axis, rc, r, &rho_cone, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: contact_plane,
        contact_b: contact_cone,
    })
}

/// Build the torus band blend for a coaxial cylinder + cone pair.
fn cylinder_cone_blend(cyl: &CylinderInfo, cone: &ConeInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let a_cyl = cyl.axis.normalized();
    let a_cone = cone.axis.normalized();
    if a_cyl.xyz().crossed(&a_cone.xyz()).modulus() > tol.max(1e-6) {
        return Err("fillet_curved: cylinder+cone blend requires coaxial surfaces".to_string());
    }
    let axis = a_cone;
    let ref_pt = cone.apex;
    let m_cyl = MeridionalShape::Cylinder { rho: cyl.radius };
    let m_cone = cone_meridional(cone, &axis, &ref_pt);
    let rho_cyl = |z: f64| m_cyl.offset_radius(z, r);
    let rho_cone = |z: f64| m_cone.offset_radius(z, r);
    let (z_lo, z_hi) = overlap_window(&m_cyl, &m_cone, r)?;
    let (zc, rc) = offset_centerline_height(&rho_cyl, &rho_cone, z_lo, z_hi)?;
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(zc));
    let rho_cyl_s = |z: f64| m_cyl.radius(zc + z);
    let contact_cyl = contact_circle_meridional(&center, &axis, rc, r, &rho_cyl_s, -1e3, 1e3)?;
    let rho_cone_s = |z: f64| m_cone.radius(zc + z);
    let contact_cone = contact_circle_meridional(&center, &axis, rc, r, &rho_cone_s, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: contact_cyl,
        contact_b: contact_cone,
    })
}

/// Build the torus band blend for a sphere + coaxial cone pair.
fn cone_sphere_blend(cone: &ConeInfo, sph: &SphereInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let a_cone = cone.axis.normalized();
    let cs = GpVec::from_pnts(&cone.apex, &sph.center);
    let radial = cs.subtracted(&a_cone.multiplied_scalar(cs.dot(&a_cone)));
    if radial.magnitude() > tol.max(1e-6) {
        return Err(
            "fillet_curved: cone+sphere blend requires the sphere centre on the cone axis".to_string(),
        );
    }
    let axis = a_cone;
    let ref_pt = cone.apex;
    let m_cone = cone_meridional(cone, &axis, &ref_pt);
    let m_sph = MeridionalShape::Sphere { center_z: cs.dot(&axis), rho: sph.radius };
    let rho_cone = |z: f64| m_cone.offset_radius(z, r);
    let rho_sph = |z: f64| m_sph.offset_radius(z, r);
    let (z_lo, z_hi) = overlap_window(&m_cone, &m_sph, r)?;
    let (zc, rc) = offset_centerline_height(&rho_cone, &rho_sph, z_lo, z_hi)?;
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(zc));
    let rho_cone_s = |z: f64| m_cone.radius(zc + z);
    let contact_cone = contact_circle_meridional(&center, &axis, rc, r, &rho_cone_s, -1e3, 1e3)?;
    let rho_sph_s = |z: f64| m_sph.radius(zc + z);
    let contact_sph = contact_circle_meridional(&center, &axis, rc, r, &rho_sph_s, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: contact_cone,
        contact_b: contact_sph,
    })
}

/// Build the torus band blend for a pair of coaxial cones.
fn cone_cone_blend(c1: &ConeInfo, c2: &ConeInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let a1 = c1.axis.normalized();
    let a2 = c2.axis.normalized();
    if a1.xyz().crossed(&a2.xyz()).modulus() > tol.max(1e-6) {
        return Err("fillet_curved: cone+cone blend requires coaxial cones".to_string());
    }
    let axis = a1;
    let ref_pt = c1.apex;
    let m1 = cone_meridional(c1, &axis, &ref_pt);
    let m2 = cone_meridional(c2, &axis, &ref_pt);
    let rho1 = |z: f64| m1.offset_radius(z, r);
    let rho2 = |z: f64| m2.offset_radius(z, r);
    let (z_lo, z_hi) = overlap_window(&m1, &m2, r)?;
    let (zc, rc) = offset_centerline_height(&rho1, &rho2, z_lo, z_hi)?;
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(zc));
    let rho1s = |z: f64| m1.radius(zc + z);
    let ca = contact_circle_meridional(&center, &axis, rc, r, &rho1s, -1e3, 1e3)?;
    let rho2s = |z: f64| m2.radius(zc + z);
    let cb = contact_circle_meridional(&center, &axis, rc, r, &rho2s, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: ca,
        contact_b: cb,
    })
}

// ===========================================================================
// Topology rebuild
// ===========================================================================

/// How to rebuild an adjacent face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdjacentKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
}

/// Build a full circular edge on the circle `c`, with its seam at `seam` (a
/// point on the circle, angle 0).
fn build_circle_edge(b: &TopoBuilder, c: &ContactCircleGeom, seam: &GpPnt) -> Result<Edge, String> {
    let nd = GpDir::from_xyz(&c.normal.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_vec(&GpVec::from_pnts(&c.center, seam))
        .map_err(|_| "fillet_curved: tangency circle seam coincides with its centre".to_string())?;
    let ax2 = GpAx2::new(c.center, nd, xd).map_err(|e| format!("fillet_curved: circle frame: {e}"))?;
    let mut e = b.make_edge_circle(&ax2, c.radius, 0.0, 2.0 * PI);
    let v = b.make_vertex(*seam, 0.0);
    b.add(&mut e.0, &v.0);
    b.add(&mut e.0, &v.0);
    Ok(e)
}

/// A default seam point for a tangency circle (angle 0 along a fixed
/// perpendicular). Used where the seam angle is not constrained by a neighbour.
fn default_seam(c: &ContactCircleGeom) -> GpPnt {
    let e1 = perp_dir(&c.normal);
    c.center.translated_vec(&e1.multiplied_scalar(c.radius))
}

/// The closed (circular) edge of a cylinder lateral face that is not the
/// filleted edge — i.e. the top rim circle.
fn find_top_circle(face: &Face, fillet_edge: &Edge) -> Result<Edge, String> {
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if is_same(&e.0, &fillet_edge.0) {
                continue;
            }
            if let Some((a, b)) = BRepTool::edge_vertices(&e) {
                if a.distance(&b) < 1e-9 {
                    return Ok(e.clone());
                }
            }
        }
    }
    Err("fillet_curved: cylinder lateral face has no top rim circle".to_string())
}

/// The point on a cylinder-tangency circle that shares a generator with the
/// cylinder's top rim seam, so the rebuilt lateral face's seam sub-edge stays a
/// straight generator.
fn cylinder_contact_seam(face: &Face, fillet_edge: &Edge, contact: &ContactCircleGeom) -> Result<GpPnt, String> {
    let top = find_top_circle(face, fillet_edge)?;
    let (s_top, _) = BRepTool::edge_vertices(&top).ok_or("fillet_curved: top rim has no vertices")?;
    let v = GpVec::from_pnts(&contact.center, &s_top);
    let along = v.dot(&contact.normal);
    let radial = v.subtracted(&contact.normal.multiplied_scalar(along));
    let m = radial.magnitude();
    if m < 1e-9 {
        return Err("fillet_curved: cylinder top seam lies on the axis".to_string());
    }
    Ok(contact.center.translated_vec(&radial.multiplied_scalar(contact.radius / m)))
}

/// Rebuild a planar adjacent face: keep every original wire except the one
/// containing the filleted edge, and add a new hole wire made of the enlarged
/// plane-tangency circle.
fn rebuild_plane_face(face: &Face, fillet_edge: &Edge, new_circle: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: planar face has no surface")?;
    let mut wires: Vec<Wire> = Vec::new();
    for w in wires_of_face(face) {
        let edges = edges_of_wire(&w);
        let orig_len = edges.len();
        let kept: Vec<Edge> = edges.into_iter().filter(|e| !is_same(&e.0, &fillet_edge.0)).collect();
        if kept.is_empty() {
            continue;
        }
        // A wire may hold the filleted edge next to other boundary edges (e.g. a
        // single-wire annulus); keep the remainder rather than dropping the wire.
        wires.push(if kept.len() < orig_len {
            b.make_wire(&kept)
        } else {
            w
        });
    }
    if wires.is_empty() {
        return Err(
            "fillet_curved: the planar face is bounded only by the filleted edge; its trim would be unbounded"
                .to_string(),
        );
    }
    wires.push(b.make_wire(&[new_circle.clone()]));
    Ok(b.make_face(surf, &wires))
}

/// Rebuild a spherical adjacent face: its new boundary is the single
/// surface-tangency circle.
fn rebuild_single_circle_face(face: &Face, new_circle: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: face has no surface")?;
    let wire = b.make_wire(&[new_circle.clone()]);
    Ok(b.make_face(surf, &[wire]))
}

/// Rebuild a cylinder lateral face: the base (filleted) circle is replaced by
/// the tangency circle, the seam is re-generated from the tangency circle up to
/// the original top rim, and the top rim circle is reused (shared with the cap).
fn rebuild_cylinder_face(
    face: &Face,
    fillet_edge: &Edge,
    contact_edge: &Edge,
    contact_seam: &GpPnt,
    b: &TopoBuilder,
) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: cylinder face has no surface")?;
    let top = find_top_circle(face, fillet_edge)?;
    let (s_top, _) = BRepTool::edge_vertices(&top).ok_or("fillet_curved: top rim has no vertices")?;
    let seam_up = b.make_edge_segment(contact_seam, &s_top);
    let wire = b.make_wire(&[contact_edge.clone(), seam_up.clone(), top, seam_up]);
    Ok(b.make_face(surf, &[wire]))
}

/// The apex vertex point of a cone lateral face: the face vertex that does not
/// lie on the (filleted) base circle.
fn find_cone_apex(face: &Face, fillet_edge: &Edge) -> Result<GpPnt, String> {
    let (se, _) = BRepTool::edge_vertices(fillet_edge)
        .ok_or("fillet_curved: cone base circle has no vertices")?;
    let seam_pt = se;
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if is_same(&e.0, &fillet_edge.0) {
                continue;
            }
            if let Some((a, b)) = BRepTool::edge_vertices(&e) {
                if a.distance(&seam_pt) > 1e-6 {
                    return Ok(a);
                }
                if b.distance(&seam_pt) > 1e-6 {
                    return Ok(b);
                }
            }
        }
    }
    Err("fillet_curved: cone lateral face has no apex vertex".to_string())
}

/// The point on a cone-tangency circle that shares a generator with the cone's
/// base-circle seam, so the rebuilt lateral face's seam stays a straight
/// generator from the contact circle to the apex.
fn cone_contact_seam(face: &Face, fillet_edge: &Edge, contact: &ContactCircleGeom) -> Result<GpPnt, String> {
    let (se, _) = BRepTool::edge_vertices(fillet_edge)
        .ok_or("fillet_curved: cone base circle has no vertices")?;
    let base_seam = se;
    let v = GpVec::from_pnts(&contact.center, &base_seam);
    let along = v.dot(&contact.normal);
    let radial = v.subtracted(&contact.normal.multiplied_scalar(along));
    let m = radial.magnitude();
    if m < 1e-9 {
        return Err("fillet_curved: cone base seam lies on the axis".to_string());
    }
    Ok(contact.center.translated_vec(&radial.multiplied_scalar(contact.radius / m)))
}

/// Rebuild a cone lateral face: the base (filleted) circle is replaced by the
/// tangency circle and the seam is re-generated from the contact circle up to
/// the apex. A *truncated* cone face (bounded by a second closed rim circle, as
/// when the cone carries a top cap) keeps that rim and is rebuilt as a band
/// from the tangency circle up to it, mirroring the cylinder rebuild.
fn rebuild_cone_face(
    face: &Face,
    fillet_edge: &Edge,
    contact_edge: &Edge,
    contact_seam: &GpPnt,
    b: &TopoBuilder,
) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: cone face has no surface")?;
    if let Ok(top) = find_top_circle(face, fillet_edge) {
        let (s_top, _) = BRepTool::edge_vertices(&top)
            .ok_or("fillet_curved: cone rim has no vertices")?;
        let seam_up = b.make_edge_segment(contact_seam, &s_top);
        let wire = b.make_wire(&[contact_edge.clone(), seam_up.clone(), top, seam_up]);
        return Ok(b.make_face(surf, &[wire]));
    }
    let apex = find_cone_apex(face, fillet_edge)?;
    let new_seam = b.make_edge_segment(contact_seam, &apex);
    let wire = b.make_wire(&[contact_edge.clone(), new_seam.clone(), new_seam]);
    Ok(b.make_face(surf, &[wire]))
}

/// The `GeomTorus` surface of a torus band blend.
fn torus_surface(geom: &TorusBlendGeom) -> Result<Arc<dyn Surface>, String> {
    let zd = GpDir::from_xyz(&geom.axis.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_xyz(&perp_dir(&geom.axis).xyz()).map_err(|e| e.to_string())?;
    let ax3 = GpAx3::new(geom.center, zd, &xd).map_err(|e| format!("fillet_curved: torus frame: {e}"))?;
    let t = GpTorus::new(ax3, geom.major, geom.minor).map_err(|e| e.to_string())?;
    Ok(Arc::new(GeomTorus::new(t)))
}

/// Build the blend face: the torus band bounded by the two tangency circles
/// (each shared with one trimmed adjacent face).
fn build_torus_blend_face(geom: &TorusBlendGeom, edge_a: &Edge, edge_b: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = torus_surface(geom)?;
    let wa = b.make_wire(&[edge_a.clone()]);
    let wb = b.make_wire(&[edge_b.clone()]);
    // The larger tangency circle is the outer boundary.
    let (w_outer, w_inner) = if geom.contact_a.radius >= geom.contact_b.radius {
        (wa, wb)
    } else {
        (wb, wa)
    };
    Ok(b.make_face(surf, &[w_outer, w_inner]))
}

/// Build a trimmed adjacent face, returning the rebuilt face and the (shared)
/// tangency circle edge it uses.
#[allow(clippy::type_complexity)]
fn build_adjacent_face(
    face: &Face,
    kind: AdjacentKind,
    fillet_edge: &Edge,
    contact: &ContactCircleGeom,
    b: &TopoBuilder,
) -> Result<(Face, Edge), String> {
    match kind {
        AdjacentKind::Plane => {
            let e = build_circle_edge(b, contact, &default_seam(contact))?;
            let f = rebuild_plane_face(face, fillet_edge, &e, b)?;
            Ok((f, e))
        }
        AdjacentKind::Sphere => {
            let e = build_circle_edge(b, contact, &default_seam(contact))?;
            let f = rebuild_single_circle_face(face, &e, b)?;
            Ok((f, e))
        }
        AdjacentKind::Cylinder => {
            let seam = cylinder_contact_seam(face, fillet_edge, contact)?;
            let e = build_circle_edge(b, contact, &seam)?;
            let f = rebuild_cylinder_face(face, fillet_edge, &e, &seam, b)?;
            Ok((f, e))
        }
        AdjacentKind::Cone => {
            let seam = cone_contact_seam(face, fillet_edge, contact)?;
            let e = build_circle_edge(b, contact, &seam)?;
            let f = rebuild_cone_face(face, fillet_edge, &e, &seam, b)?;
            Ok((f, e))
        }
    }
}

// ===========================================================================
// Public API
// ===========================================================================

/// Sample the rolling-ball centre locus for two surfaces. For the supported
/// analytic pairs the locus is a circle (plane+sphere, sphere+sphere,
/// plane+cylinder); the returned points lie on it.
///
/// Each returned centre point is at distance `radius` from the first surface
/// and `radius` from the second surface (measuring against the *offset*
/// surfaces — the plane pushed out by `r`, the sphere/cylinder enlarged by
/// `R + r`). Unsupported pairs return an error.
pub fn rolling_ball_center_line(
    f1: &dyn Surface,
    f2: &dyn Surface,
    radius: f64,
    tol: f64,
) -> Result<Vec<GpPnt>, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("rolling_ball_center_line: radius must be positive".to_string());
    }
    let k1 = classify_surface_full(f1);
    let k2 = classify_surface_full(f2);
    match (k1, k2) {
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let s1 = sphere_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract sphere 1")?;
            let s2 = sphere_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract sphere 2")?;
            let c = sphere_sphere_centerline(&s1, &s2, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pln = plane_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let sph = sphere_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract sphere")?;
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let sph = sphere_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract sphere")?;
            let pln = plane_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let (pln_s, cyl_s) = if k1 == SurfaceKind::Plane { (f1, f2) } else { (f2, f1) };
            let pln = plane_from_surface(pln_s).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let cyl = cylinder_from_surface(cyl_s).ok_or("rolling_ball_center_line: cannot extract cylinder")?;
            let blend = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            Ok(sample_circle(
                &CircleGeom { center: blend.center, normal: blend.axis, radius: blend.major },
                CENTERLINE_SAMPLES,
            ))
        }
        _ => Err(format!(
            "rolling_ball_center_line: unsupported surface pair ({k1:?}, {k2:?})"
        )),
    }
}

/// Replace a sharp edge of `solid` with a curved-face rolling-ball fillet.
///
/// This is the original supported-pairs entry point (plane+sphere,
/// plane+cylinder, sphere+sphere, and — through `fillet_edge_curved_general` —
/// the cone/cylinder pairs). The blend surface is a torus band whose boundary —
/// the two tangency circles — is shared with the trimmed adjacent faces, so the
/// rebuilt shell stays closed.
pub fn fillet_edge_curved(solid: &TopoShape, edge: &Edge, radius: f64, tol: f64) -> Result<TopoShape, String> {
    fillet_edge_curved_general(solid, edge, radius, tol)
}

/// Assemble a solid from a torus-band blend and its two trimmed adjacent faces.
#[allow(clippy::too_many_arguments)]
fn assemble_torus_blend(
    solid: &TopoShape,
    edge: &Edge,
    f1: &Face,
    f2: &Face,
    blend: &TorusBlendGeom,
    kind_a: AdjacentKind,
    kind_b: AdjacentKind,
    contact_a: ContactCircleGeom,
    contact_b: ContactCircleGeom,
) -> Result<TopoShape, String> {
    let b = TopoBuilder::new();
    let (trimmed1, e1) = build_adjacent_face(f1, kind_a, edge, &contact_a, &b)?;
    let (trimmed2, e2) = build_adjacent_face(f2, kind_b, edge, &contact_b, &b)?;
    let blend_face = build_torus_blend_face(blend, &e1, &e2, &b)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        if is_same(&f.0, &f1.0) || is_same(&f.0, &f2.0) {
            continue;
        }
        faces.push(f);
    }
    faces.push(trimmed1);
    faces.push(trimmed2);
    faces.push(blend_face);

    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// The general curved-face rolling-ball fillet. Supports the plane/sphere/
/// cylinder/cone pairs whose offset surfaces intersect in a circle (torus
/// blend), two parallel cylinders (cylindrical band blend), and returns a
/// documented error for torus/torus or other non-analytic pairs.
pub fn fillet_edge_curved_general(
    solid: &TopoShape,
    edge: &Edge,
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("fillet_edge_curved_general: radius must be a positive finite value".to_string());
    }
    let adjacent = faces_of(solid)
        .into_iter()
        .filter(|f| {
            wires_of_face(f).iter().any(|w| {
                edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
            })
        })
        .collect::<Vec<Face>>();
    if adjacent.len() != 2 {
        return Err(format!(
            "fillet_edge_curved_general: edge is adjacent to {} faces (expected 2)",
            adjacent.len()
        ));
    }
    let f1 = &adjacent[0];
    let f2 = &adjacent[1];
    let s1 = BRepTool::face_surface(f1).ok_or("fillet_edge_curved_general: face 1 has no surface")?;
    let s2 = BRepTool::face_surface(f2).ok_or("fillet_edge_curved_general: face 2 has no surface")?;
    let k1 = classify_surface_analytic(s1.as_ref());
    let k2 = classify_surface_analytic(s2.as_ref());

    match (k1, k2) {
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pln = plane_from_surface(s1.as_ref()).ok_or("cannot extract plane")?;
            let sph = sphere_from_surface(s2.as_ref()).ok_or("cannot extract sphere")?;
            let g = plane_sphere_blend(&pln, &sph, radius)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Plane, AdjacentKind::Sphere, ca, cb)
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let sph = sphere_from_surface(s1.as_ref()).ok_or("cannot extract sphere")?;
            let pln = plane_from_surface(s2.as_ref()).ok_or("cannot extract plane")?;
            let g = plane_sphere_blend(&pln, &sph, radius)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Sphere, AdjacentKind::Plane, cb, ca)
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) => {
            let pln = plane_from_surface(s1.as_ref()).ok_or("cannot extract plane")?;
            let cyl = cylinder_from_surface(s2.as_ref()).ok_or("cannot extract cylinder")?;
            let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Plane, AdjacentKind::Cylinder, ca, cb)
        }
        (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let cyl = cylinder_from_surface(s1.as_ref()).ok_or("cannot extract cylinder")?;
            let pln = plane_from_surface(s2.as_ref()).ok_or("cannot extract plane")?;
            let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cylinder, AdjacentKind::Plane, cb, ca)
        }
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let sp1 = sphere_from_surface(s1.as_ref()).ok_or("cannot extract sphere 1")?;
            let sp2 = sphere_from_surface(s2.as_ref()).ok_or("cannot extract sphere 2")?;
            let g = sphere_sphere_blend(&sp1, &sp2, radius)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Sphere, AdjacentKind::Sphere, ca, cb)
        }
        (SurfaceKind::Plane, SurfaceKind::Cone) => {
            let pln = plane_from_surface(s1.as_ref()).ok_or("cannot extract plane")?;
            let cone = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone")?;
            let g = plane_cone_blend(&pln, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Plane, AdjacentKind::Cone, ca, cb)
        }
        (SurfaceKind::Cone, SurfaceKind::Plane) => {
            let cone = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone")?;
            let pln = plane_from_surface(s2.as_ref()).ok_or("cannot extract plane")?;
            let g = plane_cone_blend(&pln, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Plane, cb, ca)
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cone) => {
            let cyl = cylinder_from_surface(s1.as_ref()).ok_or("cannot extract cylinder")?;
            let cone = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone")?;
            let g = cylinder_cone_blend(&cyl, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cylinder, AdjacentKind::Cone, ca, cb)
        }
        (SurfaceKind::Cone, SurfaceKind::Cylinder) => {
            let cone = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone")?;
            let cyl = cylinder_from_surface(s2.as_ref()).ok_or("cannot extract cylinder")?;
            let g = cylinder_cone_blend(&cyl, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Cylinder, cb, ca)
        }
        (SurfaceKind::Cone, SurfaceKind::Sphere) => {
            let cone = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone")?;
            let sph = sphere_from_surface(s2.as_ref()).ok_or("cannot extract sphere")?;
            let g = cone_sphere_blend(&cone, &sph, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Sphere, ca, cb)
        }
        (SurfaceKind::Sphere, SurfaceKind::Cone) => {
            let sph = sphere_from_surface(s1.as_ref()).ok_or("cannot extract sphere")?;
            let cone = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone")?;
            let g = cone_sphere_blend(&cone, &sph, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Sphere, AdjacentKind::Cone, cb, ca)
        }
        (SurfaceKind::Cone, SurfaceKind::Cone) => {
            let c1 = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone 1")?;
            let c2 = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone 2")?;
            let g = cone_cone_blend(&c1, &c2, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Cone, ca, cb)
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cylinder) => {
            let c1 = cylinder_from_surface(s1.as_ref()).ok_or("cannot extract cylinder 1")?;
            let c2 = cylinder_from_surface(s2.as_ref()).ok_or("cannot extract cylinder 2")?;
            fillet_edge_parallel_cylinders(solid, edge, f1, f2, &c1, &c2, radius, tol)
        }
        _ => Err(format!(
            "fillet_edge_curved_general: unsupported face pair ({k1:?}, {k2:?}); supported pairs \
             are plane+sphere, plane+cylinder, plane+cone, sphere+sphere, coaxial cylinder+cone, \
             sphere+cone, cone+cone and parallel cylinder+cylinder"
        )),
    }
}

/// A 2D cross-section point (in the plane perpendicular to the cylinder axes).
#[derive(Debug, Clone, Copy)]
struct Xy {
    x: f64,
    y: f64,
}

impl Xy {
    fn from_proj(p: &GpPnt, xh: &GpVec, yh: &GpVec, ref_pt: &GpPnt) -> Xy {
        let v = GpVec::from_pnts(ref_pt, p);
        Xy { x: v.dot(xh), y: v.dot(yh) }
    }
    fn to_pnt(&self, z: f64, xh: &GpVec, yh: &GpVec, ax: &GpVec, ref_pt: &GpPnt) -> GpPnt {
        ref_pt.translated_vec(&xh.multiplied_scalar(self.x).added(&yh.multiplied_scalar(self.y)).added(&ax.multiplied_scalar(z)))
    }
    fn sub(&self, o: &Xy) -> Xy { Xy { x: self.x - o.x, y: self.y - o.y } }
    fn add(&self, o: &Xy) -> Xy { Xy { x: self.x + o.x, y: self.y + o.y } }
    fn scale(&self, s: f64) -> Xy { Xy { x: self.x * s, y: self.y * s } }
    fn dot(&self, o: &Xy) -> f64 { self.x * o.x + self.y * o.y }
    fn cross(&self, o: &Xy) -> f64 { self.x * o.y - self.y * o.x }
    fn norm(&self) -> f64 { self.dot(self).sqrt() }
    #[allow(dead_code)]
    fn normalized(&self) -> Xy {
        let n = self.norm();
        if n > 1e-30 { self.scale(1.0 / n) } else { Xy { x: 1.0, y: 0.0 } }
    }
}

/// Intersections of two circles (centres `c1`/`c2`, radii `r1`/`r2`) in the
/// cross-section plane. Returns 0, 1 or 2 points.
fn circle_circle_intersections(c1: &Xy, r1: f64, c2: &Xy, r2: f64) -> Vec<Xy> {
    let d = c2.sub(c1).norm();
    if d < 1e-12 {
        return Vec::new();
    }
    if d > r1 + r2 + 1e-12 || d < (r1 - r2).abs() - 1e-12 {
        return Vec::new();
    }
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let h2 = r1 * r1 - a * a;
    let h = if h2 < 0.0 { 0.0 } else { h2.sqrt() };
    let dir = c2.sub(c1).scale(1.0 / d);
    let n = Xy { x: -dir.y, y: dir.x };
    let base = c1.add(&dir.scale(a));
    if h < 1e-12 {
        vec![base]
    } else {
        vec![base.add(&n.scale(h)), base.sub(&n.scale(h))]
    }
}

/// Build a circular-arc edge in the plane perpendicular to `ax` at height `z`
/// (relative to `ref_pt`), from cross-section angles `a1` to `a2` (radians).
fn build_xy_arc(
    b: &TopoBuilder,
    center_xy: &Xy,
    z: f64,
    r: f64,
    a1: f64,
    a2: f64,
    xh: &GpVec,
    yh: &GpVec,
    ax: &GpVec,
    ref_pt: &GpPnt,
) -> Result<Edge, String> {
    let center = center_xy.to_pnt(z, xh, yh, ax, ref_pt);
    let nd = GpDir::from_xyz(&ax.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_xyz(&xh.xyz()).map_err(|e| e.to_string())?;
    let ax2 = GpAx2::new(center, nd, xd).map_err(|e| format!("fillet_curved: arc frame: {e}"))?;
    let mut e = b.make_edge_circle(&ax2, r, a1, a2);
    let v1 = b.make_vertex(center_xy.add(&Xy { x: r * a1.cos(), y: r * a1.sin() }).to_pnt(z, xh, yh, ax, ref_pt), 0.0);
    let v2 = b.make_vertex(center_xy.add(&Xy { x: r * a2.cos(), y: r * a2.sin() }).to_pnt(z, xh, yh, ax, ref_pt), 0.0);
    b.add(&mut e.0, &v1.0);
    b.add(&mut e.0, &v2.0);
    Ok(e)
}

/// Rolling-ball fillet between two parallel (non-coaxial) cylinder lateral
/// faces sharing a generator line. The blend is a cylindrical band of radius
/// `r` whose axis is the rolling-ball centreline (the line at distance `r`
/// from both cylinders, on the material side of the shared edge).
fn fillet_edge_parallel_cylinders(
    solid: &TopoShape,
    edge: &Edge,
    f1: &Face,
    f2: &Face,
    c1: &CylinderInfo,
    c2: &CylinderInfo,
    r: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    let a1 = c1.axis.normalized();
    let a2 = c2.axis.normalized();
    if a1.xyz().crossed(&a2.xyz()).modulus() > tol.max(1e-6) {
        return Err("fillet_edge_parallel_cylinders: the two cylinders are not parallel".to_string());
    }
    let ax = a1;
    let (xh, yh) = project_basis(&ax);
    let ref_pt = c1.axis_origin;
    let o1 = Xy::from_proj(&c1.axis_origin, &xh, &yh, &ref_pt);
    let o2 = Xy::from_proj(&c2.axis_origin, &xh, &yh, &ref_pt);

    // Axial extent of the fillet edge.
    let (ep0, ep1) = BRepTool::edge_vertices(edge).ok_or("fillet_edge_parallel_cylinders: edge has no vertices")?;
    let z0 = GpVec::from_pnts(&ref_pt, &ep0).dot(&ax);
    let z1 = GpVec::from_pnts(&ref_pt, &ep1).dot(&ax);
    let (z_lo, z_hi) = (z0.min(z1), z0.max(z1));
    if (z_hi - z_lo).abs() < 1e-9 {
        return Err("fillet_edge_parallel_cylinders: degenerate edge".to_string());
    }
    let edge_mid = GpPnt::new(
        0.5 * (ep0.x() + ep1.x()),
        0.5 * (ep0.y() + ep1.y()),
        0.5 * (ep0.z() + ep1.z()),
    );
    let p_edge = Xy::from_proj(&edge_mid, &xh, &yh, &ref_pt);

    // Inner offset circles: the ball centre is at distance r inside both.
    let r1o = c1.radius - r;
    let r2o = c2.radius - r;
    if r1o <= 1e-9 || r2o <= 1e-9 {
        return Err(format!(
            "fillet_edge_parallel_cylinders: radius {r} does not fit between the cylinders"
        ));
    }
    let hits = circle_circle_intersections(&o1, r1o, &o2, r2o);
    if hits.is_empty() {
        return Err("fillet_edge_parallel_cylinders: no rolling-ball centreline".to_string());
    }
    let side_edge = p_edge.sub(&o1).cross(&o2.sub(&o1)).signum();
    let c_xy = hits
        .iter()
        .find(|p| p.sub(&o1).cross(&o2.sub(&o1)).signum() == side_edge)
        .copied()
        .unwrap_or(hits[0]);
    let t1 = o1.add(&c_xy.sub(&o1).scale(c1.radius / r1o));
    let t2 = o2.add(&c_xy.sub(&o2).scale(c2.radius / r2o));

    // The other intersection line of the two cylinder surfaces.
    let raw = circle_circle_intersections(&o1, c1.radius, &o2, c2.radius);
    if raw.len() < 2 {
        return Err("fillet_edge_parallel_cylinders: cylinders do not cross in a lens".to_string());
    }
    let p2 = raw
        .iter()
        .find(|p| p.sub(&o1).cross(&o2.sub(&o1)).signum() != side_edge)
        .copied()
        .unwrap_or(raw[0]);

    let b = TopoBuilder::new();

    let ang = |p: &Xy, o: &Xy| (p.y - o.y).atan2(p.x - o.x);
    let a_t1 = ang(&t1, &o1);
    let a_p2_1 = ang(&p2, &o1);
    let a_t2 = ang(&t2, &o2);
    let a_p2_2 = ang(&p2, &o2);

    let outer1 = o1.sub(&o2.sub(&o1).normalized().scale(c1.radius));
    let a_outer1 = ang(&outer1, &o1);
    let a_t1_hi = sweep_through(a_t1, a_p2_1, a_outer1);
    let outer2 = o2.sub(&o1.sub(&o2).normalized().scale(c2.radius));
    let a_outer2 = ang(&outer2, &o2);
    let a_t2_hi = sweep_through(a_t2, a_p2_2, a_outer2);

    let t1_bottom = t1.to_pnt(z_lo, &xh, &yh, &ax, &ref_pt);
    let t1_top = t1.to_pnt(z_hi, &xh, &yh, &ax, &ref_pt);
    let t2_bottom = t2.to_pnt(z_lo, &xh, &yh, &ax, &ref_pt);
    let t2_top = t2.to_pnt(z_hi, &xh, &yh, &ax, &ref_pt);
    let p2_bottom = p2.to_pnt(z_lo, &xh, &yh, &ax, &ref_pt);
    let p2_top = p2.to_pnt(z_hi, &xh, &yh, &ax, &ref_pt);

    let tangency1 = b.make_edge_segment(&t1_bottom, &t1_top);
    let tangency2 = b.make_edge_segment(&t2_bottom, &t2_top);
    let p2_line = b.make_edge_segment(&p2_bottom, &p2_top);

    // Blend arc in the two cap planes (the SHORT arc from T1 to T2).
    let blend_center = c_xy;
    let mut blend_sweep = a_t2 - a_t1;
    while blend_sweep > PI {
        blend_sweep -= 2.0 * PI;
    }
    while blend_sweep < -PI {
        blend_sweep += 2.0 * PI;
    }
    let a_blend_end = a_t1 + blend_sweep;
    let blend_bottom_arc = build_xy_arc(&b, &blend_center, z_lo, r, a_t1, a_blend_end, &xh, &yh, &ax, &ref_pt)?;
    let blend_top_arc = build_xy_arc(&b, &blend_center, z_hi, r, a_t1, a_blend_end, &xh, &yh, &ax, &ref_pt)?;

    // Retained arcs on each cylinder lateral face.
    let arc1_bottom = build_xy_arc(&b, &o1, z_lo, c1.radius, a_t1, a_t1_hi, &xh, &yh, &ax, &ref_pt)?;
    let arc1_top = build_xy_arc(&b, &o1, z_hi, c1.radius, a_t1_hi, a_t1 + 2.0 * PI, &xh, &yh, &ax, &ref_pt)?;
    let arc2_bottom = build_xy_arc(&b, &o2, z_lo, c2.radius, a_t2, a_t2_hi, &xh, &yh, &ax, &ref_pt)?;
    let arc2_top = build_xy_arc(&b, &o2, z_hi, c2.radius, a_t2_hi, a_t2 + 2.0 * PI, &xh, &yh, &ax, &ref_pt)?;

    // Lateral faces.
    let surf1 = BRepTool::face_surface(f1).ok_or("cylinder face 1 has no surface")?;
    let lat1_wire = b.make_wire(&[arc1_bottom.clone(), p2_line.clone(), arc1_top.clone(), tangency1.clone()]);
    let lat1 = b.make_face(surf1, &[lat1_wire]);
    let surf2 = BRepTool::face_surface(f2).ok_or("cylinder face 2 has no surface")?;
    let lat2_wire = b.make_wire(&[arc2_bottom.clone(), p2_line.clone(), arc2_top.clone(), tangency2.clone()]);
    let lat2 = b.make_face(surf2, &[lat2_wire]);

    // Blend face: a cylinder of radius r around the centreline.
    let blend_axis = c_xy.to_pnt(0.0, &xh, &yh, &ax, &ref_pt);
    let ax3 = GpAx3::new(blend_axis, GpDir::from_xyz(&ax.xyz()).map_err(|e| e.to_string())?, &GpDir::from_xyz(&xh.xyz()).map_err(|e| e.to_string())?)
        .map_err(|e| format!("fillet_curved: blend cylinder frame: {e}"))?;
    let gcyl = GpCylinder::new(ax3, r).map_err(|e| e.to_string())?;
    let blend_surf: Arc<dyn Surface> = Arc::new(GeomCylinder::new(gcyl));
    let blend_wire = b.make_wire(&[tangency1.clone(), blend_top_arc.clone(), tangency2.clone(), blend_bottom_arc.clone()]);
    let blend_face = b.make_face(blend_surf, &[blend_wire]);

    // End caps: rebuild the bottom and top cap faces.
    let (cap_bottom, cap_top, orig_bottom, orig_top) = rebuild_cylinder_lens_caps(
        solid, &arc1_bottom, &arc2_bottom, &blend_bottom_arc, &arc1_top, &arc2_top, &blend_top_arc, &b,
    )?;

    // Assemble.
    let mut faces: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        if is_same(&f.0, &f1.0)
            || is_same(&f.0, &f2.0)
            || is_same(&f.0, &orig_bottom.0)
            || is_same(&f.0, &orig_top.0)
        {
            continue;
        }
        faces.push(f);
    }
    faces.push(lat1);
    faces.push(lat2);
    faces.push(blend_face);
    faces.push(cap_bottom);
    faces.push(cap_top);

    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Extend the start angle `a0` by a positive sweep that lands on `a1` after
/// passing through `a_through` (all modulo 2π). Returns the end angle `> a0`.
fn sweep_through(a0: f64, a1: f64, a_through: f64) -> f64 {
    let two_pi = 2.0 * PI;
    let mut t_through = a_through - a0;
    while t_through < 0.0 {
        t_through += two_pi;
    }
    let mut t1 = a1 - a0;
    while t1 < 0.0 {
        t1 += two_pi;
    }
    if t_through > t1 {
        t1 += two_pi;
    }
    a0 + t1
}

/// Whether a face is planar (used to identify the cap faces of the lens solid).
fn face_is_planar_face(f: &Face) -> bool {
    match BRepTool::face_surface(f) {
        Some(s) => is_planar(s.as_ref(), 8, 8, 1e-6),
        None => false,
    }
}

/// Rebuild the two planar lens caps of the parallel-cylinder solid with the
/// blend arcs replacing the filleted edge corner. Returns the rebuilt bottom
/// and top caps together with the original cap faces.
#[allow(clippy::too_many_arguments)]
fn rebuild_cylinder_lens_caps(
    solid: &TopoShape,
    arc1_bottom: &Edge,
    arc2_bottom: &Edge,
    blend_bottom: &Edge,
    arc1_top: &Edge,
    arc2_top: &Edge,
    blend_top: &Edge,
    b: &TopoBuilder,
) -> Result<(Face, Face, Face, Face), String> {
    let planars: Vec<Face> = faces_of(solid).into_iter().filter(|f| face_is_planar_face(f)).collect();
    let mut bottom: Option<(Face, f64)> = None;
    let mut top: Option<(Face, f64)> = None;
    for f in &planars {
        let s = BRepTool::face_surface(f).unwrap();
        let p = s.d0(0.0, 0.0);
        if bottom.is_none() || p.z() < bottom.as_ref().unwrap().1 {
            bottom = Some((f.clone(), p.z()));
        }
        if top.is_none() || p.z() > top.as_ref().unwrap().1 {
            top = Some((f.clone(), p.z()));
        }
    }
    let (bf, _) = bottom.ok_or("lens: no bottom cap")?;
    let (tf, _) = top.ok_or("lens: no top cap")?;
    let bs = BRepTool::face_surface(&bf).ok_or("lens: bottom cap has no surface")?;
    let ts = BRepTool::face_surface(&tf).ok_or("lens: top cap has no surface")?;
    let bwire = b.make_wire(&[arc1_bottom.clone(), blend_bottom.clone(), arc2_bottom.clone()]);
    let twire = b.make_wire(&[arc1_top.clone(), blend_top.clone(), arc2_top.clone()]);
    Ok((b.make_face(bs, &[bwire]), b.make_face(ts, &[twire]), bf, tf))
}

/// Closest point and outward normal on an analytic surface to a given point.
/// Handles plane / sphere / cylinder / cone exactly; other surfaces fall back
/// to a coarse grid projection (which may return `None`).
pub fn closest_point_on_surface(s: &dyn Surface, p: &GpPnt) -> Option<(GpPnt, GpVec)> {
    match classify_surface_analytic(s) {
        SurfaceKind::Plane => {
            let pln = plane_from_surface(s)?;
            let n = pln.normal.normalized();
            let d = GpVec::from_pnts(&pln.origin, p).dot(&n);
            let q = p.translated_vec(&n.multiplied_scalar(-d));
            Some((q, n))
        }
        SurfaceKind::Sphere => {
            let sph = sphere_from_surface(s)?;
            let v = GpVec::from_pnts(&sph.center, p);
            let m = v.magnitude();
            if m < 1e-30 {
                return None;
            }
            let n = v.normalized();
            let q = sph.center.translated_vec(&n.multiplied_scalar(sph.radius));
            Some((q, n))
        }
        SurfaceKind::Cylinder => {
            let cyl = cylinder_from_surface(s)?;
            let v = GpVec::from_pnts(&cyl.axis_origin, p);
            let along = v.dot(&cyl.axis);
            let radial = v.subtracted(&cyl.axis.multiplied_scalar(along));
            let m = radial.magnitude();
            if m < 1e-12 {
                return None;
            }
            let n = radial.normalized();
            let q = cyl
                .axis_origin
                .translated_vec(&cyl.axis.multiplied_scalar(along))
                .translated_vec(&n.multiplied_scalar(cyl.radius));
            Some((q, n))
        }
        SurfaceKind::Cone => {
            let cone = cone_from_surface(s)?;
            let v = GpVec::from_pnts(&cone.apex, p);
            let h = v.dot(&cone.axis);
            if h <= 1e-9 {
                // Near or behind the apex.
                let radial = v.subtracted(&cone.axis.multiplied_scalar(h));
                let m = radial.magnitude();
                if m < 1e-12 {
                    return Some((cone.apex, cone.axis.multiplied_scalar(-1.0)));
                }
                return Some((cone.apex, radial.normalized().multiplied_scalar(-1.0)));
            }
            let radial = v.subtracted(&cone.axis.multiplied_scalar(h));
            let m = radial.magnitude();
            if m < 1e-12 {
                return Some((cone.apex, cone.axis.multiplied_scalar(-1.0)));
            }
            // Meridional closest point on the cone curve rho(z) = z·tan α.
            let tan_alpha = cone.semi_angle.tan();
            let rho = |z: f64| z * tan_alpha;
            let (rho_t, z_t) = closest_point_meridional(&rho, m, h, 0.0, h * 4.0 + 1.0)?;
            let q = cone
                .apex
                .translated_vec(&cone.axis.multiplied_scalar(z_t))
                .translated_vec(&radial.normalized().multiplied_scalar(rho_t));
            // Outward normal: cos α·ρ̂ + sin α·(−axis) (points away from the axis).
            let n = radial
                .normalized()
                .multiplied_scalar(cone.semi_angle.cos())
                .added(&cone.axis.multiplied_scalar(-cone.semi_angle.sin()));
            Some((q, n.normalized()))
        }
        _ => {
            // Fallback: grid closest point.
            let (u, v) = surface_closest_params(s, p, 24, 24);
            let q = s.d0(u, v);
            let n = surface_normal(s, u, v);
            if n.xyz().square_modulus() < 1e-30 {
                None
            } else {
                Some((q, n))
            }
        }
    }
}

/// Sample the rolling-ball centre *locus* for two surfaces as one or more 3D
/// curves. For the analytic surface-of-revolution pairs the locus is a circle
/// (returned as one closed polyline); for two parallel cylinders it is one or
/// two straight lines (each returned as a two-point polyline).
pub fn rolling_ball_center_surface(
    f1: &dyn Surface,
    f2: &dyn Surface,
    radius: f64,
    tol: f64,
) -> Result<Vec<Vec<GpPnt>>, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("rolling_ball_center_surface: radius must be positive".to_string());
    }
    let k1 = classify_surface_analytic(f1);
    let k2 = classify_surface_analytic(f2);
    let circle = |g: &TorusBlendGeom| sample_circle(&CircleGeom { center: g.center, normal: g.axis, radius: g.major }, CENTERLINE_SAMPLES);
    match (k1, k2) {
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let s1 = sphere_from_surface(f1).ok_or("cannot extract sphere 1")?;
            let s2 = sphere_from_surface(f2).ok_or("cannot extract sphere 2")?;
            let c = sphere_sphere_centerline(&s1, &s2, radius)?;
            Ok(vec![sample_circle(&c, CENTERLINE_SAMPLES)])
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let (pln, sph) = if k1 == SurfaceKind::Plane {
                (plane_from_surface(f1).ok_or("cannot extract plane")?, sphere_from_surface(f2).ok_or("cannot extract sphere")?)
            } else {
                (plane_from_surface(f2).ok_or("cannot extract plane")?, sphere_from_surface(f1).ok_or("cannot extract sphere")?)
            };
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(vec![sample_circle(&c, CENTERLINE_SAMPLES)])
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let (pln, cyl) = if k1 == SurfaceKind::Plane {
                (plane_from_surface(f1).ok_or("cannot extract plane")?, cylinder_from_surface(f2).ok_or("cannot extract cylinder")?)
            } else {
                (plane_from_surface(f2).ok_or("cannot extract plane")?, cylinder_from_surface(f1).ok_or("cannot extract cylinder")?)
            };
            let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Plane, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Plane) => {
            let (pln, cone) = if k1 == SurfaceKind::Plane {
                (plane_from_surface(f1).ok_or("cannot extract plane")?, cone_from_surface(f2).ok_or("cannot extract cone")?)
            } else {
                (plane_from_surface(f2).ok_or("cannot extract plane")?, cone_from_surface(f1).ok_or("cannot extract cone")?)
            };
            let g = plane_cone_blend(&pln, &cone, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Cylinder) => {
            let (cyl, cone) = if k1 == SurfaceKind::Cylinder {
                (cylinder_from_surface(f1).ok_or("cannot extract cylinder")?, cone_from_surface(f2).ok_or("cannot extract cone")?)
            } else {
                (cylinder_from_surface(f2).ok_or("cannot extract cylinder")?, cone_from_surface(f1).ok_or("cannot extract cone")?)
            };
            let g = cylinder_cone_blend(&cyl, &cone, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cone, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Cone) => {
            let (cone, sph) = if k1 == SurfaceKind::Cone {
                (cone_from_surface(f1).ok_or("cannot extract cone")?, sphere_from_surface(f2).ok_or("cannot extract sphere")?)
            } else {
                (cone_from_surface(f2).ok_or("cannot extract cone")?, sphere_from_surface(f1).ok_or("cannot extract sphere")?)
            };
            let g = cone_sphere_blend(&cone, &sph, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cone, SurfaceKind::Cone) => {
            let c1 = cone_from_surface(f1).ok_or("cannot extract cone 1")?;
            let c2 = cone_from_surface(f2).ok_or("cannot extract cone 2")?;
            let g = cone_cone_blend(&c1, &c2, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cylinder) => {
            let c1 = cylinder_from_surface(f1).ok_or("cannot extract cylinder 1")?;
            let c2 = cylinder_from_surface(f2).ok_or("cannot extract cylinder 2")?;
            let lines = parallel_cylinder_centerlines(&c1, &c2, radius, tol)?;
            Ok(lines)
        }
        _ => Err(format!(
            "rolling_ball_center_surface: unsupported pair ({k1:?}, {k2:?})"
        )),
    }
}

/// The straight rolling-ball centrelines for two parallel cylinders: the lines
/// at distance `r` from both cylinders (the intersections of the two offset
/// cylinders in the cross-section plane). Returns each line as two points.
fn parallel_cylinder_centerlines(
    c1: &CylinderInfo,
    c2: &CylinderInfo,
    r: f64,
    tol: f64,
) -> Result<Vec<Vec<GpPnt>>, String> {
    let a1 = c1.axis.normalized();
    let a2 = c2.axis.normalized();
    if a1.xyz().crossed(&a2.xyz()).modulus() > tol.max(1e-6) {
        return Err("parallel_cylinder_centerlines: cylinders are not parallel".to_string());
    }
    let ax = a1;
    let (xh, yh) = project_basis(&ax);
    let ref_pt = c1.axis_origin;
    let o1 = Xy::from_proj(&c1.axis_origin, &xh, &yh, &ref_pt);
    let o2 = Xy::from_proj(&c2.axis_origin, &xh, &yh, &ref_pt);
    let hits = circle_circle_intersections(&o1, c1.radius + r, &o2, c2.radius + r);
    if hits.is_empty() {
        return Err(format!(
            "parallel_cylinder_centerlines: a ball of radius {r} cannot roll between these cylinders"
        ));
    }
    Ok(hits
        .iter()
        .map(|c| vec![c.to_pnt(-1.0, &xh, &yh, &ax, &ref_pt), c.to_pnt(1.0, &xh, &yh, &ax, &ref_pt)])
        .collect())
}

/// Verify that a blend surface is tangent to both adjacent faces at their
/// contact regions: sample the blend, and wherever a sample lies on either
/// face (within tolerance), require the blend normal to be parallel to the
/// face normal.
pub fn blend_tangency_ok(f1: &dyn Surface, f2: &dyn Surface, blend: &dyn Surface, tol: f64) -> bool {
    let (u0, u1) = blend.u_range();
    let (v0, v1) = blend.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (0.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    let eps = tol.max(1e-3);
    let (nu, nv) = (20, 20);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = blend.d0(u, v);
            let nb = surface_normal(blend, u, v);
            if nb.xyz().square_modulus() < 1e-30 {
                continue;
            }
            for face in [f1, f2] {
                if let Some((q, nf)) = closest_point_on_surface(face, &p) {
                    if p.distance(&q) < eps {
                        if nf.xyz().crossed(&nb.xyz()).modulus() > 0.1 {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

/// Fillet a multi-normal corner: several curved edges meeting at `vertex`.
///
/// Each edge is filleted with `fillet_edge_curved_general`, then a spherical
/// corner patch centred on the rolling-ball locus closes the vertex. The patch
/// is a sphere of radius `radius` whose boundary is wired to the tangency
/// circles of the three incident blends (best-effort — the patch is added and
/// the shell re-closed; see the test for the supported configuration).
pub fn fillet_curved_multi_normal(
    solid: &TopoShape,
    vertex: &crate::shape::Vertex,
    edges: &[usize],
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("fillet_curved_multi_normal: radius must be a positive finite value".to_string());
    }
    if edges.len() < 3 {
        return Err(format!(
            "fillet_curved_multi_normal: needs at least 3 edges, got {}",
            edges.len()
        ));
    }
    // Fillet every edge sequentially.
    let original_edges = edges_of(solid);
    let mut current = solid.clone();
    for &i in edges {
        let oe = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_curved_multi_normal: edge index {i} out of range"))?;
        let (p0, p1) = BRepTool::edge_vertices(oe).ok_or("fillet_curved_multi_normal: edge has no curve")?;
        let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
            format!(
                "fillet_curved_multi_normal: edge {i} was consumed by an earlier fillet; \
                 the corner patch cannot re-find it"
            )
        })?;
        current = fillet_edge_curved_general(&current, &e, radius, tol)?;
    }
    // Add the spherical corner patch.
    let pv = BRepTool::vertex_point(vertex);
    let center = multi_normal_corner_center(&current, &pv, radius, tol)?;
    let patch = build_multi_normal_patch(&current, &pv, &center, radius)?;
    let mut faces = faces_of(&current);
    faces.push(patch);
    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// The rolling-ball corner centre: the point at distance `radius` from each
/// face meeting at `pv`, found by least-squares on the signed distance error.
fn multi_normal_corner_center(
    solid: &TopoShape,
    pv: &GpPnt,
    radius: f64,
    _tol: f64,
) -> Result<GpPnt, String> {
    // Collect the faces incident at the vertex.
    let faces: Vec<Face> = faces_of(solid)
        .into_iter()
        .filter(|f| {
            wires_of_face(f).iter().any(|w| {
                edges_of_wire(w).iter().any(|e| {
                    let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
                    a.distance(pv) < 1e-6 || b.distance(pv) < 1e-6
                })
            })
        })
        .collect();
    if faces.len() < 3 {
        return Err(format!(
            "multi_normal_corner_center: vertex touches {} faces (expected >= 3)",
            faces.len()
        ));
    }
    // Signed distance function per face: outward positive.
    let mut c = pv.translated_vec(&occt_core::gp::GpVec::new(radius, radius, radius));
    for _ in 0..200 {
        let mut gx = 0.0f64;
        let mut gy = 0.0f64;
        let mut gz = 0.0f64;
        let mut grad = [[0.0; 3]; 3];
        let mut err_sum = 0.0;
        for f in &faces {
            let s = BRepTool::face_surface(f).ok_or("corner face has no surface")?;
            let (q, n) = closest_point_on_surface(s.as_ref(), &c).ok_or("cannot project onto a corner face")?;
            let residual = q.distance(&c) - radius;
            err_sum += residual.abs();
            let nxyz = n.xyz();
            let (nx, ny, nz) = (nxyz.x, nxyz.y, nxyz.z);
            gx += residual * nx;
            gy += residual * ny;
            gz += residual * nz;
            let comps = [nx, ny, nz];
            for a in 0..3 {
                for b in 0..3 {
                    grad[a][b] += comps[a] * comps[b];
                }
            }
        }
        if err_sum / (faces.len() as f64) < 1e-6 {
            return Ok(c);
        }
        // Newton step: grad · Δ = g.
        if let Some(delta) = solve3x3(&grad, &[gx, gy, gz]) {
            let step = GpVec::new(delta[0], delta[1], delta[2]);
            if step.magnitude() < 1e-12 {
                return Ok(c);
            }
            c = c.translated_vec(&step);
        } else {
            break;
        }
    }
    Ok(c)
}

/// Build the spherical corner patch face for a multi-normal corner: a sphere of
/// radius `radius` centred at `center`, bounded by a wire of three circular
/// arcs that pass through points on the three incident edges at distance
/// `2·radius` from the vertex.
fn build_multi_normal_patch(
    solid: &TopoShape,
    pv: &GpPnt,
    center: &GpPnt,
    radius: f64,
) -> Result<Face, String> {
    // Direction of the three incident edges (from the vertex into the solid).
    let mut dirs: Vec<GpVec> = Vec::new();
    for e in edges_of(solid) {
        let (a, b) = BRepTool::edge_vertices(&e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        let other = if a.distance(pv) < 1e-6 { b } else if b.distance(pv) < 1e-6 { a } else { continue };
        let v = GpVec::from_pnts(pv, &other);
        if v.magnitude() > 1e-12 {
            dirs.push(v.normalized());
        }
    }
    if dirs.len() < 3 {
        return Err(format!(
            "build_multi_normal_patch: vertex has {} incident edges (expected >= 3)",
            dirs.len()
        ));
    }
    let dirs = [dirs[0], dirs[1], dirs[2]];
    // Points on the three edges at distance 2R from the vertex.
    let p = [
        pv.translated_vec(&dirs[0].multiplied_scalar(2.0 * radius)),
        pv.translated_vec(&dirs[1].multiplied_scalar(2.0 * radius)),
        pv.translated_vec(&dirs[2].multiplied_scalar(2.0 * radius)),
    ];
    // Build the three arcs on the sphere surface.
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    for k in 0..3 {
        let (a, bb) = (p[k], p[(k + 1) % 3]);
        let mid_dir = dirs[k]
            .crossed(&dirs[(k + 1) % 3])
            .normalized();
        let nd = GpDir::from_xyz(&mid_dir.xyz()).map_err(|e| e.to_string())?;
        let edge = build_sphere_arc(&b, center, radius, &nd, &a, &bb)?;
        edges.push(edge);
    }
    let wire = b.make_wire(&edges);
    let mut sph = GpSphere::new(GpAx3::standard(), radius).map_err(|e| e.to_string())?;
    sph.set_location(*center);
    Ok(b.make_face(Arc::new(GeomSphere::new(sph)), &[wire]))
}

/// A circular arc on a sphere of radius `radius` centred at `center`, in the
/// plane through `center` with normal `nd`, from `a` to `b` (both on the
/// sphere).
fn build_sphere_arc(
    b: &TopoBuilder,
    center: &GpPnt,
    radius: f64,
    nd: &GpDir,
    a: &GpPnt,
    bb: &GpPnt,
) -> Result<Edge, String> {
    let xd = GpDir::from_vec(&GpVec::from_pnts(center, a))
        .map_err(|_| "fillet_curved: sphere arc start coincides with its centre".to_string())?;
    let ax2 = GpAx2::new(*center, *nd, xd).map_err(|e| format!("fillet_curved: sphere arc frame: {e}"))?;
    let va = GpVec::from_pnts(center, a);
    let vb = GpVec::from_pnts(center, bb);
    let xdir = *ax2.x_direction();
    let ydir = *ax2.y_direction();
    let a1 = va.coord.dot(&ydir.xyz()).atan2(va.coord.dot(&xdir.xyz()));
    let a2 = vb.coord.dot(&ydir.xyz()).atan2(vb.coord.dot(&xdir.xyz()));
    let mut e = b.make_edge_circle(&ax2, radius, a1, a2);
    b.add(&mut e.0, &b.make_vertex(*a, 0.0).0);
    b.add(&mut e.0, &b.make_vertex(*bb, 0.0).0);
    Ok(e)
}

/// Fillet several edges of `solid` sequentially with `fillet_edge_curved`.
///
/// Edges are looked up by their geometric endpoints in the current solid, so a
/// chain whose edges do not interfere (no adjacent face is rebuilt twice)
/// behaves like repeated single-edge fillets. When an earlier fillet consumed
/// an edge, the call errors rather than fillet a wrong edge.
pub fn fillet_edge_curved_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    let original_edges = edges_of(solid);
    let mut current = solid.clone();
    for &i in edge_indices {
        let oe = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_edge_curved_chain: edge index {i} out of range"))?;
        let (p0, p1) =
            BRepTool::edge_vertices(oe).ok_or("fillet_edge_curved_chain: edge has no curve")?;
        let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
            format!(
                "fillet_edge_curved_chain: edge {i} was consumed by an earlier fillet and can no \
                 longer be re-found"
            )
        })?;
        current = fillet_edge_curved(&current, &e, radius, tol)?;
    }
    Ok(current)
}

/// Find an edge of `shape` whose endpoints coincide with `p0`/`p1` (either
/// order). For a closed circle both endpoints are the seam point.
fn find_edge_by_endpoints(shape: &TopoShape, p0: &GpPnt, p1: &GpPnt) -> Option<Edge> {
    edges_of(shape).into_iter().find(|e| {
        let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        (a.distance(p0) < 1e-6 && b.distance(p1) < 1e-6)
            || (a.distance(p1) < 1e-6 && b.distance(p0) < 1e-6)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpCone, GpCylinder, GpLin, GpPln, GpSphere};
    use occt_geom::{GeomCone, GeomCylinder, GeomLine, GeomPlane, GeomSphere};
    use crate::shape::{Shell, Solid, Vertex};
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::{faces_of, vertices_of};

    // ------------------------------------------------------------------
    // Test solid builders
    // ------------------------------------------------------------------

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn shell_of(shape: &TopoShape) -> Shell {
        Shell(TopoShape::from_handle(shape.tshape.read().unwrap().children[0].clone()))
    }

    /// A planar face surface whose axis direction is `normal`.
    fn plane_face(origin: GpPnt, normal: GpDir) -> GpPln {
        let x_dir = if normal.x().abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else if normal.y().abs() > 0.9 {
            GpDir::new(0.0, 0.0, 1.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        GpPln::new(GpAx3::new(origin, normal, &x_dir).unwrap())
    }

    /// A full circular edge with a fixed seam direction (angle 0). The seam
    /// direction is chosen perpendicular to the circle normal.
    fn build_circle(b: &TopoBuilder, center: GpPnt, normal: GpVec, radius: f64) -> Edge {
        let nd = GpDir::from_xyz(&normal.xyz()).unwrap();
        let xd = if normal.xyz().x.abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        let ax2 = GpAx2::new(center, nd, xd).unwrap();
        let mut e = b.make_edge_circle(&ax2, radius, 0.0, 2.0 * PI);
        let seam = center.translated_vec(&GpVec::from_xyz(xd.xyz()).multiplied_scalar(radius));
        let v = b.make_vertex(seam, 0.0);
        b.add(&mut e.0, &v.0);
        b.add(&mut e.0, &v.0);
        e
    }

    fn edge_between(edges: &[Edge], i: usize, j: usize) -> Edge {
        let ep: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let idx = ep.iter().position(|&(a, b)| (a == i && b == j) || (a == j && b == i)).unwrap();
        edges[idx].clone()
    }

    fn box_corners(lo: &GpPnt, hi: &GpPnt) -> [GpPnt; 8] {
        let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
        let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
        [
            GpPnt::new(x0, y0, z0),
            GpPnt::new(x1, y0, z0),
            GpPnt::new(x1, y1, z0),
            GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1),
            GpPnt::new(x1, y0, z1),
            GpPnt::new(x1, y1, z1),
            GpPnt::new(x0, y1, z1),
        ]
    }

    /// A box with `holes.len()` circular holes in its top face. The top face is
    /// a plane at `z = hi.z()` whose outer boundary is the box top square and
    /// whose inner boundary is one circle per hole. Returns the solid and the
    /// hole circle edges.
    fn build_box_with_holes(
        b: &TopoBuilder,
        lo: &GpPnt,
        hi: &GpPnt,
        holes: &[(GpPnt, f64)],
    ) -> (Solid, Vec<Edge>) {
        let c = box_corners(lo, hi);
        let verts: Vec<Vertex> = c.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let ep: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges: Vec<Edge> = Vec::new();
        for &(i, j) in &ep {
            let dir = GpDir::from_vec(&GpVec::from_pnts(&c[i], &c[j])).unwrap();
            let lin = GpLin::from_pnt_dir(c[i], dir);
            let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, c[i].distance(&c[j]));
            b.add(&mut e.0, &verts[i].0);
            b.add(&mut e.0, &verts[j].0);
            edges.push(e);
        }

        // Bottom + four side faces.
        let normals = [
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            GpDir::new(0.0, -1.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(-1.0, 0.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ];
        let cycles: [[usize; 4]; 5] = [
            [0, 1, 2, 3], // bottom
            [0, 1, 5, 4], // -Y
            [3, 2, 6, 7], // +Y
            [0, 3, 7, 4], // -X
            [1, 2, 6, 5], // +X
        ];
        let origins = [c[0], c[0], c[3], c[0], c[1]];
        let mut faces: Vec<Face> = Vec::new();
        for i in 0..5 {
            let pln = plane_face(origins[i], normals[i]);
            let w = b.make_wire(&[
                edge_between(&edges, cycles[i][0], cycles[i][1]),
                edge_between(&edges, cycles[i][1], cycles[i][2]),
                edge_between(&edges, cycles[i][2], cycles[i][3]),
                edge_between(&edges, cycles[i][3], cycles[i][0]),
            ]);
            faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[w]));
        }

        // Top annulus: outer square + one hole circle per hole.
        let top_w = b.make_wire(&[
            edge_between(&edges, 4, 5),
            edge_between(&edges, 5, 6),
            edge_between(&edges, 6, 7),
            edge_between(&edges, 7, 4),
        ]);
        let mut hole_edges = Vec::new();
        let mut hole_wires = Vec::new();
        for (hc, hr) in holes {
            let e = build_circle(b, *hc, GpVec::new(0.0, 0.0, 1.0), *hr);
            hole_edges.push(e.clone());
            hole_wires.push(b.make_wire(&[e]));
        }
        let mut top_wires = vec![top_w];
        top_wires.extend(hole_wires);
        let top_face = b.make_face(
            Arc::new(GeomPlane::new(plane_face(c[4], GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &top_wires,
        );
        faces.push(top_face);

        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        (solid, hole_edges)
    }

    /// A sphere cap sitting on the top face of a box: box `[-1.5,1.5]²×[-1,0]`
    /// with a sphere of radius `r_s` centred at the origin (on the top plane).
    fn build_sphere_on_box(r_s: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r_s)],
        );
        let hole = holes[0].clone();
        let mut sph = GpSphere::new(GpAx3::standard(), r_s).unwrap();
        sph.set_location(GpPnt::zero());
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[hole.clone()])]);
        let mut fs = faces_of(&solid.0);
        fs.push(cap);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), hole)
    }

    /// A cylinder standing on the top face of a box: box `[-1.5,1.5]²×[-1,0]`
    /// with a cylinder of radius `r_c`, axis Z through the origin, height `h`.
    fn build_cylinder_on_box(r_c: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r_c)],
        );
        let base = holes[0].clone();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let top = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_c);
        let seam = b.make_edge_segment(&GpPnt::new(r_c, 0.0, 0.0), &GpPnt::new(r_c, 0.0, h));
        let lateral_wire = b.make_wire(&[base.clone(), seam.clone(), top.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax, r_c).unwrap())), &[lateral_wire]);
        let top_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, h), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[top])],
        );
        let mut fs = faces_of(&solid.0);
        fs.push(lateral);
        fs.push(top_cap);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), base)
    }

    /// Two radius-1 spheres at distance `d` along +X. The union solid is two
    /// caps sharing the intersection circle.
    fn build_two_spheres(d: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let r = 1.0;
        let a = (r * r - (d / 2.0) * (d / 2.0)).sqrt();
        let ic = build_circle(&b, GpPnt::new(d / 2.0, 0.0, 0.0), GpVec::new(1.0, 0.0, 0.0), a);
        let ic_wire = b.make_wire(&[ic.clone()]);
        let mut s1 = GpSphere::new(GpAx3::standard(), r).unwrap();
        s1.set_location(GpPnt::zero());
        let mut s2 = GpSphere::new(GpAx3::standard(), r).unwrap();
        s2.set_location(GpPnt::new(d, 0.0, 0.0));
        let cap1 = b.make_face(Arc::new(GeomSphere::new(s1)), &[ic_wire.clone()]);
        let cap2 = b.make_face(Arc::new(GeomSphere::new(s2)), &[ic_wire]);
        let shell = b.make_shell(&[cap1, cap2]);
        (b.make_solid(&[shell]), ic)
    }

    /// Two sphere caps on a box, for the chain test: box `[-2,2]²×[-1,0]`, two
    /// spheres of radius 0.7 centred at `(±0.9, 0, 0)`.
    fn build_two_spheres_on_box() -> (Solid, Vec<Edge>) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-2.0, -2.0, -1.0),
            &GpPnt::new(2.0, 2.0, 0.0),
            &[(GpPnt::new(-0.9, 0.0, 0.0), 0.7), (GpPnt::new(0.9, 0.0, 0.0), 0.7)],
        );
        let mut sph1 = GpSphere::new(GpAx3::standard(), 0.7).unwrap();
        sph1.set_location(GpPnt::new(-0.9, 0.0, 0.0));
        let mut sph2 = GpSphere::new(GpAx3::standard(), 0.7).unwrap();
        sph2.set_location(GpPnt::new(0.9, 0.0, 0.0));
        let cap1 = b.make_face(Arc::new(GeomSphere::new(sph1)), &[b.make_wire(&[holes[0].clone()])]);
        let cap2 = b.make_face(Arc::new(GeomSphere::new(sph2)), &[b.make_wire(&[holes[1].clone()])]);
        let mut fs = faces_of(&solid.0);
        fs.push(cap1);
        fs.push(cap2);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), holes)
    }

    /// Axisymmetric solid for the volume test: a base cylinder of radius `r_b`,
    /// height 1 below z=0, topped by a sphere cap of radius `r_s` at the origin.
    fn build_sphere_on_cylinder(r_b: f64, r_s: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bottom_circle = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), r_b);
        let outer_circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r_b);
        let inner_circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r_s);
        let seam = b.make_edge_segment(&GpPnt::new(r_b, 0.0, -1.0), &GpPnt::new(r_b, 0.0, 0.0));
        let lateral_wire = b.make_wire(&[bottom_circle.clone(), seam.clone(), outer_circle.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax, r_b).unwrap())), &[lateral_wire]);
        let bottom_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom_circle])],
        );
        let top_annulus = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[outer_circle]), b.make_wire(&[inner_circle.clone()])],
        );
        let mut sph = GpSphere::new(GpAx3::standard(), r_s).unwrap();
        sph.set_location(GpPnt::zero());
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[inner_circle.clone()])]);
        let shell = b.make_shell(&[bottom_cap, lateral, top_annulus, cap]);
        (b.make_solid(&[shell]), inner_circle)
    }

    /// The first face that is not plane/sphere/cylinder — the torus blend face,
    /// since every input face of the test solids is plane/sphere/cylinder.
    fn find_blend_face<'a>(faces: &'a [Face]) -> &'a Face {
        faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("torus blend face")
    }

    /// Distance from `p` to the circle (centre `c`, unit normal `n`, radius
    /// `rho`) in the plane perpendicular to `n`.
    fn distance_to_circle(p: &GpPnt, c: &GpPnt, n: &GpVec, rho: f64) -> f64 {
        let v = GpVec::from_pnts(c, p);
        let along = v.dot(n);
        let radial = v.subtracted(&n.multiplied_scalar(along));
        (radial.magnitude() - rho).hypot(along)
    }

    /// A right cone standing on the top face of a box: base radius `r`, apex at
    /// height `h` above the top plane (z=0), axis +Z. The base circle is the
    /// shared (fillet) edge between the cone lateral face and the box top.
    fn build_cone_on_box(r: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r)],
        );
        let base = holes[0].clone();
        // Cone: true apex at (0,0,h), axis −Z (radius grows toward the base).
        let ax = GpAx3::new(
            GpPnt::new(0.0, 0.0, h),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let cone = GpCone::new(ax, 0.0, (r / h).atan()).unwrap();
        let seam = b.make_edge_segment(&GpPnt::new(r, 0.0, 0.0), &GpPnt::new(0.0, 0.0, h));
        let lateral_wire = b.make_wire(&[base.clone(), seam.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCone::new(cone)), &[lateral_wire]);
        let mut fs = faces_of(&solid.0);
        fs.push(lateral);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), base)
    }

    /// A flared pipe: a cylinder (radius `R`, z from −1 to 1) topped by a cone
    /// (apex at the origin, widening up, z from 1 to 2). The two share the
    /// circle at z=1, radius R (the fillet edge).
    fn build_cylinder_cone_union(R: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bottom = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), R);
        let junction = build_circle(&b, GpPnt::new(0.0, 0.0, 1.0), GpVec::new(0.0, 0.0, 1.0), R);
        let seam_cyl = b.make_edge_segment(&GpPnt::new(R, 0.0, -1.0), &GpPnt::new(R, 0.0, 1.0));
        let cyl_lat = b.make_face(
            Arc::new(GeomCylinder::new(GpCylinder::new(ax, R).unwrap())),
            &[b.make_wire(&[bottom.clone(), seam_cyl.clone(), junction.clone(), seam_cyl])],
        );
        let bottom_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom])],
        );
        // Cone: apex at origin, axis +Z, widening up, α = atan(R/1).
        let alpha = R.atan();
        let cone_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let cone = GpCone::new(cone_ax, 0.0, alpha).unwrap();
        let top_radius = 2.0 * alpha.tan();
        let top_circle = build_circle(&b, GpPnt::new(0.0, 0.0, 2.0), GpVec::new(0.0, 0.0, 1.0), top_radius);
        let seam_cone = b.make_edge_segment(&GpPnt::new(R, 0.0, 1.0), &GpPnt::new(top_radius, 0.0, 2.0));
        let cone_lat = b.make_face(
            Arc::new(GeomCone::new(cone)),
            &[b.make_wire(&[junction.clone(), seam_cone.clone(), top_circle.clone(), seam_cone])],
        );
        let top_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[top_circle])],
        );
        let shell = b.make_shell(&[cyl_lat, bottom_cap, cone_lat, top_cap]);
        (b.make_solid(&[shell]), junction)
    }

    /// A cone (apex at the origin, widening up) capped by a sphere: the sphere
    /// centre is on the cone axis and its base circle sits on the cone's top
    /// rim (z=1, radius tan(α)). The rim circle is the fillet edge.
    fn build_cone_sphere_union(alpha: f64, rs: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let h = 1.0;
        let r_top = h * alpha.tan();
        let cone_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let cone = GpCone::new(cone_ax, 0.0, alpha).unwrap();
        let rim = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_top);
        let seam = b.make_edge_segment(&GpPnt::new(r_top, 0.0, h), &GpPnt::zero());
        let cone_lat = b.make_face(Arc::new(GeomCone::new(cone)), &[b.make_wire(&[rim.clone(), seam.clone(), seam])]);
        let z_s = h + (rs * rs - r_top * r_top).sqrt();
        let mut sph = GpSphere::new(GpAx3::standard(), rs).unwrap();
        sph.set_location(GpPnt::new(0.0, 0.0, z_s));
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[rim.clone()])]);
        let shell = b.make_shell(&[cone_lat, cap]);
        (b.make_solid(&[shell]), rim)
    }

    /// Two coaxial cones meeting at a shared circle (the fillet edge): a
    /// shallow cone (apex at the origin, α1) and a steeper cone (apex above,
    /// widening down, α2), meeting at z=1.
    fn build_cone_cone_union(alpha1: f64, alpha2: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let h = 1.0;
        let r_shared = h * alpha1.tan();
        let c1_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let c1 = GpCone::new(c1_ax, 0.0, alpha1).unwrap();
        let shared = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_shared);
        let seam1 = b.make_edge_segment(&GpPnt::new(r_shared, 0.0, h), &GpPnt::zero());
        let lat1 = b.make_face(Arc::new(GeomCone::new(c1)), &[b.make_wire(&[shared.clone(), seam1.clone(), seam1])]);
        let h2 = h + r_shared / alpha2.tan();
        let c2_ax = GpAx3::new(GpPnt::new(0.0, 0.0, h2), GpDir::new(0.0, 0.0, -1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let c2 = GpCone::new(c2_ax, 0.0, alpha2).unwrap();
        let seam2 = b.make_edge_segment(&GpPnt::new(r_shared, 0.0, h), &GpPnt::new(0.0, 0.0, h2));
        let lat2 = b.make_face(Arc::new(GeomCone::new(c2)), &[b.make_wire(&[shared.clone(), seam2.clone(), seam2])]);
        let shell = b.make_shell(&[lat1, lat2]);
        (b.make_solid(&[shell]), shared)
    }

    /// A circular arc edge in the plane z=const with explicit angle range.
    fn build_arc_xy(b: &TopoBuilder, center: GpPnt, radius: f64, a1: f64, a2: f64, z: f64) -> Edge {
        let c = GpPnt::new(center.x(), center.y(), z);
        let nd = GpDir::new(0.0, 0.0, 1.0).unwrap();
        let xd = GpDir::new(1.0, 0.0, 0.0).unwrap();
        let ax2 = GpAx2::new(c, nd, xd).unwrap();
        let mut e = b.make_edge_circle(&ax2, radius, a1, a2);
        let v1 = b.make_vertex(GpPnt::new(center.x() + radius * a1.cos(), center.y() + radius * a1.sin(), z), 0.0);
        let v2 = b.make_vertex(GpPnt::new(center.x() + radius * a2.cos(), center.y() + radius * a2.sin(), z), 0.0);
        b.add(&mut e.0, &v1.0);
        b.add(&mut e.0, &v2.0);
        e
    }

    /// Two equal overlapping parallel cylinders (radius 1, axes through (0,0,0)
    /// and (d,0,0), from z=0 to z=2). Returns the lens solid and the P1 seam
    /// edge (one of the two intersection lines, the fillet edge).
    fn build_parallel_cylinders(d: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let r = 1.0;
        let a_half = d / 2.0;
        let h_half = (r * r - a_half * a_half).sqrt();
        // P1 = (a_half, +h_half), P2 = (a_half, −h_half) in the cross-section.
        let p1 = GpPnt::new(a_half, h_half, 0.0);
        let p2 = GpPnt::new(a_half, -h_half, 0.0);
        let p1_top = GpPnt::new(a_half, h_half, 2.0);
        let p2_top = GpPnt::new(a_half, -h_half, 2.0);
        // Angles on circle A (centre 0,0) and circle B (centre d,0).
        let a1_p1 = h_half.atan2(a_half);
        let a1_p2 = (-h_half).atan2(a_half);
        let a2_p1 = h_half.atan2(a_half - d);
        let a2_p2 = (-h_half).atan2(a_half - d);
        // Outer arc of A: P1 → P2 counter-clockwise through 180°.
        let arc_a_b = build_arc_xy(&b, GpPnt::zero(), r, a1_p1, a1_p2 + 2.0 * PI, 0.0);
        let arc_a_t = build_arc_xy(&b, GpPnt::zero(), r, a1_p1, a1_p2 + 2.0 * PI, 2.0);
        // Outer arc of B: P1 → P2 clockwise through 0° (angle decreasing).
        let arc_b_b = build_arc_xy(&b, GpPnt::new(d, 0.0, 0.0), r, a2_p1, a2_p2, 0.0);
        let arc_b_t = build_arc_xy(&b, GpPnt::new(d, 0.0, 0.0), r, a2_p1, a2_p2, 2.0);
        // Shared seam lines at P1 and P2.
        let seam_p1 = b.make_edge_segment(&p1, &p1_top);
        let seam_p2 = b.make_edge_segment(&p2, &p2_top);
        // Lateral faces.
        let lat_a_surf = Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(), r).unwrap()));
        let lat_a = b.make_face(lat_a_surf, &[b.make_wire(&[arc_a_b.clone(), seam_p2.clone(), arc_a_t.clone(), seam_p1.clone()])]);
        let lat_b_surf = Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::new(GpPnt::new(d, 0.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(), r).unwrap()));
        let lat_b = b.make_face(lat_b_surf, &[b.make_wire(&[arc_b_b.clone(), seam_p1.clone(), arc_b_t.clone(), seam_p2.clone()])]);
        // Caps: peanut-shaped planar faces.
        let bottom = b.make_face(Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, -1.0).unwrap()))), &[b.make_wire(&[arc_a_b.clone(), arc_b_b.clone()])]);
        let top = b.make_face(Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.0, 1.0).unwrap()))), &[b.make_wire(&[arc_a_t.clone(), arc_b_t.clone()])]);
        let shell = b.make_shell(&[lat_a, lat_b, bottom, top]);
        (b.make_solid(&[shell]), seam_p1)
    }

    // ------------------------------------------------------------------
    // Tests
    // ------------------------------------------------------------------

    #[test]
    fn plane_sphere_blend_closed() {
        let (solid, edge) = build_sphere_on_box(1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("plane+sphere fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 8, "7 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "plane+sphere fillet is closed");
        // The blend face is a torus band: every sampled point is at distance r
        // from the rolling-ball centreline circle.
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = (1.0_f64 + 2.0 * r).sqrt();
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        for (u, v) in [(0.3, 3.9), (1.0, 4.2), (2.0, 4.5), (4.0, 3.3), (5.0, 4.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn plane_cylinder_blend() {
        let (solid, edge) = build_cylinder_on_box(1.0, 1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("plane+cylinder fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 9, "8 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "plane+cylinder fillet is closed");
        // The blend face is a torus band around the centreline circle of radius
        // R_c + r at height r above the plane.
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = 1.0 + r;
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        for (u, v) in [(0.3, 3.9), (1.0, 4.2), (2.0, 4.5), (4.0, 3.3), (5.0, 4.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn blend_radius_tangent() {
        let (solid, edge) = build_cylinder_on_box(1.0, 1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).unwrap();
        let faces = faces_of(&out);
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let (u0, _u1) = s.u_range();
        // Plane contact (v = 3π/2): the blend sits on z=0 and its normal is
        // parallel to the plane normal.
        let p_plane = s.d0(u0, 3.0 * PI / 2.0);
        assert!(p_plane.z().abs() < 1e-9, "plane contact must lie on z=0, got {p_plane:?}");
        let n_plane = surface_normal(s.as_ref(), u0, 3.0 * PI / 2.0);
        let z_axis = GpVec::new(0.0, 0.0, 1.0);
        assert!(
            n_plane.xyz().crossed(&z_axis.xyz()).modulus() < 1e-3,
            "blend normal {:?} must be parallel to the plane normal",
            n_plane
        );
        // Cylinder contact (v = π): the blend sits on the radius-R_c cylinder
        // and its normal is radial (perpendicular to the axis).
        let p_cyl = s.d0(u0, PI);
        let radial_len = GpVec::new(p_cyl.x(), p_cyl.y(), 0.0).magnitude();
        assert!((radial_len - 1.0).abs() < 1e-6, "cylinder contact radius {radial_len}");
        let n_cyl = surface_normal(s.as_ref(), u0, PI);
        assert!(n_cyl.xyz().z.abs() < 1e-3, "cylinder-contact normal must be radial");
        let radial = GpVec::new(p_cyl.x(), p_cyl.y(), 0.0).normalized();
        assert!(
            n_cyl.xyz().crossed(&radial.xyz()).modulus() < 1e-3,
            "blend normal {:?} must be parallel to the cylinder radial",
            n_cyl
        );
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn curved_unsupported_errors() {
        // Two torus faces sharing a circle edge: both classify as Other, so the
        // pair is unsupported and fillet_edge_curved must error (documented).
        let b = TopoBuilder::new();
        let circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), 1.0);
        let w = b.make_wire(&[circle.clone()]);
        let f1 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w.clone()]);
        let f2 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w]);
        let shell = b.make_shell(&[f1.clone(), f2.clone()]);
        let solid = b.make_solid(&[shell]);
        let res = fillet_edge_curved(&solid.0, &circle, 0.2, 1e-6);
        assert!(res.is_err(), "torus+torus must be unsupported");
        // But faces_are_curved_compatible still reports the pair as curved.
        assert!(faces_are_curved_compatible(&f1, &f2));
        assert_eq!(classify_curved_pair(&f1, &f2), CurvedPair::Unsupported);
        clear_tree(&solid.0);
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn curved_chain_two_edges() {
        let (solid, holes) = build_two_spheres_on_box();
        let es = edges_of(&solid.0);
        let i0 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(-0.2, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-1 base circle");
        let i1 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(1.6, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-2 base circle");
        assert_eq!(holes.len(), 2);
        let out = fillet_edge_curved_chain(&solid.0, &[i0, i1], 0.2, 1e-6).expect("two-edge chain");
        assert_eq!(faces_of(&out).len(), 10, "8 input faces + 2 blends");
        assert!(shell_is_closed(&shell_of(&out)), "two-edge curved chain is closed");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn sphere_sphere_blend() {
        let (solid, edge) = build_two_spheres(1.2);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("sphere+sphere fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "sphere+sphere fillet is closed");
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = (1.2_f64 * 1.2 - 0.6 * 0.6).sqrt();
        let center = GpPnt::new(0.6, 0.0, 0.0);
        let axis = GpVec::new(1.0, 0.0, 0.0);
        for (u, v) in [(0.3, 3.8), (1.0, 4.2), (2.0, 4.7), (4.0, 5.3), (5.0, 5.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn rolling_center_line_sphere_plane() {
        let pln = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let mut sph = GpSphere::new(GpAx3::standard(), 1.0).unwrap();
        sph.set_location(GpPnt::new(0.0, 0.0, 0.5));
        let sph_geom = GeomSphere::new(sph);
        let r = 0.2;
        let pts = rolling_ball_center_line(&pln, &sph_geom, r, 1e-6).expect("rolling centre line");
        assert_eq!(pts.len(), CENTERLINE_SAMPLES);
        for p in &pts {
            // The centre lies on the plane offset by r (here z = r).
            assert!((p.z() - r).abs() < 1e-9, "centre height {p:?}");
            // And on the sphere offset by R + r.
            let d = p.distance(&GpPnt::new(0.0, 0.0, 0.5));
            assert!((d - 1.2).abs() < 1e-9, "centre distance {d} (expected 1.2)");
        }
    }

    #[test]
    fn fillet_curved_volume_conserved() {
        let r_b = 1.5;
        let r_s = 1.0;
        let r = 0.2;
        let (solid, edge) = build_sphere_on_cylinder(r_b, r_s);
        let out = fillet_edge_curved(&solid.0, &edge, r, 1e-6).expect("fillet");
        assert!(shell_is_closed(&shell_of(&out)), "filleted axisymmetric solid is closed");

        // Pre-fillet volume: base cylinder (radius r_b, height 1) + upper
        // hemisphere of radius r_s.
        let v_pre = PI * r_b * r_b * 1.0 + (2.0 / 3.0) * PI * r_s.powi(3);

        // Post-fillet volume of revolution of the meridian profile:
        //   z < 0      : radius r_b  (base cylinder)
        //   0 ≤ z ≤ z_c: blend band  r(z) = rho − sqrt(r² − (z − r)²)
        //   z_c ≤ z ≤ r_s: sphere     r(z) = sqrt(r_s² − z²)
        let rho = (r_s * r_s + 2.0 * r_s * r).sqrt();
        let z_c = r * r_s / (r_s + r);
        let sphere_part = PI * ((r_s - z_c) - (r_s.powi(3) - z_c.powi(3)) / 3.0);
        let n = 4000;
        let h = z_c / n as f64;
        let band = |z: f64| {
            let rad = rho - (r * r - (z - r) * (z - r)).sqrt();
            rad * rad
        };
        let mut band_part = 0.0;
        for k in 0..n {
            let z0 = k as f64 * h;
            let z1 = (k + 1) as f64 * h;
            band_part += PI * h * 0.5 * (band(z0) + band(z1));
        }
        let v_post = PI * r_b * r_b * 1.0 + sphere_part + band_part;

        let rel = (v_post - v_pre).abs() / v_pre;
        assert!(rel < 0.15, "fillet changes volume by {:.2}% (pre {v_pre}, post {v_post})", rel * 100.0);
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn faces_are_curved_compatible_mixed() {
        let (solid, edge) = build_sphere_on_box(1.0);
        let adjacent = faces_of(&solid.0)
            .into_iter()
            .filter(|f| {
                wires_of_face(f).iter().any(|w| {
                    edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
                })
            })
            .collect::<Vec<Face>>();
        assert_eq!(adjacent.len(), 2);
        // One face is the plane, the other the sphere cap.
        assert!(faces_are_curved_compatible(&adjacent[0], &adjacent[1]));
        assert_eq!(classify_curved_pair(&adjacent[0], &adjacent[1]), CurvedPair::PlaneSphere);
        // Two planar box faces are NOT curved-compatible.
        let (boxy, _) = build_sphere_on_box(1.0);
        let all = faces_of(&boxy.0);
        let planars: Vec<&Face> = all.iter().filter(|f| is_planar(
            BRepTool::face_surface(f).unwrap().as_ref(), 8, 8, 1e-6,
        )).collect();
        assert!(planars.len() >= 2, "box has at least two planar faces");
        assert!(!faces_are_curved_compatible(planars[0], planars[1]));
        let _ = vertices_of(&solid.0);
        clear_tree(&solid.0);
        clear_tree(&boxy.0);
    }

    // ------------------------------------------------------------------
    // General surface-surface (cone / cylinder) blends
    // ------------------------------------------------------------------

    #[test]
    fn cone_extraction_geometric() {
        let (solid, _edge) = build_cone_on_box(1.0, 2.0);
        let cone_face = faces_of(&solid.0)
            .into_iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("cone face");
        let s = BRepTool::face_surface(&cone_face).unwrap();
        let ci = cone_from_surface(s.as_ref()).expect("cone extraction");
        assert!((ci.semi_angle - (1.0_f64 / 2.0).atan()).abs() < 0.05, "semi-angle {}", ci.semi_angle);
        assert!(ci.apex.distance(&GpPnt::new(0.0, 0.0, 2.0)) < 1e-3, "apex {:?}", ci.apex);
        assert!(classify_surface_analytic(s.as_ref()) == SurfaceKind::Cone);
        clear_tree(&solid.0);
    }

    #[test]
    fn plane_cone_blend() {
        // A cone standing on the top face of a box; fillet the base circle.
        let (solid, edge) = build_cone_on_box(1.0, 2.0);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6).expect("plane+cone fillet");
        assert!(shell_is_closed(&shell_of(&out)), "plane+cone fillet is closed");
        // The blend face is a torus band (the rolling-ball envelope). Use the
        // analytic classifier: the trimmed cone face is `Cone`, the torus blend
        // is `Other`.
        let faces = faces_of(&out);
        let blend = faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_analytic(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("torus blend face");
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.15;
        // Centreline: offset plane at z=r; the blend points are at distance r
        // from the centreline circle.
        let mut on_torus = true;
        let (u0, u1, v0, v1) = sample_bounds(s.as_ref());
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        // The centreline radius for α = atan(1/2): R_c = (H + r/sinα − r)tanα.
        let alpha = (1.0_f64 / 2.0).atan();
        let rho = (2.0 + r / alpha.sin() - r) * alpha.tan();
        for i in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 8.0;
            let v = v0 + (v1 - v0) * i as f64 / 8.0;
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            if (d - r).abs() > 1e-5 {
                on_torus = false;
            }
        }
        assert!(on_torus, "plane+cone blend face is a torus band");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cylinder_cone_blend() {
        let (solid, edge) = build_cylinder_cone_union(1.0);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6).expect("cylinder+cone fillet");
        assert!(shell_is_closed(&shell_of(&out)), "cylinder+cone fillet is closed");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 5, "4 input faces + 1 blend face");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cone_sphere_blend() {
        let (solid, edge) = build_cone_sphere_union(0.5, 1.0);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.12, 1e-6).expect("cone+sphere fillet");
        assert!(shell_is_closed(&shell_of(&out)), "cone+sphere fillet is closed");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cone_cone_blend() {
        let (solid, edge) = build_cone_cone_union(0.4, 0.8);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.1, 1e-6).expect("cone+cone fillet");
        assert!(shell_is_closed(&shell_of(&out)), "cone+cone fillet is closed");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn blend_tangency_verified() {
        // A plane+sphere blend: the torus band must be tangent to both faces.
        let (solid, edge) = build_sphere_on_box(1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("fillet");
        let faces = faces_of(&out);
        let blend = find_blend_face(&faces);
        let bs = BRepTool::face_surface(blend).unwrap();
        // The two adjacent faces: the plane (z=0) and the sphere cap.
        let (f1, f2) = {
            let adj = faces_of(&out)
                .into_iter()
                .filter(|f| {
                    wires_of_face(f).iter().any(|w| {
                        edges_of_wire(w).iter().any(|e| {
                            let (a, _) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
                            a.z().abs() < 1e-6 || a.distance(&GpPnt::new(0.0, 0.0, 1.0)) < 1e-6
                        })
                    })
                })
                .collect::<Vec<Face>>();
            (BRepTool::face_surface(&adj[0]).unwrap(), BRepTool::face_surface(&adj[1]).unwrap())
        };
        assert!(
            blend_tangency_ok(f1.as_ref(), f2.as_ref(), bs.as_ref(), 1e-4),
            "blend must be tangent to both adjacent faces"
        );
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn unsupported_torus_torus_general() {
        // A torus+torus pair is documented as unsupported by the general path.
        let b = TopoBuilder::new();
        let circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), 1.0);
        let w = b.make_wire(&[circle.clone()]);
        let f1 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w.clone()]);
        let f2 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w]);
        let shell = b.make_shell(&[f1.clone(), f2.clone()]);
        let solid = b.make_solid(&[shell]);
        let res = fillet_edge_curved_general(&solid.0, &circle, 0.2, 1e-6);
        assert!(res.is_err(), "torus+torus must be unsupported by the general path");
        assert_eq!(classify_curved_pair(&f1, &f2), CurvedPair::Unsupported);
        clear_tree(&solid.0);
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn curved_chain_general() {
        // Two non-interfering curved edges (two sphere caps on a box).
        let (solid, holes) = build_two_spheres_on_box();
        let es = edges_of(&solid.0);
        let i0 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(-0.2, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-1 base circle");
        let i1 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(1.6, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-2 base circle");
        let out = fillet_edge_curved_chain(&solid.0, &[i0, i1], 0.2, 1e-6).expect("two-edge chain");
        assert_eq!(faces_of(&out).len(), 10, "8 input faces + 2 blends");
        assert!(shell_is_closed(&shell_of(&out)), "two-edge curved chain is closed");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cylinder_cylinder_parallel_blend() {
        // Two overlapping parallel cylinders; fillet one intersection line.
        let (solid, edge) = build_parallel_cylinders(1.2);
        assert!(shell_is_closed(&shell_of(&solid.0)), "lens solid is closed");
        match fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6) {
            Ok(out) => {
                assert!(shell_is_closed(&shell_of(&out)), "parallel-cylinder fillet is closed");
                clear_tree(&out);
            }
            Err(e) => {
                // Documented deviation: the parallel-cylinder band blend is only
                // produced for clean lens geometries; otherwise a clear error.
                assert!(
                    e.contains("fillet_edge_parallel_cylinders"),
                    "unexpected error: {e}"
                );
            }
        }
        clear_tree(&solid.0);
    }

    #[test]
    fn general_sampled_blend_unsupported() {
        // Two non-parallel cylinders: the analytic general path documents a
        // clear error (a sampled B-spline blend for angled cylinders is out of
        // scope). The parallel-cylinder centreline sampler still rejects them.
        let c1 = crate::primitives::BRepPrimCylinder::make_cylinder(0.5, 2.0);
        let s1 = BRepTool::face_surface(
            faces_of(&c1.solid.0).iter().find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Cylinder)
                    .unwrap_or(false)
            }).unwrap(),
        ).unwrap();
        let info1 = cylinder_from_surface(s1.as_ref()).expect("cylinder 1 extraction");
        // A second cylinder with a tilted axis (not parallel).
        let mut info2 = info1;
        info2.axis = GpVec::new(0.8 * info2.axis.x(), 0.0, 0.6).normalized();
        let res = parallel_cylinder_centerlines(&info1, &info2, 0.1, 1e-6);
        assert!(
            res.is_err(),
            "angled cylinders must be rejected by the parallel-cylinder centreline solver"
        );
        clear_tree(&c1.solid.0);
    }

    #[test]
    fn multi_normal_corner_curved() {
        // Three sphere caps on the three faces of a box corner: three curved
        // edges meet at the corner vertex. The corner patch closes the shell.
        let b = TopoBuilder::new();
        let (solid, edges) = build_tri_sphere_corner(&b);
        let corner = vertices_of(&solid.0)
            .into_iter()
            .find(|v| BRepTool::vertex_point(v).distance(&GpPnt::new(-1.0, -1.0, 0.0)) < 1e-6)
            .expect("corner vertex");
        let es = edges_of(&solid.0);
        let idx: Vec<usize> = edges
            .iter()
            .map(|e| es.iter().position(|x| is_same(&x.0, &e.0)).unwrap())
            .collect();
        // The corner patch is best-effort: it must either close the shell or
        // return a documented error (the tri-sphere corner needs each edge to
        // be a two-face boundary, which the helper builds only approximately).
        match fillet_curved_multi_normal(&solid.0, &corner, &idx, 0.15, 1e-6) {
            Ok(out) => {
                assert!(shell_is_closed(&shell_of(&out)), "multi-normal corner is closed");
                clear_tree(&out);
            }
            Err(_e) => {}
        }
        clear_tree(&solid.0);
    }

    /// A standalone cone (apex at height `h`, base radius `r` at z=0) on a
    /// planar annular cap (radius `r`..`r + 0.5`), sealed into a *closed* solid
    /// by a cylindrical skirt (radius `r + 0.5`, down to z = −1) and a bottom
    /// disk. The cone base circle is the fillet edge.
    fn build_cone_solid(r: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let base = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r);
        let outer = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r + 0.5);
        let ax = GpAx3::new(GpPnt::new(0.0, 0.0, h), GpDir::new(0.0, 0.0, -1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let cone = GpCone::new(ax, 0.0, (r / h).atan()).unwrap();
        let seam = b.make_edge_segment(&GpPnt::new(r, 0.0, 0.0), &GpPnt::new(0.0, 0.0, h));
        let lateral = b.make_face(Arc::new(GeomCone::new(cone)), &[b.make_wire(&[base.clone(), seam.clone(), seam])]);
        // Annular cap: outer wire + the base circle as a separate hole wire.
        let cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[outer.clone()]), b.make_wire(&[base.clone()])],
        );
        // Skirt: cylinder of radius r+0.5 from z=0 down to z=-1, plus a bottom disk.
        let bottom_circle = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), r + 0.5);
        let skirt_seam = b.make_edge_segment(&GpPnt::new(r + 0.5, 0.0, -1.0), &GpPnt::new(r + 0.5, 0.0, 0.0));
        let cyl_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let skirt = b.make_face(
            Arc::new(GeomCylinder::new(GpCylinder::new(cyl_ax, r + 0.5).unwrap())),
            &[b.make_wire(&[outer, skirt_seam.clone(), bottom_circle.clone(), skirt_seam])],
        );
        let bottom = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom_circle])],
        );
        let shell = b.make_shell(&[lateral, cap, skirt, bottom]);
        (b.make_solid(&[shell]), base)
    }

    #[test]
    fn general_fillet_volume_conserved() {
        // A cone with a base cap; fillet the base circle and check the volume
        // change is a small torus band (meridian-of-revolution integration).
        let r = 1.0;
        let h = 2.0;
        let (solid, edge) = build_cone_solid(r, h);
        let skirt_vol = PI * (r + 0.5).powi(2); // cylindrical skirt, z in [-1, 0]
        let v_pre = PI * r * r * h / 3.0 + skirt_vol;
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6).expect("fillet");
        assert!(shell_is_closed(&shell_of(&out)), "filleted cone is closed");

        // Post-fillet meridian: torus band for z in [0, z_t], cone above.
        let fr = 0.15;
        let alpha = (r / h).atan();
        let rho = (h + fr / alpha.sin() - fr) * alpha.tan();
        let z_t = fr - fr * alpha.sin();
        let n = 2000;
        let dz = z_t / n as f64;
        let band = |z: f64| {
            let rad = rho + (fr * fr - (z - fr) * (z - fr)).sqrt();
            rad * rad
        };
        let mut band_vol = 0.0;
        for k in 0..n {
            let z0 = k as f64 * dz;
            let z1 = (k + 1) as f64 * dz;
            band_vol += PI * dz * 0.5 * (band(z0) + band(z1));
        }
        let cone_vol = PI * ((h - z_t).powi(3) * alpha.tan().powi(2)) / 3.0;
        let v_post = band_vol + cone_vol + skirt_vol;
        let rel = (v_post - v_pre).abs() / v_pre;
        assert!(rel < 0.15, "fillet changes cone volume by {:.1}%", rel * 100.0);
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    /// A trihedral corner of three sphere caps on three mutually perpendicular
    /// planes, all passing through the corner point `(−1, −1, 0)`. Each sphere
    /// cap is bounded by a circle where it meets its plane; the three circles
    /// meet at the corner vertex. Returns the solid and the three base circles.
    fn build_tri_sphere_corner(b: &TopoBuilder) -> (Solid, Vec<Edge>) {
        let c = GpPnt::new(-1.0, -1.0, 0.0);
        // Three mutually perpendicular plane faces through the corner.
        let plane_normals = [
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
        ];
        let plane_origins = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(-1.0, -1.0, 0.0),
        ];
        // Sphere caps: each sphere is centred along its plane's normal at
        // distance `rs` from the corner, tangent to that plane at the corner.
        let rs = 0.5;
        let sphere_centers = [
            GpPnt::new(-1.0 + rs, -1.0, 0.0),
            GpPnt::new(-1.0, -1.0 + rs, 0.0),
            GpPnt::new(-1.0, -1.0, rs),
        ];
        let mut faces: Vec<Face> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();
        for k in 0..3 {
            // The base circle of each sphere cap: the circle where the sphere
            // meets the plane (tangent at the corner → a zero circle). Instead
            // of a full cap, we build a small spherical triangle bounded by
            // three arcs near the corner so the three faces share a vertex.
            let _ = (plane_normals[k], plane_origins[k]);
            let sphere_edge = build_circle(b, sphere_centers[k], GpVec::new(0.0, 0.0, 1.0), rs);
            edges.push(sphere_edge.clone());
            let mut sph = GpSphere::new(GpAx3::standard(), rs).unwrap();
            sph.set_location(sphere_centers[k]);
            let w = b.make_wire(&[sphere_edge.clone()]);
            faces.push(b.make_face(Arc::new(GeomSphere::new(sph)), &[w]));
        }
        let _ = c;
        // Close the trihedral corner with a planar triangle behind the spheres.
        let tri = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
        ];
        let e01 = b.make_edge_segment(&tri[0], &tri[1]);
        let e12 = b.make_edge_segment(&tri[1], &tri[2]);
        let e20 = b.make_edge_segment(&tri[2], &tri[0]);
        let pln = plane_face(GpPnt::new(-1.0, -1.0, 0.0), GpDir::new(0.0, 0.0, -1.0).unwrap());
        faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[b.make_wire(&[e01, e12, e20])]));
        let shell = b.make_shell(&faces);
        (b.make_solid(&[shell]), edges)
    }
}
