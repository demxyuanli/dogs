use super::prelude::*;
use super::*;

/// TEMPORARY BUILD FIX (concurrent agent): component `i` of a `GpXyz`.

pub(super) fn xyz_component(c: &occt_core::gp::xyz::GpXyz, i: usize) -> f64 {
    match i {
        0 => c.x(),
        1 => c.y(),
        _ => c.z(),
    }
}

/// Number of sampled points on a rolling-ball centreline circle.
pub(super) const CENTERLINE_SAMPLES: usize = 48;

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
pub(super) struct CircleGeom {
    pub(super) center: GpPnt,
    /// Unit circle-plane normal.
    pub(super) normal: GpVec,
    pub(super) radius: f64,
}

/// A tangency circle of the blend on one adjacent face.
#[derive(Debug, Clone, Copy)]
pub(super) struct ContactCircleGeom {
    pub(super) center: GpPnt,
    /// Unit circle-plane normal.
    pub(super) normal: GpVec,
    pub(super) radius: f64,
}

/// Full description of a torus band blend.
#[derive(Debug, Clone)]
pub(super) struct TorusBlendGeom {
    /// Torus centre (on the centreline circle's axis).
    pub(super) center: GpPnt,
    /// Unit torus axis.
    pub(super) axis: GpVec,
    /// Centreline (major) radius.
    pub(super) major: f64,
    /// Fillet (minor) radius.
    pub(super) minor: f64,
    pub(super) contact_a: ContactCircleGeom,
    pub(super) contact_b: ContactCircleGeom,
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
pub(super) fn is_curved_surface(s: &dyn Surface) -> bool {
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
pub(super) fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
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
pub(super) fn smallest_eigenvector(m: &[[f64; 3]; 3]) -> Option<[f64; 3]> {
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
pub(super) fn solve3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
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
pub(super) fn cone_from_surface_nappe(s: &dyn Surface, sign: f64) -> Option<ConeInfo> {
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
pub(super) fn project_basis(axis: &GpVec) -> (GpVec, GpVec) {
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
pub(super) fn perp_dir(z: &GpVec) -> GpVec {
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
pub(super) fn plane_sphere_centerline(pln: &PlaneInfo, sph: &SphereInfo, r: f64) -> Result<CircleGeom, String> {
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
pub(super) fn sphere_sphere_centerline(s1: &SphereInfo, s2: &SphereInfo, r: f64) -> Result<CircleGeom, String> {
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
pub(super) fn sample_circle(c: &CircleGeom, n: usize) -> Vec<GpPnt> {
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
