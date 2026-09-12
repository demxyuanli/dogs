use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Surface classification
// ---------------------------------------------------------------------------

/// Classify a surface, subdividing the analytic revolution types beyond the
/// plane/sphere detection of [`crate::brep_surface::classify_surface`].
///
/// The `cylinder` / `cone` / `torus` checks reuse the geometric-invariant
/// reconstruction from `gprop_analytic.rs` (a trait object cannot be downcast,
/// so the parameters are recovered by sampling).

pub fn classify_surface_kind(s: &dyn Surface) -> SurfaceKind {
    let base = classify_surface(s);
    if base == SurfaceKind::Plane || base == SurfaceKind::Sphere {
        return base;
    }
    if cylinder_params(s).is_some() {
        return SurfaceKind::Cylinder;
    }
    if cone_params(s).is_some() {
        return SurfaceKind::Cone;
    }
    if torus_params(s).is_some() {
        return SurfaceKind::Torus;
    }
    SurfaceKind::Other
}

/// Cylinder parameters `(center, unit axis, radius)` from sampled invariants:
/// the axis is along `∂S/∂v` and two opposite points of a `u`-ring give the
/// radius (their midpoint lies on the axis).
pub(super) fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axis.divided(m);
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let r = q0.distance(&q1) / 2.0;
    if r <= 1e-9 {
        return None;
    }
    let center = mid(&q0, &q1);
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

/// Cone parameters `(apex, unit axis, semi-angle)`: two fixed-`u` generatrices
/// intersect at the apex; every sampled point makes the same angle with the axis.
pub(crate) fn cone_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let g0a = s.d0(0.0, 0.0);
    let g0b = s.d0(0.0, 1.0);
    let d0 = GpVec::from_pnts(&g0a, &g0b);
    if d0.magnitude() < 1e-9 {
        return None;
    }
    let g1a = s.d0(FRAC_PI_2, 0.0);
    let g1b = s.d0(FRAC_PI_2, 1.0);
    let d1 = GpVec::from_pnts(&g1a, &g1b);
    if d1.magnitude() < 1e-9 {
        return None;
    }
    let apex = closest_point_lines(g0a, d0, g1a, d1)?;
    // The cone axis runs through the apex and the midpoint of two opposite
    // points on a u-ring. Sample the ring at v=1 (not v=0 — every u collapses
    // to the apex there, making the axis a zero vector).
    let m0 = mid(&s.d0(0.0, 1.0), &s.d0(PI, 1.0));
    let axis = GpVec::from_pnts(&apex, &m0);
    let am = axis.magnitude();
    if am < 1e-9 {
        return None;
    }
    let ax = axis.divided(am);
    // Use the v=1 point for the apex angle: `g0a` (v=0) may coincide with the
    // apex itself (a cone whose parameter v=0 sits at the vertex), making the
    // vector from apex to g0a a zero vector and the angle undefined.
    let alpha = GpVec::from_pnts(&apex, &g0b).angle(&ax);
    if alpha < 1e-6 || alpha > FRAC_PI_2 - 1e-6 {
        return None;
    }
    let (nu, nv) = (6, 5);
    for i in 0..=nu {
        for j in 1..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            // Sample the valid cone half only (v > 0): the surface's v-range is
            // unbounded but d0 is only well-defined on the cone sheet away from
            // the apex; negative v can land on the opposite (invalid) sheet.
            let v = j as f64 / nv as f64;
            let p = s.d0(u, v);
            if p.distance(&apex) < 1e-9 {
                continue;
            }
            let ang = GpVec::from_pnts(&apex, &p).angle(&ax);
            if (ang - alpha).abs() > 1e-3 * alpha.max(0.01) {
                return None;
            }
        }
    }
    Some((apex, ax, alpha))
}

/// Torus parameters `(center, axis, major radius, minor radius)`.
pub(crate) fn torus_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64, f64)> {
    let top = s.d0(0.0, FRAC_PI_2);
    let bot = s.d0(0.0, -FRAC_PI_2);
    let axv = GpVec::from_pnts(&bot, &top);
    let m = axv.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axv.divided(m);
    let e0 = s.d0(0.0, 0.0);
    let e1 = s.d0(PI, 0.0);
    let center = mid(&e0, &e1);
    let rel = GpVec::from_pnts(&center, &top);
    let rel_axis = rel.dot(&ax);
    let in_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let big_r = in_plane.magnitude();
    let small_r = rel_axis.abs();
    if big_r <= 1e-9 || small_r <= 1e-9 {
        return None;
    }
    let (nu, nv) = (6, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = 2.0 * PI * j as f64 / nv as f64;
            let p = s.d0(u, v);
            let rr = GpVec::from_pnts(&center, &p);
            let ra = rr.dot(&ax);
            // Signed radial distance `R + r·cos v` — not the unsigned magnitude.
            // A spindle torus (R < r) folds through itself, so its inner equator
            // (v ≈ π) has a negative rho; the unsigned |R + r·cos v| breaks the
            // tube equation there. Recover the sign from the torus identity
            // |rr|² = rho² + ra² ⇒ rho = (|rr|² + R² − r²) / (2R), which is exact
            // for both ring (R > r) and spindle (R < r) tori (ElSLib convention).
            let rho = (rr.square_magnitude() + big_r.powi(2) - small_r.powi(2)) / (2.0 * big_r);
            let err = ((rho - big_r).powi(2) + ra.powi(2)).sqrt();
            if (err - small_r).abs() > 1e-3 * small_r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, big_r, small_r))
}

// ---------------------------------------------------------------------------
// Small vector helpers
// ---------------------------------------------------------------------------

pub(super) fn mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

pub(super) fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

pub(super) fn pnt_add_vec(p: GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::from_xyz(&p.coord.added(&v.coord))
}

/// Closest point between two (possibly skew) lines; `None` when parallel or
/// too far apart.
pub(super) fn closest_point_lines(p1: GpPnt, d1: GpVec, p2: GpPnt, d2: GpVec) -> Option<GpPnt> {
    let r = GpVec::from_pnts(&p1, &p2);
    let a = d1.dot(&d1);
    let b = d1.dot(&d2);
    let c = d2.dot(&d2);
    let e = d1.dot(&r);
    let f = d2.dot(&r);
    let denom = a * c - b * b;
    if denom.abs() < 1e-12 {
        return None;
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

// ---------------------------------------------------------------------------
// Analytic surface projector
// ---------------------------------------------------------------------------

/// Exact UV projection of a 3D point onto an analytic surface, recovered from
/// the surface's own parameterization so the results agree with
/// `edge_pcurve_on_face`.
pub(super) enum Projector {
    Plane { o: GpPnt, x: GpVec, y: GpVec },
    Cylinder { o: GpPnt, x: GpVec, y: GpVec, z: GpVec },
    Cone { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, alpha: f64 },
    Sphere { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r: f64 },
    Torus { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r_big: f64, r_small: f64 },
}

impl Projector {
    pub(super) fn from_surface(s: &dyn Surface, kind: SurfaceKind) -> Option<Projector> {
        match kind {
            SurfaceKind::Plane => {
                let (o, x, y) = plane_axes(s)?;
                Some(Projector::Plane { o, x, y })
            }
            SurfaceKind::Cylinder => {
                let (o, x, y, z) = cylinder_axes(s)?;
                Some(Projector::Cylinder { o, x, y, z })
            }
            SurfaceKind::Cone => {
                let (apex, ax, alpha) = cone_params(s)?;
                let (x, y, z) = cone_frame(s, &apex, &ax)?;
                // The projection origin must be the surface's placement plane
                // point (v=0), not its apex: `cone_value` is parameterized from
                // the placement (ElSLib::ConeValue), so an apex-based origin
                // would offset the projected v by RefRadius/sin(α) and the UV
                // window would not be the exact inverse of `surface.d0`.
                let o = cone_placement(s, &apex, &z)?;
                Some(Projector::Cone { o, x, y, z, alpha })
            }
            SurfaceKind::Sphere => {
                let center = sphere_center(s)?;
                let r = s.d0(0.0, 0.0).distance(&center);
                if r < 1e-12 {
                    return None;
                }
                let x = GpVec::from_pnts(&center, &s.d0(0.0, 0.0)).normalized();
                let y = GpVec::from_pnts(&center, &s.d0(FRAC_PI_2, 0.0)).normalized();
                let z = x.crossed(&y).normalized();
                Some(Projector::Sphere { o: center, x, y, z, r })
            }
            SurfaceKind::Torus => {
                let (center, ax, r_big, r_small) = torus_params(s)?;
                let x = GpVec::from_pnts(&center, &s.d0(0.0, 0.0)).normalized();
                let y = GpVec::from_pnts(&center, &s.d0(FRAC_PI_2, 0.0)).normalized();
                // The axis must be the surface's *actual* revolution axis (the
                // `ElSLib::TorusParameters` `gp_Ax3` direction), recovered by
                // `torus_params` from `top − bot` — not `x × y`, which flips for
                // a left-handed STEP torus and would make `v = atan2(z·d, …)` land
                // on the wrong half of the tube (the screw's thread got v < 0 for
                // points above its centre).
                let z = ax;
                Some(Projector::Torus { o: center, x, y, z, r_big, r_small })
            }
            SurfaceKind::Other => None,
        }
    }

    /// Whether the surface's `u` parameter is periodic (revolution surfaces).
    pub(super) fn is_u_periodic(&self) -> bool {
        matches!(
            self,
            Projector::Cylinder { .. } | Projector::Cone { .. } | Projector::Sphere { .. } | Projector::Torus { .. }
        )
    }

    /// Whether the surface's `v` parameter is periodic (torus tube direction).
    pub(super) fn is_v_periodic(&self) -> bool {
        matches!(self, Projector::Torus { .. })
    }

    /// Project a single 3D point onto the face `(u, v)` domain. `u` is wrapped
    /// into the principal angle range `(-π, π]`; use [`Self::project_seq`] to
    /// unwrap a continuous sequence across the seam.
    pub(super) fn project(&self, p: &GpPnt) -> GpPnt2d {
        match self {
            Projector::Plane { o, x, y } => {
                let z = x.crossed(y);
                if let (Ok(zd), Ok(xd)) = (GpDir::from_vec(&z), GpDir::from_vec(x)) {
                    if let Ok(ax) = GpAx3::new(*o, zd, &xd) {
                        return occt_geom::projlib::project_pln_pnt(&GpPln::new(ax), p);
                    }
                }
                let d = p.coord.subtracted(&o.coord);
                let xx = x.xyz().dot(x.xyz());
                let yy = y.xyz().dot(y.xyz());
                let xy = x.xyz().dot(y.xyz());
                let dx = d.dot(x.xyz());
                let dy = d.dot(y.xyz());
                let det = xx * yy - xy * xy;
                if det.abs() < 1e-24 {
                    GpPnt2d::new(dx, dy)
                } else {
                    GpPnt2d::new((dx * yy - dy * xy) / det, (dy * xx - dx * xy) / det)
                }
            }
            Projector::Cylinder { o, x, y, z } => {
                let d = GpVec::from_pnts(o, p);
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), d.dot(z))
            }
            Projector::Cone { o, x, y, z, alpha } => {
                let d = GpVec::from_pnts(o, p);
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), d.dot(z) / alpha.cos())
            }
            Projector::Sphere { o, x, y, z, r } => {
                let d = GpVec::from_pnts(o, p);
                let v = (d.dot(z) / r).clamp(-1.0, 1.0).asin();
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), v)
            }
            Projector::Torus { o, x, y, z, r_big, r_small } => {
                let d = GpVec::from_pnts(o, p);
                let dz = d.dot(z);
                // Signed radial distance `R + r·cos v`. The unsigned |R + r·cos v|
                // folds a spindle torus (R < r) through its inner equator (v ≈ π),
                // where the radial distance is negative; `atan2(dz, |rho| − R)`
                // then returns v ≈ 0 instead of v ≈ π. Recover the sign from the
                // torus identity |d|² = rho² + dz² ⇒ rho = (|d|² + R² − r²)/(2R).
                let rho = (d.square_magnitude() + r_big.powi(2) - r_small.powi(2)) / (2.0 * r_big);
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), dz.atan2(rho - r_big))
            }
        }
    }

    /// Project a sequence of edge points, unwrapping the `u` (and, for a torus,
    /// the `v`) coordinate across the seam so the result stays continuous and on
    /// the same sheet of a self-intersecting surface.
    pub(super) fn project_seq(&self, curve: &dyn Curve, t_vals: &[f64]) -> Vec<GpPnt2d> {
        let mut out = Vec::with_capacity(t_vals.len());
        let mut prev_u = f64::NAN;
        let mut prev_v = f64::NAN;
        for &t in t_vals {
            let q = self.project(&curve.d0(t));
            let mut u = q.x();
            let mut v = q.y();
            if self.is_u_periodic() && prev_u.is_finite() {
                while u - prev_u > PI {
                    u -= 2.0 * PI;
                }
                while u - prev_u < -PI {
                    u += 2.0 * PI;
                }
            }
            if self.is_v_periodic() && prev_v.is_finite() {
                while v - prev_v > PI {
                    v -= 2.0 * PI;
                }
                while v - prev_v < -PI {
                    v += 2.0 * PI;
                }
            }
            prev_u = u;
            prev_v = v;
            out.push(GpPnt2d::new(u, v));
        }
        out
    }
}

/// Analytic projection of a 3D point onto an analytic surface
/// (plane/cylinder/cone/sphere/torus) — the OCCT `ProjLib` inverse. Returns
/// `None` for non-analytic surfaces (B-spline, …), where the caller must fall
/// back to a numeric inversion.
///
/// Unlike a Newton iterate seeded from the surface range centre, this is the
/// exact inverse of `d0` even on unbounded-v surfaces (a cylinder/cone point far
/// from v=0), where Newton from the centre cannot converge.
pub fn project_point_on_surface(s: &dyn Surface, p: &GpPnt) -> Option<GpPnt2d> {
    let kind = classify_surface_kind(s);
    let proj = Projector::from_surface(s, kind)?;
    Some(proj.project(p))
}

/// Continuous projection of a 3D curve onto an analytic surface — the OCCT
/// `ProjLib` `projected curve` route. Each sample is projected and then
/// unwrapped across the `u`/`v` seam so the polyline stays on the *same sheet*
/// of a self-intersecting surface (a spindle torus), which a point-wise
/// projection cannot guarantee. Restricted to the torus: for the other analytic
/// surfaces the point-wise grid search is already exact (their `d0` inverse has
/// no sheet ambiguity), so the caller keeps that path. Returns `None` otherwise.
/// `ShapeConstruct_ProjectCurveOnSurface.cxx:71`.
const THE_NCONTROL: usize = 23;

/// `generateCurvePoints` count (`cxx:531-556`).
fn generate_curve_point_count(curve: &dyn Curve, first: f64, last: f64) -> usize {
    let mut n = THE_NCONTROL;
    if curve.bspline_poles().is_some() {
        let knots = curve.parameter_intervals(0);
        let mut used = 0usize;
        for w in knots.windows(2) {
            if w[1] > first && w[0] < last {
                used += 1;
            }
        }
        let deg = curve.nurbs_degree().unwrap_or(1);
        let min_pnt = used * (deg + 1);
        while n < min_pnt {
            n += THE_NCONTROL - 1;
        }
    }
    n.max(2)
}

/// `ShapeConstruct_ProjectCurveOnSurface::Perform` after the plane-only
/// `projectAnalytic` arm (`cxx:708-774`): sample `generateCurvePoints` then
/// project. Isoline shortcuts belong to `BOPTools` MakePCurve, not this path.
pub fn project_curve_on_surface_perform(
    curve: &dyn Curve,
    surf: &dyn Surface,
    first: f64,
    last: f64,
) -> Option<Arc<dyn Curve2d>> {
    if !first.is_finite() || !last.is_finite() || last - first < 1e-15 {
        return None;
    }
    let n = generate_curve_point_count(curve, first, last);
    let t_vals: Vec<f64> = (0..n)
        .map(|i| {
            if i == 0 {
                first
            } else if i + 1 == n {
                last
            } else {
                first + (last - first) * (i as f64) / ((n - 1) as f64)
            }
        })
        .collect();
    let pts = if let Some(pts) = project_curve_on_surface(surf, curve, &t_vals) {
        pts
    } else {
        let kind = classify_surface_kind(surf);
        if let Some(proj) = Projector::from_surface(surf, kind) {
            proj.project_seq(curve, &t_vals)
        } else {
            t_vals
                .iter()
                .map(|&t| {
                    let p = curve.d0(t);
                    let (u, v) = crate::brep_surface::surface_closest_params(surf, &p, 32, 32);
                    GpPnt2d::new(u, v)
                })
                .collect()
        }
    };
    if pts.len() < 2 {
        return None;
    }
    bspline_from_samples(&pts, first, last, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
}

pub fn project_curve_on_surface(
    s: &dyn Surface,
    curve: &dyn Curve,
    t_vals: &[f64],
) -> Option<Vec<GpPnt2d>> {
    let kind = classify_surface_kind(s);
    if kind != SurfaceKind::Torus {
        return None;
    }
    let proj = Projector::from_surface(s, kind)?;
    Some(proj.project_seq(curve, t_vals))
}

/// Recover the plane's location and in-plane axes from its surface `d0`.
pub(super) fn plane_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let o = surf.d0(0.0, 0.0);
    let x = GpVec::from_pnts(&o, &surf.d0(1.0, 0.0));
    let y = GpVec::from_pnts(&o, &surf.d0(0.0, 1.0));
    if x.square_magnitude() < 1e-24 || y.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y))
}

/// Recover a cylinder's axis point, unit axes and radius from its surface `d0`.
pub(super) fn cylinder_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
    let a = surf.d0(0.0, 0.0);
    let b = surf.d0(PI, 0.0);
    let c = surf.d0(FRAC_PI_2, 0.0);
    let o = mid(&a, &b);
    let r = 0.5 * a.distance(&b);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::new((a.x() - o.x()) / r, (a.y() - o.y()) / r, (a.z() - o.z()) / r);
    let y = GpVec::new((c.x() - o.x()) / r, (c.y() - o.y()) / r, (c.z() - o.z()) / r);
    let z = GpVec::from_pnts(&a, &surf.d0(0.0, 1.0));
    if z.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y, z))
}

/// Orthonormal radial frame of a cone: `x` is the `u = 0` radial direction
/// (the generatrix direction projected perpendicular to the axis).
pub(super) fn cone_frame(s: &dyn Surface, _apex: &GpPnt, z: &GpVec) -> Option<(GpVec, GpVec, GpVec)> {
    let dz = GpVec::from_pnts(&s.d0(0.0, 0.0), &s.d0(0.0, 1.0));
    let x = dz.subtracted(&z.multiplied_scalar(dz.dot(z))).normalized();
    if x.square_magnitude() < 1e-20 {
        return None;
    }
    let y = z.crossed(&x).normalized();
    Some((x, y, *z))
}

/// The cone's placement point (the v=0 ring plane, `ElSLib::ConeValue` base).
///
/// `s.d0(0,0)` lies on the v=0 ring: `placement + RefRadius·x`. Subtracting
/// its radial component about the axis (from the apex) recovers the placement
/// point on the axis.
pub(super) fn cone_placement(s: &dyn Surface, apex: &GpPnt, z: &GpVec) -> Option<GpPnt> {
    let p0 = s.d0(0.0, 0.0);
    let d = GpVec::from_pnts(apex, &p0);
    let radial = d.subtracted(&z.multiplied_scalar(d.dot(z)));
    Some(GpPnt::from_xyz(&p0.coord.subtracted(&radial.coord)))
}

// ---------------------------------------------------------------------------
// Analytic isoparametric pcurves
// ---------------------------------------------------------------------------

/// A straight edge (unbounded parameter range).
pub(super) fn is_line(curve: &dyn Curve) -> bool {
    !curve.first_parameter().is_finite() || !curve.last_parameter().is_finite()
}

/// A circle edge (periodic with period 2π).
pub(super) fn is_circle(curve: &dyn Curve) -> bool {
    curve.is_periodic() && (curve.period() - 2.0 * PI).abs() < 1e-9
}

/// Fit a clamped B-spline through `pts` with the exact parameter range [a, b].
/// Degree 1 is a polyline through the sample points (endpoints interpolated).
pub(super) fn bspline_from_samples(pts: &[GpPnt2d], a: f64, b: f64, degree: usize) -> Result<Geom2dBSplineCurve, String> {
    if pts.len() < 2 {
        return Err("pcurve_full: need at least 2 points".into());
    }
    if degree >= pts.len() {
        return Err("pcurve_full: degree must be below pole count".into());
    }
    let mut knots = occt_core::bspl::knots::build_uniform_knots(pts.len(), degree);
    for k in knots.iter_mut() {
        *k = a + (b - a) * *k;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots, degree).map_err(|e| e.to_string())
}

/// Build a `Geom2dLine` (unit-speed) or a degree-1 B-spline (otherwise) through
/// the first and last projected points over the parameter range `[a, b]`.
pub(super) fn line_from_pts(pts: &[GpPnt2d], a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    let q0 = pts[0];
    let q1 = pts[pts.len() - 1];
    let d = GpVec2d::new(q1.x() - q0.x(), q1.y() - q0.y());
    let mag = d.magnitude();
    if mag < 1e-15 {
        return None;
    }
    let dir = GpDir2d::from_vec2d(&d).ok()?;
    let len = b - a;
    if (mag - len).abs() <= 1e-6 * len.max(1.0) {
        let loc = GpPnt2d::new(q0.x() - a * dir.x, q0.y() - a * dir.y);
        Some(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir))))
    } else {
        let bs = bspline_from_samples(pts, a, b, 1).ok()?;
        Some(Arc::new(bs))
    }
}

/// Analytic pcurve of a line edge on an analytic revolution surface: a
/// generatrix projects to a `u`-isoline (varying `v`) or `v`-isoline.
/// `None` when the projection is not an isoparametric line (the caller then
/// samples instead).
pub(super) fn iso_line_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    if b - a < 1e-15 {
        return None;
    }
    let t_vals = [a, 0.5 * (a + b), b];
    let pts = proj.project_seq(curve, &t_vals);
    let tol = 1e-6 * (1.0 + pts[0].x().abs() + pts[0].y().abs());
    let u_const = pts.iter().all(|q| (q.x() - pts[0].x()).abs() < tol);
    let v_const = pts.iter().all(|q| (q.y() - pts[0].y()).abs() < tol);
    if !u_const && !v_const {
        return None;
    }
    line_from_pts(&pts, a, b)
}

/// Analytic pcurve of a circle edge on an analytic revolution surface: a
/// latitude (constant `v`) or meridian (constant `u`) projects to a straight
/// UV line; a U-periodic full circle unwraps across the seam. `None` when the
/// projection is not an isoparametric line.
pub(super) fn iso_circle_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    if b - a < 1e-15 {
        return None;
    }
    let n = 9;
    let t_vals: Vec<f64> = (0..n).map(|i| a + (b - a) * i as f64 / (n - 1) as f64).collect();
    let pts = proj.project_seq(curve, &t_vals);
    let tol = 1e-6 * (1.0 + pts[0].x().abs() + pts[0].y().abs());
    let u_const = pts.iter().all(|q| (q.x() - pts[0].x()).abs() < tol);
    let v_const = pts.iter().all(|q| (q.y() - pts[0].y()).abs() < tol);
    if u_const == v_const {
        // Degenerate (both) or non-isoparametric (neither): let the sampling
        // path handle it.
        return None;
    }
    line_from_pts(&pts, a, b)
}

// ---------------------------------------------------------------------------
// Main entry points
// ---------------------------------------------------------------------------

/// A 2D curve linearly reparameterized from `[p0, p1]` onto `[a, b]`: a parameter
/// `t ∈ [a, b]` maps to `inner.d0(p0 + (t − a)·(p1 − p0)/(b − a))`. Aligns a
/// STEP-imported pcurve (parameterized over its natural range) with the edge's
/// 3D range.
#[derive(Clone)]
pub struct ReparamCurve2d {
    pub(super) inner: Arc<dyn Curve2d>,
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) p0: f64,
    pub(super) p1: f64,
}

impl ReparamCurve2d {
    pub fn new(inner: Arc<dyn Curve2d>, a: f64, b: f64, p0: f64, p1: f64) -> Self {
        Self { inner, a, b, p0, p1 }
    }

    pub(super) fn scale(&self) -> f64 {
        if (self.b - self.a).abs() < 1e-30 {
            1.0
        } else {
            (self.p1 - self.p0) / (self.b - self.a)
        }
    }

    pub(super) fn map(&self, t: f64) -> f64 {
        self.p0 + (t - self.a) * self.scale()
    }
}

impl Curve2d for ReparamCurve2d {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.inner.d0(self.map(u))
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p, d) = self.inner.d1(self.map(u));
        let s = self.scale();
        (p, GpVec2d::new(d.x() * s, d.y() * s))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1, d2) = self.inner.d2(self.map(u));
        let s = self.scale();
        (
            p,
            GpVec2d::new(d1.x() * s, d1.y() * s),
            GpVec2d::new(d2.x() * s * s, d2.y() * s * s),
        )
    }
    fn first_parameter(&self) -> f64 {
        self.a
    }
    fn last_parameter(&self) -> f64 {
        self.b
    }
    fn continuity(&self) -> u8 {
        self.inner.continuity()
    }
    fn transform(&mut self, t: &GpTrsf2d) {
        let mut inner = self.inner.clone_dyn();
        inner.transform(t);
        self.inner = Arc::from(inner);
    }
    fn reverse(&mut self) {
        let mut inner = self.inner.clone_dyn();
        inner.reverse();
        self.inner = Arc::from(inner);
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
}

/// Wrap a 2D curve in a linear reparameterization from `[p0, p1]` onto `[a, b]`.
pub fn reparam_curve2d(
    inner: Arc<dyn Curve2d>,
    a: f64,
    b: f64,
    p0: f64,
    p1: f64,
) -> Arc<dyn Curve2d> {
    Arc::new(ReparamCurve2d::new(inner, a, b, p0, p1))
}

/// Parameter of the 2D curve closest to `uv` (the 1-D projection of a surface
/// point onto the pcurve). A bounded curve is searched on `[first, last]`; an
/// unbounded line is projected in closed form.
pub fn project_uv_on_curve2d(curve: &dyn Curve2d, uv: GpPnt2d) -> f64 {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    if f.is_finite() && l.is_finite() {
        let mut best = f;
        let mut best_d = f64::INFINITY;
        for i in 0..=64 {
            let u = f + (l - f) * i as f64 / 64.0;
            let d = curve.d0(u).distance(&uv);
            if d < best_d {
                best_d = d;
                best = u;
            }
        }
        let mut step = (l - f) / 64.0;
        for _ in 0..32 {
            step *= 0.5;
            let mut u = best;
            let mut dd = best_d;
            for &du in &[-step, step] {
                let nu = (u + du).clamp(f, l);
                let d = curve.d0(nu).distance(&uv);
                if d < dd {
                    dd = d;
                    u = nu;
                }
            }
            best = u;
            best_d = dd;
        }
        best
    } else {
        let p0 = curve.d0(0.0);
        let dir = curve.d0(1.0).coord.subtracted(&p0.coord);
        let len2 = dir.square_modulus();
        if len2 < 1e-30 {
            return 0.0;
        }
        uv.coord.subtracted(&p0.coord).dot(&dir) / len2
    }
}
