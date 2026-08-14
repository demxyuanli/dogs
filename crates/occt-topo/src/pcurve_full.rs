//! Full edge→face UV curve (pcurve) construction — the Phase 16b superset.
//!
//! Extends [`crate::pcurve::make_pcurve_on_face`] from the plane/cylinder pair
//! to every analytic surface (plane, cylinder, cone, sphere, torus) and to
//! general (B-spline) faces:
//!
//! * **Plane** — the Phase 16a implementation is reused verbatim, so
//!   [`make_pcurve_full`] and [`crate::pcurve::make_pcurve_on_face`] agree
//!   exactly: a line edge becomes a `Geom2dLine`, a circle edge a
//!   `Geom2dCircle`, anything else a sampled degree-1 B-spline.
//! * **Cylinder / cone / sphere / torus** — the surface's own parameterization
//!   is recovered from sampled geometry invariants and used to project the edge
//!   exactly.  A generatrix (line edge) or latitude/meridian (circle edge) maps
//!   to a `u`- or `v`-isoline, returned as a `Geom2dLine` (unit-speed) or a
//!   degree-1 B-spline otherwise.  Non-isoparametric edges fall through to the
//!   sampling path.
//! * **Anything else (B-spline faces)** — the edge is sampled at 48 parameters
//!   and the projected `(u, v)` points are fitted with a clamped degree-1
//!   B-spline whose parameter range matches the edge range.
//!
//! Source: `BOPTools_AlgoTools2D::MakePCurveOnFace` / `AdjustPCurveOnSurf`
//! (TKBO) with the `ProjLib` analytic projection rules for the revolution
//! surfaces.

use std::f64::consts::{FRAC_PI_2, PI};
use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt, GpPnt2d, GpTrsf2d, GpVec, GpVec2d};
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::{Geom2dBSplineCurve, Geom2dLine};

use crate::brep_surface::{
    classify_surface, edge_pcurve_on_face, face_uv_bounds, sphere_center, SurfaceKind,
};
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

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
fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
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

fn mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

fn pnt_add_vec(p: GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::from_xyz(&p.coord.added(&v.coord))
}

/// Closest point between two (possibly skew) lines; `None` when parallel or
/// too far apart.
fn closest_point_lines(p1: GpPnt, d1: GpVec, p2: GpPnt, d2: GpVec) -> Option<GpPnt> {
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
enum Projector {
    Plane { o: GpPnt, x: GpVec, y: GpVec },
    Cylinder { o: GpPnt, x: GpVec, y: GpVec, z: GpVec },
    Cone { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, alpha: f64 },
    Sphere { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r: f64 },
    Torus { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r_big: f64, r_small: f64 },
}

impl Projector {
    fn from_surface(s: &dyn Surface, kind: SurfaceKind) -> Option<Projector> {
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
    fn is_u_periodic(&self) -> bool {
        matches!(
            self,
            Projector::Cylinder { .. } | Projector::Cone { .. } | Projector::Sphere { .. } | Projector::Torus { .. }
        )
    }

    /// Whether the surface's `v` parameter is periodic (torus tube direction).
    fn is_v_periodic(&self) -> bool {
        matches!(self, Projector::Torus { .. })
    }

    /// Project a single 3D point onto the face `(u, v)` domain. `u` is wrapped
    /// into the principal angle range `(-π, π]`; use [`Self::project_seq`] to
    /// unwrap a continuous sequence across the seam.
    fn project(&self, p: &GpPnt) -> GpPnt2d {
        match self {
            Projector::Plane { o, x, y } => {
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
    fn project_seq(&self, curve: &dyn Curve, t_vals: &[f64]) -> Vec<GpPnt2d> {
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
fn plane_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let o = surf.d0(0.0, 0.0);
    let x = GpVec::from_pnts(&o, &surf.d0(1.0, 0.0));
    let y = GpVec::from_pnts(&o, &surf.d0(0.0, 1.0));
    if x.square_magnitude() < 1e-24 || y.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y))
}

/// Recover a cylinder's axis point, unit axes and radius from its surface `d0`.
fn cylinder_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
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
fn cone_frame(s: &dyn Surface, _apex: &GpPnt, z: &GpVec) -> Option<(GpVec, GpVec, GpVec)> {
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
fn cone_placement(s: &dyn Surface, apex: &GpPnt, z: &GpVec) -> Option<GpPnt> {
    let p0 = s.d0(0.0, 0.0);
    let d = GpVec::from_pnts(apex, &p0);
    let radial = d.subtracted(&z.multiplied_scalar(d.dot(z)));
    Some(GpPnt::from_xyz(&p0.coord.subtracted(&radial.coord)))
}

// ---------------------------------------------------------------------------
// Analytic isoparametric pcurves
// ---------------------------------------------------------------------------

/// A straight edge (unbounded parameter range).
fn is_line(curve: &dyn Curve) -> bool {
    !curve.first_parameter().is_finite() || !curve.last_parameter().is_finite()
}

/// A circle edge (periodic with period 2π).
fn is_circle(curve: &dyn Curve) -> bool {
    curve.is_periodic() && (curve.period() - 2.0 * PI).abs() < 1e-9
}

/// Fit a clamped B-spline through `pts` with the exact parameter range [a, b].
/// Degree 1 is a polyline through the sample points (endpoints interpolated).
fn bspline_from_samples(pts: &[GpPnt2d], a: f64, b: f64, degree: usize) -> Result<Geom2dBSplineCurve, String> {
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
fn line_from_pts(pts: &[GpPnt2d], a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
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
fn iso_line_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
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
fn iso_circle_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
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

/// Construct the pcurve of `edge` on `face` — the edge's 3D curve projected
/// into the face surface's `(u, v)` parameter domain, with analytic cases for
/// every analytic surface and a sampling fallback for B-spline faces.
///
/// The returned pcurve is parameterized over the edge's range `[a, b]` so
/// `d0(a)` / `d0(b)` land on the projected endpoints.
pub fn make_pcurve_full(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
    let Some(curve) = GeometryRegistry::global().edge_curve(&edge.0) else {
        return Err("make_pcurve_full: edge has no 3D curve".into());
    };
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Err("make_pcurve_full: edge range is empty or unbounded".into());
    }
    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
        return Err("make_pcurve_full: face has no surface".into());
    };
    let kind = classify_surface_kind(surf.as_ref());

    // Plane faces: reuse the Phase-16a implementation verbatim, so
    // make_pcurve_full and pcurve::make_pcurve_on_face agree exactly.
    if kind == SurfaceKind::Plane {
        return crate::pcurve::make_pcurve_on_face(edge, face);
    }

    // Analytic revolution surfaces: exact isoparametric pcurves.
    if let Some(proj) = Projector::from_surface(surf.as_ref(), kind) {
        if is_line(curve.as_ref()) {
            if let Some(pc) = iso_line_pcurve(curve.as_ref(), &proj, a, b) {
                return Ok(pc);
            }
        }
        if is_circle(curve.as_ref()) {
            if let Some(pc) = iso_circle_pcurve(curve.as_ref(), &proj, a, b) {
                return Ok(pc);
            }
        }
    }

    // General B-spline surface (or a non-isoparametric edge): sample the
    // projection and fit a clamped degree-1 B-spline over the edge range.
    let pts = edge_pcurve_on_face(edge, face, 48);
    if pts.len() < 2 {
        return Err("make_pcurve_full: too few projected samples".into());
    }
    let bs = bspline_from_samples(&pts, a, b, 1)?;
    Ok(Arc::new(bs))
}

/// Direction correction of a pcurve: `+1` when the pcurve parameter follows the
/// edge's 3D direction, `-1` when it must be reversed.
///
/// The surface-mapped pcurve tangent `∂S/∂u·u' + ∂S/∂v·v'` is compared against
/// the edge's 3D tangent (flipped for a REVERSED edge). Surfaces whose `d1`
/// returns zero vectors fall back to comparing endpoint displacements.
pub fn pc_curve_orientation(edge: &Edge, face: &Face, curve: &dyn Curve2d) -> f64 {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return 1.0;
    }
    let tm = 0.5 * (a + b);
    let (Some(c3d), Some(surf)) = (
        GeometryRegistry::global().edge_curve(&edge.0),
        GeometryRegistry::global().face_surface(&face.0),
    ) else {
        return 1.0;
    };

    let (_, mut d3d) = c3d.d1(tm);
    if edge.orientation().is_reversed() {
        d3d = d3d.reversed();
    }
    let q = curve.d0(tm);
    let (_, du, dv) = surf.d1(q.x(), q.y());
    let d2d = curve.d1(tm).1;
    let mapped = du.multiplied_scalar(d2d.x()).add(&dv.multiplied_scalar(d2d.y()));
    let dot = d3d.dot(&mapped);
    if dot.abs() < 1e-12 {
        // Degenerate surface tangent (analytic `d1` may be zero): compare the
        // endpoint displacement of the edge curve against the UV-mapped one.
        let p0 = c3d.d0(a);
        let p1 = c3d.d0(b);
        let mut disp3 = GpVec::from_pnts(&p0, &p1);
        if edge.orientation().is_reversed() {
            disp3 = disp3.reversed();
        }
        let q0 = curve.d0(a);
        let q1 = curve.d0(b);
        let suv0 = surf.d0(q0.x(), q0.y());
        let suv1 = surf.d0(q1.x(), q1.y());
        let disp_uv = GpVec::from_pnts(&suv0, &suv1);
        return if disp3.dot(&disp_uv) >= 0.0 { 1.0 } else { -1.0 };
    }
    if dot > 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// Wrap `x` into the half-open interval `[lo, hi)` by adding/subtracting the
/// period `hi − lo` (a no-op when the interval is degenerate).
fn wrap_into(x: f64, lo: f64, hi: f64) -> f64 {
    let p = hi - lo;
    if p <= 0.0 || !p.is_finite() {
        return x;
    }
    let mut r = (x - lo) % p;
    if r < 0.0 {
        r += p;
    }
    lo + r
}

/// The `u`-period of a periodic surface (its `u`-range span).
fn u_period(s: &dyn Surface) -> f64 {
    let (u0, u1) = s.u_range();
    if u0.is_finite() && u1.is_finite() {
        u1 - u0
    } else {
        2.0 * PI
    }
}

/// The `v`-period of a periodic surface (its `v`-range span).
fn v_period(s: &dyn Surface) -> f64 {
    let (v0, v1) = s.v_range();
    if v0.is_finite() && v1.is_finite() {
        v1 - v0
    } else {
        2.0 * PI
    }
}

/// Bring a pcurve inside the face's UV bounds.
///
/// For periodic surface dimensions (the `u` of cylinder/cone/sphere/torus and
/// the `v` of a torus) the curve is translated by whole periods so its midpoint
/// lies in range. If part of the curve still leaves a non-periodic bounded
/// dimension, it is trimmed to the in-bounds parameter subrange.
pub fn trim_pcurve_to_face(curve: &Arc<dyn Curve2d>, face: &Face, _tol: f64) -> Result<Arc<dyn Curve2d>, String> {
    let (umin, umax, vmin, vmax) = face_uv_bounds(face);
    let surf = GeometryRegistry::global()
        .face_surface(&face.0)
        .ok_or("trim_pcurve_to_face: face has no surface")?;
    let u_periodic = surf.is_u_periodic();
    let v_periodic = surf.is_v_periodic();
    let up = if u_periodic { u_period(surf.as_ref()) } else { 0.0 };
    let vp = if v_periodic { v_period(surf.as_ref()) } else { 0.0 };

    let (mut c0, mut c1) = (curve.first_parameter(), curve.last_parameter());
    if !c0.is_finite() || !c1.is_finite() || c1 - c0 < 1e-15 {
        // Unbounded 2D curve (a Geom2dLine): evaluate a unit window around 0.
        c0 = 0.0;
        c1 = 1.0;
    }

    // Sample the curve and compute the whole-period shift that brings the
    // midpoint inside the periodic bounds.
    let n = 33;
    let t_vals: Vec<f64> = (0..n).map(|i| c0 + (c1 - c0) * i as f64 / (n - 1) as f64).collect();
    let qm = curve.d0(0.5 * (c0 + c1));
    let mut du = 0.0;
    let mut dv = 0.0;
    if u_periodic && up > 0.0 && umin.is_finite() {
        du = wrap_into(qm.x(), umin, umin + up) - qm.x();
    }
    if v_periodic && vp > 0.0 && vmin.is_finite() {
        dv = wrap_into(qm.y(), vmin, vmin + vp) - qm.y();
    }

    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(&GpVec2d::new(du, dv));
    let shifted: Arc<dyn Curve2d> = Arc::from(curve.transformed(&tr));

    let inside = |q: &GpPnt2d| {
        (if umin.is_finite() { q.x() >= umin - 1e-7 } else { true })
            && (if umax.is_finite() { q.x() <= umax + 1e-7 } else { true })
            && (if vmin.is_finite() { q.y() >= vmin - 1e-7 } else { true })
            && (if vmax.is_finite() { q.y() <= vmax + 1e-7 } else { true })
    };

    let shifted_pts: Vec<GpPnt2d> = t_vals.iter().map(|&t| shifted.d0(t)).collect();
    if shifted_pts.iter().all(inside) {
        return Ok(shifted);
    }

    // Some samples still violate a non-periodic bounded dimension: trim to the
    // in-bounds parameter subrange around the midpoint.
    let in_bounds: Vec<bool> = shifted_pts.iter().map(inside).collect();
    let Some(first) = in_bounds.iter().position(|&b| b) else {
        return Err("trim_pcurve_to_face: pcurve lies entirely outside the face UV bounds".into());
    };
    let last = in_bounds.iter().rposition(|&b| b).unwrap();
    let t0 = t_vals[first];
    let t1 = t_vals[last];
    let kept: Vec<GpPnt2d> = shifted_pts[first..=last].to_vec();
    let bs = bspline_from_samples(&kept, t0, t1, 1)?;
    Ok(Arc::new(bs))
}

/// Evaluate a pcurve at `param`, validating that the parameter is inside the
/// curve's range (a `Geom2dLine` accepts any parameter).
pub fn pcurve_point_on_face(curve: &dyn Curve2d, param: f64) -> Result<GpPnt2d, String> {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if a.is_finite() && b.is_finite() && (param < a - 1e-9 || param > b + 1e-9) {
        return Err(format!(
            "pcurve_point_on_face: parameter {param} outside range [{a}, {b}]"
        ));
    }
    Ok(curve.d0(param))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use occt_core::gp::{GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus};
    use occt_geom::{GeomBSplineCurve, GeomCone, GeomCylinder, GeomSphere, GeomTorus};

    fn cyl_face(radius: f64) -> Face {
        let b = TopoBuilder::new();
        let cyl = GpCylinder::new(GpAx3::standard(), radius).unwrap();
        b.make_face(Arc::new(GeomCylinder::new(cyl)), &[])
    }

    #[test]
    fn classify_surface_kinds() {
        let b = TopoBuilder::new();
        let pln = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let s = GeometryRegistry::global().face_surface(&pln.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Plane);

        let fc = cyl_face(1.0);
        let s = GeometryRegistry::global().face_surface(&fc.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Cylinder);

        let sph = GpSphere::new(GpAx3::standard(), 2.0).unwrap();
        let fs = b.make_face(Arc::new(GeomSphere::new(sph)), &[]);
        let s = GeometryRegistry::global().face_surface(&fs.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Sphere);

        let cone = GpCone::new(GpAx3::standard(), 1.0, 0.5f64.atan()).unwrap();
        let fcn = b.make_face(Arc::new(GeomCone::new(cone)), &[]);
        let s = GeometryRegistry::global().face_surface(&fcn.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Cone);

        let tor = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let ft = b.make_face(Arc::new(GeomTorus::new(tor)), &[]);
        let s = GeometryRegistry::global().face_surface(&ft.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Torus);
    }

    #[test]
    fn unit_box_line_edge_is_line_pcurve_with_boundary_endpoints() {
        let b = unit_box();
        let edge = &b.edges[4]; // (0,0,1) → (1,0,1)
        let face = &b.faces[1]; // top (z = 1)
        let pc = make_pcurve_full(edge, face).expect("pcurve");
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::Line);
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 1.0).abs() < 1e-6 && (p1.y() - 0.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn cylinder_generatrix_is_v_line() {
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0));
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::Line);
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        // Generatrix at angle 0: u = atan2(0, 1) = 0, v = axial coordinate.
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 0.0).abs() < 1e-6 && (p1.y() - 1.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn cylinder_top_circle_is_u_line_across_the_seam() {
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let ax2 = GpAx2::new(
            GpPnt::new(0.0, 0.0, 1.0),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let edge = b.make_edge_circle(&ax2, 1.0, 0.0, 2.0 * PI);
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        // The cap circle is a v = 1 isoline in the U-periodic direction.
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::Line);
        let p0 = pc.d0(0.0);
        let pm = pc.d0(PI);
        let p2 = pc.d0(2.0 * PI);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 1.0).abs() < 1e-6, "start {:?}", p0);
        assert!((pm.x() - PI).abs() < 1e-6 && (pm.y() - 1.0).abs() < 1e-6, "mid {:?}", pm);
        assert!((p2.x() - 2.0 * PI).abs() < 1e-6 && (p2.y() - 1.0).abs() < 1e-6, "end {:?}", p2);
    }

    #[test]
    fn general_bspline_edge_on_plane_face_samples_boundary() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.5, 0.4, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let curve = Arc::new(GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap());
        let edge = b.make_edge(curve, 0.0, 1.0);
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::BSpline);
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 1.0).abs() < 1e-6 && (p1.y() - 0.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn agrees_with_make_pcurve_on_face_for_plane_line() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        let full = make_pcurve_full(edge, face).expect("full pcurve");
        let base = crate::pcurve::make_pcurve_on_face(edge, face).expect("base pcurve");
        let (a, z) = GeometryRegistry::global().edge_parameters(&edge.0);
        for i in 0..=8 {
            let t = a + (z - a) * i as f64 / 8.0;
            let q1 = full.d0(t);
            let q2 = base.d0(t);
            assert!((q1.x() - q2.x()).abs() < 1e-9 && (q1.y() - q2.y()).abs() < 1e-9, "t {t}");
        }
    }

    #[test]
    fn trim_shifts_periodic_pcurve_into_face_bounds() {
        let face = cyl_face(1.0);
        // A degree-1 B-spline with u in [6.5, 7.5] at v = 0: u is outside the
        // face's [0, 2π] range.
        let pts = [
            GpPnt2d::new(6.5, 0.0),
            GpPnt2d::new(7.0, 0.0),
            GpPnt2d::new(7.5, 0.0),
        ];
        let bs = bspline_from_samples(&pts, 6.5, 7.5, 1).unwrap();
        let pc: Arc<dyn Curve2d> = Arc::new(bs);
        let trimmed = trim_pcurve_to_face(&pc, &face, 1e-7).expect("trimmed");
        let (umin, umax, vmin, vmax) = face_uv_bounds(&face);
        for i in 0..=16 {
            let t = 6.5 + (7.5 - 6.5) * i as f64 / 16.0;
            let q = trimmed.d0(t);
            assert!(q.x() >= umin - 1e-6 && q.x() <= umax + 1e-6, "u {} at t {}", q.x(), t);
            assert!(q.y() >= vmin - 1e-6 && q.y() <= vmax + 1e-6, "v {} at t {}", q.y(), t);
        }
    }

    #[test]
    fn orientation_matches_edge_direction() {
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0));
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        assert_eq!(pc_curve_orientation(&edge, &face, pc.as_ref()), 1.0);
        // A reversed edge runs against the stored curve direction.
        let mut rev = edge.clone();
        rev.0.set_orientation(crate::abs::Orientation::Reversed);
        assert_eq!(pc_curve_orientation(&rev, &face, pc.as_ref()), -1.0);
    }

    #[test]
    fn pcurve_point_evaluation() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let edge = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        let q = pcurve_point_on_face(pc.as_ref(), 0.5).expect("point");
        // The unit-speed Geom2dLine pcurve maps edge parameter → (u, 0).
        assert!((q.x() - 0.5).abs() < 1e-6 && (q.y() - 0.0).abs() < 1e-6, "q {:?}", q);

        // A finite-range B-spline pcurve rejects out-of-range parameters.
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.5, 0.4, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let curve = Arc::new(GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap());
        let be = b.make_edge(curve, 0.0, 1.0);
        let bpc = make_pcurve_full(&be, &face).expect("bspline pcurve");
        assert!(pcurve_point_on_face(bpc.as_ref(), 5.0).is_err(), "out of range must fail");
    }

    #[test]
    fn non_isoparametric_circle_on_cylinder_falls_back_to_sampling() {
        // A circle in the X-Z plane (not a cap circle, and only touching the
        // cylinder at the two axis points) is not an isoparametric curve, so
        // the pcurve is built by the sampling path. Its endpoints still land on
        // the projected surface parameters.
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let ax2 = GpAx2::new(
            GpPnt::zero(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let edge = b.make_edge_circle(&ax2, 1.0, 0.0, PI);
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(PI);
        // Endpoint (1, 0, 0) → (u=0, v=0); endpoint (−1, 0, 0) → (u=π, v=0).
        assert!((p0.x() - 0.0).abs() < 1e-5 && (p0.y() - 0.0).abs() < 1e-5, "start {:?}", p0);
        assert!((p1.x() - PI).abs() < 1e-5 && (p1.y() - 0.0).abs() < 1e-5, "end {:?}", p1);
    }
}
