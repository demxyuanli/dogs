//! [`Projector`]: the analytic per-type UV projection of a 3D point or curve, with the [`project_point_on_surface`] entry point.

use super::*;

/// Exact UV projection of a 3D point onto an analytic surface, recovered from
/// the surface's own parameterization so the results agree with
/// `edge_pcurve_on_face`.
pub(crate) enum Projector {
    Plane { o: GpPnt, x: GpVec, y: GpVec },
    Cylinder { o: GpPnt, x: GpVec, y: GpVec, z: GpVec },
    Cone { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, alpha: f64 },
    Sphere { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r: f64 },
    Torus { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r_big: f64, r_small: f64 },
}

impl Projector {
    pub(crate) fn from_surface(s: &dyn Surface, kind: SurfaceKind) -> Option<Projector> {
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
    pub(crate) fn is_u_periodic(&self) -> bool {
        matches!(
            self,
            Projector::Cylinder { .. } | Projector::Cone { .. } | Projector::Sphere { .. } | Projector::Torus { .. }
        )
    }

    /// Whether the surface's `v` parameter is periodic (torus tube direction).
    pub(crate) fn is_v_periodic(&self) -> bool {
        matches!(self, Projector::Torus { .. })
    }

    /// Project a single 3D point onto the face `(u, v)` domain. `u` is wrapped
    /// into the principal angle range `(-π, π]`; use [`Self::project_seq`] to
    /// unwrap a continuous sequence across the seam.
    pub(crate) fn project(&self, p: &GpPnt) -> GpPnt2d {
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
    pub(crate) fn project_seq(&self, curve: &dyn Curve, t_vals: &[f64]) -> Vec<GpPnt2d> {
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
pub fn surface_project_point(s: &dyn Surface, p: &GpPnt) -> Option<GpPnt2d> {
    let kind = classify_surface_kind(s);
    let proj = Projector::from_surface(s, kind)?;
    Some(proj.project(p))
}

