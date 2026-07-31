//! Surface parameterization utilities: isoparametric curves, reparameterization,
//! normal grids, curvature, periodic adjustment.
//! Source: `GeomAdaptor_Surface.hxx`, `GeomAPI_ProjectPointOnSurf.hxx`.

use std::sync::Arc;
use crate::curve::Curve;
use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpVec, GpTrsf};

/// An isoparametric curve of a surface with its metadata.
///
/// `direction` is `'u'` when the curve varies the `u` parameter (fixed `v`) and
/// `'v'` when it varies `v` (fixed `u`); `fixed` holds the constant parameter.
#[derive(Clone)]
pub struct IsoparametricCurve {
    pub curve: Arc<dyn Curve>,
    pub fixed: f64,
    pub direction: char,
}

impl Curve for IsoparametricCurve {
    fn d0(&self, u: f64) -> GpPnt { self.curve.d0(u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { self.curve.d1(u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { self.curve.d2(u) }
    fn first_parameter(&self) -> f64 { self.curve.first_parameter() }
    fn last_parameter(&self) -> f64 { self.curve.last_parameter() }
    fn is_periodic(&self) -> bool { self.curve.is_periodic() }
    fn period(&self) -> f64 { self.curve.period() }
    fn continuity(&self) -> u8 { self.curve.continuity() }
    // ponytail: adapter over an immutable Arc; transforms are no-ops like
    // GeomTrimmedCurve. Apply transforms to the underlying surface instead.
    fn transform(&mut self, _t: &GpTrsf) {}
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}

/// Surface-calling isoparametric curve adapter. `c(t) = s(t, v)` (direction
/// `'u'`) or `c(t) = s(u, t)` (direction `'v'`), evaluated directly from the
/// surface rather than a pre-computed curve.
#[derive(Clone)]
struct IsoCurve {
    surface: Arc<dyn Surface>,
    fixed: f64,
    direction: char,
}

impl IsoCurve {
    fn new(surface: Arc<dyn Surface>, fixed: f64, direction: char) -> Self {
        Self { surface, fixed, direction }
    }

    fn eval(&self, t: f64) -> GpPnt {
        match self.direction {
            'u' => self.surface.d0(t, self.fixed),
            _ => self.surface.d0(self.fixed, t),
        }
    }

    fn range(&self) -> (f64, f64) {
        match self.direction {
            'u' => self.surface.u_range(),
            _ => self.surface.v_range(),
        }
    }

    fn periodic(&self) -> bool {
        match self.direction {
            'u' => self.surface.is_u_periodic(),
            _ => self.surface.is_v_periodic(),
        }
    }
}

impl Curve for IsoCurve {
    fn d0(&self, t: f64) -> GpPnt { self.eval(t) }

    fn d1(&self, t: f64) -> (GpPnt, GpVec) {
        let (p, d) = match self.direction {
            'u' => {
                let (p, du, _) = self.surface.d1(t, self.fixed);
                (p, du)
            }
            _ => {
                let (p, _, dv) = self.surface.d1(self.fixed, t);
                (p, dv)
            }
        };
        if d.xyz().square_modulus() > 1e-30 {
            return (p, d);
        }
        // Fallback: central finite differences (some surfaces return zero d1).
        let h = 1e-6;
        let p1 = self.eval(t + h);
        let p2 = self.eval(t - h);
        (
            p,
            GpVec::new(
                (p1.x() - p2.x()) / (2.0 * h),
                (p1.y() - p2.y()) / (2.0 * h),
                (p1.z() - p2.z()) / (2.0 * h),
            ),
        )
    }

    fn d2(&self, t: f64) -> (GpPnt, GpVec, GpVec) {
        // ponytail: finite-difference second derivative; surfaces expose no d2.
        let h = 1e-6;
        let p = self.eval(t);
        let (_, d1) = self.d1(t);
        let pp = self.eval(t + h);
        let pm = self.eval(t - h);
        let d2 = GpVec::new(
            (pp.x() - 2.0 * p.x() + pm.x()) / (h * h),
            (pp.y() - 2.0 * p.y() + pm.y()) / (h * h),
            (pp.z() - 2.0 * p.z() + pm.z()) / (h * h),
        );
        (p, d1, d2)
    }

    fn first_parameter(&self) -> f64 { self.range().0 }
    fn last_parameter(&self) -> f64 { self.range().1 }
    fn is_periodic(&self) -> bool { self.periodic() }
    fn period(&self) -> f64 {
        if self.periodic() {
            let (a, b) = self.range();
            b - a
        } else {
            0.0
        }
    }
    fn continuity(&self) -> u8 { self.surface.continuity() }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}

/// Curve `c(u) = s(u, v)` over the surface's `u` range.
pub fn iso_u(s: &dyn Surface, v: f64) -> Arc<dyn Curve> {
    Arc::new(IsoparametricCurve {
        curve: Arc::new(IsoCurve::new(Arc::from(s.clone_dyn()), v, 'u')),
        fixed: v,
        direction: 'u',
    })
}

/// Curve `c(v) = s(u, v)` over the surface's `v` range.
pub fn iso_v(s: &dyn Surface, u: f64) -> Arc<dyn Curve> {
    Arc::new(IsoparametricCurve {
        curve: Arc::new(IsoCurve::new(Arc::from(s.clone_dyn()), u, 'v')),
        fixed: u,
        direction: 'v',
    })
}

/// Grid of isoparametric curves: `nu` curves of constant `u` (varying `v`) plus
/// `nv` curves of constant `v` (varying `u`), sampled across the parameter
/// ranges. Unbounded ranges clamp to `[-1, 1]`.
pub fn iso_curves_grid(s: &dyn Surface, nu: usize, nv: usize) -> (Vec<Arc<dyn Curve>>, Vec<Arc<dyn Curve>>) {
    let (u0, u1) = finite_range(s.u_range());
    let (v0, v1) = finite_range(s.v_range());
    let u_const: Vec<Arc<dyn Curve>> = (0..nu.max(1))
        .map(|i| {
            let u = if nu <= 1 { 0.5 * (u0 + u1) } else { u0 + (u1 - u0) * i as f64 / (nu - 1) as f64 };
            iso_v(s, u)
        })
        .collect();
    let v_const: Vec<Arc<dyn Curve>> = (0..nv.max(1))
        .map(|j| {
            let v = if nv <= 1 { 0.5 * (v0 + v1) } else { v0 + (v1 - v0) * j as f64 / (nv - 1) as f64 };
            iso_u(s, v)
        })
        .collect();
    (u_const, v_const)
}

/// Surface adapter that maps a new parameter range `[new_u0, new_u1] ×
/// [new_v0, new_v1]` linearly onto the original surface's range. The reported
/// range of the returned surface is the new one.
#[derive(Clone)]
struct ReparamSurface {
    surface: Arc<dyn Surface>,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    new_u0: f64,
    new_u1: f64,
    new_v0: f64,
    new_v1: f64,
}

impl ReparamSurface {
    fn map_u(&self, u: f64) -> f64 {
        self.u0 + (u - self.new_u0) * (self.u1 - self.u0) / (self.new_u1 - self.new_u0)
    }
    fn map_v(&self, v: f64) -> f64 {
        self.v0 + (v - self.new_v0) * (self.v1 - self.v0) / (self.new_v1 - self.new_v0)
    }
}

impl Surface for ReparamSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.surface.d0(self.map_u(u), self.map_v(v))
    }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let su = (self.u1 - self.u0) / (self.new_u1 - self.new_u0);
        let sv = (self.v1 - self.v0) / (self.new_v1 - self.new_v0);
        let (p, du, dv) = self.surface.d1(self.map_u(u), self.map_v(v));
        (p, du.multiplied_scalar(su), dv.multiplied_scalar(sv))
    }
    fn u_range(&self) -> (f64, f64) { (self.new_u0, self.new_u1) }
    fn v_range(&self) -> (f64, f64) { (self.new_v0, self.new_v1) }
    fn is_u_periodic(&self) -> bool { self.surface.is_u_periodic() }
    fn is_v_periodic(&self) -> bool { self.surface.is_v_periodic() }
    fn continuity(&self) -> u8 { self.surface.continuity() }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}

/// Reparameterize `s` so the returned surface's range is `[new_u0, new_u1] ×
/// [new_v0, new_v1]` while preserving geometry: new `(u', v')` maps onto the
/// original range linearly.
pub fn reparameterize_surface(
    s: &dyn Surface,
    new_u0: f64,
    new_u1: f64,
    new_v0: f64,
    new_v1: f64,
) -> Arc<dyn Surface> {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    Arc::new(ReparamSurface {
        surface: Arc::from(s.clone_dyn()),
        u0,
        u1,
        v0,
        v1,
        new_u0,
        new_u1,
        new_v0,
        new_v1,
    })
}

/// Unit normals of `s` on a `nu × nv` grid, robust to surfaces whose `d1`
/// returns zero vectors (fallback: finite differences of `d0`).
pub fn surface_normal_grid(s: &dyn Surface, nu: usize, nv: usize) -> Vec<Vec<GpVec>> {
    let (u0, u1) = finite_range(s.u_range());
    let (v0, v1) = finite_range(s.v_range());
    let nu = nu.max(2);
    let nv = nv.max(2);
    (0..nu)
        .map(|i| {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            (0..nv)
                .map(|j| {
                    let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
                    unit_normal(s, u, v)
                })
                .collect()
        })
        .collect()
}

/// Principal curvatures `(k_min, k_max)` at `(u, v)`.
///
/// Approximation: finite-difference shape operator. Unit normal derivatives
/// `dN/du`, `dN/dv` and surface tangents are obtained numerically, and the
/// Weingarten matrix `-dN` is solved in the tangent basis; the eigenvalues are
/// the principal curvatures. Exact where `d1` is analytic, approximate where it
/// is zero (sphere, cylinder ports).
pub fn surface_curvature_at(s: &dyn Surface, u: f64, v: f64) -> (f64, f64) {
    let eps = 1e-6;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-5 } else { eps };
    let hv = if v1 > v0 { (v1 - v0) * 1e-5 } else { eps };

    let (su, sv) = surface_tangents(s, u, v);
    if su.xyz().square_modulus() < 1e-30 || sv.xyz().square_modulus() < 1e-30 {
        return (0.0, 0.0);
    }
    let n0 = unit_normal(s, u, v);
    let du_n = unit_normal(s, u + hu, v).subtracted(&n0).divided(hu);
    let dv_n = unit_normal(s, u, v + hv).subtracted(&n0).divided(hv);

    // Shape operator W = -dN in the (Su, Sv) basis:
    // -du_n = w00 Su + w10 Sv ; -dv_n = w01 Su + w11 Sv.
    let e = su.dot(&su);
    let f = su.dot(&sv);
    let g = sv.dot(&sv);
    let det = e * g - f * f;
    if det.abs() < 1e-30 {
        return (0.0, 0.0);
    }
    let a0 = du_n.reversed().dot(&su);
    let a1 = du_n.reversed().dot(&sv);
    let b0 = dv_n.reversed().dot(&su);
    let b1 = dv_n.reversed().dot(&sv);
    let w00 = (a0 * g - a1 * f) / det;
    let w10 = (e * a1 - f * a0) / det;
    let w01 = (b0 * g - b1 * f) / det;
    let w11 = (e * b1 - f * b0) / det;

    // Eigenvalues of [[w00, w01], [w10, w11]].
    let tr = w00 + w11;
    let disc = (tr * tr - 4.0 * (w00 * w11 - w01 * w10)).max(0.0).sqrt();
    let k1 = 0.5 * (tr - disc);
    let k2 = 0.5 * (tr + disc);
    (k1, k2)
}

/// Wrap `u` into `[u0, u1)` when the surface is u-periodic (and `v` likewise).
pub fn surface_periodic_adjust(s: &dyn Surface, u: f64, v: f64) -> (f64, f64) {
    let wrap = |val: f64, a: f64, b: f64| -> f64 {
        let p = b - a;
        if p > 0.0 && p.is_finite() {
            a + (val - a).rem_euclid(p)
        } else {
            val
        }
    };
    let (mut u, mut v) = (u, v);
    if s.is_u_periodic() {
        let (a, b) = s.u_range();
        u = wrap(u, a, b);
    }
    if s.is_v_periodic() {
        let (a, b) = s.v_range();
        v = wrap(v, a, b);
    }
    (u, v)
}

// --- internals -------------------------------------------------------------

fn finite_range((a, b): (f64, f64)) -> (f64, f64) {
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

/// Unit normal of `s` at `(u, v)`; falls back to central finite differences of
/// `d0` when the `d1` cross product vanishes.
fn unit_normal(s: &dyn Surface, u: f64, v: f64) -> GpVec {
    let (_, du, dv) = s.d1(u, v);
    let n = du.xyz().crossed(dv.xyz());
    if n.square_modulus() > 1e-30 {
        let m = n.modulus();
        return GpVec::new(n.x / m, n.y / m, n.z / m);
    }
    let eps = 1e-6;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { eps };
    let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { eps };
    let p0 = s.d0(u, v);
    let duv = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
    let dvv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
    let n = duv.xyz().crossed(dvv.xyz());
    let m = n.modulus();
    if m > 1e-30 {
        GpVec::new(n.x / m, n.y / m, n.z / m)
    } else {
        GpVec::zero()
    }
}

/// Surface tangents at `(u, v)`; finite-difference fallback when `d1` is zero.
fn surface_tangents(s: &dyn Surface, u: f64, v: f64) -> (GpVec, GpVec) {
    let (_, du, dv) = s.d1(u, v);
    if du.xyz().square_modulus() > 1e-30 && dv.xyz().square_modulus() > 1e-30 {
        return (du, dv);
    }
    let eps = 1e-6;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-5 } else { eps };
    let hv = if v1 > v0 { (v1 - v0) * 1e-5 } else { eps };
    let p0 = s.d0(u, v);
    let pu = s.d0(u + hu, v);
    let pv = s.d0(u, v + hv);
    (
        GpVec::new(
            (pu.x() - p0.x()) / hu,
            (pu.y() - p0.y()) / hu,
            (pu.z() - p0.z()) / hu,
        ),
        GpVec::new(
            (pv.x() - p0.x()) / hv,
            (pv.y() - p0.y()) / hv,
            (pv.z() - p0.z()) / hv,
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeomCylinder, GeomPlane, GeomSphere};
    use occt_core::gp::{GpAx3, GpCylinder as GpCyl, GpPln, GpPnt, GpSphere as GpSph};

    #[test]
    fn iso_u_of_plane_is_line() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let c = iso_u(&plane, 0.5);
        // Plane z=0: (u, v, 0), fixed v=0.5 → line along +x.
        assert!(c.d0(0.0).distance(&GpPnt::new(0.0, 0.5, 0.0)) < 1e-12);
        assert!(c.d0(1.0).distance(&GpPnt::new(1.0, 0.5, 0.0)) < 1e-12);
        assert!(c.d0(2.0).distance(&GpPnt::new(2.0, 0.5, 0.0)) < 1e-12);
        let (_, d0) = c.d1(0.0);
        let (_, d1) = c.d1(1.0);
        assert!(d0.subtracted(&d1).magnitude() < 1e-12, "derivative not constant");
    }

    #[test]
    fn iso_u_of_cylinder_is_circle() {
        let cyl = GeomCylinder::new(GpCyl::new(GpAx3::standard(), 2.0).unwrap());
        let c = iso_u(&cyl, 3.0);
        assert!(c.first_parameter().abs() < 1e-12);
        assert!((c.last_parameter() - 2.0 * std::f64::consts::PI).abs() < 1e-12);
        // (r cos u, r sin u, v) at r=2, v=3.
        assert!(c.d0(0.0).distance(&GpPnt::new(2.0, 0.0, 3.0)) < 1e-9);
        assert!(c.d0(std::f64::consts::FRAC_PI_2).distance(&GpPnt::new(0.0, 2.0, 3.0)) < 1e-9);
    }

    #[test]
    fn iso_curves_grid_counts() {
        let s = GeomSphere::new(GpSph::new(GpAx3::standard(), 1.0).unwrap());
        let (u_const, v_const) = iso_curves_grid(&s, 4, 5);
        assert_eq!(u_const.len(), 4);
        assert_eq!(v_const.len(), 5);
        // u-constant curves fix u and vary v over the v-range.
        assert!((u_const[0].first_parameter() - (-std::f64::consts::FRAC_PI_2)).abs() < 1e-12);
        // v-constant curves fix v and vary u over the u-range.
        assert!(v_const[0].first_parameter().abs() < 1e-12);
    }

    #[test]
    fn reparameterize_surface_maps_new_range_onto_old() {
        let s = GeomSphere::new(GpSph::new(GpAx3::standard(), 1.0).unwrap());
        let rp = reparameterize_surface(&s, 0.0, 1.0, 0.0, 1.0);
        let (ru0, ru1) = rp.u_range();
        let (rv0, rv1) = rp.v_range();
        assert!((ru0 - 0.0).abs() < 1e-12 && (ru1 - 1.0).abs() < 1e-12);
        assert!((rv0 - 0.0).abs() < 1e-12 && (rv1 - 1.0).abs() < 1e-12);
        // New u=0.5 → old u=π; new v=0.5 → old v=0.
        let p = rp.d0(0.5, 0.5);
        let expected = s.d0(std::f64::consts::PI, 0.0);
        assert!(p.distance(&expected) < 1e-9, "p={p:?} expected={expected:?}");
    }

    #[test]
    fn normal_grid_of_plane_is_constant() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let grid = surface_normal_grid(&plane, 4, 4);
        for row in &grid {
            for n in row {
                assert!(n.x().abs() < 1e-12 && n.y().abs() < 1e-12, "n={n:?}");
                assert!((n.z() - 1.0).abs() < 1e-12, "n={n:?}");
            }
        }
    }

    #[test]
    fn curvature_of_sphere() {
        let s = GeomSphere::new(GpSph::new(GpAx3::standard(), 2.0).unwrap());
        let (k1, k2) = surface_curvature_at(&s, 0.0, 0.0);
        // Radius-2 sphere: both principal curvatures have magnitude 1/2.
        assert!((k1.abs() - 0.5).abs() < 0.05, "k1={k1}");
        assert!((k2.abs() - 0.5).abs() < 0.05, "k2={k2}");
    }

    #[test]
    fn periodic_adjust_wraps_u() {
        let cyl = GeomCylinder::new(GpCyl::new(GpAx3::standard(), 1.0).unwrap());
        let (u, v) = surface_periodic_adjust(&cyl, 3.0 * std::f64::consts::PI, 5.0);
        assert!((u - std::f64::consts::PI).abs() < 1e-12, "u={u}");
        // v is not periodic and stays untouched.
        assert!((v - 5.0).abs() < 1e-12);
    }
}
