//! `Surface` trait implementation of `Geom_BSplineSurface` (OCCT
//! `Geom_BSplineSurface` / `Geom_Surface` overrides) plus the `Rational` /
//! `Resolution` flag helpers.

use super::*;
use crate::curve::Curve;
use crate::surface::Surface;
use occt_core::gp::{GpTrsf, GpVec};
use std::sync::Arc;

impl GeomBSplineSurface {
    /// `Epsilon` of a non negative value (`Standard_Real.hxx:242-246`): the
    /// distance from the value to the next double toward `RealLast()`, i.e. one
    /// ULP. `Rational` always calls it with `std::abs(...)`.
    fn occt_epsilon(x: f64) -> f64 {
        if x > 0.0 {
            if x.is_infinite() {
                return f64::NEG_INFINITY;
            }
            f64::from_bits(x.to_bits() + 1) - x
        } else if x < 0.0 {
            if x.is_infinite() {
                return f64::INFINITY;
            }
            f64::from_bits(x.to_bits() - 1) - x
        } else if x == 0.0 {
            f64::MIN_POSITIVE
        } else {
            x
        }
    }

    /// `Geom_BSplineSurface::Rational` (`Geom_BSplineSurface.cxx:110-137`).
    ///
    /// Returns the two flags in the order `Geom_BSplineSurface::Resolution`
    /// forwards them to `BSplSLib::Resolution(..., URational, VRational, ...)`
    /// (`Geom_BSplineSurface_1.cxx:2203-2217`). The stored flags are swapped
    /// relative to their names: the scan over adjacent U poles sets the flag
    /// which `Resolution` forwards as `VRational`, and the scan over adjacent V
    /// poles sets the one forwarded as `URational`
    /// (`Geom_BSplineSurface.hxx:1312-1313`). Weights are indexed [U][V].
    fn occt_rational_flags(&self) -> (bool, bool) {
        let Some(w) = self.weights.as_ref() else {
            return (false, false);
        };
        let nu = w.len();
        if nu == 0 {
            return (false, false);
        }
        let nv = w[0].len();

        let mut v_rational = false;
        'u_poles: for j in 0..nv {
            for i in 0..nu.saturating_sub(1) {
                if (w[i][j] - w[i + 1][j]).abs() > Self::occt_epsilon(w[i][j].abs()) {
                    v_rational = true;
                    break 'u_poles;
                }
            }
        }

        let mut u_rational = false;
        'v_poles: for i in 0..nu {
            for j in 0..nv.saturating_sub(1) {
                if (w[i][j] - w[i][j + 1]).abs() > Self::occt_epsilon(w[i][j].abs()) {
                    u_rational = true;
                    break 'v_poles;
                }
            }
        }

        (u_rational, v_rational)
    }
}

impl Surface for GeomBSplineSurface {
    /// `Geom_BSplineSurface::EvalD0` / `D0` (`Geom_BSplineSurface_1.cxx:131-159`)
    /// via `BSplSLib::D0`.
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.eval_d0_bspl(u, v)
    }

    /// `Geom_BSplineSurface::EvalD1` / `D1` (`Geom_BSplineSurface_1.cxx:163-195`)
    /// via `BSplSLib::D1`.
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        self.eval_d1_bspl(u, v)
    }

    /// `Geom_BSplineSurface::EvalD2` / `D2` (`Geom_BSplineSurface_1.cxx:199-235`)
    /// via `BSplSLib::D2`.
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        self.eval_d2_bspl(u, v)
    }

    fn u_range(&self) -> (f64, f64) {
        knot_window(&self.knots_u, self.deg_u)
    }

    fn v_range(&self) -> (f64, f64) {
        knot_window(&self.knots_v, self.deg_v)
    }

    fn continuity(&self) -> u8 {
        self.deg_u.min(self.deg_v).min(u8::MAX as usize) as u8
    }

    /// `Geom_BSplineSurface::IsCNu` (`Geom_BSplineSurface_1.cxx:50-76`).
    fn is_cnu(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_BSplineSurface::IsCNu: Geom_BSplineSurface_1.cxx:52".into());
        }
        let (a, b) = self.u_range();
        let shape = occt_core::bspl::local_continuity(
            &self.knots_u, self.deg_u, self.u_periodic, a, b,
        );
        Ok(crate::surface::smooth_is_cn(shape, n))
    }

    /// `Geom_BSplineSurface::IsCNv` (`Geom_BSplineSurface_1.cxx:80-106`).
    fn is_cnv(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_BSplineSurface::IsCNv: Geom_BSplineSurface_1.cxx:82".into());
        }
        let (a, b) = self.v_range();
        let shape = occt_core::bspl::local_continuity(
            &self.knots_v, self.deg_v, self.v_periodic, a, b,
        );
        Ok(crate::surface::smooth_is_cn(shape, n))
    }

    fn transform(&mut self, t: &GpTrsf) {
        // `Geom_BSplineSurface::Transform` clears `myMaxDerivInvOk`
        // (`Geom_BSplineSurface_1.cxx:935`).
        self.res_latch.invalidate();
        for row in self.poles.iter_mut() {
            for p in row.iter_mut() {
                *p = p.transformed(t);
            }
        }
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }

    fn is_bspline_surface(&self) -> bool {
        true
    }

    fn osculating_bspline(&self) -> Option<GeomBSplineSurface> {
        Some(self.clone())
    }

    fn eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> GpVec {
        self.eval_dn_bspl(u, v, nu, nv)
    }

    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        Ok(self.eval_dn(u, v, nu, nv))
    }

    /// `Geom_BSplineSurface::EvalD3` (`Geom_BSplineSurface_1.cxx:239-274`):
    /// `D3UUV` is `(2,1)`, `D3UVV` is `(1,2)`.
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        let d3 = self.eval_d3_bspl(u, v);
        Ok(match (nu, nv) {
            (2, 1) => d3.d3_uuv,
            (1, 2) => d3.d3_uvv,
            (3, 0) => d3.d3u,
            (0, 3) => d3.d3v,
            _ => GpVec::zero(),
        })
    }

    fn is_u_periodic(&self) -> bool {
        self.u_periodic
    }

    /// `Geom_Surface::UReversed` (`Geom_Surface.cxx:33-38`) = `Copy()` +
    /// [`Self::u_reverse`].
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        if self.knots_u.len() < 2 || self.poles.is_empty() {
            return None;
        }
        let mut out = self.clone();
        out.u_reverse();
        Some(Arc::new(out) as Arc<dyn Surface>)
    }

    /// `Geom_Surface::VReversed` (`Geom_Surface.cxx:43-48`) = `Copy()` +
    /// [`Self::v_reverse`].
    fn v_reversed(&self) -> Option<Arc<dyn Surface>> {
        if self.knots_v.len() < 2 || self.poles.is_empty() {
            return None;
        }
        let mut out = self.clone();
        out.v_reverse();
        Some(Arc::new(out) as Arc<dyn Surface>)
    }

    /// `Geom_BSplineSurface::UReversedParameter`
    /// (`Geom_BSplineSurface_1.cxx:1594-1598`):
    /// `FlatU(1) + FlatU(Upper) - U`.
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        match (self.knots_u.first(), self.knots_u.last()) {
            (Some(a), Some(b)) => a + b - u,
            _ => u,
        }
    }

    /// `Geom_BSplineSurface::VReversedParameter`
    /// (`Geom_BSplineSurface_1.cxx:1625-1629`):
    /// `FlatV(1) + FlatV(Upper) - V`.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        match (self.knots_v.first(), self.knots_v.last()) {
            (Some(a), Some(b)) => a + b - v,
            _ => v,
        }
    }

    fn is_v_periodic(&self) -> bool {
        self.v_periodic
    }

    fn is_u_closed(&self) -> bool {
        let nu = self.poles.len();
        if nu < 2 {
            return false;
        }
        let a = &self.poles[0];
        let b = &self.poles[nu - 1];
        if a.len() != b.len() {
            return false;
        }
        a.iter()
            .zip(b.iter())
            .all(|(p, q)| p.square_distance(q) <= 1e-14)
    }

    fn is_v_closed(&self) -> bool {
        if self.poles.is_empty() || self.poles[0].len() < 2 {
            return false;
        }
        let nv = self.poles[0].len();
        self.poles.iter().all(|row| {
            row.len() == nv && row[0].square_distance(&row[nv - 1]) <= 1e-14
        })
    }

    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        self.u_iso(u).ok().map(|c| Arc::new(c) as Arc<dyn Curve>)
    }

    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        self.v_iso(v).ok().map(|c| Arc::new(c) as Arc<dyn Curve>)
    }

    fn u_degree(&self) -> i32 {
        self.deg_u as i32
    }
    fn v_degree(&self) -> i32 {
        self.deg_v as i32
    }
    fn bspline_surface_poles(&self) -> Option<&[Vec<GpPnt>]> { Some(&self.poles) }

    fn bspline_surface_uknots(&self) -> Option<&[f64]> { Some(&self.knots_u) }

    fn bspline_surface_vknots(&self) -> Option<&[f64]> { Some(&self.knots_v) }

    fn bspline_surface_weights(&self) -> Option<&[Vec<f64>]> { self.weights.as_deref() }

    fn nb_u_poles(&self) -> i32 {
        self.nb_poles_u() as i32
    }
    fn nb_v_poles(&self) -> i32 {
        self.nb_poles_v() as i32
    }
    fn nb_u_intervals(&self, _continuity: u8) -> i32 {
        Self::unique_knots(&self.knots_u).len().saturating_sub(1) as i32
    }
    fn nb_v_intervals(&self, _continuity: u8) -> i32 {
        Self::unique_knots(&self.knots_v).len().saturating_sub(1) as i32
    }
    fn u_intervals(&self, _continuity: u8) -> Vec<f64> {
        Self::unique_knots(&self.knots_u)
    }
    fn v_intervals(&self, _continuity: u8) -> Vec<f64> {
        Self::unique_knots(&self.knots_v)
    }

    /// `Geom_BSplineSurface::Resolution` (`Geom_BSplineSurface_1.cxx:2197-2221`):
    /// the `BSplSLib::Resolution` bound is latched behind `myMaxDerivInvOk` and
    /// each query divides it by `Tolerance3D`.
    fn uv_resolution(&self, r3d: f64) -> Option<(f64, f64)> {
        let (u_rational, v_rational) = self.occt_rational_flags();
        let (max_deriv_u, max_deriv_v) = self.res_latch.or_compute(|| {
            occt_core::bspl::bspline_surface_max_derivative(
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                self.deg_u as i32,
                self.deg_v as i32,
                u_rational,
                v_rational,
            )
        });
        if max_deriv_u != 0.0 && max_deriv_v != 0.0 {
            Some((r3d / max_deriv_u, r3d / max_deriv_v))
        } else {
            Some((0.0, 0.0))
        }
    }
}
