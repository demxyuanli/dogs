//! `Surface` for `GeomBezierSurface`.

use std::sync::Arc;

use super::{de_casteljau, GeomBezierSurface};
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;
use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

impl Surface for GeomBezierSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        if let Some(bs) = self.rational_bspline() {
            return bs.d0(u, v);
        }
        let mut rows = Vec::with_capacity(self.n_u);
        for i in 0..self.n_u {
            let start = i * self.n_v;
            rows.push(de_casteljau(&self.poles[start..start + self.n_v], v));
        }
        de_casteljau(&rows, u)
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        if let Some(bs) = self.rational_bspline() {
            return bs.d1(u, v);
        }
        (self.d0(u, v), self.du(u, v), self.dv(u, v))
    }

    /// `Geom_BezierSurface::D2`. Exact through the equivalent B-spline
    /// representation (a Bezier is a B-spline on `[0,1]` with full knot
    /// multiplicities, so all derivatives coincide), instead of the trait's
    /// `h = 1e-6` central difference (audit A15). The only case without that
    /// representation is a degenerate surface with a single pole row
    /// (`n_u < 2 || n_v < 2`), where `d1` is still analytic and the second
    /// derivatives vanish.
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        if let Some(bs) = self.rational_bspline() {
            return bs.d2(u, v);
        }
        match self.osculating_bspline() {
            Some(bs) => bs.d2(u, v),
            None => {
                let (p, du, dv) = self.d1(u, v);
                (p, du, dv, GpVec::zero(), GpVec::zero(), GpVec::zero())
            }
        }
    }

    /// `Geom_BezierSurface::EvalDN` (`Geom_BezierSurface.cxx:1666-1724`).
    /// Both rational and polynomial arms call `BSplSLib::DN` on the Bezier
    /// knot sequence. Polynomial uses the equivalent non-rational B-spline.
    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if let Some(bs) = self.rational_bspline() {
            return Ok(bs.eval_dn(u, v, nu, nv));
        }
        if let Some(bs) = self.osculating_bspline() {
            return Ok(bs.eval_dn(u, v, nu, nv));
        }
        Err("UNPORTED Geom_Surface::EvalDN".into())
    }

    /// `Geom_BezierSurface::EvalD3` (`Geom_BezierSurface.cxx:1594-1661`) via
    /// `BSplSLib::D3` on the equivalent B-spline knot sequence.
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if let Some(bs) = self.rational_bspline() {
            return bs.try_eval_d3_uv(u, v, nu, nv);
        }
        if let Some(bs) = self.osculating_bspline() {
            return bs.try_eval_d3_uv(u, v, nu, nv);
        }
        Err("UNPORTED Geom_BezierSurface::EvalD3: Geom_BezierSurface.cxx:1594".into())
    }

    fn osculating_bspline(&self) -> Option<GeomBSplineSurface> {
        if self.n_u < 2 || self.n_v < 2 {
            return None;
        }
        let deg_u = self.n_u - 1;
        let deg_v = self.n_v - 1;
        let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); self.n_v]; self.n_u];
        for i in 0..self.n_u {
            for j in 0..self.n_v {
                poles[i][j] = self.poles[i * self.n_v + j];
            }
        }
        GeomBSplineSurface::from_poles_knots_mults(
            poles,
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            vec![(deg_u + 1) as i32, (deg_u + 1) as i32],
            vec![(deg_v + 1) as i32, (deg_v + 1) as i32],
            deg_u,
            deg_v,
        )
        .ok()
    }

    fn u_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    /// `IsKind(STANDARD_TYPE(Geom_BezierSurface))` (`GeomAdaptor_Surface.cxx:480`).
    fn is_bezier_surface(&self) -> bool {
        true
    }

    fn bezier_control_grid(&self) -> Option<(Vec<GpPnt>, usize, usize)> {
        Some((self.poles.clone(), self.n_u, self.n_v))
    }

    fn bezier_weights(&self) -> Option<(Vec<f64>, bool, bool)> {
        let weights = self.weights.as_ref()?;
        if !(self.u_rational || self.v_rational) || weights.len() != self.n_u * self.n_v {
            return None;
        }
        Some((weights.clone(), self.u_rational, self.v_rational))
    }

    fn v_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    /// `Geom_BezierSurface::Continuity` (`Geom_BezierSurface.cxx:1409-1412`).
    fn continuity(&self) -> u8 {
        occt_core::kernel::geomabs::Shape::CN as u8
    }

    /// `Geom_BezierSurface::IsUClosed` (`Geom_BezierSurface.cxx:1944-1957`).
    fn is_u_closed(&self) -> bool {
        let upper = self.n_u - 1;
        for j in 0..self.n_v {
            let a = self.poles[j];
            let b = self.poles[upper * self.n_v + j];
            if a.distance(&b) > occt_core::precision::CONFUSION {
                return false;
            }
        }
        true
    }

    /// `Geom_BezierSurface::IsVClosed` (`Geom_BezierSurface.cxx:1961-1973`).
    fn is_v_closed(&self) -> bool {
        let upper = self.n_v - 1;
        for i in 0..self.n_u {
            let a = self.poles[i * self.n_v];
            let b = self.poles[i * self.n_v + upper];
            if a.distance(&b) > occt_core::precision::CONFUSION {
                return false;
            }
        }
        true
    }

    /// `Geom_BezierSurface::IsUPeriodic` (`Geom_BezierSurface.cxx:1977-1980`).
    fn is_u_periodic(&self) -> bool {
        false
    }

    /// `Geom_BezierSurface::IsVPeriodic` (`Geom_BezierSurface.cxx:1984-1987`).
    fn is_v_periodic(&self) -> bool {
        false
    }

    /// `Geom_BezierSurface::IsCNu` (`Geom_BezierSurface.cxx:1901-1904`).
    fn is_cnu(&self, _n: i32) -> Result<bool, String> {
        Ok(true)
    }

    /// `Geom_BezierSurface::IsCNv` (`Geom_BezierSurface.cxx:1908`).
    fn is_cnv(&self, _n: i32) -> Result<bool, String> {
        Ok(true)
    }

    /// `Geom_BezierSurface::Transform` (`Geom_BezierSurface.cxx:1929-1940`).
    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.poles {
            *p = p.transformed(t);
        }
    }

    fn uv_resolution(&self, r3d: f64) -> Option<(f64, f64)> {
        Some(self.resolution(r3d))
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }

    fn u_degree(&self) -> i32 {
        GeomBezierSurface::u_degree(self) as i32
    }
    fn v_degree(&self) -> i32 {
        GeomBezierSurface::v_degree(self) as i32
    }

    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        self.u_iso(u).ok().map(|c| Arc::new(c) as Arc<dyn Curve>)
    }

    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        self.v_iso(v).ok().map(|c| Arc::new(c) as Arc<dyn Curve>)
    }

    fn nb_u_poles(&self) -> i32 {
        self.n_u as i32
    }
    fn nb_v_poles(&self) -> i32 {
        self.n_v as i32
    }

    /// `Geom_Surface::UReversed` = `Copy()` + `UReverse` (`Geom_Surface.cxx:33-38`).
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut copy = self.clone();
        copy.u_reverse();
        Some(Arc::new(copy))
    }
}
