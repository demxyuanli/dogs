use crate::curve::Curve;
use occt_core::bspl::{curve_tools, knots};
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

use super::curve::GeomBSplineCurve;
use super::distinct_knots_and_mults_of;

impl GeomBSplineCurve {
    /// `Geom_BSplineCurve::EvalD0/D1/D2/D3` preamble
    /// (`Geom_BSplineCurve_1.cxx:182-190`).
    ///
    /// OCCT hands `BSplCLib::D0/D1/D2/D3` the *distinct* knots plus `&myMults`
    /// (not `myFlatKnots` / `NoMults`), together with a span index located in
    /// that distinct-knot array and post-corrected by `aSpanIndex--` when
    /// `aNewU` sits below `myKnots(aSpanIndex)`. The flat-knot arm is only used
    /// by `LocalD0..LocalD3` and `EvalDN`.
    fn eval_span(&self, u: f64) -> (f64, i32, &[f64], &[i32]) {
        let uknots = self.distinct.as_slice();
        let umults = self.mults.as_slice();
        let (mut span, new_u) = occt_core::bspl::locate::locate_parameter(
            self.degree as i32,
            uknots,
            Some(umults),
            u,
            self.periodic,
            0,
        );
        if span >= 1 {
            if let Some(&k) = uknots.get(span as usize - 1) {
                if new_u < k {
                    span -= 1;
                }
            }
        }
        // Debug builds re-derive the pair and compare, so a mutation path that
        // forgets to refresh the representation shows up in the existing parity
        // tests instead of silently shifting geometry.
        debug_assert!(
            {
                let (d, m) = distinct_knots_and_mults_of(&self.knots, self.degree, self.periodic);
                d.as_slice() == uknots && m.as_slice() == umults
            },
            "GeomBSplineCurve: stored distinct knots out of date"
        );
        (new_u, span, uknots, umults)
    }

    /// Finite-difference first derivative: (f(u+h) - f(u-h)) / (2h).
    fn fd_d1(&self, u: f64) -> GpVec {
        let h = 1e-6;
        let p1 = self.d0(u + h);
        let p2 = self.d0(u - h);
        GpVec::new(
            (p1.x() - p2.x()) / (2.0 * h),
            (p1.y() - p2.y()) / (2.0 * h),
            (p1.z() - p2.z()) / (2.0 * h),
        )
    }

    /// Finite-difference second derivative: (f(u+h) - 2f(u) + f(u-h)) / h^2.
    fn fd_d2(&self, u: f64) -> GpVec {
        let h = 1e-6;
        let p1 = self.d0(u + h);
        let p0 = self.d0(u);
        let p2 = self.d0(u - h);
        GpVec::new(
            (p1.x() - 2.0 * p0.x() + p2.x()) / (h * h),
            (p1.y() - 2.0 * p0.y() + p2.y()) / (h * h),
            (p1.z() - 2.0 * p0.z() + p2.z()) / (h * h),
        )
    }

    /// Private de Boor evaluation. The public d0 delegates to `eval`; this is
    /// kept as a self-contained reference.
    #[allow(dead_code)]
    fn de_boor(&self, u: f64) -> GpPnt {
        let n = self.poles.len();
        if n == 0 { return GpPnt::zero(); }
        let idx = knots::hunt(&self.knots, u).max(self.degree).min(n - 1);
        let mut pts = vec![GpPnt::zero(); self.degree + 1];
        for k in 0..=self.degree {
            let pi = idx - self.degree + k;
            pts[k] = if pi < n { self.poles[pi] } else { self.poles[n - 1] };
        }
        for r in 1..=self.degree {
            for i in (r..=self.degree).rev() {
                let k0 = idx + i - self.degree;
                let k1 = k0 + self.degree + 1 - r;
                let alpha = (u - self.knots[k0]) / (self.knots[k1] - self.knots[k0]);
                if alpha.is_finite() {
                    pts[i] = GpPnt::new(
                        (1.0 - alpha) * pts[i - 1].x() + alpha * pts[i].x(),
                        (1.0 - alpha) * pts[i - 1].y() + alpha * pts[i].y(),
                        (1.0 - alpha) * pts[i - 1].z() + alpha * pts[i].z(),
                    );
                }
            }
        }
        pts[self.degree]
    }
}

impl Curve for GeomBSplineCurve {
    fn d0(&self, u: f64) -> GpPnt {
        // `Geom_BSplineCurve::EvalD0` (`Geom_BSplineCurve_1.cxx:174-194`) ->
        // `BSplCLib::D0` over `myKnots` + `&myMults` (distinct knots).
        let (new_u, span, uknots, umults) = self.eval_span(u);
        occt_core::bspl::curve_eval::curve_d0(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        )
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        // `Geom_BSplineCurve::EvalD1` (`:198-227`) -> `BSplCLib::D1`.
        let (new_u, span, uknots, umults) = self.eval_span(u);
        occt_core::bspl::curve_eval::curve_d1(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        )
    }

    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        // `Geom_BSplineCurve::EvalDN` (`Geom_BSplineCurve_1.cxx:300-316`).
        // Illegal N<1 returns zero instead of throw. Eval-rep is empty.
        if n < 1 {
            return GpVec::zero();
        }
        occt_core::bspl::curve_dn::curve_dn(
            u,
            n,
            0,
            self.degree as i32,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.knots,
            None,
        )
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        // `Geom_BSplineCurve::EvalD2` (`:231-261`) -> `BSplCLib::D2`.
        let (new_u, span, uknots, umults) = self.eval_span(u);
        let r = occt_core::bspl::curve_eval::curve_d2(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        );
        r
    }

    /// `Geom_BSplineCurve::EvalD3` (`:265-296`) -> `BSplCLib::D3`, which
    /// `Geom_OffsetCurve` calls through `basisCurve->D3(U, …)` in
    /// `CalculateD2`. Without this override the trait default returned a
    /// **zero** third derivative for every B-spline basis.
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (new_u, span, uknots, umults) = self.eval_span(u);
        occt_core::bspl::curve_eval::curve_d3(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        )
    }

    fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }
    fn is_periodic(&self) -> bool { self.periodic }
    fn continuity(&self) -> u8 {
        occt_core::bspl::local_continuity(
            &self.knots,
            self.degree,
            self.periodic,
            self.first_parameter(),
            self.last_parameter(),
        )
    }

    fn transform(&mut self, t: &GpTrsf) {
        for p in self.poles.iter_mut() {
            *p = p.transformed(t);
        }
    }

    fn reverse(&mut self) {
        // `Geom_BSplineCurve::Reverse` (`cxx:496-511`): `BSplCLib::Reverse` on
        // the distinct knots and on the multiplicities, then the poles
        // (periodic window `myFlatKnots.Upper() - Deg - 1`, `cxx:502-510`),
        // then `updateKnots()`.
        let (mut uknots, mut umults) = self.distinct_knots_and_mults();
        curve_tools::reverse_distinct_knots(&mut uknots);
        umults.reverse();
        if self.periodic {
            curve_tools::reverse_periodic_span(&mut self.poles, self.knots.len(), self.degree);
            if let Some(w) = self.weights.as_mut() {
                curve_tools::reverse_periodic_span(w, self.knots.len(), self.degree);
            }
            self.knots = knots::knot_sequence_periodic(&uknots, &umults, self.degree as i32);
        } else {
            self.poles.reverse();
            if let Some(w) = self.weights.as_mut() {
                w.reverse();
            }
            self.knots =
                occt_core::bspl::banded_interp::knot_sequence(&uknots, &umults, self.degree as i32);
        }
        self.refresh_knot_rep();
    }

    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
    fn bspline_poles(&self) -> Option<&[GpPnt]> { Some(&self.poles) }
    fn bspline_weights(&self) -> Option<&[f64]> { self.weights.as_deref() }
    fn bspline_knots(&self) -> Option<&[f64]> { Some(&self.knots) }
    fn nurbs_degree(&self) -> Option<usize> { Some(self.degree) }
    fn resolution(&self, r3d: f64) -> f64 {
        occt_core::bspl::bspline_curve_resolution(
            &self.poles,
            self.weights.as_deref(),
            &self.knots,
            self.degree as i32,
            r3d,
        )
    }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        let eps = self
            .resolution(occt_core::precision::CONFUSION)
            .min(occt_core::precision::PCONFUSION);
        occt_core::bspl::adaptor_intervals(
            &self.knots,
            self.degree,
            self.periodic,
            continuity,
            self.first_parameter(),
            self.last_parameter(),
            eps,
        )
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity).len().saturating_sub(1).max(1) as i32
    }
}
