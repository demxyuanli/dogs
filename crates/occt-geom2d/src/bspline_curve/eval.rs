use crate::curve::Curve2d;
use occt_core::bspl::{curve_tools, eval, knots};
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};

use super::curve::Geom2dBSplineCurve;

impl Curve2d for Geom2dBSplineCurve {
    fn d0(&self, u: f64) -> GpPnt2d {
        // `Geom2d_BSplineCurve::D0` (`Geom2d_BSplineCurve.cxx`) ->
        // `BSplCLib::D0` over `myKnots` + `&myMults`, i.e. the exact path the
        // 3D `GeomBSplineCurve::d0` already takes. A hand-rolled De Boor
        // evaluation differs from `BSplCLib` in the last bits, which is enough
        // to move a pcurve node across a range boundary (e.g. a seam arc whose
        // `v` sits on the surface bound) and change the range splitter seeds.
        if self.xs.is_empty() {
            return GpPnt2d::zero();
        }
        let (new_u, span, uknots, umults) = self.eval_span(u);
        let p = occt_core::bspl::curve_eval::curve_d0(
            new_u,
            span,
            self.degree as i32,
            self.periodic,
            &self.poles_3d(),
            self.weights.as_deref(),
            &uknots,
            Some(&umults),
        );
        GpPnt2d::new(p.x(), p.y())
    }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D1` / `BSplCLib::D1`.
        let poles = self.poles_3d();
        if self.periodic {
            let d = occt_core::bspl::curve_dn::curve_dn(
                u, 1, 0, self.degree as i32, true, &poles, self.weights.as_deref(), &self.knots,
                None,
            );
            return (self.d0(u), GpVec2d::new(d.x(), d.y()));
        }
        if let Some(w) = self.weights.as_deref() {
            // Rational uses the homogeneous quotient already computed by
            // `eval_curve_rational_d2`.
            let (p, d1, _) = eval::eval_curve_rational_d2(&poles, w, &self.knots, self.degree, u);
            return (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(d1.x(), d1.y()));
        }
        let (p, d) = eval::eval_curve_d1(&poles, &self.knots, self.degree, u);
        (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(d.x(), d.y()))
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D2` / `BSplCLib::D2`.
        let poles = self.poles_3d();
        if self.periodic {
            let d1 = occt_core::bspl::curve_dn::curve_dn(
                u, 1, 0, self.degree as i32, true, &poles, self.weights.as_deref(), &self.knots,
                None,
            );
            let d2 = occt_core::bspl::curve_dn::curve_dn(
                u, 2, 0, self.degree as i32, true, &poles, self.weights.as_deref(), &self.knots,
                None,
            );
            return (
                self.d0(u),
                GpVec2d::new(d1.x(), d1.y()),
                GpVec2d::new(d2.x(), d2.y()),
            );
        }
        if let Some(w) = self.weights.as_deref() {
            let (p, d1, d2) = eval::eval_curve_rational_d2(&poles, w, &self.knots, self.degree, u);
            return (
                GpPnt2d::new(p.x(), p.y()),
                GpVec2d::new(d1.x(), d1.y()),
                GpVec2d::new(d2.x(), d2.y()),
            );
        }
        let (p, d1, d2) = eval::eval_curve_d2(&poles, &self.knots, self.degree, u);
        (
            GpPnt2d::new(p.x(), p.y()),
            GpVec2d::new(d1.x(), d1.y()),
            GpVec2d::new(d2.x(), d2.y()),
        )
    }

    fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }
    fn is_periodic(&self) -> bool { self.periodic }
    fn period(&self) -> f64 {
        if self.periodic {
            self.last_parameter() - self.first_parameter()
        } else {
            0.0
        }
    }
    fn continuity(&self) -> u8 {
        // `Geom2d_BSplineCurve::Continuity` / `GeomAdaptor` LocalContinuity.
        occt_core::bspl::local_continuity(
            &self.knots,
            self.degree,
            self.periodic,
            self.first_parameter(),
            self.last_parameter(),
        )
    }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        occt_core::bspl::adaptor_intervals(
            &self.knots,
            self.degree,
            self.periodic,
            continuity,
            self.first_parameter(),
            self.last_parameter(),
            occt_core::precision::PCONFUSION,
        )
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }

    fn transform(&mut self, t: &GpTrsf2d) {
        for i in 0..self.xs.len() {
            let mut p = GpPnt2d::new(self.xs[i], self.ys[i]);
            p.transform(t);
            self.xs[i] = p.x();
            self.ys[i] = p.y();
        }
    }

    fn reverse(&mut self) {
        // `Geom2d_BSplineCurve::Reverse` (`Geom2d_BSplineCurve.cxx:677-696`):
        // `BSplCLib::Reverse(myKnots)` + `BSplCLib::Reverse(myMults)` + reverse
        // the poles (+ weights) + `updateKnots()`. The knot walk
        // (`BSplCLib.cxx:802-824`) keeps `Knots(Lower)` and rebuilds the gaps,
        // so a first knot below half an ulp of the last knot is not flushed to
        // 0. Reflecting the flat vector (`K -> K_first + K_last - K`) does that
        // flush. Periodic poles use the same window as the 3D curve
        // (`Geom_BSplineCurve.cxx:502-510`).
        let (mut uknots, mut umults) = self.distinct_knots_and_mults();
        curve_tools::reverse_distinct_knots(&mut uknots);
        umults.reverse();
        if self.periodic {
            curve_tools::reverse_periodic_span(&mut self.xs, self.knots.len(), self.degree);
            curve_tools::reverse_periodic_span(&mut self.ys, self.knots.len(), self.degree);
            if let Some(w) = self.weights.as_mut() {
                curve_tools::reverse_periodic_span(w, self.knots.len(), self.degree);
            }
            self.knots = knots::knot_sequence_periodic(&uknots, &umults, self.degree as i32);
        } else {
            self.xs.reverse();
            self.ys.reverse();
            if let Some(w) = self.weights.as_mut() {
                w.reverse();
            }
            self.knots =
                occt_core::bspl::banded_interp::knot_sequence(&uknots, &umults, self.degree as i32);
        }
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
    fn is_bspline2d(&self) -> bool {
        true
    }
    /// `Geom2dAdaptor_Curve::Resolution`'s `GeomAbs_BSplineCurve` arm
    /// (`Geom2dAdaptor_Curve.cxx:1211-1215`) -> `Geom2d_BSplineCurve::Resolution`
    /// (`Geom2d_BSplineCurve_1.cxx:764-806`), implemented as the inherent
    /// [`Geom2dBSplineCurve::resolution`].
    fn resolution(&self, ruv: f64) -> f64 {
        Geom2dBSplineCurve::resolution(self, ruv)
    }
    /// `Geom2d_BSplineCurve::IsRational()` (`Geom2dAdaptor_Curve.cxx:1296-1298`).
    fn is_rational(&self) -> bool {
        self.weights.is_some()
    }
    /// `Geom2d_BSplineCurve::Weights()` (`GeomLib.cxx:598-620`).
    fn bspline_weights2d(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }
    fn poles2d(&self) -> Option<Vec<GpPnt2d>> {
        Some(self.xs.iter().zip(self.ys.iter()).map(|(x, y)| GpPnt2d::new(*x, *y)).collect())
    }
    fn set_poles2d(&mut self, poles: &[GpPnt2d]) {
        for (i, p) in poles.iter().enumerate() {
            self.xs[i] = p.x();
            self.ys[i] = p.y();
        }
    }

    /// `Geom2d_BSplineCurve::NbKnots()` (`Geom2d_BSplineCurve_1.cxx:598-601`):
    /// the number of distinct knots. The curve stores the expanded knot
    /// sequence, so compress it as `Geom_BSplineCurve::Knots` does.
    fn bspline_nb_knots(&self) -> Option<usize> {
        Some(knots::unique_knots_mults(&self.knots).0.len())
    }

    /// `Geom2d_BSplineCurve::Knots()` flat knot sequence (`Knot(j)` in
    /// `ShapeAnalysis_TransferParametersProj::CorrectParameter`,
    /// `Proj.cxx:268-279`).
    fn bspline_knots2d(&self) -> Option<&[f64]> {
        Some(&self.knots)
    }

    /// `Geom2d_BSplineCurve::Degree()` (`Geom2d_BSplineCurve_1.cxx:168-171`).
    fn bspline_degree(&self) -> Option<usize> {
        Some(self.degree)
    }

    fn bspline_poles2d(&self) -> Option<(&[f64], &[f64])> {
        Some((&self.xs, &self.ys))
    }

    /// `Geom2d_BSplineCurve::Copy()` (`Geom2d_BSplineCurve.cxx:109-136`): a deep
    /// copy of the stored poles/knots. The port's clone is that copy.
    fn bspline_copy2d(&self) -> Option<Geom2dBSplineCurve> {
        Some(self.clone())
    }

    /// `Geom2d_BSplineCurve::Knots()` / `Multiplicities()` as the distinct
    /// knots + multiplicities pair that the `Geom_BSplineCurve` constructor
    /// takes (`GeomLib.cxx:610-618`). The port stores no weights, so the
    /// rational branch of `GeomLib::To3d` cannot be reached from here.
    fn bspline_distinct_knots_mults(&self) -> Option<(Vec<f64>, Vec<i32>)> {
        Some(self.distinct_knots_and_mults())
    }
}
