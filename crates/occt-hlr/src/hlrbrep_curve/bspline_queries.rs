//! `HLRBRep_Curve` B-spline / Bezier queries: `IsRational` (`lxx:165-168`),
//! `Degree` (`lxx:145-148`), `NbPoles` (`lxx:152-155`), `Poles`
//! (`cxx:476-508`), `PolesAndWeights` (`cxx:533-603`), `NbKnots`
//! (`lxx:158-161`), `Knots` (`cxx:606-619`) and `Multiplicities`
//! (`cxx:621-633`).
//!
//! OCCT fills caller-owned `NCollection_Array1` buffers with the `Lower()` ..
//! `Upper()` bounds of the destination and indexes the source poles with the
//! same 1-based indices; the port answers owned `Vec`s in source order.

use occt_core::gp::{GpPnt, GpPnt2d};

use super::HlrBRepCurve;

impl HlrBRepCurve {
    /// Source poles of the current curve: `Geom_BSplineCurve::Poles()` when
    /// the adaptor type is `GeomAbs_BSplineCurve`, `Geom_BezierCurve::Poles()`
    /// otherwise (`cxx:481-500`).
    fn source_poles(&self) -> Option<&[GpPnt]> {
        let c = self.curve.as_ref()?;
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => c.bspline_poles(),
            _ => c.bezier_poles(),
        }
    }

    /// Paired poles / weights of the current curve (`cxx:540-566`).
    fn source_poles_and_weights(&self) -> Option<(&[GpPnt], Option<&[f64]>)> {
        let c = self.curve.as_ref()?;
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => {
                Some((c.bspline_poles()?, c.bspline_weights()))
            }
            _ => Some((c.bezier_poles()?, c.bezier_weights())),
        }
    }

    /// `IsRational()` (`lxx:165-168`).
    pub fn is_rational(&self) -> bool {
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => self
                .curve
                .as_ref()
                .and_then(|c| c.as_bspline_curve())
                .map(|b| b.is_rational())
                .unwrap_or(false),
            _ => self
                .curve
                .as_ref()
                .and_then(|c| c.bezier_weights())
                .is_some(),
        }
    }

    /// `Degree()` (`lxx:145-148`): the B-spline degree or `NbPoles() - 1`.
    pub fn degree(&self) -> i32 {
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => self
                .curve
                .as_ref()
                .and_then(|c| c.as_bspline_curve())
                .map(|b| b.degree() as i32)
                .unwrap_or(0),
            _ => self.source_poles().map(|p| p.len() as i32 - 1).unwrap_or(0),
        }
    }

    /// `NbPoles()` (`lxx:152-155`).
    pub fn nb_poles(&self) -> i32 {
        self.source_poles().map(|p| p.len() as i32).unwrap_or(0)
    }

    /// `Poles(TP)` (`cxx:476-508`): every pole transformed by the projector,
    /// keeping only `X` and `Y`.
    pub fn poles(&self) -> Vec<GpPnt2d> {
        let Some(src) = self.source_poles() else {
            return Vec::new();
        };
        let mut out = Vec::with_capacity(src.len());
        for p in src {
            let mut p3 = *p;
            if let Some(proj) = &self.proj {
                proj.transform_pnt(&mut p3);
            }
            out.push(GpPnt2d::new(p3.x(), p3.y()));
        }
        out
    }

    /// `PolesAndWeights(TP, TW)` (`cxx:533-603`).
    pub fn poles_and_weights(&self) -> (Vec<GpPnt2d>, Vec<f64>) {
        let Some((src, weights)) = self.source_poles_and_weights() else {
            return (Vec::new(), Vec::new());
        };
        let mut tp = Vec::with_capacity(src.len());
        for p in src {
            let mut p3 = *p;
            if let Some(proj) = &self.proj {
                proj.transform_pnt(&mut p3);
            }
            tp.push(GpPnt2d::new(p3.x(), p3.y()));
        }
        let tw = match weights {
            Some(w) => w.to_vec(),
            None => vec![1.0; src.len()],
        };
        (tp, tw)
    }

    /// `NbKnots()` (`lxx:158-161`).
    pub fn nb_knots(&self) -> i32 {
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => self
                .curve
                .as_ref()
                .and_then(|c| c.as_bspline_curve())
                .map(|b| b.nb_knots() as i32)
                .unwrap_or(0),
            _ => 0,
        }
    }

    /// `Knots(kn)` (`cxx:606-619`): only the B-spline branch writes, and it
    /// writes the distinct knots.
    pub fn knots(&self) -> Vec<f64> {
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => self
                .curve
                .as_ref()
                .and_then(|c| c.as_bspline_curve())
                .map(|b| b.distinct_knots().to_vec())
                .unwrap_or_default(),
            _ => Vec::new(),
        }
    }

    /// `Multiplicities(mu)` (`cxx:621-633`).
    pub fn multiplicities(&self) -> Vec<i32> {
        match self.curve_kind_3d() {
            occt_topo::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve => self
                .curve
                .as_ref()
                .and_then(|c| c.as_bspline_curve())
                .map(|b| b.multiplicities().to_vec())
                .unwrap_or_default(),
            _ => Vec::new(),
        }
    }
}
