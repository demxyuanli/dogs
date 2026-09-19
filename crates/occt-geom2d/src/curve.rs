//! Abstract 2D parametric curve trait. Source: `Geom2d_Curve.hxx`
use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d, GpVec2d, GpTrsf2d};

/// Parametric 2D curve. Replaces OCCT Geom2d_Curve.
pub trait Curve2d: Send + Sync {
    fn d0(&self, u: f64) -> GpPnt2d;
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d);
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d);
    fn value(&self, u: f64) -> GpPnt2d { self.d0(u) }
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    fn is_periodic(&self) -> bool { false }
    fn period(&self) -> f64 { 0.0 }
    fn continuity(&self) -> u8;
    fn transform(&mut self, t: &GpTrsf2d);
    fn reverse(&mut self);

    fn clone_dyn(&self) -> Box<dyn Curve2d>;

    /// `Geom2dAdaptor_Curve::Intervals` break points for the requested shape.
    /// Default is a single span `[First, Last]`.
    fn parameter_intervals(&self, _continuity: u8) -> Vec<f64> {
        vec![self.first_parameter(), self.last_parameter()]
    }

    /// `Geom2dAdaptor_Curve::NbIntervals`.
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }

    /// `Geom2dAdaptor_Curve::GetType() == GeomAbs_Line`.
    fn is_line(&self) -> bool {
        false
    }
    /// `Geom2dAdaptor_Curve::Line`. `None` otherwise.
    fn gp_lin2d(&self) -> Option<GpLin2d> {
        None
    }
    /// `Geom2dAdaptor_Curve::GetType() == GeomAbs_Circle` then `Circle()`.
    fn gp_circ2d(&self) -> Option<GpCirc2d> {
        None
    }

    /// Basis of a `Geom2d_TrimmedCurve`; `None` for every other type.
    fn trimmed_basis(&self) -> Option<&dyn Curve2d> {
        None
    }

    /// Basis of a `Geom2d_OffsetCurve` (`Geom2d_OffsetCurve::BasisCurve()`,
    /// `Geom2d_OffsetCurve.cxx:174-177`); `None` for every other type. Mirrors
    /// `IsKind(STANDARD_TYPE(Geom2d_OffsetCurve))` then `BasisCurve()` in
    /// `Geom2dAdaptor_Curve.cxx:1373-1377`.
    fn offset_basis(&self) -> Option<&dyn Curve2d> {
        None
    }

    /// `Geom2d_BezierCurve::NbPoles()` (`Geom2d_BezierCurve.cxx:600`); `Some`
    /// only when this curve is a `Geom2d_BezierCurve`, mirroring
    /// `IsKind(STANDARD_TYPE(Geom2d_BezierCurve))` (`Geom2dAdaptor_Curve.cxx:1360-1363`).
    fn bezier_nb_poles(&self) -> Option<usize> {
        None
    }

    /// `Geom2d_BSplineCurve::NbKnots()` (`Geom2d_BSplineCurve_1.cxx:598`); `Some`
    /// only when this curve is a `Geom2d_BSplineCurve`, mirroring
    /// `IsKind(STANDARD_TYPE(Geom2d_BSplineCurve))` (`Geom2dAdaptor_Curve.cxx:1364-1368`).
    fn bspline_nb_knots(&self) -> Option<usize> {
        None
    }

    /// `Geom2d_BSplineCurve::Degree()` (`Geom2d_BSplineCurve_1.cxx:168`); `Some`
    /// only when this curve is a `Geom2d_BSplineCurve`
    /// (`Geom2dAdaptor_Curve.cxx:1367`).
    fn bspline_degree(&self) -> Option<usize> {
        None
    }

    /// (x, y) poles of a non-rational `Geom2d_BSplineCurve`; `None` otherwise.
    /// Feeds the control-net branch of `GeomBndLib_BSplineCurve2d::Box`
    /// (`GeomBndLib_BSplineCurve2d.cxx:52-57`).
    fn bspline_poles2d(&self) -> Option<(&[f64], &[f64])> {
        None
    }

    fn transformed(&self, t: &GpTrsf2d) -> Box<dyn Curve2d> {
        let mut c = self.clone_dyn();
        c.transform(t);
        c
    }
    fn reversed(&self) -> Box<dyn Curve2d> {
        let mut c = self.clone_dyn();
        c.reverse();
        c
    }
}

// ponytail: Arc<dyn Curve2d> works for cloning without custom Clone impls
