//! Abstract 2D parametric curve trait. Source: `Geom2d_Curve.hxx`
use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d, GpTrsf2d};

use crate::Geom2dBSplineCurve;

/// Parametric 2D curve. Replaces OCCT Geom2d_Curve.
pub trait Curve2d: Send + Sync {
    fn d0(&self, u: f64) -> GpPnt2d;
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d);
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d);
    /// `Geom2d_Curve::EvalD3`. Default returns a zero third derivative, mirroring
    /// the 3D `Curve::d3` default; overridden by `Geom2dCircle`/`Geom2dEllipse`/
    /// `Geom2dHyperbola` (`clib2d.rs:294-328`) for `Geom2d_OffsetCurve::d2`
    /// (task T-63). The remaining implementors keep the zero default — see the
    /// `eval_dn` note below.
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        let (p, d1, d2) = self.d2(u);
        (p, d1, d2, GpVec2d::zero())
    }
    /// `Geom2d_Curve::EvalDN` (`Geom2d_Curve.hxx:210`, pure virtual).
    /// The default covers `N = 1..3` through `d1`/`d2`/`d3`; `N < 1` (OCCT throws
    /// `Geom2d_UndefinedDerivative`) and `N > 3` return a zero vector — the same
    /// convention as the 3D `Curve::eval_dn` (`occt-geom/src/curve.rs:12-25`).
    /// UNPORTED for `N > 3`: the 2D B-spline/Bezier/trimmed/offset implementors
    /// do not override this yet (3D `GeomBSplineCurve` uses `BSplCLib::DN`; the
    /// 2D `Geom2d_BSplineCurve::EvalDN` / `Geom2d_BezierCurve::EvalDN` have no
    /// counterpart here). The elementary curves *do* override it via
    /// `ElCLib::*DN` (`clib2d.rs:329-432`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        if n < 1 {
            return GpVec2d::zero();
        }
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            _ => GpVec2d::zero(),
        }
    }
    fn value(&self, u: f64) -> GpPnt2d { self.d0(u) }
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    fn is_periodic(&self) -> bool { false }
    fn period(&self) -> f64 { 0.0 }
    fn continuity(&self) -> u8;
    fn transform(&mut self, t: &GpTrsf2d);
    fn reverse(&mut self);

    /// `Geom2d_Curve::ReversedParameter` — the parameter on the *reversed*
    /// curve of the point at parameter `U`. OCCT implements it per class:
    /// `Geom2d_Line` → `-U` (`Geom2d_Line.cxx:135`),
    /// `Geom2d_Circle` → `2*pi - U` (`Geom2d_Circle.cxx:122`),
    /// `Geom2d_BSplineCurve` → `first + last - U`
    /// (`Geom2d_BSplineCurve.cxx:700`), `Geom2d_TrimmedCurve` → the basis
    /// curve's (`Geom2d_TrimmedCurve.cxx:91`). The default here is the affine
    /// form used by the polynomial classes.
    fn reversed_parameter(&self, u: f64) -> f64 {
        self.first_parameter() + self.last_parameter() - u
    }

    /// `Geom2d_Curve::TransformedParameter` (`Geom2d_Curve.cxx:41-44`): the
    /// parameter on `me->Transformed(T)` of the point at `U`. Base returns `U`;
    /// `Geom2d_Line`/`Geom2d_Parabola` scale it by `|T.ScaleFactor()|`, and
    /// `Geom2d_TrimmedCurve`/`Geom2d_OffsetCurve` delegate to the basis curve.
    fn transformed_parameter(&self, u: f64, _t: &GpTrsf2d) -> f64 {
        u
    }

    /// `Geom2d_BSplineCurve::Pole` / `Geom2d_BezierCurve::Pole`: the control
    /// points, for `ShapeBuild_Edge::TransformPCurve`'s affinity branch
    /// (`ShapeBuild_Edge.cxx:642-698`). `None` for curves that are not
    /// pole-based.
    fn poles2d(&self) -> Option<Vec<GpPnt2d>> {
        None
    }

    /// Setter counterpart of [`Curve2d::poles2d`].
    fn set_poles2d(&mut self, _poles: &[GpPnt2d]) {}

    /// `STANDARD_TYPE(Geom2d_BSplineCurve)` — for `Geom2dConvert::CurveToBSplineCurve`.
    fn is_bspline2d(&self) -> bool {
        false
    }

    /// `STANDARD_TYPE(Geom2d_BezierCurve)`.
    fn is_bezier2d(&self) -> bool {
        false
    }

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
    /// `Geom2dAdaptor_Curve::GetType() == GeomAbs_Ellipse` then `Ellipse()`.
    /// Mirrors the 3D `Curve::gp_ellipse` (`GeomAdaptor_Curve.cxx:266-270`).
    fn gp_elips2d(&self) -> Option<GpElips2d> {
        None
    }
    /// `Geom2dAdaptor_Curve::GetType() == GeomAbs_Parabola` then `Parabola()`.
    fn gp_parab2d(&self) -> Option<GpParab2d> {
        None
    }
    /// `Geom2dAdaptor_Curve::GetType() == GeomAbs_Hyperbola` then `Hyperbola()`.
    fn gp_hypr2d(&self) -> Option<GpHypr2d> {
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

    /// `Geom2d_OffsetCurve::Offset()` (`Geom2d_OffsetCurve.cxx:169-172`); the
    /// companion of [`Curve2d::offset_basis`]. Needed by `GeomLib::To3d`'s
    /// `Geom2d_OffsetCurve` arm (`GeomLib.cxx:575-581`).
    fn offset_value(&self) -> Option<f64> {
        None
    }

    /// `Geom2d_BSplineCurve::Weights()`. `None` when the curve is not a
    /// B-spline or is not rational. `Geom2dBSplineCurve` reports its stored
    /// weights; it exists so that `GeomLib::To3d`'s B-spline arm
    /// (`GeomLib.cxx:598-620`) can be written the way OCCT writes it.
    fn bspline_weights2d(&self) -> Option<&[f64]> {
        None
    }

    /// `Geom2d_BezierCurve::Weights()` (`GeomLib.cxx:582-597`).
    fn bezier_weights2d(&self) -> Option<&[f64]> {
        None
    }

    /// `Geom2d_BSplineCurve::Knots()` / `Multiplicities()`: the **distinct**
    /// knots and their multiplicities, i.e. what `GeomLib::To3d`'s B-spline arm
    /// feeds to the `Geom_BSplineCurve(Poles, Knots, Mults, Degree, IsPeriodic)`
    /// constructor (`GeomLib.cxx:610-618`). Mirror of
    /// `GeomBSplineCurve::distinct_knots_and_mults`; `None` for every other type.
    fn bspline_distinct_knots_mults(&self) -> Option<(Vec<f64>, Vec<i32>)> {
        None
    }

    /// `Geom2d_BezierCurve::NbPoles()` (`Geom2d_BezierCurve.cxx:600`); `Some`
    /// only when this curve is a `Geom2d_BezierCurve`, mirroring
    /// `IsKind(STANDARD_TYPE(Geom2d_BezierCurve))` (`Geom2dAdaptor_Curve.cxx:1360-1363`).
    fn bezier_nb_poles(&self) -> Option<usize> {
        None
    }

    /// `Geom2dAdaptor_Curve::IsRational()` (`Geom2dAdaptor_Curve.cxx:1291-1302`):
    /// the Bezier and B-spline arms return the stored curve's `IsRational()`,
    /// and every other `GeomAbs_CurveType` returns `false` (`:1299-1300`).
    ///
    /// UNPORTED (outside this crate): the `occt-topo` wrappers
    /// `pcurve::PlaneKeepParam2d` and `pcurve_full::ReparamCurve2d` delegate
    /// `bezier_nb_poles` but keep this default, so a rational curve behind them
    /// still reports `false`.
    fn is_rational(&self) -> bool {
        false
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

    /// `Geom2d_BSplineCurve::Copy()` as fed to
    /// `GeomBndLib_BSplineCurve2d::Box`'s `Segment` arm
    /// (`GeomBndLib_BSplineCurve2d.cxx:45-49`): the owned B-spline the evaluator
    /// was built from (`GeomBndLib_Curve2d.cxx:137-143`), which for a
    /// `Geom2d_TrimmedCurve` is its `BasisCurve()` (`:150-157`). `None` for
    /// every non-`Geom2d_BSplineCurve` — OCCT then picks a different evaluator.
    fn bspline_copy2d(&self) -> Option<Geom2dBSplineCurve> {
        None
    }

    /// The flat knot sequence of a `Geom2d_BSplineCurve` (`Geom2d_BSplineCurve::
    /// Knots()`, the `Knot(j)` values `ShapeAnalysis_TransferParametersProj::
    /// CorrectParameter` snaps onto, `Proj.cxx:268-279`); `None` otherwise.
    fn bspline_knots2d(&self) -> Option<&[f64]> {
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
