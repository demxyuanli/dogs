//! `Curve2d` view of `HLRBRep_Curve` - the bridge the `IntCurve` kernels need.
//!
//! `HLRBRep_CInter` is `IntCurve_IntCurveCurveGen<HLRBRep_CurveTool,
//! Geom2dIntGInter>` (`HLRBRep_CInter.hxx:26-34`); its `TheCurveTool` is a
//! static forwarder over `HLRBRep_Curve` (`HLRBRep_CurveTool.lxx:31-257`, ported
//! as [`crate::hlrbrep_curve_tool`]) and its `TheImpTool` is `Geom2dInt_GInter`.
//! The port's `Geom2dInt_GInter` takes `&dyn Curve2d`
//! (`crates/occt-geom2d/src/geom2d_int/ginter.rs:82-107`) rather than a tool
//! type parameter, so the two halves are joined by implementing `Curve2d` on
//! `HLRBRep_Curve` itself.
//!
//! Only the five `GeomAbs_CurveType` values the `IntCurve` dispatcher branches
//! on are mapped (`IntCurve_IntCurveCurveGen.gxx:249-812`): `Line`, `Circle`,
//! `Ellipse`, `Parabola`, `Hyperbola`. Every other type falls through to
//! `GeomAbsCurveType2d::OtherCurve`, which is the same `default:` arm
//! `HLRBRep_CurveTool::GetType` (`.lxx:173-176`) sends `BezierCurve`,
//! `BSplineCurve` and `OtherCurve` to.
//!
//! The type is read from `HLRBRep_Curve::GetType()` - the cached *projected*
//! type (`HLRBRep_Curve.lxx:138-141`) - not recomputed from the 2D samples, so
//! `Line()` / `Circle()` answer `None` whenever the projected type is something
//! else, mirroring `HLRBRep_CurveTool::Line` (`lxx:186-189`) which forwards to
//! `HLRBRep_Curve::Line()` unconditionally.
//!
//! UNREACHABLE in `IntCurve`: `reverse` and `transform`. `HLRBRep_Curve` has no
//! `Reverse` / `Transform` at all (`HLRBRep_Curve.hxx:44-213`), and the
//! `IntCurve` kernels never ask a `Curve2d` for either - they reverse the
//! intersection *result* (`IntRes2d_Intersection::SetReversedParameters`) and
//! transform `gp` conics, never the adaptor. They panic rather than silently
//! answer something OCCT cannot compute.

use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::hlrbrep_curve::HlrBRepCurve;

impl Curve2d for HlrBRepCurve {
    /// `HLRBRep_Curve::D0(U, P)` (`.lxx:145-151` through
    /// `HLRBRep_CurveTool::D0`, `.lxx:117-120`).
    fn d0(&self, u: f64) -> GpPnt2d {
        HlrBRepCurve::d0(self, u)
    }

    /// `HLRBRep_Curve::D1(U, P, T)` (`HLRBRep_CurveTool::D1`, `.lxx:124-131`).
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        HlrBRepCurve::d1(self, u)
    }

    /// `HLRBRep_Curve::D2(U, P, T, N)` (`HLRBRep_CurveTool::D2`, `.lxx:135-144`).
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        HlrBRepCurve::d2(self, u)
    }

    /// `HLRBRep_Curve::D3(U, P, V1, V2, V3)` (`HLRBRep_CurveTool::D3`,
    /// `.lxx:148-155`).
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        HlrBRepCurve::d3(self, u)
    }

    /// `HLRBRep_Curve::DN(U, N)` (`HLRBRep_CurveTool::DN`, `.lxx:159-162`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        HlrBRepCurve::dn(self, u, n)
    }

    /// `HLRBRep_CurveTool::Value` (`.lxx:110-113`).
    fn value(&self, u: f64) -> GpPnt2d {
        HlrBRepCurve::value(self, u)
    }

    /// `HLRBRep_Curve::FirstParameter` (`.lxx:64-67`).
    fn first_parameter(&self) -> f64 {
        HlrBRepCurve::first_parameter(self)
    }

    /// `HLRBRep_Curve::LastParameter` (`.lxx:71-74`).
    fn last_parameter(&self) -> f64 {
        HlrBRepCurve::last_parameter(self)
    }

    /// `HLRBRep_Curve::IsPeriodic` (`.lxx:109-112`).
    fn is_periodic(&self) -> bool {
        HlrBRepCurve::is_periodic(self)
    }

    /// `HLRBRep_Curve::Period` (`.lxx:116-119`).
    fn period(&self) -> f64 {
        HlrBRepCurve::period(self)
    }

    /// `HLRBRep_Curve::Continuity` (`.lxx:80-83`).
    fn continuity(&self) -> u8 {
        HlrBRepCurve::continuity(self)
    }

    /// `HLRBRep_Curve` has no `Transform` (`HLRBRep_Curve.hxx:44-213`).
    fn transform(&mut self, _t: &GpTrsf2d) {
        unreachable!(
            "HLRBRep_Curve has no Transform (HLRBRep_Curve.hxx:44-213); \
             IntCurve never transforms an adaptor"
        )
    }

    /// `HLRBRep_Curve` has no `Reverse` (`HLRBRep_Curve.hxx:44-213`).
    fn reverse(&mut self) {
        unreachable!(
            "HLRBRep_Curve has no Reverse (HLRBRep_Curve.hxx:44-213); \
             IntCurve never reverses an adaptor"
        )
    }

    /// `HLRBRep_CurveTool::NbIntervals` (`.lxx:56-63`): `NbIntervals(GeomAbs_C1)`.
    fn nb_intervals(&self, continuity: u8) -> i32 {
        HlrBRepCurve::nb_intervals(self, continuity)
    }

    /// `HLRBRep_Curve::Intervals(T, S)` (`.lxx:95-98`).
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        HlrBRepCurve::intervals(self, continuity)
    }

    /// `HLRBRep_CurveTool::GetType() == GeomAbs_Line` (`.lxx:173-176`).
    fn is_line(&self) -> bool {
        matches!(self.get_type(), GeomAbsCurveType::Line)
    }

    /// `HLRBRep_CurveTool::Line` (`.lxx:186-189`).
    fn gp_lin2d(&self) -> Option<GpLin2d> {
        match self.get_type() {
            GeomAbsCurveType::Line => HlrBRepCurve::line(self),
            _ => None,
        }
    }

    /// `HLRBRep_CurveTool::Circle` (`.lxx:192-195`).
    fn gp_circ2d(&self) -> Option<GpCirc2d> {
        match self.get_type() {
            GeomAbsCurveType::Circle => HlrBRepCurve::circle(self),
            _ => None,
        }
    }

    /// `HLRBRep_CurveTool::Ellipse` (`.lxx:198-201`).
    fn gp_elips2d(&self) -> Option<GpElips2d> {
        match self.get_type() {
            GeomAbsCurveType::Ellipse => HlrBRepCurve::ellipse(self),
            _ => None,
        }
    }

    /// `HLRBRep_CurveTool::Parabola` (`.lxx:204-207`).
    fn gp_parab2d(&self) -> Option<GpParab2d> {
        match self.get_type() {
            GeomAbsCurveType::Parabola => Some(HlrBRepCurve::parabola(self)),
            _ => None,
        }
    }

    /// `HLRBRep_CurveTool::Hyperbola` (`.lxx:211-214`).
    fn gp_hypr2d(&self) -> Option<GpHypr2d> {
        match self.get_type() {
            GeomAbsCurveType::Hyperbola => Some(HlrBRepCurve::hyperbola(self)),
            _ => None,
        }
    }

    /// `HLRBRep_CurveTool::NbSamples(C, u1, u2)` (`.lxx` forwards to
    /// `HLRBRep_CurveTool.cxx:53-83`): the HLR instantiation samples with its
    /// own cap of 50 and no circle-radius arm, so it overrides
    /// `Geom2dInt_Geom2dCurveTool`'s count
    /// (`Geom2dInt_Geom2dCurveTool.cxx:23-70`).
    fn nb_samples_2d(&self, u0: f64, u1: f64) -> usize {
        let n = crate::hlrbrep_curve_tool::nb_samples_range(self, u0, u1);
        if n < 1 {
            1
        } else {
            n as usize
        }
    }

    /// `HLRBRep_CurveTool::NbSamples(C)` (`HLRBRep_CurveTool.cxx:22-49`).
    fn nb_samples_2d_full(&self) -> usize {
        let n = crate::hlrbrep_curve_tool::nb_samples(self);
        if n < 1 {
            1
        } else {
            n as usize
        }
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
}
