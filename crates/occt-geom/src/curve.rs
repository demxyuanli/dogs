//! Abstract 3D parametric curve. Source: `Geom_Curve.hxx`
use occt_core::gp::{GpPnt, GpVec, GpTrsf};

pub trait Curve: Send + Sync {
    fn d0(&self, u: f64) -> GpPnt;
    fn d1(&self, u: f64) -> (GpPnt, GpVec);
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec);
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (p, d1, d2) = self.d2(u);
        (p, d1, d2, GpVec::zero())
    }
    /// `Geom_Curve::EvalDN` (`Geom_Curve.hxx:210`, pure virtual).
    /// Illegal `N < 1` returns zero instead of throw. `GeomBSplineCurve`
    /// overrides via `BSplCLib::DN`.
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        if n < 1 {
            return GpVec::zero();
        }
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            _ => GpVec::zero(),
        }
    }
    fn value(&self, u: f64) -> GpPnt { self.d0(u) }
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    fn is_periodic(&self) -> bool { false }
    fn period(&self) -> f64 { 0.0 }
    fn continuity(&self) -> u8;
    /// Radius when this is a `Geom_Circle`; `None` for every other type.
    /// Source: `Adaptor3d_Curve::GetType() == GeomAbs_Circle`.
    fn circle_radius(&self) -> Option<f64> { None }
    /// `Geom_Line::Lin` / `Adaptor3d_Curve::Line`. `None` otherwise.
    fn gp_line(&self) -> Option<occt_core::gp::GpLin> { None }
    /// `Geom_Circle::Circ` / `Adaptor3d_Curve::Circle`. `None` otherwise.
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> { None }
    /// `Geom_Ellipse::Elips` / `Adaptor3d_Curve::Ellipse`. `None` otherwise.
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> { None }
    /// `Geom_Hyperbola::Hypr` / `Adaptor3d_Curve::Hyperbola`. `None` otherwise.
    fn gp_hyperbola(&self) -> Option<occt_core::gp::GpHypr> { None }
    /// `Geom_Parabola::Parab` / `Adaptor3d_Curve::Parabola`. `None` otherwise.
    fn gp_parabola(&self) -> Option<occt_core::gp::GpParab> { None }
    /// `GeomAdaptor_Curve::GetType() == GeomAbs_OffsetCurve` companion:
    /// the basis curve of a `Geom_OffsetCurve` and its offset value
    /// (`Geom_OffsetCurve::BasisCurve()` / `Offset()`). `None` for every other
    /// curve type; used by `IntTools_EdgeEdge::ResolutionCoeff`/`Resolution`.
    fn offset_curve(&self) -> Option<(std::sync::Arc<dyn Curve>, f64)> { None }
    /// `Geom_TrimmedCurve`. OCCT `IsKind(STANDARD_TYPE(Geom_TrimmedCurve))`.
    fn is_geom_trimmed(&self) -> bool { false }
    /// Basis `[First, Last]` of a `Geom_TrimmedCurve` before the `[0, 1]` remap.
    /// `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps to the basis curve;
    /// a `Geom_Circle` keeps radian parameters on that interval.
    fn trimmed_basis_range(&self) -> Option<(f64, f64)> { None }
    /// Basis curve of a `Geom_TrimmedCurve` plus its basis parameter range.
    /// `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps the trim and keeps the
    /// basis curve, so `GeomBndLib_Curve` sees the basis with the basis range.
    fn untrimmed_basis(&self) -> Option<(std::sync::Arc<dyn Curve>, f64, f64)> { None }
    /// `Adaptor3d_Curve::GetType() == GeomAbs_Line`.
    fn is_line(&self) -> bool { false }
    /// `Adaptor3d_Curve::Degree` for Bezier / BSpline.
    fn nurbs_degree(&self) -> Option<usize> { None }
    /// `Adaptor3d_Curve::Intervals` including the range ends (`GeomAbs_CN`).
    fn parameter_intervals(&self, _continuity: u8) -> Vec<f64> {
        vec![self.first_parameter(), self.last_parameter()]
    }
    /// Flat knot sequence of a `Geom_BSplineCurve`; `None` otherwise.
    fn bspline_knots(&self) -> Option<&[f64]> {
        None
    }
    /// Poles when this is a `Geom_BSplineCurve`; `None` otherwise.
    /// Source: `Adaptor3d_Curve::GetType() == GeomAbs_BSplineCurve`.
    fn bspline_poles(&self) -> Option<&[GpPnt]> { None }
    /// `Geom_BSplineCurve::Weights()`; `None` when the curve is not rational or
    /// not a B-spline. Needed by the STEP writer, whose
    /// `GeomToStep_MakeBSplineCurveWithKnotsAndRationalBSplineCurve`
    /// (`GeomToStep_MakeBoundedCurve.cxx:52-62`) branches on `IsRational`.
    fn bspline_weights(&self) -> Option<&[f64]> { None }
    /// Poles when this is a `Geom_BezierCurve`; `None` otherwise.
    /// Source: `Adaptor3d_Curve::GetType() == GeomAbs_BezierCurve`.
    fn bezier_poles(&self) -> Option<&[GpPnt]> { None }
    /// `Adaptor3d_Curve::NbIntervals`. Default one span.
    fn nb_intervals(&self, _continuity: u8) -> i32 { 1 }
    /// `GeomAdaptor_Curve::Resolution` (`cxx:1116-1148`).
    /// Default `Precision::Parametric(r3d) = r3d * PConfusion / Confusion`.
    fn resolution(&self, r3d: f64) -> f64 {
        r3d * occt_core::precision::PCONFUSION / occt_core::precision::CONFUSION
    }
    fn transform(&mut self, t: &GpTrsf);
    fn reverse(&mut self);

    /// `Geom_Curve::ReversedParameter` — the parameter on the *reversed* curve
    /// of the point at parameter `U`. OCCT implements it per class:
    /// `Geom_Line`/`Geom_Hyperbola`/`Geom_Parabola` -> `-U`
    /// (`Geom_Line.cxx:72`, `Geom_Hyperbola.cxx:167`, `Geom_Parabola.cxx:84`),
    /// `Geom_Circle`/`Geom_Ellipse` -> `2*pi - U` (`Geom_Circle.cxx:87`,
    /// `Geom_Ellipse.cxx:164`), which equals the affine form on their
    /// `[0, 2*pi]` range, `Geom_BSplineCurve` -> `first + last - U`
    /// (`Geom_BSplineCurve.cxx:520`), `Geom_BezierCurve` -> `1 - U`
    /// (`Geom_BezierCurve.cxx:381`, on the always-`[0, 1]` range), and
    /// `Geom_TrimmedCurve` / `Geom_OffsetCurve` -> the basis curve's
    /// (`Geom_TrimmedCurve.cxx:88`, `Geom_OffsetCurve.cxx:104`). The default
    /// here is that shared affine form.
    fn reversed_parameter(&self, u: f64) -> f64 {
        self.first_parameter() + self.last_parameter() - u
    }

    fn clone_dyn(&self) -> Box<dyn Curve>;
    fn transformed(&self, t: &GpTrsf) -> Box<dyn Curve> { let mut c = self.clone_dyn(); c.transform(t); c }
    fn reversed(&self) -> Box<dyn Curve> { let mut c = self.clone_dyn(); c.reverse(); c }
    /// `Geom_Curve::Translate`.
    fn translated(&self, v: &GpVec) -> Box<dyn Curve> {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(v);
        self.transformed(&t)
    }
}
