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
    /// `Geom_Curve::IsCN`. `N < 0` is the range error.
    fn is_cn(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_Curve::IsCN: negative order".into());
        }
        Ok(curve_smooth_is_cn(self.continuity(), n))
    }
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
    /// Basis of a trimmed curve for `IsKind(Geom_TrimmedCurve)` walks
    /// (`ShapeFix_Wire.cxx:1686-1697`). Default is the curve from
    /// [`Curve::untrimmed_basis`]. A trim that keeps basis parameters and must
    /// not advertise the `[0, 1]` remap overrides this alone.
    fn trim_basis(&self) -> Option<std::sync::Arc<dyn Curve>> {
        self.untrimmed_basis().map(|(c, _, _)| c)
    }
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
    /// `Geom_BezierCurve::Weights()` / `IsRational()`; `None` when polynomial.
    fn bezier_weights(&self) -> Option<&[f64]> { None }
    /// `Adaptor3d_Curve::NbIntervals`. Default one span.
    fn nb_intervals(&self, _continuity: u8) -> i32 { 1 }
    /// `AdvApprox_EvaluatorFunction::Evaluate` body for a 3D curve: ignore
    /// `StartEnd` and call `D0`/`D1`/`D2`. `CurveOnSurface` overrides this to
    /// re-trim when `StartEnd` changes (`GeomLib.cxx:1009-1016`).
    fn approx_eval(&self, _first: f64, _last: f64, u: f64, order: i32, out: &mut [f64]) -> i32 {
        match order {
            0 => {
                let p = self.d0(u);
                out[0] = p.x();
                out[1] = p.y();
                out[2] = p.z();
            }
            1 => {
                let (_, v) = self.d1(u);
                out[0] = v.x();
                out[1] = v.y();
                out[2] = v.z();
            }
            2 => {
                let (_, _, v) = self.d2(u);
                out[0] = v.x();
                out[1] = v.y();
                out[2] = v.z();
            }
            _ => {
                out[0] = 0.0;
                out[1] = 0.0;
                out[2] = 0.0;
                return 3;
            }
        }
        0
    }
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

    /// `Geom_Curve::TransformedParameter` (`Geom_Curve.cxx:49-52`): default
    /// identity. `Geom_Line` scales by `|ScaleFactor|` (`Geom_Line.cxx:225-232`).
    fn transformed_parameter(&self, u: f64, _t: &GpTrsf) -> f64 {
        u
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

fn curve_smooth_is_cn(shape: u8, n: i32) -> bool {
    match shape {
        6 => true,
        0 | 1 => n <= 0,
        2 | 3 => n <= 1,
        4 => n <= 2,
        5 => n <= 3,
        _ => false,
    }
}
