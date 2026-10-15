//! [`ReparamCurve2d`] and the two conversion tails `GeomLib::SameRange`
//! (`GeomLib.cxx:842-970`) is built from.

use super::*;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

/// A 2D curve linearly reparameterized from `[p0, p1]` onto `[a, b]`: a parameter
/// `t ∈ [a, b]` maps to `inner.d0(p0 + (t − a)·(p1 − p0)/(b − a))`. Aligns a
/// STEP-imported pcurve (parameterized over its natural range) with the edge's
/// 3D range.
///
/// This is a stand-in, not an OCCT class. `GeomLib::SameRange` only ever stores
/// the input curve unchanged, a translated `Geom2d_Line` (`cxx:864-871`), a
/// rotated `Geom2d_Circle` (`:872-889`), a re-trimmed `Geom2d_TrimmedCurve`
/// (`:890-901`) or a concrete `Geom2d_BSplineCurve` (`:912-921`, `:962-968`),
/// and throws `Standard_DomainError` for a kind
/// `Geom2dConvert::CurveToBSplineCurve` cannot convert
/// (`Geom2dConvert.cxx:442-444`). [`reparam_curve2d`] builds one only where the
/// faithful conversion declined, so the callers that consume `GeomLib::SameRange`
/// (`shhealing::wire_fix::pcurve_range::geom_lib_same_range`) hand back OCCT's
/// concrete curve and this type stays off the shared pipeline.
#[derive(Clone)]
pub struct ReparamCurve2d {
    pub(crate) inner: Arc<dyn Curve2d>,
    pub(crate) a: f64,
    pub(crate) b: f64,
    pub(crate) p0: f64,
    pub(crate) p1: f64,
}

impl ReparamCurve2d {
    pub fn new(inner: Arc<dyn Curve2d>, a: f64, b: f64, p0: f64, p1: f64) -> Self {
        Self { inner, a, b, p0, p1 }
    }

    pub(crate) fn scale(&self) -> f64 {
        if (self.b - self.a).abs() < 1e-30 {
            1.0
        } else {
            (self.p1 - self.p0) / (self.b - self.a)
        }
    }

    pub(crate) fn map(&self, t: f64) -> f64 {
        self.p0 + (t - self.a) * self.scale()
    }
}

impl Curve2d for ReparamCurve2d {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.inner.d0(self.map(u))
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p, d) = self.inner.d1(self.map(u));
        let s = self.scale();
        (p, GpVec2d::new(d.x() * s, d.y() * s))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1, d2) = self.inner.d2(self.map(u));
        let s = self.scale();
        (
            p,
            GpVec2d::new(d1.x() * s, d1.y() * s),
            GpVec2d::new(d2.x() * s * s, d2.y() * s * s),
        )
    }
    fn first_parameter(&self) -> f64 {
        self.a
    }
    fn last_parameter(&self) -> f64 {
        self.b
    }
    fn continuity(&self) -> u8 {
        self.inner.continuity()
    }
    fn transform(&mut self, t: &GpTrsf2d) {
        let mut inner = self.inner.clone_dyn();
        inner.transform(t);
        self.inner = Arc::from(inner);
    }
    fn reverse(&mut self) {
        let mut inner = self.inner.clone_dyn();
        inner.reverse();
        self.inner = Arc::from(inner);
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
    // Preserve the inner `GeomAbs_Line` classification through the wrapper so
    // `Adaptor3d_CurveOnSurface::EvalKPart` (`cxx:1552-1732`) and the
    // `GeomBndLib_Curve2d` kind dispatch keep iso/line detection. A bare
    // default `is_line=false` densifies Offset seams.
    fn is_line(&self) -> bool {
        self.inner.is_line()
    }
    fn gp_lin2d(&self) -> Option<occt_core::gp::GpLin2d> {
        self.inner.gp_lin2d()
    }
    // Kind delegation. `Geom2dAdaptor_Curve::load` unwraps the stored curve
    // before classifying it (`Geom2dAdaptor_Curve.cxx:226-230`), so the
    // wrapper's sample count (`Geom2dAdaptor_Curve::NbSamples`,
    // `Geom2dAdaptor_Curve.cxx:1351-1394`) and bounding-box sample budget
    // (`Geom2dAdaptor_Curve::GetType`, mirrored by
    // `geom_bnd_lib_sample2d::sample_kind_of`) must come from the inner curve.
    // A bare trait default makes a reparameterized B-spline pcurve report
    // `GeomAbs_OtherCurve` and fall back to 20 samples.
    fn gp_circ2d(&self) -> Option<occt_core::gp::GpCirc2d> {
        self.inner.gp_circ2d()
    }
    fn bezier_nb_poles(&self) -> Option<usize> {
        self.inner.bezier_nb_poles()
    }
    fn bspline_nb_knots(&self) -> Option<usize> {
        self.inner.bspline_nb_knots()
    }
    fn bspline_degree(&self) -> Option<usize> {
        self.inner.bspline_degree()
    }
    // No `bspline_copy2d` override: `GeomBndLib_Curve2d.cxx:137-143` builds the
    // `GeomBndLib_BSplineCurve2d` evaluator from a concrete
    // `Geom2d_BSplineCurve`, and a wrapper only exists where the faithful
    // `Geom2dConvert::CurveToBSplineCurve` declined, so the trait default
    // (`None`, which selects a different evaluator) is the correct answer.
    fn offset_basis(&self) -> Option<&dyn Curve2d> {
        self.inner.offset_basis()
    }
    /// `Geom2dAdaptor_Curve::Resolution` in this wrapper's parameter space.
    ///
    /// For the kinds `GeomLib::SameRange` stores unconverted the parameterization
    /// is kept -- a line is only translated (`GeomLib.cxx:864-871`) and a circle
    /// only rotated (`:872-889`), so `scale = ±1` and the inner answer carries
    /// over unchanged. Every other kind is converted to a `Geom2d_BSplineCurve`
    /// and `BSplCLib::Reparametrize`d onto `[a, b]` (`:912-921`, `:962-968`),
    /// which divides a parameter increment by the `scale` factor of this map.
    /// A wrapper exists only where that conversion declined, so this override is
    /// what keeps its `Resolution` -- the `Precision::Parametric` default would
    /// under-step by orders of magnitude on a reparameterized conic or Bezier
    /// pcurve.
    fn resolution(&self, ruv: f64) -> f64 {
        let s = self.scale();
        let res = self.inner.resolution(ruv);
        if s.abs() <= 1e-30 {
            res
        } else {
            res / s.abs()
        }
    }
}

/// The `[Udeb, Ufin]` window the unequal-span arm trims to
/// (`GeomLib.cxx:934-959`).
///
/// The periodicity test reads `BasisCurve()` when the input is a
/// `Geom2d_TrimmedCurve` (`cxx:929-932`, `aCCheck`): a periodic basis keeps the
/// requested `[FirstOnCurve, LastOnCurve]` window (`cxx:936-938`) and falls back
/// to the curve's own range when it collapses (`cxx:940-943`), every other kind
/// takes `[max(FirstParameter, FirstOnCurve), min(LastParameter, LastOnCurve)]`
/// (`cxx:948-949`) and falls back the same way (`cxx:950-957`, both under
/// `Precision::PConfusion()`).
fn same_range_window(inner: &dyn Curve2d, first_on: f64, last_on: f64) -> (f64, f64) {
    let pconfusion = occt_core::precision::PCONFUSION;
    let a_check: &dyn Curve2d = inner.trimmed_basis().unwrap_or(inner);
    if a_check.is_periodic() {
        if (last_on - first_on).abs() > pconfusion {
            (first_on, last_on)
        } else {
            (inner.first_parameter(), inner.last_parameter())
        }
    } else {
        let (cf, cl) = (inner.first_parameter(), inner.last_parameter());
        // `std::max`/`std::min` with `+-inf` yield the finite operand, so the
        // bounded-curve guard only keeps the port off `inf - inf`.
        let udeb = if cf.is_finite() { cf.max(first_on) } else { first_on };
        let ufin = if cl.is_finite() { cl.min(last_on) } else { last_on };
        if (ufin - udeb).abs() > pconfusion {
            (udeb, ufin)
        } else {
            (cf, cl)
        }
    }
}

/// The conversion tail both `SameRange` arms share: `Geom2d_TrimmedCurve(CurvePtr,
/// U1, U2)` (`cxx:912-913`, `:952`), `Geom2dConvert::CurveToBSplineCurve`
/// (`cxx:915`, `:962`), then `BSplCLib::Reparametrize(RequestedFirst,
/// RequestedLast, Knots)` over `BS->Knots()` and `BS->SetKnots(Knots)`
/// (`cxx:918-920`, `:965-967`).
///
/// `Knots()` is the **distinct** knot array: `Geom2d_BSplineCurve::Knots()`
/// returns `myKnots`, which holds one entry per `Multiplicities()` slot
/// (`Geom2d_BSplineCurve_1.cxx`), and `BSplCLib::Reparametrize` maps
/// `Knots(Knots.Lower())`/`Knots(Knots.Upper())` onto
/// `[RequestedFirst, RequestedLast]` (`BSplCLib.cxx:756-798`) -- the ported
/// `reparametrize` is that same function on that same array.
///
/// `None` where the port has no equivalent branch: OCCT's
/// `CurveToBSplineCurve` throws `Standard_DomainError` there
/// (`Geom2dConvert.cxx:442-444`), and an infinite bound is only reachable
/// through the degenerate `SetTrim` fallback above (`cxx:954-958`), where
/// OCCT builds a trimmed curve with infinite bounds.
pub(crate) fn trimmed_reparam_bspline(
    curve: Arc<dyn Curve2d>,
    u1: f64,
    u2: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Geom2dBSplineCurve> {
    if !u1.is_finite() || !u2.is_finite() {
        return None;
    }
    let tc = Geom2dTrimmedCurve::new(curve, u1, u2);
    let mut bs = occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
        &tc,
        occt_geom2d::geom2d_convert::DEFAULT_PARAMETERISATION,
    )?;
    let (mut distinct, _) = bs.distinct_knots_and_mults();
    occt_core::bspl::knots::reparametrize(req_first, req_last, &mut distinct);
    bs.set_knots(&distinct).ok()?;
    Some(bs)
}

/// The concrete `Geom2d_BSplineCurve` of the unequal-span arm
/// (`GeomLib.cxx:924-969`): resolve the window ([`same_range_window`]) then run
/// the shared tail ([`trimmed_reparam_bspline`]). OCCT dispatches on no curve
/// kind here -- a line takes the same path as a conic or a B-spline.
fn same_range_bspline2d(
    inner: &dyn Curve2d,
    a: f64,
    b: f64,
    p0: f64,
    p1: f64,
) -> Option<Geom2dBSplineCurve> {
    let (u_deb, u_fin) = same_range_window(inner, p0, p1);
    trimmed_reparam_bspline(Arc::from(inner.clone_dyn()), u_deb, u_fin, a, b)
}

/// The curve `GeomLib::SameRange` (`GeomLib.cxx:842-970`) hands back for a
/// `[p0, p1]`-parameterized curve requested on `[a, b]`, for its unequal-span
/// arm (`:924-969`).
///
/// Every kind `Geom2dConvert::CurveToBSplineCurve` can convert -- line,
/// trimmed line, circle, ellipse, hyperbola, parabola, Bezier, B-spline and
/// `Geom2d_OffsetCurve` (`Geom2dConvert.cxx:181-449`, including the `RationalC1`
/// circle/ellipse split through `CompCurveToBSplineCurve` at `:237-263`/`:276-302`
/// and the offset approximation at `:353-368`/`:425-440`) -- yields the concrete
/// reparameterized `Geom2d_BSplineCurve` OCCT stores (`cxx:962-968`), so the
/// dispatched callers see the class OCCT tests for with
/// `STANDARD_TYPE(Geom2d_BSplineCurve)`: `Geom2dConvert::CurveToBSplineCurve`
/// (`Geom2dConvert.cxx:420-423`), `GeomLib::To3d` (`GeomLib.cxx:603-632`),
/// `GeomLib::isIsoLine` (`GeomLib.cxx:3010-3029`, on the `GeomAbs_BSplineCurve`
/// `GetType()`) and `ShapeBuild_Edge::TransformPCurve`
/// (`ShapeBuild_Edge.cxx:679-698`).
///
/// Only a kind the conversion declines keeps the linear [`ReparamCurve2d`]
/// stand-in, which reproduces the OCCT parameterization for sampling but not the
/// concrete class. OCCT cannot reach that state from `GeomLib::SameRange`: its
/// conversion throws and leaves the out-parameter null, which
/// `BRepLib::SameRange` then stores as the pcurve (`BRepLib.cxx:236-254`); the
/// port has no null pcurve, so that throw maps to this wrapper.
pub fn reparam_curve2d(
    inner: Arc<dyn Curve2d>,
    a: f64,
    b: f64,
    p0: f64,
    p1: f64,
) -> Arc<dyn Curve2d> {
    if let Some(bs) = same_range_bspline2d(inner.as_ref(), a, b, p0, p1) {
        return Arc::new(bs);
    }
    Arc::new(ReparamCurve2d::new(inner, a, b, p0, p1))
}

/// Parameter of the 2D curve closest to `uv` (the 1-D projection of a surface
/// point onto the pcurve). A bounded curve is searched on `[first, last]`; an
/// unbounded line is projected in closed form.
pub fn project_uv_on_curve2d(curve: &dyn Curve2d, uv: GpPnt2d) -> f64 {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    if f.is_finite() && l.is_finite() {
        let mut best = f;
        let mut best_d = f64::INFINITY;
        for i in 0..=64 {
            let u = f + (l - f) * i as f64 / 64.0;
            let d = curve.d0(u).distance(&uv);
            if d < best_d {
                best_d = d;
                best = u;
            }
        }
        let mut step = (l - f) / 64.0;
        for _ in 0..32 {
            step *= 0.5;
            let mut u = best;
            let mut dd = best_d;
            for &du in &[-step, step] {
                let nu = (u + du).clamp(f, l);
                let d = curve.d0(nu).distance(&uv);
                if d < dd {
                    dd = d;
                    u = nu;
                }
            }
            best = u;
            best_d = dd;
        }
        best
    } else {
        let p0 = curve.d0(0.0);
        let dir = curve.d0(1.0).coord.subtracted(&p0.coord);
        let len2 = dir.square_modulus();
        if len2 < 1e-30 {
            return 0.0;
        }
        uv.coord.subtracted(&p0.coord).dot(&dir) / len2
    }
}

