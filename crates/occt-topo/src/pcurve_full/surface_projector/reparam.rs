//! [`ReparamCurve2d`]: the linear 2D reparameterization that aligns a STEP pcurve natural range with the edge 3D range, and its helpers.

use super::*;

/// A 2D curve linearly reparameterized from `[p0, p1]` onto `[a, b]`: a parameter
/// `t ∈ [a, b]` maps to `inner.d0(p0 + (t − a)·(p1 − p0)/(b − a))`. Aligns a
/// STEP-imported pcurve (parameterized over its natural range) with the edge's
/// 3D range.
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
    // Preserve Geom2d_Line kind through the wrapper so
    // `Adaptor3d_CurveOnSurface::EvalKPart` (`cxx:1552-1732`) and
    // `GeomLib::SameRange` Line translate (`cxx:862-870`) keep iso/line
    // detection. A bare default `is_line=false` densifies Offset seams.
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
    // The concrete `Geom2d_BSplineCurve` `GeomLib::SameRange` would have
    // stored; see [`same_range_bspline2d`]. `GeomBndLib_Curve2d.cxx:137-143`
    // builds the `GeomBndLib_BSplineCurve2d` evaluator from exactly that curve,
    // so without this the wrapper falls through to the `OtherCurve2d` sampler
    // (`GeomBndLib_BSplineCurve2d.cxx:31-59` never runs).
    fn bspline_copy2d(&self) -> Option<Geom2dBSplineCurve> {
        same_range_bspline2d(self.inner.as_ref(), self.a, self.b, self.p0, self.p1)
    }
    fn offset_basis(&self) -> Option<&dyn Curve2d> {
        self.inner.offset_basis()
    }
    /// `Geom2dAdaptor_Curve::Resolution` in this wrapper's parameter space.
    ///
    /// `GeomLib::SameRange` keeps the parameterization for a line (it only
    /// translates the line, `GeomLib.cxx:866-870`) and for a circle (it only
    /// rotates it, `:871-889`), so for those `scale = ±1` and the inner answer
    /// carries over unchanged. Every other kind is converted to a
    /// `Geom2d_BSplineCurve` and `BSplCLib::Reparametrize`d onto `[a, b]`
    /// (`:908-921`, `:924-969`), which divides a parameter increment by the
    /// `scale` factor of this map. Without this override the wrapper would fall
    /// back to the `Precision::Parametric` default and under-step by orders of
    /// magnitude on a reparameterized conic or Bezier pcurve.
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

/// The `Geom2d_BSplineCurve` `GeomLib::SameRange` returns for a B-spline input
/// whose `[p0, p1]` range is remapped onto `[a, b]`.
///
/// OCCT's `SameRange` stores a concrete curve, never a wrapper. Both the
/// equal-span arm (`GeomLib.cxx:908-921`) and the unequal-span arm
/// (`GeomLib.cxx:924-969`) end the same way:
///
/// * clamp the kept window to the curve's own range — `:947-949`
///   `Udeb = max(FirstParameter(), FirstOnCurve)`,
///   `Ufin = min(LastParameter(), LastOnCurve)`, and keep the whole curve when
///   `|Ufin - Udeb| <= PConfusion` (`:951-957`); the same clamp sits in
///   `Geom2dConvert.cxx:199-209`;
/// * `Geom2dConvert::CurveToBSplineCurve` (`Geom2dConvert.cxx:347-351` for a
///   trimmed input, `:420-423` for a plain one): copy the basis and `Segment`
///   the copy to `[U1, U2]`;
/// * `BSplCLib::Reparametrize(a, b, Knots)` + `SetKnots` (`:961-969`).
///
/// `None` when the port has no equivalent branch. The conic inputs are handled
/// through [`occt_geom2d::geom2d_convert::trimmed_curve_to_bspline_curve`]
/// (`Convert_CircleToBSplineCurve` / `Convert_EllipseToBSplineCurve` /
/// `Convert_HyperbolaToBSplineCurve` / `Convert_ParabolaToBSplineCurve`,
/// `Geom2dConvert.cxx:228-321`); an offset curve needs
/// `Geom2dConvert_ApproxCurve` (`:353-368`, `:425-440`) and a `RationalC1`
/// circle/ellipse longer than 6 rad needs
/// `Geom2dConvert_CompCurveToBSplineCurve` (`:237-263`, `:276-302`) — those
/// stay on the linear [`ReparamCurve2d`] stand-in.
pub(crate) fn same_range_bspline2d(
    inner: &dyn Curve2d,
    a: f64,
    b: f64,
    p0: f64,
    p1: f64,
) -> Option<Geom2dBSplineCurve> {
    let pconfusion = occt_core::precision::PCONFUSION;
    // A conic basis (`Geom2d_Circle` / `Geom2d_Ellipse` / hyperbola / parabola)
    // exhausts the conic arm of `Geom2dConvert::CurveToBSplineCurve`
    // (`Geom2dConvert.cxx:228-321`): `GeomLib::SameRange` wraps the curve in a
    // `Geom2d_TrimmedCurve` (`GeomLib.cxx:912-915`, `:926-962`), converts it
    // (default `Convert_TgtThetaOver2`) and reparameterizes the knots
    // (`:961-969`).
    if inner.gp_circ2d().is_some()
        || inner.gp_elips2d().is_some()
        || inner.gp_hypr2d().is_some()
        || inner.gp_parab2d().is_some()
    {
        // `Geom2dConvert.cxx:199-209`: clamp to the basis range only when the
        // basis is not periodic (`GeomLib.cxx:934-958`).
        let (mut u_deb, mut u_fin) = if inner.is_periodic() {
            (p0, p1)
        } else {
            (
                if inner.first_parameter().is_finite() {
                    inner.first_parameter().max(p0)
                } else {
                    p0
                },
                if inner.last_parameter().is_finite() {
                    inner.last_parameter().min(p1)
                } else {
                    p1
                },
            )
        };
        if (u_fin - u_deb).abs() <= pconfusion {
            u_deb = inner.first_parameter();
            u_fin = inner.last_parameter();
        }
        let mut bs = occt_geom2d::geom2d_convert::trimmed_curve_to_bspline_curve(
            inner,
            u_deb,
            u_fin,
            occt_geom2d::geom2d_convert::DEFAULT_PARAMETERISATION,
        )?;
        occt_core::bspl::knots::reparametrize(a, b, &mut bs.knots);
        return Some(bs);
    }
    let mut a_curve = inner.bspline_copy2d()?;
    let (a_u1, a_u2) = (a_curve.first_parameter(), a_curve.last_parameter());
    let (mut u_deb, mut u_fin) = (a_u1.max(p0), a_u2.min(p1));
    if (u_fin - u_deb).abs() <= pconfusion {
        u_deb = a_u1;
        u_fin = a_u2;
    }
    if (u_deb - a_u1).abs() > pconfusion || (u_fin - a_u2).abs() > pconfusion {
        a_curve.segment(u_deb, u_fin, pconfusion).ok()?;
    }
    occt_core::bspl::knots::reparametrize(a, b, &mut a_curve.knots);
    Some(a_curve)
}

/// The curve `GeomLib::SameRange` (`GeomLib.cxx:842-969`) hands back for a
/// `[p0, p1]`-parameterized curve requested on `[a, b]`.
///
/// A B-spline input yields the concrete segmented/reparameterized
/// `Geom2d_BSplineCurve` OCCT stores (`:961-969`), so the port's dispatched
/// callers see the same class OCCT tests for with
/// `STANDARD_TYPE(Geom2d_BSplineCurve)`:
/// `Geom2dConvert::CurveToBSplineCurve` (`Geom2dConvert.cxx:420-423`),
/// `GeomLib::To3d` (`GeomLib.cxx:603-632`),
/// `GeomLib::isIsoLine` (`GeomLib.cxx:3010-3029`, on the
/// `GeomAbs_BSplineCurve` `GetType()`) and
/// `ShapeBuild_Edge::TransformPCurve` (`ShapeBuild_Edge.cxx:679-698`).
///
/// Every other kind keeps the linear [`ReparamCurve2d`] stand-in, which
/// reproduces the OCCT parameterization for sampling but not the concrete
/// class.
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

