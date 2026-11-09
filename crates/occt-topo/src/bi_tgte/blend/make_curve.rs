//! `MakeCurve_Function` and `MakeCurve` (`BiTgte_Blend.cxx:381-470`): turn a
//! `BiTgte_CurveOnEdge` into the 3D curve the rolling ball traces.
//!
//! `MakeCurve` (`cxx:414-470`) has two arms:
//!
//! * a zero-radius `GeomAbs_Circle` adaptor is used analytically
//!   (`cxx:424-428`): `Geom_Circle` + `Geom_TrimmedCurve` over
//!   `[FirstParameter, LastParameter]`;
//! * otherwise the centre curve is approximated (`cxx:429-466`) with
//!   `Approx_FitAndDivide(F, 8, 8, Precision::Approximation(),
//!   Precision::Approximation(), true)` - the default
//!   `AppParCurves_TangencyPoint` constraints at both ends - then the resulting
//!   Bezier pieces are glued into one C1 B-spline by
//!   `Convert_CompBezierCurvesToBSplineCurve` and reparametrised from
//!   `[0, 1]` onto `[FirstParameter, LastParameter]`
//!   (`BSplCLib::Reparametrize`, `cxx:460`).
//!
//! `Approx_FitAndDivide` and `Approx_FitAndDivide2d` are the **same**
//! `Approx_ComputeCLine.gxx`, instantiated by
//! `Approx_FitAndDivide_0.cxx:20-26` and `Approx_FitAndDivide2d_0.cxx` with
//! `MultiLine = AppCont_Function`; the 3D class simply leaves the 2D tolerance
//! unused, which is what passing the same value to
//! [`FitAndDivide2d`] does here.
//!
//! `MakeCurve_Function` (`cxx:383-412`) reports one 3D point and no 2D point,
//! `Value` is `BiTgte_CurveOnEdge::Value` and `D1` answers `false`.

use std::sync::Arc;

use occt_core::bspl::knots::reparametrize;
use occt_core::convert::CompBezierCurves;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::APPROXIMATION;
use occt_geom::app_cont::function::Function;
use occt_geom::app_par_curves::constraint::Constraint;
use occt_geom::approx::fit_and_divide::FitAndDivide2d;
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom::{Curve, GeomBSplineCurve, GeomCircle};

use crate::bi_tgte::CurveOnEdge;
use crate::geom_abs_curve_type::GeomAbsCurveType;

/// `MakeCurve_Function` (`cxx:383-412`): the `AppCont_Function` view of a
/// `BiTgte_CurveOnEdge` handed to `Approx_FitAndDivide`.
struct MakeCurveFunction {
    /// `myCurve` (`cxx:386`).
    curve: CurveOnEdge,
}

impl MakeCurveFunction {
    /// `MakeCurve_Function(const BiTgte_CurveOnEdge&)` (`cxx:388-392`):
    /// `myNbPnt = 1`, `myNbPnt2d = 0`.
    fn new(curve: &CurveOnEdge) -> Self {
        Self {
            curve: curve.clone(),
        }
    }
}

impl Function for MakeCurveFunction {
    /// `GetNumberOfPoints` (`hxx:37-41`) with `myNbPnt = 1`, `myNbPnt2d = 0`
    /// (`cxx:390-391`).
    fn get_number_of_points(&self) -> (i32, i32) {
        (1, 0)
    }

    /// `FirstParameter()` (`cxx:394`).
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    /// `LastParameter()` (`cxx:396`).
    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    /// `Value(theT, thePnt2d, thePnt)` (`cxx:398-403`): `thePnt(1) =
    /// myCurve.Value(theT)`, always `true`.
    fn value(&self, the_u: f64, _the_pnt2d: &mut [GpPnt2d], the_pnt: &mut [GpPnt]) -> bool {
        the_pnt[0] = self.curve.d0(the_u);
        true
    }

    /// `D1(...)` (`cxx:405-411`): always `false`.
    fn d1(&self, _the_u: f64, _the_vec2d: &mut [GpVec2d], _the_vec: &mut [GpVec]) -> bool {
        false
    }
}

/// `MakeCurve(const BiTgte_CurveOnEdge&)` (`cxx:414-470`).
pub(super) fn make_curve(hc: &CurveOnEdge) -> Arc<dyn Curve> {
    if hc.get_type() == GeomAbsCurveType::Circle {
        // `cxx:424-428`: `Geom_Circle(HC.Circle())` trimmed on the adaptor's
        // range.
        let circle: Arc<dyn Curve> = Arc::new(GeomCircle::new(hc.circle()));
        Arc::new(GeomTrimmedCurve::new(
            circle,
            hc.first_parameter(),
            hc.last_parameter(),
        ))
    } else {
        // `cxx:429-466`.
        let f = MakeCurveFunction::new(hc);
        let fit = FitAndDivide2d::with_line(
            &f,
            8,
            8,
            APPROXIMATION,
            APPROXIMATION,
            true,
            Constraint::TangencyPoint,
            Constraint::TangencyPoint,
        );

        let nb_curves = fit.nb_multi_curves() as i32;
        // `cxx:437`: `Convert_CompBezierCurvesToBSplineCurve Conv` with the
        // header default angular tolerance `1.0e-4`.
        let mut conv = CompBezierCurves::new(1.0e-4);
        for i in 1..=nb_curves {
            let mc = fit.value(i);
            conv.add_curve(&mc.curve3d(1));
        }
        conv.perform();

        // `cxx:451-460`.
        let (knots, mults) = conv.knots_and_mults();
        let mut new_knots = knots.to_vec();
        let new_mults = mults.to_vec();
        reparametrize(hc.first_parameter(), hc.last_parameter(), &mut new_knots);
        let new_poles = conv.poles().to_vec();
        let degree = conv.degree() as usize;

        // `cxx:462`: `new Geom_BSplineCurve(NewPoles, NewKnots, NewMults,
        // Conv.Degree())`.
        Arc::new(
            GeomBSplineCurve::from_poles_knots_mults(new_poles, new_knots, new_mults, degree)
                .expect("Geom_BSplineCurve: invalid pole/knot/multiplicity data"),
        )
    }
}
