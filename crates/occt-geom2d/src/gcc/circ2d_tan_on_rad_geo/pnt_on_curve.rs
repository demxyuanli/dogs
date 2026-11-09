//! `Geom2dGcc_Circ2dTanOnRadGeo(gp_Pnt2d, Adaptor, Radius, Tolerance)`
//! (`Geom2dGcc_Circ2dTanOnRadGeo.cxx:736-816`): the centre locus is the circle
//! of the requested radius around the operand point, intersected with the
//! carrier.
//!
//! Translation notes:
//! - no qualifier reaches this constructor, so `qualifier1` is
//!   `GccEnt_noqualifier` and `pnttg1sol` is the operand point itself
//!   (`cxx:805`, `:809`),
//! - `par1sol` is computed before `pnttg1sol` is written (`cxx:808-809`),
//! - OCCT builds the intersector with the six-argument constructor
//!   `Geom2dInt_TheIntConicCurveOfGInter Intp(Circ, D1, OnCurv, D2, Tol, Tol);`
//!   (`cxx:793`), the same call as `Perform`.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::GpPnt2d;
use occt_core::intres2d::IntRes2dDomain;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::geom2d_int::IntConicCurveGen;

use super::{
    circle_at, circle_domain, Geom2dGccCirc2dTanOnRadGeo, MAX_SOL, THEFIRST, THELAST,
};

/// `Geom2dGcc_Circ2dTanOnRadGeo(gp_Pnt2d, Geom2dAdaptor_Curve, double, double)`
/// (`cxx:736-816`).
pub(super) fn new(
    point1: &GpPnt2d,
    on_curve: &Geom2dAdaptorCurve,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2dTanOnRadGeo {
    let mut s = Geom2dGccCirc2dTanOnRadGeo::empty();
    let tol = tolerance.abs();
    if radius < 0.0 {
        panic!("Standard_NegativeValue in Geom2dGcc_Circ2dTanOnRadGeo");
    }
    let circ = circle_at(*point1, radius);
    let d1 = circle_domain(&circ, tol);
    let carrier = on_curve
        .curve()
        .expect("Geom2dAdaptor_Curve::Curve on the carrier")
        .clone();
    let first = curve_tool::first_parameter(on_curve).max(THEFIRST);
    let last = curve_tool::last_parameter(on_curve).min(THELAST);
    let p_first =
        curve_tool::value(on_curve, first).expect("Geom2dAdaptor_Curve::Value on the carrier");
    let p_last =
        curve_tool::value(on_curve, last).expect("Geom2dAdaptor_Curve::Value on the carrier");
    let d2 = IntRes2dDomain::bounded(&p_first, first, tol, &p_last, last, tol);
    let mut inter = IntConicCurveGen::new();
    inter.perform_circle(&circ, &d1, carrier.as_ref(), &d2, tol, tol);
    if inter.is_done() {
        if !inter.result().is_empty() {
            for i in 1..=inter.result().nb_points() {
                if s.nbr_sol >= MAX_SOL {
                    break;
                }
                s.nbr_sol += 1;
                let n = s.nbr_sol;
                let point = inter.result().point(i);
                let center = *point.value();
                let par_on_first = point.param_on_first();
                let par_on_second = point.param_on_second();
                s.cirsol[n - 1] = circle_at(center, radius);
                s.qualifier1[n - 1] = GccEntPosition::NoQualifier;
                s.the_same1[n - 1] = 0;
                s.pararg1[n - 1] = par_on_first;
                s.parcen3[n - 1] = par_on_second;
                s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg1sol[n - 1]);
                s.pnttg1sol[n - 1] = *point1;
                s.pntcen3[n - 1] = center;
            }
        }
        s.well_done = true;
    }
    s
}
