//! `Geom2dGcc_Circ2d2TanOnGeo(gp_Pnt2d, gp_Pnt2d, OnCurv, Tolerance)`
//! (`Geom2dGcc_Circ2d2TanOnGeo.cxx:1064-1143`): the bisector of the two points
//! (the perpendicular line returned by `GccAna_Pnt2dBisec`) is intersected with
//! the carrier; every intersection point is a centre and the radius is the
//! distance to the second point.
//!
//! Translation notes:
//! - `Bis.HasSolution()` guards the intersection (`cxx:1111`); the port keeps
//!   that `if` so a coincident point pair answers `!IsDone()` with no solution,
//! - `pnttg1sol = Point1` and `pnttg2sol = Point2` (`cxx:1131-1132`) with both
//!   `pararg` zero (`cxx:1133-1134`),
//! - `WellDone = true` (`cxx:1140`) is inside `if (Intp.IsDone())` and outside
//!   `if (!Intp.IsEmpty())`.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::GpPnt2d;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::gcc_ana_bisec::GccAnaPnt2dBisec;
use crate::geom2d_int::IntConicCurveGen;

use super::{carrier_domain, circle_at, empty_domain, Geom2dGccCirc2d2TanOnGeo, MAX_SOL};

/// `Geom2dGcc_Circ2d2TanOnGeo::Geom2dGcc_Circ2d2TanOnGeo(gp_Pnt2d, gp_Pnt2d,
/// Geom2dAdaptor_Curve, double)` (`cxx:1064-1143`).
pub(super) fn new(
    point1: &GpPnt2d,
    point2: &GpPnt2d,
    on_curve: &Geom2dAdaptorCurve,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanOnGeo {
    let mut s = Geom2dGccCirc2d2TanOnGeo::empty();
    let tol = tolerance.abs();
    let bis = GccAnaPnt2dBisec::new(point1, point2);
    if bis.is_done() {
        let tol1 = tol;
        let tol2 = tol1;
        let mut inter = IntConicCurveGen::new();
        let (d2, curve) = carrier_domain(on_curve, tol);
        let d1 = empty_domain();
        if bis.has_solution() {
            let line = bis
                .this_solution()
                .expect("GccAna_Pnt2dBisec::ThisSolution");
            inter.perform_line(&line, &d1, curve.as_ref(), &d2, tol1, tol2);
            if inter.is_done() {
                if !inter.result().is_empty() {
                    for j in 1..=inter.result().nb_points() {
                        if s.nbr_sol >= MAX_SOL {
                            break;
                        }
                        let center = *inter.result().point(j).value();
                        let radius = point2.distance(&center);
                        s.nbr_sol += 1;
                        let n = s.nbr_sol;
                        s.cirsol[n - 1] = circle_at(center, radius);
                        s.qualifier1[n - 1] = GccEntPosition::NoQualifier;
                        s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                        s.the_same1[n - 1] = 0;
                        s.the_same2[n - 1] = 0;
                        s.pntcen[n - 1] = center;
                        s.pnttg1sol[n - 1] = *point1;
                        s.pnttg2sol[n - 1] = *point2;
                        s.pararg1[n - 1] = 0.0;
                        s.pararg2[n - 1] = 0.0;
                        s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], point1);
                        s.par2sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], point2);
                        s.parcen3[n - 1] = inter.result().point(j).param_on_second();
                    }
                }
                s.well_done = true;
            }
        }
    }
    s
}
