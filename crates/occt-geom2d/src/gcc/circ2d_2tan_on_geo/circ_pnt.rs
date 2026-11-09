//! `Geom2dGcc_Circ2d2TanOnGeo(QualifiedCirc, gp_Pnt2d, OnCurv, Tolerance)`
//! (`Geom2dGcc_Circ2d2TanOnGeo.cxx:716-899`): the bisector of the circle and
//! the point (a circle, a line, an ellipse or a hyperbola branch, `cxx:774-816`)
//! is intersected with the carrier; the radius is the distance from the centre
//! to the point.
//!
//! Translation notes:
//! - the qualifier tests use `<= Tol` here (`cxx:831`, `:840`) where the
//!   circle/circle constructor uses `< Tol`,
//! - `par1sol` is written twice (`cxx:886-887`, `0.` then
//!   `ElCLib::Parameter(cirsol, pnttg1sol)`); only the second survives, so the
//!   port writes once,
//! - `par2sol` is assigned after `pnttg2sol = Point2` (`cxx:893`), and
//!   `pararg2` is zero (`cxx:891`),
//! - `WellDone = true` (`cxx:897`) is inside `if (Intp.IsDone())` and outside
//!   `if (!Intp.IsEmpty())`.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::GpPnt2d;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc_ana_bisec::GccAnaCircPnt2dBisec;
use crate::gcc_int::GccIntBisec;
use crate::geom2d_int::IntConicCurveGen;

use super::{
    carrier_domain, circle_at, circle_domain, dir_from, ellipse_domain, empty_domain,
    hyperbola_domain, par_on_circle, Geom2dGccCirc2d2TanOnGeo, MAX_SOL,
};

/// `Geom2dGcc_Circ2d2TanOnGeo::Geom2dGcc_Circ2d2TanOnGeo(QualifiedCirc,
/// gp_Pnt2d, Geom2dAdaptor_Curve, double)` (`cxx:716-899`).
pub(super) fn new(
    qualified1: &GccEntQualifiedCirc,
    point2: &GpPnt2d,
    on_curve: &Geom2dAdaptorCurve,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanOnGeo {
    let mut s = Geom2dGccCirc2d2TanOnGeo::empty();
    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnGeo");
    }
    let tol = tolerance.abs();
    let c1 = qualified1.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let bis = GccAnaCircPnt2dBisec::new(&c1, point2);
    if bis.is_done() {
        let tol1 = tol;
        let tol2 = tol1;
        let mut inter = IntConicCurveGen::new();
        let nbsolution = bis.nb_solutions();
        let (d2, curve) = carrier_domain(on_curve, tol);
        for i in 1..=nbsolution {
            let sol = bis
                .this_solution(i)
                .expect("GccAna_CircPnt2dBisec::ThisSolution");
            match sol {
                GccIntBisec::Circle(circ) => {
                    let d1 = circle_domain(&circ, tol1, tol2);
                    inter.perform_circle(&circ, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Line(line) => {
                    let d1 = empty_domain();
                    inter.perform_line(&line, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Ellipse(elips) => {
                    let d1 = ellipse_domain(&elips, tol1, tol2);
                    inter.perform_ellipse(&elips, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Hyperbola(hypr) => {
                    let d1 = hyperbola_domain(&hypr, tol1, tol2);
                    inter.perform_hyperbola(&hypr, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Parabola(_) | GccIntBisec::Point(_) => {
                    panic!("Standard_ConstructionError in Geom2dGcc_Circ2d2TanOnGeo");
                }
            }
            if inter.is_done() {
                if !inter.result().is_empty() {
                    for j in 1..=inter.result().nb_points() {
                        if s.nbr_sol >= MAX_SOL {
                            break;
                        }
                        let center = *inter.result().point(j).value();
                        let radius = center.distance(point2);
                        let dist1 = center1.distance(&center);
                        let mut ok = false;
                        if qualified1.is_enclosed() {
                            if dist1 - r1 <= tol {
                                ok = true;
                            }
                        } else if qualified1.is_outside() {
                            if r1 - dist1 <= tol {
                                ok = true;
                            }
                        } else if qualified1.is_enclosing() {
                            ok = true;
                        } else if qualified1.is_unqualified() {
                            ok = true;
                        }
                        if ok {
                            s.nbr_sol += 1;
                            let n = s.nbr_sol;
                            s.cirsol[n - 1] = circle_at(center, radius);
                            let distcc1 = center.distance(&center1);
                            s.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                qualified1.qualifier()
                            } else if (distcc1 + radius - r1).abs() < tol {
                                GccEntPosition::Enclosed
                            } else if (distcc1 - r1 - radius).abs() < tol {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosing
                            };
                            s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                            if dist1 <= tol && (radius - r1).abs() <= tol {
                                s.the_same1[n - 1] = 1;
                            } else {
                                s.the_same1[n - 1] = 0;
                                let dc1 = dir_from(&center1, &center);
                                let p = GpPnt2d::new(
                                    center.x() + radius * dc1.x(),
                                    center.y() + radius * dc1.y(),
                                );
                                let par1 = par_on_circle(&s.cirsol[n - 1], &p);
                                let arg1 = parameter_circ2d(&c1, &p);
                                s.pnttg1sol[n - 1] = p;
                                s.par1sol[n - 1] = par1;
                                s.pararg1[n - 1] = arg1;
                            }
                            s.the_same2[n - 1] = 0;
                            s.pnttg2sol[n - 1] = *point2;
                            s.pntcen[n - 1] = center;
                            s.parcen3[n - 1] = inter.result().point(j).param_on_second();
                            s.pararg2[n - 1] = 0.0;
                            s.par2sol[n - 1] = par_on_circle(&s.cirsol[n - 1], point2);
                        }
                    }
                }
                s.well_done = true;
            }
        }
    }
    s
}
