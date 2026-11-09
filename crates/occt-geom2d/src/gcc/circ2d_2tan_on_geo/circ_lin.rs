//! `Geom2dGcc_Circ2d2TanOnGeo(QualifiedCirc, QualifiedLin, OnCurv, Tolerance)`
//! (`Geom2dGcc_Circ2d2TanOnGeo.cxx:325-522`): the bisector of the circle and
//! the line (a line or a parabola, `cxx:386-403`) is intersected with the
//! carrier; the radius is the distance from the centre to the line.
//!
//! Translation notes:
//! - `if (Qualified1.IsEnclosing() && dist1 > Radius) ok = false;`
//!   (`cxx:456-459`) runs after the `Qualified2` test and is the only place the
//!   enclosing qualifier rejects a candidate,
//! - `WellDone = true` (`cxx:519`) is inside `if (Intp.IsDone())` but outside
//!   `if (!Intp.IsEmpty())`, so an empty intersection still marks the solver
//!   done with no solution,
//! - `TheSame2` is unconditionally zero (`cxx:507`).

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::{GpDir2d, GpPnt2d};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc_ana_bisec::GccAnaCircLin2dBisec;
use crate::gcc_int::GccIntBisec;
use crate::geom2d_int::IntConicCurveGen;

use super::{
    carrier_domain, circle_at, dir_from, empty_domain, par_on_circle, par_on_line, parabola_domain,
    signed_normal, Geom2dGccCirc2d2TanOnGeo, MAX_SOL,
};

/// `Geom2dGcc_Circ2d2TanOnGeo::Geom2dGcc_Circ2d2TanOnGeo(QualifiedCirc,
/// QualifiedLin, Geom2dAdaptor_Curve, double)` (`cxx:325-522`).
pub(super) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedLin,
    on_curve: &Geom2dAdaptorCurve,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanOnGeo {
    let mut s = Geom2dGccCirc2d2TanOnGeo::empty();
    let tol = tolerance.abs();
    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
        || !(qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnGeo");
    }
    let c1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let origin2 = l2.location();
    let dir2 = *l2.direction();
    let norm_l2 = GpDir2d::new(-dir2.y(), dir2.x()).expect("gp_Dir2d of the line normal");
    let bis = GccAnaCircLin2dBisec::new(&c1, &l2);
    if bis.is_done() {
        let tol1 = tol;
        let tol2 = tol1;
        let mut inter = IntConicCurveGen::new();
        let nbsolution = bis.nb_solutions();
        let (d2, curve) = carrier_domain(on_curve, tol);
        for i in 1..=nbsolution {
            let sol = bis
                .this_solution(i)
                .expect("GccAna_CircLin2dBisec::ThisSolution");
            match sol {
                GccIntBisec::Line(line) => {
                    let d1 = empty_domain();
                    inter.perform_line(&line, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Parabola(parab) => {
                    let d1 = parabola_domain(&parab, tol1);
                    inter.perform_parabola(&parab, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Circle(_)
                | GccIntBisec::Ellipse(_)
                | GccIntBisec::Hyperbola(_)
                | GccIntBisec::Point(_) => {
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
                        let dist1 = center.distance(&center1);
                        let mut ok = false;
                        if qualified1.is_enclosed() {
                            if dist1 - r1 < tol {
                                ok = true;
                            }
                        } else if qualified1.is_outside() {
                            if r1 - dist1 < tol {
                                ok = true;
                            }
                        } else if qualified1.is_enclosing() || qualified1.is_unqualified() {
                            ok = true;
                        }
                        let radius = l2.distance(&center);
                        if qualified2.is_enclosed() && ok {
                            ok = false;
                            if (origin2.x() - center.x()) * -dir2.y()
                                + (origin2.y() - center.y()) * dir2.x()
                                <= 0.0
                            {
                                ok = true;
                            }
                        } else if qualified2.is_outside() && ok {
                            ok = false;
                            if (origin2.x() - center.x()) * -dir2.y()
                                + (origin2.y() - center.y()) * dir2.x()
                                >= 0.0
                            {
                                ok = true;
                            }
                        }
                        if qualified1.is_enclosing() && dist1 > radius {
                            ok = false;
                        }
                        if ok {
                            s.nbr_sol += 1;
                            let n = s.nbr_sol;
                            s.cirsol[n - 1] = circle_at(center, radius);
                            let dc2 = dir_from(&origin2, &center);
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
                            s.qualifier2[n - 1] = if !qualified2.is_unqualified() {
                                qualified2.qualifier()
                            } else if dc2.dot(&norm_l2) > 0.0 {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosed
                            };
                            if dist1 <= tol && (radius - c1.radius()).abs() <= tol {
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
                            let sign = dc2.dot(&GpDir2d::new(-dir2.y(), dir2.x())
                                .expect("gp_Dir2d of the line normal"));
                            let dc2 = signed_normal(sign, &dir2);
                            let p = GpPnt2d::new(
                                center.x() + radius * dc2.x(),
                                center.y() + radius * dc2.y(),
                            );
                            let par2 = par_on_circle(&s.cirsol[n - 1], &p);
                            let arg2 = par_on_line(&l2, &p);
                            s.pnttg2sol[n - 1] = p;
                            s.par2sol[n - 1] = par2;
                            s.pararg2[n - 1] = arg2;
                            s.pntcen[n - 1] = center;
                            s.parcen3[n - 1] = inter.result().point(j).param_on_second();
                        }
                    }
                }
                s.well_done = true;
            }
        }
    }
    s
}


