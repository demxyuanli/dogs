//! `Geom2dGcc_Circ2d2TanOnGeo(QualifiedLin, QualifiedLin, OnCurv, Tolerance)`
//! (`Geom2dGcc_Circ2d2TanOnGeo.cxx:537-705`): the bisector of the two lines is
//! intersected with the carrier; every accepted centre gives the radius
//! `(dist1 + dist2) / 2`.
//!
//! Translation notes:
//! - `Radius` is declared before the loops (`cxx:568`) and only assigned inside
//!   the three `Qualified2` arms (`cxx:636`, `:645`, `:650`); a candidate that
//!   reaches `if (ok)` through the `Qualified2.IsEnclosing()` / unqualified
//!   path keeps the previous radius, which the port reproduces,
//! - `TheSame1` / `TheSame2` are unconditionally zero (`cxx:684-685`),
//! - `WellDone = true` (`cxx:702`) is inside `if (Intp.IsDone())` and outside
//!   `if (!Intp.IsEmpty())`.

use occt_core::gp::{GpDir2d, GpPnt2d};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use crate::gcc_ana_bisec::GccAnaLin2dBisec;
use crate::geom2d_int::IntConicCurveGen;

use super::{
    carrier_domain, circle_at, dir_from, empty_domain, par_on_circle, par_on_line, signed_normal,
    Geom2dGccCirc2d2TanOnGeo, MAX_SOL,
};

/// `Geom2dGcc_Circ2d2TanOnGeo::Geom2dGcc_Circ2d2TanOnGeo(QualifiedLin,
/// QualifiedLin, Geom2dAdaptor_Curve, double)` (`cxx:537-705`).
pub(super) fn new(
    qualified1: &GccEntQualifiedLin,
    qualified2: &GccEntQualifiedLin,
    on_curve: &Geom2dAdaptorCurve,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanOnGeo {
    let mut s = Geom2dGccCirc2d2TanOnGeo::empty();
    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
        || !(qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnGeo");
    }
    let tol = tolerance.abs();
    let mut radius = 0.0_f64;
    let l1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let dir1 = *l1.direction();
    let dir2 = *l2.direction();
    let dnor1 = GpDir2d::new(-dir1.y(), dir1.x()).expect("gp_Dir2d of the line normal");
    let dnor2 = GpDir2d::new(-dir2.y(), dir2.x()).expect("gp_Dir2d of the line normal");
    let origin1 = l1.location();
    let origin2 = l2.location();
    let bis = GccAnaLin2dBisec::new(&l1, &l2);
    if bis.is_done() {
        let tol1 = tol;
        let tol2 = tol1;
        let mut inter = IntConicCurveGen::new();
        let nbsolution = bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions");
        let (d2, curve) = carrier_domain(on_curve, tol);
        let d1 = empty_domain();
        for i in 1..=nbsolution {
            let line = bis
                .this_solution(i)
                .expect("GccAna_Lin2dBisec::ThisSolution");
            inter.perform_line(&line, &d1, curve.as_ref(), &d2, tol1, tol2);
            if inter.is_done() {
                if !inter.result().is_empty() {
                    for j in 1..=inter.result().nb_points() {
                        if s.nbr_sol >= MAX_SOL {
                            break;
                        }
                        let center = *inter.result().point(j).value();
                        let dist1 = l1.distance(&center);
                        let dist2 = l2.distance(&center);
                        let mut ok = false;
                        if qualified1.is_enclosed() {
                            if (origin1.x() - center.x()) * -dir1.y()
                                + (origin1.y() - center.y()) * dir1.x()
                                <= 0.0
                            {
                                ok = true;
                            }
                        } else if qualified1.is_outside() {
                            if (origin1.x() - center.x()) * -dir1.y()
                                + (origin1.y() - center.y()) * dir1.x()
                                >= 0.0
                            {
                                ok = true;
                            }
                        } else if qualified1.is_unqualified() {
                            ok = true;
                        }
                        if qualified2.is_enclosed() && ok {
                            ok = false;
                            if (origin2.x() - center.x()) * -dir2.y()
                                + (origin2.y() - center.y()) * dir2.x()
                                <= 0.0
                            {
                                ok = true;
                                radius = (dist1 + dist2) / 2.0;
                            }
                        } else if qualified2.is_outside() && ok {
                            ok = false;
                            if (origin2.x() - center.x()) * -dir2.y()
                                + (origin2.y() - center.y()) * dir2.x()
                                >= 0.0
                            {
                                ok = true;
                                radius = (dist1 + dist2) / 2.0;
                            }
                        } else if qualified2.is_unqualified() && ok {
                            radius = (dist1 + dist2) / 2.0;
                        }
                        if ok {
                            s.nbr_sol += 1;
                            let n = s.nbr_sol;
                            s.cirsol[n - 1] = circle_at(center, radius);
                            let dc1 = dir_from(&origin1, &center);
                            let dc2 = dir_from(&origin2, &center);
                            s.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                qualified1.qualifier()
                            } else if dc1.dot(&dnor1) > 0.0 {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosed
                            };
                            s.qualifier2[n - 1] = if !qualified2.is_unqualified() {
                                qualified2.qualifier()
                            } else if dc2.dot(&dnor2) > 0.0 {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosed
                            };
                            s.the_same1[n - 1] = 0;
                            s.the_same2[n - 1] = 0;
                            let sign = dc1.dot(&dnor1);
                            let dc1 = signed_normal(sign, &dir1);
                            let p1 = GpPnt2d::new(
                                center.x() + radius * dc1.x(),
                                center.y() + radius * dc1.y(),
                            );
                            let par1 = par_on_circle(&s.cirsol[n - 1], &p1);
                            let arg1 = par_on_line(&l1, &p1);
                            s.pnttg1sol[n - 1] = p1;
                            s.par1sol[n - 1] = par1;
                            s.pararg1[n - 1] = arg1;
                            let sign = dc2.dot(&dnor2);
                            let dc2 = signed_normal(sign, &dir2);
                            let p2 = GpPnt2d::new(
                                center.x() + radius * dc2.x(),
                                center.y() + radius * dc2.y(),
                            );
                            let par2 = par_on_circle(&s.cirsol[n - 1], &p2);
                            let arg2 = par_on_line(&l2, &p2);
                            s.pnttg2sol[n - 1] = p2;
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
