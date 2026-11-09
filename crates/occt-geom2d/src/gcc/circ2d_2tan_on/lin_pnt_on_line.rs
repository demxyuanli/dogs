//! `GccAna_Circ2d2TanOn(QualifiedLin, gp_Pnt2d, gp_Lin2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_4.cxx:51-241`): circles tangent to a line and passing
//! through a point, with the centre on a straight line.
//!
//! `cxx:89-96` returns no solution at all when `L1` and `OnLine` are the same
//! line. The boundary block (`cxx:98-150`) projects `Point2` on `OnLine`. The
//! general case (`cxx:157-240`) intersects the single
//! `GccAna_LinPnt2dBisec` branch of `(L1, Point2)` with `OnLine`.
//!
//! Translation notes:
//! - `cxx:89-90` uses `Precision::Confusion()` in both tests; the port uses
//!   [`occt_core::precision::CONFUSION`],
//! - `cxx:89-96` leaves `WellDone = false` and `NbrSol = 0` and returns, so the
//!   object reports "not done"; the port returns the empty constructor result,
//! - `cxx:223` re-projects `dc1` on the same normal it just dotted, so `sign`
//!   and `dc1` duplicate `dc2`/`normal`; the port keeps the separate locals.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::precision::CONFUSION;

use crate::gcc_ana_bisec::GccAnaLinPnt2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use super::{
    axis_at, dir_between, left_normal, offset_xy, qualifier_of_line, side, GccAnaCirc2d2TanOn,
    MAX_SOL,
};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    point2: &GpPnt2d,
    on_line: &GpLin2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();

    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let tol = tolerance.abs();
    let l1 = qualified1.qualified();
    let origin_l1 = l1.location();
    let dir_l1 = *l1.direction();
    let normal = left_normal(&dir_l1);

    // `cxx:89-96`: identical lines give no solution and no `WellDone`.
    if dir_l1.is_equal(on_line.direction(), CONFUSION)
        && on_line.distance(&origin_l1) < CONFUSION
    {
        return me;
    }

    // `cxx:98-150`: boundary case.
    let dp2l = on_line.distance(point2);
    let on_norm = left_normal(on_line.direction());
    let mut pinterm = offset_xy(point2, dp2l, &on_norm);
    if on_line.distance(&pinterm) > tol {
        pinterm = offset_xy(point2, -dp2l, &on_norm);
    }
    let dist = l1.distance(&pinterm);
    if (dist - dp2l).abs() <= tol {
        let dirbid = dir_between(&pinterm, &origin_l1);
        if qualified1.is_enclosed() && dirbid.dot(&normal) < 0.0 {
            me.well_done = true;
        } else if qualified1.is_outside() && dirbid.dot(&normal) > 0.0 {
            me.well_done = true;
        } else if qualified1.is_unqualified() {
            me.well_done = true;
        }
        if me.well_done {
            me.nbr_sol = 1;
            let n = me.nbr_sol;
            me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), dp2l);
            me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
            let dc2 = dir_between(&pinterm, &origin_l1);
            me.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                qualified1.qualifier()
            } else if dc2.dot(&normal) > 0.0 {
                GccEntPosition::Outside
            } else {
                GccEntPosition::Enclosed
            };
            let sign = dc2.dot(&normal);
            let dc2 = signed_normal(sign, &normal);
            me.pnttg1sol[n - 1] = offset_xy(&pinterm, dp2l, &dc2);
            me.pnttg2sol[n - 1] = *point2;
            me.pntcen[n - 1] = pinterm;
            me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
            me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
            me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
            me.pararg2[n - 1] = 0.0;
            me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
            return me;
        }
    }

    // `cxx:157-240`: general case.
    let bis = GccAnaLinPnt2dBisec::new(&l1, point2);
    if bis.is_done() {
        let sol = bis
            .this_solution()
            .expect("GccAna_LinPnt2dBisec::ThisSolution");
        let typ = sol.arc_type();
        let mut intp = IntAna2dAnaIntersection::new();
        if typ == GccIntIType::Lin {
            intp.perform_lin_lin(on_line, &sol.line().unwrap());
        }
        if typ == GccIntIType::Par {
            intp.perform_lin_conic(
                on_line,
                &IntAna2dConic::from_parab2d(&sol.parabola().unwrap()),
            );
        }
        if intp.is_done() {
            if !intp.is_empty() {
                for j in 1..=intp.nb_points() {
                    let center = *intp.point(j).value();
                    let radius = l1.distance(&center);
                    let mut ok = false;
                    if qualified1.is_enclosed() {
                        if side(&origin_l1, &center, &normal) <= 0.0 {
                            ok = true;
                        }
                    } else if qualified1.is_outside() {
                        if side(&origin_l1, &center, &normal) >= 0.0 {
                            ok = true;
                        }
                    } else if qualified1.is_unqualified() {
                        ok = true;
                    }
                    if ok {
                        if me.nbr_sol == MAX_SOL {
                            break;
                        }
                        me.nbr_sol += 1;
                        let n = me.nbr_sol;
                        me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                        me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                        let dc2 = dir_between(&center, &origin_l1);
                        me.qualifier1[n - 1] =
                            qualifier_of_line(qualified1, dc2.dot(&normal));
                        me.the_same1[n - 1] = 0;
                        me.the_same2[n - 1] = 0;
                        let dc1 = dir_between(&center, &origin_l1);
                        let sign = dc1.dot(&normal);
                        let dc1 = signed_normal(sign, &normal);
                        me.pnttg1sol[n - 1] = offset_xy(&center, radius, &dc1);
                        me.pnttg2sol[n - 1] = *point2;
                        me.pntcen[n - 1] = center;
                        me.par1sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                        me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
                        me.par2sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                        me.pararg2[n - 1] = 0.0;
                        me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
                    }
                }
            }
            me.well_done = true;
        }
    }
    me
}

/// `gp_Dir2d(sign * normal.XY())`.
fn signed_normal(sign: f64, normal: &GpDir2d) -> GpDir2d {
    GpDir2d::new(sign * normal.x, sign * normal.y).expect("gp_Dir2d from the signed normal")
}
