//! `GccAna_Circ2d2TanOn(QualifiedCirc, gp_Pnt2d, gp_Lin2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_3.cxx:48-276`): circles tangent to a circle and passing
//! through a point, with the centre on a straight line.
//!
//! The boundary block (`cxx:88-146`) projects `Point2` on `OnLine` and tests the
//! three circle qualifiers against `dist`. The general case (`cxx:152-276`)
//! intersects each `GccAna_CircPnt2dBisec` branch of `(C1, Point2)` with
//! `OnLine`.
//!
//! Translation notes:
//! - `cxx:90-93` recomputes `pinterm` with the opposite sign only when the first
//!   projection is farther than `Tol`, so both arms are meaningful here,
//! - the `IsUnqualified` arm (`cxx:214-222`) tests `std::abs(popradius - dist1)`
//!   as a truthiness condition on a `double`, i.e. "not exactly zero"; the port
//!   keeps `!= 0.0`,
//! - `cxx:260` writes `pararg1(i)` with the outer bisec index `i`, not the
//!   solution index `NbrSol`; the port keeps the stale index,
//! - the `nbsol = 2` of the commented-out `IsUnqualified` block (`cxx:207-213`)
//!   never executes in OCCT, so every unqualified solution stores `Radius(1)`
//!   only.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::GccAnaCircPnt2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use super::{
    axis_at, dir_between, left_normal, offset_xy, qualifier_of_circle, GccAnaCirc2d2TanOn, MAX_SOL,
};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    point2: &GpPnt2d,
    on_line: &GpLin2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();
    let tol = tolerance.abs();

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let mut radius = [0.0f64; MAX_SOL];
    let c1 = qualified1.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();

    // `cxx:88-146`: boundary case.
    let dp2l = on_line.distance(point2);
    let on_norm = left_normal(on_line.direction());
    let mut pinterm = offset_xy(point2, dp2l, &on_norm);
    if on_line.distance(&pinterm) > tol {
        pinterm = offset_xy(point2, -dp2l, &on_norm);
    }
    let dist = pinterm.distance(&center1);
    if qualified1.is_enclosed() && (r1 - dist - dp2l).abs() <= tol {
        me.well_done = true;
    } else if qualified1.is_enclosing() && (r1 + dist - dp2l).abs() <= tol {
        me.well_done = true;
    } else if qualified1.is_outside() && (dist - dp2l).abs() <= tol {
        me.well_done = true;
    } else if qualified1.is_unqualified()
        && ((dist - dp2l).abs() <= tol
            || (r1 - dist - dp2l).abs() <= tol
            || (r1 + dist - dp2l).abs() <= tol)
    {
        me.well_done = true;
    }
    if me.well_done {
        me.nbr_sol = 1;
        let n = me.nbr_sol;
        me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), dp2l);

        let dc1 = dir_between(&pinterm, &center1);
        let distcc1 = pinterm.distance(&center1);
        me.qualifier1[n - 1] = qualifier_of_circle(qualified1, distcc1, dp2l, r1, tol);
        me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
        me.pnttg1sol[n - 1] = offset_xy(&pinterm, dp2l, &dc1);
        me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
        me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
        me.pnttg2sol[n - 1] = *point2;
        me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
        me.pntcen[n - 1] = me.cirsol[n - 1].location();
        me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
        return me;
    }

    // `cxx:152-276`: general case.
    let bis = GccAnaCircPnt2dBisec::new(&c1, point2);
    if bis.is_done() {
        let nbsolution = bis.nb_solutions();
        for i in 1..=nbsolution {
            let sol = bis
                .this_solution(i)
                .expect("GccAna_CircPnt2dBisec::ThisSolution");
            let typ = sol.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
            if typ == GccIntIType::Lin {
                intp.perform_lin_lin(on_line, &sol.line().unwrap());
            } else if typ == GccIntIType::Cir {
                intp.perform_lin_circ(on_line, &sol.circle().unwrap());
            } else if typ == GccIntIType::Ell {
                intp.perform_lin_conic(
                    on_line,
                    &IntAna2dConic::from_elips2d(&sol.ellipse().unwrap()),
                );
            } else if typ == GccIntIType::Hpr {
                intp.perform_lin_conic(
                    on_line,
                    &IntAna2dConic::from_hypr2d(&sol.hyperbola().unwrap()),
                );
            }
            if intp.is_done() {
                if !intp.is_empty() {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = center1.distance(&center);
                        let nbsol = 1usize;
                        let mut ok = false;
                        if qualified1.is_enclosed() {
                            if dist1 - c1.radius() <= tolerance {
                                ok = true;
                                radius[0] = (c1.radius() - dist1).abs();
                            }
                        } else if qualified1.is_outside() {
                            if c1.radius() - dist1 <= tolerance {
                                ok = true;
                                radius[0] = (c1.radius() - dist1).abs();
                            }
                        } else if qualified1.is_enclosing() {
                            ok = true;
                            radius[0] = c1.radius() + dist1;
                        } else if qualified1.is_unqualified() {
                            let popradius = center.distance(point2);
                            if (popradius - dist1) != 0.0 {
                                ok = true;
                                radius[0] = popradius;
                            }
                        }
                        if ok {
                            for k in 1..=nbsol {
                                if me.nbr_sol == MAX_SOL {
                                    break;
                                }
                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
                                me.cirsol[n - 1] =
                                    GpCirc2d::new(axis_at(center), radius[k - 1]);
                                let distcc1 = center.distance(&center1);
                                me.qualifier1[n - 1] = qualifier_of_circle(
                                    qualified1,
                                    distcc1,
                                    radius[k - 1],
                                    r1,
                                    tol,
                                );
                                me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                                if center.distance(&center1) <= tolerance
                                    && (radius[k - 1] - c1.radius()).abs() <= tolerance
                                {
                                    me.the_same1[n - 1] = 1;
                                } else {
                                    me.the_same1[n - 1] = 0;
                                    let dc1 = dir_between(&center, &center1);
                                    me.pnttg1sol[n - 1] =
                                        offset_xy(&center, radius[k - 1], &dc1);
                                    me.par1sol[n - 1] = parameter_circ2d(
                                        &me.cirsol[n - 1],
                                        &me.pnttg1sol[n - 1],
                                    );
                                    // `cxx:260`: the stale `i`, not `NbrSol`.
                                    me.pararg1[i as usize - 1] =
                                        parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                                }
                                me.the_same2[n - 1] = 0;
                                me.pnttg2sol[n - 1] = *point2;
                                me.par2sol[n - 1] = parameter_circ2d(
                                    &me.cirsol[n - 1],
                                    &me.pnttg2sol[n - 1],
                                );
                                me.pararg2[n - 1] = 0.0;
                                me.pntcen[n - 1] = center;
                                me.parcen3[n - 1] =
                                    parameter_lin2d(on_line, &me.pntcen[n - 1]);
                            }
                        }
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}
