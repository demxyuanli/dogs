//! `GccAna_Circ2d2TanOn(QualifiedCirc, gp_Pnt2d, gp_Circ2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_8.cxx:42-315`): circles tangent to a circle and passing
//! through a point, with the centre on a circle.
//!
//! The boundary block (`cxx:80-181`) is skipped entirely when `OnCirc` and `C1`
//! share a centre (`cxx:88`); otherwise it tests the two centres
//! `center1 +/- (distcco -/+ Ron) * dircc`. The general case (`cxx:186-314`)
//! intersects each `GccAna_CircPnt2dBisec` branch of `(C1, Point2)` with
//! `OnCirc`.
//!
//! Translation notes:
//! - `SameCenter` is initialised to `false` (`cxx:88`) and never set, so the
//!   `if (!SameCenter)` guards of `cxx:99`, `:136`, `:139` and `:159` always
//!   hold and the `else` arms are dead; the port keeps the same branches,
//! - the general case skips a candidate whose centre coincides with `C1`'s
//!   (`cxx:258-262`) or that is not on `OnCirc` within `Tol` (`cxx:263-266`),
//!   and stops at four solutions (`cxx:267-270`),
//! - the boundary `continue`-free loop of `cxx:129-176` keeps `WellDone = true`
//!   even when the `i != 1 || !SameCenter` guard skips the store.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::precision::CONFUSION;

use crate::gcc_ana_bisec::GccAnaCircPnt2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use super::{axis_at, dir_between, offset_xy, qualifier_of_circle, GccAnaCirc2d2TanOn, MAX_SOL};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    point2: &GpPnt2d,
    on_circ: &GpCirc2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let tol = tolerance.abs();
    let mut radius = [0.0f64; MAX_SOL];
    let c1 = qualified1.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();

    // `cxx:80-181`: boundary case, skipped when the centres coincide.
    let mut nbsol1 = 1i32;
    let same_center = false;
    if !on_circ.location().is_equal(&center1, CONFUSION) {
        let ron = on_circ.radius();
        let distcco = on_circ.location().distance(&center1);
        let dircc = dir_between(&center1, &on_circ.location());
        let mut pinterm = offset_xy(&center1, distcco - ron, &dircc);
        let mut distcc2 = pinterm.distance(point2);
        let mut distcc1 = pinterm.distance(&center1);
        let mut d1 = (distcc2 - (distcc1 - r1).abs()).abs();
        let mut d2 = (distcc2 - (distcc1 + r1)).abs();
        if d1 > tol || d2 > tol {
            if !same_center {
                pinterm = offset_xy(&center1, distcco + ron, &dircc);
            }
            distcc2 = pinterm.distance(point2);
            distcc1 = pinterm.distance(&center1);
            d1 = (distcc2 - (distcc1 - r1).abs()).abs();
            d2 = (distcc2 - (distcc1 + r1)).abs();
            if d1 > tol || d2 > tol {
                nbsol1 = 0;
            }
        }
        if nbsol1 > 0 {
            if qualified1.is_enclosed() || qualified1.is_outside() {
                nbsol1 = 1;
                radius[0] = (distcc1 - r1).abs();
            } else if qualified1.is_enclosing() {
                nbsol1 = 1;
                radius[0] = r1 + distcc1;
            } else if qualified1.is_unqualified() {
                nbsol1 = 2;
                radius[0] = (distcc1 - r1).abs();
                radius[1] = r1 + distcc1;
            }
            for i in 1..=nbsol1 as usize {
                me.well_done = true;
                if i != 1 || !same_center {
                    if me.nbr_sol == MAX_SOL {
                        break;
                    }
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), radius[i - 1]);
                    let dc1 = dir_between(&pinterm, &center1);
                    me.qualifier1[n - 1] =
                        qualifier_of_circle(qualified1, distcc1, radius[i - 1], r1, tol);
                    me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                    me.pnttg1sol[n - 1] = offset_xy(&pinterm, radius[i - 1], &dc1);
                    me.pnttg2sol[n - 1] = *point2;
                    me.pntcen[n - 1] = me.cirsol[n - 1].location();
                    me.par1sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                    me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                    me.par2sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                    me.pararg2[n - 1] = 0.0;
                    me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
                }
            }
            if me.well_done {
                return me;
            }
        }
    }

    // `cxx:186-314`: general case.
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
                intp.perform_lin_circ(&sol.line().unwrap(), on_circ);
            } else if typ == GccIntIType::Cir {
                intp.perform_circ_circ(on_circ, &sol.circle().unwrap());
            } else if typ == GccIntIType::Hpr {
                intp.perform_circ_conic(
                    on_circ,
                    &IntAna2dConic::from_hypr2d(&sol.hyperbola().unwrap()),
                );
            } else if typ == GccIntIType::Ell {
                intp.perform_circ_conic(
                    on_circ,
                    &IntAna2dConic::from_elips2d(&sol.ellipse().unwrap()),
                );
            }
            if intp.is_done() {
                if !intp.is_empty()
                    && !intp.parallel_elements()
                    && !intp.identical_elements()
                {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = center1.distance(&center);
                        let mut nbsol = 1usize;
                        let mut ok = false;
                        if qualified1.is_enclosed() {
                            if dist1 - c1.radius() <= tol {
                                ok = true;
                                radius[0] = (c1.radius() - dist1).abs();
                            }
                        } else if qualified1.is_outside() {
                            if c1.radius() - dist1 <= tol {
                                ok = true;
                                radius[0] = (c1.radius() - dist1).abs();
                            }
                        } else if qualified1.is_enclosing() {
                            ok = true;
                            radius[0] = c1.radius() + dist1;
                        } else if qualified1.is_unqualified() {
                            ok = true;
                            nbsol = 2;
                            radius[0] = (c1.radius() - dist1).abs();
                            radius[1] = c1.radius() + dist1;
                        }
                        if ok {
                            for k in 1..=nbsol {
                                // `cxx:258-266`: two guards OCCT applies here only.
                                if center1.is_equal(&center, CONFUSION) {
                                    continue;
                                }
                                if on_circ.distance(&center) > tol {
                                    continue;
                                }
                                if me.nbr_sol == 4 {
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
                                if distcc1 <= tol
                                    && (radius[k - 1] - c1.radius()).abs() <= tol
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
                                    me.pararg1[n - 1] =
                                        parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                                }
                                me.the_same2[n - 1] = 0;
                                me.pntcen[n - 1] = center;
                                me.pnttg2sol[n - 1] = *point2;
                                me.pararg2[n - 1] = 0.0;
                                me.par2sol[n - 1] = parameter_circ2d(
                                    &me.cirsol[n - 1],
                                    &me.pnttg2sol[n - 1],
                                );
                                me.parcen3[n - 1] =
                                    parameter_circ2d(on_circ, &me.pntcen[n - 1]);
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
