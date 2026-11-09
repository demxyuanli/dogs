//! `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedCirc, gp_Circ2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_6.cxx:45-401`): circles tangent to two circles with the
//! centre on a circle.
//!
//! The boundary block (`cxx:90-211`) walks the two candidate centres
//! `center1 +/- (distcco -/+ Ron) * dircc`, keeps the one whose radius matches
//! `C2` and collects the qualifying radii. The general case (`cxx:216-400`)
//! intersects each `GccAna_Circ2dBisec` branch of `(C1, C2)` with `OnCirc`.
//!
//! Translation notes:
//! - `Rbid(2)` is never assigned (`cxx:278-282` writes `Rbid(1)` twice) and
//!   `Rbid`/`RBid` are otherwise uninitialised while the `isol`/`jsol` loops
//!   read `Rbid(isol)` for `isol <= 2`; the port's zero-filled arrays are the
//!   only defined reading of those slots,
//! - the `nnsol` loop writes up to `Radius(4)` (`cxx:320`) into an OCCT
//!   `(1, 2)` array; the port sizes it at [`MAX_SOL`],
//! - unlike `GccAna_Circ2d2TanOn.cxx:315` this file qualifies with `Radius(k)`
//!   (`cxx:337`, `:342`, `:353`, `:358`), the solution index; the port keeps
//!   the difference,
//! - `cxx:227` and `:237` call `Perform` with the arguments in the other order
//!   (`Sol->Line()` first); both map to [`IntAna2dAnaIntersection::perform_lin_circ`].

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::GccAnaCirc2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::GccEntQualifiedCirc;
use super::{axis_at, dir_between, offset_xy, qualifier_of_circle, GccAnaCirc2d2TanOn, MAX_SOL};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedCirc,
    on_circ: &GpCirc2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();
    let tol = tolerance.abs();

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
        || !(qualified2.is_enclosed()
            || qualified2.is_enclosing()
            || qualified2.is_outside()
            || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let c1 = qualified1.qualified();
    let c2 = qualified2.qualified();
    let mut r1 = c1.radius();
    let mut r2 = c2.radius();
    let center1 = c1.location();
    let center2 = c2.location();
    let mut radius = [0.0f64; MAX_SOL];
    let mut rradius = [0.0f64; MAX_SOL];

    // `cxx:90-211`: boundary case.
    let mut nbsol1 = 1i32;
    let mut nbsol2 = 0i32;
    let ron = on_circ.radius();
    let distcco = on_circ.location().distance(&center1);
    let dircc = dir_between(&center1, &on_circ.location());
    let mut pinterm = offset_xy(&center1, distcco - ron, &dircc);
    let mut distcc2 = pinterm.distance(&center2);
    let mut distcc1 = pinterm.distance(&center1);
    let mut d1 = (distcc2 - r2 - (distcc1 - r1).abs()).abs();
    let mut d2 = (distcc2 + r2 - (distcc1 - r1).abs()).abs();
    let mut d3 = (distcc2 - r2 - (distcc1 + r1)).abs();
    let mut d4 = (distcc2 + r2 - (distcc1 + r1)).abs();
    if d1 > tol || d2 > tol || d3 > tol || d4 > tol {
        pinterm = offset_xy(&center1, distcco + ron, &dircc);
        distcc2 = pinterm.distance(&center2);
        distcc1 = pinterm.distance(&center1);
        d1 = (distcc2 - r2 - (distcc1 - r1).abs()).abs();
        d2 = (distcc2 + r2 - (distcc1 - r1).abs()).abs();
        d3 = (distcc2 - r2 - (distcc1 + r1)).abs();
        d4 = (distcc2 + r2 - (distcc1 + r1)).abs();
        if d1 > tol || d2 > tol || d3 > tol || d4 > tol {
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
        if qualified2.is_enclosed() || qualified2.is_outside() {
            nbsol2 = 1;
            rradius[0] = (distcc2 - r2).abs();
        } else if qualified2.is_enclosing() {
            nbsol2 = 1;
            rradius[0] = r2 + distcc2;
        } else if qualified2.is_unqualified() {
            nbsol2 = 2;
            rradius[0] = (distcc2 - r2).abs();
            rradius[1] = r2 + distcc2;
        }
        for i in 1..=nbsol1 as usize {
            for j in 1..=nbsol2 as usize {
                if (radius[i - 1] - rradius[j - 1]).abs() <= tol {
                    me.well_done = true;
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), radius[i - 1]);

                    let dc1 = dir_between(&pinterm, &center1);
                    let dc2 = dir_between(&pinterm, &center2);
                    distcc1 = pinterm.distance(&center1);
                    distcc2 = pinterm.distance(&center2);
                    me.qualifier1[n - 1] =
                        qualifier_of_circle(qualified1, distcc1, radius[i - 1], r1, tol);
                    me.qualifier2[n - 1] =
                        qualifier_of_circle(qualified2, distcc2, radius[i - 1], r2, tol);
                    me.pnttg1sol[n - 1] = offset_xy(&pinterm, radius[i - 1], &dc1);
                    me.pnttg2sol[n - 1] = offset_xy(&pinterm, radius[i - 1], &dc2);
                    me.pntcen[n - 1] = me.cirsol[n - 1].location();
                    me.par1sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                    me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                    me.par2sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                    me.pararg2[n - 1] = parameter_circ2d(&c2, &me.pnttg2sol[n - 1]);
                    me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
                }
            }
        }
        if me.well_done {
            return me;
        }
    }

    // `cxx:216-400`: general case.
    let bis = GccAnaCirc2dBisec::new(&c1, &c2);
    if bis.is_done() {
        let mut rbid = [0.0f64; MAX_SOL];
        let mut rbid2 = [0.0f64; MAX_SOL];
        let nbsolution = bis
            .nb_solutions()
            .expect("GccAna_Circ2dBisec::NbSolutions");
        for i in 1..=nbsolution as usize {
            let sol = bis
                .this_solution(i as i32)
                .expect("GccAna_Circ2dBisec::ThisSolution");
            let typ = sol.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
            if typ == GccIntIType::Cir {
                intp.perform_circ_circ(on_circ, &sol.circle().unwrap());
            } else if typ == GccIntIType::Lin {
                intp.perform_lin_circ(&sol.line().unwrap(), on_circ);
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
                        let dist1 = center.distance(&center1);
                        let dist2 = center.distance(&center2);
                        let mut nbsol = 0usize;
                        let mut nsol = 0usize;
                        let mut nnsol = 0usize;
                        r1 = c1.radius();
                        r2 = c2.radius();
                        if qualified1.is_enclosed() {
                            if dist1 - r1 < tol {
                                nbsol = 1;
                                rbid[0] = (r1 - dist1).abs();
                            }
                        } else if qualified1.is_outside() {
                            if r1 - dist1 < tol {
                                nbsol = 1;
                                rbid[0] = (dist1 - r1).abs();
                            }
                        } else if qualified1.is_enclosing() {
                            nbsol = 1;
                            rbid[0] = dist1 + r1;
                        } else if qualified1.is_unqualified() {
                            nbsol = 2;
                            rbid[0] = dist1 + r1;
                            rbid[0] = (dist1 - r1).abs();
                        }
                        if qualified2.is_enclosed() && nbsol != 0 {
                            if dist2 - r2 < tol {
                                nsol = 1;
                                rbid2[0] = (r2 - dist2).abs();
                            }
                        } else if qualified2.is_outside() && nbsol != 0 {
                            if r2 - dist2 < tol {
                                nsol = 1;
                                rbid2[0] = (r2 - dist2).abs();
                            }
                        } else if qualified2.is_enclosing() && nbsol != 0 {
                            nsol = 1;
                            rbid2[0] = dist2 + r2;
                        } else if qualified2.is_unqualified() && nbsol != 0 {
                            nsol = 2;
                            rbid2[0] = dist2 + r2;
                            rbid2[1] = (r2 - dist2).abs();
                        }
                        for isol in 1..=nbsol {
                            for jsol in 1..=nsol {
                                if (rbid[isol - 1] - rbid2[jsol - 1]).abs() <= tol {
                                    nnsol += 1;
                                    radius[nnsol - 1] =
                                        (rbid2[jsol - 1] + rbid[isol - 1]) / 2.0;
                                }
                            }
                        }
                        if nnsol > 0 {
                            for k in 1..=nnsol {
                                if me.nbr_sol == MAX_SOL {
                                    break;
                                }
                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
                                me.cirsol[n - 1] =
                                    GpCirc2d::new(axis_at(center), radius[k - 1]);
                                distcc1 = center.distance(&center1);
                                distcc2 = center.distance(&center2);
                                me.qualifier1[n - 1] = qualifier_of_circle(
                                    qualified1,
                                    distcc1,
                                    radius[k - 1],
                                    r1,
                                    tol,
                                );
                                me.qualifier2[n - 1] = qualifier_of_circle(
                                    qualified2,
                                    distcc2,
                                    radius[k - 1],
                                    r2,
                                    tol,
                                );
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
                                    me.pararg1[n - 1] =
                                        parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                                }
                                if center.distance(&center2) <= tolerance
                                    && (radius[k - 1] - c2.radius()).abs() <= tolerance
                                {
                                    me.the_same2[n - 1] = 1;
                                } else {
                                    me.the_same2[n - 1] = 0;
                                    let dc2 = dir_between(&center, &center2);
                                    me.pnttg2sol[n - 1] =
                                        offset_xy(&center, radius[k - 1], &dc2);
                                    me.par2sol[n - 1] = parameter_circ2d(
                                        &me.cirsol[n - 1],
                                        &me.pnttg2sol[n - 1],
                                    );
                                    me.pararg2[n - 1] =
                                        parameter_circ2d(&c2, &me.pnttg2sol[n - 1]);
                                }
                                me.pntcen[n - 1] = center;
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
