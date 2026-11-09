//! `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedCirc, gp_Lin2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn.cxx:42-378`): circles tangent to two circles with the
//! centre on a straight line.
//!
//! The four boundary cases (`cxx:99-126`) compare the four combinations
//! `d1 = dist1 + R1`, `d2 = dist2 + R2`, `d3 = dist1 - R1`, `d4 = dist2 - R2`
//! and collect the corresponding radii; the centre is then the intersection of
//! `OnLine` with the line joining the two circle centres (`cxx:127-189`).
//! The general case (`cxx:194-377`) intersects each `GccAna_Circ2dBisec` branch
//! of `(C1, C2)` with `OnLine`.
//!
//! Translation notes:
//! - `Radius`, `Rbid` and `RBid` are `NCollection_Array1<double>(1, 2)` in
//!   OCCT while the boundary block writes up to `Radius(4)` (`cxx:99-126`) and
//!   the `nnsol` loop up to `Radius(4)` (`cxx:293-302`); the port sizes all
//!   three at [`MAX_SOL`],
//! - `Rbid(2)` is never assigned (`cxx:259-262` writes `Rbid(1)` twice) and
//!   `Rbid`/`RBid` are uninitialised otherwise; the port's zero-filled arrays
//!   are the only defined reading of those slots,
//! - the general case compares `Radius(i)` (`cxx:315`, `:320`, `:332`, `:337`)
//!   where `i` is the outer bisec index, not the solution index `k`; the port
//!   keeps the stale index,
//! - `pararg1` of the boundary block is evaluated from `pnttg2sol`
//!   (`cxx:356`); the general case uses `pnttg1sol` (`cxx:354`).

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::GccAnaCirc2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::GccEntQualifiedCirc;
use super::{axis_at, dir_between, qualifier_of_circle, GccAnaCirc2d2TanOn, MAX_SOL};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedCirc,
    on_line: &GpLin2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();
    let tol = tolerance.abs();
    let mut nbsol = 0usize;

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

    let mut dist1 = on_line.distance(&center1);
    let mut dist2 = on_line.distance(&center2);
    let d1 = dist1 + r1;
    let d2 = dist2 + r2;
    let d3 = dist1 - r1;
    let d4 = dist2 - r2;

    // `cxx:99-126`: four boundary cases.
    if (d3 - d4).abs() < tol
        && (qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        nbsol += 1;
        radius[nbsol - 1] = d3.abs();
        me.well_done = true;
    }
    if (d1 - d2).abs() < tol
        && (qualified1.is_enclosing() || qualified1.is_unqualified())
        && (qualified2.is_enclosing() || qualified2.is_unqualified())
    {
        nbsol += 1;
        radius[nbsol - 1] = d1.abs();
        me.well_done = true;
    }
    if (d1 - d4).abs() < tol
        && (qualified1.is_enclosing() || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        nbsol += 1;
        radius[nbsol - 1] = d1.abs();
        me.well_done = true;
    }
    if (d3 - d2).abs() < tol
        && (qualified2.is_enclosing() || qualified2.is_unqualified())
        && (qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
    {
        nbsol += 1;
        radius[nbsol - 1] = d3.abs();
        me.well_done = true;
    }

    // `cxx:127-189`: the centre is `OnLine ^ (center1 -> center2)`.
    let l = GpLin2d::from_pnt_dir(center1, dir_between(&center1, &center2));
    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_lin(on_line, &l);
    if intp.is_done() {
        if !intp.is_empty() {
            for j in 1..=intp.nb_points() {
                let center = *intp.point(j).value();
                for i in 1..=nbsol {
                    if me.nbr_sol == MAX_SOL {
                        break;
                    }
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    let rad = radius[i - 1];
                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), rad);
                    me.well_done = true;

                    let distcc1 = center.distance(&center1);
                    let distcc2 = center.distance(&center2);
                                    me.qualifier1[n - 1] = qualifier_of_circle(
                                        qualified1, distcc1, rad, r1, tol,
                                    );
                                    me.qualifier2[n - 1] = qualifier_of_circle(
                                        qualified2, distcc2, rad, r2, tol,
                                    );

                    let dc1 = dir_between(&center, &center1);
                    let dc2 = dir_between(&center, &center2);
                    me.pnttg1sol[n - 1] =
                        GpPnt2d::new(center.x() + rad * dc1.x, center.y() + rad * dc1.y);
                    me.pnttg2sol[n - 1] =
                        GpPnt2d::new(center.x() + rad * dc2.x, center.y() + rad * dc2.y);
                    me.par1sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                    me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                    me.par2sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                    me.pararg2[n - 1] = parameter_circ2d(&c2, &me.pnttg2sol[n - 1]);
                    me.pntcen[n - 1] = me.cirsol[n - 1].location();
                    me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
                }
            }
        }
    }

    // `cxx:194-377`: general case.
    if !me.well_done {
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
                if typ == GccIntIType::Cir {
                    intp.perform_lin_circ(on_line, &sol.circle().unwrap());
                } else if typ == GccIntIType::Lin {
                    intp.perform_lin_lin(on_line, &sol.line().unwrap());
                } else if typ == GccIntIType::Hpr {
                    intp.perform_lin_conic(
                        on_line,
                        &IntAna2dConic::from_hypr2d(&sol.hyperbola().unwrap()),
                    );
                } else if typ == GccIntIType::Ell {
                    intp.perform_lin_conic(
                        on_line,
                        &IntAna2dConic::from_elips2d(&sol.ellipse().unwrap()),
                    );
                }
                if intp.is_done() {
                    if !intp.is_empty() {
                        for j in 1..=intp.nb_points() {
                            let center = *intp.point(j).value();
                            dist1 = center.distance(&center1);
                            dist2 = center.distance(&center2);
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
                                    let distcc1 = center.distance(&center1);
                                    let distcc2 = center.distance(&center2);
                                    // `cxx:311-341`: the stale `Radius(i)`.
                                    let stale = radius[i - 1];
                                    me.qualifier1[n - 1] =
                                        qualifier_of_circle(qualified1, distcc1, stale, r1, tol);
                                    me.qualifier2[n - 1] =
                                        qualifier_of_circle(qualified2, distcc2, stale, r2, tol);
                                    let rad = radius[k - 1];
                                    if center.distance(&center1) <= tol
                                        && (rad - c1.radius()).abs() <= tol
                                    {
                                        me.the_same1[n - 1] = 1;
                                    } else {
                                        me.the_same1[n - 1] = 0;
                                        let dc1 = dir_between(&center, &center1);
                                        me.pnttg1sol[n - 1] = GpPnt2d::new(
                                            center.x() + rad * dc1.x,
                                            center.y() + rad * dc1.y,
                                        );
                                        me.par1sol[n - 1] = parameter_circ2d(
                                            &me.cirsol[n - 1],
                                            &me.pnttg1sol[n - 1],
                                        );
                                        // `cxx:356`: `pnttg2sol`, not `pnttg1sol`.
                                        me.pararg1[n - 1] =
                                            parameter_circ2d(&c1, &me.pnttg2sol[n - 1]);
                                    }
                                    if center.distance(&center2) <= tol
                                        && (rad - c2.radius()).abs() <= tol
                                    {
                                        me.the_same2[n - 1] = 1;
                                    } else {
                                        me.the_same2[n - 1] = 0;
                                        let dc2 = dir_between(&center, &center2);
                                        me.pnttg2sol[n - 1] = GpPnt2d::new(
                                            center.x() + rad * dc2.x,
                                            center.y() + rad * dc2.y,
                                        );
                                        me.par2sol[n - 1] = parameter_circ2d(
                                            &me.cirsol[n - 1],
                                            &me.pnttg2sol[n - 1],
                                        );
                                        me.pararg2[n - 1] =
                                            parameter_circ2d(&c2, &me.pnttg2sol[n - 1]);
                                    }
                                    me.pntcen[n - 1] = center;
                                    me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
                                }
                            }
                        }
                    }
                    me.well_done = true;
                }
            }
        }
    }
    me
}
