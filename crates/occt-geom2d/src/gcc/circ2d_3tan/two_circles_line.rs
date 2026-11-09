//! `GccAna_Circ2d3Tan(GccEnt_QualifiedCirc, GccEnt_QualifiedCirc,
//! GccEnt_QualifiedLin, Tolerance)` (`GccAna_Circ2d3Tan_1.cxx:38-415`): circles
//! tangent to two circles and a straight line.
//!
//! Centers come from intersecting the `GccAna_CircLin2dBisec` of `(C1, L3)`
//! with the one of `(C2, L3)` (`cxx:96-133`); the radius is then selected by the
//! qualifier of all three arguments (`cxx:159-290`).
//!
//! Translation notes:
//! - `cxx:271` and `cxx:281` test `Qualified2.IsOutside()` /
//!   `Qualified2.IsUnqualified()` inside the *third* qualifier's chain, i.e. the
//!   `L3` qualifier is never consulted there. The port keeps that slip,
//! - `cxx:489-497`: `distcc2` is `Center.Distance(center1)`, not `center2`; kept,
//! - `cxx:387-390` sets `par1sol` (not `par2sol`) inside the `TheSame2` block
//!   when the solution center equals `pnttg1sol`; kept,
//! - the `Center.X() > 1e10` / `dist3 > 1e10` early `break`s of `cxx:144-159`
//!   are kept as written,
//! - comparisons of this constructor use `< Tol` where the circ/circ/circ one
//!   uses `<= Tol` (`cxx:299`, `:303`, ...); the port keeps each operator.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::precision::CONFUSION;

use crate::gcc_ana_bisec::GccAnaCircLin2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use super::{axis_at, GccAnaCirc2d3Tan};

/// `nbsol2 = 1` at `cxx:203/218/231` is a dead store in OCCT (only the
/// `IsUnqualified` arm reads `nbsol2`); the port keeps it for branch parity.
#[allow(unused_assignments)]
pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedCirc,
    qualified3: &GccEntQualifiedLin,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !((qualified1.is_enclosed() || qualified1.is_enclosing() || qualified1.is_outside()
        || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_enclosing() || qualified2.is_outside()
            || qualified2.is_unqualified())
        && (qualified3.is_enclosed() || qualified3.is_outside() || qualified3.is_unqualified()))
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let c1 = qualified1.qualified();
    let c2 = qualified2.qualified();
    let l3 = qualified3.qualified();
    let r1 = c1.radius();
    let r2 = c2.radius();
    let center1 = c1.location();

    let origin3 = l3.location();
    let dir3 = *l3.direction();
    let norm_l3 = GpDir2d::new(-dir3.y, dir3.x).expect("gp_Dir2d(-dir.Y(), dir.X())");

    let mut radius = [0.0f64; 3];

    let bis1 = GccAnaCircLin2dBisec::new(&c1, &l3);
    let bis2 = GccAnaCircLin2dBisec::new(&c2, &l3);
    if bis1.is_done() && bis2.is_done() {
        let nbsolution1 = bis1.nb_solutions();
        let nbsolution2 = bis2.nb_solutions();
        for i in 1..=nbsolution1 {
            let sol1 = bis1.this_solution(i).expect("GccAna_CircLin2dBisec::ThisSolution");
            let typ1 = sol1.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
            for k in 1..=nbsolution2 {
                let sol2 = bis2.this_solution(k).expect("GccAna_CircLin2dBisec::ThisSolution");
                let typ2 = sol2.arc_type();
                if typ1 == GccIntIType::Lin {
                    if typ2 == GccIntIType::Lin {
                        intp.perform_lin_lin(&sol1.line().unwrap(), &sol2.line().unwrap());
                    } else if typ2 == GccIntIType::Par {
                        intp.perform_lin_conic(
                            &sol1.line().unwrap(),
                            &IntAna2dConic::from_parab2d(&sol2.parabola().unwrap()),
                        );
                    }
                } else if typ1 == GccIntIType::Par {
                    if typ2 == GccIntIType::Lin {
                        intp.perform_lin_conic(
                            &sol2.line().unwrap(),
                            &IntAna2dConic::from_parab2d(&sol1.parabola().unwrap()),
                        );
                    } else if typ2 == GccIntIType::Par {
                        intp.perform_parab_conic(
                            &sol1.parabola().unwrap(),
                            &IntAna2dConic::from_parab2d(&sol2.parabola().unwrap()),
                        );
                    }
                }
                if intp.is_done() {
                    if !intp.is_empty() {
                        for j in 1..=intp.nb_points() {
                            let mut rradius = 0.0;
                            let center = *intp.point(j).value();

                            // `cxx:143-147`.
                            if center.x() > 1e10 || center.y() > 1e10 {
                                break;
                            }

                            let dist1 = center.distance(&c1.location());
                            let dist2 = center.distance(&c2.location());
                            let dist3 = l3.distance(&center);

                            if dist3 > 1e10 {
                                break;
                            }

                            let mut nbsol1 = 0i32;
                            let mut nbsol2 = 0i32;
                            let mut nbsol3 = 0i32;
                            let mut ok = false;
                            if qualified1.is_enclosed() {
                                if dist1 - r1 < tolerance {
                                    radius[0] = (r1 - dist1).abs();
                                    nbsol1 = 1;
                                    ok = true;
                                }
                            } else if qualified1.is_outside() {
                                if r1 - dist1 < tolerance {
                                    radius[0] = (r1 - dist1).abs();
                                    nbsol1 = 1;
                                    ok = true;
                                }
                            } else if qualified1.is_enclosing() {
                                ok = true;
                                nbsol1 = 1;
                                radius[0] = (r1 - dist1).abs();
                            } else if qualified1.is_unqualified() {
                                ok = true;
                                nbsol1 = 2;
                                radius[0] = (r1 - dist1).abs();
                                radius[1] = r1 + dist1;
                            }
                            if qualified2.is_enclosed() && ok {
                                if dist2 - r2 < tolerance {
                                    for ii in 1..=nbsol1 {
                                        if (radius[ii as usize - 1] - (r2 - dist2).abs()) < tol {
                                            radius[0] = (r2 - dist2).abs();
                                            ok = true;
                                            nbsol2 = 1;
                                        }
                                    }
                                }
                            } else if qualified2.is_outside() && ok {
                                if r2 - dist2 < tolerance {
                                    for ii in 1..=nbsol1 {
                                        if (radius[ii as usize - 1] - (r2 - dist2).abs()) < tol {
                                            radius[0] = (r2 - dist2).abs();
                                            ok = true;
                                            nbsol2 = 1;
                                        }
                                    }
                                }
                            } else if qualified2.is_enclosing() && ok {
                                for ii in 1..=nbsol1 {
                                    if (radius[ii as usize - 1] - r2 - dist2).abs() < tol {
                                        radius[0] = r2 + dist2;
                                        ok = true;
                                        nbsol2 = 1;
                                    }
                                }
                            } else if qualified2.is_unqualified() && ok {
                                for ii in 1..=nbsol1 {
                                    if (radius[ii as usize - 1] - (r2 - dist2).abs()) < tol {
                                        rradius = (r2 - dist2).abs();
                                        ok = true;
                                        nbsol2 += 1;
                                    } else if (radius[ii as usize - 1] - r2 - dist2).abs() < tol {
                                        rradius = r2 + dist2;
                                        ok = true;
                                        nbsol2 += 1;
                                    }
                                }
                                if nbsol2 == 1 {
                                    radius[0] = rradius;
                                } else if nbsol2 == 2 {
                                    radius[0] = (r2 - dist2).abs();
                                    radius[1] = r2 + dist2;
                                }
                            }
                            if qualified3.is_enclosed() && ok {
                                if (((origin3.x() - center.x()) * (-dir3.y))
                                    + ((origin3.y() - center.y()) * (dir3.x)))
                                    <= 0.0
                                {
                                    ok = true;
                                    nbsol3 = 1;
                                }
                            } else if qualified2.is_outside() && ok {
                                // `cxx:271`: `Qualified2`, not `Qualified3`.
                                if (((origin3.x() - center.x()) * (-dir3.y))
                                    + ((origin3.y() - center.y()) * (dir3.x)))
                                    >= 0.0
                                {
                                    ok = true;
                                    nbsol3 = 1;
                                }
                            } else if qualified2.is_unqualified() && ok {
                                // `cxx:281`: `Qualified2`, not `Qualified3`.
                                ok = true;
                                nbsol3 = 1;
                            }
                            if ok {
                                for ind3 in 1..=nbsol3 {
                                    if me.nbr_sol >= super::MAX_SOL {
                                        break;
                                    }
                                    me.nbr_sol += 1;
                                    let n = me.nbr_sol;
                                    let rad = radius[ind3 as usize - 1];
                                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), rad);

                                    let distcc1 = center.distance(&center1);
                                    me.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                        qualified1.qualifier()
                                    } else if (distcc1 + rad - r1).abs() < tol {
                                        GccEntPosition::Enclosed
                                    } else if (distcc1 - r1 - rad).abs() < tol {
                                        GccEntPosition::Outside
                                    } else {
                                        GccEntPosition::Enclosing
                                    };

                                    // `cxx:312`: `center1`, not `center2`.
                                    let distcc2 = center.distance(&center1);
                                    me.qualifier2[n - 1] = if !qualified2.is_unqualified() {
                                        qualified2.qualifier()
                                    } else if (distcc2 + rad - r2).abs() < tol {
                                        GccEntPosition::Enclosed
                                    } else if (distcc2 - r2 - rad).abs() < tol {
                                        GccEntPosition::Outside
                                    } else {
                                        GccEntPosition::Enclosing
                                    };

                                    let dc3 = GpDir2d::new(origin3.x() - center.x(), origin3.y() - center.y())
                                        .expect("gp_Dir2d(origin3.XY() - Center.XY())");
                                    me.qualifier3[n - 1] = if !qualified3.is_unqualified() {
                                        qualified3.qualifier()
                                    } else if dc3.dot(&norm_l3) > 0.0 {
                                        GccEntPosition::Outside
                                    } else {
                                        GccEntPosition::Enclosed
                                    };

                                    if center.distance(&c1.location()) <= tolerance
                                        && (rad - r1).abs() <= tolerance
                                    {
                                        me.the_same1[n - 1] = 1;
                                    } else {
                                        me.the_same1[n - 1] = 0;
                                        let dc = GpDir2d::new(
                                            c1.location().x() - center.x(),
                                            c1.location().y() - center.y(),
                                        )
                                        .expect("gp_Dir2d(C1.Location().XY() - Center.XY())");
                                        me.pnttg1sol[n - 1] =
                                            GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                        if me.cirsol[n - 1].location().is_equal(
                                            &me.pnttg1sol[n - 1],
                                            CONFUSION,
                                        ) {
                                            me.par1sol[n - 1] = 1.0;
                                        } else {
                                            me.par1sol[n - 1] =
                                                parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                                        }
                                        if c1.location().is_equal(&me.pnttg1sol[n - 1], CONFUSION) {
                                            me.pararg1[n - 1] = 1.0;
                                        } else {
                                            me.pararg1[n - 1] =
                                                parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                                        }
                                    }

                                    if center.distance(&c2.location()) <= tolerance
                                        && (rad - r2).abs() <= tolerance
                                    {
                                        me.the_same2[n - 1] = 1;
                                    } else {
                                        me.the_same2[n - 1] = 0;
                                        let dc = GpDir2d::new(
                                            c2.location().x() - center.x(),
                                            c2.location().y() - center.y(),
                                        )
                                        .expect("gp_Dir2d(C2.Location().XY() - Center.XY())");
                                        me.pnttg2sol[n - 1] =
                                            GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                        // `cxx:383-390`: OCCT compares against
                                        // `pnttg1sol` and writes `par1sol`.
                                        if me.cirsol[n - 1].location().is_equal(
                                            &me.pnttg1sol[n - 1],
                                            CONFUSION,
                                        ) {
                                            me.par1sol[n - 1] = 1.0;
                                        } else {
                                            me.par2sol[n - 1] =
                                                parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                                        }
                                        if c2.location().is_equal(&me.pnttg2sol[n - 1], CONFUSION) {
                                            me.pararg2[n - 1] = 1.0;
                                        } else {
                                            me.pararg2[n - 1] =
                                                parameter_circ2d(&c2, &me.pnttg2sol[n - 1]);
                                        }
                                    }

                                    me.the_same3[n - 1] = 0;
                                    let dc = GpDir2d::new(
                                        l3.location().x() - center.x(),
                                        l3.location().y() - center.y(),
                                    )
                                    .expect("gp_Dir2d(L3.Location().XY() - Center.XY())");
                                    let sign = dc.dot(&norm_l3);
                                    let dc = GpDir2d::new(sign * norm_l3.x, sign * norm_l3.y)
                                        .expect("gp_Dir2d(sign * normL3.XY())");
                                    me.pnttg3sol[n - 1] =
                                        GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                    me.par3sol[n - 1] =
                                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg3sol[n - 1]);
                                    me.pararg3[n - 1] =
                                        parameter_lin2d(&l3, &me.pnttg3sol[n - 1]);
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
