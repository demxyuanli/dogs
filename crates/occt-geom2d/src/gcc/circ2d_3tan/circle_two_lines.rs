//! `GccAna_Circ2d3Tan(GccEnt_QualifiedCirc, GccEnt_QualifiedLin,
//! GccEnt_QualifiedLin, Tolerance)` (`GccAna_Circ2d3Tan_2.cxx:38-323`): circles
//! tangent to one circle and two straight lines.
//!
//! Centers come from intersecting the `GccAna_CircLin2dBisec` of `(C1, L2)`
//! with the `GccAna_Lin2dBisec` of `(L2, L3)` (`cxx:95-119`); the radius is
//! fixed by the circle qualifier and then validated against both lines
//! (`cxx:160-237`).
//!
//! Translation notes:
//! - `cxx:297` is the one place in this constructor that writes `pnttg1sol`
//!   without a `Precision::Confusion()` guard, while `_1.cxx` has one; the port
//!   keeps each file's own variant,
//! - comparisons here are `< Tol` (`cxx:247-259`), not `<= Tol`,
//! - `NCollection_Array1<...>(1, 8)` becomes [`super::MAX_SOL`], and the
//!   `NbrSol` cap of this constructor is 8 (the loops cannot exceed it).

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::{GccAnaCircLin2dBisec, GccAnaLin2dBisec};
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedLin,
    qualified3: &GccEntQualifiedLin,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !((qualified1.is_enclosed() || qualified1.is_enclosing() || qualified1.is_outside()
        || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
        && (qualified3.is_enclosed() || qualified3.is_outside() || qualified3.is_unqualified()))
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let c1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let l3 = qualified3.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let origin2 = l2.location();
    let dir2 = *l2.direction();
    let norm_l2 = GpDir2d::new(-dir2.y, dir2.x).expect("gp_Dir2d(-dir.Y(), dir.X())");
    let origin3 = l3.location();
    let dir3 = *l3.direction();
    let norm_l3 = GpDir2d::new(-dir3.y, dir3.x).expect("gp_Dir2d(-dir.Y(), dir.X())");

    let mut radius = [0.0f64; 3];

    let bis1 = GccAnaCircLin2dBisec::new(&c1, &l2);
    let bis2 = GccAnaLin2dBisec::new(&l2, &l3);
    if bis1.is_done() && bis2.is_done() {
        let nbsolution1 = bis1.nb_solutions();
        let nbsolution2 = bis2.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions");
        for i in 1..=nbsolution1 {
            let sol1 = bis1.this_solution(i).expect("GccAna_CircLin2dBisec::ThisSolution");
            let typ1 = sol1.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
            for k in 1..=nbsolution2 {
                let sol2 = bis2.this_solution(k as usize).expect("GccAna_Lin2dBisec::ThisSolution");
                if typ1 == GccIntIType::Lin {
                    intp.perform_lin_lin(&sol1.line().unwrap(), &sol2);
                } else if typ1 == GccIntIType::Par {
                    intp.perform_lin_conic(
                        &sol2,
                        &IntAna2dConic::from_parab2d(&sol1.parabola().unwrap()),
                    );
                }
                if intp.is_done() && !intp.is_empty() && !intp.parallel_elements()
                    && !intp.identical_elements()
                {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = center.distance(&center1);
                        let dist2 = l2.distance(&center);
                        let dist3 = l3.distance(&center);
                        let mut nbsol1 = 0i32;
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
                            if (((origin2.x() - center.x()) * (-dir2.y))
                                + ((origin2.y() - center.y()) * (dir2.x)))
                                <= 0.0
                            {
                                for ii in 1..=nbsol1 {
                                    if (dist2 - radius[ii as usize - 1]).abs() < tol {
                                        ok = true;
                                        radius[0] = radius[ii as usize - 1];
                                    }
                                }
                            }
                        } else if qualified2.is_outside() && ok {
                            if (((origin2.x() - center.x()) * (-dir2.y))
                                + ((origin2.y() - center.y()) * (dir2.x)))
                                >= 0.0
                            {
                                for ii in 1..=nbsol1 {
                                    if (dist2 - radius[ii as usize - 1]).abs() < tol {
                                        ok = true;
                                        radius[0] = radius[ii as usize - 1];
                                    }
                                }
                            }
                        } else if qualified2.is_unqualified() && ok {
                            for ii in 1..=nbsol1 {
                                if (dist2 - radius[ii as usize - 1]).abs() < tol {
                                    ok = true;
                                    radius[0] = radius[ii as usize - 1];
                                }
                            }
                        }
                        if qualified3.is_enclosed() && ok {
                            if (((origin3.x() - center.x()) * (-dir3.y))
                                + ((origin3.y() - center.y()) * (dir3.x)))
                                <= 0.0
                                && (dist3 - radius[0]).abs() < tol
                            {
                                ok = true;
                                nbsol3 = 1;
                            }
                        } else if qualified3.is_outside() && ok {
                            if (((origin3.x() - center.x()) * (-dir3.y))
                                + ((origin3.y() - center.y()) * (dir3.x)))
                                >= 0.0
                                && (dist3 - radius[0]).abs() < tol
                            {
                                ok = true;
                                nbsol3 = 1;
                            }
                        } else if qualified3.is_unqualified() && ok {
                            if (dist3 - radius[0]).abs() < tol {
                                ok = true;
                                nbsol3 = 1;
                            }
                        }
                        if ok {
                            for m in 1..=nbsol3 {
                                if me.nbr_sol >= super::MAX_SOL {
                                    break;
                                }
                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
                                let rad = radius[m as usize - 1];
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

                                let dc2 = GpDir2d::new(
                                    origin2.x() - center.x(),
                                    origin2.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin2.XY() - Center.XY())");
                                me.qualifier2[n - 1] = if !qualified2.is_unqualified() {
                                    qualified2.qualifier()
                                } else if dc2.dot(&norm_l2) > 0.0 {
                                    GccEntPosition::Outside
                                } else {
                                    GccEntPosition::Enclosed
                                };

                                let dc3 = GpDir2d::new(
                                    origin3.x() - center.x(),
                                    origin3.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin3.XY() - Center.XY())");
                                me.qualifier3[n - 1] = if !qualified3.is_unqualified() {
                                    qualified3.qualifier()
                                } else if dc3.dot(&norm_l3) > 0.0 {
                                    GccEntPosition::Outside
                                } else {
                                    GccEntPosition::Enclosed
                                };

                                if center.distance(&center1) <= tolerance
                                    && (rad - r1).abs() <= tolerance
                                {
                                    me.the_same1[n - 1] = 1;
                                } else {
                                    me.the_same1[n - 1] = 0;
                                    let dc = GpDir2d::new(
                                        center1.x() - center.x(),
                                        center1.y() - center.y(),
                                    )
                                    .expect("gp_Dir2d(center1.XY() - Center.XY())");
                                    me.pnttg1sol[n - 1] =
                                        GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                    me.par1sol[n - 1] =
                                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                                    me.pararg1[n - 1] =
                                        parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                                }

                                me.the_same2[n - 1] = 0;
                                me.the_same3[n - 1] = 0;
                                let dc = GpDir2d::new(
                                    origin2.x() - center.x(),
                                    origin2.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin2.XY() - Center.XY())");
                                let sign = dc.dot(&norm_l2);
                                let dc = GpDir2d::new(sign * norm_l2.x, sign * norm_l2.y)
                                    .expect("gp_Dir2d(sign * normL2.XY())");
                                me.pnttg2sol[n - 1] =
                                    GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                me.par2sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                                me.pararg2[n - 1] =
                                    parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);

                                let dc = GpDir2d::new(
                                    origin3.x() - center.x(),
                                    origin3.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin3.XY() - Center.XY())");
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
                    me.well_done = true;
                }
            }
        }
    }
    me
}
