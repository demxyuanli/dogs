//! `GccAna_Circ2d3Tan(GccEnt_QualifiedLin, GccEnt_QualifiedLin, gp_Pnt2d,
//! Tolerance)` (`GccAna_Circ2d3Tan_6.cxx:34-251`): circles tangent to two
//! straight lines and one point.
//!
//! Centers come from intersecting the `GccAna_Lin2dBisec` of `(L1, L2)` with
//! the `GccAna_LinPnt2dBisec` of `(L1, Point3)` (`cxx:84-102`); the radius is
//! the distance to `L1` and is validated against `L2` (`cxx:141-181`) and the
//! point (`cxx:182-191`).
//!
//! Translation notes:
//! - `pnttg3sol.Init(Point3)` (`cxx:69`) fills every tangency-3 point before
//!   the loop; `par3sol` is evaluated from that field (`cxx:241`),
//! - the `Qualified2.IsEnclosed()` branch (`cxx:141-155`) tests the signed
//!   distance without first checking `ok`, unlike the two other branches; the
//!   port keeps both that asymmetry and the sign-only rejection that leaves
//!   `ok` untouched,
//! - the `(1, 2)` OCCT arrays (`cxx:40-56`) are not clamped by the source, the
//!   port only stops at [`super::MAX_SOL`] to stay memory safe.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::{GccAnaLin2dBisec, GccAnaLinPnt2dBisec};
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    qualified2: &GccEntQualifiedLin,
    point3: &GpPnt2d,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
        || !(qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    for slot in me.pnttg3sol.iter_mut() {
        *slot = *point3;
    }

    let l1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let norm_l1 = GpDir2d::new(-dir1.y, dir1.x).expect("gp_Dir2d(-dir1.Y(), dir1.X())");
    let origin2 = l2.location();
    let dir2 = *l2.direction();
    let norm_l2 = GpDir2d::new(-dir2.y, dir2.x).expect("gp_Dir2d(-dir2.Y(), dir2.X())");

    let bis1 = GccAnaLin2dBisec::new(&l1, &l2);
    let bis2 = GccAnaLinPnt2dBisec::new(&l1, point3);
    if bis1.is_done() && bis2.is_done() {
        let nbsolution1 = bis1.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions");
        let sol2 = bis2
            .this_solution()
            .expect("GccAna_LinPnt2dBisec::ThisSolution");
        let typ2 = sol2.arc_type();
        for i in 1..=nbsolution1 {
            let sol1 = bis1
                .this_solution(i)
                .expect("GccAna_Lin2dBisec::ThisSolution");
            let mut intp = IntAna2dAnaIntersection::new();
            if typ2 == GccIntIType::Lin {
                intp.perform_lin_lin(&sol1, &sol2.line().unwrap());
            } else if typ2 == GccIntIType::Par {
                intp.perform_lin_conic(
                    &sol1,
                    &IntAna2dConic::from_parab2d(&sol2.parabola().unwrap()),
                );
            }
            if intp.is_done() {
                if !intp.is_empty() {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = l1.distance(&center);
                        let dist2 = l2.distance(&center);
                        let dist3 = center.distance(point3);
                        let mut radius = 0.0f64;
                        let mut nbsol3 = 0i32;
                        let mut ok = false;
                        if qualified1.is_enclosed() {
                            if ((origin1.x() - center.x()) * (-dir1.y))
                                + ((origin1.y() - center.y()) * dir1.x)
                                <= 0.0
                            {
                                ok = true;
                                radius = dist1;
                            }
                        } else if qualified1.is_outside() {
                            if ((origin1.x() - center.x()) * (-dir1.y))
                                + ((origin1.y() - center.y()) * dir1.x)
                                >= 0.0
                            {
                                ok = true;
                                radius = dist1;
                            }
                        } else if qualified1.is_unqualified() {
                            ok = true;
                            radius = dist1;
                        }
                        if qualified2.is_enclosed() {
                            if ((origin2.x() - center.x()) * (-dir2.y))
                                + ((origin2.y() - center.y()) * dir2.x)
                                <= 0.0
                            {
                                if (dist2 - radius).abs() < tol {
                                } else {
                                    ok = false;
                                }
                            }
                        } else if qualified2.is_outside() && ok {
                            if ((origin2.x() - center.x()) * (-dir2.y))
                                + ((origin2.y() - center.y()) * dir2.x)
                                >= 0.0
                            {
                                if (dist2 - radius).abs() < tol {
                                } else {
                                    ok = false;
                                }
                            }
                        } else if qualified2.is_unqualified() && ok {
                            if (dist2 - radius).abs() < tol {
                            } else {
                                ok = false;
                            }
                        }
                        if ok {
                            if (dist3 - radius).abs() < tol {
                                nbsol3 = 1;
                            } else {
                                ok = false;
                            }
                        }
                        if ok {
                            for _k in 1..=nbsol3 {
                                if me.nbr_sol >= super::MAX_SOL {
                                    break;
                                }
                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
                                me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);

                                let dc1 = GpDir2d::new(
                                    origin1.x() - center.x(),
                                    origin1.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin1.XY() - Center.XY())");
                                me.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                    qualified1.qualifier()
                                } else if dc1.dot(&norm_l1) > 0.0 {
                                    GccEntPosition::Outside
                                } else {
                                    GccEntPosition::Enclosed
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
                                me.qualifier3[n - 1] = GccEntPosition::NoQualifier;
                                me.the_same1[n - 1] = 0;
                                let dc = GpDir2d::new(
                                    origin1.x() - center.x(),
                                    origin1.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin1.XY() - Center.XY())");
                                let sign = dc.dot(&norm_l1);
                                let dc = GpDir2d::new(sign * norm_l1.x, sign * norm_l1.y)
                                    .expect("gp_Dir2d(sign * gp_XY(-dir1.Y(), dir1.X()))");
                                me.pnttg1sol[n - 1] =
                                    GpPnt2d::new(center.x() + radius * dc.x, center.y() + radius * dc.y);
                                me.par1sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                                me.pararg1[n - 1] =
                                    parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
                                me.the_same2[n - 1] = 0;
                                let dc = GpDir2d::new(
                                    origin2.x() - center.x(),
                                    origin2.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin2.XY() - Center.XY())");
                                let sign = dc.dot(&norm_l2);
                                let dc = GpDir2d::new(sign * norm_l2.x, sign * norm_l2.y)
                                    .expect("gp_Dir2d(sign * gp_XY(-dir2.Y(), dir2.X()))");
                                me.pnttg2sol[n - 1] =
                                    GpPnt2d::new(center.x() + radius * dc.x, center.y() + radius * dc.y);
                                me.par2sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                                me.pararg2[n - 1] =
                                    parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);
                                me.the_same3[n - 1] = 0;
                                me.par3sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg3sol[n - 1]);
                                me.pararg3[n - 1] = 0.0;
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
