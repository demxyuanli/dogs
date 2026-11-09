//! `GccAna_Circ2d3Tan(GccEnt_QualifiedCirc, GccEnt_QualifiedLin, gp_Pnt2d,
//! Tolerance)` (`GccAna_Circ2d3Tan_5.cxx:38-309`): circles tangent to one
//! circle, one straight line and one point.
//!
//! Centers come from intersecting the `GccAna_CircLin2dBisec` of `(C1, L2)`
//! with the `GccAna_LinPnt2dBisec` of `(L2, Point3)` (`cxx:94-125`); the radius
//! is fixed by the circle qualifier, checked against the line (`cxx:170-211`)
//! and finally against the point distance (`cxx:213-217`).
//!
//! Translation notes:
//! - `MaxRad = 1e10` / `MinRad = 1e-6` (`cxx:68`) reject too large and too
//!   small radii,
//! - the `NbrSol == 4` caps (`cxx:222-225`, `cxx:303-306`) match the `(1, 4)`
//!   arrays of this constructor; [`super::MAX_SOL`] only sizes storage,
//! - when the `Qualified2` branch rejects the signed distance the OCCT `ok`
//!   flag is left untouched (`cxx:170-211`), so a solution can still be stored;
//!   the port keeps that behaviour.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::{GccAnaCircLin2dBisec, GccAnaLinPnt2dBisec};
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedLin,
    point3: &GpPnt2d,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();
    let max_rad = 1e10f64;
    let min_rad = 1e-6f64;

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
        || !(qualified2.is_enclosed()
            || qualified2.is_outside()
            || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let c1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let origin2 = l2.location();
    let dir2 = *l2.direction();
    let norm_l2 = GpDir2d::new(-dir2.y, dir2.x).expect("gp_Dir2d(-dir.Y(), dir.X())");

    let mut radius = [0.0f64; 2];

    let bis1 = GccAnaCircLin2dBisec::new(&c1, &l2);
    let bis2 = GccAnaLinPnt2dBisec::new(&l2, point3);
    if bis1.is_done() && bis2.is_done() {
        let nbsolution1 = bis1.nb_solutions();
        let sol2 = bis2
            .this_solution()
            .expect("GccAna_LinPnt2dBisec::ThisSolution");
        let typ2 = sol2.arc_type();
        for i in 1..=nbsolution1 {
            let sol1 = bis1
                .this_solution(i)
                .expect("GccAna_CircLin2dBisec::ThisSolution");
            let typ1 = sol1.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
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
                        let center = *intp.point(j).value();
                        let dist1 = center.distance(&c1.location());
                        let dist2 = l2.distance(&center);
                        let dist3 = center.distance(point3);
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
                        let side = ((l2.location().x() - center.x()) * (-l2.direction().y()))
                            + ((l2.location().y() - center.y()) * l2.direction().x());
                        if qualified2.is_enclosed() && ok {
                            if side <= 0.0 {
                                for ii in 1..=nbsol1 {
                                    if (dist2 - radius[ii as usize - 1]).abs() < tol {
                                        ok = true;
                                        radius[0] = radius[ii as usize - 1];
                                    }
                                }
                            }
                        } else if qualified2.is_outside() && ok {
                            if side >= 0.0 {
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
                        if (dist3 - radius[0]).abs() <= tol && ok {
                            ok = true;
                            nbsol3 = 1;
                        }
                        if ok {
                            for k in 1..=nbsol3 {
                                if me.nbr_sol == 4 {
                                    break;
                                }
                                let rad = radius[k as usize - 1];
                                // `cxx:227-233`: too large / too small radius.
                                if rad > max_rad {
                                    break;
                                }
                                if rad.abs() < min_rad {
                                    break;
                                }

                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
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
                                me.qualifier3[n - 1] = GccEntPosition::NoQualifier;

                                if center.distance(&c1.location()) <= tolerance
                                    && (rad - r1).abs() <= tolerance
                                {
                                    me.the_same1[n - 1] = 1;
                                } else {
                                    me.the_same1[n - 1] = 0;
                                    // `cxx:280`: `dc = Center.XY() - C1.Location().XY()`.
                                    let dc = GpDir2d::new(
                                        center.x() - c1.location().x(),
                                        center.y() - c1.location().y(),
                                    )
                                    .expect("gp_Dir2d(Center.XY() - C1.Location().XY())");
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
                                    l2.location().x() - center.x(),
                                    l2.location().y() - center.y(),
                                )
                                .expect("gp_Dir2d(L2.Location().XY() - Center.XY())");
                                let sign = dc.dot(&GpDir2d::new(-l2.direction().y(), l2.direction().x())
                                    .expect("gp_Dir2d(-L2.Direction().Y(), L2.Direction().X())"));
                                let dc = GpDir2d::new(
                                    sign * -l2.direction().y(),
                                    sign * l2.direction().x(),
                                )
                                .expect("gp_Dir2d(sign * gp_XY(-L2.Direction().Y(), L2.Direction().X()))");
                                me.pnttg2sol[n - 1] =
                                    GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                me.par2sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                                me.pararg2[n - 1] =
                                    parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);
                                me.pnttg3sol[n - 1] = *point3;
                                me.par3sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg3sol[n - 1]);
                                me.pararg3[n - 1] = 0.0;
                            }
                        }
                    }
                }
                me.well_done = true;
            }
            if me.nbr_sol == 4 {
                break;
            }
        }
    }
    me
}
