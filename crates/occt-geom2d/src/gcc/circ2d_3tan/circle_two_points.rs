//! `GccAna_Circ2d3Tan(GccEnt_QualifiedCirc, gp_Pnt2d, gp_Pnt2d, Tolerance)`
//! (`GccAna_Circ2d3Tan_7.cxx:36-241`): circles tangent to one circle and two
//! points.
//!
//! Centers come from intersecting the `GccAna_Pnt2dBisec` of `(Point2, Point3)`
//! with each `GccAna_CircPnt2dBisec` branch of `(C1, Point2)` (`cxx:84-112`);
//! the radius is fixed by the circle qualifier and validated against both
//! distances (`cxx:160-173`).
//!
//! Translation notes:
//! - coincident points return `WellDone = false` (`cxx:80-84`), unlike the
//!   collinear guard of `_9.cxx`,
//! - `nbsol3` keeps the *last* matching radius index (`cxx:166-172`), not the
//!   first, and the `NbrSol >= 2` break (`cxx:179-182`) leaves `WellDone`
//!   untouched only until the outer block reaches it,
//! - `dc.Reverse()` when the stored qualifier is `GccEnt_enclosed`
//!   (`cxx:215-219`) is kept,
//! - the `(1, 2)` OCCT arrays (`cxx:42-58`) are sized at [`super::MAX_SOL`].

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::precision::CONFUSION;

use crate::gcc_ana_bisec::{GccAnaCircPnt2dBisec, GccAnaPnt2dBisec};
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    point2: &GpPnt2d,
    point3: &GpPnt2d,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let c1 = qualified1.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let mut radius = [0.0f64; 2];

    // `cxx:80-84`: `Point2.IsEqual(Point3, Precision::Confusion())`.
    if point2.is_equal(point3, CONFUSION) {
        return me;
    }

    let bis1 = GccAnaPnt2dBisec::new(point2, point3);
    let bis2 = GccAnaCircPnt2dBisec::new(&c1, point2);

    if bis1.is_done() && bis2.is_done() {
        let nbsolution2 = bis2.nb_solutions();
        for i in 1..=nbsolution2 {
            let sol2 = bis2
                .this_solution(i)
                .expect("GccAna_CircPnt2dBisec::ThisSolution");
            let typ2 = sol2.arc_type();
            let sol1 = bis1
                .this_solution()
                .expect("GccAna_Pnt2dBisec::ThisSolution");
            let mut intp = IntAna2dAnaIntersection::new();
            if typ2 == GccIntIType::Cir {
                intp.perform_lin_circ(&sol1, &sol2.circle().unwrap());
            } else if typ2 == GccIntIType::Lin {
                intp.perform_lin_lin(&sol1, &sol2.line().unwrap());
            } else if typ2 == GccIntIType::Hpr {
                intp.perform_lin_conic(
                    &sol1,
                    &IntAna2dConic::from_hypr2d(&sol2.hyperbola().unwrap()),
                );
            } else if typ2 == GccIntIType::Ell {
                intp.perform_lin_conic(
                    &sol1,
                    &IntAna2dConic::from_elips2d(&sol2.ellipse().unwrap()),
                );
            }

            if intp.is_done() {
                if !intp.is_empty() {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = center.distance(&center1);
                        let dist2 = center.distance(point2);
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
                            radius[0] = r1 + dist1;
                        } else if qualified1.is_unqualified() {
                            ok = true;
                            nbsol1 = 2;
                            radius[0] = (r1 - dist1).abs();
                            radius[1] = r1 + dist1;
                        }
                        if ok {
                            ok = false;
                            for ii in 1..=nbsol1 {
                                let rad = radius[ii as usize - 1];
                                if (dist2 - rad).abs() <= tol && (dist3 - rad).abs() <= tol {
                                    nbsol3 = ii;
                                    ok = true;
                                }
                            }
                        }

                        if ok {
                            if me.nbr_sol >= 2 {
                                break;
                            }
                            me.nbr_sol += 1;
                            let n = me.nbr_sol;
                            let rad = radius[nbsol3 as usize - 1];
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
                            me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                            me.qualifier3[n - 1] = GccEntPosition::NoQualifier;
                            if center.distance(&center1) <= tolerance
                                && (rad - r1).abs() <= tolerance
                            {
                                me.the_same1[n - 1] = 1;
                            } else {
                                me.the_same1[n - 1] = 0;
                                let mut dc = GpDir2d::new(
                                    center1.x() - center.x(),
                                    center1.y() - center.y(),
                                )
                                .expect("gp_Dir2d(center1.XY() - Center.XY())");
                                // `cxx:215-219`: moving to the edge of the source circle.
                                if me.qualifier1[n - 1] == GccEntPosition::Enclosed {
                                    dc.reverse();
                                }
                                me.pnttg1sol[n - 1] =
                                    GpPnt2d::new(center.x() + rad * dc.x, center.y() + rad * dc.y);
                                me.par1sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                                me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                            }

                            me.the_same2[n - 1] = 0;
                            me.pnttg2sol[n - 1] = *point2;
                            me.par2sol[n - 1] =
                                parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                            me.pararg2[n - 1] = 0.0;
                            me.the_same3[n - 1] = 0;
                            me.pnttg3sol[n - 1] = *point3;
                            me.par3sol[n - 1] =
                                parameter_circ2d(&me.cirsol[n - 1], &me.pnttg3sol[n - 1]);
                            me.pararg3[n - 1] = 0.0;
                        }
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}
