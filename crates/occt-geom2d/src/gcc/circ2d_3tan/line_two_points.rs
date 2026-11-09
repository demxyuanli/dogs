//! `GccAna_Circ2d3Tan(GccEnt_QualifiedLin, gp_Pnt2d, gp_Pnt2d, Tolerance)`
//! (`GccAna_Circ2d3Tan_8.cxx:34-177`): circles tangent to one straight line and
//! two points.
//!
//! The center is the intersection of the `GccAna_Pnt2dBisec` of the two points
//! with the `GccAna_LinPnt2dBisec` of `(L1, Point2)` (`cxx:87-101`); the radius
//! is the distance to `L1`, validated against the first point (`cxx:122-133`).
//!
//! Translation notes:
//! - coincident points return `WellDone = false` (`cxx:78-82`),
//! - unlike `_6.cxx`, the line qualifier here is a single `if` over the signed
//!   distance and leaves `ok` untouched when it fails (`cxx:112-120`),
//! - the `Intp` parallel/identical filters (`cxx:104`) are applied, and the
//!   `(1, 2)` OCCT arrays (`cxx:41-57`) are sized at [`super::MAX_SOL`].

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::precision::CONFUSION;

use crate::gcc_ana_bisec::{GccAnaLinPnt2dBisec, GccAnaPnt2dBisec};
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    point2: &GpPnt2d,
    point3: &GpPnt2d,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let l1 = qualified1.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let norm_l1 = GpDir2d::new(-dir1.y, dir1.x).expect("gp_Dir2d(-dir1.Y(), dir1.X())");

    if point2.is_equal(point3, CONFUSION) {
        return me;
    }

    let bis1 = GccAnaPnt2dBisec::new(point2, point3);
    let bis2 = GccAnaLinPnt2dBisec::new(&l1, point2);
    if bis1.is_done() && bis2.is_done() {
        let linint1 = bis1
            .this_solution()
            .expect("GccAna_Pnt2dBisec::ThisSolution");
        let sol2 = bis2
            .this_solution()
            .expect("GccAna_LinPnt2dBisec::ThisSolution");
        let typ2 = sol2.arc_type();
        let mut intp = IntAna2dAnaIntersection::new();
        if typ2 == GccIntIType::Lin {
            intp.perform_lin_lin(&linint1, &sol2.line().unwrap());
        } else if typ2 == GccIntIType::Par {
            intp.perform_lin_conic(
                &linint1,
                &IntAna2dConic::from_parab2d(&sol2.parabola().unwrap()),
            );
        }
        if intp.is_done() {
            if !intp.is_empty() && !intp.parallel_elements() && !intp.identical_elements() {
                for j in 1..=intp.nb_points() {
                    let center = *intp.point(j).value();
                    let dist1 = l1.distance(&center);
                    let dist2 = center.distance(point2);

                    let mut radius = 0.0f64;
                    let nbsol3;
                    let mut ok = false;
                    let difference = ((origin1.x() - center.x()) * (-dir1.y))
                        + ((origin1.y() - center.y()) * dir1.x);
                    if (qualified1.is_enclosed() && difference <= 0.0)
                        || (qualified1.is_outside() && difference >= 0.0)
                        || qualified1.is_unqualified()
                    {
                        ok = true;
                        radius = dist1;
                    }
                    if ok {
                        if (dist2 - radius).abs() <= tol {
                            nbsol3 = 1i32;
                        } else {
                            ok = false;
                            nbsol3 = 0;
                        }
                    } else {
                        nbsol3 = 0;
                    }
                    if ok {
                        for _k in 1..=nbsol3 {
                            if me.nbr_sol >= super::MAX_SOL {
                                break;
                            }
                            me.nbr_sol += 1;
                            let n = me.nbr_sol;
                            me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);

                            let dc1 = GpDir2d::new(origin1.x() - center.x(), origin1.y() - center.y())
                                .expect("gp_Dir2d(origin1.XY() - Center.XY())");
                            me.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                qualified1.qualifier()
                            } else if dc1.dot(&norm_l1) > 0.0 {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosed
                            };
                            me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                            me.qualifier3[n - 1] = GccEntPosition::NoQualifier;
                            me.the_same1[n - 1] = 0;
                            let dc = GpDir2d::new(origin1.x() - center.x(), origin1.y() - center.y())
                                .expect("gp_Dir2d(origin1.XY() - Center.XY())");
                            let sign = dc.dot(&norm_l1);
                            let dc = GpDir2d::new(sign * norm_l1.x, sign * norm_l1.y)
                                .expect("gp_Dir2d(sign * gp_XY(-dir1.Y(), dir1.X()))");
                            me.pnttg1sol[n - 1] =
                                GpPnt2d::new(center.x() + radius * dc.x, center.y() + radius * dc.y);
                            me.par1sol[n - 1] =
                                parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                            me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
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
            }
            me.well_done = true;
        }
    }
    me
}
