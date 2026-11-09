//! `GccAna_Circ2d2TanOn(QualifiedLin, QualifiedLin, gp_Circ2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_9.cxx:31-379`): circles tangent to two lines with the
//! centre on a circle. The fourth constructor argument is unnamed and unused in
//! OCCT (`cxx:34`).
//!
//! `cxx:79-149` is the parallel-lines case (`D1.Angle(D2) <= gp::Resolution()`)
//! where the single bisector is intersected with `OnCirc`. `cxx:150-295` selects
//! which bisector branches and which sign the two qualifiers require, by
//! quadrants of `D1.Angle(D2)`. `cxx:296-356` intersects the selected branches
//! and filters the points by `sgn`/`s` against the bisector direction.
//! `cxx:357-378` fills the tangency points and parameters of every solution.
//!
//! Translation notes:
//! - the tolerance is not a constructor argument here: the qualifier tests of
//!   `cxx:83-99` and the acceptance test of `cxx:312-315` are sign tests only,
//! - the parallel arm overwrites `ok` from the `Qualified2` test instead of
//!   combining it with `&&` (`cxx:92-101`); the port keeps the overwrite,
//! - `cxx:317` takes the solution centre from `Intp.Point(i)` (the outer branch
//!   index) while the radius comes from `Intp.Point(l)` (`cxx:321`); the port
//!   keeps the mismatched indices,
//! - `cxx:289-294` sets `nbsol = 4` without guarding on `Bis.IsDone()`; the port
//!   keeps the assignment outside the guard,
//! - `D1.Angle(D2)` is the signed angle of `gp_Dir2d` in `(-PI, PI]`, which the
//!   `>= 0.0`/`<= 0.0` quadrant tests depend on.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::intana2d::IntAna2dAnaIntersection;
use occt_core::precision::REAL_SMALL;

use crate::gcc_ana_bisec::GccAnaLin2dBisec;

use super::super::gcc_ent::GccEntQualifiedLin;
use super::{axis_at, left_normal, offset_xy, qualifier_of_line, side, GccAnaCirc2d2TanOn, MAX_SOL};

/// `GccAna_Circ2d2TanOn(QualifiedLin, QualifiedLin, gp_Circ2d)` (`_9.cxx:31-379`).
///
/// `GccAna_Circ2d2TanOn_9.cxx:92-101` overwrites `ok` from the second qualifier
/// instead of combining it with `&&`, so the first qualifier's value is a dead
/// store; `cxx:71`/`:73`/`:76`/`:78`/`:80` are the overwritten assignments.
#[allow(unused_assignments)]
pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    qualified2: &GccEntQualifiedLin,
    on_circ: &GpCirc2d,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();

    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
        || !(qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let l1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let bis = GccAnaLin2dBisec::new(&l1, &l2);
    let mut i = 0usize;
    let mut j = 0usize;
    let mut nbsol = 0usize;
    let mut sgn = 1.0f64;
    let mut s = 1.0f64;
    let mut ok = false;
    let d1 = *l1.direction();
    let d2 = *l2.direction();
    let dnor1 = left_normal(&d1);
    let dnor2 = left_normal(&d2);
    let origin_l1 = l1.location();
    let origin_l2 = l2.location();

    // `cxx:79-149`: parallel lines.
    if d1.angle(&d2) <= REAL_SMALL {
        if qualified1.is_enclosed() {
            ok = side(&origin_l1, &origin_l2, &dnor1) <= 0.0;
        } else if qualified1.is_outside() {
            ok = side(&origin_l1, &origin_l2, &dnor1) >= 0.0;
        } else {
            ok = true;
        }
        if qualified2.is_enclosed() {
            ok = side(&origin_l1, &origin_l2, &dnor2) >= 0.0;
        } else if qualified2.is_outside() {
            ok = side(&origin_l1, &origin_l2, &dnor2) <= 0.0;
        } else {
            ok = true;
        }
        if ok {
            let branch = bis.this_solution(1).expect("GccAna_Lin2dBisec::ThisSolution");
            let mut intp = IntAna2dAnaIntersection::new();
            intp.perform_lin_circ(&branch, on_circ);
            if intp.is_done() {
                me.well_done = true;
                if !intp.is_empty() {
                    for l in 1..=intp.nb_points() {
                        if me.nbr_sol == MAX_SOL {
                            break;
                        }
                        me.nbr_sol += 1;
                        let n = me.nbr_sol;
                        let pt = *intp.point(l).value();
                        me.cirsol[n - 1] = GpCirc2d::new(axis_at(pt), l1.distance(&pt));
                        let dc1 = GpDir2d::new(origin_l1.x() - pt.x(), origin_l1.y() - pt.y())
                            .expect("gp_Dir2d from the centre to originL1");
                        let dc2 = GpDir2d::new(origin_l2.x() - pt.x(), origin_l2.y() - pt.y())
                            .expect("gp_Dir2d from the centre to originL2");
                        me.qualifier1[n - 1] = qualifier_of_line(qualified1, dc1.dot(&dnor1));
                        me.qualifier2[n - 1] = qualifier_of_line(qualified2, dc2.dot(&dnor2));
                    }
                }
            }
        }
    } else if qualified1.is_enclosed() && qualified2.is_enclosed() {
        if bis.is_done() {
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") == 2 {
                nbsol = 1;
                i = 2;
                j = 1;
                sgn = -1.0;
            }
        }
    } else if qualified1.is_enclosed() && qualified2.is_outside() {
        if bis.is_done() {
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                nbsol = 1;
                i = 1;
                j = 1;
                if d1.angle(&d2) >= 0.0 {
                    sgn = -1.0;
                }
            }
        }
    } else if qualified1.is_outside() && qualified2.is_enclosed() {
        if bis.is_done() {
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                nbsol = 1;
                i = 1;
                j = 1;
                if d1.angle(&d2) <= 0.0 {
                    sgn = -1.0;
                }
            }
        }
    } else if qualified1.is_outside() && qualified2.is_outside() {
        if bis.is_done() {
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                nbsol = 1;
                i = 2;
                j = 1;
            }
        }
    } else if qualified1.is_unqualified() && qualified2.is_enclosed() {
        if bis.is_done() {
            nbsol = 2;
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                i = 1;
                j = 2;
            }
            if d1.angle(&d2) >= 0.0 {
                s = -1.0;
            } else {
                sgn = -1.0;
            }
        }
    } else if qualified1.is_unqualified() && qualified2.is_outside() {
        if bis.is_done() {
            nbsol = 2;
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                i = 1;
                j = 2;
            }
            if d1.angle(&d2) >= 0.0 {
                s = -1.0;
                sgn = -1.0;
            }
        }
    } else if qualified1.is_enclosed() && qualified2.is_unqualified() {
        if bis.is_done() {
            nbsol = 2;
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                i = 1;
                j = 2;
            }
            if d1.angle(&d2) >= 0.0 {
                sgn = -1.0;
            } else {
                s = -1.0;
            }
        }
    } else if qualified1.is_outside() && qualified2.is_unqualified() {
        if bis.is_done() {
            nbsol = 2;
            if bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions") >= 1 {
                i = 1;
                j = 2;
            }
            if d1.angle(&d2) <= 0.0 {
                s = -1.0;
                sgn = -1.0;
            }
        }
    } else if qualified1.is_unqualified() && qualified2.is_unqualified() {
        nbsol = 4;
        i = 1;
        j = 2;
    }

    // `cxx:296-356`: intersect the selected branches and filter the points.
    if nbsol >= 1 {
        if bis.is_done() {
            let mut kk = 0i32;
            for k in i..=(i + j - 1) {
                kk += 1;
                let branch = bis.this_solution(k).expect("GccAna_Lin2dBisec::ThisSolution");
                let branch_loc = branch.location();
                let branch_dir = *branch.direction();
                let mut intp = IntAna2dAnaIntersection::new();
                intp.perform_lin_circ(&branch, on_circ);
                if intp.is_done() {
                    if !intp.is_empty() {
                        for l in 1..=intp.nb_points() {
                            let pt_l = *intp.point(l).value();
                            // `gp_Vec2d V(Intp.Point(l).Value(), Bis....Location())`.
                            let v_dot_dir = side(&branch_loc, &pt_l, &branch_dir);
                            if (kk == 1 && sgn * v_dot_dir >= 0.0)
                                || (kk == 2 && sgn * s * v_dot_dir >= 0.0)
                                || nbsol == 4
                            {
                                if me.nbr_sol == MAX_SOL {
                                    break;
                                }
                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
                                // `cxx:317`: the centre comes from `Intp.Point(i)`.
                                let pt = *intp.point(i as i32).value();
                                me.cirsol[n - 1] = GpCirc2d::new(axis_at(pt), l1.distance(&pt_l));
                                let dc1 = GpDir2d::new(origin_l1.x() - pt.x(), origin_l1.y() - pt.y())
                                    .expect("gp_Dir2d from the centre to originL1");
                                let dc2 = GpDir2d::new(origin_l2.x() - pt.x(), origin_l2.y() - pt.y())
                                    .expect("gp_Dir2d from the centre to originL2");
                                me.qualifier1[n - 1] =
                                    qualifier_of_line(qualified1, dc1.dot(&dnor1));
                                me.qualifier2[n - 1] =
                                    qualifier_of_line(qualified2, dc2.dot(&dnor2));
                            }
                        }
                    }
                    me.well_done = true;
                }
            }
        }
    }

    // `cxx:357-378`: tangency points and parameters.
    let mut i = 1usize;
    while i <= me.nbr_sol {
        let pbid = me.cirsol[i - 1].location();
        let radius = me.cirsol[i - 1].radius();
        let dc2 = GpDir2d::new(origin_l1.x() - pbid.x(), origin_l1.y() - pbid.y())
            .expect("gp_Dir2d from the centre to originL1");
        let sign = dc2.dot(&left_normal(&d1));
        let dc2 = GpDir2d::new(sign * left_normal(&d1).x, sign * left_normal(&d1).y)
            .expect("gp_Dir2d from the signed normal");
        me.pnttg1sol[i - 1] = offset_xy(&pbid, radius, &dc2);
        let dc2 = GpDir2d::new(origin_l2.x() - pbid.x(), origin_l2.y() - pbid.y())
            .expect("gp_Dir2d from the centre to originL2");
        let sign = dc2.dot(&left_normal(&d2));
        let dc2 = GpDir2d::new(sign * left_normal(&d2).x, sign * left_normal(&d2).y)
            .expect("gp_Dir2d from the signed normal");
        me.pnttg2sol[i - 1] = offset_xy(&pbid, radius, &dc2);
        me.pntcen[i - 1] = pbid;
        me.par1sol[i - 1] = parameter_circ2d(&me.cirsol[i - 1], &me.pnttg1sol[i - 1]);
        me.pararg1[i - 1] = parameter_lin2d(&l1, &me.pnttg1sol[i - 1]);
        me.par2sol[i - 1] = parameter_circ2d(&me.cirsol[i - 1], &me.pnttg2sol[i - 1]);
        me.pararg2[i - 1] = parameter_lin2d(&l2, &me.pnttg2sol[i - 1]);
        me.parcen3[i - 1] = parameter_circ2d(on_circ, &me.pntcen[i - 1]);
        i += 1;
    }
    me
}
