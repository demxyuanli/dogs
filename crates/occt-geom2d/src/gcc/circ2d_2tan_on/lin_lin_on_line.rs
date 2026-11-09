//! `GccAna_Circ2d2TanOn(QualifiedLin, QualifiedLin, gp_Lin2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_2.cxx:30-189`): circles tangent to two lines with the
//! centre on a straight line.
//!
//! `GccAna_Lin2dBisec` gives one bisector when `L1` and `L2` are parallel and
//! two when they intersect (`cxx:73-79`); each is intersected with `OnLine` and
//! the intersection point is kept when it is not on `L1` and its distance to
//! `L1` is below `1 / Tol` (`cxx:92-101`). The qualifier pass then rejects or
//! keeps the collected centres (`cxx:110-188`).
//!
//! Translation notes:
//! - the qualifier pass removes a rejected solution by shifting the remaining
//!   `cirsol` entries down and decrementing the loop index (`cxx:171-186`); the
//!   port keeps the same in-place compaction and only compacts `cirsol`, which
//!   is all OCCT moves,
//! - `Radius < 1.0 / Tol` (`cxx:96`) is an unbounded radius guard, not a
//!   distance test; the port keeps it verbatim,
//! - `sign1`/`sign2` are `dc.Dot(gp_Dir2d(-dir.Y(), dir.X()))`, so at
//!   `sign == 0` the `gp_Dir2d(sign * gp_XY(...))` of `cxx:144`/`:154` raises
//!   `Standard_ConstructionError`; the port's `expect` panics instead.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use crate::gcc_ana_bisec::GccAnaLin2dBisec;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use super::{axis_at, dir_between, left_normal, offset_xy, GccAnaCirc2d2TanOn};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    qualified2: &GccEntQualifiedLin,
    on_line: &GpLin2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();

    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
        || !(qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let tol = tolerance.abs();
    let l1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let dir_l1 = *l1.direction();
    let dir_l2 = *l2.direction();
    let origin_l1 = l1.location();
    let origin_l2 = l2.location();
    let bis = GccAnaLin2dBisec::new(&l1, &l2);

    if bis.is_done() {
        let nb = bis.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions");
        if nb == 1 || nb == 2 {
            for k in 1..=nb {
                let branch = bis.this_solution(k).expect("GccAna_Lin2dBisec::ThisSolution");
                let mut intp = IntAna2dAnaIntersection::new();
                intp.perform_lin_lin(&branch, on_line);
                if intp.is_done() {
                    me.well_done = true;
                    if !intp.identical_elements() && !intp.parallel_elements() && !intp.is_empty() {
                        for l in 1..=intp.nb_points() {
                            let pt = *intp.point(l).value();
                            let radius = l1.distance(&pt);
                            if !l1.contains(&pt, tol) && radius < 1.0 / tol && me.nbr_sol < 2 {
                                me.nbr_sol += 1;
                                let n = me.nbr_sol;
                                me.cirsol[n - 1] = GpCirc2d::new(axis_at(pt), radius);
                            }
                        }
                    }
                }
            }
        }
    }

    // `cxx:110-188`: qualifier pass.
    let mut i = 1usize;
    while i <= me.nbr_sol {
        let pbid = me.cirsol[i - 1].location();
        let radius = me.cirsol[i - 1].radius();
        let mut ok;

        let mut dc1 = dir_between(&origin_l1, &pbid);
        let sign1 = dc1.dot(&left_normal(&dir_l1));
        if sign1 > 0.0 {
            ok = qualified1.is_unqualified() || qualified1.is_outside();
        } else {
            ok = qualified1.is_unqualified() || qualified1.is_enclosed();
        }

        let dc2 = dir_between(&origin_l2, &pbid);
        let sign2 = dc2.dot(&left_normal(&dir_l2));
        if sign2 > 0.0 {
            ok = ok && (qualified2.is_unqualified() || qualified2.is_outside());
        } else {
            ok = ok && (qualified2.is_unqualified() || qualified2.is_enclosed());
        }

        if ok {
            dc1 = signed_normal(sign1, &left_normal(&dir_l1));
            me.pnttg1sol[i - 1] = offset_xy(&pbid, radius, &dc1);
            me.qualifier1[i - 1] = if sign1 > 0.0 {
                GccEntPosition::Outside
            } else {
                GccEntPosition::Enclosed
            };
            let dc2 = signed_normal(sign2, &left_normal(&dir_l2));
            me.pnttg2sol[i - 1] = offset_xy(&pbid, radius, &dc2);
            me.qualifier2[i - 1] = if sign2 > 0.0 {
                GccEntPosition::Outside
            } else {
                GccEntPosition::Enclosed
            };
            me.pntcen[i - 1] = pbid;
            me.par1sol[i - 1] = parameter_circ2d(&me.cirsol[i - 1], &me.pnttg1sol[i - 1]);
            me.pararg1[i - 1] = parameter_lin2d(&l1, &me.pnttg1sol[i - 1]);
            me.par2sol[i - 1] = parameter_circ2d(&me.cirsol[i - 1], &me.pnttg2sol[i - 1]);
            me.pararg2[i - 1] = parameter_lin2d(&l2, &me.pnttg2sol[i - 1]);
            me.parcen3[i - 1] = parameter_lin2d(on_line, &me.pntcen[i - 1]);
        } else if i == me.nbr_sol {
            me.nbr_sol -= 1;
        } else {
            for k in (i + 1)..=me.nbr_sol {
                me.cirsol[k - 2] = me.cirsol[k - 1];
            }
            me.nbr_sol -= 1;
            i -= 1;
        }
        i += 1;
    }
    me
}

/// `gp_Dir2d(sign * gp_XY(-dir.Y(), dir.X()))`.
fn signed_normal(sign: f64, normal: &GpDir2d) -> GpDir2d {
    GpDir2d::new(sign * normal.x, sign * normal.y).expect("gp_Dir2d from the signed normal")
}
