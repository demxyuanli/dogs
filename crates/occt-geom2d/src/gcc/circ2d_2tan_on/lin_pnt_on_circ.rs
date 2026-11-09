//! `GccAna_Circ2d2TanOn(QualifiedLin, gp_Pnt2d, gp_Circ2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_10.cxx:47-235`): circles tangent to a line and passing
//! through a point, with the centre on a circle.
//!
//! The boundary block (`cxx:84-146`) puts the centre at `Point2 + (distpc + Ron)
//! * dir` and accepts it when the distance to `L1` matches. The general case
//! (`cxx:152-235`) intersects the single `GccAna_LinPnt2dBisec` branch of
//! `(L1, Point2)` with `OnCirc`.
//!
//! Translation notes:
//! - `cxx:88-92` builds `dir` from `OnCirc.Location() - Point2` and then
//!   reassigns it to the `L1` normal at `cxx:94`, so the later side tests use
//!   the normal; the port rebinds the same way,
//! - the general case uses `dist1` (`cxx:222`), which still holds
//!   `L1.Distance(pinterm)` from the boundary block (`cxx:90`), instead of
//!   `Radius`; the port keeps the stale value,
//! - `cxx:139-140` assigns `parcen3` from `OnCirc` and immediately overwrites it
//!   with `0.`; both assignments are kept,
//! - `qualifier2` is written twice in the general case (`cxx:219`, `:229`);
//!   both are `GccEnt_noqualifier`.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::GccAnaLinPnt2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use super::{
    axis_at, dir_between, left_normal, offset_xy, qualifier_of_line, side, GccAnaCirc2d2TanOn,
    MAX_SOL,
};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    point2: &GpPnt2d,
    on_circ: &GpCirc2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();

    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let tol = tolerance.abs();
    let l1 = qualified1.qualified();
    let origin_l1 = l1.location();
    let dir_l1 = *l1.direction();
    let norm_l1 = left_normal(&dir_l1);

    // `cxx:84-146`: boundary case.
    let ron = on_circ.radius();
    let distpc = on_circ.location().distance(point2);
    let dir = dir_between(point2, &on_circ.location());
    let pinterm = offset_xy(point2, distpc + ron, &dir);
    let dist1 = l1.distance(&pinterm);
    if (dist1 - distpc + ron).abs() <= tol {
        let direc = dir_between(&pinterm, &origin_l1);
        if qualified1.is_outside() {
            if direc.dot(&norm_l1) >= 0.0 {
                me.well_done = true;
            }
        } else if qualified1.is_enclosed() {
            if direc.dot(&norm_l1) <= 0.0 {
                me.well_done = true;
            }
        } else {
            me.well_done = true;
        }
        if me.well_done {
            me.nbr_sol = 1;
            let n = me.nbr_sol;
            me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), dist1);

            let dc1 = dir_between(&pinterm, &origin_l1);
            let sign = dc1.dot(&norm_l1);
            me.qualifier1[n - 1] = qualifier_of_line(qualified1, dc1.dot(&norm_l1));
            me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
            let dc1 = signed_normal(sign, &norm_l1);
            me.pnttg1sol[n - 1] = offset_xy(&pinterm, dist1, &dc1);
            me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
            me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
            me.pntcen[n - 1] = pinterm;
            me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
            // `cxx:140`: overwritten straight after.
            me.parcen3[n - 1] = 0.0;
            me.pnttg2sol[n - 1] = *point2;
            me.pararg2[n - 1] = 0.0;
            me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
            return me;
        }
    }

    // `cxx:152-235`: general case.
    let bis = GccAnaLinPnt2dBisec::new(&l1, point2);
    if bis.is_done() {
        let sol = bis
            .this_solution()
            .expect("GccAna_LinPnt2dBisec::ThisSolution");
        let typ = sol.arc_type();
        let mut intp = IntAna2dAnaIntersection::new();
        if typ == GccIntIType::Lin {
            intp.perform_lin_circ(&sol.line().unwrap(), on_circ);
        }
        if typ == GccIntIType::Par {
            intp.perform_circ_conic(
                on_circ,
                &IntAna2dConic::from_parab2d(&sol.parabola().unwrap()),
            );
        }
        if intp.is_done() {
            if !intp.is_empty() {
                for j in 1..=intp.nb_points() {
                    let center = *intp.point(j).value();
                    let radius = l1.distance(&center);
                    let mut ok = false;
                    if qualified1.is_enclosed() {
                        if side(&origin_l1, &center, &norm_l1) <= 0.0 {
                            ok = true;
                        }
                    } else if qualified1.is_outside() {
                        if side(&origin_l1, &center, &norm_l1) >= 0.0 {
                            ok = true;
                        }
                    } else if qualified1.is_unqualified() {
                        ok = true;
                    }
                    if ok {
                        if me.nbr_sol == MAX_SOL {
                            break;
                        }
                        me.nbr_sol += 1;
                        let n = me.nbr_sol;
                        me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                        me.the_same1[n - 1] = 0;
                        me.the_same2[n - 1] = 0;
                        let dc1 = dir_between(&center, &origin_l1);
                        let sign = dc1.dot(&norm_l1);
                        me.qualifier1[n - 1] = qualifier_of_line(qualified1, dc1.dot(&norm_l1));
                        me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                        let dc1 = signed_normal(sign, &norm_l1);
                        // `cxx:222`: the stale `dist1` of the boundary block.
                        me.pnttg1sol[n - 1] = offset_xy(&center, dist1, &dc1);
                        me.par1sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                        me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
                        me.pntcen[n - 1] = center;
                        me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
                        me.pnttg2sol[n - 1] = *point2;
                        me.pararg2[n - 1] = 0.0;
                        me.par2sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                        me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                    }
                }
            }
            me.well_done = true;
        }
    }
    me
}

/// `gp_Dir2d(sign * gp_XY(-dirL1.Y(), dirL1.X()))`.
fn signed_normal(sign: f64, normal: &GpDir2d) -> GpDir2d {
    GpDir2d::new(sign * normal.x, sign * normal.y).expect("gp_Dir2d from the signed normal")
}
