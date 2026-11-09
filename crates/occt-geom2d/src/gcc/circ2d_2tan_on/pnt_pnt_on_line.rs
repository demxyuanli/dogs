//! `GccAna_Circ2d2TanOn(gp_Pnt2d, gp_Pnt2d, gp_Lin2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_5.cxx:38-99`): circles through two points with the
//! centre on a straight line.
//!
//! The line of points equidistant to `Point1` and `Point2` (`cxx:70-71`) is
//! intersected with `OnLine`; every intersection point is a centre and the
//! radius is its distance to `Point1` (`cxx:78-83`). Coincident points report
//! `WellDone = true` with no solution (`cxx:65-68`).
//!
//! Translation notes:
//! - `gp_Dir2d(P1.Y() - P2.Y(), P2.X() - P1.X())` (`cxx:71`) is the direction
//!   orthogonal to `P1P2`; the port keeps the swapped operands,
//! - the two `gp_Pnt2d` arguments are stored as the tangency points and their
//!   parameters are fixed at `0.` (`cxx:85-89`).

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use super::super::gcc_ent::GccEntPosition;
use super::{axis_at, GccAnaCirc2d2TanOn, MAX_SOL};

pub(crate) fn new(
    point1: &GpPnt2d,
    point2: &GpPnt2d,
    on_line: &GpLin2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();
    let dist = point1.distance(point2);
    if dist < tolerance.abs() {
        me.well_done = true;
        return me;
    }

    let mid = GpPnt2d::new(
        (point1.x() + point2.x()) / 2.0,
        (point1.y() + point2.y()) / 2.0,
    );
    let axis_dir = GpDir2d::new(point1.y() - point2.y(), point2.x() - point1.x())
        .expect("gp_Dir2d of the perpendicular bisector");
    let l1 = GpLin2d::from_pnt_dir(mid, axis_dir);
    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_lin(&l1, on_line);
    if intp.is_done() {
        if !intp.is_empty() {
            for i in 1..=intp.nb_points() {
                if me.nbr_sol == MAX_SOL {
                    break;
                }
                me.nbr_sol += 1;
                let n = me.nbr_sol;
                let center = *intp.point(i).value();
                me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), point1.distance(&center));
                me.qualifier1[n - 1] = GccEntPosition::NoQualifier;
                me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                me.pnttg1sol[n - 1] = *point1;
                me.pnttg2sol[n - 1] = *point2;
                me.pntcen[n - 1] = me.cirsol[n - 1].location();
                me.pararg1[n - 1] = 0.0;
                me.pararg2[n - 1] = 0.0;
                me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
            }
        }
        me.well_done = true;
    }
    me
}
