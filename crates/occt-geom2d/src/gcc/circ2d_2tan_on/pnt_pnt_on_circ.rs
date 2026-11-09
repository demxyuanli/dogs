//! `GccAna_Circ2d2TanOn(gp_Pnt2d, gp_Pnt2d, gp_Circ2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_11.cxx:47-123`): circles through two points with the
//! centre on a circle.
//!
//! `cxx:70-81` reports `WellDone = true` with no solution when the two points
//! coincide or when one of them lies two diameters away from the other with
//! respect to `OnCirc`. Otherwise the perpendicular bisector of `Point1`/`Point2`
//! is intersected with `OnCirc` and every intersection point is a centre
//! (`cxx:82-122`).
//!
//! Translation notes:
//! - the `C1 = gp_Circ2d(...)` block of `cxx:83-95` is unreachable: it sits
//!   behind the very conditions that made `cxx:78-81` return, and its second
//!   `else if` repeats the first condition. The port keeps the block verbatim
//!   behind the same (never taken) guard,
//! - `Tol = std::abs(Tolerance)` (`cxx:71`) is the only tolerance used, the raw
//!   `Tolerance` argument elsewhere in the package is not.

use occt_core::elib::clib2d::parameter_circ2d;
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
    on_circ: &GpCirc2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();
    let tol = tolerance.abs();
    let dist = point1.distance(point2);
    let dp1cen = point1.distance(&on_circ.location());
    let dp2cen = point2.distance(&on_circ.location());
    let r = on_circ.radius();
    let mut c1 = *on_circ;
    if dist < tol
        || (dp1cen + 2.0 * r - dp2cen).abs() < tol
        || (dp2cen + 2.0 * r - dp1cen).abs() < tol
    {
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

    // `cxx:83-95`: unreachable, kept for parity.
    if (dp1cen + 2.0 * r - dp2cen).abs() < tol || (dp2cen + 2.0 * r - dp1cen).abs() < tol {
        if (dp1cen + 2.0 * r - dp2cen).abs() < tol {
            c1 = GpCirc2d::new(
                axis_at(on_circ.location()),
                on_circ.radius() + (dp2cen - dp1cen - 2.0 * r).abs(),
            );
        } else if (dp1cen + 2.0 * r - dp2cen).abs() < tol {
            c1 = GpCirc2d::new(
                axis_at(on_circ.location()),
                on_circ.radius() + (dp2cen - dp1cen - 2.0 * r).abs(),
            );
        }
    }

    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_circ(&l1, &c1);
    if intp.is_done() {
        if !intp.is_empty() {
            for i in 1..=intp.nb_points() {
                if me.nbr_sol == MAX_SOL {
                    break;
                }
                me.nbr_sol += 1;
                let n = me.nbr_sol;
                let pt = *intp.point(i).value();
                me.cirsol[n - 1] = GpCirc2d::new(axis_at(pt), point1.distance(&pt));
                me.qualifier1[n - 1] = GccEntPosition::NoQualifier;
                me.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                me.pnttg1sol[n - 1] = *point1;
                me.pararg1[n - 1] = 0.0;
                me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                me.pnttg2sol[n - 1] = *point2;
                me.pararg2[n - 1] = 0.0;
                me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                me.pntcen[n - 1] = me.cirsol[n - 1].location();
                me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
            }
        }
        me.well_done = true;
    }
    me
}
