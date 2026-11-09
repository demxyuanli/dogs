//! `GccAna_Circ2dTanOnRad(gp_Pnt2d, gp_Circ2d, Radius, Tolerance)`
//! (`GccAna_Circ2dTanOnRad_5.cxx:37-124`): circles of the given radius through
//! a point with the centre on a circle.
//!
//! `cxx:64-68` rejects the radii that cannot reach `OnCirc` from `Point1`.
//! `cxx:70-98` handles the two exact tangency radii (one centre, inside or
//! outside `OnCirc`); `cxx:99-123` intersects the circle of centre `Point1` and
//! radius `Radius` with `OnCirc`.
//!
//! Translation notes:
//! - the branch of `cxx:83-88` builds the centre from
//!   `gp_Dir2d(Loc.X() - signe*P.X(), Loc.Y() - signe*P.Y())` and then subtracts
//!   `Roncirc` times it from `Loc`; the port keeps the same two steps,
//! - `cxx:112-118` stores `qualifier1`, `TheSame1`, `pnttg1sol`, `pntcen3`,
//!   `par1sol` and `parcen3` at the **hard-coded** index `1`, not at `NbrSol`;
//!   with two intersections the second solution therefore overwrites indices
//!   `1` with the fields of `cirsol(1)`. The port keeps the fixed index,
//! - `par1sol(1)`/`parcen3(1)` are recomputed on every intersection from
//!   `cirsol(1)`/`pntcen3(1)`, so they always describe the first solution.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use super::super::circ2d_2tan_on::axis_at;
use super::super::gcc_ent::GccEntPosition;
use super::{GccAnaCirc2dTanOnRad, MAX_SOL};

pub(crate) fn new(
    point1: &GpPnt2d,
    on_circ: &GpCirc2d,
    radius: f64,
    tolerance: f64,
) -> GccAnaCirc2dTanOnRad {
    let mut me = GccAnaCirc2dTanOnRad::empty();
    let tol = tolerance.abs();
    let roncirc = on_circ.radius();
    let dist1 = point1.distance(&on_circ.location()) - roncirc;
    let dist2 = point1.distance(&on_circ.location()) + roncirc;

    if radius < 0.0 {
        panic!("Standard_NegativeValue in GccAna_Circ2dTanOnRad");
    }

    if dist1 - radius > tol || tol < radius - dist2 {
        me.well_done = true;
    } else {
        let mut signe = 0i32;
        if (dist1 - radius).abs() < tol {
            signe = 1;
        } else if (dist2 - radius).abs() < tol {
            signe = -1;
        }
        if signe != 0 {
            let loc = on_circ.location();
            let d = GpDir2d::new(
                loc.x() - signe as f64 * point1.x(),
                loc.y() - signe as f64 * point1.y(),
            )
            .expect("gp_Dir2d of the tangency radius");
            let center = GpPnt2d::new(loc.x() - roncirc * d.x, loc.y() - roncirc * d.y);
            me.cirsol[0] = GpCirc2d::new(axis_at(center), radius);
            me.qualifier1[0] = GccEntPosition::NoQualifier;
            me.the_same1[0] = 0;
            me.pnttg1sol[0] = *point1;
            me.pntcen3[0] = me.cirsol[0].location();
            me.pararg1[0] = 0.0;
            me.par1sol[0] = parameter_circ2d(&me.cirsol[0], &me.pnttg1sol[0]);
            me.parcen3[0] = parameter_circ2d(on_circ, &me.pntcen3[0]);
            me.well_done = true;
            me.nbr_sol = 1;
        } else {
            let mut intp = IntAna2dAnaIntersection::new();
            intp.perform_circ_circ(
                on_circ,
                &GpCirc2d::new(axis_at(*point1), radius),
            );
            if intp.is_done() {
                if !intp.is_empty() {
                    for i in 1..=intp.nb_points() {
                        if me.nbr_sol == MAX_SOL {
                            break;
                        }
                        me.nbr_sol += 1;
                        let center = *intp.point(i).value();
                        me.cirsol[me.nbr_sol - 1] = GpCirc2d::new(axis_at(center), radius);
                        me.qualifier1[0] = GccEntPosition::NoQualifier;
                        me.the_same1[0] = 0;
                        me.pnttg1sol[0] = *point1;
                        me.pntcen3[0] = me.cirsol[0].location();
                        me.par1sol[0] = parameter_circ2d(&me.cirsol[0], &me.pnttg1sol[0]);
                        me.parcen3[0] = parameter_circ2d(on_circ, &me.pntcen3[0]);
                        me.pararg1[0] = 0.0;
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}
