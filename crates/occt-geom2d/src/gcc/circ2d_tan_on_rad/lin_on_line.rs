//! `GccAna_Circ2dTanOnRad(QualifiedLin, gp_Lin2d, Radius, Tolerance)`
//! (`GccAna_Circ2dTanOnRad_1.cxx:42-153`): circles of the given radius tangent
//! to a line with the centre on another line.
//!
//! `cxx:79-82` reports `WellDone = true` with no solution when the two lines are
//! parallel. Otherwise `eps` selects the side(s) of `L1` at distance `Radius`
//! (`cxx:83-101`), each side builds the parallel line `L1para` (`cxx:108`) and
//! every intersection point of `L1para` with `OnLine` is a centre
//! (`cxx:109-148`).
//!
//! Translation notes:
//! - `gp_Vec2d(Center, origin1)` of `cxx:134` is `origin1 - Center`, i.e. the
//!   shared `side(origin1, center, normal)` helper with `normal` the left normal
//!   of `L1`,
//! - `cxx:138`/`:141` store `pntcen3(NbrSol)` from `cirsol(1).Location()` (the
//!   first solution), not from `cirsol(NbrSol)`; kept verbatim,
//! - `TheSame1(NbrSol) = 0` of `cxx:133` matches the zero initialisation the
//!   port already has, and is kept as the assignment.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use super::super::circ2d_2tan_on::{
    axis_at, dir_between, left_normal, qualifier_of_line, side,
};
use super::super::gcc_ent::GccEntQualifiedLin;
use super::{GccAnaCirc2dTanOnRad, MAX_SOL};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    on_line: &GpLin2d,
    radius: f64,
    tolerance: f64,
) -> GccAnaCirc2dTanOnRad {
    let mut me = GccAnaCirc2dTanOnRad::empty();
    let tol = tolerance.abs();
    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in GccAna_Circ2dTanOnRad");
    }
    if radius < 0.0 {
        panic!("Standard_NegativeValue in GccAna_Circ2dTanOnRad");
    }

    let l1 = qualified1.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let (dx1, dy1) = (dir1.x, dir1.y);
    let norm_l1 = left_normal(&dir1);

    if on_line.direction().is_parallel(&dir1, tol) {
        me.well_done = true;
    } else {
        let eps: [f64; 2];
        let nbsol;
        if qualified1.is_enclosed() {
            eps = [-1.0, 0.0];
            nbsol = 1usize;
        } else if qualified1.is_outside() {
            eps = [1.0, 0.0];
            nbsol = 1usize;
        } else {
            eps = [1.0, -1.0];
            nbsol = 2usize;
        }
        let (lx1, ly1) = (origin1.x(), origin1.y());
        for j in 1..=nbsol {
            let l1para = GpLin2d::from_pnt_dir(
                GpPnt2d::new(
                    lx1 + eps[j - 1] * radius * dy1,
                    ly1 - eps[j - 1] * radius * dx1,
                ),
                dir1,
            );
            let mut intp = IntAna2dAnaIntersection::new();
            intp.perform_lin_lin(on_line, &l1para);
            if intp.is_done() {
                if !intp.is_empty() {
                    for i in 1..=intp.nb_points() {
                        if me.nbr_sol == MAX_SOL {
                            break;
                        }
                        me.nbr_sol += 1;
                        let n = me.nbr_sol;
                        let center = *intp.point(i).value();
                        me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                        let dc1 = dir_between(&center, &origin1);
                        me.qualifier1[n - 1] = qualifier_of_line(qualified1, dc1.dot(&norm_l1));
                        me.the_same1[n - 1] = 0;
                        if side(&origin1, &center, &norm_l1) > 0.0 {
                            me.pnttg1sol[n - 1] = GpPnt2d::new(
                                center.x() + radius * -dy1,
                                center.y() + radius * dx1,
                            );
                        } else {
                            me.pnttg1sol[n - 1] = GpPnt2d::new(
                                center.x() - radius * -dy1,
                                center.y() - radius * dx1,
                            );
                        }
                        me.pntcen3[n - 1] = me.cirsol[0].location();
                        me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                        me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
                        me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen3[n - 1]);
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}
