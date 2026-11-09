//! `GccAna_Circ2dTanOnRad(QualifiedLin, gp_Circ2d, Radius, Tolerance)`
//! (`GccAna_Circ2dTanOnRad_4.cxx:38-175`): circles of the given radius tangent
//! to a line with the centre on a circle.
//!
//! `cxx:93-97` rejects the configurations where no centre of distance `Radius`
//! from `L1` can reach `OnCirc`. `cxx:102-110` clamps `dist1`/`dist2` towards
//! `Radius`; both values are dead afterwards (kept for parity). `cxx:112-127`
//! selects the side(s) of `L1`, `cxx:132-171` intersects the parallel of `L1` at
//! distance `Radius` with `OnCirc`.
//!
//! Translation notes:
//! - the `dist1`/`dist2` stores of `cxx:102-110` are never read again; they are
//!   kept because they are part of the OCCT body, hence the
//!   `unused_assignments` allowance below,
//! - `cxx:149` overwrites the loop's `sign` with `(int)dc1.Dot(normL1)`, a
//!   truncation towards zero, so the next iteration of the `j` loop starts from
//!   that truncated value; the port keeps the overwrite,
//! - `cxx:164` builds `pnttg1sol` with `(+sign*Radius*dir1.Y(), -sign*Radius*dir1.X())`
//!   while `cxx:134` builds the parallel line with the **opposite** signs
//!   (`(-sign*Radius*dir1.Y(), +sign*Radius*dir1.X())`); both are preserved,
//! - `L1 = Qualified1.Qualified()` of `cxx:94` is a redundant re-read of
//!   `cxx:78`; dropped, since it assigns the same value to the same variable.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use super::super::circ2d_2tan_on::{axis_at, dir_between, left_normal, qualifier_of_line};
use super::super::gcc_ent::GccEntQualifiedLin;
use super::{GccAnaCirc2dTanOnRad, MAX_SOL};

#[allow(unused_assignments)]
pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    on_circ: &GpCirc2d,
    radius: f64,
    tolerance: f64,
) -> GccAnaCirc2dTanOnRad {
    let mut me = GccAnaCirc2dTanOnRad::empty();
    let tol = tolerance.abs();
    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in GccAna_Circ2dTanOnRad");
    }

    let l1 = qualified1.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let norm_l1 = left_normal(&dir1);
    let mut dist1 = l1.distance(&on_circ.location()) - on_circ.radius();
    let mut dist2 = l1.distance(&on_circ.location()) + on_circ.radius();

    if radius < 0.0 {
        panic!("Standard_NegativeValue in GccAna_Circ2dTanOnRad");
    }

    if dist1 - radius > tol || tol < radius - dist2 {
        me.well_done = true;
    } else {
        if dist1 - radius > 0.0 {
            dist1 = radius;
        } else if dist2 - radius < 0.0 {
            dist2 = radius;
        }

        let nbsol;
        let mut sign;
        if qualified1.is_enclosed() {
            nbsol = 1usize;
            sign = -1i32;
        } else if qualified1.is_outside() {
            nbsol = 1usize;
            sign = 1i32;
        } else {
            nbsol = 2usize;
            sign = 1i32;
        }
        for _j in 1..=nbsol {
            sign = -sign;
            let l = GpLin2d::from_pnt_dir(
                GpPnt2d::new(
                    origin1.x() - sign as f64 * radius * dir1.y,
                    origin1.y() + sign as f64 * radius * dir1.x,
                ),
                dir1,
            );
            let mut intp = IntAna2dAnaIntersection::new();
            intp.perform_lin_circ(&l, on_circ);
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
                        sign = dc1.dot(&norm_l1) as i32;
                        me.qualifier1[n - 1] = qualifier_of_line(qualified1, dc1.dot(&norm_l1));
                        me.pntcen3[n - 1] = me.cirsol[n - 1].location();
                        me.pnttg1sol[n - 1] = GpPnt2d::new(
                            me.pntcen3[n - 1].x() + sign as f64 * radius * dir1.y,
                            me.pntcen3[n - 1].y() - sign as f64 * radius * dir1.x,
                        );
                        me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);
                        me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                        me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen3[n - 1]);
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}
