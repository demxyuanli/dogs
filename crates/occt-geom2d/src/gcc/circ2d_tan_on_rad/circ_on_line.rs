//! `GccAna_Circ2dTanOnRad(QualifiedCirc, gp_Lin2d, Radius, Tolerance)`
//! (`GccAna_Circ2dTanOnRad.cxx:55-322`): circles of the given radius tangent to
//! a circle with the centre on a line.
//!
//! `cxx:106-201` eliminates the configurations the qualifier rules out and, in
//! the three exact boundary cases, places the single centre on the diameter of
//! `C1` normal to `OnLine` (`signe != 0`). Otherwise `Coef` holds the squared
//! offsets `|Radius -/+ R1|` and `cxx:242-318` solves the second-degree equation
//! for a centre at distance `Radius` from `C1` and on `OnLine`.
//!
//! Translation notes:
//! - `Coef(1, 2)` is 1-based in OCCT; the port keeps the 1-based `j` and
//!   indexes `coef[j - 1]`,
//! - `cxx:266-268` swaps `A`/`B` and the pivot point when `|A| <= Tol`
//!   (`D` keeps the original `A`), so the general case then solves for `xc`
//!   instead of `yc` while `C` stays unswapped; kept verbatim,
//! - `cxx:278-281` drops a root whose centre is further than `Tol` from
//!   `OnLine`; the `continue` skips the whole storage block, as in OCCT,
//! - `cxx:317` stores `pntcen3` from `cirsol(NbrSol)`, i.e. the centre of the
//!   solution just built,
//! - the qualifier of `cxx:210-225` and `:285-303` is the shared three-branch
//!   test of [`super::super::circ2d_2tan_on::qualifier_of_circle`].
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::math_direct_poly_roots::DirectPolynomialRoots;

use super::super::circ2d_2tan_on::{axis_at, dir_between, offset_xy, qualifier_of_circle};
use super::super::gcc_ent::GccEntQualifiedCirc;
use super::{GccAnaCirc2dTanOnRad, MAX_SOL};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    on_line: &GpLin2d,
    radius: f64,
    tolerance: f64,
) -> GccAnaCirc2dTanOnRad {
    let mut me = GccAnaCirc2dTanOnRad::empty();
    let tol = tolerance.abs();
    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2dTanOnRad");
    }

    let c1 = qualified1.qualified();
    if radius < 0.0 {
        panic!("Standard_NegativeValue in GccAna_Circ2dTanOnRad");
    }

    let mut coef = [0.0f64; 2];
    let mut nbsol = 0usize;
    let mut signe = 0i32;
    let mut center = GpPnt2d::zero();

    let r1 = c1.radius();
    let dist = on_line.distance(&c1.location());
    let dir = *on_line.direction();
    let (xdir, ydir) = (dir.x, dir.y);
    let loc = on_line.location();
    let (lxloc, lyloc) = (loc.x(), loc.y());
    let center1 = c1.location();
    let (x1, y1) = (center1.x(), center1.y());

    if qualified1.is_enclosed() {
        if tol < radius - r1 + dist {
            me.well_done = true;
        } else if (radius - r1 + dist).abs() < tol {
            me.well_done = true;
            me.nbr_sol = 1;
            center = boundary_center(x1, y1, lxloc, lyloc, xdir, ydir, dist);
            signe = 1;
        } else {
            coef[0] = (r1 - radius) * (r1 - radius);
            nbsol = 1;
        }
    } else if qualified1.is_enclosing() {
        if r1 + dist - radius > tol {
            me.well_done = true;
        } else if r1 + dist - radius > 0.0 {
            me.well_done = true;
            me.nbr_sol = 1;
            center = boundary_center(x1, y1, lxloc, lyloc, xdir, ydir, dist);
            signe = -1;
        } else {
            coef[0] = (radius - r1) * (radius - r1);
            nbsol = 1;
        }
    } else if dist - r1 - radius > tol {
        me.well_done = false;
    } else if (dist - r1 - radius).abs() < tol {
        me.well_done = true;
        me.nbr_sol = 1;
        center = boundary_center(x1, y1, lxloc, lyloc, xdir, ydir, dist);
        signe = -1;
    } else if qualified1.is_outside() {
        coef[0] = (radius + r1) * (radius + r1);
        nbsol = 1;
    } else {
        coef[0] = (radius - r1) * (radius - r1);
        coef[1] = (radius + r1) * (radius + r1);
        nbsol = 2;
    }

    if signe != 0 {
        let n = me.nbr_sol;
        me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
        let distcc1 = center.distance(&center1);
        me.qualifier1[n - 1] = qualifier_of_circle(qualified1, distcc1, radius, r1, tol);
        if (radius - r1).abs() <= tol {
            me.the_same1[n - 1] = 1;
        } else {
            let dir1cen = dir_between(&center1, &center);
            me.pnttg1sol[n - 1] = offset_xy(&center, signe as f64 * radius, &dir1cen);
            me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
            me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
        }
        me.pntcen3[n - 1] = me.cirsol[me.nbr_sol - 1].location();
        me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen3[n - 1]);
    } else if nbsol > 0 {
        for j in 1..=nbsol {
            let (a0, b0, c0) = on_line.coefficients();
            let d = a0;
            let (a, b, x0, y0) = if d.abs() <= tol {
                (b0, d, y1, x1)
            } else {
                (a0, b0, x1, y1)
            };
            let sol = DirectPolynomialRoots::new3(
                (b * b + a * a) / (a * a),
                2.0 * (b * c0 / (a * a) + (b / a) * x0 - y0),
                x0 * x0 + y0 * y0 + c0 * c0 / (a * a) - coef[j - 1] + 2.0 * c0 * x0 / a,
            );
            if sol.is_done() {
                for i in 1..=sol.nb_solutions() {
                    let center = if d.abs() > tol {
                        let yc = sol.value(i);
                        GpPnt2d::new(-(b / a) * yc - c0 / a, yc)
                    } else {
                        let xc = sol.value(i);
                        GpPnt2d::new(xc, -(b / a) * xc - c0 / a)
                    };
                    if on_line.distance(&center) > tol {
                        continue;
                    }
                    if me.nbr_sol == MAX_SOL {
                        break;
                    }
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                    let distcc1 = center.distance(&center1);
                    me.qualifier1[n - 1] = qualifier_of_circle(qualified1, distcc1, radius, r1, tol);
                    let dir1cen = dir_between(&center1, &center);
                    if radius > r1 || center.distance(&center1) > r1 {
                        me.pnttg1sol[n - 1] = offset_xy(&center, radius, &dir1cen);
                    } else {
                        me.pnttg1sol[n - 1] = offset_xy(&center, -radius, &dir1cen);
                    }
                    me.pntcen3[n - 1] = me.cirsol[n - 1].location();
                    me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                    me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                    me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen3[n - 1]);
                }
                me.well_done = true;
            }
        }
    }
    me
}

/// The single boundary centre of `cxx:113-182`: the point of the line through
/// `C1` normal to `OnLine`, on the side given by the sign of
/// `-ydir*(x1-lxloc) + xdir*(y1-lyloc)`.
fn boundary_center(
    x1: f64,
    y1: f64,
    lxloc: f64,
    lyloc: f64,
    xdir: f64,
    ydir: f64,
    dist: f64,
) -> GpPnt2d {
    if -ydir * (x1 - lxloc) + xdir * (y1 - lyloc) < 0.0 {
        GpPnt2d::new(x1 - ydir * dist, y1 + xdir * dist)
    } else {
        GpPnt2d::new(x1 + ydir * dist, y1 - xdir * dist)
    }
}
