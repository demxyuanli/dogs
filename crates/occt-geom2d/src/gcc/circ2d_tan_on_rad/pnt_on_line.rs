//! `GccAna_Circ2dTanOnRad(gp_Pnt2d, gp_Lin2d, Radius, Tolerance)`
//! (`GccAna_Circ2dTanOnRad_2.cxx:43-179`): circles of the given radius through
//! a point with the centre on a line.
//!
//! `cxx:70-72` reports `WellDone = true` with no solution when the point is
//! further than `Radius + Tol` from the line, and the code still falls through
//! to the three-way dispatch below. `cxx:83-108` is the exact tangency case, one
//! centre. `cxx:109-129` is the point-on-the-line case: two centres, one on each
//! side along the *line* direction. `cxx:130-178` solves the general quadratic
//! `|centre - Point1| = Radius` intersected with the line.
//!
//! Translation notes:
//! - `cxx:71-72` sets `WellDone = true` without `return`; the port keeps the
//!   fall-through, so a point further than `Radius + Tol` still reaches the
//!   general branch (whose discriminant is negative there),
//! - the two solutions of `cxx:109-129` use `gp_Ax2d(center, OnLine.Direction())`
//!   (`cxx:115`, `:117`), not the `dirx` every other branch uses,
//! - `cxx:136`/`:143` compare `A == 0.0` / `A != 0.0` exactly, without
//!   tolerance; the port keeps the exact comparison,
//! - `xbid` of `cxx:78` only exists to swap `x1`/`y1`; the port uses a direct
//!   tuple swap.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::{line_value_ax2d, parameter_circ2d, parameter_lin2d};
use occt_core::gp::ax22d::GpAx22d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::math_direct_poly_roots::DirectPolynomialRoots;

use super::super::circ2d_2tan_on::axis_at;
use super::super::gcc_ent::GccEntPosition;
use super::{GccAnaCirc2dTanOnRad, MAX_SOL};

pub(crate) fn new(
    point1: &GpPnt2d,
    on_line: &GpLin2d,
    radius: f64,
    tolerance: f64,
) -> GccAnaCirc2dTanOnRad {
    let mut me = GccAnaCirc2dTanOnRad::empty();
    let tol = tolerance.abs();
    if radius < 0.0 {
        panic!("Standard_NegativeValue in GccAna_Circ2dTanOnRad");
    }

    let dp1lin = on_line.distance(point1);
    if dp1lin > radius + tol {
        me.well_done = true;
    }

    let (x1, y1) = (point1.x(), point1.y());
    let dir = *on_line.direction();
    let (xdir, ydir) = (dir.x, dir.y);
    let loc = on_line.location();
    let (lxloc, lyloc) = (loc.x(), loc.y());

    if (dp1lin - radius).abs() < tol {
        me.well_done = true;
        me.nbr_sol = 1;
        let center = if -ydir * (x1 - lxloc) + xdir * (y1 - lyloc) < 0.0 {
            GpPnt2d::new(x1 - ydir * dp1lin, y1 + xdir * dp1lin)
        } else {
            GpPnt2d::new(x1 + ydir * dp1lin, y1 - xdir * dp1lin)
        };
        let n = me.nbr_sol;
        me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
        me.qualifier1[n - 1] = GccEntPosition::NoQualifier;
        me.the_same1[n - 1] = 0;
        me.pnttg1sol[n - 1] = *point1;
        me.pntcen3[n - 1] = me.cirsol[n - 1].location();
        me.pararg1[n - 1] = 0.0;
        me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
        me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen3[n - 1]);
    } else if dp1lin < tol {
        me.pntcen3[0] = GpPnt2d::new(x1 + radius * xdir, y1 + radius * ydir);
        me.pntcen3[1] = GpPnt2d::new(x1 - radius * xdir, y1 - radius * ydir);
        me.pntcen3[0] = line_value_ax2d(parameter_lin2d(on_line, &me.pntcen3[0]), on_line.position());
        me.pntcen3[1] = line_value_ax2d(parameter_lin2d(on_line, &me.pntcen3[1]), on_line.position());
        me.cirsol[0] = GpCirc2d::new(GpAx22d::from_xdir(me.pntcen3[0], dir), radius);
        me.cirsol[1] = GpCirc2d::new(GpAx22d::from_xdir(me.pntcen3[1], dir), radius);
        me.the_same1[0] = 0;
        me.pnttg1sol[0] = *point1;
        me.pararg1[0] = 0.0;
        me.par1sol[0] = parameter_circ2d(&me.cirsol[0], &me.pnttg1sol[0]);
        me.parcen3[0] = parameter_lin2d(on_line, &me.pntcen3[0]);
        me.the_same1[1] = 0;
        me.pnttg1sol[1] = *point1;
        me.pararg1[1] = 0.0;
        me.par1sol[1] = parameter_circ2d(&me.cirsol[1], &me.pnttg1sol[1]);
        me.parcen3[1] = parameter_lin2d(on_line, &me.pntcen3[1]);
        me.nbr_sol = 2;
    } else {
        let (a0, b0, c0) = on_line.coefficients();
        let d = a0;
        let (a, b, x0, y0) = if a0 == 0.0 {
            (b0, d, y1, x1)
        } else {
            (a0, b0, x1, y1)
        };
        if a != 0.0 {
            let sol = DirectPolynomialRoots::new3(
                (b * b + a * a) / (a * a),
                2.0 * (b * c0 / (a * a) + (b / a) * x0 - y0),
                x0 * x0 + y0 * y0 + c0 * c0 / (a * a) - radius * radius + 2.0 * c0 * x0 / a,
            );
            if sol.is_done() {
                for i in 1..=sol.nb_solutions() {
                    let center = if d != 0.0 {
                        let yc = sol.value(i);
                        GpPnt2d::new(-(b / a) * yc - c0 / a, yc)
                    } else {
                        let xc = sol.value(i);
                        GpPnt2d::new(xc, -(b / a) * xc - c0 / a)
                    };
                    if me.nbr_sol == MAX_SOL {
                        break;
                    }
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                    me.qualifier1[n - 1] = GccEntPosition::NoQualifier;
                    me.the_same1[n - 1] = 0;
                    me.pnttg1sol[n - 1] = *point1;
                    me.pntcen3[n - 1] = me.cirsol[n - 1].location();
                    me.pararg1[n - 1] = 0.0;
                    me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                    me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen3[n - 1]);
                }
                me.well_done = true;
            }
        }
    }
    me
}
