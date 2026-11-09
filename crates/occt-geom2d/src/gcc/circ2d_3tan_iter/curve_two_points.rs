//! `Geom2dGcc_Circ2d3TanIter(Qualified1, Point1, Point2, Param1, Tolerance)`
//! (`cxx:912-1034`): one `Geom2dGcc_QCurve` and two points, dispatched to the
//! `CiCiCu` system with both points as radius-0 circles.
//!
//! `Qualified2` / `Qualified3` are reported as `GccEnt_noqualifier` and both
//! `pararg2` / `pararg3` are `0.` (`cxx:1020-1029`).

use std::f64::consts::PI;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_math::MathVector;

use crate::gcc::circ2d_3tan::GccAnaCirc2d3Tan;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::gcc::qualified_curve::Geom2dGccQCurve;

use super::{
    accept_curve, angle_or_zero, dot_ratio, point_circle, vec_between, Geom2dGccCirc2d3TanIter,
    ACCEPT_TOL, ANGULAR_TOL,
};

/// `Geom2dGcc_Circ2d3TanIter::Geom2dGcc_Circ2d3TanIter` (`cxx:912-1034`).
pub(super) fn new(
    q1: &Geom2dGccQCurve,
    point2: &GpPnt2d,
    point3: &GpPnt2d,
    param1: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d3TanIter {
    let mut s = Geom2dGccCirc2d3TanIter::not_done();
    if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d3TanIter");
    }
    let tol_val = tolerance.abs();
    let c1 = point_circle(point2);
    let c2 = point_circle(point3);
    let cu1 = q1.qualified().clone();
    let mut func = Geom2dGccFunctionTanCuCuCu::new_ci_ci_cu(&c1, &c2, &cu1);
    let mut umin = MathVector::with_init(1, 3, 0.0);
    let mut umax = MathVector::with_init(1, 3, 0.0);
    let mut ufirst = MathVector::with_init(1, 3, 0.0);
    let mut tol = MathVector::with_init(1, 3, 0.0);
    umin.set_value(1, 0.0);
    umin.set_value(2, 0.0);
    umin.set_value(3, curve_tool::first_parameter(&cu1));
    umax.set_value(1, 2.0 * PI);
    umax.set_value(2, 2.0 * PI);
    umax.set_value(3, curve_tool::last_parameter(&cu1));
    ufirst.set_value(1, PI);
    ufirst.set_value(2, PI);
    ufirst.set_value(3, param1);
    tol.set_value(1, ANGULAR_TOL);
    tol.set_value(2, ANGULAR_TOL);
    tol.set_value(3, curve_tool::eps_x(&cu1, tol_val));
    let x = match Geom2dGccCirc2d3TanIter::solve(&mut func, &tol, &ufirst, &umin, &umax) {
        Some(x) => x,
        None => return s,
    };
    let (point3new, tan3) = curve_tool::d1(&cu1, x.value(3)).expect("undefined curve value");
    let circ = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(point2, point3, &point3new, tol_val);
    if circ.is_done() {
        let sol = circ.this_solution(1);
        let centre = sol.location();
        let tan2 = GpVec2d::new(-x.value(2).sin(), x.value(2).cos());
        let tan1 = GpVec2d::new(-x.value(1).sin(), x.value(1).cos());
        let vec1 = vec_between(point2, &centre);
        let vec2 = vec_between(point3, &centre);
        let vec3 = vec_between(&point3new, &centre);
        let dot1 = dot_ratio(&vec1, &tan1);
        let dot2 = dot_ratio(&vec2, &tan2);
        let dot3 = dot_ratio(&vec3, &tan3);
        if dot1 <= ACCEPT_TOL && dot2 <= ACCEPT_TOL && dot3 <= ACCEPT_TOL {
            if accept_curve(q1.qualifier(), angle_or_zero(&vec1, &tan1)) {
                s.commit(
                    sol,
                    [q1.qualifier(), GccEntPosition::NoQualifier, GccEntPosition::NoQualifier],
                    point3new,
                    x.value(3),
                    0.0,
                    *point2,
                    0.0,
                    0.0,
                    *point3,
                    0.0,
                    0.0,
                );
            }
        }
    }
    s
}
