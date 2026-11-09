//! `Geom2dGcc_Circ2d3TanIter(Qualified1, Qualified2, Point3, Param1, Param2,
//! Tolerance)` (`cxx:778-910`): two `Geom2dGcc_QCurve`s and a point,
//! dispatched to the `CiCuCu` system with the point as a radius-0 circle.
//!
//! `cxx:887-889` re-tests `Qualified1` for the second curve, so `Qualified2`
//! never gates the result; the port keeps that.

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

/// `Geom2dGcc_Circ2d3TanIter::Geom2dGcc_Circ2d3TanIter` (`cxx:778-910`).
#[allow(clippy::too_many_arguments)]
pub(super) fn new(
    q1: &Geom2dGccQCurve,
    q2: &Geom2dGccQCurve,
    point3: &GpPnt2d,
    param1: f64,
    param2: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d3TanIter {
    let mut s = Geom2dGccCirc2d3TanIter::not_done();
    if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified())
        || !(q2.is_enclosed() || q2.is_enclosing() || q2.is_outside() || q2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d3TanIter");
    }
    let tol_val = tolerance.abs();
    let c1 = point_circle(point3);
    let cu1 = q1.qualified().clone();
    let cu2 = q2.qualified().clone();
    let mut func = Geom2dGccFunctionTanCuCuCu::new_ci_cu_cu(&c1, &cu1, &cu2);
    let mut umin = MathVector::with_init(1, 3, 0.0);
    let mut umax = MathVector::with_init(1, 3, 0.0);
    let mut ufirst = MathVector::with_init(1, 3, 0.0);
    let mut tol = MathVector::with_init(1, 3, 0.0);
    umin.set_value(1, 0.0);
    umin.set_value(2, curve_tool::first_parameter(&cu1));
    umin.set_value(3, curve_tool::first_parameter(&cu2));
    umax.set_value(1, 2.0 * PI);
    umax.set_value(2, curve_tool::last_parameter(&cu1));
    umax.set_value(3, curve_tool::last_parameter(&cu2));
    ufirst.set_value(1, PI);
    ufirst.set_value(2, param1);
    ufirst.set_value(3, param2);
    tol.set_value(1, ANGULAR_TOL);
    tol.set_value(2, curve_tool::eps_x(&cu1, tol_val));
    tol.set_value(3, curve_tool::eps_x(&cu2, tol_val));
    let x = match Geom2dGccCirc2d3TanIter::solve(&mut func, &tol, &ufirst, &umin, &umax) {
        Some(x) => x,
        None => return s,
    };
    let (point1, tan1) = curve_tool::d1(&cu1, x.value(2)).expect("undefined curve value");
    let (point2, tan2) = curve_tool::d1(&cu2, x.value(3)).expect("undefined curve value");
    let circ = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(point3, &point1, &point2, tol_val);
    if circ.is_done() {
        let sol = circ.this_solution(1);
        let centre = sol.location();
        let tan3 = GpVec2d::new(-x.value(1).sin(), x.value(1).cos());
        let vec1 = vec_between(&point1, &centre);
        let vec2 = vec_between(&point2, &centre);
        let vec3 = vec_between(point3, &centre);
        let dot1 = dot_ratio(&vec1, &tan1);
        let dot2 = dot_ratio(&vec2, &tan2);
        let dot3 = dot_ratio(&vec3, &tan3);
        if dot1 <= ACCEPT_TOL && dot2 <= ACCEPT_TOL && dot3 <= ACCEPT_TOL {
            if accept_curve(q1.qualifier(), angle_or_zero(&vec1, &tan1))
                // `cxx:888`: `Qualified1` is tested again for the second curve.
                && accept_curve(q1.qualifier(), angle_or_zero(&vec2, &tan2))
            {
                s.commit(
                    sol,
                    [q1.qualifier(), q2.qualifier(), GccEntPosition::NoQualifier],
                    point1,
                    x.value(2),
                    0.0,
                    point2,
                    x.value(3),
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
