//! `Geom2dGcc_Circ2d3TanIter(Qualified1, Qualified2, Point3, Param1, Param2,
//! Tolerance)` (`cxx:1318-1453`): one circle, one `Geom2dGcc_QCurve` and a
//! point, dispatched to the `CiCiCu` system with the point as a radius-0
//! circle.
//!
//! Two OCCT quirks are kept: the curve is evaluated at `Ufirst(2)`, the
//! point-circle parameter (`cxx:1381`), and `dot1` / `dot3` omit the tangent
//! norm from the denominator (`cxx:1398-1418`).

use std::f64::consts::PI;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::REAL_SMALL;
use occt_math::MathVector;

use crate::gcc::circ2d_3tan::GccAnaCirc2d3Tan;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc::qualified_curve::Geom2dGccQCurve;

use super::{
    accept_circ, accept_curve, angle_or_zero, dot_ratio, point_circle, vec_between,
    Geom2dGccCirc2d3TanIter, ACCEPT_TOL, ANGULAR_TOL,
};

/// `Geom2dGcc_Circ2d3TanIter::Geom2dGcc_Circ2d3TanIter` (`cxx:1318-1453`).
#[allow(clippy::too_many_arguments)]
pub(super) fn new(
    q1: &GccEntQualifiedCirc,
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
    let c1 = q1.qualified();
    let cu2 = q2.qualified().clone();
    let c3 = point_circle(point3);
    let mut func = Geom2dGccFunctionTanCuCuCu::new_ci_ci_cu(&c1, &c3, &cu2);
    let mut umin = MathVector::with_init(1, 3, 0.0);
    let mut umax = MathVector::with_init(1, 3, 0.0);
    let mut ufirst = MathVector::with_init(1, 3, 0.0);
    let mut tol = MathVector::with_init(1, 3, 0.0);
    umin.set_value(1, 0.0);
    umin.set_value(3, curve_tool::first_parameter(&cu2));
    umin.set_value(2, 0.0);
    umax.set_value(1, 2.0 * PI);
    umax.set_value(3, curve_tool::last_parameter(&cu2));
    umax.set_value(2, 2.0 * PI);
    ufirst.set_value(1, param1);
    ufirst.set_value(2, PI);
    ufirst.set_value(3, param2);
    tol.set_value(1, ANGULAR_TOL);
    tol.set_value(2, ANGULAR_TOL);
    tol.set_value(3, curve_tool::eps_x(&cu2, tol_val));
    let x = match Geom2dGccCirc2d3TanIter::solve(&mut func, &tol, &ufirst, &umin, &umax) {
        Some(x) => x,
        None => return s,
    };
    let centre1 = c1.location();
    let r1 = c1.radius();
    let point1 = GpPnt2d::new(
        centre1.x() + r1 * x.value(1).cos(),
        centre1.y() + r1 * x.value(1).sin(),
    );
    // `cxx:1381`: `Ufirst(2)`, the point-circle parameter.
    let (point2, tan2) = curve_tool::d1(&cu2, x.value(2)).expect("undefined curve value");
    let circ = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(&point1, &point2, point3, tol_val);
    if circ.is_done() {
        let sol = circ.this_solution(1);
        let centre = sol.location();
        let tan1 = GpVec2d::new(-x.value(1).sin(), x.value(1).cos());
        let tan3 = GpVec2d::new(-x.value(3).sin(), x.value(3).cos());
        let vec1 = vec_between(&point1, &centre);
        let vec2 = vec_between(&point2, &centre);
        let vec3 = vec_between(point3, &centre);
        // `cxx:1398-1418`: no tangent norm in the denominators of 1 and 3.
        let n1 = vec1.magnitude();
        let dot1 = if n1 >= REAL_SMALL { vec1.dot(&tan1) / n1 } else { 0.0 };
        let dot2 = dot_ratio(&vec2, &tan2);
        let n3 = vec3.magnitude();
        let dot3 = if n3 >= REAL_SMALL { vec3.dot(&tan3) / n3 } else { 0.0 };
        if dot1 <= ACCEPT_TOL && dot2 <= ACCEPT_TOL && dot3 <= ACCEPT_TOL {
            let dist = centre1.distance(&centre);
            let r_sol = sol.radius();
            if accept_circ(q1.qualifier(), r1, r_sol, dist)
                && accept_curve(q2.qualifier(), angle_or_zero(&vec2, &tan2))
            {
                s.commit(
                    sol,
                    [q1.qualifier(), q2.qualifier(), GccEntPosition::NoQualifier],
                    point1,
                    x.value(1),
                    0.0,
                    point2,
                    x.value(2),
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
