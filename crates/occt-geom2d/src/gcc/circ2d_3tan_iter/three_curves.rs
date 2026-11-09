//! `Geom2dGcc_Circ2d3TanIter(Qualified1, Qualified2, Qualified3, Param1,
//! Param2, Param3, Tolerance)` (`cxx:44-183`): three `Geom2dGcc_QCurve`
//! operands, dispatched to the `CuCuCu` system.
//!
//! `par2sol` / `par3sol` are chord lengths from the first tangency point
//! (`cxx:174`, `:178`), the only constructor of the class that reports them.

use occt_math::MathVector;

use crate::gcc::circ2d_3tan::GccAnaCirc2d3Tan;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
use crate::gcc::qualified_curve::Geom2dGccQCurve;

use super::{accept_curve, angle_or_zero, dot_ratio, vec_between, Geom2dGccCirc2d3TanIter, ACCEPT_TOL};

/// `Geom2dGcc_Circ2d3TanIter::Geom2dGcc_Circ2d3TanIter` (`cxx:44-183`).
pub(super) fn new(
    q1: &Geom2dGccQCurve,
    q2: &Geom2dGccQCurve,
    q3: &Geom2dGccQCurve,
    param1: f64,
    param2: f64,
    param3: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d3TanIter {
    let mut s = Geom2dGccCirc2d3TanIter::not_done();
    if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified())
        || !(q2.is_enclosed() || q2.is_enclosing() || q2.is_outside() || q2.is_unqualified())
        || !(q3.is_enclosed() || q3.is_enclosing() || q3.is_outside() || q3.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d3TanIter");
    }
    let tol_val = tolerance.abs();
    let cu1 = q1.qualified().clone();
    let cu2 = q2.qualified().clone();
    let cu3 = q3.qualified().clone();
    let mut func = Geom2dGccFunctionTanCuCuCu::new_cu_cu_cu(&cu1, &cu2, &cu3);
    let mut umin = MathVector::with_init(1, 3, 0.0);
    let mut umax = MathVector::with_init(1, 3, 0.0);
    let mut ufirst = MathVector::with_init(1, 3, 0.0);
    let mut tol = MathVector::with_init(1, 3, 0.0);
    umin.set_value(1, curve_tool::first_parameter(&cu1));
    umin.set_value(2, curve_tool::first_parameter(&cu2));
    umin.set_value(3, curve_tool::first_parameter(&cu3));
    umax.set_value(1, curve_tool::last_parameter(&cu1));
    umax.set_value(2, curve_tool::last_parameter(&cu2));
    umax.set_value(3, curve_tool::last_parameter(&cu3));
    ufirst.set_value(1, param1);
    ufirst.set_value(2, param2);
    ufirst.set_value(3, param3);
    tol.set_value(1, curve_tool::eps_x(&cu1, tol_val));
    tol.set_value(2, curve_tool::eps_x(&cu2, tol_val));
    tol.set_value(3, curve_tool::eps_x(&cu3, tol_val));
    let x = match Geom2dGccCirc2d3TanIter::solve(&mut func, &tol, &ufirst, &umin, &umax) {
        Some(x) => x,
        None => return s,
    };
    let (point1, tan1) = curve_tool::d1(&cu1, x.value(1)).expect("undefined curve value");
    let (point2, tan2) = curve_tool::d1(&cu2, x.value(2)).expect("undefined curve value");
    let (point3, tan3) = curve_tool::d1(&cu3, x.value(3)).expect("undefined curve value");
    let circ = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(&point1, &point2, &point3, tol_val);
    if circ.is_done() {
        let sol = circ.this_solution(1);
        let centre = sol.location();
        let vec1 = vec_between(&point1, &centre);
        let vec2 = vec_between(&point2, &centre);
        let vec3 = vec_between(&point3, &centre);
        let dot1 = dot_ratio(&vec1, &tan1);
        let dot2 = dot_ratio(&vec2, &tan2);
        let dot3 = dot_ratio(&vec3, &tan3);
        if dot1 <= ACCEPT_TOL && dot2 <= ACCEPT_TOL && dot3 <= ACCEPT_TOL {
            if accept_curve(q1.qualifier(), angle_or_zero(&vec1, &tan1))
                && accept_curve(q2.qualifier(), angle_or_zero(&vec2, &tan2))
                && accept_curve(q3.qualifier(), angle_or_zero(&vec3, &tan3))
            {
                let par2sol = point2.distance(&point1);
                let par3sol = point3.distance(&point1);
                s.commit(
                    sol,
                    [q1.qualifier(), q2.qualifier(), q3.qualifier()],
                    point1,
                    x.value(1),
                    0.0,
                    point2,
                    x.value(2),
                    par2sol,
                    point3,
                    x.value(3),
                    par3sol,
                );
            }
        }
    }
    s
}
