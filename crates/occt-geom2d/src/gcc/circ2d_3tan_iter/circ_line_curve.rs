//! `Geom2dGcc_Circ2d3TanIter(Qualified1, Qualified2, Qualified3, Param1,
//! Param2, Param3, Tolerance)` (`cxx:1171-1315`): one circle, one line and one
//! `Geom2dGcc_QCurve`, dispatched to the `CiLiCu` system.
//!
//! `cxx:1228-1229` calls `Func.Value` before `Root.Root`; harmless because the
//! system is stateless and the port drops the call.

use std::f64::consts::PI;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_math::MathVector;

use crate::gcc::circ2d_3tan::GccAnaCirc2d3Tan;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
use crate::gcc::gcc_ent::{GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::Geom2dGccQCurve;

use super::{
    accept_circ, accept_curve, accept_lin, angle_or_zero, dot_ratio, vec_between,
    Geom2dGccCirc2d3TanIter, ACCEPT_TOL, ANGULAR_TOL, LINEAR_TOL, REAL_FIRST, REAL_LAST,
};

/// `Geom2dGcc_Circ2d3TanIter::Geom2dGcc_Circ2d3TanIter` (`cxx:1171-1315`).
#[allow(clippy::too_many_arguments)]
pub(super) fn new(
    q1: &GccEntQualifiedCirc,
    q2: &GccEntQualifiedLin,
    q3: &Geom2dGccQCurve,
    param1: f64,
    param2: f64,
    param3: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d3TanIter {
    let mut s = Geom2dGccCirc2d3TanIter::not_done();
    if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified())
        || !(q2.is_enclosed() || q2.is_outside() || q2.is_unqualified())
        || !(q3.is_enclosed() || q3.is_enclosing() || q3.is_outside() || q3.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d3TanIter");
    }
    let tol_val = tolerance.abs();
    let c1 = q1.qualified();
    let l2 = q2.qualified();
    let cu3 = q3.qualified().clone();
    let mut func = Geom2dGccFunctionTanCuCuCu::new_ci_li_cu(&c1, &l2, &cu3);
    let mut umin = MathVector::with_init(1, 3, 0.0);
    let mut umax = MathVector::with_init(1, 3, 0.0);
    let mut ufirst = MathVector::with_init(1, 3, 0.0);
    let mut tol = MathVector::with_init(1, 3, 0.0);
    umin.set_value(1, 0.0);
    umin.set_value(2, REAL_FIRST);
    umin.set_value(3, curve_tool::first_parameter(&cu3));
    umax.set_value(1, 2.0 * PI);
    umax.set_value(2, REAL_LAST);
    umax.set_value(3, curve_tool::last_parameter(&cu3));
    ufirst.set_value(1, param1);
    ufirst.set_value(2, param2);
    ufirst.set_value(3, param3);
    tol.set_value(1, ANGULAR_TOL);
    tol.set_value(2, LINEAR_TOL);
    tol.set_value(3, curve_tool::eps_x(&cu3, tol_val));
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
    let loc2 = l2.location();
    let dir2 = l2.direction();
    let point2 = GpPnt2d::new(loc2.x() + x.value(2) * dir2.x(), loc2.y() + x.value(2) * dir2.y());
    let (point3, tan3) = curve_tool::d1(&cu3, x.value(3)).expect("undefined curve value");
    let circ = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(&point1, &point2, &point3, tol_val);
    if circ.is_done() {
        let sol = circ.this_solution(1);
        let centre = sol.location();
        let tan1 = GpVec2d::new(-x.value(1).sin(), x.value(1).cos());
        let tan2 = GpVec2d::new(dir2.x(), dir2.y());
        let vec1 = vec_between(&point1, &centre);
        let vec2 = vec_between(&point2, &centre);
        let vec3 = vec_between(&point3, &centre);
        let dot1 = dot_ratio(&vec1, &tan1);
        let dot2 = dot_ratio(&vec2, &tan2);
        let dot3 = dot_ratio(&vec3, &tan3);
        if dot1 <= ACCEPT_TOL && dot2 <= ACCEPT_TOL && dot3 <= ACCEPT_TOL {
            let dist = centre1.distance(&centre);
            let r_sol = sol.radius();
            if accept_circ(q1.qualifier(), r1, r_sol, dist) {
                // `cxx:1289`: absolute centre coordinates.
                let pscal = centre.x() * (-dir2.y()) + centre.y() * dir2.x();
                if accept_lin(q2.qualifier(), pscal)
                    && accept_curve(q3.qualifier(), angle_or_zero(&vec3, &tan3))
                {
                    s.commit(
                        sol,
                        [q1.qualifier(), q2.qualifier(), q3.qualifier()],
                        point1,
                        x.value(1),
                        0.0,
                        point2,
                        x.value(2),
                        0.0,
                        point3,
                        x.value(3),
                        0.0,
                    );
                }
            }
        }
    }
    s
}
