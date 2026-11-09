//! `Geom2dGcc_Circ2d3TanIter(Qualified1, Qualified2, Qualified3, Param1,
//! Param2, Param3, Tolerance)` (`cxx:636-776`): two lines and one
//! `Geom2dGcc_QCurve`, dispatched to the `LiLiCu` system.
//!
//! `cxx:707` takes `pscal` from the **absolute** centre coordinates and reuses
//! the same value for `Qualified2` (`cxx:749-750`); the port keeps both.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_math::MathVector;

use crate::gcc::circ2d_3tan::GccAnaCirc2d3Tan;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
use crate::gcc::gcc_ent::GccEntQualifiedLin;
use crate::gcc::qualified_curve::Geom2dGccQCurve;

use super::{
    accept_curve, accept_lin, angle_or_zero, dot_ratio, vec_between, Geom2dGccCirc2d3TanIter,
    ACCEPT_TOL, LINEAR_TOL, REAL_FIRST, REAL_LAST,
};

/// `Geom2dGcc_Circ2d3TanIter::Geom2dGcc_Circ2d3TanIter` (`cxx:636-776`).
#[allow(clippy::too_many_arguments)]
pub(super) fn new(
    q1: &GccEntQualifiedLin,
    q2: &GccEntQualifiedLin,
    q3: &Geom2dGccQCurve,
    param1: f64,
    param2: f64,
    param3: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d3TanIter {
    let mut s = Geom2dGccCirc2d3TanIter::not_done();
    if !(q1.is_enclosed() || q1.is_outside() || q1.is_unqualified())
        || !(q2.is_enclosed() || q2.is_outside() || q2.is_unqualified())
        || !(q3.is_enclosed() || q3.is_enclosing() || q3.is_outside() || q3.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d3TanIter");
    }
    let tol_val = tolerance.abs();
    let l1 = q1.qualified();
    let l2 = q2.qualified();
    let cu3 = q3.qualified().clone();
    let mut func = Geom2dGccFunctionTanCuCuCu::new_li_li_cu(&l1, &l2, &cu3);
    let mut umin = MathVector::with_init(1, 3, 0.0);
    let mut umax = MathVector::with_init(1, 3, 0.0);
    let mut ufirst = MathVector::with_init(1, 3, 0.0);
    let mut tol = MathVector::with_init(1, 3, 0.0);
    umin.set_value(1, REAL_FIRST);
    umin.set_value(2, REAL_FIRST);
    umin.set_value(3, curve_tool::first_parameter(&cu3));
    umax.set_value(1, REAL_LAST);
    umax.set_value(2, REAL_LAST);
    umax.set_value(3, curve_tool::last_parameter(&cu3));
    ufirst.set_value(1, param1);
    ufirst.set_value(2, param2);
    ufirst.set_value(3, param3);
    tol.set_value(1, LINEAR_TOL);
    tol.set_value(2, LINEAR_TOL);
    tol.set_value(3, curve_tool::eps_x(&cu3, tol_val));
    let x = match Geom2dGccCirc2d3TanIter::solve(&mut func, &tol, &ufirst, &umin, &umax) {
        Some(x) => x,
        None => return s,
    };
    let loc1 = l1.location();
    let dir1 = l1.direction();
    let loc2 = l2.location();
    let dir2 = l2.direction();
    let point1 = GpPnt2d::new(loc1.x() + x.value(1) * dir1.x(), loc1.y() + x.value(1) * dir1.y());
    let point2 = GpPnt2d::new(loc2.x() + x.value(2) * dir2.x(), loc2.y() + x.value(2) * dir2.y());
    let (point3, tan3) = curve_tool::d1(&cu3, x.value(3)).expect("undefined curve value");
    let circ = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(&point1, &point2, &point3, tol_val);
    if circ.is_done() {
        let sol = circ.this_solution(1);
        let centre = sol.location();
        // `cxx:707`: absolute centre coordinates, no line-location offset.
        let pscal = centre.x() * (-dir1.y()) + centre.y() * dir1.x();
        if accept_lin(q1.qualifier(), pscal) {
            let tan1 = GpVec2d::new(dir1.x(), dir1.y());
            let tan2 = GpVec2d::new(dir2.x(), dir2.y());
            let vec1 = vec_between(&point1, &centre);
            let vec2 = vec_between(&point2, &centre);
            let vec3 = vec_between(&point3, &centre);
            let dot1 = dot_ratio(&vec1, &tan1);
            let dot2 = dot_ratio(&vec2, &tan2);
            let dot3 = dot_ratio(&vec3, &tan3);
            if dot1 <= ACCEPT_TOL && dot2 <= ACCEPT_TOL && dot3 <= ACCEPT_TOL {
                // `cxx:749-750` reuses the `pscal` of the first line.
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
