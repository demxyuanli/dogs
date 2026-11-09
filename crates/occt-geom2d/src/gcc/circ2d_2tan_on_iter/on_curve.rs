//! Constructors whose carrier is a `Geom2dAdaptor_Curve`
//! (`Geom2dGcc_Circ2d2TanOnIter.cxx:920-1344`).

use super::*;

impl Geom2dGccCirc2d2TanOnIter {
    /// `(QCurve, QCurve, Geom2dAdaptor_Curve)` (`cxx:920-1035`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_cu_on_cu(
        q1: &Geom2dGccQCurve,
        q2: &Geom2dGccQCurve,
        on_curv: &Geom2dAdaptorCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if !(q1.is_enclosed()
            || q1.is_enclosing()
            || q1.is_outside()
            || q1.is_unqualified())
            || !(q2.is_enclosed()
                || q2.is_enclosing()
                || q2.is_outside()
                || q2.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnIter");
        }
        let tol_val = tolerance.abs();
        let cu1 = q1.qualified();
        let cu2 = q2.qualified();
        let mut umin = MathVector::with_init(1, 4, 0.0);
        let mut umax = MathVector::with_init(1, 4, 0.0);
        let mut ufirst = MathVector::with_init(1, 4, 0.0);
        let mut tol = MathVector::with_init(1, 4, 0.0);
        umin.set_value(1, curve_tool::first_parameter(cu1));
        umin.set_value(2, curve_tool::first_parameter(cu2));
        umin.set_value(3, curve_tool::first_parameter(on_curv));
        umin.set_value(4, 0.0);
        umax.set_value(1, curve_tool::last_parameter(cu1));
        umax.set_value(2, curve_tool::last_parameter(cu2));
        umax.set_value(3, curve_tool::last_parameter(on_curv));
        umax.set_value(4, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        ufirst.set_value(3, param3);
        tol.set_value(1, curve_tool::eps_x(cu1, tolerance.abs()));
        tol.set_value(2, curve_tool::eps_x(cu2, tolerance.abs()));
        tol.set_value(3, curve_tool::eps_x(on_curv, tolerance.abs()));
        tol.set_value(4, tol_val / 10.);
        let point1 = curve_tool::value(cu1, param1).expect("undefined curve value");
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = curve_tool::value(on_curv, param3).expect("undefined curve value");
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_cu_cu_on_cu(
            cu1,
            cu2,
            on_curv,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1, tan1) = curve_tool::d1(cu1, x.value(1)).expect("undefined curve value");
        let (point2, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let (point3, tan3) = curve_tool::d1(on_curv, x.value(3)).expect("undefined curve value");
        let _ = tan3;
        let dist1 = point3.distance(&point1);
        let dist2 = point3.distance(&point2);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3), (dist1 + dist2) / 2.);
            let angle1 = angle_or_zero(&vec_between(&point1, &point3), &tan1);
            let angle2 = angle_or_zero(&vec_between(&point2, &point3), &tan2);
            if accept_curve(q1.qualifier(), angle1) && accept_curve(q2.qualifier(), angle2) {
                s.commit(
                    sol,
                    q1.qualifier(),
                    q2.qualifier(),
                    point1,
                    x.value(1),
                    0.0,
                    point2,
                    x.value(2),
                    point2.distance(&point1),
                    point3,
                    x.value(3),
                );
            }
        }
        s
    }

    /// `(QualifiedCirc, QCurve, Geom2dAdaptor_Curve)` (`cxx:1037-1146`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_cu_on_cu(
        q1: &GccEntQualifiedCirc,
        q2: &Geom2dGccQCurve,
        on_curv: &Geom2dAdaptorCurve,
        param1: f64,
        param2: f64,
        param_on: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if !(q1.is_enclosed()
            || q1.is_enclosing()
            || q1.is_outside()
            || q1.is_unqualified())
            || !(q2.is_enclosed()
                || q2.is_enclosing()
                || q2.is_outside()
                || q2.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnIter");
        }
        let tol_val = tolerance.abs();
        let c1 = q1.qualified();
        let r1 = c1.radius();
        let cu2 = q2.qualified();
        let mut umin = MathVector::with_init(1, 4, 0.0);
        let mut umax = MathVector::with_init(1, 4, 0.0);
        let mut ufirst = MathVector::with_init(1, 4, 0.0);
        let mut tol = MathVector::with_init(1, 4, 0.0);
        umin.set_value(1, REAL_FIRST);
        umin.set_value(2, curve_tool::first_parameter(cu2));
        umin.set_value(3, curve_tool::first_parameter(on_curv));
        umin.set_value(4, 0.0);
        umax.set_value(1, REAL_LAST);
        umax.set_value(2, curve_tool::last_parameter(cu2));
        umax.set_value(3, curve_tool::last_parameter(on_curv));
        umax.set_value(4, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        ufirst.set_value(3, param_on);
        tol.set_value(1, ANGULAR_TOL);
        tol.set_value(2, curve_tool::eps_x(cu2, tolerance.abs()));
        tol.set_value(3, curve_tool::eps_x(on_curv, tolerance.abs()));
        tol.set_value(4, tol_val / 10.);
        let point1 = circle_value_ax22d(param1, c1.position(), r1);
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = curve_tool::value(on_curv, param_on).expect("undefined curve value");
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_ci_cu_on_cu(
            &c1,
            cu2,
            on_curv,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point2, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let (point3, tan3) = curve_tool::d1(on_curv, x.value(3)).expect("undefined curve value");
        let (point1, tan1) = circle_d1_ax22d(x.value(1), c1.position(), r1);
        let _ = tan1;
        let _ = tan3;
        let dist1 = point3.distance(&point1);
        let dist2 = point3.distance(&point2);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3), (dist1 + dist2) / 2.);
            let angle2 = angle_or_zero(&vec_between(&point2, &point3), &tan2);
            let dist = c1.location().distance(&point3);
            let rsol = sol.radius();
            if accept_circ(q1.qualifier(), r1, rsol, dist) && accept_curve(q2.qualifier(), angle2) {
                let par1 = par_on_solution(&sol, &point1);
                let par2 = par_on_solution(&sol, &point2);
                s.commit(
                    sol,
                    q1.qualifier(),
                    q2.qualifier(),
                    point1,
                    x.value(1),
                    par1,
                    point2,
                    x.value(2),
                    par2,
                    point3,
                    x.value(3),
                );
            }
        }
        s
    }

    /// `(QualifiedLin, QCurve, Geom2dAdaptor_Curve)` (`cxx:1148-1253`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_li_cu_on_cu(
        q1: &GccEntQualifiedLin,
        q2: &Geom2dGccQCurve,
        on_curv: &Geom2dAdaptorCurve,
        param1: f64,
        param2: f64,
        param_on: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if !(q1.is_enclosed() || q1.is_outside() || q1.is_unqualified())
            || !(q2.is_enclosed()
                || q2.is_enclosing()
                || q2.is_outside()
                || q2.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnIter");
        }
        let tol_val = tolerance.abs();
        let l1 = q1.qualified();
        let cu2 = q2.qualified();
        let mut umin = MathVector::with_init(1, 4, 0.0);
        let mut umax = MathVector::with_init(1, 4, 0.0);
        let mut ufirst = MathVector::with_init(1, 4, 0.0);
        let mut tol = MathVector::with_init(1, 4, 0.0);
        umin.set_value(1, REAL_FIRST);
        umin.set_value(2, curve_tool::first_parameter(cu2));
        umin.set_value(3, curve_tool::first_parameter(on_curv));
        umin.set_value(4, 0.0);
        umax.set_value(1, REAL_LAST);
        umax.set_value(2, curve_tool::last_parameter(cu2));
        umax.set_value(3, curve_tool::last_parameter(on_curv));
        umax.set_value(4, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        ufirst.set_value(3, param_on);
        tol.set_value(1, LINEAR_TOL);
        tol.set_value(2, curve_tool::eps_x(cu2, tolerance.abs()));
        tol.set_value(3, curve_tool::eps_x(on_curv, tolerance.abs()));
        tol.set_value(4, tol_val / 10.);
        let point1 = line_value_ax2d(param1, l1.position());
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = curve_tool::value(on_curv, param_on).expect("undefined curve value");
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_li_cu_on_cu(
            &l1,
            cu2,
            on_curv,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1, tan1) = line_d1_ax2d(x.value(1), l1.position());
        let (point2, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let (point3, tan3) = curve_tool::d1(on_curv, x.value(3)).expect("undefined curve value");
        let _ = tan1;
        let _ = tan3;
        let dist1 = point3.distance(&point1);
        let dist2 = point3.distance(&point2);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3), (dist1 + dist2) / 2.);
            let angle2 = angle_or_zero(&vec_between(&point2, &point3), &tan2);
            let pscal = signed_side(&point3, &l1);
            if accept_lin(q1.qualifier(), pscal) && accept_curve(q2.qualifier(), angle2) {
                let par1 = par_on_solution(&sol, &point1);
                let par2 = par_on_solution(&sol, &point2);
                s.commit(
                    sol,
                    q1.qualifier(),
                    q2.qualifier(),
                    point1,
                    x.value(1),
                    par1,
                    point2,
                    x.value(2),
                    par2,
                    point3,
                    x.value(3),
                );
            }
        }
        s
    }

    /// `(QCurve, gp_Pnt2d, Geom2dAdaptor_Curve)` (`cxx:1255-1344`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_pt_on_cu(
        q1: &Geom2dGccQCurve,
        point2: &GpPnt2d,
        on_curv: &Geom2dAdaptorCurve,
        param1: f64,
        param_on: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if !(q1.is_enclosed()
            || q1.is_enclosing()
            || q1.is_outside()
            || q1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnIter");
        }
        let tol_val = tolerance.abs();
        let cu1 = q1.qualified();
        let mut umin = MathVector::with_init(1, 3, 0.0);
        let mut umax = MathVector::with_init(1, 3, 0.0);
        let mut ufirst = MathVector::with_init(1, 3, 0.0);
        let mut tol = MathVector::with_init(1, 3, 0.0);
        umin.set_value(1, curve_tool::first_parameter(cu1));
        umin.set_value(2, REAL_FIRST);
        umin.set_value(3, curve_tool::first_parameter(on_curv));
        umax.set_value(1, curve_tool::last_parameter(cu1));
        umax.set_value(2, REAL_LAST);
        umax.set_value(3, curve_tool::last_parameter(on_curv));
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param_on);
        tol.set_value(1, curve_tool::eps_x(cu1, tolerance.abs()));
        tol.set_value(2, curve_tool::eps_x(on_curv, tolerance.abs()));
        tol.set_value(3, tol_val / 10.);
        let point1 = curve_tool::value(cu1, param1).expect("undefined curve value");
        let point3 = curve_tool::value(on_curv, param_on).expect("undefined curve value");
        ufirst.set_value(3, (point3.distance(point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_cu_pt_on_cu(
            cu1,
            point2,
            on_curv,
            ufirst.value(3).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1, tan1) = curve_tool::d1(cu1, x.value(1)).expect("undefined curve value");
        let (point3, tan3) = curve_tool::d1(on_curv, x.value(3)).expect("undefined curve value");
        let _ = tan3;
        let dist1 = point3.distance(&point1);
        let dist2 = point3.distance(point2);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3), (dist1 + dist2) / 2.);
            let angle1 = angle_or_zero(&vec_between(&point1, &point3), &tan1);
            if accept_curve(q1.qualifier(), angle1) {
                let par1 = par_on_solution(&sol, &point1);
                let par2 = par_on_solution(&sol, point2);
                s.commit(
                    sol,
                    q1.qualifier(),
                    GccEntPosition::NoQualifier,
                    point1,
                    x.value(1),
                    par1,
                    *point2,
                    0.0,
                    par2,
                    point3,
                    x.value(3),
                );
            }
        }
        s
    }
}
