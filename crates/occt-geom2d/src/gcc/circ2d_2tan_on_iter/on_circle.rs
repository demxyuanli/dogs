//! Constructors whose carrier is a `gp_Circ2d` (`Geom2dGcc_Circ2d2TanOnIter.cxx:482-918`).

use super::*;

impl Geom2dGccCirc2d2TanOnIter {
    /// `(QualifiedCirc, QCurve, gp_Circ2d)` (`cxx:482-595`).
    ///
    /// The carrier is dispatched to the `(gp_Circ2d, curve, gp_Circ2d)`
    /// constructor, whose `TheType` is `CuCuOnCi` (`cxx:538`), so the circle
    /// operand is not seen by the system.
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_cu_on_ci(
        q1: &GccEntQualifiedCirc,
        q2: &Geom2dGccQCurve,
        on_circ: &GpCirc2d,
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
        let c1 = q1.qualified();
        let r1 = c1.radius();
        let cu2 = q2.qualified();
        let mut umin = MathVector::with_init(1, 4, 0.0);
        let mut umax = MathVector::with_init(1, 4, 0.0);
        let mut ufirst = MathVector::with_init(1, 4, 0.0);
        let mut tol = MathVector::with_init(1, 4, 0.0);
        umin.set_value(1, REAL_FIRST);
        umin.set_value(2, curve_tool::first_parameter(cu2));
        umin.set_value(3, REAL_FIRST);
        umin.set_value(4, 0.0);
        umax.set_value(1, REAL_LAST);
        umax.set_value(2, curve_tool::last_parameter(cu2));
        umax.set_value(3, REAL_LAST);
        umax.set_value(4, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        ufirst.set_value(3, param3);
        tol.set_value(1, ANGULAR_TOL);
        tol.set_value(2, curve_tool::eps_x(cu2, tolerance.abs()));
        tol.set_value(3, ANGULAR_TOL);
        tol.set_value(4, tol_val / 10.);
        let point1 = circle_value_ax22d(param1, c1.position(), r1);
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = circle_value_ax22d(param3, on_circ.position(), on_circ.radius());
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_ci_cu_on_ci(
            &c1,
            cu2,
            on_circ,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1, tan1) = circle_d1_ax22d(x.value(1), c1.position(), r1);
        let (point2, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let point3 = circle_value_ax22d(x.value(3), on_circ.position(), on_circ.radius());
        let _ = tan1;
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

    /// `(QualifiedLin, QCurve, gp_Circ2d)` (`cxx:597-704`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_li_cu_on_ci(
        q1: &GccEntQualifiedLin,
        q2: &Geom2dGccQCurve,
        on_circ: &GpCirc2d,
        param1: f64,
        param2: f64,
        param3: f64,
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
        umin.set_value(3, REAL_FIRST);
        umin.set_value(4, 0.0);
        umax.set_value(1, REAL_LAST);
        umax.set_value(2, curve_tool::last_parameter(cu2));
        umax.set_value(3, REAL_LAST);
        umax.set_value(4, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        ufirst.set_value(3, param3);
        tol.set_value(1, LINEAR_TOL);
        tol.set_value(2, curve_tool::eps_x(cu2, tolerance.abs()));
        tol.set_value(3, ANGULAR_TOL);
        tol.set_value(4, tol_val / 10.);
        let point1 = line_value_ax2d(param1, l1.position());
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = circle_value_ax22d(param3, on_circ.position(), on_circ.radius());
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_li_cu_on_ci(
            &l1,
            cu2,
            on_circ,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1new, tan1) = line_d1_ax2d(x.value(1), l1.position());
        let (point2new, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let (point3new, tan3) = circle_d1_ax22d(x.value(3), on_circ.position(), on_circ.radius());
        let _ = tan1;
        let _ = tan3;
        let dist1 = point3new.distance(&point1new);
        let dist2 = point3new.distance(&point2new);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3new), (dist1 + dist2) / 2.);
            let angle2 = angle_or_zero(&vec_between(&point2new, &point3new), &tan2);
            let pscal = signed_side(&point3new, &l1);
            if accept_lin(q1.qualifier(), pscal) && accept_curve(q2.qualifier(), angle2) {
                let par1 = par_on_solution(&sol, &point1new);
                let par2 = par_on_solution(&sol, &point2new);
                s.commit(
                    sol,
                    q1.qualifier(),
                    q2.qualifier(),
                    point1new,
                    x.value(1),
                    par1,
                    point2new,
                    x.value(2),
                    par2,
                    point3new,
                    x.value(3),
                );
            }
        }
        s
    }

    /// `(QCurve, QCurve, gp_Circ2d)` (`cxx:706-827`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_cu_on_ci(
        q1: &Geom2dGccQCurve,
        q2: &Geom2dGccQCurve,
        on_circ: &GpCirc2d,
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
        umin.set_value(3, REAL_FIRST);
        umin.set_value(4, 0.0);
        umax.set_value(1, curve_tool::last_parameter(cu1));
        umax.set_value(2, curve_tool::last_parameter(cu2));
        umax.set_value(3, REAL_LAST);
        umax.set_value(4, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        ufirst.set_value(3, param3);
        tol.set_value(1, curve_tool::eps_x(cu1, tolerance.abs()));
        tol.set_value(2, curve_tool::eps_x(cu2, tolerance.abs()));
        tol.set_value(3, ANGULAR_TOL);
        tol.set_value(4, tol_val / 10.);
        let point1 = curve_tool::value(cu1, param1).expect("undefined curve value");
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let r1 = on_circ.radius();
        let loc = on_circ.location();
        let point3 = GpPnt2d::new(loc.x() + r1 * param3.cos(), loc.y() + r1 * param3.sin());
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_cu_cu_on_ci(
            cu1,
            cu2,
            on_circ,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1, tan1) = curve_tool::d1(cu1, x.value(1)).expect("undefined curve value");
        let (point2, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let u3 = x.value(3);
        let point3 = GpPnt2d::new(loc.x() + r1 * u3.cos(), loc.y() + r1 * u3.sin());
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
                    u3,
                );
            }
        }
        s
    }

    /// `(QCurve, gp_Pnt2d, gp_Circ2d)` (`cxx:829-918`).
    pub fn new_cu_pt_on_ci(
        q1: &Geom2dGccQCurve,
        point2: &GpPnt2d,
        on_circ: &GpCirc2d,
        param1: f64,
        param2: f64,
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
        umin.set_value(3, 0.0);
        umax.set_value(1, curve_tool::last_parameter(cu1));
        umax.set_value(2, REAL_LAST);
        umax.set_value(3, REAL_LAST);
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        tol.set_value(1, curve_tool::eps_x(cu1, tolerance.abs()));
        tol.set_value(2, ANGULAR_TOL);
        tol.set_value(3, tol_val / 10.);
        let point1 = curve_tool::value(cu1, param1).expect("undefined curve value");
        let point3 = circle_value_ax22d(param2, on_circ.position(), on_circ.radius());
        ufirst.set_value(3, (point3.distance(point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_cu_pt_on_ci(
            cu1,
            point2,
            on_circ,
            ufirst.value(3).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1new, tan1) = curve_tool::d1(cu1, x.value(1)).expect("undefined curve value");
        let (point3new, tan3) = circle_d1_ax22d(x.value(2), on_circ.position(), on_circ.radius());
        let _ = tan3;
        let dist1 = point3new.distance(&point1new);
        let dist2 = point3new.distance(point2);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3new), (dist1 + dist2) / 2.);
            let angle1 = angle_or_zero(&vec_between(&point1new, &point3new), &tan1);
            if accept_curve(q1.qualifier(), angle1) {
                let par1 = par_on_solution(&sol, &point1new);
                let par2 = par_on_solution(&sol, point2);
                s.commit(
                    sol,
                    q1.qualifier(),
                    GccEntPosition::NoQualifier,
                    point1new,
                    x.value(1),
                    par1,
                    *point2,
                    0.0,
                    par2,
                    point3new,
                    x.value(3),
                );
            }
        }
        s
    }
}
