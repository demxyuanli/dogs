//! Constructors whose carrier is a `gp_Lin2d` (`Geom2dGcc_Circ2d2TanOnIter.cxx:45-480`).

use super::*;

impl Geom2dGccCirc2d2TanOnIter {
    /// `(QualifiedLin, QCurve, gp_Lin2d)` (`cxx:45-151`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_li_cu_on_li(
        q1: &GccEntQualifiedLin,
        q2: &Geom2dGccQCurve,
        on_line: &GpLin2d,
        param1: f64,
        param2: f64,
        param3: f64,
        tolang: f64,
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
        let tol_val = tolang.abs();
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
        // `cxx:90-93`: all four tolerances are `1.e-15`.
        tol.set_value(1, LINEAR_TOL);
        tol.set_value(2, curve_tool::eps_x(cu2, tolang));
        tol.set_value(3, LINEAR_TOL);
        tol.set_value(4, LINEAR_TOL);
        let point1 = line_value_ax2d(param1, l1.position());
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = line_value_ax2d(param3, on_line.position());
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_li_cu_on_li(
            &l1,
            cu2,
            on_line,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1new, _tan1) = line_d1_ax2d(x.value(1), l1.position());
        let (point2new, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let tan3 = GpVec2d::from_dir2d(on_line.direction());
        let point3new = line_value_ax2d(x.value(3), on_line.position());
        let dist1 = point3new.distance(&point1new);
        let dist2 = point3new.distance(&point2new);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3new), (dist1 + dist2) / 2.);
            let angle2 = angle_or_zero(&vec_between(&point2new, &point3new), &tan2);
            let _ = tan3;
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

    /// `(QCurve, QCurve, gp_Lin2d)` (`cxx:153-272`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_cu_on_li(
        q1: &Geom2dGccQCurve,
        q2: &Geom2dGccQCurve,
        on_line: &GpLin2d,
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
        tol.set_value(3, LINEAR_TOL);
        tol.set_value(4, tol_val / 10.);
        let point1 = curve_tool::value(cu1, param1).expect("undefined curve value");
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = line_value_ax2d(param3, on_line.position());
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_cu_cu_on_li(
            cu1,
            cu2,
            on_line,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1new, tan1) = curve_tool::d1(cu1, x.value(1)).expect("undefined curve value");
        let (point2new, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        let tan3 = GpVec2d::from_dir2d(on_line.direction());
        let point3new = line_value_ax2d(x.value(3), on_line.position());
        let _ = tan3;
        let dist1 = point3new.distance(&point1new);
        let dist2 = point3new.distance(&point2new);
        // `cxx:220-223`: the degenerate solution abandons the constructor.
        if (dist1 + dist2) / 2. < tol_val {
            return s;
        }
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3new), (dist1 + dist2) / 2.);
            let angle1 = angle_or_zero(&vec_between(&point1new, &point3new), &tan1);
            let angle2 = angle_or_zero(&vec_between(&point2new, &point3new), &tan2);
            if accept_curve(q1.qualifier(), angle1) && accept_curve(q2.qualifier(), angle2) {
                s.commit(
                    sol,
                    q1.qualifier(),
                    q2.qualifier(),
                    point1new,
                    x.value(1),
                    0.0,
                    point2new,
                    x.value(2),
                    point2new.distance(&point1new),
                    point3new,
                    x.value(3),
                );
            }
        }
        s
    }

    /// `(QCurve, gp_Pnt2d, gp_Lin2d)` (`cxx:274-364`).
    pub fn new_cu_pt_on_li(
        q1: &Geom2dGccQCurve,
        point2: &GpPnt2d,
        on_line: &GpLin2d,
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
        tol.set_value(2, LINEAR_TOL);
        tol.set_value(3, tol_val / 10.);
        let point1 = curve_tool::value(cu1, param1).expect("undefined curve value");
        let point3 = line_value_ax2d(param2, on_line.position());
        ufirst.set_value(3, (point3.distance(point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_cu_pt_on_li(
            cu1,
            point2,
            on_line,
            ufirst.value(3).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1new, tan1) = curve_tool::d1(cu1, x.value(1)).expect("undefined curve value");
        let (point3new, tan3) = line_d1_ax2d(x.value(2), on_line.position());
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
                    x.value(2),
                    par2,
                    point3new,
                    x.value(3),
                );
            }
        }
        s
    }

    /// `(QualifiedCirc, QCurve, gp_Lin2d)` (`cxx:366-480`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_cu_on_li(
        q1: &GccEntQualifiedCirc,
        q2: &Geom2dGccQCurve,
        on_line: &GpLin2d,
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
        tol.set_value(3, LINEAR_TOL);
        tol.set_value(4, tol_val / 10.);
        let point1 = circle_value_ax22d(param1, c1.position(), r1);
        let point2 = curve_tool::value(cu2, param2).expect("undefined curve value");
        let point3 = line_value_ax2d(param3, on_line.position());
        ufirst.set_value(4, (point3.distance(&point2) + point3.distance(&point1)) / 2.);
        let mut func = Geom2dGccFunctionTanCuCuOnCu::new_ci_cu_on_li(
            &c1,
            cu2,
            on_line,
            ufirst.value(4).max(tol_val),
        );
        let x = match Self::solve(&mut func, &tol, &ufirst, &umin, &umax) {
            Some(x) => x,
            None => return s,
        };
        let (point1, tan1, _nor2) = circle_d2_ax22d(x.value(1), c1.position(), r1);
        let (point2, tan2) = curve_tool::d1(cu2, x.value(2)).expect("undefined curve value");
        // `cxx:436`: `ElCLib::Value(Ufirst(1), OnLine)`.
        let point3 = line_value_ax2d(x.value(1), on_line.position());
        let dist1 = point3.distance(&point1);
        let dist2 = point3.distance(&point2);
        if (dist1 - dist2).abs() / 2. <= tol_val {
            let sol = GpCirc2d::new(axis_at(point3), (dist1 + dist2) / 2.);
            let angle2 = angle_or_zero(&vec_between(&point2, &point3), &tan2);
            let _ = tan1;
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
}
