//! `AppParCurves_Function::Value` (`AppParCurves_Function.gxx:189-293`).
//!
//! Translation notes:
//! - The unconstrained arm delegates to
//!   `AppParCurves_LeastSquare::Error` (`gxx:202-206`), which already reports
//!   the sum of squares and the largest 2d / 3d distances.
//! - The constrained arm first refreshes the Bezier poles, then hands the
//!   curve to `ResolConstraint` (which moves the constrained poles onto the
//!   points) and recomputes `F` from the moved poles (`gxx:214-292`).
//! - `myMultiCurve` is replaced from `BezierValue()` (`gxx:219`), so the
//!   poles `ResolConstraint` writes are the ones measured below.

use occt_math::MathVector;

use crate::app_par_curves::resol_constraint::ResolConstraint;

use super::ParFunction;
use crate::approx::line_tool::LineTool;

/// `Value(X, F)` (`gxx:189-293`).
pub(super) fn value<L: LineTool>(this: &mut ParFunction<L>, x: &MathVector, f: &mut f64) -> bool {
    this.my_parameters = x.clone();

    // Resolution moindres carres (`gxx:194-201`).
    this.my_least_square.perform_parameters(&this.my_parameters);
    if !this.my_least_square.is_done() {
        this.done = false;
        return false;
    }

    if !this.contraintes {
        // `MyLeastSquare.Error(FVal, ERR3d, ERR2d); F = FVal;` (`gxx:202-206`).
        let (fval, err3d, err2d) = this.my_least_square.error();
        this.f_val = fval;
        this.err3d = err3d;
        this.err2d = err2d;
        *f = fval;
        return true;
    }

    // Resolution avec contraintes (`gxx:210-292`).
    let npol = this.degre + 1;
    this.err3d = 0.0;
    this.err2d = 0.0;

    this.my_multi_curve = this.my_least_square.bezier_value();
    this.a = this.my_least_square.function_matrix().clone();

    let resol = ResolConstraint::new(
        &this.my_multi_line,
        &mut this.my_multi_curve,
        this.first_p,
        this.last_p,
        &this.my_constraints,
        &this.a,
        this.my_least_square.derivative_function_matrix(),
        1.0e-10,
    );
    if !resol.is_done() {
        this.done = false;
        return false;
    }

    // Calcul de F = Sum||C(ui)-Ptli||2 sur toutes les courbes (`gxx:237-291`).
    let mut fval = 0.0;
    let mut ptcxci = MathVector::new(1, npol.max(1) as usize);
    let mut ptcyci = MathVector::new(1, npol.max(1) as usize);
    let mut ptczci = MathVector::new(1, npol.max(1) as usize);

    for ci in 1..=this.nb_cu {
        let dimen = this.tabdim[(ci - 1) as usize];
        for j in 1..=npol {
            let mp = this.my_multi_curve.value(j);
            if dimen == 3 {
                let p = mp.point(ci);
                ptcxci.set_value(j as usize, p.x());
                ptcyci.set_value(j as usize, p.y());
                ptczci.set_value(j as usize, p.z());
            } else {
                let p = mp.point2d(ci);
                ptcxci.set_value(j as usize, p.x());
                ptcyci.set_value(j as usize, p.y());
            }
        }

        for i in this.adeb..=this.afin {
            let mut aa = 0.0;
            let mut bb = 0.0;
            let mut cc = 0.0;
            for j in 1..=npol {
                let aij = this.a.value(i as usize, j as usize);
                aa += aij * ptcxci.value(j as usize);
                bb += aij * ptcyci.value(j as usize);
                if dimen == 3 {
                    cc += aij * ptczci.value(j as usize);
                }
            }
            let fx = aa - this.ptlx.value(i as usize, ci as usize);
            let fy = bb - this.ptly.value(i as usize, ci as usize);
            let mut myf = fx * fx + fy * fy;
            if dimen == 3 {
                let fz = cc - this.ptlz.value(i as usize, ci as usize);
                myf += fz * fz;
                this.my_f.set_value(i as usize, ci as usize, myf);
                let fi = myf;
                if fi.sqrt() > this.err3d {
                    this.err3d = fi.sqrt();
                }
            } else {
                this.my_f.set_value(i as usize, ci as usize, myf);
                let fi = myf;
                if fi.sqrt() > this.err2d {
                    this.err2d = fi.sqrt();
                }
            }
            fval += myf;
        }
    }

    this.f_val = fval;
    *f = fval;
    true
}


