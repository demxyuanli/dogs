//! `math_BracketedRoot` Brent root on a bracket.
//! Source: `math_BracketedRoot.cxx:22-121`, `.hxx:42-47`.

use crate::math_fn::MathFunction;

/// Default `NbIterations` (`math_BracketedRoot.hxx:46`).
pub const BRACKETED_ROOT_NB_ITERATIONS: i32 = 100;
/// Default `ZEPS` (`math_BracketedRoot.hxx:47`).
pub const BRACKETED_ROOT_ZEPS: f64 = 1.0e-12;

/// Result of `math_BracketedRoot`.
#[derive(Debug, Clone, Copy)]
pub struct BracketedRoot {
    pub done: bool,
    pub root: f64,
    pub value: f64,
    pub nb_iterations: i32,
}

/// `math_BracketedRoot::math_BracketedRoot` (`cxx:22-121`).
pub fn bracketed_root<F: MathFunction>(
    f: &mut F,
    bound1: f64,
    bound2: f64,
    tolerance: f64,
    nb_iterations: i32,
    zeps: f64,
) -> BracketedRoot {
    let mut fa = 0.0;
    let mut a = bound1;
    let mut c = 0.0;
    let mut d = 0.0;
    let mut e = 0.0;
    let mut root_val = bound2;
    let mut error = 0.0;
    f.value(a, &mut fa);
    f.value(root_val, &mut error);
    if fa * error > 0.0 {
        return BracketedRoot {
            done: false,
            root: root_val,
            value: error,
            nb_iterations: 0,
        };
    }
    let mut fc = error;
    for nb_iter in 1..=nb_iterations {
        if error * fc > 0.0 {
            c = a;
            fc = fa;
            d = root_val - a;
            e = d;
        }
        if fc.abs() < fa.abs() {
            a = root_val;
            root_val = c;
            c = a;
            fa = error;
            error = fc;
            fc = fa;
        }
        let tol1 = 2.0 * zeps * root_val.abs() + 0.5 * tolerance;
        let xm = 0.5 * (c - root_val);
        if xm.abs() <= tol1 || error == 0.0 {
            return BracketedRoot {
                done: true,
                root: root_val,
                value: error,
                nb_iterations: nb_iter,
            };
        }
        if e.abs() >= tol1 && fa.abs() > error.abs() {
            let s = error / fa;
            let (p, mut q) = if a == c {
                (2.0 * xm * s, 1.0 - s)
            } else {
                let qv = fa / fc;
                let r = error / fc;
                (
                    s * (2.0 * xm * qv * (qv - r) - (root_val - a) * (r - 1.0)),
                    (qv - 1.0) * (r - 1.0) * (s - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            }
            let p = p.abs();
            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if 2.0 * p < if min1 < min2 { min1 } else { min2 } {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = root_val;
        fa = error;
        if d.abs() > tol1 {
            root_val += d;
        } else {
            root_val += if xm > 0.0 { tol1.abs() } else { -tol1.abs() };
        }
        f.value(root_val, &mut error);
    }
    BracketedRoot {
        done: false,
        root: root_val,
        value: error,
        nb_iterations,
    }
}
