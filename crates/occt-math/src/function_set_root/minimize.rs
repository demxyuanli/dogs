//! `MinimizeDirection` of `math_FunctionSetRoot.cxx` (`:198-436`).

use super::*;
use super::direction::DirFunction;
use crate::BrentMinimum;

// ---------------------------------------------------------------------------
// MinimizeDirection (math_FunctionSetRoot.cxx:198-436)
// ---------------------------------------------------------------------------

/// `MinimizeDirection(P0, P1, P2, F1, Delta, Tol, F)`
/// (`math_FunctionSetRoot.cxx:198-264`): minimisation from three points.
pub(super) fn minimize_direction_3pt<F: MathFunctionSetWithDerivatives>(
    p0: &MathVector,
    p1: &MathVector,
    p2: &MathVector,
    f1: f64,
    delta: &mut MathVector,
    tol: &MathVector,
    dir: &mut DirFunction,
    f: &mut F,
) -> bool {
    // (1) 1D parametric tolerance.
    let mut tol1d: f64 = 2.1;
    let eps = 1.0e-16;
    for ii in 1..=tol.len() {
        let invnorme = delta.value(ii).abs();
        if invnorme > eps {
            tol1d = tol1d.min(tol.value(ii) / invnorme);
        }
    }
    if tol1d > 1.9 {
        return false; // Pas la peine de se fatiguer
    }
    tol1d /= 3.0;

    *delta = p1.subtracted(p0);
    let mut invnorme = delta.norm();
    if invnorme <= eps {
        return false;
    }
    invnorme = 1.0 / invnorme;

    dir.initialize(p1, delta);

    // (2) On minimise.
    let ax = -1.0;
    let bx = 0.0;
    let cx = p2.subtracted(p1).norm() * invnorme;
    if cx < 1.0e-2 {
        return false;
    }

    let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
    let mut closure = |x: f64| dir.value_scalar(f, x);
    sol.perform_fallible(&mut closure, ax, bx, cx);

    if sol.is_done() {
        let tsol = sol.location().unwrap();
        if sol.minimum().unwrap() < f1 {
            delta.multiply_scalar(tsol);
            return true;
        }
    }
    false
}

/// `MinimizeDirection(P, Dir, PValue, PDirValue, Gradient, DGradient, Tol, F)`
/// (`math_FunctionSetRoot.cxx:267-436`): minimisation from two points and a
/// derivative.
pub(super) fn minimize_direction_2pt<F: MathFunctionSetWithDerivatives>(
    p: &MathVector,
    dir_vec: &mut MathVector,
    p_value: f64,
    p_dir_value: f64,
    gradient: &MathVector,
    d_gradient: &MathVector,
    tol: &MathVector,
    dir: &mut DirFunction,
    f: &mut F,
) -> bool {
    if is_infinite(p_value) || is_infinite(p_dir_value) {
        return false;
    }
    // (0) 1D parametric tolerance.
    let mut good = false;
    let eps = 1.0e-20;
    let mut tol1d: f64 = 1.1;
    let mut result = p_value;
    for ii in 1..=tol.len() {
        let absdir = dir_vec.value(ii).abs();
        if absdir > eps {
            tol1d = tol1d.min(tol.value(ii) / absdir);
        }
    }
    if tol1d > 0.9 {
        return false;
    }

    // (1) Premiere interpolation quadratique.
    let df1 = gradient.dot(dir_vec);
    let df2 = d_gradient.dot(dir_vec);
    let mut tsol;
    if df1 < -eps && df2 > eps {
        // cuvette
        tsol = -df1 / (df2 - df1);
    } else {
        let cx = p_value;
        let bx = df1;
        let ax = p_dir_value - (bx + cx);
        if ax.abs() <= eps {
            // cas lineaire
            if bx.abs() >= eps {
                tsol = -cx / bx;
            } else {
                tsol = 0.0;
            }
        } else {
            // cas quadratique
            let mut disc = bx * bx - 4.0 * ax * cx;
            if disc > 1.0e-9 {
                disc = disc.sqrt();
                tsol = -(bx + disc);
                let tsolbis = disc - bx;
                if tsolbis.abs() < tsol.abs() {
                    tsol = tsolbis;
                }
                tsol /= 2.0 * ax;
            } else {
                // pas ou peu de racine : on "extremise"
                tsol = -(0.5 * bx) / ax;
            }
        }
    }

    if tsol.abs() >= 1.0 {
        return false; // resultat sans interet
    }

    dir.initialize(p, dir_vec);
    // math_FunctionSetRoot.cxx:358 ignores the bool returned by
    // MyDirFunction::Value, leaving fsol indeterminate on failure. We model the
    // failure as "no value": the progress test below is skipped and the Brent
    // search is entered.
    let fsol = dir.value_scalar(f, tsol);
    if let Some(v) = fsol {
        if v < p_value {
            good = true;
            result = v;
        }
    }

    // (2) Recherche en bonne et due forme si l'on n'a pas assez progresse.
    let need_brent = fsol.map_or(true, |v| v > 0.2 * p_value) && tol1d < 0.5;
    if need_brent {
        let (ax, bx, cx);
        if tsol < 0.0 {
            ax = tsol;
            bx = 0.0;
            cx = 1.0;
        } else {
            ax = 0.0;
            bx = tsol;
            cx = 1.0;
        }
        let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
        let mut closure = |x: f64| dir.value_scalar(f, x);
        sol.perform_fallible(&mut closure, ax, bx, cx);
        if sol.is_done() {
            if sol.minimum().unwrap() <= result {
                tsol = sol.location().unwrap();
                good = true;
                result = sol.minimum().unwrap();

                // Objective function changes too fast -> extra computations.
                if gradient.norm2() > 1.0 / SQUARE_CONFUSION && tsol > ax && tsol < cx {
                    sol.perform_fallible(&mut closure, ax, (ax + tsol) / 2.0, tsol);
                    if sol.is_done() && sol.minimum().unwrap() <= result {
                        tsol = sol.location().unwrap();
                        good = true;
                        result = sol.minimum().unwrap();
                    }
                    sol.perform_fallible(&mut closure, tsol, (cx + tsol) / 2.0, cx);
                    if sol.is_done() && sol.minimum().unwrap() <= result {
                        tsol = sol.location().unwrap();
                        good = true;
                        result = sol.minimum().unwrap();
                    }
                }
            }
        }
    }

    if good {
        // mise a jour du Delta
        dir_vec.multiply_scalar(tsol);
    }
    good
}
