//! Fletcher-Reeves / Polak-Ribiere conjugate-gradient minimization.
//! Source: `math_FRPR.cxx` / `math_FRPR.hxx`.
//!
//! Implements the FRPR conjugate-gradient algorithm for unconstrained
//! minimization of a function of several variables. The gradient is provided
//! by the caller; a 1D bracket + Brent line search is used along each
//! conjugate direction (mirroring `math_BracketMinimum` / `math_BrentMinimum`).

use crate::{MathVector, MathStatus};

/// Result of a successful FRPR run.
#[derive(Debug, Clone)]
pub struct FrprResult {
    /// Minimizing point.
    pub location: MathVector,
    /// Objective value at `location`.
    pub minimum: f64,
    /// Number of outer (conjugate-direction) iterations.
    pub nb_iter: usize,
    /// Gradient at `location`.
    pub gradient: MathVector,
    /// Always `true` for a `Ok` result; kept for API symmetry with OCCT.
    pub done: bool,
}

/// FRPR conjugate-gradient optimizer.
///
/// Defaults follow the OCCT constructor: tolerance 1e-6, at most 200
/// iterations, ZEPS 1e-12.
pub struct Frpr {
    /// Convergence tolerance on the relative change of the objective.
    pub tolerance: f64,
    /// Maximum number of conjugate-direction iterations.
    pub max_iter: usize,
    /// Small absolute guard used inside the convergence test.
    pub zeps: f64,
}

impl Default for Frpr {
    fn default() -> Self {
        Self::new()
    }
}

impl Frpr {
    /// Create with OCCT default parameters.
    pub fn new() -> Self {
        Self { tolerance: 1e-6, max_iter: 200, zeps: 1e-12 }
    }

    /// Minimize `f` starting from `start`; `grad` supplies the gradient of `f`.
    ///
    /// The solution is found when
    /// `2*|F(i) - F(i-1)| <= tolerance*(|F(i)| + |F(i-1)| + zeps)`.
    pub fn perform<F, G>(&self, f: F, grad: G, start: &MathVector) -> Result<FrprResult, MathStatus>
    where
        F: Fn(&MathVector) -> f64,
        G: Fn(&MathVector) -> MathVector,
    {
        let n = start.len();
        if n == 0 {
            return Err(MathStatus::FunctionError);
        }
        let mut location = start.clone();
        let mut the_grad = grad(&location);
        if has_non_finite(&the_grad) {
            return Err(MathStatus::FunctionError);
        }
        let mut prev_min = f(&location);
        if !prev_min.is_finite() {
            return Err(MathStatus::FunctionError);
        }

        let mut g = the_grad.opposite();
        let mut h = g.clone();
        the_grad = g.clone(); // first search direction = -gradient

        for its in 1..=self.max_iter {
            // Minimize f along direction `the_grad` from `location`.
            let Some((scale, the_min)) = minimize_direction(&f, &location, &the_grad) else {
                return Err(MathStatus::DirectionSearchError);
            };
            if !the_min.is_finite() {
                return Err(MathStatus::FunctionError);
            }
            // Apply the step: location += scale * direction.
            for i in 1..=n {
                the_grad.set_value(i, the_grad.value(i) * scale);
            }
            for i in 1..=n {
                location.set_value(i, location.value(i) + the_grad.value(i));
            }
            if 2.0 * (the_min - prev_min).abs()
                <= self.tolerance * (the_min.abs() + prev_min.abs() + self.zeps)
            {
                let gradient = grad(&location);
                return Ok(FrprResult {
                    location,
                    minimum: the_min,
                    nb_iter: its,
                    gradient,
                    done: true,
                });
            }
            // Re-evaluate objective and gradient at the new location.
            prev_min = f(&location);
            the_grad = grad(&location);
            if !prev_min.is_finite() || has_non_finite(&the_grad) {
                return Err(MathStatus::FunctionError);
            }

            // Polak-Ribiere conjugate gradient update.
            let mut dgg = 0.0;
            let mut gg = 0.0;
            for j in 1..=n {
                gg += g.value(j) * g.value(j);
                dgg += (the_grad.value(j) + g.value(j)) * the_grad.value(j);
            }
            if gg == 0.0 {
                // Exact stationary point reached.
                return Err(MathStatus::FunctionError);
            }
            let gam = dgg / gg;
            g = the_grad.opposite();
            for j in 1..=n {
                the_grad.set_value(j, g.value(j) + gam * h.value(j));
            }
            h = the_grad.clone();
        }
        Err(MathStatus::TooManyIterations)
    }
}

/// 1D line search along `dir` from `p0`: returns `(step, f(p0 + step*dir))`.
fn minimize_direction<F: Fn(&MathVector) -> f64>(
    f: &F,
    p0: &MathVector,
    dir: &MathVector,
) -> Option<(f64, f64)> {
    let n = p0.len();
    let dirf = |t: f64| -> f64 {
        let mut p = MathVector::new(1, n);
        for i in 1..=n {
            p.set_value(i, p0.value(i) + t * dir.value(i));
        }
        f(&p)
    };
    let (ax, bx, cx) = bracket_minimum(&dirf, 0.0, 1.0)?;
    brent_minimum(&dirf, ax, bx, cx, 1e-10, 1e-12, 200)
}

/// Bracket a local minimum of `f` starting from points `a0`, `b0`
/// (Numerical Recipes `mnbrak`, as in `math_BracketMinimum`).
fn bracket_minimum<F: Fn(f64) -> f64>(f: &F, a0: f64, b0: f64) -> Option<(f64, f64, f64)> {
    const GOLD: f64 = 1.618_034;
    const GLIMIT: f64 = 100.0;
    const TINY: f64 = 1.0e-20;

    let mut ax = a0;
    let mut bx = b0;
    let mut fa = f(ax);
    let mut fb = f(bx);
    if fb > fa {
        std::mem::swap(&mut ax, &mut bx);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut cx = bx + GOLD * (bx - ax);
    let mut fc = f(cx);
    while fb > fc {
        let r = (bx - ax) * (fb - fc);
        let q = (bx - cx) * (fb - fa);
        let mut u = bx - ((bx - cx) * q - (bx - ax) * r)
            / (2.0 * sign((q - r).abs().max(TINY), q - r));
        let ulim = bx + GLIMIT * (cx - bx);
        let fu;
        if (bx - u) * (u - cx) > 0.0 {
            // u lies between b and c.
            let fu0 = f(u);
            if fu0 < fc {
                return Some((bx, u, cx));
            }
            if fu0 > fb {
                return Some((ax, bx, u));
            }
            // Take the next probe after (b, c).
            u = cx + GOLD * (cx - bx);
            fu = f(u);
        } else if (cx - u) * (u - ulim) > 0.0 {
            // u beyond c but before the limit.
            fu = f(u);
        } else if (u - ulim) * (ulim - cx) >= 0.0 {
            // u at or beyond the limit.
            u = ulim;
            fu = f(u);
        } else {
            // Reset u to the next probe after (b, c).
            u = cx + GOLD * (cx - bx);
            fu = f(u);
        }
        ax = bx;
        bx = cx;
        cx = u;
        fa = fb;
        fb = fc;
        fc = fu;
    }
    Some((ax, bx, cx))
}

/// Brent 1D minimization on bracket `(ax, bx, cx)` (as `math_BrentMinimum`).
fn brent_minimum<F: Fn(f64) -> f64>(
    f: &F,
    ax: f64,
    bx: f64,
    cx: f64,
    xtol: f64,
    zeps: f64,
    itermax: usize,
) -> Option<(f64, f64)> {
    const CGOLD: f64 = 0.381_966;

    let mut a = ax.min(cx);
    let mut b = ax.max(cx);
    let mut x = bx;
    let mut w = bx;
    let mut v = bx;
    let mut fx = f(x);
    let mut fw = fx;
    let mut fv = fx;
    let mut e: f64 = 0.0;
    let mut d: f64 = f64::MAX;
    for _ in 0..itermax {
        let xm = 0.5 * (a + b);
        let tol1 = xtol * x.abs() + zeps;
        let tol2 = 2.0 * tol1;
        if x <= tol2 + a && x >= b - tol2 {
            return Some((x, fx));
        }
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            let mut q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let etemp = e;
            e = d;
            if p.abs() >= (0.5 * q * etemp).abs() || p <= q * (a - x) || p >= q * (b - x) {
                e = if x >= xm { a - x } else { b - x };
                d = CGOLD * e;
            } else {
                d = p / q;
                let u = x + d;
                if u - a < tol2 || b - u < tol2 {
                    d = sign(tol1, xm - x);
                }
            }
        } else {
            e = if x >= xm { a - x } else { b - x };
            d = CGOLD * e;
        }
        let u = if d.abs() >= tol1 { x + d } else { x + sign(tol1, d) };
        let fu = f(u);
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    None
}

/// `sign(a, b)`: magnitude of `a`, sign of `b` (as the C `SIGN` macro).
fn sign(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

fn has_non_finite(v: &MathVector) -> bool {
    v.as_slice().iter().any(|x| !x.is_finite())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BFGS;

    fn rosenbrock(x: &MathVector) -> f64 {
        (1.0 - x.value(1)).powi(2) + 100.0 * (x.value(2) - x.value(1) * x.value(1)).powi(2)
    }

    fn rosenbrock_grad(x: &MathVector) -> MathVector {
        let x1 = x.value(1);
        let x2 = x.value(2);
        MathVector::from_slice(&[
            2.0 * (x1 - 1.0) - 400.0 * x1 * (x2 - x1 * x1),
            200.0 * (x2 - x1 * x1),
        ])
    }

    #[test]
    fn frpr_rosenbrock() {
        let frpr = Frpr::new();
        let start = MathVector::from_slice(&[-1.0, 1.0]);
        let r = frpr.perform(rosenbrock, rosenbrock_grad, &start).expect("FRPR converges");
        assert!(r.minimum < 1e-4, "min={}", r.minimum);
        assert!((r.location.value(1) - 1.0).abs() < 1e-2, "x={}", r.location.value(1));
        assert!((r.location.value(2) - 1.0).abs() < 1e-2, "y={}", r.location.value(2));
        assert!(r.done);
        assert!(r.nb_iter > 0 && r.nb_iter < 200);
    }

    #[test]
    fn frpr_from_origin() {
        // Starting from (0,0) the Rosenbrock valley is harder; allow more iterations.
        let frpr = Frpr { tolerance: 1e-6, max_iter: 500, zeps: 1e-12 };
        let start = MathVector::from_slice(&[0.0, 0.0]);
        let r = frpr.perform(rosenbrock, rosenbrock_grad, &start).expect("FRPR converges");
        assert!(r.minimum < 1e-3, "min={}", r.minimum);
        assert!((r.location.value(1) - 1.0).abs() < 5e-2);
        assert!((r.location.value(2) - 1.0).abs() < 5e-2);
    }

    #[test]
    fn frpr_matches_bfgs() {
        let frpr = Frpr::new();
        let bfgs = BFGS::new();
        let start = MathVector::from_slice(&[-1.0, 1.0]);
        let fr = frpr.perform(rosenbrock, rosenbrock_grad, &start).unwrap();
        let bf = bfgs.minimize(rosenbrock, rosenbrock_grad, &start).unwrap();
        assert!((fr.minimum - rosenbrock(&bf)).abs() < 1e-4);
    }

    #[test]
    fn frpr_quadratic() {
        let f = |x: &MathVector| (x.value(1) - 3.0).powi(2) + (x.value(2) + 2.0).powi(2);
        let g = |x: &MathVector| {
            MathVector::from_slice(&[2.0 * (x.value(1) - 3.0), 2.0 * (x.value(2) + 2.0)])
        };
        let frpr = Frpr::new();
        let r = frpr.perform(f, g, &MathVector::from_slice(&[0.0, 0.0])).unwrap();
        assert!((r.location.value(1) - 3.0).abs() < 1e-4);
        assert!((r.location.value(2) + 2.0).abs() < 1e-4);
        assert!(r.minimum < 1e-8);
    }

    #[test]
    fn frpr_too_many_iterations() {
        let f = |x: &MathVector| (x.value(1) - 3.0).powi(2);
        let g = |x: &MathVector| MathVector::from_slice(&[2.0 * (x.value(1) - 3.0)]);
        let frpr = Frpr { tolerance: 1e-12, max_iter: 1, zeps: 1e-12 };
        let err = frpr.perform(f, g, &MathVector::from_slice(&[0.0])).unwrap_err();
        assert_eq!(err, MathStatus::TooManyIterations);
    }

    #[test]
    fn frpr_nan_returns_error() {
        let f = |_: &MathVector| f64::NAN;
        let g = |_: &MathVector| MathVector::from_slice(&[0.0]);
        let frpr = Frpr::new();
        assert_eq!(frpr.perform(f, g, &MathVector::from_slice(&[1.0])).unwrap_err(),
                   MathStatus::FunctionError);
    }
}
