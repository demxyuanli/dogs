//! Scalar root-finding and 1-D optimization.
//! Source: `math_BissecNewton`, `math_FunctionRoot`, `math_BrentMinimum`.

/// Bisection root finder on `[a, b]`, which must bracket a root (`f(a)·f(b) < 0`).
pub fn bisection<F>(f: &F, a: f64, b: f64, tol: f64) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let (mut lo, mut hi) = (a, b);
    let mut flo = f(lo);
    let fhi = f(hi);
    if flo * fhi > 0.0 {
        return Err("bisection: root not bracketed".to_string());
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let fm = f(mid);
        if (hi - lo).abs() < tol || fm.abs() < tol {
            return Ok(mid);
        }
        if flo * fm <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = fm;
        }
    }
    Err("bisection: max iterations reached".to_string())
}

/// Secant root finder starting from two initial guesses `x0`, `x1`.
pub fn secant<F>(f: &F, x0: f64, x1: f64, tol: f64) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let (mut x_prev, mut x_cur) = (x0, x1);
    let mut f_prev = f(x_prev);
    for _ in 0..200 {
        let f_cur = f(x_cur);
        if f_cur.abs() < tol {
            return Ok(x_cur);
        }
        let denom = f_cur - f_prev;
        if denom.abs() < 1e-300 {
            return Err("secant: zero slope".to_string());
        }
        let x_next = x_cur - f_cur * (x_cur - x_prev) / denom;
        if (x_next - x_cur).abs() < tol {
            return Ok(x_next);
        }
        x_prev = x_cur;
        x_cur = x_next;
        f_prev = f_cur;
    }
    Err("secant: max iterations reached".to_string())
}

/// Golden-section minimization of `f` on `[a, b]`.
pub fn golden_section<F>(f: &F, a: f64, b: f64, tol: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mut a, mut b) = (a, b);
    let (mut c, mut d) = (b - (b - a) / phi, a + (b - a) / phi);
    let (mut fc, mut fd) = (f(c), f(d));
    while (b - a).abs() > tol {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) / phi;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) / phi;
            fd = f(d);
        }
    }
    0.5 * (a + b)
}

/// Brent's method for a root on `[a, b]`, which must bracket a root.
pub fn brent<F>(f: &F, a: f64, b: f64, tol: f64) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let (mut a, mut b) = (a, b);
    let mut c = b;
    let mut fa = f(a);
    let mut fb = f(b);
    let mut fc = fb;
    if fa * fb > 0.0 {
        return Err("brent: root not bracketed".to_string());
    }
    let eps = f64::EPSILON;
    let mut d = b - a;
    let mut e = d;
    for _ in 0..200 {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = 2.0 * eps * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let p;
            let mut q;
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let rr = fb / fc;
                p = s * (2.0 * xm * qq * (qq - rr) - (b - a) * (rr - 1.0));
                q = (qq - 1.0) * (rr - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            let p = p.abs();
            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min(e.abs() * q) {
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
        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else if xm > 0.0 {
            b += tol1;
        } else {
            b -= tol1;
        }
        fb = f(b);
    }
    Err("brent: max iterations reached".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bisection_square() {
        // x² - 4 = 0 on [0, 5]
        let r = bisection(&|x| x * x - 4.0, 0.0, 5.0, 1e-12).unwrap();
        assert!((r - 2.0).abs() < 1e-10);
    }

    #[test]
    fn secant_cos() {
        // cos(x) - x = 0
        let r = secant(&|x| x.cos() - x, 0.0, 1.0, 1e-12).unwrap();
        assert!((r.cos() - r).abs() < 1e-10);
    }

    #[test]
    fn golden_section_minimum() {
        // (x-2)² + 3 on [0, 5]
        let x = golden_section(&|x| (x - 2.0) * (x - 2.0) + 3.0, 0.0, 5.0, 1e-8);
        assert!((x - 2.0).abs() < 1e-6);
    }

    #[test]
    fn brent_cuberoot() {
        // x³ - 2 = 0
        let r = brent(&|x| x * x * x - 2.0, 0.0, 2.0, 1e-12).unwrap();
        assert!((r - 2.0_f64.cbrt()).abs() < 1e-10);
    }

    #[test]
    fn bisection_unbracketed() {
        assert!(bisection(&|x| x * x + 1.0, 0.0, 1.0, 1e-12).is_err());
    }
}
