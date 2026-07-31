//! Numerical helpers — root finding, interpolation, geometry numeric utilities.
use crate::gp::GpPnt;

/// Bisection root finding for continuous f on [a,b] where f(a)*f(b) <= 0.
/// Returns root or None if no sign change. max_iter guards convergence.
pub fn bisection<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, max_iter: usize) -> Option<f64> {
    let mut lo = a; let mut hi = b;
    let flo = f(lo); let fhi = f(hi);
    if flo == 0.0 { return Some(lo); }
    if fhi == 0.0 { return Some(hi); }
    if flo.signum() == fhi.signum() { return None; }

    for _ in 0..max_iter {
        let mid = 0.5 * (lo + hi);
        let fm = f(mid);
        if fm == 0.0 || (hi - lo) * 0.5 < tol { return Some(mid); }
        if flo.signum() == fm.signum() { lo = mid; }
        else { hi = mid; }
    }
    Some(0.5 * (lo + hi))
}

/// False position (regula falsi) root finding.
pub fn regula_falsi<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, max_iter: usize) -> Option<f64> {
    let mut lo = a; let mut hi = b;
    let mut flo = f(lo); let mut fhi = f(hi);
    if flo == 0.0 { return Some(lo); }
    if fhi == 0.0 { return Some(hi); }
    if flo.signum() == fhi.signum() { return None; }

    for _ in 0..max_iter {
        let x = (lo * fhi - hi * flo) / (fhi - flo);
        let fx = f(x);
        if fx == 0.0 || (hi - lo).abs() < tol { return Some(x); }
        if flo.signum() == fx.signum() { lo = x; flo = fx; }
        else { hi = x; fhi = fx; }
    }
    Some(0.5 * (lo + hi))
}

/// Secant method root finding (no bracket needed, may diverge).
pub fn secant<F: Fn(f64) -> f64>(f: &F, x0: f64, x1: f64, tol: f64, max_iter: usize) -> Option<f64> {
    let mut xn1 = x0; let mut xn = x1;
    let mut fn1 = f(xn1);
    for _ in 0..max_iter {
        let fn_val = f(xn);
        let denom = fn_val - fn1;
        if denom.abs() < 1e-300 { return None; }
        let x_next = xn - fn_val * (xn - xn1) / denom;
        if (x_next - xn).abs() < tol { return Some(x_next); }
        xn1 = xn; fn1 = fn_val; xn = x_next;
    }
    Some(xn)
}

/// Catmull-Rom spline interpolation through points. Returns y at t in [0,1]
/// within segment i (P0,P1,P2,P3 control points).
pub fn catmull_rom(p0: f64, p1: f64, p2: f64, p3: f64, t: f64) -> f64 {
    0.5 * ((2.0*p1) + (-p0+p2)*t + (2.0*p0-5.0*p1+4.0*p2-p3)*t*t + (-p0+3.0*p1-3.0*p2+p3)*t*t*t)
}

/// Cubic Hermite interpolation between p0,p1 with tangents m0,m1.
pub fn hermite(p0: f64, p1: f64, m0: f64, m1: f64, t: f64) -> f64 {
    let t2 = t*t; let t3 = t2*t;
    let h00 = 2.0*t3 - 3.0*t2 + 1.0;
    let h10 = t3 - 2.0*t2 + t;
    let h01 = -2.0*t3 + 3.0*t2;
    let h11 = t3 - t2;
    h00*p0 + h10*m0 + h01*p1 + h11*m1
}

/// Solve quadratic ax² + bx + c = 0. Returns up to 2 sorted roots.
pub fn quadratic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < 1e-300 {
        if b.abs() < 1e-300 { return vec![]; }
        return vec![-c / b];
    }
    let disc = b*b - 4.0*a*c;
    if disc < 0.0 { return vec![]; }
    let sq = disc.sqrt();
    let q = -0.5 * (b + b.signum() * sq);
    let mut roots = vec![q / a, c / q];
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
    roots
}

/// Solve cubic x³ + a x² + b x + c = 0 (depressed via Cardano). Returns real roots.
pub fn cubic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    let p = b - a*a/3.0;
    let q = 2.0*a*a*a/27.0 - a*b/3.0 + c;
    let disc = q*q/4.0 + p*p*p/27.0;
    let mut roots = Vec::new();
    if disc > 0.0 {
        let sq = disc.sqrt();
        let u = (-q/2.0 + sq).cbrt();
        let v = (-q/2.0 - sq).cbrt();
        roots.push(u + v - a/3.0);
    } else if disc == 0.0 {
        let u = (-q/2.0).cbrt();
        roots.push(2.0*u - a/3.0);
        roots.push(-u - a/3.0);
    } else {
        let r = (-p*p*p/27.0).sqrt();
        let theta = (-q / (2.0*r)).clamp(-1.0, 1.0).acos();
        for k in 0..3 {
            let ang = (theta + 2.0*std::f64::consts::PI*k as f64) / 3.0;
            roots.push(2.0*r.cbrt()*ang.cos() - a/3.0);
        }
    }
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
    roots
}

/// Distance between two points.
pub fn distance(a: &GpPnt, b: &GpPnt) -> f64 { a.coord.subtracted(&b.coord).modulus() }

/// Squared distance between two points (avoids sqrt).
pub fn square_distance(a: &GpPnt, b: &GpPnt) -> f64 { a.coord.subtracted(&b.coord).square_modulus() }

/// Clamp value to [lo, hi].
pub fn clamp(v: f64, lo: f64, hi: f64) -> f64 { v.max(lo).min(hi) }

/// Smoothstep interpolation.
pub fn smoothstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Solve linear system a*x + b*y = e, c*x + d*y = f. Returns (x, y).
pub fn solve_2x2(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> Option<(f64, f64)> {
    let det = a*d - b*c;
    if det.abs() < 1e-300 { return None; }
    let x = (e*d - b*f) / det;
    let y = (a*f - e*c) / det;
    Some((x, y))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bisection_linear() {
        let f = |x: f64| x - 2.0;
        let r = bisection(&f, 0.0, 5.0, 1e-12, 100).unwrap();
        assert!((r - 2.0).abs() < 1e-12);
    }

    #[test]
    fn quadratic() {
        // x² - 3x + 2 = 0 → x=1,2
        let roots = quadratic_roots(1.0, -3.0, 2.0);
        assert_eq!(roots.len(), 2);
        assert!((roots[0] - 1.0).abs() < 1e-12);
        assert!((roots[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn cubic_real_root() {
        // x³ - 6x² + 11x - 6 = 0 → x=1,2,3
        let roots = cubic_roots(-6.0, 11.0, -6.0);
        assert_eq!(roots.len(), 3);
        assert!((roots[0] - 1.0).abs() < 1e-10);
        assert!((roots[1] - 2.0).abs() < 1e-10);
        assert!((roots[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn hermite_matches_endpoints() {
        assert!((hermite(1.0, 3.0, 2.0, 4.0, 0.0) - 1.0).abs() < 1e-12);
        assert!((hermite(1.0, 3.0, 2.0, 4.0, 1.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn solve_2x2_basic() {
        // 2x + y = 5, x + 3y = 6 → x=1.8, y=1.4
        let (x, y) = solve_2x2(2.0, 1.0, 1.0, 3.0, 5.0, 6.0).unwrap();
        assert!((x - 1.8).abs() < 1e-12);
        assert!((y - 1.4).abs() < 1e-12);
    }
}
