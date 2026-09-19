//! `math_DirectPolynomialRoots` - all real roots of a real polynomial of
//! degree <= 4 by direct algebraic methods.
//! Source: `src/FoundationClasses/TKMath/math/math_DirectPolynomialRoots.cxx`
//! and `math_DirectPolynomialRoots.hxx`.
//!
//! Consumers in the 2D intersection chain: `IntAna2d_Outils`'s
//! `MyDirectPolynomialRoots` (`IntAna2d_Outils.cxx:21-215`) and
//! `math_TrigonometricFunctionRoots` (`math_TrigonometricFunctionRoots.cxx:218`
//! and `:367`).

use crate::precision::epsilon;

/// `ZERO_THRESHOLD` (`cxx:44`).
const ZERO_THRESHOLD: f64 = 1.0e-30;
/// `MACHINE_EPSILON = RealEpsilon()` (`cxx:46`).
const MACHINE_EPSILON: f64 = f64::EPSILON;
/// `FLOATING_RADIX` (`cxx:51`).
const FLOATING_RADIX: f64 = 2.0;
/// `INV_LOG_RADIX = 1.0 / std::log(2.0)` (`cxx:54`).
const INV_LOG_RADIX: f64 = 1.4426950408889634; // 1.0 / ln(2)
/// `MAX_NEWTON_ITERATIONS` (`cxx:57`).
const MAX_NEWTON_ITERATIONS: i32 = 10;
/// `OVERFLOW_LIMIT` (`cxx:60`).
const OVERFLOW_LIMIT: f64 = 1.0e+80;

/// `EvaluatePolynomial` (`cxx:66-74`): Horner, `thePoly[0]` is the highest power.
fn evaluate_polynomial(n: i32, poly: &[f64], x: f64) -> f64 {
    let mut result = poly[0];
    for i in 1..n {
        result = result * x + poly[i as usize];
    }
    result
}

/// `EvaluatePolynomialWithDerivative` (`cxx:77-89`).
fn evaluate_polynomial_with_derivative(n: i32, poly: &[f64], x: f64) -> (f64, f64) {
    let mut value = poly[0] * x + poly[1];
    let mut derivative = poly[0];
    for i in 2..n {
        derivative = derivative * x + value;
        value = value * x + poly[i as usize];
    }
    (value, derivative)
}

/// `RefineRoot` (`cxx:92-120`).
fn refine_root(n: i32, poly: &[f64], initial_guess: f64) -> f64 {
    let mut value = 0.0;
    let mut solution = initial_guess;
    let initial_value = evaluate_polynomial(n, poly, initial_guess);

    for _iter in 1..MAX_NEWTON_ITERATIONS {
        let (v, d) = evaluate_polynomial_with_derivative(n, poly, solution);
        value = v;

        if d.abs() <= ZERO_THRESHOLD {
            break;
        }

        let delta = -value / d;

        if delta.abs() <= MACHINE_EPSILON * solution.abs() {
            break;
        }

        solution += delta;
    }

    // Return improved solution only if it is better.
    if value.abs() <= initial_value.abs() {
        solution
    } else {
        initial_guess
    }
}

/// `ScaleAndRefineAllRoots` (`cxx:135-144`).
fn scale_and_refine_all_roots(roots: &mut [f64], nb_roots: i32, scale_factor: f64, coeffs: &[f64]) {
    let n = coeffs.len() as i32;
    for i in 0..nb_roots as usize {
        roots[i] *= scale_factor;
        roots[i] = refine_root(n, coeffs, roots[i]);
    }
}

/// `ComputeBaseExponent` (`cxx:147-158`).
fn compute_base_exponent(value: f64) -> i32 {
    if value > 1.0 {
        (value.ln() * INV_LOG_RADIX) as i32
    } else if value < -1.0 {
        (-((-value).ln()) * INV_LOG_RADIX) as i32
    } else {
        0
    }
}

/// `ScaledCoefficients` (`cxx:162-195`). `e` of `ScaleQuartic` and `d` of
/// `ScaleCubic` are stored but never read afterwards, exactly as in OCCT.
#[allow(dead_code)]
struct ScaledCoefficients {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    scale_factor: f64,
}

impl ScaledCoefficients {
    /// `ScaleQuartic` (`cxx:167-179`).
    fn scale_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        let exp = compute_base_exponent(e) / 4;
        let scale_factor = FLOATING_RADIX.powf(exp as f64);
        let sf2 = scale_factor * scale_factor;
        Self {
            a: a / scale_factor,
            b: b / sf2,
            c: c / (sf2 * scale_factor),
            d: d / (sf2 * sf2),
            e: e / (sf2 * sf2),
            scale_factor,
        }
    }

    /// `ScaleCubic` (`cxx:181-193`).
    fn scale_cubic(a: f64, b: f64, c: f64, d: f64) -> Self {
        let exp = compute_base_exponent(d) / 3;
        let scale_factor = FLOATING_RADIX.powf(exp as f64);
        let sf2 = scale_factor * scale_factor;
        Self {
            a: a / scale_factor,
            b: b / sf2,
            c: c / (sf2 * scale_factor),
            d,
            e: 0.0,
            scale_factor,
        }
    }
}

/// `ComputeSpecialDiscriminant` (`cxx:197-222`).
fn compute_special_discriminant(beta: f64, gamma: f64, del: f64, a1: f64) -> f64 {
    let sigma = beta * gamma / 3.0 - 2.0 * beta * beta * beta / 27.0;
    let psi = gamma * gamma * (4.0 * gamma - beta * beta) / 27.0;

    let d1 = if sigma >= 0.0 {
        sigma + 2.0 * (-a1).sqrt()
    } else {
        sigma - 2.0 * (-a1).sqrt()
    };

    let d2 = psi / d1;

    if (del - d1).abs() >= 18.0 * MACHINE_EPSILON * (del.abs() + d1.abs())
        && (del - d2).abs() >= 24.0 * MACHINE_EPSILON * (del.abs() + d2.abs())
    {
        return (del - d1) * (del - d2) / 4.0;
    }

    0.0
}

/// `SolveCubicThreeRealRoots` (`cxx:225-267`).
fn solve_cubic_three_real_roots(
    beta: f64,
    gamma: f64,
    del: f64,
    p: f64,
    q: f64,
    discr: f64,
    roots: &mut [f64; 4],
) {
    if beta == 0.0 && q == 0.0 {
        // Special case: x^3 + Px = 0
        roots[0] = (-p).sqrt();
        roots[1] = -roots[0];
        roots[2] = 0.0;
    } else {
        let sb = if beta >= 0.0 { 1.0 } else { -1.0 };
        let omega = (0.5 * q / (-discr).sqrt()).atan();
        let sp3 = (-p / 3.0).sqrt();
        let y1 = -2.0 * sb * sp3 * (std::f64::consts::PI / 6.0 - sb * omega / 3.0).cos();

        roots[0] = -beta / 3.0 + y1;

        if beta * q <= 0.0 {
            roots[1] = -beta / 3.0 + 2.0 * sp3 * (omega / 3.0).sin();
        } else {
            // Alternative formula for better accuracy.
            let dbg = del - beta * gamma;
            let sdbg = if dbg >= 0.0 { 1.0 } else { -1.0 };
            let den1 = 8.0 * beta * beta / 9.0 - 4.0 * beta * y1 / 3.0 - 2.0 * q / y1;
            let den2 = 2.0 * y1 * y1 - q / y1;
            roots[1] = dbg / den1 + sdbg * (-27.0 * discr).sqrt() / den2;
        }

        // Use Vieta's formula for the third root.
        roots[2] = -del / (roots[0] * roots[1]);
    }
}

/// `SolveCubicOneRealRoot` (`cxx:270-306`).
fn solve_cubic_one_real_root(
    beta: f64,
    del: f64,
    p: f64,
    q: f64,
    discr: f64,
    roots: &mut [f64; 4],
) {
    let mut u = discr.sqrt() + (q / 2.0).abs();
    u = if u >= 0.0 {
        u.powf(1.0 / 3.0)
    } else {
        -u.abs().powf(1.0 / 3.0)
    };

    let h = if p >= 0.0 {
        u * u + p / 3.0 + (p / u) * (p / u) / 9.0
    } else {
        u * q.abs() / (u * u - p / 3.0)
    };

    if beta * q >= 0.0 {
        if h.abs() <= f64::MIN_POSITIVE && q.abs() <= f64::MIN_POSITIVE {
            roots[0] = -beta / 3.0 - u + p / (3.0 * u);
        } else {
            roots[0] = -beta / 3.0 - q / h;
        }
    } else {
        roots[0] = -del / (beta * beta / 9.0 + h - beta * q / (3.0 * h));
    }
}

/// `SolveCubicMultipleRoots` (`cxx:309-341`).
fn solve_cubic_multiple_roots(
    beta: f64,
    gamma: f64,
    del: f64,
    p: f64,
    q: f64,
    roots: &mut [f64; 4],
    nb_roots: &mut i32,
) {
    *nb_roots = 3;
    let sq = if q >= 0.0 { 1.0 } else { -1.0 };
    let sp3 = (-p / 3.0).sqrt();

    if beta * q <= 0.0 {
        roots[0] = -beta / 3.0 + sq * sp3;
        roots[1] = roots[0];

        if beta * q == 0.0 {
            roots[2] = -beta / 3.0 - 2.0 * sq * sp3;
        } else {
            roots[2] = -del / (roots[0] * roots[1]);
        }
    } else {
        roots[0] = -gamma / (beta + 3.0 * sq * sp3);
        roots[1] = roots[0];
        roots[2] = -beta / 3.0 - 2.0 * sq * sp3;
    }
}

/// `ShouldReduceDegreeQuartic` (`cxx:344-392`).
fn should_reduce_degree_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> bool {
    if a.abs() <= ZERO_THRESHOLD {
        return true;
    }

    // Modified by jgv, 22.01.09 for numerical stability.
    let mut max_coeff = ZERO_THRESHOLD;
    max_coeff = max_coeff.max(b.abs());
    max_coeff = max_coeff.max(c.abs());
    max_coeff = max_coeff.max(d.abs());
    max_coeff = max_coeff.max(e.abs());

    if max_coeff > ZERO_THRESHOLD {
        max_coeff = epsilon(100.0 * max_coeff);
    }

    if a.abs() <= max_coeff {
        let max_coeff1000 = 1000.0 * max_coeff;
        let mut with_a = false;

        if b.abs() > ZERO_THRESHOLD && b.abs() <= max_coeff1000 {
            with_a = true;
        }
        if c.abs() > ZERO_THRESHOLD && c.abs() <= max_coeff1000 {
            with_a = true;
        }
        if d.abs() > ZERO_THRESHOLD && d.abs() <= max_coeff1000 {
            with_a = true;
        }
        if e.abs() > ZERO_THRESHOLD && e.abs() <= max_coeff1000 {
            with_a = true;
        }

        return !with_a;
    }

    false
}

/// `SolveFerrariResolvent` (`cxx:395-428`).
fn solve_ferrari_resolvent(a: f64, b: f64, c: f64, d: f64) -> Option<f64> {
    // Construct resolvent cubic: Y^3 + R3*Y^2 + S3*Y + T3 = 0
    let r3 = -b;
    let s3 = a * c - 4.0 * d;
    let t3 = d * (4.0 * b - a * a) - c * c;

    let cubic_solver = DirectPolynomialRoots::new4(1.0, r3, s3, t3);

    if !cubic_solver.is_done() {
        return None;
    }

    // Choose the largest root for numerical stability.
    let mut y0 = cubic_solver.value(1);
    for i in 2..=cubic_solver.nb_solutions() {
        if cubic_solver.value(i) > y0 {
            y0 = cubic_solver.value(i);
        }
    }

    Some(y0)
}

/// `QuarticFactorization` (`cxx:431-444`).
struct QuarticFactorization {
    p1: f64,
    q1: f64,
    p2: f64,
    q2: f64,
}

/// `FactorQuarticViaFerrari` (`cxx:447-506`).
fn factor_quartic_via_ferrari(a: f64, b: f64, c: f64, d: f64, y0: f64) -> QuarticFactorization {
    // Compute discriminant and parameters.
    let discr = a * y0 * 0.5 - c;
    let sdiscr = if discr >= 0.0 { 1.0 } else { -1.0 };

    // Compute P0 and Q0 for the quadratic factors.
    let mut p0 = a * a * 0.25 - b + y0;
    p0 = if p0 < 0.0 { 0.0 } else { p0.sqrt() };

    let mut q0 = y0 * y0 * 0.25 - d;

    // Handle the case where Q0^2 is very close to zero more robustly.
    if q0.abs() < 10.0 * MACHINE_EPSILON {
        q0 = 0.0;
    } else {
        q0 = if q0 < 0.0 { 0.0 } else { q0.sqrt() };
    }

    let a_demi = a * 0.5;
    let y_demi = y0 * 0.5;
    let sdiscr_q0 = sdiscr * q0;

    let mut p1 = a_demi + p0;
    let mut q1 = y_demi + sdiscr_q0;
    let mut p2 = a_demi - p0;
    let mut q2 = y_demi - sdiscr_q0;

    // Clean up near-zero coefficients.
    let eps = 100.0 * MACHINE_EPSILON;

    if p1.abs() <= eps {
        p1 = 0.0;
    }
    if p2.abs() <= eps {
        p2 = 0.0;
    }
    if q1.abs() <= eps {
        q1 = 0.0;
    }
    if q2.abs() <= eps {
        q2 = 0.0;
    }

    QuarticFactorization { p1, q1, p2, q2 }
}

/// `math_DirectPolynomialRoots` (`math_DirectPolynomialRoots.hxx:41-218`).
#[derive(Debug, Clone, Copy)]
pub struct DirectPolynomialRoots {
    done: bool,
    infinite_status: bool,
    nb_sol: i32,
    roots: [f64; 4],
}

impl DirectPolynomialRoots {
    /// Ctor for `A*x^4 + B*x^3 + C*x^2 + D*x + E = 0` (`cxx:511-521`).
    pub fn new5(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        let mut r = Self {
            infinite_status: false,
            done: true,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve5(a, b, c, d, e);
        r
    }

    /// Ctor for `A*x^3 + B*x^2 + C*x + D = 0` (`cxx:525-533`).
    pub fn new4(a: f64, b: f64, c: f64, d: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve4(a, b, c, d);
        r
    }

    /// Ctor for `A*x^2 + B*x + C = 0` (`cxx:537-545`).
    pub fn new3(a: f64, b: f64, c: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve3(a, b, c);
        r
    }

    /// Ctor for `A*x + B = 0` (`cxx:549-555`).
    pub fn new2(a: f64, b: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve2(a, b);
        r
    }

    /// `IsDone` (`hxx:220-223`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `InfiniteRoots` (`hxx:225-228`).
    pub fn infinite_roots(&self) -> bool {
        self.infinite_status
    }

    /// `NbSolutions` (`hxx:230-234`); `StdFail_InfiniteSolutions` is raised
    /// when there is an infinity of roots.
    pub fn nb_solutions(&self) -> i32 {
        assert!(!self.infinite_status, "StdFail_InfiniteSolutions");
        self.nb_sol
    }

    /// `Value` (`hxx:236-242`); 1-based, with `Standard_RangeError` outside
    /// `[1, NbSolutions]`.
    pub fn value(&self, index: i32) -> f64 {
        assert!(!self.infinite_status, "StdFail_InfiniteSolutions");
        assert!(
            index >= 1 && index <= self.nb_sol,
            "Standard_RangeError in math_DirectPolynomialRoots::Value"
        );
        self.roots[(index - 1) as usize]
    }

    /// `Solve` quartic (`cxx:558-625`).
    fn solve5(&mut self, a: f64, b: f64, c: f64, d: f64, e: f64) {
        // Check for degree reduction.
        if should_reduce_degree_quartic(a, b, c, d, e) {
            self.solve4(b, c, d, e);
            return;
        }

        // Normalize coefficients.
        let na = b / a;
        let nb = c / a;
        let nc = d / a;
        let nd = e / a;

        // Scale coefficients to avoid overflow/underflow.
        // OCCT passes `nd` twice (`cxx:577`); the stored `e` is never read.
        let scaled = ScaledCoefficients::scale_quartic(na, nb, nc, nd, nd);

        // Solve Ferrari's resolvent cubic.
        let y0 = match solve_ferrari_resolvent(scaled.a, scaled.b, scaled.c, scaled.d) {
            Some(y) => y,
            None => {
                self.done = false;
                return;
            }
        };

        // Factor into two quadratics.
        let factors = factor_quartic_via_ferrari(scaled.a, scaled.b, scaled.c, scaled.d, y0);

        // Solve first quadratic.
        let quadratic1 = DirectPolynomialRoots::new3(1.0, factors.p1, factors.q1);
        if !quadratic1.is_done() {
            self.done = false;
            return;
        }

        // Solve second quadratic.
        let quadratic2 = DirectPolynomialRoots::new3(1.0, factors.p2, factors.q2);
        if !quadratic2.is_done() {
            self.done = false;
            return;
        }

        // Collect all roots.
        self.nb_sol = quadratic1.nb_sol + quadratic2.nb_sol;
        let mut index = 0usize;
        for i in 0..quadratic1.nb_sol as usize {
            self.roots[index] = quadratic1.roots[i];
            index += 1;
        }
        for i in 0..quadratic2.nb_sol as usize {
            self.roots[index] = quadratic2.roots[i];
            index += 1;
        }

        // Apply inverse scaling and Newton-Raphson refinement to all roots.
        let coeffs = [a, b, c, d, e];
        scale_and_refine_all_roots(
            &mut self.roots,
            self.nb_sol,
            scaled.scale_factor,
            &coeffs,
        );
    }

    /// `Solve` cubic (`cxx:627-708`).
    fn solve4(&mut self, a: f64, b: f64, c: f64, d: f64) {
        // Check for degree reduction.
        if a.abs() <= ZERO_THRESHOLD {
            self.solve3(b, c, d);
            return;
        }

        // Normalize coefficients.
        let beta = b / a;
        let gamma = c / a;
        let del = d / a;

        // Scale to avoid overflow/underflow. OCCT passes `del` twice
        // (`cxx:648`); the stored `d` is never read.
        let scaled = ScaledCoefficients::scale_cubic(beta, gamma, del, del);

        // Transform to depressed cubic: t^3 + Pt + Q = 0
        let p1 = scaled.b;
        let p2 = -(scaled.a * scaled.a) / 3.0;
        let mut p = p1 + p2;
        let ep = 5.0 * MACHINE_EPSILON * (p1.abs() + p2.abs());
        if p.abs() <= ep {
            p = 0.0;
        }

        let q1 = scaled.c;
        let q2 = -scaled.a * scaled.b / 3.0;
        let q3 = 2.0 * (scaled.a * scaled.a * scaled.a) / 27.0;
        let mut q = q1 + q2 + q3;
        let eq = 10.0 * MACHINE_EPSILON * (q1.abs() + q2.abs() + q3.abs());
        if q.abs() <= eq {
            q = 0.0;
        }

        // Check for overflow.
        if p.abs() > OVERFLOW_LIMIT {
            self.done = false;
            return;
        }

        // Compute discriminant.
        let a1 = (p * p * p) / 27.0;
        let a2 = (q * q) / 4.0;
        let mut discr = a1 + a2;

        // Special handling for P < 0.
        if p < 0.0 {
            discr = compute_special_discriminant(scaled.a, scaled.b, scaled.c, a1);
        }

        // Solve based on discriminant.
        if discr < 0.0 {
            // Three distinct real roots.
            self.nb_sol = 3;
            solve_cubic_three_real_roots(
                scaled.a,
                scaled.b,
                scaled.c,
                p,
                q,
                discr,
                &mut self.roots,
            );
        } else if discr > 0.0 {
            // One real root.
            self.nb_sol = 1;
            solve_cubic_one_real_root(scaled.a, scaled.c, p, q, discr, &mut self.roots);
        } else {
            // Multiple roots.
            let mut nb = 0;
            solve_cubic_multiple_roots(
                scaled.a,
                scaled.b,
                scaled.c,
                p,
                q,
                &mut self.roots,
                &mut nb,
            );
            self.nb_sol = nb;
        }

        // Apply inverse scaling and Newton-Raphson refinement to all roots.
        let coeffs = [a, b, c, d];
        scale_and_refine_all_roots(
            &mut self.roots,
            self.nb_sol,
            scaled.scale_factor,
            &coeffs,
        );
    }

    /// `Solve` quadratic (`cxx:712-764`).
    fn solve3(&mut self, a: f64, b: f64, c: f64) {
        // Check for degree reduction.
        if a.abs() <= ZERO_THRESHOLD {
            self.solve2(b, c);
            return;
        }

        // Solve normalized quadratic x^2 + P*x + Q = 0.
        let p = b / a;
        let q = c / a;

        // Compute discriminant with error bounds.
        let eps_d = 3.0 * MACHINE_EPSILON * (p * p + (4.0 * q).abs());
        let mut discrim = p * p - 4.0 * q;

        if discrim.abs() <= eps_d {
            discrim = 0.0;
        }

        if discrim < 0.0 {
            // No real roots.
            self.nb_sol = 0;
        } else if discrim == 0.0 {
            // Double root.
            self.nb_sol = 2;
            self.roots[0] = -0.5 * p;
            self.roots[0] = refine_root(3, &[1.0, p, q], self.roots[0]);
            self.roots[1] = self.roots[0];
        } else {
            // Two distinct real roots - use numerically stable formula.
            self.nb_sol = 2;
            if p > 0.0 {
                self.roots[0] = -(p + discrim.sqrt()) / 2.0;
            } else {
                self.roots[0] = -(p - discrim.sqrt()) / 2.0;
            }
            self.roots[0] = refine_root(3, &[1.0, p, q], self.roots[0]);
            self.roots[1] = q / self.roots[0];
            self.roots[1] = refine_root(3, &[1.0, p, q], self.roots[1]);
        }
    }

    /// `Solve` linear (`cxx:768-786`).
    fn solve2(&mut self, a: f64, b: f64) {
        if a.abs() <= ZERO_THRESHOLD {
            if b.abs() <= ZERO_THRESHOLD {
                // 0 = 0: infinite solutions.
                self.infinite_status = true;
                return;
            }
            // 0*x + B = 0: no solution.
            self.nb_sol = 0;
            return;
        }

        // Normal case: unique solution.
        self.nb_sol = 1;
        self.roots[0] = -b / a;
    }
}
