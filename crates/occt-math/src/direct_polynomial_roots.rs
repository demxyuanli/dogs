//! `math_DirectPolynomialRoots` — all real roots of a polynomial of degree
//! `<= 4` by direct algebraic methods.
//!
//! Source: `FoundationClasses/TKMath/math/math_DirectPolynomialRoots.cxx:1-811`
//! and `.hxx:29-226` (OCCT 8.0.0). The class body was rewritten for 8.0.0
//! (Ferrari + Cardano + coefficient scaling + Newton polish); the port follows
//! that file line by line:
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | anonymous-namespace helpers | `cxx:41-508` | module-private `fn`s |
//! | `Solve(A,B,C,D,E)` quartic | `cxx:556-626` | [`DirectPolynomialRoots::solve_quartic`] |
//! | `Solve(A,B,C,D)` cubic | `cxx:628-709` | [`DirectPolynomialRoots::solve_cubic`] |
//! | `Solve(A,B,C)` quadratic | `cxx:711-765` | [`DirectPolynomialRoots::solve_quadratic`] |
//! | `Solve(A,B)` linear | `cxx:767-786` | [`DirectPolynomialRoots::solve_linear`] |
//! | `IsDone` / `InfiniteRoots` / `NbSolutions` / `Value` | `hxx:198-220` | same names, snake_case |
//!
//! Translation notes:
//! - `RealEpsilon()` is [`occt_core::precision::REAL_EPSILON`], `RealSmall()` is
//!   [`occt_core::precision::REAL_SMALL`], and `Epsilon(v)` (`Standard_Real.hxx:242-246`,
//!   the ULP of `v`) is [`occt_core::precision::epsilon`],
//! - the four root slots stay a fixed `[f64; 4]`: `myRoots` is only ever read
//!   through `Value` for indices `<= myNbSol` (`hxx:214-219`),
//! - `NbSolutions` / `Value` raise `StdFail_InfiniteSolutions` and
//!   `Standard_RangeError` (`hxx:203-219`); the port returns `0` / `f64::NAN`
//!   instead of throwing. No caller in the ported tree indexes past
//!   `nb_solutions()`,
//! - `Dump` (`cxx:790-811`) only formats diagnostics to a stream; it is
//!   `UNPORTED` because the port has no `Standard_OStream` equivalent.

use occt_core::precision::{epsilon, REAL_EPSILON, REAL_SMALL};

/// Threshold below which coefficients are considered effectively zero
/// (`cxx:45`).
const ZERO_THRESHOLD: f64 = 1.0e-30;

/// Machine epsilon for floating-point arithmetic precision (`cxx:48`).
const MACHINE_EPSILON: f64 = REAL_EPSILON;

/// Floating-point radix (base) for coefficient scaling (`cxx:51`).
const FLOATING_RADIX: f64 = 2.0;

/// Inverse of natural logarithm of radix for exponent calculations (`cxx:54`).
const INV_LOG_RADIX: f64 = 1.0 / std::f64::consts::LN_2;

/// Maximum number of Newton-Raphson refinement iterations (`cxx:57`).
const MAX_NEWTON_ITERATIONS: i32 = 10;

/// Upper limit to prevent numerical overflow in computations (`cxx:60`).
const OVERFLOW_LIMIT: f64 = 1.0e+80;

/// `EvaluatePolynomial` (`cxx:66-73`): Horner over `poly[0..n]`.
fn evaluate_polynomial(n: usize, poly: &[f64], x: f64) -> f64 {
    let mut result = poly[0];
    for i in 1..n {
        result = result * x + poly[i];
    }
    result
}

/// `EvaluatePolynomialWithDerivative` (`cxx:76-89`).
fn evaluate_polynomial_with_derivative(
    n: usize,
    poly: &[f64],
    x: f64,
    value: &mut f64,
    derivative: &mut f64,
) {
    *value = poly[0] * x + poly[1];
    *derivative = poly[0];

    for i in 2..n {
        *derivative = *derivative * x + *value;
        *value = *value * x + poly[i];
    }
}

/// `RefineRoot` (`cxx:92-120`): Newton-Raphson polish, kept only when the
/// residual does not grow.
fn refine_root(n: usize, poly: &[f64], initial_guess: f64) -> f64 {
    let mut value = 0.0;
    let mut derivative = 0.0;
    let mut solution = initial_guess;
    let initial_value = evaluate_polynomial(n, poly, initial_guess);

    for _iter in 1..MAX_NEWTON_ITERATIONS {
        evaluate_polynomial_with_derivative(n, poly, solution, &mut value, &mut derivative);

        if derivative.abs() <= ZERO_THRESHOLD {
            break;
        }

        let delta = -value / derivative;

        if delta.abs() <= MACHINE_EPSILON * solution.abs() {
            break;
        }

        solution += delta;
    }

    // `cxx:119`: return improved solution only if it's better.
    if value.abs() <= initial_value.abs() {
        solution
    } else {
        initial_guess
    }
}

/// `ComputeBaseExponent` (`cxx:147-158`).
fn compute_base_exponent(value: f64) -> i32 {
    if value > 1.0 {
        (value.ln() * INV_LOG_RADIX) as i32
    } else if value < -1.0 {
        (-(-value).ln() * INV_LOG_RADIX) as i32
    } else {
        0
    }
}

/// `ScaledCoefficients` (`cxx:163-196`).
struct ScaledCoefficients {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    scale_factor: f64,
}

impl ScaledCoefficients {
    /// `ScaleQuartic` (`cxx:168-181`).
    ///
    /// `cxx:578` calls it as `ScaleQuartic(A, B, C, D, D)`, so the fifth
    /// coefficient the body divides is the fourth: the `E` member the struct
    /// also holds is never read back by the quartic path. The port keeps the
    /// same argument list and drops the unread member.
    fn scale_quartic(a: f64, b: f64, c: f64, d: f64) -> Self {
        let exp = compute_base_exponent(d) / 4;
        let scale_factor = FLOATING_RADIX.powf(exp as f64);
        let s2 = scale_factor * scale_factor;
        Self {
            a: a / scale_factor,
            b: b / s2,
            c: c / (s2 * scale_factor),
            d: d / (s2 * s2),
            scale_factor,
        }
    }

    /// `ScaleCubic` (`cxx:183-195`). `cxx:648` passes `theDel` twice; the
    /// fourth argument is stored but never read.
    fn scale_cubic(a: f64, b: f64, c: f64, d: f64) -> Self {
        let exp = compute_base_exponent(d) / 3;
        let scale_factor = FLOATING_RADIX.powf(exp as f64);
        let s2 = scale_factor * scale_factor;
        Self {
            a: a / scale_factor,
            b: b / s2,
            c: c / (s2 * scale_factor),
            d,
            scale_factor,
        }
    }
}

/// `ComputeSpecialDiscriminant` (`cxx:199-222`).
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

/// `SolveCubicThreeRealRoots` (`cxx:225-268`).
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
        // `cxx:235-240`: x^3 + Px = 0.
        roots[0] = (-p).sqrt();
        roots[1] = -roots[0];
        roots[2] = 0.0;
        return;
    }

    let sb = if beta >= 0.0 { 1.0 } else { -1.0 };
    let omega = (0.5 * q / (-discr).sqrt()).atan();
    let sp3 = (-p / 3.0).sqrt();
    let y1 = -2.0 * sb * sp3 * (std::f64::consts::PI / 6.0 - sb * omega / 3.0).cos();

    roots[0] = -beta / 3.0 + y1;

    // `cxx:250-265`.
    if beta * q <= 0.0 {
        roots[1] = -beta / 3.0 + 2.0 * sp3 * (omega / 3.0).sin();
    } else {
        let dbg = del - beta * gamma;
        let sdbg = if dbg >= 0.0 { 1.0 } else { -1.0 };
        let den1 = 8.0 * beta * beta / 9.0 - 4.0 * beta * y1 / 3.0 - 2.0 * q / y1;
        let den2 = 2.0 * y1 * y1 - q / y1;
        roots[1] = dbg / den1 + sdbg * (-27.0 * discr).sqrt() / den2;
    }

    // `cxx:267`: Vieta's formula for the third root.
    roots[2] = -del / (roots[0] * roots[1]);
}

/// `SolveCubicOneRealRoot` (`cxx:271-306`).
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
        -(u.abs().powf(1.0 / 3.0))
    };

    let h = if p >= 0.0 {
        u * u + p / 3.0 + (p / u) * (p / u) / 9.0
    } else {
        u * q.abs() / (u * u - p / 3.0)
    };

    // `cxx:292-305`.
    if beta * q >= 0.0 {
        if h.abs() <= REAL_SMALL && q.abs() <= REAL_SMALL {
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
    nb_roots: &mut usize,
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

/// `ShouldReduceDegreeQuartic` (`cxx:344-393`).
fn should_reduce_degree_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> bool {
    if a.abs() <= ZERO_THRESHOLD {
        return true;
    }

    // `cxx:356-362`: max of the four lower coefficients.
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

/// `SolveFerrariResolvent` (`cxx:396-428`): the largest root of
/// `Y^3 + R3*Y^2 + S3*Y + T3 = 0`, `None` for `!IsDone()`.
fn solve_ferrari_resolvent(a: f64, b: f64, c: f64, d: f64) -> Option<f64> {
    let r3 = -b;
    let s3 = a * c - 4.0 * d;
    let t3 = d * (4.0 * b - a * a) - c * c;

    let cubic = DirectPolynomialRoots::cubic(1.0, r3, s3, t3);

    if !cubic.is_done() {
        return None;
    }

    // `cxx:418-426`: choose the largest root for numerical stability.
    let mut y0 = cubic.value(1);
    for i in 2..=cubic.nb_solutions() {
        if cubic.value(i) > y0 {
            y0 = cubic.value(i);
        }
    }

    Some(y0)
}

/// `QuarticFactorization` (`cxx:430-445`).
struct QuarticFactorization {
    p1: f64,
    q1: f64,
    p2: f64,
    q2: f64,
}

/// `FactorQuarticViaFerrari` (`cxx:448-506`).
fn factor_quartic_via_ferrari(a: f64, b: f64, c: f64, d: f64, y0: f64) -> QuarticFactorization {
    let discr = a * y0 * 0.5 - c;
    let sdiscr = if discr >= 0.0 { 1.0 } else { -1.0 };

    // `cxx:457-464`.
    let mut p0 = a * a * 0.25 - b + y0;
    p0 = if p0 < 0.0 { 0.0 } else { p0.sqrt() };

    // `cxx:466-473`.
    let mut q0 = y0 * y0 * 0.25 - d;
    if q0.abs() < 10.0 * MACHINE_EPSILON {
        q0 = 0.0;
    } else {
        q0 = if q0 < 0.0 { 0.0 } else { q0.sqrt() };
    }

    let a_demi = a * 0.5;
    let y_demi = y0 * 0.5;
    let sdiscr_q0 = sdiscr * q0;

    let mut factors = QuarticFactorization {
        p1: a_demi + p0,
        q1: y_demi + sdiscr_q0,
        p2: a_demi - p0,
        q2: y_demi - sdiscr_q0,
    };

    // `cxx:487-503`: clean up near-zero coefficients.
    let eps = 100.0 * MACHINE_EPSILON;
    if factors.p1.abs() <= eps {
        factors.p1 = 0.0;
    }
    if factors.p2.abs() <= eps {
        factors.p2 = 0.0;
    }
    if factors.q1.abs() <= eps {
        factors.q1 = 0.0;
    }
    if factors.q2.abs() <= eps {
        factors.q2 = 0.0;
    }

    factors
}

/// `math_DirectPolynomialRoots` (`hxx:47-226`).
#[derive(Debug, Clone)]
pub struct DirectPolynomialRoots {
    /// `myDone` (`hxx:216`).
    done: bool,
    /// `myInfiniteStatus` (`hxx:217`).
    infinite_status: bool,
    /// `myNbSol` (`hxx:218`).
    nb_sol: usize,
    /// `double myRoots[4]` (`hxx:219`).
    roots: [f64; 4],
}

impl DirectPolynomialRoots {
    /// `math_DirectPolynomialRoots(A, B, C, D, E)` (`cxx:511-521`).
    pub fn quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve_quartic(a, b, c, d, e);
        r
    }

    /// `math_DirectPolynomialRoots(A, B, C, D)` (`cxx:525-533`).
    pub fn cubic(a: f64, b: f64, c: f64, d: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve_cubic(a, b, c, d);
        r
    }

    /// `math_DirectPolynomialRoots(A, B, C)` (`cxx:537-544`).
    pub fn quadratic(a: f64, b: f64, c: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve_quadratic(a, b, c);
        r
    }

    /// `math_DirectPolynomialRoots(A, B)` (`cxx:548-554`).
    pub fn linear(a: f64, b: f64) -> Self {
        let mut r = Self {
            done: true,
            infinite_status: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        r.solve_linear(a, b);
        r
    }

    /// `IsDone` (`hxx:198-201`).
    pub fn is_done(&self) -> bool { self.done }

    /// `InfiniteRoots` (`hxx:203-206`).
    pub fn infinite_roots(&self) -> bool { self.infinite_status }

    /// `NbSolutions` (`hxx:208-212`): OCCT raises `StdFail_InfiniteSolutions`
    /// when the roots are infinite; the port answers `0`.
    pub fn nb_solutions(&self) -> usize {
        if self.infinite_status {
            return 0;
        }
        self.nb_sol
    }

    /// `Value(theIndex)` (`hxx:214-219`), 1-based. OCCT raises
    /// `Standard_RangeError` outside `1..=NbSolutions`; the port answers
    /// `f64::NAN`.
    pub fn value(&self, index: usize) -> f64 {
        if self.infinite_status || index < 1 || index > self.nb_sol {
            return f64::NAN;
        }
        self.roots[index - 1]
    }

    /// The roots `1..=NbSolutions`, in OCCT order.
    pub fn roots(&self) -> &[f64] {
        &self.roots[..self.nb_solutions()]
    }

    /// `Solve(A, B, C, D, E)` (`cxx:556-626`).
    fn solve_quartic(&mut self, a: f64, b: f64, c: f64, d: f64, e: f64) {
        // `cxx:561-565`: check for degree reduction.
        if should_reduce_degree_quartic(a, b, c, d, e) {
            self.solve_cubic(b, c, d, e);
            return;
        }

        // `cxx:569-573`: normalize coefficients.
        let na = b / a;
        let nb = c / a;
        let nc = d / a;
        let nd = e / a;

        // `cxx:576-578`: scale coefficients to avoid overflow/underflow.
        let scaled = ScaledCoefficients::scale_quartic(na, nb, nc, nd);

        // `cxx:580-587`: solve Ferrari's resolvent cubic.
        let Some(y0) = solve_ferrari_resolvent(scaled.a, scaled.b, scaled.c, scaled.d) else {
            self.done = false;
            return;
        };

        // `cxx:590-591`: factor into two quadratics.
        let factors = factor_quartic_via_ferrari(scaled.a, scaled.b, scaled.c, scaled.d, y0);

        // `cxx:594-598`, `cxx:601-605`.
        let q1 = DirectPolynomialRoots::quadratic(1.0, factors.p1, factors.q1);
        if !q1.is_done() {
            self.done = false;
            return;
        }
        let q2 = DirectPolynomialRoots::quadratic(1.0, factors.p2, factors.q2);
        if !q2.is_done() {
            self.done = false;
            return;
        }

        // `cxx:610-624`: collect all roots, then inverse-scale and refine.
        self.nb_sol = q1.nb_solutions() + q2.nb_solutions();
        let mut index = 0;
        for i in 0..q1.nb_solutions() {
            self.roots[index] = q1.roots[i];
            index += 1;
        }
        for i in 0..q2.nb_solutions() {
            self.roots[index] = q2.roots[i];
            index += 1;
        }

        let poly = [a, b, c, d, e];
        for i in 0..self.nb_sol {
            self.roots[i] *= scaled.scale_factor;
            self.roots[i] = refine_root(poly.len(), &poly, self.roots[i]);
        }
    }

    /// `Solve(A, B, C, D)` (`cxx:628-709`).
    fn solve_cubic(&mut self, a: f64, b: f64, c: f64, d: f64) {
        // `cxx:634-638`: check for degree reduction.
        if a.abs() <= ZERO_THRESHOLD {
            self.solve_quadratic(b, c, d);
            return;
        }

        // `cxx:641-644`: normalize coefficients.
        let beta = b / a;
        let gamma = c / a;
        let del = d / a;

        // `cxx:647-649`: scale to avoid overflow/underflow.
        let scaled = ScaledCoefficients::scale_cubic(beta, gamma, del, del);

        // `cxx:651-659`: transform to the depressed cubic t^3 + Pt + Q = 0.
        let p1 = scaled.b;
        let p2 = -(scaled.a * scaled.a) / 3.0;
        let mut p = p1 + p2;
        let ep = 5.0 * MACHINE_EPSILON * (p1.abs() + p2.abs());
        if p.abs() <= ep {
            p = 0.0;
        }

        // `cxx:660-668`.
        let q1 = scaled.c;
        let q2 = -scaled.a * scaled.b / 3.0;
        let q3 = 2.0 * (scaled.a * scaled.a * scaled.a) / 27.0;
        let mut q = q1 + q2 + q3;
        let eq = 10.0 * MACHINE_EPSILON * (q1.abs() + q2.abs() + q3.abs());
        if q.abs() <= eq {
            q = 0.0;
        }

        // `cxx:670-675`: check for overflow.
        if p.abs() > OVERFLOW_LIMIT {
            self.done = false;
            return;
        }

        // `cxx:677-688`: compute the discriminant.
        let a1 = (p * p * p) / 27.0;
        let a2 = (q * q) / 4.0;
        let mut discr = a1 + a2;

        // `cxx:681-685`: special handling for P < 0.
        if p < 0.0 {
            discr = compute_special_discriminant(scaled.a, scaled.b, scaled.c, a1);
        }

        // `cxx:688-708`.
        if discr < 0.0 {
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
            self.nb_sol = 1;
            solve_cubic_one_real_root(scaled.a, scaled.c, p, q, discr, &mut self.roots);
        } else {
            solve_cubic_multiple_roots(
                scaled.a,
                scaled.b,
                scaled.c,
                p,
                q,
                &mut self.roots,
                &mut self.nb_sol,
            );
        }

        // `cxx:707-708`: inverse scaling and Newton-Raphson refinement.
        let poly = [a, b, c, d];
        for i in 0..self.nb_sol {
            self.roots[i] *= scaled.scale_factor;
            self.roots[i] = refine_root(poly.len(), &poly, self.roots[i]);
        }
    }

    /// `Solve(A, B, C)` (`cxx:711-765`).
    fn solve_quadratic(&mut self, a: f64, b: f64, c: f64) {
        // `cxx:716-720`: check for degree reduction.
        if a.abs() <= ZERO_THRESHOLD {
            self.solve_linear(b, c);
            return;
        }

        // `cxx:722-724`: solve the normalized quadratic x^2 + P*x + Q = 0.
        let p = b / a;
        let q = c / a;

        // `cxx:726-734`: compute the discriminant with error bounds.
        let eps_d = 3.0 * MACHINE_EPSILON * (p * p + (4.0 * q).abs());
        let mut discrim = p * p - 4.0 * q;

        if discrim.abs() <= eps_d {
            discrim = 0.0;
        }

        if discrim < 0.0 {
            // `cxx:736-739`: no real roots.
            self.nb_sol = 0;
        } else if discrim == 0.0 {
            // `cxx:740-748`: double root.
            self.nb_sol = 2;
            self.roots[0] = -0.5 * p;
            self.roots[0] = refine_root(3, &[1.0, p, q], self.roots[0]);
            self.roots[1] = self.roots[0];
        } else {
            // `cxx:749-764`: two distinct real roots, stable formula.
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

    /// `Solve(A, B)` (`cxx:767-786`).
    fn solve_linear(&mut self, a: f64, b: f64) {
        if a.abs() <= ZERO_THRESHOLD {
            if b.abs() <= ZERO_THRESHOLD {
                // `cxx:773-777`: 0 = 0, infinite solutions.
                self.infinite_status = true;
                return;
            }
            // `cxx:778-780`: 0*x + B = 0, no solution.
            self.nb_sol = 0;
            return;
        }

        // `cxx:782-785`: unique solution.
        self.nb_sol = 1;
        self.roots[0] = -b / a;
    }
}
