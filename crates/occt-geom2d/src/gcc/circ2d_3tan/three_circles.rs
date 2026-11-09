//! `GccAna_Circ2d3Tan(GccEnt_QualifiedCirc, GccEnt_QualifiedCirc,
//! GccEnt_QualifiedCirc, Tolerance)` (`GccAna_Circ2d3Tan.cxx:31-1017`): circles
//! tangent to three circles.
//!
//! The solver walks the eight `(R +- Ri)` sign systems of `cxx:180-198`:
//!
//! 1. first step (`cxx:229-312`) eliminates `R` from every pair of tangency
//!    equations and stores the coefficients `A2..F2` / `A3..F3` plus the
//!    `Beta/Gamma/Delta` of the second and third circle,
//! 2. second step (`cxx:321-671`) solves the resulting quadratic system per
//!    sign, handling the `a2 == 0` degenerate case and the quartic special case
//!    when two circles touch (`IsTouch`, `cxx:211-218`),
//! 3. third step (`cxx:676-775`) recovers `R` for every `(X, Y)` pair,
//! 4. fourth step (`cxx:777-864`) filters the triplets against the third circle,
//! 5. fifth step (`cxx:866-1017`) applies the qualifiers and computes the
//!    tangency points and parameters.
//!
//! Translation notes:
//! - OCCT stores the intermediate solutions in `NCollection_Array1` of fixed
//!   size (`xSol/ySol(1, 64)`, `xSol1/ySol1/rSol1(1, 32)`, `rSol(1, 16)`); the
//!   port pushes into `Vec` so the same index arithmetic (`FirstSol`/
//!   `FirstSol1`, both 1-based) is safe even where OCCT would write past the
//!   upper bound,
//! - `FirstSol` is overwritten by the fourth step, exactly as OCCT does, so the
//!   fifth step reads the fourth step's boundaries and `FirstSol1` is only ever
//!   read by the fourth step,
//! - the `distcc2`/`distcc3` of `cxx:915` and `cxx:933` are
//!   `Center.Distance(center1)` in OCCT (not `center2`/`center3`); the port
//!   keeps that, since it only feeds the unqualified-qualifier decision.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::math_direct_poly_roots::DirectPolynomialRoots;

use super::super::gcc_ent::GccEntQualifiedCirc;
use super::{axis_at, qualifier_of_circle, GccAnaCirc2d3Tan};

/// `cirsol.Upper()` of the circ/circ/circ constructor (`cxx:42-59`).
const UPPER: usize = 16;

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedCirc,
    qualified3: &GccEntQualifiedCirc,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !((qualified1.is_enclosed() || qualified1.is_enclosing() || qualified1.is_outside()
        || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_enclosing() || qualified2.is_outside()
            || qualified2.is_unqualified())
        && (qualified3.is_enclosed() || qualified3.is_enclosing() || qualified3.is_outside()
            || qualified3.is_unqualified()))
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let cir1 = qualified1.qualified();
    let cir2 = qualified2.qualified();
    let cir3 = qualified3.qualified();
    let r1 = cir1.radius();
    let r2 = cir2.radius();
    let r3 = cir3.radius();
    let center1 = cir1.location();
    let center2 = cir2.location();
    let center3 = cir3.location();

    let x1 = center1.x();
    let x2 = center2.x();
    let x3 = center3.x();
    let y1 = center1.y();
    let y2 = center2.y();
    let y3 = center3.y();

    let dir2x = x1 - x2;
    let dir2y = y1 - y2;
    let dir3x = x1 - x3;
    let dir3y = y1 - y3;

    // `cxx:100-144`: the degenerate configurations return without a solution.
    if ((r1 - r2).abs() <= tolerance && center1.is_equal(&center2, tolerance))
        || ((r1 - r3).abs() <= tolerance && center1.is_equal(&center3, tolerance))
        || ((r2 - r3).abs() <= tolerance && center2.is_equal(&center3, tolerance))
    {
        return me;
    }
    if (dir2x * dir3y - dir2y * dir3x).abs() <= tolerance {
        let dist1 = center1.distance(&center2);
        let dist2 = center1.distance(&center3);
        let dist3 = center2.distance(&center3);
        if ((r1 - r2).abs() - dist1).abs() <= tolerance {
            if ((r1 - r3).abs() - dist2).abs() <= tolerance {
                if ((r2 - r3).abs() - dist3).abs() <= tolerance {
                    return me;
                }
            } else if (r1 + r3 - dist2).abs() <= tolerance {
                if (r2 + r3 - dist3).abs() <= tolerance {
                    return me;
                }
            }
        } else if (r1 + r2 - dist1).abs() <= tolerance {
            if ((r1 - r3).abs() - dist2).abs() <= tolerance
                && (r2 + r3 - dist3).abs() <= tolerance
            {
                // `cxx:131-137`: both conditions hold, nothing to reject.
            } else if ((r2 - r3).abs() - dist3).abs() <= tolerance
                && (r1 + r3 - dist2).abs() <= tolerance
            {
                return me;
            }
        }
    }

    let mut a2 = [0.0f64; 9];
    let mut b2 = [0.0f64; 9];
    let mut c2 = [0.0f64; 9];
    let mut d2 = [0.0f64; 9];
    let mut e2 = [0.0f64; 9];
    let mut f2 = [0.0f64; 9];
    let mut a3 = [0.0f64; 9];
    let mut b3 = [0.0f64; 9];
    let mut c3 = [0.0f64; 9];
    let mut d3 = [0.0f64; 9];
    let mut e3 = [0.0f64; 9];
    let mut f3 = [0.0f64; 9];
    let mut beta2 = [0.0f64; 9];
    let mut gamma2 = [0.0f64; 9];
    let mut delta2 = [0.0f64; 9];
    let mut beta3 = [0.0f64; 9];
    let mut gamma3 = [0.0f64; 9];
    let mut delta3 = [0.0f64; 9];

    // `cxx:211-218`.
    let is_touch = ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) - (r1 - r2) * (r1 - r2)).abs()
        <= tolerance
        || ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) - (r1 + r2) * (r1 + r2)).abs()
            <= tolerance
        || ((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) - (r1 - r3) * (r1 - r3)).abs()
            <= tolerance
        || ((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) - (r1 + r3) * (r1 + r3)).abs()
            <= tolerance
        || ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) - (r2 - r3) * (r2 - r3)).abs()
            <= tolerance
        || ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) - (r2 + r3) * (r2 + r3)).abs()
            <= tolerance;

    // First step (`cxx:229-312`).
    for i in 1..=8 {
        let minus2 = i == 1 || i == 4 || i == 5 || i == 8;
        if minus2 {
            if (r1 - r2).abs() > tolerance {
                beta2[i] = (x1 - x2) / (r1 - r2);
                gamma2[i] = (y1 - y2) / (r1 - r2);
                delta2[i] = (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + (r1 - r2) * (r1 - r2))
                    / (2.0 * (r1 - r2));
            }
        } else {
            beta2[i] = (x1 - x2) / (r1 + r2);
            gamma2[i] = (y1 - y2) / (r1 + r2);
            delta2[i] = (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + (r1 + r2) * (r1 + r2))
                / (2.0 * (r1 + r2));
        }
        if minus2 && (r1 - r2).abs() <= tolerance {
            a2[i] = 0.0;
            b2[i] = 0.0;
            c2[i] = 0.0;
            d2[i] = x2 - x1;
            e2[i] = y2 - y1;
            f2[i] = x1 * x1 - x2 * x2 + y1 * y1 - y2 * y2;
        } else {
            a2[i] = beta2[i] * beta2[i] - 1.0;
            b2[i] = beta2[i] * gamma2[i];
            c2[i] = gamma2[i] * gamma2[i] - 1.0;
            d2[i] = beta2[i] * delta2[i] + x2;
            e2[i] = gamma2[i] * delta2[i] + y2;
            f2[i] = delta2[i] * delta2[i] - x2 * x2 - y2 * y2;
        }

        let minus3 = i == 1 || i == 3 || i == 6 || i == 8;
        if minus3 {
            if (r1 - r3).abs() > tolerance {
                beta3[i] = (x1 - x3) / (r1 - r3);
                gamma3[i] = (y1 - y3) / (r1 - r3);
                delta3[i] = (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + (r1 - r3) * (r1 - r3))
                    / (2.0 * (r1 - r3));
            }
        } else {
            beta3[i] = (x1 - x3) / (r1 + r3);
            gamma3[i] = (y1 - y3) / (r1 + r3);
            delta3[i] = (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + (r1 + r3) * (r1 + r3))
                / (2.0 * (r1 + r3));
        }
        if minus3 && (r1 - r3).abs() <= tolerance {
            a3[i] = 0.0;
            b3[i] = 0.0;
            c3[i] = 0.0;
            d3[i] = x3 - x1;
            e3[i] = y3 - y1;
            f3[i] = x1 * x1 - x3 * x3 + y1 * y1 - y3 * y3;
        } else {
            a3[i] = beta3[i] * beta3[i] - 1.0;
            b3[i] = beta3[i] * gamma3[i];
            c3[i] = gamma3[i] * gamma3[i] - 1.0;
            d3[i] = beta3[i] * delta3[i] + x3;
            e3[i] = gamma3[i] * delta3[i] + y3;
            f3[i] = delta3[i] * delta3[i] - x3 * x3 - y3 * y3;
        }
    }

    // Second step (`cxx:321-671`). `xsol`/`ysol` are OCCT's `xSol`/`ySol`,
    // `first_sol` its `FirstSol`; both keep the 1-based indexing.
    let mut xsol: Vec<f64> = Vec::new();
    let mut ysol: Vec<f64> = Vec::new();
    let mut first_sol = [0usize; 10];

    for i in 1..=8 {
        let a2i = a2[i];
        let a3i = a3[i];
        let b2i = b2[i];
        let b3i = b3[i];
        let c2i = c2[i];
        let c3i = c3[i];
        let d2i = d2[i];
        let d3i = d3[i];
        let e2i = e2[i];
        let e3i = e3[i];
        let f2i = f2[i];
        let f3i = f3[i];

        first_sol[i] = xsol.len() + 1;

        if ((i == 2 || i == 5 || i == 6 || i == 8)
            && (qualified1.is_enclosed() || qualified1.is_enclosing()))
            || ((i == 1 || i == 3 || i == 4 || i == 7) && qualified1.is_outside())
        {
            continue;
        }
        if ((i == 3 || i == 5 || i == 7 || i == 8)
            && (qualified2.is_enclosed() || qualified2.is_enclosing()))
            || ((i == 1 || i == 2 || i == 4 || i == 6) && qualified2.is_outside())
        {
            continue;
        }
        if ((i == 4 || i == 6 || i == 7 || i == 8)
            && (qualified3.is_enclosed() || qualified3.is_enclosing()))
            || ((i == 1 || i == 2 || i == 3 || i == 5) && qualified3.is_outside())
        {
            continue;
        }

        // `cxx:361-372`: the first circle itself solves the system.
        if (a2i - a3i).abs() <= tolerance
            && (b2i - b3i).abs() <= tolerance
            && (c2i - c3i).abs() <= tolerance
            && (d2i - d3i).abs() <= tolerance
            && (e2i - e3i).abs() <= tolerance
            && (f2i - f3i).abs() <= tolerance
        {
            xsol.push(x1);
            ysol.push(y1);
            continue;
        }

        if a2i.abs() <= tolerance {
            // 1) a2 = 0 (`cxx:373-488`).
            // 1.1) b2*y + d2 = 0 (`cxx:374-411`).
            let y_roots = DirectPolynomialRoots::new3(c2i, 2.0 * e2i, f2i);
            if y_roots.is_done() && !y_roots.infinite_roots() {
                for k in 1..=y_roots.nb_solutions() {
                    let y = y_roots.value(k);
                    if (k != 2 || (y - y_roots.value(1)).abs() > 10.0 * tolerance)
                        && (b2i * y + d2i).abs() <= b2i * tolerance
                    {
                        let x_roots = DirectPolynomialRoots::new3(
                            a3i,
                            2.0 * (b3i * y + d3i),
                            c3i * (y * y) + 2.0 * e3i * y + f3i,
                        );
                        if x_roots.is_done() && !x_roots.infinite_roots() {
                            for j in 1..=x_roots.nb_solutions() {
                                let x = x_roots.value(j);
                                if j != 2 || (x - x_roots.value(1)).abs() > 10.0 * tolerance {
                                    xsol.push(x);
                                    ysol.push(y);
                                }
                            }
                        }
                    }
                }
            }

            // 1.2) b2*y + d2 != 0 (`cxx:413-488`).
            let a = a3i * c2i * c2i - 4.0 * b2i * (b3i * c2i - b2i * c3i);
            let b = 4.0 * a3i * c2i * e2i - 4.0 * b3i * (c2i * d2i + 2.0 * b2i * e2i)
                + 4.0 * b2i * (2.0 * c3i * d2i - c2i * d3i + 2.0 * b2i * e3i);
            let c = 2.0 * a3i * (c2i * f2i + 2.0 * e2i * e2i)
                - 4.0 * b3i * (b2i * f2i + 2.0 * e2i * d2i)
                + 4.0 * c3i * d2i * d2i
                - 4.0 * d3i * (c2i * d2i + 2.0 * b2i * e2i)
                + 16.0 * b2i * e3i * d2i
                + 4.0 * b2i * b2i * f3i;
            let d = 4.0 * a3i * e2i * f2i - 4.0 * b3i * d2i * f2i
                - 4.0 * d3i * (b2i * f2i + 2.0 * d2i * e2i)
                + 8.0 * d2i * d2i * e3i
                + 8.0 * b2i * d2i * f3i;
            let e = a3i * f2i * f2i - 4.0 * d2i * d3i * f2i + 4.0 * d2i * d2i * f3i;

            if is_touch {
                // Derivation of `A*y^4 + ... + E` (`cxx:427-461`).
                let y_roots1 = DirectPolynomialRoots::new4(4.0 * a, 3.0 * b, 2.0 * c, d);
                if y_roots1.is_done() && !y_roots1.infinite_roots() {
                    for k in 1..=y_roots1.nb_solutions() {
                        let y = y_roots1.value(k);
                        let mut is_same = false;
                        for l in 1..k {
                            if (y - y_roots1.value(l)).abs() <= 10.0 * tolerance {
                                is_same = true;
                            }
                        }
                        let epsilon = (((4.0 * a * y).abs() + (3.0 * b).abs()) * y).abs()
                            + (2.0 * c).abs();
                        let epsilon = (epsilon * y).abs() + d.abs();
                        if ((((a * y + b) * y + c) * y + d) * y + e).abs()
                            <= epsilon * tolerance
                        {
                            if !is_same && (b2i * y + d2i).abs() > b2i * tolerance {
                                let x = -(c2i * (y * y) + 2.0 * e2i * y + f2i) / (2.0 * (b2i * y + d2i));
                                xsol.push(x);
                                ysol.push(y);
                            }
                        }
                    }
                }
            }

            let y_roots1 = DirectPolynomialRoots::new5(a, b, c, d, e);
            if y_roots1.is_done() && !y_roots1.infinite_roots() {
                for k in 1..=y_roots1.nb_solutions() {
                    let y = y_roots1.value(k);
                    let mut is_same = false;
                    let first_index = if i == 1 { 1 } else { first_sol[i] };
                    for l in first_index..=xsol.len() {
                        if (y - ysol[l - 1]).abs() <= 10.0 * tolerance {
                            is_same = true;
                        }
                    }
                    if !is_same && (b2i * y + d2i).abs() > b2i * tolerance {
                        let x = -(c2i * (y * y) + 2.0 * e2i * y + f2i) / (2.0 * (b2i * y + d2i));
                        xsol.push(x);
                        ysol.push(y);
                    }
                }
            }
        } else {
            // 2) a2 != 0 (`cxx:489-671`).
            let m = 2.0 * a3i * b2i * b2i / (a2i * a2i) - 2.0 * b2i * b3i / a2i
                - a3i * c2i / a2i
                + c3i;
            let n = 4.0 * a3i * b2i * d2i / (a2i * a2i) - 2.0 * b3i * d2i / a2i
                - 2.0 * b2i * d3i / a2i
                - 2.0 * a3i * e2i / a2i
                + 2.0 * e3i;
            let t = 2.0 * a3i * d2i * d2i / (a2i * a2i) - 2.0 * d2i * d3i / a2i
                - a3i * f2i / a2i
                + f3i;
            let s = 2.0 * b3i - 2.0 * a3i * b2i / a2i;
            let v = 2.0 * d3i - 2.0 * d2i * a3i / a2i;

            if s.abs() <= tolerance && v.abs() <= tolerance {
                let y_roots = DirectPolynomialRoots::new3(m, n, t);
                if y_roots.is_done() && !y_roots.infinite_roots() {
                    for k in 1..=y_roots.nb_solutions() {
                        let y = y_roots.value(k);
                        let p = -(b2i * y + d2i) / a2i;
                        let q = (c2i * (y * y) + 2.0 * e2i * y + f2i) / a2i;
                        let epsilon = pq_epsilon(a2i, b2i, c2i, d2i, e2i, y);
                        if (k != 2 || (y - y_roots.value(1)).abs() > 10.0 * tolerance)
                            && p * p - q >= -epsilon * tolerance
                        {
                            let x_roots = DirectPolynomialRoots::new3(
                                a2i,
                                2.0 * (b2i * y + d2i),
                                c2i * y * y + 2.0 * e2i * y + f2i,
                            );
                            if x_roots.is_done() && !x_roots.infinite_roots() {
                                for l in 1..=x_roots.nb_solutions() {
                                    let x = x_roots.value(l);
                                    if l != 2 || (x - x_roots.value(1)).abs() > 10.0 * tolerance {
                                        xsol.push(x);
                                        ysol.push(y);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                let a = s * s * (b2i * b2i - a2i * c2i) - m * m * a2i * a2i;
                let b = 2.0 * s * v * (b2i * b2i - a2i * c2i)
                    + 2.0 * s * s * (b2i * d2i - a2i * e2i)
                    - 2.0 * m * n * a2i * a2i;
                let c = v * v * (b2i * b2i - a2i * c2i)
                    + 4.0 * s * v * (b2i * d2i - a2i * e2i)
                    + s * s * (d2i * d2i - a2i * f2i)
                    - a2i * a2i * (2.0 * m * t + n * n);
                let d = 2.0 * v * v * (b2i * d2i - a2i * e2i)
                    + 2.0 * s * v * (d2i * d2i - a2i * f2i)
                    - 2.0 * n * t * a2i * a2i;
                let e = v * v * (d2i * d2i - a2i * f2i) - t * t * a2i * a2i;

                if is_touch {
                    let y_roots1 = DirectPolynomialRoots::new4(4.0 * a, 3.0 * b, 2.0 * c, d);
                    if y_roots1.is_done() && !y_roots1.infinite_roots() {
                        for k in 1..=y_roots1.nb_solutions() {
                            let y = y_roots1.value(k);
                            let p = -(b2i * y + d2i) / a2i;
                            let q = (c2i * (y * y) + 2.0 * e2i * y + f2i) / a2i;
                            let mut is_same = false;
                            let first_index = if i == 1 { 1 } else { first_sol[i] };
                            for l in first_index..=xsol.len() {
                                if (y - ysol[l - 1]).abs() <= 10.0 * tolerance {
                                    is_same = true;
                                }
                            }
                            let epsilon_q = (((4.0 * a * y).abs() + (3.0 * b).abs()) * y).abs()
                                + (2.0 * c).abs();
                            let epsilon_q = (epsilon_q * y).abs() + d.abs();
                            if ((((a * y + b) * y + c) * y + d) * y + e).abs()
                                <= epsilon_q * tolerance
                            {
                                let epsilon = pq_epsilon(a2i, b2i, c2i, d2i, e2i, y);
                                if !is_same && p * p - q >= -epsilon * tolerance {
                                    let x_roots = DirectPolynomialRoots::new3(
                                        a2i,
                                        2.0 * (b2i * y + d2i),
                                        c2i * y * y + 2.0 * e2i * y + f2i,
                                    );
                                    if x_roots.is_done() && !x_roots.infinite_roots() {
                                        for l in 1..=x_roots.nb_solutions() {
                                            let x = x_roots.value(l);
                                            if l != 2
                                                || (x - x_roots.value(1)).abs()
                                                    > 10.0 * tolerance
                                            {
                                                xsol.push(x);
                                                ysol.push(y);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                let y_roots = DirectPolynomialRoots::new5(a, b, c, d, e);
                if y_roots.is_done() && !y_roots.infinite_roots() {
                    for k in 1..=y_roots.nb_solutions() {
                        let y = y_roots.value(k);
                        let p = -(b2i * y + d2i) / a2i;
                        let q = (c2i * (y * y) + 2.0 * e2i * y + f2i) / a2i;
                        let mut is_same = false;
                        for l in 1..k {
                            if (y - y_roots.value(l)).abs() <= 10.0 * tolerance {
                                is_same = true;
                            }
                        }
                        let epsilon = pq_epsilon(a2i, b2i, c2i, d2i, e2i, y);
                        if !is_same && p * p - q >= -epsilon * tolerance {
                            let x_roots = DirectPolynomialRoots::new3(
                                a2i,
                                2.0 * (b2i * y + d2i),
                                c2i * y * y + 2.0 * e2i * y + f2i,
                            );
                            if x_roots.is_done() && !x_roots.infinite_roots() {
                                for l in 1..=x_roots.nb_solutions() {
                                    let x = x_roots.value(l);
                                    if l != 2 || (x - x_roots.value(1)).abs() > 10.0 * tolerance {
                                        xsol.push(x);
                                        ysol.push(y);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    first_sol[9] = xsol.len() + 1;

    // Third step (`cxx:676-775`): recover `R`, `r > 0`.
    let mut xsol1: Vec<f64> = Vec::new();
    let mut ysol1: Vec<f64> = Vec::new();
    let mut rsol1: Vec<f64> = Vec::new();
    let mut first_sol1 = [0usize; 10];
    for i in 1..=8 {
        first_sol1[i] = xsol1.len() + 1;
        for j in first_sol[i]..first_sol[i + 1] {
            let x = xsol[j - 1];
            let y = ysol[j - 1];
            if (i == 1 || i == 4 || i == 5 || i == 8) && (r1 - r2).abs() <= tolerance {
                if i == 1 || i == 4 {
                    let mut r = r1 + ((x - x1) * (x - x1) + (y - y1) * (y - y1)).sqrt();
                    let mut epsilon =
                        10.0 * (2.0 * (r - r2).abs() + (x - x2).abs() + (y - y2).abs());
                    if ((r - r2) * (r - r2) - (x - x2) * (x - x2) - (y - y2) * (y - y2)).abs()
                        <= epsilon * tolerance
                    {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                    r = r1 - ((x - x1) * (x - x1) + (y - y1) * (y - y1)).sqrt();
                    epsilon = 10.0 * (2.0 * (r - r2).abs() + (x - x2).abs() + (y - y2).abs());
                    if r > tolerance
                        && ((r - r2) * (r - r2) - (x - x2) * (x - x2) - (y - y2) * (y - y2))
                            .abs()
                            <= epsilon * tolerance
                    {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                } else {
                    // `i == 5 || i == 8`.
                    let r = -r1 + ((x - x1) * (x - x1) + (y - y1) * (y - y1)).sqrt();
                    if r > tolerance {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                }
            } else {
                if i == 1 || i == 4 {
                    let r = r2 + beta2[i] * x + gamma2[i] * y + delta2[i];
                    if r > tolerance {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                }
                if i == 5 || i == 8 {
                    let r = -r2 - beta2[i] * x - gamma2[i] * y - delta2[i];
                    if r > tolerance {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                }
                if i == 3 || i == 7 {
                    let r = -r2 + beta2[i] * x + gamma2[i] * y + delta2[i];
                    if r > tolerance {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                }
                if i == 2 || i == 6 {
                    let r = r2 - beta2[i] * x - gamma2[i] * y - delta2[i];
                    if r > tolerance {
                        xsol1.push(x);
                        ysol1.push(y);
                        rsol1.push(r);
                    }
                }
            }
        }
    }
    first_sol1[9] = xsol1.len() + 1;

    // Fourth step (`cxx:777-864`): keep the triplets consistent with circle 3.
    let mut xsol2: Vec<f64> = Vec::new();
    let mut ysol2: Vec<f64> = Vec::new();
    let mut rsol2: Vec<f64> = Vec::new();
    for i in 1..=8 {
        first_sol[i] = xsol2.len() + 1;
        for j in first_sol1[i]..first_sol1[i + 1] {
            let x = xsol1[j - 1];
            let y = ysol1[j - 1];
            let r = rsol1[j - 1];
            if (i == 1 || i == 3 || i == 6 || i == 8) && (r1 - r3).abs() <= tolerance {
                if i == 1 || i == 3 {
                    let epsilon =
                        10.0 * (2.0 * (r - r3).abs() + (x - x3).abs() + (y - y3).abs());
                    if ((r - r3) * (r - r3) - (x - x3) * (x - x3) - (y - y3) * (y - y3)).abs()
                        <= epsilon * tolerance
                    {
                        xsol2.push(x);
                        ysol2.push(y);
                        rsol2.push(r);
                    }
                } else {
                    // `i == 6 || i == 8`.
                    let epsilon = 10.0 * (2.0 * (r + r3) + (x - x3).abs() + (y - y3).abs());
                    if ((r + r3) * (r + r3) - (x - x3) * (x - x3) - (y - y3) * (y - y3)).abs()
                        <= epsilon * tolerance
                    {
                        xsol2.push(x);
                        ysol2.push(y);
                        rsol2.push(r);
                    }
                }
            } else {
                let epsilon = 10.0 * (beta3[i].abs() + gamma3[i].abs() + 1.0);
                if i == 1 || i == 3 {
                    if (r3 + beta3[i] * x + gamma3[i] * y + delta3[i] - r).abs()
                        <= epsilon * tolerance
                    {
                        xsol2.push(x);
                        ysol2.push(y);
                        rsol2.push(r);
                    }
                }
                if i == 6 || i == 8 {
                    if (r3 + beta3[i] * x + gamma3[i] * y + delta3[i] + r).abs()
                        <= epsilon * tolerance
                    {
                        xsol2.push(x);
                        ysol2.push(y);
                        rsol2.push(r);
                    }
                }
                if i == 4 || i == 7 {
                    if (beta3[i] * x + gamma3[i] * y + delta3[i] - r - r3).abs()
                        <= epsilon * tolerance
                    {
                        xsol2.push(x);
                        ysol2.push(y);
                        rsol2.push(r);
                    }
                }
                if i == 2 || i == 5 {
                    if (r - r3 + beta3[i] * x + gamma3[i] * y + delta3[i]).abs()
                        <= epsilon * tolerance
                    {
                        xsol2.push(x);
                        ysol2.push(y);
                        rsol2.push(r);
                    }
                }
            }
        }
    }
    first_sol[9] = xsol2.len() + 1;

    // Fifth step (`cxx:866-1017`): qualifiers, tangency points, parameters.
    for i in 1..=8 {
        for j in first_sol[i]..first_sol[i + 1] {
            let r = rsol2[j - 1];
            if (qualified1.is_enclosed() && r > r1) || (qualified1.is_enclosing() && r < r1) {
                continue;
            }
            if (qualified2.is_enclosed() && r > r2) || (qualified2.is_enclosing() && r < r2) {
                continue;
            }
            if (qualified3.is_enclosed() && r > r3) || (qualified3.is_enclosing() && r < r3) {
                continue;
            }

            let mut n = me.nbr_sol + 1;
            if n > UPPER {
                n = UPPER;
            }
            me.nbr_sol = n;
            let k = n - 1;

            let center = GpPnt2d::new(xsol2[j - 1], ysol2[j - 1]);
            me.cirsol[k] = GpCirc2d::new(axis_at(center), r);

            let distcc1 = center.distance(&center1);
            me.qualifier1[k] = qualifier_of_circle(qualified1, distcc1, r, r1, tol);
            // `cxx:915`: OCCT measures against `center1` here.
            let distcc2 = center.distance(&center1);
            me.qualifier2[k] = qualifier_of_circle(qualified2, distcc2, r, r2, tol);
            // `cxx:933`: same slip, against `center1`.
            let distcc3 = center.distance(&center1);
            me.qualifier3[k] = qualifier_of_circle(qualified3, distcc3, r, r3, tol);

            if center.distance(&cir1.location()) <= tolerance {
                me.the_same1[k] = 1;
            } else {
                me.the_same1[k] = 0;
                let dc = if (i == 2 || i == 5 || i == 6 || i == 8) || r > r1 {
                    GpDir2d::new(cir1.location().x() - center.x(), cir1.location().y() - center.y())
                        .expect("gp_Dir2d from circle centre to solution centre")
                } else {
                    GpDir2d::new(center.x() - cir1.location().x(), center.y() - cir1.location().y())
                        .expect("gp_Dir2d from solution centre to circle centre")
                };
                me.pnttg1sol[k] = GpPnt2d::new(center.x() + r * dc.x, center.y() + r * dc.y);
                me.par1sol[k] = parameter_circ2d(&me.cirsol[k], &me.pnttg1sol[k]);
                me.pararg1[k] = parameter_circ2d(&cir1, &me.pnttg1sol[k]);
            }

            if center.distance(&cir2.location()) <= tolerance {
                me.the_same2[k] = 1;
            } else {
                me.the_same2[k] = 0;
                let dc = if (i == 3 || i == 5 || i == 7 || i == 8) || r > r2 {
                    GpDir2d::new(cir2.location().x() - center.x(), cir2.location().y() - center.y())
                        .expect("gp_Dir2d from circle centre to solution centre")
                } else {
                    GpDir2d::new(center.x() - cir2.location().x(), center.y() - cir2.location().y())
                        .expect("gp_Dir2d from solution centre to circle centre")
                };
                me.pnttg2sol[k] = GpPnt2d::new(center.x() + r * dc.x, center.y() + r * dc.y);
                me.par2sol[k] = parameter_circ2d(&me.cirsol[k], &me.pnttg2sol[k]);
                me.pararg2[k] = parameter_circ2d(&cir2, &me.pnttg2sol[k]);
            }

            if center.distance(&cir3.location()) <= tolerance {
                me.the_same3[k] = 1;
            } else {
                me.the_same3[k] = 0;
                let dc = if (i == 4 || i == 6 || i == 7 || i == 8) || r > r3 {
                    GpDir2d::new(cir3.location().x() - center.x(), cir3.location().y() - center.y())
                        .expect("gp_Dir2d from circle centre to solution centre")
                } else {
                    GpDir2d::new(center.x() - cir3.location().x(), center.y() - cir3.location().y())
                        .expect("gp_Dir2d from solution centre to circle centre")
                };
                me.pnttg3sol[k] = GpPnt2d::new(center.x() + r * dc.x, center.y() + r * dc.y);
                me.par3sol[k] = parameter_circ2d(&me.cirsol[k], &me.pnttg3sol[k]);
                me.pararg3[k] = parameter_circ2d(&cir3, &me.pnttg3sol[k]);
            }
        }
    }
    me.well_done = true;
    me
}

/// `Epsilon = 2. * (|(b2*b2 + |a2*c2|) * y| + |b2*d2| + |a2*e2|) / (a2*a2)`
/// (`GccAna_Circ2d3Tan.cxx:513-517`, `:589-592`, `:641-645`).
fn pq_epsilon(a2: f64, b2: f64, c2: f64, d2: f64, e2: f64, y: f64) -> f64 {
    2.0 * (((b2 * b2 + (a2 * c2).abs()) * y).abs() + (b2 * d2).abs() + (a2 * e2).abs())
        / (a2 * a2)
}
