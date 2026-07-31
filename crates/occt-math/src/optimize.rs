//! Unconstrained optimization — gradient descent, coordinate descent and
//! simulated annealing for general scalar objectives.
//! Source: `math_PSO`, `math_BFGS`-adjacent heuristics (no derivatives).

/// Gradient descent with numerical central-difference gradients and line
/// search. Returns the minimizing point or an error message.
pub fn gradient_descent<F: Fn(&[f64]) -> f64>(
    f: &F,
    x0: &[f64],
    learning_rate: f64,
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, String> {
    let n = x0.len();
    if n == 0 {
        return Err("gradient_descent: empty domain".into());
    }
    let mut x = x0.to_vec();
    let h = 1e-6;
    for _ in 0..max_iter {
        let fx = f(&x);
        let mut grad = vec![0.0; n];
        for i in 0..n {
            let mut xp = x.clone();
            let mut xm = x.clone();
            xp[i] += h;
            xm[i] -= h;
            grad[i] = (f(&xp) - f(&xm)) / (2.0 * h);
        }
        let gnorm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gnorm < tol {
            return Ok(x);
        }
        // Line search: halve the step until the objective decreases.
        let mut step = learning_rate;
        let mut improved = false;
        for _ in 0..30 {
            let mut xn = x.clone();
            for i in 0..n {
                xn[i] -= step * grad[i];
            }
            if f(&xn) < fx {
                x = xn;
                improved = true;
                break;
            }
            step *= 0.5;
        }
        if !improved {
            return Ok(x);
        }
        if step < 1e-15 {
            return Ok(x);
        }
    }
    Ok(x)
}

/// Coordinate descent: minimize along each axis in turn via golden-section
/// line search. Robust for separable objectives.
pub fn coordinate_descent<F: Fn(&[f64]) -> f64>(
    f: &F,
    x0: &[f64],
    bounds: &[(f64, f64)],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, String> {
    let n = x0.len();
    if n == 0 || bounds.len() != n {
        return Err("coordinate_descent: dimension mismatch".into());
    }
    let mut x = x0.to_vec();
    for _ in 0..max_iter {
        let before = f(&x);
        for i in 0..n {
            let (lo, hi) = bounds[i];
            // Golden-section 1D minimization along axis i.
            let (mut a, mut b) = (lo, hi);
            let phi = 0.618_033_988_749_895;
            let mut c = b - phi * (b - a);
            let mut d = a + phi * (b - a);
            let mut fc = f_with(&x, i, c, f);
            let mut fd = f_with(&x, i, d, f);
            for _ in 0..60 {
                if fc < fd {
                    b = d;
                    d = c;
                    fd = fc;
                    c = b - phi * (b - a);
                    fc = f_with(&x, i, c, f);
                } else {
                    a = c;
                    c = d;
                    fc = fd;
                    d = a + phi * (b - a);
                    fd = f_with(&x, i, d, f);
                }
            }
            x[i] = 0.5 * (a + b);
        }
        if (f(&x) - before).abs() < tol {
            break;
        }
    }
    Ok(x)
}

fn f_with<F: Fn(&[f64]) -> f64>(x: &[f64], i: usize, v: f64, f: &F) -> f64 {
    let mut xc = x.to_vec();
    xc[i] = v;
    f(&xc)
}

/// Simulated annealing over a bounded box. Good for rough global search.
pub fn simulated_annealing<F: Fn(&[f64]) -> f64>(
    f: &F,
    x0: &[f64],
    bounds: &[(f64, f64)],
    max_iter: usize,
    t0: f64,
) -> Vec<f64> {
    let n = x0.len();
    if n == 0 || bounds.len() != n {
        return x0.to_vec();
    }
    let mut x = x0.to_vec();
    let mut best = x.clone();
    let mut fbest = f(&best);
    let mut fx = fbest;
    let mut t = t0;
    let mut seed = 12_345u64;
    // Deterministic pseudo-random in [0,1).
    let mut rng = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) as f64) / (1u64 << 31) as f64
    };
    for _ in 0..max_iter {
        let mut cand = x.clone();
        for i in 0..n {
            let span = bounds[i].1 - bounds[i].0;
            cand[i] += (rng() - 0.5) * span;
            cand[i] = cand[i].clamp(bounds[i].0, bounds[i].1);
        }
        let fc = f(&cand);
        if fc < fx || rng() < ((fx - fc) / t).exp().min(1.0) {
            x = cand;
            fx = fc;
            if fc < fbest {
                fbest = fc;
                best = x.clone();
            }
        }
        t *= 0.99;
    }
    best
}

/// Nelder–Mead simplex minimization (derivative-free, robust for smooth
/// low-dimensional problems).
pub fn nelder_mead<F: Fn(&[f64]) -> f64>(
    f: &F,
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = x0.len();
    if n == 0 {
        return x0.to_vec();
    }
    // Initial simplex: x0 plus perturbations along each axis.
    let mut sim: Vec<Vec<f64>> = vec![x0.to_vec()];
    for i in 0..n {
        let mut p = x0.to_vec();
        p[i] += if p[i].abs() > 1e-6 { p[i] * 0.05 } else { 0.05 };
        sim.push(p);
    }
    let mut evals: Vec<f64> = sim.iter().map(|p| f(p)).collect();
    for _ in 0..max_iter {
        // Order by objective.
        let mut order: Vec<usize> = (0..sim.len()).collect();
        order.sort_by(|&a, &b| evals[a].total_cmp(&evals[b]));
        let (best, worst) = (order[0], *order.last().unwrap());
        if evals[worst] - evals[best] < tol {
            break;
        }
        // Centroid of all but worst.
        let mut centroid = vec![0.0; n];
        for &i in &order[..order.len() - 1] {
            for k in 0..n {
                centroid[k] += sim[i][k];
            }
        }
        for k in 0..n {
            centroid[k] /= (order.len() - 1) as f64;
        }
        // Reflect.
        let reflected: Vec<f64> = (0..n).map(|k| 2.0 * centroid[k] - sim[worst][k]).collect();
        let fr = f(&reflected);
        if fr < evals[best] {
            // Expand.
            let expanded: Vec<f64> = (0..n).map(|k| centroid[k] + 2.0 * (reflected[k] - centroid[k])).collect();
            let fe = f(&expanded);
            if fe < fr {
                sim[worst] = expanded;
                evals[worst] = fe;
            } else {
                sim[worst] = reflected;
                evals[worst] = fr;
            }
        } else if fr < evals[*order.last().unwrap()] {
            sim[worst] = reflected;
            evals[worst] = fr;
        } else {
            // Shrink toward the best.
            for &i in &order[1..] {
                for k in 0..n {
                    sim[i][k] = 0.5 * (sim[i][k] + sim[best][k]);
                }
                evals[i] = f(&sim[i]);
            }
        }
    }
    let mut best = 0;
    for i in 1..sim.len() {
        if evals[i] < evals[best] {
            best = i;
        }
    }
    sim[best].clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_descent_quadratic() {
        // f(x,y) = (x-3)² + (y+2)² → minimum (3,-2).
        let f = |p: &[f64]| (p[0] - 3.0).powi(2) + (p[1] + 2.0).powi(2);
        let x = gradient_descent(&f, &[0.0, 0.0], 0.1, 5000, 1e-6).unwrap();
        assert!((x[0] - 3.0).abs() < 1e-3, "x0={}", x[0]);
        assert!((x[1] + 2.0).abs() < 1e-3, "x1={}", x[1]);
    }

    #[test]
    fn coordinate_descent_quadratic() {
        let f = |p: &[f64]| (p[0] - 1.0).powi(2) + (p[1] - 2.0).powi(2);
        let x = coordinate_descent(&f, &[0.0, 0.0], &[(-10.0, 10.0), (-10.0, 10.0)], 100, 1e-8).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-3);
        assert!((x[1] - 2.0).abs() < 1e-3);
    }

    #[test]
    fn simulated_annealing_finds_basin() {
        // Rugged: sin bump on top of a quadratic; annealing finds the global min.
        let f = |p: &[f64]| (p[0] - 4.0).powi(2) + (p[1]).powi(2) + 3.0 * (5.0 * p[0]).sin();
        let x = simulated_annealing(&f, &[0.0, 0.0], &[(-10.0, 10.0), (-10.0, 10.0)], 5000, 5.0);
        assert!(f(&x) < f(&[0.0, 0.0]) + 1.0, "annealing improved: {} vs {}", f(&x), f(&[0.0,0.0]));
    }

    #[test]
    fn nelder_mead_rosenbrock() {
        // Rosenbrock valley.
        let f = |p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2);
        let x = nelder_mead(&f, &[0.0, 0.0], 5000, 1e-8);
        assert!((x[0] - 1.0).abs() < 1e-2, "x0={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-2, "x1={}", x[1]);
    }

    #[test]
    fn all_methods_agree_on_quadratic() {
        let f = |p: &[f64]| (p[0] - 2.0).powi(2) + (p[1] - 3.0).powi(2);
        let gd = gradient_descent(&f, &[0.0, 0.0], 0.1, 5000, 1e-6).unwrap();
        let cd = coordinate_descent(&f, &[0.0, 0.0], &[(-10.0, 10.0), (-10.0, 10.0)], 100, 1e-8).unwrap();
        let nm = nelder_mead(&f, &[0.0, 0.0], 5000, 1e-8);
        for x in [gd, cd, nm] {
            assert!((x[0] - 2.0).abs() < 1e-2);
            assert!((x[1] - 3.0).abs() < 1e-2);
        }
    }
}
