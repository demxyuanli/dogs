//! Gauss integration points and weights. Source: `math_GaussSingleIntegration.cxx`
//! Generates Gauss-Legendre quadrature points for numerical integration.

/// Gauss-Legendre quadrature rule for interval [a, b].
/// Returns (points, weights) arrays of length n.
pub fn gauss_legendre(a: f64, b: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0f64; n];
    let mut w = vec![0.0f64; n];
    let m = (n + 1) / 2;
    let xm = 0.5 * (b + a);
    let xl = 0.5 * (b - a);

    for i in 1..=m {
        let mut z = (std::f64::consts::PI * (i as f64 - 0.25) / (n as f64 + 0.5)).cos();
        let mut pp = 0.0;
        loop {
            let mut p1 = 1.0;
            let mut p2 = 0.0;
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j as f64 - 1.0) * z * p2 - (j as f64 - 1.0) * p3) / j as f64;
            }
            pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / pp;
            if (z - z1).abs() < 1e-15 { break; }
        }
        x[i - 1] = xm - xl * z;
        x[n - i] = xm + xl * z;
        w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - i] = w[i - 1];
    }
    (x, w)
}

/// Kronrod extension for error estimation.
/// Returns Gauss-Kronrod (n=15) points and weights.
pub fn gauss_kronrod15(a: f64, b: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // GK15 nodes and weights for [-1, 1]
    let nodes: [f64; 15] = [
        -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993945,
        -0.5860872354676911, -0.4058451513773972, -0.20778495500789848, 0.0,
        0.20778495500789848, 0.4058451513773972, 0.5860872354676911, 0.7415311855993945,
        0.8648644233597691, 0.9491079123427585, 0.9914553711208126,
    ];
    let gauss_w: [f64; 7] = [
        0.1294849661688697, 0.27970539148927664, 0.3818300505051189, 0.4179591836734694,
        0.3818300505051189, 0.27970539148927664, 0.1294849661688697,
    ];
    let kronrod_w: [f64; 15] = [
        0.022935322010529224, 0.06309209262997856, 0.10479001032225019, 0.14065325971552592,
        0.1690047266392679, 0.19035057806478542, 0.2044329400752989, 0.20948214108472782,
        0.2044329400752989, 0.19035057806478542, 0.1690047266392679, 0.14065325971552592,
        0.10479001032225019, 0.06309209262997856, 0.022935322010529224,
    ];

    let xm = 0.5 * (b + a);
    let xl = 0.5 * (b - a);
    let pts: Vec<f64> = nodes.iter().map(|&n| xm + xl * n).collect();
    let w_gauss: Vec<f64> = (0..7).map(|i| {
        let j = 2 * i + 1;
        xl * gauss_w[i]
    }).collect();
    let w_kronrod: Vec<f64> = kronrod_w.iter().map(|&w| xl * w).collect();
    (pts, w_gauss, w_kronrod)
}

/// Integrate f over [a, b] using n-point Gauss-Legendre quadrature.
pub fn integrate<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    let (pts, w) = gauss_legendre(a, b, n);
    pts.iter().zip(w.iter()).map(|(&x, &w)| f(x) * w).sum()
}

/// Integrate f over [a, b] with error estimate using Gauss-Kronrod (G7/K15).
pub fn integrate_with_error<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64) -> (f64, f64) {
    let (pts, w_gauss, w_kronrod) = gauss_kronrod15(a, b);
    let g7: f64 = pts.iter().enumerate()
        .filter(|(i, _)| i % 2 == 1)
        .map(|(i, &x)| f(x) * w_gauss[i / 2]).sum();
    let k15: f64 = pts.iter().enumerate()
        .map(|(i, &x)| f(x) * w_kronrod[i]).sum();
    (g7, (200.0 * (g7 - k15).abs()).cbrt()) // crude error estimate
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrate_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let v = integrate(&|x| x * x, 0.0, 1.0, 5);
        assert!((v - 1.0 / 3.0).abs() < 1e-14);
    }

    #[test]
    fn integrate_sin() {
        // ∫₀^π sin(x) dx = 2
        let v = integrate(&|x| x.sin(), 0.0, std::f64::consts::PI, 10);
        assert!((v - 2.0).abs() < 1e-14);
    }
}
