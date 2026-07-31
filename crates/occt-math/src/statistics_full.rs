//! Extended statistics — moments, percentiles, covariance, correlation and
//! regression. Source: `math_Recipes`, `math_Statistics`.

use std::collections::HashMap;

/// Mean of a sample.
pub fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().sum::<f64>() / xs.len() as f64
}

/// Sample variance (n−1 denominator). Returns 0 for < 2 samples.
pub fn variance(xs: &[f64]) -> f64 {
    if xs.len() < 2 {
        return 0.0;
    }
    let m = mean(xs);
    xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (xs.len() - 1) as f64
}

/// Standard deviation.
pub fn std_dev(xs: &[f64]) -> f64 {
    variance(xs).sqrt()
}

/// Central moments up to order 4 (μ2=var, μ3, μ4).
pub fn central_moments(xs: &[f64]) -> [f64; 5] {
    let m = mean(xs);
    let mut acc = [0.0; 5];
    for &x in xs {
        let d = x - m;
        acc[1] += d;
        acc[2] += d * d;
        acc[3] += d * d * d;
        acc[4] += d * d * d * d;
    }
    let n = xs.len().max(1) as f64;
    for i in 1..5 {
        acc[i] /= n;
    }
    acc
}

/// Skewness (Fisher–Pearson). 0 for symmetric distributions.
pub fn skewness(xs: &[f64]) -> f64 {
    if xs.len() < 3 {
        return 0.0;
    }
    let m = central_moments(xs);
    let s = m[2].sqrt();
    if s.abs() < 1e-30 {
        return 0.0;
    }
    m[3] / (s * s * s)
}

/// Excess kurtosis (0 for a normal distribution).
pub fn kurtosis(xs: &[f64]) -> f64 {
    if xs.len() < 4 {
        return 0.0;
    }
    let m = central_moments(xs);
    let v = m[2];
    if v.abs() < 1e-30 {
        return 0.0;
    }
    m[4] / (v * v) - 3.0
}

/// Quantile via the sorted sample (linear interpolation between order
/// statistics). `q` in [0,1].
pub fn quantile(xs: &[f64], q: f64) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mut s = xs.to_vec();
    s.sort_by(f64::total_cmp);
    if q <= 0.0 {
        return s[0];
    }
    if q >= 1.0 {
        return *s.last().unwrap();
    }
    let pos = q * (s.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    s[lo] + frac * (s[hi] - s[lo])
}

/// Median, first and third quartiles.
pub fn quartiles(xs: &[f64]) -> (f64, f64, f64) {
    (quantile(xs, 0.25), quantile(xs, 0.5), quantile(xs, 0.75))
}

/// Mode (most frequent value, first encountered on ties). Uses rounded keys.
pub fn mode(xs: &[f64], round_to: i64) -> Option<f64> {
    let mut counts: HashMap<i64, (f64, usize)> = HashMap::new();
    for &x in xs {
        let key = (x * round_to as f64).round() as i64;
        let e = counts.entry(key).or_insert((x, 0));
        e.1 += 1;
    }
    counts.into_values().max_by(|a, b| a.1.cmp(&b.1)).map(|(v, _)| v)
}

/// Covariance of two paired samples.
pub fn covariance(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().min(ys.len());
    if n < 2 {
        return 0.0;
    }
    let mx = mean(&xs[..n]);
    let my = mean(&ys[..n]);
    xs[..n].iter().zip(&ys[..n]).map(|(&x, &y)| (x - mx) * (y - my)).sum::<f64>() / (n - 1) as f64
}

/// Pearson correlation coefficient in [-1, 1].
pub fn correlation(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().min(ys.len());
    if n < 2 {
        return 0.0;
    }
    let mx = mean(&xs[..n]);
    let my = mean(&ys[..n]);
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..n {
        let dx = xs[i] - mx;
        let dy = ys[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    let denom = (sxx * syy).sqrt();
    if denom < 1e-30 {
        0.0
    } else {
        sxy / denom
    }
}

/// Simple linear regression y = a + b·x. Returns (slope b, intercept a).
pub fn linear_regression(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len().min(ys.len());
    if n < 2 {
        return (0.0, 0.0);
    }
    let mx = mean(&xs[..n]);
    let my = mean(&ys[..n]);
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for i in 0..n {
        let dx = xs[i] - mx;
        sxx += dx * dx;
        sxy += dx * (ys[i] - my);
    }
    if sxx.abs() < 1e-30 {
        return (0.0, my);
    }
    let b = sxy / sxx;
    (b, my - b * mx)
}

/// R² of a linear fit (coefficient of determination).
pub fn r_squared(xs: &[f64], ys: &[f64]) -> f64 {
    let c = correlation(xs, ys);
    c * c
}

/// Frequency counts over integer bins (rounded values).
pub fn frequency_counts(xs: &[f64], bin_width: f64) -> Vec<(f64, usize)> {
    let mut map: HashMap<i64, (f64, usize)> = HashMap::new();
    for &x in xs {
        let key = (x / bin_width).floor() as i64;
        let base = key as f64 * bin_width;
        let e = map.entry(key).or_insert((base, 0));
        e.1 += 1;
    }
    let mut out: Vec<(f64, usize)> = map.into_values().collect();
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_variance() {
        assert!((mean(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-12);
        assert!((variance(&[1.0, 2.0, 3.0]) - 1.0).abs() < 1e-12);
        assert!((std_dev(&[0.0, 0.0, 0.0])).abs() < 1e-12);
    }

    #[test]
    fn moments_skew_kurtosis() {
        // Normal-ish samples → skew ≈ 0, kurtosis ≈ 0.
        let xs: Vec<f64> = (0..40).map(|i| (i as f64 / 40.0 - 0.5)).collect();
        assert!(skewness(&xs).abs() < 0.1);
        // Symmetric uniform → kurtosis ≈ −1.2 (excess).
        let k = kurtosis(&xs);
        assert!(k.abs() < 1.5);
        // Right-skewed: squares.
        let sq: Vec<f64> = (1..=10).map(|i| (i * i) as f64).collect();
        assert!(skewness(&sq) > 0.0);
    }

    #[test]
    fn quantiles() {
        let xs = vec![1.0, 2.0, 3.0, 4.0];
        assert!((quantile(&xs, 0.0) - 1.0).abs() < 1e-12);
        assert!((quantile(&xs, 0.5) - 2.5).abs() < 1e-12);
        assert!((quantile(&xs, 1.0) - 4.0).abs() < 1e-12);
        let (q1, med, q3) = quartiles(&xs);
        assert!(q1 < med && med < q3);
    }

    #[test]
    fn correlation_and_regression() {
        let xs: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
        assert!((correlation(&xs, &ys) - 1.0).abs() < 1e-9);
        assert!((r_squared(&xs, &ys) - 1.0).abs() < 1e-9);
        let (b, a) = linear_regression(&xs, &ys);
        assert!((b - 2.0).abs() < 1e-9);
        assert!((a - 1.0).abs() < 1e-9);
    }

    #[test]
    fn covariance_and_mode() {
        let xs = vec![1.0, 2.0, 3.0, 4.0];
        let ys = vec![2.0, 4.0, 6.0, 8.0];
        assert!(covariance(&xs, &ys) > 0.0);
        let m = mode(&[1.0, 2.0, 2.0, 3.0], 1).expect("mode");
        assert!((m - 2.0).abs() < 1e-9);
    }

    #[test]
    fn frequency_counts_works() {
        let counts = frequency_counts(&[0.1, 0.4, 1.2, 1.8, 2.0, 2.9], 1.0);
        let total: usize = counts.iter().map(|(_, c)| c).sum();
        assert_eq!(total, 6);
    }
}
