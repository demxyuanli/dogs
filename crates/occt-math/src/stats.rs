//! Statistical utilities for one- and two-dimensional data sets.

/// Arithmetic mean of `data`. Returns 0.0 for an empty slice.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    sum(data) / data.len() as f64
}

/// Population variance (divide by n).
pub fn variance(data: &[f64]) -> f64 {
    let n = data.len();
    if n == 0 {
        return 0.0;
    }
    let m = mean(data);
    data.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / n as f64
}

/// Sample variance (divide by n - 1). Returns 0.0 when n < 2.
pub fn sample_variance(data: &[f64]) -> f64 {
    let n = data.len();
    if n < 2 {
        return 0.0;
    }
    let m = mean(data);
    data.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (n - 1) as f64
}

/// Population standard deviation.
pub fn std_dev(data: &[f64]) -> f64 {
    variance(data).sqrt()
}

/// Median of `data` (averages the two middle values for even lengths).
pub fn median(data: &[f64]) -> f64 {
    let n = data.len();
    if n == 0 {
        return 0.0;
    }
    let mut s: Vec<f64> = data.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

/// Minimum and maximum of `data`. Returns (NaN, NaN) for an empty slice.
pub fn min_max(data: &[f64]) -> (f64, f64) {
    if data.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let mut min = data[0];
    let mut max = data[0];
    for &v in &data[1..] {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    (min, max)
}

/// Sum of all values.
pub fn sum(data: &[f64]) -> f64 {
    data.iter().sum()
}

/// Quantile via linear interpolation; `sorted` is assumed ascending.
/// `q` is clamped to [0, 1].
pub fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let q = q.clamp(0.0, 1.0);
    let pos = q * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        return sorted[lo];
    }
    let frac = pos - lo as f64;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

/// Sample covariance; `None` when the slices differ in length or are empty.
pub fn covariance(x: &[f64], y: &[f64]) -> Option<f64> {
    let n = x.len();
    if n == 0 || n != y.len() {
        return None;
    }
    let mx = mean(x);
    let my = mean(y);
    Some(x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum::<f64>() / n as f64)
}

/// Pearson correlation coefficient; `None` when either series has zero spread.
pub fn correlation(x: &[f64], y: &[f64]) -> Option<f64> {
    let cov = covariance(x, y)?;
    let sx = std_dev(x);
    let sy = std_dev(y);
    if sx == 0.0 || sy == 0.0 {
        return None;
    }
    Some(cov / (sx * sy))
}

/// Histogram of `data`. Returns `(bin_edges, counts)` where `bin_edges`
/// has `bins + 1` entries spanning [min, max].
pub fn histogram(data: &[f64], bins: usize) -> (Vec<f64>, Vec<usize>) {
    let bins = bins.max(1);
    if data.is_empty() {
        return (vec![0.0; bins + 1], vec![0; bins]);
    }
    let (mut lo, mut hi) = min_max(data);
    if (hi - lo).abs() < f64::EPSILON {
        // Degenerate range: widen it so every value lands in a single bin.
        let pad = hi.abs().max(1.0) * 0.5;
        lo -= pad;
        hi += pad;
    }
    let width = (hi - lo) / bins as f64;
    let edges: Vec<f64> = (0..=bins).map(|i| lo + i as f64 * width).collect();
    let mut counts = vec![0usize; bins];
    for &v in data {
        let mut idx = ((v - lo) / width) as usize;
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1;
    }
    (edges, counts)
}

/// Z-score of `value` given the distribution's `mean` and `sd`.
pub fn z_score(value: f64, mean: f64, sd: f64) -> f64 {
    (value - mean) / sd
}

/// Percentile rank (0-100): the percentage of values strictly below `value`.
pub fn percentile_rank(data: &[f64], value: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let below = data.iter().filter(|&&v| v < value).count();
    below as f64 / data.len() as f64 * 100.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_and_median() {
        let d = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((mean(&d) - 3.0).abs() < 1e-12);
        assert!((median(&d) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn population_variance() {
        let d = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        assert!((variance(&d) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn covariance_and_correlation() {
        let x = [1.0, 2.0, 3.0];
        let y = [2.0, 4.0, 6.0];
        assert!((covariance(&x, &y).unwrap() - 4.0 / 3.0).abs() < 1e-12);
        assert!((correlation(&x, &y).unwrap() - 1.0).abs() < 1e-12);
    }
}
