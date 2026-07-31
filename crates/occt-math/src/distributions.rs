//! Probability distributions and distribution sampling helpers.
//! Source: `math_PSO` sampling utilities and standard distributions.

use std::f64::consts::{PI, SQRT_2};

/// Error function. Standard 5-term rational approximation (Abramowitz &
/// Stegun 7.1.26), absolute error below ~1.5e-7. `f64::erf` is not stable
/// on this toolchain, hence the local implementation.
fn erf(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x < 0.0 {
        return -erf(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = t * (0.254_829_592
        + t * (-0.284_496_736 + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    1.0 - y * (-x * x).exp()
}

/// Inverse error function on `(-1, 1)`, via a Winitzki initial guess plus
/// Newton refinement. Returns `±∞` at `±1` and `NaN` outside the domain.
pub fn erf_inv(x: f64) -> f64 {
    if x <= -1.0 || x >= 1.0 {
        return if x == 1.0 {
            f64::INFINITY
        } else if x == -1.0 {
            f64::NEG_INFINITY
        } else {
            f64::NAN
        };
    }
    if x == 0.0 {
        return 0.0;
    }
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let t = x.abs();
    // Winitzki's approximation: accurate to ~1e-4 over the whole domain.
    let a = 0.147;
    let u = (1.0 - t * t).ln();
    let inner = 2.0 / (PI * a) + u / 2.0;
    let mut z = ((inner * inner - u / a).sqrt() - inner).sqrt();
    // Newton: z -= (erf(z) - t) / erf'(z), erf'(z) = 2/√π e^(-z²).
    for _ in 0..12 {
        let e = erf(z) - t;
        let d = (2.0 / PI.sqrt()) * (-(z * z)).exp();
        let dz = e / d;
        z -= dz;
        if dz.abs() < 1e-16 {
            break;
        }
    }
    sign * z
}

/// Gaussian probability density at `x`.
pub fn normal_pdf(x: f64, mean: f64, std: f64) -> f64 {
    let z = (x - mean) / std;
    (-0.5 * z * z).exp() / (std * (2.0 * PI).sqrt())
}

/// Gaussian cumulative distribution at `x`.
pub fn normal_cdf(x: f64, mean: f64, std: f64) -> f64 {
    0.5 * (1.0 + erf((x - mean) / (std * SQRT_2)))
}

/// Mean and sample variance of `samples`. Returns `(0, 0)` when empty and
/// `(mean, 0)` for a single sample (variance is undefined there).
pub fn mean_var(samples: &[f64]) -> (f64, f64) {
    let n = samples.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let m = samples.iter().sum::<f64>() / n as f64;
    if n < 2 {
        return (m, 0.0);
    }
    let v = samples.iter().map(|s| (s - m) * (s - m)).sum::<f64>() / (n - 1) as f64;
    (m, v)
}

/// Bin counts over `[lo, hi)` split into `bins` equal-width buckets. Values
/// outside the range are clamped into the nearest end bin.
pub fn histogram(samples: &[f64], bins: usize, lo: f64, hi: f64) -> Vec<usize> {
    if bins == 0 {
        return Vec::new();
    }
    if samples.is_empty() {
        return vec![0; bins];
    }
    let width = (hi - lo) / bins as f64;
    if width <= 0.0 {
        // Degenerate range: every sample lands in the first bin.
        let mut counts = vec![0usize; bins];
        counts[0] = samples.len();
        return counts;
    }
    let mut counts = vec![0usize; bins];
    for &s in samples {
        let mut idx = ((s - lo) / width) as isize;
        if idx < 0 {
            idx = 0;
        }
        if idx >= bins as isize {
            idx = bins as isize - 1;
        }
        counts[idx as usize] += 1;
    }
    counts
}

/// Standard normal from two uniforms `u1, u2 in (0, 1]` via Box-Muller.
pub fn box_muller(u1: f64, u2: f64) -> f64 {
    let u1 = u1.max(f64::MIN_POSITIVE);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// Gaussian distribution.
pub struct NormalDistribution {
    pub mean: f64,
    pub std: f64,
}

impl NormalDistribution {
    pub fn new(mean: f64, std: f64) -> Result<Self, String> {
        if std <= 0.0 {
            return Err("NormalDistribution::new: std must be > 0".to_string());
        }
        Ok(Self { mean, std })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        normal_pdf(x, self.mean, self.std)
    }

    pub fn cdf(&self, x: f64) -> f64 {
        normal_cdf(x, self.mean, self.std)
    }

    /// Deterministic inverse-CDF sample: `u` is a uniform draw in `[0, 1]`.
    pub fn sample(&self, u: f64) -> f64 {
        self.mean + self.std * SQRT_2 * erf_inv(2.0 * u - 1.0)
    }
}

/// Uniform distribution on `[lo, hi]`.
pub struct UniformDistribution {
    pub lo: f64,
    pub hi: f64,
}

impl UniformDistribution {
    pub fn new(lo: f64, hi: f64) -> Result<Self, String> {
        if hi <= lo {
            return Err("UniformDistribution::new: hi must be > lo".to_string());
        }
        Ok(Self { lo, hi })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        if x >= self.lo && x <= self.hi {
            1.0 / (self.hi - self.lo)
        } else {
            0.0
        }
    }

    pub fn cdf(&self, x: f64) -> f64 {
        ((x - self.lo) / (self.hi - self.lo)).clamp(0.0, 1.0)
    }

    /// Deterministic inverse-CDF sample: `u` in `[0, 1]`.
    pub fn sample(&self, u: f64) -> f64 {
        self.lo + (self.hi - self.lo) * u.clamp(0.0, 1.0)
    }
}

/// Exponential distribution with rate `rate > 0`.
pub struct ExponentialDistribution {
    pub rate: f64,
}

impl ExponentialDistribution {
    pub fn new(rate: f64) -> Result<Self, String> {
        if rate <= 0.0 {
            return Err("ExponentialDistribution::new: rate must be > 0".to_string());
        }
        Ok(Self { rate })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            0.0
        } else {
            self.rate * (-self.rate * x).exp()
        }
    }

    pub fn cdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            0.0
        } else {
            1.0 - (-self.rate * x).exp()
        }
    }

    /// Deterministic inverse-CDF sample: `u` in `[0, 1]`.
    pub fn sample(&self, u: f64) -> f64 {
        -((1.0 - u.clamp(0.0, 1.0)).ln()) / self.rate
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SplitMix64;

    fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
        let n = (n / 2) * 2;
        let h = (b - a) / n as f64;
        let mut s = f(a) + f(b);
        for i in 1..n {
            s += if i % 2 == 1 { 4.0 } else { 2.0 } * f(a + i as f64 * h);
        }
        s * h / 3.0
    }

    #[test]
    fn normal_pdf_integrates_to_one() {
        let d = NormalDistribution::new(0.0, 1.0).unwrap();
        // Wide, mean-centered range: the ±10-sigma tails are ~1e-23.
        let area = simpson(&|x| d.pdf(x), -10.0, 10.0, 40_000);
        assert!((area - 1.0).abs() < 1e-9, "area = {area}");
    }

    #[test]
    fn normal_cdf_mean_is_half() {
        let d = NormalDistribution::new(3.0, 1.5).unwrap();
        assert!((d.cdf(3.0) - 0.5).abs() < 1e-12);
        // Symmetry sanity: cdf(mean+std) = 1 - cdf(mean-std).
        let a = d.cdf(4.5);
        let b = d.cdf(1.5);
        assert!((a + b - 1.0).abs() < 1e-6, "a={a} b={b}");
    }

    #[test]
    fn normal_sample_half_is_mean() {
        let d = NormalDistribution::new(2.0, 0.7).unwrap();
        assert!((d.sample(0.5) - 2.0).abs() < 1e-12);
        // Deterministic: same input, same output.
        assert_eq!(d.sample(0.3), d.sample(0.3));
    }

    #[test]
    fn erf_inv_roundtrip() {
        for x in [-0.9, -0.5, 0.0, 0.3, 0.7, 0.95] {
            let y = erf_inv(x);
            assert!((erf(y) - x).abs() < 1e-6, "x={x} y={y}");
        }
        assert_eq!(erf_inv(0.0), 0.0);
        assert_eq!(erf_inv(1.0), f64::INFINITY);
    }

    #[test]
    fn uniform_roundtrip() {
        let d = UniformDistribution::new(2.0, 5.0).unwrap();
        for u in [0.0, 0.25, 0.5, 0.75, 1.0] {
            assert!((d.cdf(d.sample(u)) - u).abs() < 1e-12, "u = {u}");
        }
        assert_eq!(d.pdf(1.0), 0.0);
        assert!((d.pdf(3.5) - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn exponential_pdf_cdf_sample() {
        let d = ExponentialDistribution::new(0.5).unwrap();
        assert!((d.cdf(1.0) - (1.0 - (-0.5f64).exp())).abs() < 1e-12);
        assert!((d.sample(0.5) - 2.0 * 2.0f64.ln()).abs() < 1e-12);
        assert_eq!(d.pdf(-1.0), 0.0);
        assert!((d.pdf(1.0) - 0.5 * (-0.5f64).exp()).abs() < 1e-12);
    }

    #[test]
    fn mean_var_of_known_data() {
        let (m, v) = mean_var(&[1.0, 2.0, 3.0]);
        assert!((m - 2.0).abs() < 1e-12);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn histogram_counts_sum_to_samples() {
        let samples = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 0.0];
        let counts = histogram(&samples, 5, 0.0, 1.0);
        assert_eq!(counts.len(), 5);
        assert_eq!(counts.iter().sum::<usize>(), samples.len());
        // Degenerate range clamps everything into one bin.
        let counts = histogram(&[1.0, 2.0, 3.0], 4, 2.0, 2.0);
        assert_eq!(counts[0], 3);
    }

    #[test]
    fn box_muller_standard_normal_sanity() {
        let mut rng = SplitMix64::new(99);
        let mut sum = 0.0;
        let mut sq = 0.0;
        let n = 50_000;
        for _ in 0..n {
            let z = box_muller(rng.next_f64(), rng.next_f64());
            sum += z;
            sq += z * z;
        }
        let mean = sum / n as f64;
        let var = sq / n as f64 - mean * mean;
        assert!(mean.abs() < 0.02, "mean = {mean}");
        assert!((var - 1.0).abs() < 0.05, "var = {var}");
    }
}
