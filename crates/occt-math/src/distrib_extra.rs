//! Additional probability distributions — chi-square, Student-t, F, gamma,
//! beta and Poisson. Source: standard statistical functions
//! (`math_Recipes`-style algorithms).
//!
//! The distributions are built on numerically stable special functions: a
//! Lanczos/Numerical-Recipes log-gamma, the regularized incomplete gamma via
//! series/continued fraction, and the regularized incomplete beta via the
//! continued-fraction representation.

use std::f64::consts::{PI, SQRT_2};

// ---------------------------------------------------------------------------
// Special functions
// ---------------------------------------------------------------------------

/// Log-gamma `ln Γ(x)` for `x > 0` (Numerical Recipes `gammln`). For `x ≤ 0`
/// the reflection formula is used, which is valid away from the poles.
fn lgamma(x: f64) -> f64 {
    if x < 0.5 {
        let pi = PI;
        // Γ(x)Γ(1−x) = π / sin(πx)
        (pi / (pi * x).sin()).ln() - lgamma(1.0 - x)
    } else {
        let cof = [
            76.180_091_729_471_46,
            -86.505_320_329_416_77,
            24.014_098_240_830_91,
            -1.231_739_572_450_155,
            0.120_865_097_386_617_9e-2,
            -0.539_523_938_495_3e-5,
        ];
        let mut y = x;
        let tmp = x + 5.5;
        let tmp = tmp - (x + 0.5) * tmp.ln();
        let mut ser = 1.000_000_000_190_015;
        for j in 0..6 {
            y += 1.0;
            ser += cof[j] / y;
        }
        -tmp + (2.506_628_274_631_000_5 * ser / x).ln()
    }
}

/// Complementary error function `erfc(x)` (Numerical Recipes approximation,
/// accurate to ~1.2e-7).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.265_512_23
            + t
                * (1.000_023_68
                    + t
                        * (0.374_091_96
                            + t
                                * (0.096_784_18
                                    + t
                                        * (-0.186_288_06
                                            + t
                                                * (0.278_868_07
                                                    + t
                                                        * (-1.135_203_98
                                                            + t
                                                                * (1.488_515_87
                                                                    + t
                                                                        * (-0.822_152_23
                                                                            + t * 0.170_872_77)))))))))
        .exp();
    ans
}

/// Error function `erf(x)`.
fn erf(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x >= 0.0 {
        1.0 - erfc(x)
    } else {
        erfc(-x) - 1.0
    }
}

/// Cumulative distribution function of the standard normal distribution.
pub fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / SQRT_2))
}

/// Regularized incomplete gamma `P(a, x) = γ(a, x)/Γ(a)`.
///
/// Uses the series representation for `x < a + 1` and the continued-fraction
/// representation of the complement `Q(a, x) = 1 − P(a, x)` otherwise.
fn regularized_gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if a <= 0.0 {
        // Limit: P(0⁺, x) = 1 for x > 0.
        return 1.0;
    }
    if x < a + 1.0 {
        // Series for P(a, x).
        let gln = lgamma(a);
        let mut ap = a;
        let mut sum = 1.0 / a;
        let mut del = sum;
        for _ in 0..200 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-14 {
                break;
            }
        }
        (sum * (-x + a * x.ln() - gln).exp()).clamp(0.0, 1.0)
    } else {
        // Continued fraction for Q(a, x).
        (1.0 - gamma_q_cf(a, x)).clamp(0.0, 1.0)
    }
}

/// Continued-fraction evaluation of `Q(a, x) = 1 − P(a, x)`.
fn gamma_q_cf(a: f64, x: f64) -> f64 {
    let gln = lgamma(a);
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / f64::MIN_POSITIVE;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..200 {
        let an = -i as f64 * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < f64::MIN_POSITIVE {
            d = f64::MIN_POSITIVE;
        }
        c = b + an / c;
        if c.abs() < f64::MIN_POSITIVE {
            c = f64::MIN_POSITIVE;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-14 {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

/// Regularized incomplete beta `I_x(a, b)`.
fn regularized_beta_inc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        (bt * beta_cf(a, b, x) / a).clamp(0.0, 1.0)
    } else {
        (1.0 - bt * beta_cf(b, a, 1.0 - x) / b).clamp(0.0, 1.0)
    }
}

/// Continued fraction used by [`regularized_beta_inc`] (Numerical Recipes
/// `betacf`).
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    let fpmin = f64::MIN_POSITIVE;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < fpmin {
        d = fpmin;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..200 {
        let m2 = 2 * m;
        // Even step.
        let mut aa = m as f64 * (b - m as f64) * x / ((qam + m2 as f64) * (a + m2 as f64));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        h *= d * c;
        // Odd step.
        aa = -(a + m as f64) * (qab + m as f64) * x / ((a + m2 as f64) * (qap + m2 as f64));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-14 {
            break;
        }
    }
    h
}

// ---------------------------------------------------------------------------
// Distributions
// ---------------------------------------------------------------------------

/// Chi-square distribution with `k` degrees of freedom (a gamma distribution
/// with shape `k/2` and scale 2).
pub struct ChiSquareDistribution {
    pub k: f64,
}

impl ChiSquareDistribution {
    pub fn new(k: f64) -> Result<Self, String> {
        if k <= 0.0 {
            return Err("ChiSquareDistribution::new: k must be > 0".to_string());
        }
        Ok(Self { k })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        gamma_pdf_impl(x, self.k / 2.0, 2.0)
    }

    pub fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            0.0
        } else {
            regularized_gamma_p(self.k / 2.0, x / 2.0)
        }
    }

    pub fn mean(&self) -> f64 {
        self.k
    }

    pub fn variance(&self) -> f64 {
        2.0 * self.k
    }
}

/// Student-t distribution with `nu` degrees of freedom.
pub struct StudentTDistribution {
    pub nu: f64,
}

impl StudentTDistribution {
    pub fn new(nu: f64) -> Result<Self, String> {
        if nu <= 0.0 {
            return Err("StudentTDistribution::new: nu must be > 0".to_string());
        }
        Ok(Self { nu })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        let nu = self.nu;
        let logp = lgamma((nu + 1.0) / 2.0) - lgamma(nu / 2.0) - 0.5 * (nu * PI).ln()
            - ((nu + 1.0) / 2.0) * (1.0 + x * x / nu).ln();
        logp.exp()
    }

    pub fn cdf(&self, x: f64) -> f64 {
        let nu = self.nu;
        let y = nu / (nu + x * x);
        if x >= 0.0 {
            1.0 - 0.5 * regularized_beta_inc(nu / 2.0, 0.5, y)
        } else {
            0.5 * regularized_beta_inc(nu / 2.0, 0.5, y)
        }
    }

    pub fn mean(&self) -> f64 {
        0.0
    }

    pub fn variance(&self) -> f64 {
        if self.nu > 2.0 {
            self.nu / (self.nu - 2.0)
        } else {
            f64::NAN
        }
    }
}

/// F-distribution with numerator/denominator degrees of freedom `d1`, `d2`.
pub struct FDistribution {
    pub d1: f64,
    pub d2: f64,
}

impl FDistribution {
    pub fn new(d1: f64, d2: f64) -> Result<Self, String> {
        if d1 <= 0.0 || d2 <= 0.0 {
            return Err("FDistribution::new: d1 and d2 must be > 0".to_string());
        }
        Ok(Self { d1, d2 })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let (d1, d2) = (self.d1, self.d2);
        let logp = lgamma((d1 + d2) / 2.0) - lgamma(d1 / 2.0) - lgamma(d2 / 2.0)
            + (d1 / 2.0) * (d1 / d2).ln()
            + (d1 / 2.0 - 1.0) * x.ln()
            - ((d1 + d2) / 2.0) * (1.0 + d1 * x / d2).ln();
        logp.exp()
    }

    pub fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            0.0
        } else {
            let (d1, d2) = (self.d1, self.d2);
            regularized_beta_inc(d1 / 2.0, d2 / 2.0, d1 * x / (d1 * x + d2))
        }
    }

    pub fn mean(&self) -> f64 {
        if self.d2 > 2.0 {
            self.d2 / (self.d2 - 2.0)
        } else {
            f64::NAN
        }
    }

    pub fn variance(&self) -> f64 {
        let (d1, d2) = (self.d1, self.d2);
        if d2 > 4.0 {
            2.0 * d2 * d2 * (d1 + d2 - 2.0) / (d1 * (d2 - 2.0) * (d2 - 2.0) * (d2 - 4.0))
        } else {
            f64::NAN
        }
    }
}

/// Gamma distribution with shape `k` and scale `theta`.
pub struct GammaDistribution {
    pub k: f64,
    pub theta: f64,
}

impl GammaDistribution {
    pub fn new(k: f64, theta: f64) -> Result<Self, String> {
        if k <= 0.0 || theta <= 0.0 {
            return Err("GammaDistribution::new: k and theta must be > 0".to_string());
        }
        Ok(Self { k, theta })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        gamma_pdf_impl(x, self.k, self.theta)
    }

    pub fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            0.0
        } else {
            regularized_gamma_p(self.k, x / self.theta)
        }
    }

    pub fn mean(&self) -> f64 {
        self.k * self.theta
    }

    pub fn variance(&self) -> f64 {
        self.k * self.theta * self.theta
    }
}

/// Beta distribution on `[0, 1]` with shape parameters `a`, `b`.
pub struct BetaDistribution {
    pub a: f64,
    pub b: f64,
}

impl BetaDistribution {
    pub fn new(a: f64, b: f64) -> Result<Self, String> {
        if a <= 0.0 || b <= 0.0 {
            return Err("BetaDistribution::new: a and b must be > 0".to_string());
        }
        Ok(Self { a, b })
    }

    pub fn pdf(&self, x: f64) -> f64 {
        beta_pdf_impl(x, self.a, self.b)
    }

    pub fn cdf(&self, x: f64) -> f64 {
        regularized_beta_inc(self.a, self.b, x)
    }

    pub fn mean(&self) -> f64 {
        self.a / (self.a + self.b)
    }

    pub fn variance(&self) -> f64 {
        self.a * self.b / ((self.a + self.b) * (self.a + self.b) * (self.a + self.b + 1.0))
    }
}

/// Poisson distribution with mean (rate) `lambda`.
pub struct PoissonDistribution {
    pub lambda: f64,
}

impl PoissonDistribution {
    pub fn new(lambda: f64) -> Result<Self, String> {
        if lambda <= 0.0 {
            return Err("PoissonDistribution::new: lambda must be > 0".to_string());
        }
        Ok(Self { lambda })
    }

    pub fn pmf(&self, k: i64) -> f64 {
        poisson_pmf_impl(k, self.lambda)
    }

    pub fn cdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            return 0.0;
        }
        let kmax = x.floor() as i64;
        let mut sum = 0.0;
        for k in 0..=kmax {
            sum += self.pmf(k);
        }
        sum.min(1.0)
    }

    pub fn mean(&self) -> f64 {
        self.lambda
    }

    pub fn variance(&self) -> f64 {
        self.lambda
    }
}

// ---------------------------------------------------------------------------
// Internal pdf helpers (log-space for stability)
// ---------------------------------------------------------------------------

/// Gamma pdf `x^(k−1) e^(−x/θ) / (Γ(k) θ^k)`.
fn gamma_pdf_impl(x: f64, k: f64, theta: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return if k > 1.0 {
            0.0
        } else if k == 1.0 {
            1.0 / theta
        } else {
            f64::INFINITY
        };
    }
    let logp = (k - 1.0) * x.ln() - x / theta - lgamma(k) - k * theta.ln();
    logp.exp()
}

/// Beta pdf `x^(a−1) (1−x)^(b−1) / B(a, b)`.
fn beta_pdf_impl(x: f64, a: f64, b: f64) -> f64 {
    if x < 0.0 || x > 1.0 {
        return 0.0;
    }
    if x == 0.0 {
        return if a > 1.0 {
            0.0
        } else if a == 1.0 {
            b
        } else {
            f64::INFINITY
        };
    }
    if x == 1.0 {
        return if b > 1.0 {
            0.0
        } else if b == 1.0 {
            a
        } else {
            f64::INFINITY
        };
    }
    let logp = (a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() - (lgamma(a) + lgamma(b) - lgamma(a + b));
    logp.exp()
}

/// Poisson pmf `e^(−λ) λ^k / k!`.
fn poisson_pmf_impl(k: i64, lambda: f64) -> f64 {
    if k < 0 {
        return 0.0;
    }
    let logp = -lambda + k as f64 * lambda.ln() - lgamma(k as f64 + 1.0);
    logp.exp()
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn normal_cdf_known_values() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 1e-6);
        assert!((normal_cdf(1.96) - 0.975).abs() < 1e-4, "{}", normal_cdf(1.96));
        assert!((normal_cdf(-1.96) - 0.025).abs() < 1e-4);
        assert!(normal_cdf(8.0) > 1.0 - 1e-12);
    }

    #[test]
    fn chi_square_known_values() {
        let d = ChiSquareDistribution::new(1.0).unwrap();
        // χ²_1 cdf(x) = erf(√(x/2)).
        assert!((d.cdf(3.84) - 0.95).abs() < 1e-3, "cdf(3.84) = {}", d.cdf(3.84));
        assert!((d.cdf(1.0) - 0.6827).abs() < 1e-3, "cdf(1) = {}", d.cdf(1.0));
        assert!((d.mean() - 1.0).abs() < 1e-12);
        assert!((d.variance() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn chi_square_two_equals_exponential_half() {
        // χ²_2 cdf = 1 − e^(−x/2), the same as Exp(rate 0.5).
        let chi = ChiSquareDistribution::new(2.0).unwrap();
        let exp = crate::distributions::ExponentialDistribution::new(0.5).unwrap();
        for x in [0.1, 0.5, 1.0, 2.0, 5.0] {
            assert!((chi.cdf(x) - exp.cdf(x)).abs() < 1e-10, "x={x}");
        }
    }

    #[test]
    fn student_t_cdf_half_at_zero() {
        for nu in [1.0, 2.0, 5.0, 30.0] {
            let d = StudentTDistribution::new(nu).unwrap();
            assert!((d.cdf(0.0) - 0.5).abs() < 1e-10, "nu={nu}");
        }
    }

    #[test]
    fn student_t_one_is_cauchy() {
        let d = StudentTDistribution::new(1.0).unwrap();
        assert!((d.cdf(1.0) - 0.75).abs() < 1e-6, "cdf(1) = {}", d.cdf(1.0));
        // Cauchy cdf(x) = 0.5 + arctan(x)/π.
        assert!((d.cdf(2.0) - (0.5 + 2.0f64.atan() / PI)).abs() < 1e-6, "cdf(2) = {}", d.cdf(2.0));
    }

    #[test]
    fn f_distribution_known() {
        let d = FDistribution::new(1.0, 1.0).unwrap();
        // I_{x/(1+x)}(1/2, 1/2) at x=1 → I_{1/2}(1/2, 1/2) = 0.5.
        assert!((d.cdf(1.0) - 0.5).abs() < 1e-6, "F(1,1) cdf(1) = {}", d.cdf(1.0));
        let d = FDistribution::new(5.0, 10.0).unwrap();
        assert!((d.mean() - 10.0 / 8.0).abs() < 1e-10);
        assert!(d.pdf(1.0) > 0.0);
    }

    #[test]
    fn gamma_pdf_integrates_to_one() {
        let d = GammaDistribution::new(2.0, 1.5).unwrap();
        let area = simpson(&|x| d.pdf(x), 0.0, 30.0, 40_000);
        assert!((area - 1.0).abs() < 1e-6, "gamma area = {area}");
        assert!((d.mean() - 3.0).abs() < 1e-10);
        assert!((d.variance() - 4.5).abs() < 1e-9);
    }

    #[test]
    fn beta_pdf_integrates_and_symmetric() {
        let d = BetaDistribution::new(2.0, 2.0).unwrap();
        let area = simpson(&|x| d.pdf(x), 0.0, 1.0, 20_000);
        assert!((area - 1.0).abs() < 1e-6, "beta area = {area}");
        // Beta(2,2) is symmetric about 0.5.
        assert!((d.cdf(0.5) - 0.5).abs() < 1e-8, "beta cdf(0.5) = {}", d.cdf(0.5));
        assert!((d.mean() - 0.5).abs() < 1e-12);
        assert!(d.pdf(0.0) == 0.0);
        assert!(d.pdf(1.0) == 0.0);
    }

    #[test]
    fn poisson_pmf_sums_to_one() {
        let d = PoissonDistribution::new(4.0).unwrap();
        let sum: f64 = (0..=40).map(|k| d.pmf(k)).sum();
        assert!((sum - 1.0).abs() < 1e-10, "pmf sum = {sum}");
        assert!((d.mean() - 4.0).abs() < 1e-12);
        assert!((d.variance() - 4.0).abs() < 1e-12);
        // cdf is monotone and reaches 1.
        let mut last = -1.0;
        for x in (0..=20).map(|i| i as f64) {
            let c = d.cdf(x);
            assert!(c >= last, "cdf not monotone at {x}");
            last = c;
        }
        assert!((d.cdf(20.0) - 1.0).abs() < 1e-8);
    }

    #[test]
    fn invalid_params_error() {
        assert!(ChiSquareDistribution::new(0.0).is_err());
        assert!(StudentTDistribution::new(-1.0).is_err());
        assert!(FDistribution::new(0.0, 1.0).is_err());
        assert!(FDistribution::new(1.0, 0.0).is_err());
        assert!(GammaDistribution::new(0.0, 1.0).is_err());
        assert!(GammaDistribution::new(1.0, -2.0).is_err());
        assert!(BetaDistribution::new(0.0, 1.0).is_err());
        assert!(BetaDistribution::new(1.0, 0.0).is_err());
        assert!(PoissonDistribution::new(0.0).is_err());
    }

    #[test]
    fn valid_params_ok() {
        assert!(ChiSquareDistribution::new(3.0).is_ok());
        assert!(StudentTDistribution::new(5.0).is_ok());
        assert!(FDistribution::new(2.0, 3.0).is_ok());
        assert!(GammaDistribution::new(1.0, 1.0).is_ok());
        assert!(BetaDistribution::new(1.0, 1.0).is_ok());
        assert!(PoissonDistribution::new(2.5).is_ok());
    }

    #[test]
    fn cdf_monotone_and_pdf_nonnegative() {
        let dists: Vec<Box<dyn Fn(f64) -> f64>> = vec![
            Box::new(|x| ChiSquareDistribution::new(4.0).unwrap().cdf(x)),
            Box::new(|x| GammaDistribution::new(2.0, 1.0).unwrap().cdf(x)),
            Box::new(|x| FDistribution::new(3.0, 5.0).unwrap().cdf(x)),
            Box::new(|x| BetaDistribution::new(2.0, 3.0).unwrap().cdf(x)),
        ];
        for cdf in &dists {
            let mut last = -1.0;
            for i in 0..100 {
                let x = i as f64 / 100.0;
                let c = cdf(x);
                assert!(c >= last, "cdf decreased at {x}");
                last = c;
                assert!((0.0..=1.0).contains(&c));
            }
        }
        // pdfs non-negative.
        let g = GammaDistribution::new(3.0, 2.0).unwrap();
        let b = BetaDistribution::new(2.0, 5.0).unwrap();
        let t = StudentTDistribution::new(4.0).unwrap();
        let f = FDistribution::new(2.0, 4.0).unwrap();
        for i in 0..200 {
            let x = -0.5 + 1.5 * i as f64 / 199.0;
            assert!(g.pdf(x) >= 0.0);
            assert!(t.pdf(x) >= 0.0);
            let xb = i as f64 / 199.0;
            assert!(b.pdf(xb) >= 0.0);
            let xf = 0.01 + 5.0 * i as f64 / 199.0;
            assert!(f.pdf(xf) >= 0.0);
        }
    }

    #[test]
    fn student_t_variance_matches_formula() {
        let d = StudentTDistribution::new(5.0).unwrap();
        assert!((d.variance() - 5.0 / 3.0).abs() < 1e-12);
        let d = StudentTDistribution::new(2.0).unwrap();
        assert!(d.variance().is_nan()); // undefined for ν ≤ 2
    }
}
