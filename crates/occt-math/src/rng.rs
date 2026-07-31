//! Deterministic pseudo-random number generators and sampling helpers.

/// xorshift64* PRNG.
pub struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                seed
            },
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in [lo, hi).
    pub fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    /// Uniform in [0, n).
    pub fn next_usize(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        (self.next_u64() % n as u64) as usize
    }
}

/// 32-bit linear congruential generator.
pub struct Lcg32 {
    state: u32,
}

impl Lcg32 {
    pub fn new(seed: u32) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    pub fn next_u32(&mut self) -> u32 {
        self.state = self.state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        self.state
    }

    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        self.next_u32() as f64 / (u32::MAX as f64 + 1.0)
    }
}

/// SplitMix64 PRNG (tiny state, good bit mixing).
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// PCG-XSH-RR 32-bit generator (O'Neill's PCG).
pub struct Pcg32 {
    state: u64,
    inc: u64,
}

impl Pcg32 {
    pub fn new(seed: u64, stream: u64) -> Self {
        let mut pcg = Self {
            state: 0,
            inc: (stream << 1) | 1,
        };
        pcg.next_u32();
        pcg.state = pcg.state.wrapping_add(seed);
        pcg.next_u32();
        pcg
    }

    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(self.inc);
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        self.next_u32() as f64 / (u32::MAX as f64 + 1.0)
    }
}

/// In-place Fisher-Yates shuffle. `rng(i)` returns an index in `[0, i)`.
pub fn shuffle<T>(items: &mut [T], rng: &mut dyn FnMut(usize) -> usize) {
    for i in (1..items.len()).rev() {
        let j = rng(i + 1).min(i);
        items.swap(i, j);
    }
}

/// Draw a normal(mean, std) sample via Box-Muller.
pub fn normal_sample(mean: f64, std: f64, rng: &mut SplitMix64) -> f64 {
    let u1 = rng.next_f64().max(f64::MIN_POSITIVE);
    let u2 = rng.next_f64();
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
    mean + std * z
}

/// Uniform direction on the unit sphere (rejection sampling on a cube).
pub fn random_unit_vector(rng: &mut SplitMix64) -> [f64; 3] {
    loop {
        let x = rng.next_f64() * 2.0 - 1.0;
        let y = rng.next_f64() * 2.0 - 1.0;
        let z = rng.next_f64() * 2.0 - 1.0;
        let sq = x * x + y * y + z * z;
        if sq > 1e-12 && sq <= 1.0 {
            let inv = 1.0 / sq.sqrt();
            return [x * inv, y * inv, z * inv];
        }
    }
}

/// `n` independent normal(mean, std) draws.
pub fn gaussian_noise(mean: f64, std: f64, n: usize, rng: &mut SplitMix64) -> Vec<f64> {
    (0..n).map(|_| normal_sample(mean, std, rng)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_sequence() {
        let mut a = SplitMix64::new(12_345);
        let mut b = SplitMix64::new(12_345);
        for _ in 0..5 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn f64_in_unit_range() {
        let mut rng = SplitMix64::new(7);
        for _ in 0..1_000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v));
        }
    }

    #[test]
    fn shuffle_preserves_elements() {
        let mut rng = XorShift64::new(42);
        let mut v: Vec<i32> = (0..10).collect();
        let mut expected = v.clone();
        shuffle(&mut v, &mut |m| rng.next_usize(m));
        expected.sort_unstable();
        v.sort_unstable();
        assert_eq!(v, expected);
    }

    #[test]
    fn unit_vector_length() {
        let mut rng = SplitMix64::new(9);
        for _ in 0..100 {
            let v = random_unit_vector(&mut rng);
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-9);
        }
    }
}
