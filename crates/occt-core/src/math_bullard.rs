//! `math_BullardGenerator` (`math_BullardGenerator.hxx:27-56`): the xorshift
//! style generator proposed by Ian C. Bullard. Default seed is 1.

/// Random number generator used by the OCCT ports that draw probe parameters.
#[derive(Debug, Clone)]
pub struct BullardGenerator {
    hi: u32,
    lo: u32,
}

impl Default for BullardGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl BullardGenerator {
    /// `math_BullardGenerator()` with the default seed 1 (`SetSeed(1)`).
    pub fn new() -> Self {
        Self {
            hi: 1,
            lo: 1 ^ 0x49616E42,
        }
    }

    /// `NextInt()`.
    pub fn next_int(&mut self) -> u32 {
        self.hi = self.hi.wrapping_shr(2).wrapping_add(self.hi.wrapping_shl(2));
        self.hi = self.hi.wrapping_add(self.lo);
        self.lo = self.lo.wrapping_add(self.hi);
        self.hi
    }

    /// `NextReal()`: a value in [0, 1].
    pub fn next_real(&mut self) -> f64 {
        self.next_int() as f64 / u32::MAX as f64
    }
}
