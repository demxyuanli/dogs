//! Cutting strategies for `AdvApprox_ApproxAFunction`.
//!
//! Source:
//! - `AdvApprox_Cutting.hxx:28-38` (abstract `Value`)
//! - `AdvApprox_DichoCutting.cxx:20-28` (bisection)
//! - `AdvApprox_PrefCutting.cxx:20-44` (preferred points)
//! - `AdvApprox_PrefAndRec.cxx:21-75` (recommended + preferred points)

use occt_core::precision::PCONFUSION;

/// `AdvApprox_Cutting` (`AdvApprox_Cutting.hxx:28-38`): "to choose the way of
/// cutting in approximation".
pub trait Cutting {
    /// `AdvApprox_Cutting::Value` (`AdvApprox_Cutting.hxx:35-37`): writes the
    /// cutting parameter into `cutting_value` and returns true when `[a, b]` is
    /// large enough to be split.
    fn value(&self, a: f64, b: f64, cutting_value: &mut f64) -> bool;
}

/// `AdvApprox_DichoCutting` (`AdvApprox_DichoCutting.hxx:29-37`): "if Cutting
/// is necessary in `[a,b]`, we cut at `(a+b)/2`".
#[derive(Clone, Copy, Debug, Default)]
pub struct DichoCutting;

impl Cutting for DichoCutting {
    /// `AdvApprox_DichoCutting::Value` (`AdvApprox_DichoCutting.cxx:22-28`).
    fn value(&self, a: f64, b: f64, cutting_value: &mut f64) -> bool {
        // Minimum length of an interval for F(U,V): EPS1=1.e-9 (cf. MEPS1)
        // (`cxx:24-25`).
        const LG_MIN: f64 = 10.0 * PCONFUSION;
        *cutting_value = (a + b) / 2.0;
        (b - a).abs() >= 2.0 * LG_MIN
    }
}

/// `AdvApprox_PrefCutting` (`AdvApprox_PrefCutting.hxx:30-41`): "contains a
/// list of preferential points (di)i; if Cutting is necessary in `[a,b]`, we
/// cut at the di nearest from `(a+b)/2`".
#[derive(Clone, Debug)]
pub struct PrefCutting {
    points: Vec<f64>,
}

impl PrefCutting {
    /// `AdvApprox_PrefCutting(CutPnts)` (`AdvApprox_PrefCutting.cxx:20-24`).
    pub fn new(points: Vec<f64>) -> Self {
        Self { points }
    }
}

impl Cutting for PrefCutting {
    /// `AdvApprox_PrefCutting::Value` (`AdvApprox_PrefCutting.cxx:26-44`).
    fn value(&self, a: f64, b: f64, cutting_value: &mut f64) -> bool {
        // Minimum length of a parametric interval: PConfusion() (`cxx:29-30`).
        const LG_MIN: f64 = 10.0 * PCONFUSION;
        let mil = (a + b) / 2.0;
        let mut cut = mil;
        let mut dist = ((a - b) / 2.0).abs();
        for &p in &self.points {
            if (dist - LG_MIN) > (mil - p).abs() {
                cut = p;
                dist = (mil - p).abs();
            }
        }
        *cutting_value = cut;
        (cut - a).abs() >= LG_MIN && (b - cut).abs() >= LG_MIN
    }
}

/// `AdvApprox_PrefAndRec` (`AdvApprox_PrefAndRec.hxx:31-52`): "contains a list
/// of preferential points (pi)i and a list of Recommended points used in
/// cutting management. if Cutting is necessary in `[a,b]`, we cut at the di
/// nearest from `(a+b)/2`".
#[derive(Clone, Debug)]
pub struct PrefAndRec {
    recommended: Vec<f64>,
    preferred: Vec<f64>,
    weight: f64,
}

impl PrefAndRec {
    /// `AdvApprox_PrefAndRec(RecCut, PrefCut, Weight=5)`
    /// (`AdvApprox_PrefAndRec.cxx:21-34`). A `Weight <= 1` raises
    /// `Standard_DomainError("PrefAndRec : Weight is too small")`
    /// (`cxx:30-33`), reported here as `Err(())`.
    pub fn new(recommended: Vec<f64>, preferred: Vec<f64>, weight: f64) -> Result<Self, ()> {
        if weight <= 1.0 {
            return Err(());
        }
        Ok(Self {
            recommended,
            preferred,
            weight,
        })
    }
}

impl Cutting for PrefAndRec {
    /// `AdvApprox_PrefAndRec::Value` (`AdvApprox_PrefAndRec.cxx:36-75`):
    /// first look for a **preferential** point closer to `mil = (a+b)/2` than the
    /// weighted window boundary `(r*a+b)/(r+1)` (`cxx:46-56`); if none was found,
    /// look for a **recommended** point that beats `|a-b|/2 - lgmin`
    /// (`cxx:58-70`); otherwise cut at `mil`.
    fn value(&self, a: f64, b: f64, cutting_value: &mut f64) -> bool {
        // Minimum length of a parametric interval: 10*PConfusion()
        // (`cxx:38-39`).
        const LG_MIN: f64 = 10.0 * PCONFUSION;
        let mil = (a + b) / 2.0;
        let mut cut = mil;
        let mut isfound = false;

        // Search for a preferred cutting point (`cxx:46-56`).
        let mut dist = ((a * self.weight + b) / (1.0 + self.weight) - mil).abs();
        for &p in &self.preferred {
            if dist > (mil - p).abs() {
                cut = p;
                dist = (mil - cut).abs();
                isfound = true;
            }
        }

        // Search for a recommended cutting point (`cxx:58-70`).
        if !isfound {
            dist = ((a - b) / 2.0).abs();
            for &p in &self.recommended {
                if (dist - LG_MIN) > (mil - p).abs() {
                    cut = p;
                    dist = (mil - cut).abs();
                }
            }
        }

        // Result (`cxx:72-74`).
        *cutting_value = cut;
        (cut - a).abs() >= LG_MIN && (b - cut).abs() >= LG_MIN
    }
}
