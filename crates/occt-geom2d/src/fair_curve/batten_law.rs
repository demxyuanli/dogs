//! `FairCurve_BattenLaw` (`FairCurve_BattenLaw.hxx:27-59`, `.cxx:16-23`,
//! `.lxx:17-36`): the height of a batten as a linear law of the parameter.
//!
//! Rust note: the setters are `const` in OCCT (`.lxx:17-29`) and are called on a
//! batten law shared by the energy and its Gauss-integration sub-functions, so
//! the three values live in `Cell`s.

use std::cell::Cell;

/// `FairCurve_BattenLaw` (`FairCurve_BattenLaw.hxx:27-59`).
#[derive(Debug, Clone, Default)]
pub struct BattenLaw {
    middle_heigth: Cell<f64>,
    geometric_slope: Cell<f64>,
    length_sliding: Cell<f64>,
}

impl BattenLaw {
    /// `FairCurve_BattenLaw(Heigth, Slope, Sliding)`
    /// (`FairCurve_BattenLaw.cxx:16-23`).
    pub fn new(heigth: f64, slope: f64, sliding: f64) -> Self {
        Self {
            middle_heigth: Cell::new(heigth),
            geometric_slope: Cell::new(slope),
            length_sliding: Cell::new(sliding),
        }
    }

    /// `SetSliding(Sliding)` (`FairCurve_BattenLaw.lxx:17-20`).
    pub fn set_sliding(&self, sliding: f64) {
        self.length_sliding.set(sliding);
    }

    /// `SetHeigth(Heigth)` (`FairCurve_BattenLaw.lxx:22-25`).
    pub fn set_heigth(&self, heigth: f64) {
        self.middle_heigth.set(heigth);
    }

    /// `SetSlope(Slope)` (`FairCurve_BattenLaw.lxx:27-30`).
    pub fn set_slope(&self, slope: f64) {
        self.geometric_slope.set(slope);
    }

    /// `Value(T, THeigth)` (`FairCurve_BattenLaw.lxx:32-36`): the height of the
    /// neutral fibre at `T`, `None` where OCCT's `Value` returns `false`
    /// (`THeigth <= 0`).
    pub fn value(&self, t: f64) -> Option<f64> {
        let t_heigth = self.middle_heigth.get()
            + (t - 0.5) * self.geometric_slope.get() * self.length_sliding.get();
        if t_heigth > 0.0 {
            Some(t_heigth)
        } else {
            None
        }
    }
}
