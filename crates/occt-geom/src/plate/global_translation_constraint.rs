//! `Plate_GlobalTranslationConstraint` — force a set of UV points to translate
//! without deformation.
//!
//! Source: `Plate_GlobalTranslationConstraint.hxx:24-42`,
//! `Plate_GlobalTranslationConstraint.cxx:19-41`,
//! `Plate_GlobalTranslationConstraint.lxx`.

use occt_core::gp::{GpXyz, GpXY};

use crate::plate::{LinearXyzConstraint, PinpointConstraint};

/// `Plate_GlobalTranslationConstraint` (`hxx:28-42`).
#[derive(Debug, Clone)]
pub struct GlobalTranslationConstraint {
    lxyzc: LinearXyzConstraint,
}

impl GlobalTranslationConstraint {
    /// `Plate_GlobalTranslationConstraint(SOfXY)` (`cxx:19-41`).
    ///
    /// ```text
    /// myLXYZC(SOfXY.Length() - 1, SOfXY.Length())
    /// ```
    /// so `ColLen = Length - 1` rows and `RowLen = Length` columns. Every
    /// `SetPPC(i, ...)` writes `RowLen` entries; the coefficient loop writes
    /// the `Length - 1` rows of the `Length` columns.
    pub fn new(so_f_xy: &[GpXY]) -> Self {
        let len = so_f_xy.len();
        let mut lxyzc = LinearXyzConstraint::with_size(len.saturating_sub(1), len);

        for (i, xy) in so_f_xy.iter().enumerate() {
            lxyzc.set_ppc(i + 1, PinpointConstraint::new(*xy, GpXyz::zero(), 0, 0));
        }
        for i in 1..len {
            lxyzc.set_coeff(i, 1, -1.0);
            for j in 2..=len {
                if j == i + 1 {
                    lxyzc.set_coeff(i, j, 1.0);
                } else {
                    lxyzc.set_coeff(i, j, 0.0);
                }
            }
        }
        Self { lxyzc }
    }

    /// `LXYZC()` (`hxx:39`).
    pub fn lxyzc(&self) -> &LinearXyzConstraint {
        &self.lxyzc
    }
}
