//! `Plate_D1` — the order-1 derivatives of a 3d-valued function of a 2d
//! variable.
//!
//! Source: `Plate_D1.hxx:24-38`, `Plate_D1.cxx:20-28`, `Plate_D1.lxx`.

use occt_core::gp::GpXyz;

/// `Plate_D1` (`hxx:29-38`): `Du`, `Dv` with the accessors of `hxx:33-35`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateD1 {
    du: GpXyz,
    dv: GpXyz,
}

impl PlateD1 {
    /// `Plate_D1(Du, Dv)` (`cxx:20-24`).
    pub fn new(du: GpXyz, dv: GpXyz) -> Self {
        Self { du, dv }
    }

    /// `DU()` (`hxx:33`).
    pub fn du(&self) -> &GpXyz {
        &self.du
    }

    /// `DV()` (`hxx:35`).
    pub fn dv(&self) -> &GpXyz {
        &self.dv
    }
}
