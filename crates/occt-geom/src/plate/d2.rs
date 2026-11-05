//! `Plate_D2` — the order-2 derivatives of a 3d-valued function of a 2d
//! variable.
//!
//! Source: `Plate_D2.hxx:24-37`, `Plate_D2.cxx:20-29`.

use occt_core::gp::GpXyz;

/// `Plate_D2` (`hxx:28-37`): `Duu`, `Duv`, `Dvv`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateD2 {
    duu: GpXyz,
    duv: GpXyz,
    dvv: GpXyz,
}

impl PlateD2 {
    /// `Plate_D2(Duu, Duv, Dvv)` (`cxx:20-25`).
    pub fn new(duu: GpXyz, duv: GpXyz, dvv: GpXyz) -> Self {
        Self { duu, duv, dvv }
    }

    /// `Duu`.
    pub fn duu(&self) -> &GpXyz {
        &self.duu
    }

    /// `Duv`.
    pub fn duv(&self) -> &GpXyz {
        &self.duv
    }

    /// `Dvv`.
    pub fn dvv(&self) -> &GpXyz {
        &self.dvv
    }
}
