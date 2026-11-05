//! `Plate_D3` — the order-3 derivatives of a 3d-valued function of a 2d
//! variable.
//!
//! Source: `Plate_D3.hxx:24-40`, `Plate_D3.cxx:20-31`.

use occt_core::gp::GpXyz;

/// `Plate_D3` (`hxx:28-40`): `Duuu`, `Duuv`, `Duvv`, `Dvvv`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateD3 {
    duuu: GpXyz,
    duuv: GpXyz,
    duvv: GpXyz,
    dvvv: GpXyz,
}

impl PlateD3 {
    /// `Plate_D3(Duuu, Duuv, Duvv, Dvvv)` (`cxx:20-27`).
    pub fn new(duuu: GpXyz, duuv: GpXyz, duvv: GpXyz, dvvv: GpXyz) -> Self {
        Self { duuu, duuv, duvv, dvvv }
    }

    /// `Duuu`.
    pub fn duuu(&self) -> &GpXyz {
        &self.duuu
    }

    /// `Duuv`.
    pub fn duuv(&self) -> &GpXyz {
        &self.duuv
    }

    /// `Duvv`.
    pub fn duvv(&self) -> &GpXyz {
        &self.duvv
    }

    /// `Dvvv`.
    pub fn dvvv(&self) -> &GpXyz {
        &self.dvvv
    }
}
