//! `Plate_PinpointConstraint` — a point and the value / derivative order
//! imposed on the plate function there.
//!
//! Source: `Plate_PinpointConstraint.hxx:24-45`,
//! `Plate_PinpointConstraint.cxx:19-38`, `Plate_PinpointConstraint.lxx`.

use occt_core::gp::{GpXY, GpXyz};

/// `Plate_PinpointConstraint` (`hxx:29-45`).
///
/// The members keep the OCCT declaration order (`value`, `pnt2d`, `idu`, `idv`)
/// but are all `i32`/value types, so no default initialisation is lost.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PinpointConstraint {
    value: GpXyz,
    pnt2d: GpXY,
    idu: i32,
    idv: i32,
}

impl Default for PinpointConstraint {
    /// `Plate_PinpointConstraint()` (`cxx:19-26`): null point, null value,
    /// both orders zero.
    fn default() -> Self {
        Self { value: GpXyz::zero(), pnt2d: GpXY::new(0.0, 0.0), idu: 0, idv: 0 }
    }
}

impl PinpointConstraint {
    /// `Plate_PinpointConstraint(point2d, ImposedValue, iu, iv)`
    /// (`cxx:28-38`), default orders `0`.
    pub fn new(point2d: GpXY, imposed_value: GpXyz, iu: i32, iv: i32) -> Self {
        Self { value: imposed_value, pnt2d: point2d, idu: iu, idv: iv }
    }

    /// `Pnt2d()` (`lxx:20`).
    pub fn pnt2d(&self) -> &GpXY {
        &self.pnt2d
    }

    /// `Idu()` (`lxx:26`).
    pub fn idu(&self) -> i32 {
        self.idu
    }

    /// `Idv()` (`lxx:32`).
    pub fn idv(&self) -> i32 {
        self.idv
    }

    /// `Value()` (`lxx:38`).
    pub fn value(&self) -> &GpXyz {
        &self.value
    }
}
