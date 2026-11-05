//! `Plate_PlaneConstraint` — constrain a point to belong to a plane, expressed
//! as one `LinearScalarConstraint` row.
//!
//! Source: `Plate_PlaneConstraint.hxx:24-43`,
//! `Plate_PlaneConstraint.cxx:19-33`, `Plate_PlaneConstraint.lxx`.
//!
//! Translation note: `gp_XYZ::Normalize()` raises `Standard_ConstructionError`
//! on a zero vector (`gp_XYZ.hxx:598-604`); the port returns `Err` instead.

use occt_core::gp::{GpPln, GpXY};

use crate::plate::{LinearScalarConstraint, PinpointConstraint};

/// `Plate_PlaneConstraint` (`hxx:28-43`).
#[derive(Debug, Clone)]
pub struct PlaneConstraint {
    lsc: LinearScalarConstraint,
}

impl PlaneConstraint {
    /// `Plate_PlaneConstraint(point2d, pln, iu, iv)` (`cxx:19-33`).
    pub fn new(point2d: GpXY, pln: &GpPln, iu: i32, iv: i32) -> Result<Self, String> {
        let mut lsc = LinearScalarConstraint::with_size(1, 1);

        let point = pln.location().coord;
        lsc.set_ppc(1, PinpointConstraint::new(point2d, point, iu, iv));
        let mut dir = *pln.axis().direction().xyz();
        if !dir.normalize() {
            return Err("Plate_PlaneConstraint: gp_XYZ::Normalize on a null direction".to_string());
        }
        lsc.set_coeff(1, 1, dir);
        Ok(Self { lsc })
    }

    /// `LSC()` (`hxx:40`).
    pub fn lsc(&self) -> &LinearScalarConstraint {
        &self.lsc
    }
}
