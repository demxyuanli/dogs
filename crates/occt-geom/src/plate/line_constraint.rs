//! `Plate_LineConstraint` — constrain a point to belong to a straight line,
//! expressed as two `LinearScalarConstraint` rows.
//!
//! Source: `Plate_LineConstraint.hxx:24-44`, `Plate_LineConstraint.cxx:19-47`,
//! `Plate_LineConstraint.lxx`.
//!
//! Translation note: `gp_XYZ::Normalize()` raises `Standard_ConstructionError`
//! on a zero vector (`gp_XYZ.hxx:598-604`); the port returns `Err` from
//! [`LineConstraint::new`] instead.

use occt_core::gp::{GpLin, GpXyz, GpXY};

use crate::plate::{LinearScalarConstraint, PinpointConstraint};

/// `Plate_LineConstraint` (`hxx:29-44`).
#[derive(Debug, Clone)]
pub struct LineConstraint {
    lsc: LinearScalarConstraint,
}

impl LineConstraint {
    /// `Plate_LineConstraint(point2d, lin, iu, iv)` (`cxx:19-47`).
    pub fn new(point2d: GpXY, lin: &GpLin, iu: i32, iv: i32) -> Result<Self, String> {
        let mut lsc = LinearScalarConstraint::with_size(2, 1);

        let point = lin.location().coord;
        lsc.set_ppc(1, PinpointConstraint::new(point2d, point, iu, iv));

        let dir = *lin.direction().xyz();
        // one builds two directions orthogonal to dir
        let d_x = GpXyz::new(1.0, 0.0, 0.0);
        let d_y = GpXyz::new(0.0, 1.0, 0.0);

        let mut d1 = d_x.crossed(&dir);
        let d2 = d_y.crossed(&dir);
        if d2.square_modulus() > d1.square_modulus() {
            d1 = d2;
        }
        if !d1.normalize() {
            return Err("Plate_LineConstraint: gp_XYZ::Normalize on a null direction".to_string());
        }
        let mut d2 = dir.crossed(&d1);
        if !d2.normalize() {
            return Err("Plate_LineConstraint: gp_XYZ::Normalize on a null direction".to_string());
        }
        lsc.set_coeff(1, 1, d1);
        lsc.set_coeff(2, 1, d2);

        Ok(Self { lsc })
    }

    /// `LSC()` (`hxx:41`).
    pub fn lsc(&self) -> &LinearScalarConstraint {
        &self.lsc
    }
}
