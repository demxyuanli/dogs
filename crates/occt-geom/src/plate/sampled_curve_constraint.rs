//! `Plate_SampledCurveConstraint` — `m` pinpoint constraints driven by `n`
//! unknowns through the "tent" blending function `B0`.
//!
//! Source: `Plate_SampledCurveConstraint.hxx:24-44`,
//! `Plate_SampledCurveConstraint.cxx:19-57`,
//! `Plate_SampledCurveConstraint.lxx`.
//!
//! Translation note: `Standard_DimensionMismatch` (`cxx:38-41`) becomes the
//! `Err` return of [`SampledCurveConstraint::new`].

use crate::plate::{LinearXyzConstraint, PinpointConstraint};

/// `B0(t)` (`cxx:22-33`): `max(0, 1 - |t|)`.
fn b0(t: f64) -> f64 {
    let mut s = t;
    if s < 0.0 {
        s = -s;
    }
    s = 1.0 - s;
    if s < 0.0 {
        s = 0.0;
    }
    s
}

/// `Plate_SampledCurveConstraint` (`hxx:28-44`).
#[derive(Debug, Clone)]
pub struct SampledCurveConstraint {
    lxyzc: LinearXyzConstraint,
}

impl SampledCurveConstraint {
    /// `Plate_SampledCurveConstraint(SOPPC, n)` (`cxx:35-57`).
    pub fn new(soppc: &[PinpointConstraint], n: usize) -> Result<Self, String> {
        let m = soppc.len();
        let mut lxyzc = LinearXyzConstraint::with_size(n, m);

        if n > m {
            return Err("Plate_SampledCurveConstraint: dimension mismatch".to_string());
        }
        for (index, ppc) in soppc.iter().enumerate() {
            lxyzc.set_ppc(index + 1, *ppc);
        }

        let ratio = (n + 1) as f64 / (m + 1) as f64;
        for i in 1..=n {
            for j in 1..=m {
                lxyzc.set_coeff(i, j, b0(ratio * j as f64 - i as f64));
            }
        }
        Ok(Self { lxyzc })
    }

    /// `LXYZC()` (`hxx:41`).
    pub fn lxyzc(&self) -> &LinearXyzConstraint {
        &self.lxyzc
    }
}
