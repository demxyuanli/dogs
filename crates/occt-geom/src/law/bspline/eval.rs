//! `Law_BSpline` evaluation — the `Local*` window evaluators and `Resolution`.
//!
//! Source: `Law_BSpline.cxx:1316-1450` (`LocalValue`/`LocalD0`/`LocalD1`/
//! `LocalD2`/`LocalD3`/`LocalDN`) and `:1705-1766` (`Resolution`), OCCT 8.0.0.
//!
//! Translation notes: see the module header of [`super`].

use occt_core::bspl::curve_1d::{curve_1d_d0, curve_1d_d1, curve_1d_d2, curve_1d_d3, curve_1d_dn};
use occt_core::bspl::curve_resolution::bspline_curve_resolution_1d;

use super::LawBSpline;

impl LawBSpline {
    /// `Law_BSpline::LocalValue` (`Law_BSpline.cxx:1405-1411`).
    pub fn local_value(&self, u: f64, from_k1: i32, to_k2: i32) -> Result<f64, String> {
        self.local_d0(u, from_k1, to_k2)
    }

    /// `Law_BSpline::LocalD0` (`Law_BSpline.cxx:1415-1430`).
    pub fn local_d0(&self, u: f64, from_k1: i32, to_k2: i32) -> Result<f64, String> {
        if from_k1 == to_k2 {
            return Err("Standard_DomainError: Law_BSpline::LocalValue".into());
        }
        let (u, index) = self.local_index(u, from_k1, to_k2);
        Ok(curve_1d_d0(
            u,
            index,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        ))
    }

    /// `Law_BSpline::LocalD1` (`Law_BSpline.cxx:1434-1453`).
    pub fn local_d1(
        &self,
        u: f64,
        from_k1: i32,
        to_k2: i32,
    ) -> Result<(f64, f64), String> {
        if from_k1 == to_k2 {
            return Err("Standard_DomainError: Law_BSpline::LocalD1".into());
        }
        let (u, index) = self.local_index(u, from_k1, to_k2);
        Ok(curve_1d_d1(
            u,
            index,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        ))
    }

    /// `Law_BSpline::LocalD2` (`Law_BSpline.cxx:1457-1478`).
    pub fn local_d2(
        &self,
        u: f64,
        from_k1: i32,
        to_k2: i32,
    ) -> Result<(f64, f64, f64), String> {
        if from_k1 == to_k2 {
            return Err("Standard_DomainError: Law_BSpline::LocalD2".into());
        }
        let (u, index) = self.local_index(u, from_k1, to_k2);
        Ok(curve_1d_d2(
            u,
            index,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        ))
    }

    /// `Law_BSpline::LocalD3` (`Law_BSpline.cxx:1483-1516`).
    pub fn local_d3(
        &self,
        u: f64,
        from_k1: i32,
        to_k2: i32,
    ) -> Result<(f64, f64, f64, f64), String> {
        if from_k1 == to_k2 {
            return Err("Standard_DomainError: Law_BSpline::LocalD3".into());
        }
        let (u, index) = self.local_index(u, from_k1, to_k2);
        Ok(curve_1d_d3(
            u,
            index,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        ))
    }

    /// `Law_BSpline::LocalDN` (`Law_BSpline.cxx:1519-1536`).
    pub fn local_dn(
        &self,
        u: f64,
        from_k1: i32,
        to_k2: i32,
        n: i32,
    ) -> Result<f64, String> {
        if from_k1 == to_k2 {
            return Err("Standard_DomainError: Law_BSpline::LocalDN".into());
        }
        let (u, index) = self.local_index(u, from_k1, to_k2);
        Ok(curve_1d_dn(
            u,
            n,
            index,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        ))
    }

    /// `Law_BSpline::Resolution` (`Law_BSpline.cxx:1705-1766`).
    pub fn resolution(&self, tolerance_3d: f64) -> f64 {
        bspline_curve_resolution_1d(
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            self.deg,
            tolerance_3d,
        )
    }
}
