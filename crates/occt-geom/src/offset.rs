//! Offset 3D curve. Source: `Geom_OffsetCurve.hxx`
//!
//! Point and derivative formulas are a port of `Geom_OffsetCurveUtils.pxx`:
//! `CalculateD0` (`:47-63`), `CalculateD1` (`:74-116`), `CalculateD2` (`:129-201`).
//! The offset point moves along the **local normal** `Ndir = D1 ^ Direction`
//! (not along `Direction`): `P(u) = p(u) + Offset * Ndir / ||Ndir||`.
//! `EvalD0/EvalD1/EvalD2` wrappers mirror `Geom_OffsetCurve.cxx:262-334`.
//!
//! UNPORTED (documented deviations, no OCCT invention added):
//! - `CalculateD2` needs the basis **third** derivative (`EvalD3`). No curve in
//!   this crate overrides `Curve::d3`, so the trait default (zero) is used and
//!   the `D2Ndir` term of `D2` is zero-valued until `EvalD3` is ported per
//!   curve type (task T-63; OCCT `Geom_Curve::EvalD3`). D0/D1 are unaffected.
//! - `EvalD2`'s singular arm (`Geom_OffsetCurve.cxx:311-330`) computes
//!   `isDirectionChange` through `Geom_OffsetCurveUtils::AdjustDerivative`
//!   (`pxx:313-...`) when the basis `D1` magnitude is `<= gp::Resolution()`;
//!   the port passes `false` (same task T-63).
//! - OCCT throws `Geom_UndefinedValue`/`Geom_UndefinedDerivative` when a
//!   `CalculateD*` call returns false (`cxx:271-274`, `:291-294`, `:327-331`).
//!   The `Curve` trait has no failure channel, so the port returns the basis
//!   value, matching the existing convention in `bspline_surface.rs:182-185`
//!   and `offset_surface.rs:160-165`.
//! - Constructor checks (`Geom_OffsetCurve.cxx:187-217`: C0 rejection, G1
//!   upgrade of a C0 B-spline basis, direction-magnitude folding at `:181-183`)
//!   are not ported; this crate constructs the offset directly.
use std::sync::Arc;
use occt_core::gp::{GpDir, GpPnt, GpTrsf, GpVec, GpXyz};
use crate::curve::Curve;

/// `gp::Resolution()` (`gp.hxx:60`) is `RealSmall()` = `DBL_MIN`.
/// Note: `occt_core::precision::RESOLUTION` is `1e-12`, which is **not** OCCT's
/// `gp::Resolution()`; this file uses the faithful constant.
const GP_RESOLUTION: f64 = occt_core::precision::REAL_SMALL;

/// `Geom_OffsetCurve::Continuity` (`Geom_OffsetCurve.cxx:229-257`).
/// The offset of a C1 basis is only C0; C2 -> C1; C3 -> C2; G1/G2/CN kept.
/// Numeric codes follow `GeomAbs_Shape` (C0=0, G1=1, C1=2, G2=3, C2=4, C3=5, CN=6).
fn offset_continuity(basis_continuity: u8) -> u8 {
    use occt_core::kernel::geomabs::Shape;
    let shape = match basis_continuity {
        0 => Shape::C0,
        1 => Shape::G1,
        2 => Shape::C1,
        3 => Shape::G2,
        4 => Shape::C2,
        5 => Shape::C3,
        _ => Shape::CN,
    };
    let offset_shape = match shape {
        Shape::C0 => Shape::C0,
        Shape::G1 => Shape::G1,
        Shape::C1 => Shape::C0,
        Shape::G2 => Shape::G2,
        Shape::C2 => Shape::C1,
        Shape::C3 => Shape::C2,
        Shape::CN => Shape::CN,
    };
    offset_shape as u8
}

#[derive(Clone)]
pub struct GeomOffsetCurve { basis: Arc<dyn Curve>, offset: f64, direction: GpDir }

impl GeomOffsetCurve {
    pub fn new(curve: Arc<dyn Curve>, offset: f64, dir: GpDir) -> Self { Self { basis: curve, offset, direction: dir } }

    /// `Geom_OffsetCurveUtils::CalculateD0` (`Geom_OffsetCurveUtils.pxx:47-63`).
    fn calculate_d0(&self, value: &mut GpPnt, d1: &GpVec, tolerance: f64) -> bool {
        let mut ndir = d1.xyz().crossed(self.direction.xyz());
        let r = ndir.modulus();
        if r <= tolerance {
            return false;
        }
        ndir = ndir.multiplied(self.offset / r);
        *value = GpPnt::from_xyz(&value.coord.added(&ndir));
        true
    }

    /// `Geom_OffsetCurveUtils::CalculateD1` (`Geom_OffsetCurveUtils.pxx:74-116`).
    fn calculate_d1(&self, value: &mut GpPnt, d1: &mut GpVec, d2: &GpVec, tolerance: f64) -> bool {
        let mut ndir = d1.xyz().crossed(self.direction.xyz());
        let mut dndir = d2.xyz().crossed(self.direction.xyz());
        let r2 = ndir.square_modulus();
        let r = r2.sqrt();
        let r3 = r * r2;
        let dr = ndir.dot(&dndir);
        if r3 <= tolerance {
            if r2 <= tolerance {
                return false;
            }
            // We try another computation but the stability is not very good.
            dndir = dndir.multiplied(r);
            dndir = dndir.subtracted(&ndir.multiplied(dr / r));
            dndir = dndir.multiplied(self.offset / r2);
        } else {
            // Same computation as IICURV in EUCLID-IS because the stability is better.
            dndir = dndir.multiplied(self.offset / r);
            dndir = dndir.subtracted(&ndir.multiplied(self.offset * dr / r3));
        }

        ndir = ndir.multiplied(self.offset / r);
        // P(u)
        *value = GpPnt::from_xyz(&value.coord.added(&ndir));
        // P'(u)
        *d1 = d1.added(&GpVec::from_xyz(&dndir));
        true
    }

    /// `Geom_OffsetCurveUtils::CalculateD2` (`Geom_OffsetCurveUtils.pxx:129-201`).
    fn calculate_d2(
        &self,
        value: &mut GpPnt,
        d1: &mut GpVec,
        d2: &mut GpVec,
        d3: &GpVec,
        is_dir_change: bool,
        tolerance: f64,
    ) -> bool {
        let mut ndir = d1.xyz().crossed(self.direction.xyz());
        let mut dndir = d2.xyz().crossed(self.direction.xyz());
        let mut d2ndir = d3.xyz().crossed(self.direction.xyz());
        let r2 = ndir.square_modulus();
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r5 = r3 * r2;
        let dr = ndir.dot(&dndir);
        let d2r = ndir.dot(&d2ndir) + dndir.dot(&dndir);

        if r5 <= tolerance {
            if r2 * r2 <= tolerance {
                return false;
            }
            // We try another computation but the stability is not very good
            // dixit ISG.
            //  V2 = P" (U) :
            let r4 = r2 * r2;
            d2ndir = d2ndir.subtracted(&dndir.multiplied(2.0 * dr / r2));
            d2ndir = d2ndir.added(&ndir.multiplied(((3.0 * dr * dr) / r4) - (d2r / r2)));
            d2ndir = d2ndir.multiplied(self.offset / r);

            // V1 = P' (U) :
            dndir = dndir.multiplied(r);
            dndir = dndir.subtracted(&ndir.multiplied(dr / r));
            dndir = dndir.multiplied(self.offset / r2);
        } else {
            // Same computation as IICURV in EUCLID-IS because the stability is better.
            // V2 = P" (U) :
            d2ndir = d2ndir.multiplied(self.offset / r);
            d2ndir = d2ndir.subtracted(&dndir.multiplied(2.0 * self.offset * dr / r3));
            d2ndir = d2ndir.added(&ndir.multiplied(self.offset * (((3.0 * dr * dr) / r5) - (d2r / r3))));

            // V1 = P' (U) :
            dndir = dndir.multiplied(self.offset / r);
            dndir = dndir.subtracted(&ndir.multiplied(self.offset * dr / r3));
        }

        ndir = ndir.multiplied(self.offset / r);
        // P(u)
        *value = GpPnt::from_xyz(&value.coord.added(&ndir));
        // P'(u) :
        *d1 = d1.added(&GpVec::from_xyz(&dndir));
        // P"(u) :
        if is_dir_change {
            *d2 = d2.reversed();
        }
        *d2 = d2.added(&GpVec::from_xyz(&d2ndir));
        true
    }
}

impl Curve for GeomOffsetCurve {
    /// `Geom_OffsetCurve::EvalD0` (`Geom_OffsetCurve.cxx:262-276`): the basis is
    /// evaluated through `EvalD1`, so the point comes from the D1 evaluation.
    fn d0(&self, u: f64) -> GpPnt {
        let (p, d1) = self.basis.d1(u);
        let mut value = p;
        if !self.calculate_d0(&mut value, &d1, GP_RESOLUTION) {
            // `cxx:273` throws `Geom_UndefinedValue`.
            return p;
        }
        value
    }

    /// `Geom_OffsetCurve::EvalD1` (`Geom_OffsetCurve.cxx:280-296`).
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let (p, d1, d2) = self.basis.d2(u);
        let mut value = p;
        let mut a_d1 = d1;
        if !self.calculate_d1(&mut value, &mut a_d1, &d2, GP_RESOLUTION) {
            // `cxx:293` throws `Geom_UndefinedDerivative`.
            return (p, d1);
        }
        (value, a_d1)
    }

    /// `Geom_OffsetCurve::EvalD2` (`Geom_OffsetCurve.cxx:300-334`).
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, d1, d2, d3) = self.basis.d3(u);
        let mut value = p;
        let mut a_d1 = d1;
        let mut a_d2 = d2;
        // `cxx:311-330` derives `isDirectionChange` from `AdjustDerivative` when
        // the basis D1 is singular; not ported (see header), passed as false.
        if !self.calculate_d2(&mut value, &mut a_d1, &mut a_d2, &d3, false, GP_RESOLUTION) {
            // `cxx:328` throws `Geom_UndefinedDerivative`.
            return (p, d1, d2);
        }
        (value, a_d1, a_d2)
    }

    fn first_parameter(&self) -> f64 { self.basis.first_parameter() }
    fn last_parameter(&self) -> f64 { self.basis.last_parameter() }
    /// `Geom_OffsetCurve::Continuity` (`Geom_OffsetCurve.cxx:229-257`).
    fn continuity(&self) -> u8 { offset_continuity(self.basis.continuity()) }

    /// `Geom_OffsetCurve::Transform` (`Geom_OffsetCurve.cxx:454-460`):
    /// basis, direction and offset are all transformed; the offset is scaled by
    /// the *signed* `ScaleFactor()` (no `abs`, unlike the 2D class).
    fn transform(&mut self, t: &GpTrsf) {
        let mut basis = self.basis.clone_dyn();
        basis.transform(t);
        self.basis = Arc::from(basis);

        let mut dir_xyz = GpXyz::new(self.direction.x(), self.direction.y(), self.direction.z());
        t.transforms_xyz_dir(&mut dir_xyz);
        if let Ok(d) = GpDir::from_xyz(&dir_xyz) {
            self.direction = d;
        }
        self.offset *= t.scale_factor();
    }

    /// `Geom_OffsetCurve::Reverse` (`Geom_OffsetCurve.cxx:95-100`): the basis is
    /// reversed **and** the offset is negated.
    fn reverse(&mut self) {
        let mut basis = self.basis.clone_dyn();
        basis.reverse();
        self.basis = Arc::from(basis);
        self.offset = -self.offset;
    }

    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
