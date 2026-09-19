//! Offset 2D curve. Source: `Geom2d_OffsetCurve.hxx`
//!
//! Point and derivative formulas are a port of `Geom2d_OffsetCurveUtils.pxx`:
//! `CalculateD0` (`:43-53`), `CalculateD1` (`:61-101`), `CalculateD2` (`:111-177`).
//! The normal is `Ndir = (D1.Y, -D1.X)` (the tangent rotated by -90 degrees, see
//! `pxx:34`) and the offset point is `P(u) = p(u) + Offset * Ndir / ||Ndir||`.
//! `EvalD0/EvalD1/EvalD2` wrappers mirror `Geom2d_OffsetCurve.cxx:216-285`.
//!
//! UNPORTED (documented deviations, no OCCT invention added):
//! - `CalculateD2` needs the basis **third** derivative (`Geom2d_Curve::EvalD3`).
//!   No curve in this crate overrides `Curve2d::d3`, so the trait default (zero)
//!   is used and the `D2Ndir` term is zero-valued until `EvalD3` is ported per
//!   curve type (task T-63). D0/D1 are unaffected.
//! - `EvalD2`'s singular arm (`Geom2d_OffsetCurve.cxx:265-280`) derives
//!   `isDirectionChange` from `AdjustDerivative` when the basis `D1` is
//!   singular; the port passes `false` (task T-63).
//! - OCCT throws `Geom2d_UndefinedValue`/`Geom2d_UndefinedDerivative` when a
//!   `CalculateD*` call returns false (`cxx:224-227`, `:244-247`, `:280-283`).
//!   The `Curve2d` trait has no failure channel, so the port returns the basis
//!   value, matching the existing convention in `occt-geom/src/offset_surface.rs`.
//! - Constructor checks (`Geom2d_OffsetCurve.cxx:150`: C0 rejection and the G1
//!   upgrade path) are not ported.
use std::sync::Arc;
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use crate::curve::Curve2d;

/// `gp::Resolution()` (`gp.hxx:60`) is `RealSmall()` = `DBL_MIN`.
/// `occt_core::precision::RESOLUTION` is `1e-12`, which is **not** OCCT's value;
/// this file uses the faithful constant.
const GP_RESOLUTION: f64 = occt_core::precision::REAL_SMALL;

/// `Geom2d_OffsetCurve::Continuity` (`Geom2d_OffsetCurve.cxx:181-210`). Same table
/// as the 3D class: C1 -> C0, C2 -> C1, C3 -> C2, G1/G2/CN kept.
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
pub struct Geom2dOffsetCurve {
    basis: Arc<dyn Curve2d>,
    offset: f64,
}

impl Geom2dOffsetCurve {
    pub fn new(curve: Arc<dyn Curve2d>, offset: f64) -> Self { Self { basis: curve, offset } }
    pub fn basis_curve(&self) -> &Arc<dyn Curve2d> { &self.basis }
    pub fn offset_value(&self) -> f64 { self.offset }

    /// `Geom2d_OffsetCurveUtils::CalculateD0` (`Geom2d_OffsetCurveUtils.pxx:43-53`).
    fn calculate_d0(&self, value: &mut GpPnt2d, d1: &GpVec2d) -> bool {
        if d1.square_magnitude() <= GP_RESOLUTION {
            return false;
        }
        // `gp_Dir2d aNormal(theD1.Y(), -theD1.X())` — a `gp_Dir2d` is normalized.
        // The null-vector case is already excluded by the check above, exactly as
        // in OCCT (where only the `gp_Dir2d` ctor could throw).
        let normal = GpVec2d::new(d1.y(), -d1.x());
        let n = normal.multiplied_scalar(self.offset / normal.magnitude());
        *value = GpPnt2d::new(value.x() + n.x(), value.y() + n.y());
        true
    }

    /// `Geom2d_OffsetCurveUtils::CalculateD1` (`Geom2d_OffsetCurveUtils.pxx:61-101`).
    fn calculate_d1(&self, value: &mut GpPnt2d, d1: &mut GpVec2d, d2: &GpVec2d) -> bool {
        let ndir = GpVec2d::new(d1.y(), -d1.x());
        let mut dndir = GpVec2d::new(d2.y(), -d2.x());
        let r2 = ndir.square_magnitude();
        let r = r2.sqrt();
        let r3 = r * r2;
        let dr = ndir.dot(&dndir);
        if r3 <= GP_RESOLUTION {
            if r2 <= GP_RESOLUTION {
                return false;
            }
            // We try another computation but the stability is not very good.
            dndir = dndir.multiplied_scalar(r);
            dndir = dndir.subtracted(&ndir.multiplied_scalar(dr / r));
            dndir = dndir.multiplied_scalar(self.offset / r2);
        } else {
            // Same computation as IICURV in EUCLID-IS because the stability is better.
            dndir = dndir.multiplied_scalar(self.offset / r);
            dndir = dndir.subtracted(&ndir.multiplied_scalar(self.offset * dr / r3));
        }

        let n = ndir.multiplied_scalar(self.offset / r);
        // P(u)
        *value = GpPnt2d::new(value.x() + n.x(), value.y() + n.y());
        // P'(u)
        *d1 = d1.added(&dndir);
        true
    }

    /// `Geom2d_OffsetCurveUtils::CalculateD2` (`Geom2d_OffsetCurveUtils.pxx:111-177`).
    fn calculate_d2(
        &self,
        value: &mut GpPnt2d,
        d1: &mut GpVec2d,
        d2: &mut GpVec2d,
        d3: &GpVec2d,
        is_dir_change: bool,
    ) -> bool {
        let ndir = GpVec2d::new(d1.y(), -d1.x());
        let mut dndir = GpVec2d::new(d2.y(), -d2.x());
        let mut d2ndir = GpVec2d::new(d3.y(), -d3.x());
        let r2 = ndir.square_magnitude();
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r5 = r3 * r2;
        let dr = ndir.dot(&dndir);
        let d2r = ndir.dot(&d2ndir) + dndir.dot(&dndir);

        if r5 <= GP_RESOLUTION {
            if r2 * r2 <= GP_RESOLUTION {
                return false;
            }
            // We try another computation but the stability is not very good dixit ISG.
            let r4 = r2 * r2;
            //  V2 = P" (U) :
            d2ndir = d2ndir.subtracted(&dndir.multiplied_scalar(2.0 * dr / r2));
            d2ndir = d2ndir.added(&ndir.multiplied_scalar(((3.0 * dr * dr) / r4) - (d2r / r2)));
            d2ndir = d2ndir.multiplied_scalar(self.offset / r);

            // V1 = P' (U) :
            dndir = dndir.multiplied_scalar(r);
            dndir = dndir.subtracted(&ndir.multiplied_scalar(dr / r));
            dndir = dndir.multiplied_scalar(self.offset / r2);
        } else {
            // Same computation as IICURV in EUCLID-IS because the stability is better.
            // V2 = P" (U) :
            d2ndir = d2ndir.multiplied_scalar(self.offset / r);
            d2ndir = d2ndir.subtracted(&dndir.multiplied_scalar(2.0 * self.offset * dr / r3));
            d2ndir = d2ndir.added(&ndir.multiplied_scalar(self.offset * (((3.0 * dr * dr) / r5) - (d2r / r3))));

            // V1 = P' (U)
            dndir = dndir.multiplied_scalar(self.offset / r);
            dndir = dndir.subtracted(&ndir.multiplied_scalar(self.offset * dr / r3));
        }

        let n = ndir.multiplied_scalar(self.offset / r);
        // P(u)
        *value = GpPnt2d::new(value.x() + n.x(), value.y() + n.y());
        // P'(u) :
        *d1 = d1.added(&dndir);
        // P"(u) :
        if is_dir_change {
            *d2 = d2.reversed();
        }
        *d2 = d2.added(&d2ndir);
        true
    }
}

impl Curve2d for Geom2dOffsetCurve {
    /// `Geom2d_OffsetCurve::EvalD0` (`Geom2d_OffsetCurve.cxx:216-228`).
    fn d0(&self, u: f64) -> GpPnt2d {
        let (p, d1) = self.basis.d1(u);
        let mut value = p;
        if !self.calculate_d0(&mut value, &d1) {
            // `cxx:226` throws `Geom2d_UndefinedValue`.
            return p;
        }
        value
    }

    /// `Geom2d_OffsetCurve::EvalD1` (`Geom2d_OffsetCurve.cxx:232-248`).
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p, d1, d2) = self.basis.d2(u);
        let mut value = p;
        let mut a_d1 = d1;
        if !self.calculate_d1(&mut value, &mut a_d1, &d2) {
            // `cxx:246` throws `Geom2d_UndefinedDerivative`.
            return (p, d1);
        }
        (value, a_d1)
    }

    /// `Geom2d_OffsetCurve::EvalD2` (`Geom2d_OffsetCurve.cxx:252-285`).
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1, d2, d3) = self.basis.d3(u);
        let mut value = p;
        let mut a_d1 = d1;
        let mut a_d2 = d2;
        // `cxx:258-280` derives `isDirectionChange` from `AdjustDerivative` at a
        // singular basis D1; not ported (see header), passed as false.
        if !self.calculate_d2(&mut value, &mut a_d1, &mut a_d2, &d3, false) {
            // `cxx:282` throws `Geom2d_UndefinedDerivative`.
            return (p, d1, d2);
        }
        (value, a_d1, a_d2)
    }

    fn first_parameter(&self) -> f64 { self.basis.first_parameter() }
    fn last_parameter(&self) -> f64 { self.basis.last_parameter() }
    fn is_periodic(&self) -> bool { self.basis.is_periodic() }
    fn period(&self) -> f64 { self.basis.period() }
    /// `Geom2d_OffsetCurve::Continuity` (`Geom2d_OffsetCurve.cxx:181-210`).
    fn continuity(&self) -> u8 { offset_continuity(self.basis.continuity()) }

    /// `Geom2d_OffsetCurve::Transform` (`Geom2d_OffsetCurve.cxx:414-419`):
    /// the basis is transformed and the offset is scaled by `abs(ScaleFactor())`
    /// (the 2D class *does* use `abs`, unlike the 3D one).
    fn transform(&mut self, t: &GpTrsf2d) {
        let mut basis = self.basis.clone_dyn();
        basis.transform(t);
        self.basis = Arc::from(basis);
        self.offset *= t.scale_factor().abs();
    }

    /// `Geom2d_OffsetCurve::Reverse` (`Geom2d_OffsetCurve.cxx:90-95`): the basis is
    /// reversed **and** the offset is negated.
    fn reverse(&mut self) {
        let mut basis = self.basis.clone_dyn();
        basis.reverse();
        self.basis = Arc::from(basis);
        self.offset = -self.offset;
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }

    /// `Geom2d_OffsetCurve::BasisCurve()` (`Geom2d_OffsetCurve.cxx:174-177`).
    fn offset_basis(&self) -> Option<&dyn Curve2d> { Some(&*self.basis) }
}
