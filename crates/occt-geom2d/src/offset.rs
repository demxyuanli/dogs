//! Offset 2D curve. Source: `Geom2d_OffsetCurve.hxx`
//!
//! Point and derivative formulas are a port of `Geom2d_OffsetCurveUtils.pxx`:
//! `CalculateD0` (`:43-53`), `CalculateD1` (`:61-101`), `CalculateD2` (`:111-177`),
//! `CalculateD3` (`:189-279`), `AdjustDerivative` (`:296-364`).
//! The normal is `Ndir = (D1.Y, -D1.X)` (the tangent rotated by -90 degrees, see
//! `pxx:34`) and the offset point is `P(u) = p(u) + Offset * Ndir / ||Ndir||`.
//! `EvalD0/EvalD1/EvalD2/EvalD3/EvalDN` wrappers mirror
//! `Geom2d_OffsetCurve.cxx:216-356`.
//!
//! UNPORTED (documented deviations, no OCCT invention added):
//! - `AdjustDerivative` consumes basis `EvalDN` orders up to 5 (`EvalD3` feeds it
//!   `EvalDN(U, 4)`, `pxx:474`). The elementary bases are faithful
//!   (`ElCLib::*DN` 2d, `clib2d.rs:329-432`) and `Geom2dTrimmedCurve` delegates to
//!   its basis (`cxx:273-283`); the remaining bases fall back to the trait default
//!   (zero above order 3) — `Geom2d_BSplineCurve::EvalDN` /
//!   `Geom2d_BezierCurve::EvalDN` have no counterpart here
//!   (`Geom2d_Curve::eval_dn`, `curve.rs`).
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

    /// `Geom2d_OffsetCurveUtils::CalculateD3` (`Geom2d_OffsetCurveUtils.pxx:189-279`).
    /// `d4` is the basis fourth derivative (`EvaluateD3` passes `EvalDN(U, 4)`,
    /// `pxx:474`).
    #[allow(clippy::too_many_arguments)]
    fn calculate_d3(
        &self,
        value: &mut GpPnt2d,
        d1: &mut GpVec2d,
        d2: &mut GpVec2d,
        d3: &mut GpVec2d,
        d4: &GpVec2d,
        is_dir_change: bool,
    ) -> bool {
        let ndir = GpVec2d::new(d1.y(), -d1.x());
        let mut dndir = GpVec2d::new(d2.y(), -d2.x());
        let mut d2ndir = GpVec2d::new(d3.y(), -d3.x());
        let mut d3ndir = GpVec2d::new(d4.y(), -d4.x());
        let r2 = ndir.square_magnitude();
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r4 = r2 * r2;
        let r5 = r3 * r2;
        let r6 = r3 * r3;
        let r7 = r5 * r2;
        let dr = ndir.dot(&dndir);
        let d2r = ndir.dot(&d2ndir) + dndir.dot(&dndir);
        let d3r = ndir.dot(&d3ndir) + 3.0 * dndir.dot(&d2ndir);

        if r7 <= GP_RESOLUTION {
            if r6 <= GP_RESOLUTION {
                return false;
            }
            // We try another computation but the stability is not very good dixit ISG.
            // V3 = P"' (U) :
            d3ndir = d3ndir.subtracted(&d2ndir.multiplied_scalar(3.0 * dr / r2));
            d3ndir = d3ndir.subtracted(&dndir.multiplied_scalar(3.0 * ((d2r / r2) + (dr * dr / r4))));
            d3ndir = d3ndir.added(&ndir.multiplied_scalar(
                6.0 * dr * dr / r4 + 6.0 * dr * d2r / r4 - 15.0 * dr * dr * dr / r6 - d3r,
            ));
            d3ndir = d3ndir.multiplied_scalar(self.offset / r);

            // V2 = P" (U) :
            d2ndir = d2ndir.subtracted(&dndir.multiplied_scalar(2.0 * dr / r2));
            d2ndir = d2ndir.subtracted(&ndir.multiplied_scalar((3.0 * dr * dr / r4) - (d2r / r2)));
            d2ndir = d2ndir.multiplied_scalar(self.offset / r);

            // V1 = P' (U) :
            dndir = dndir.multiplied_scalar(r);
            dndir = dndir.subtracted(&ndir.multiplied_scalar(dr / r));
            dndir = dndir.multiplied_scalar(self.offset / r2);
        } else {
            // Same computation as IICURV in EUCLID-IS because the stability is better.
            // V3 = P"' (U) :
            d3ndir = d3ndir.multiplied_scalar(self.offset / r);
            d3ndir = d3ndir.subtracted(&d2ndir.multiplied_scalar(3.0 * self.offset * dr / r3));
            d3ndir = d3ndir.subtracted(&dndir.multiplied_scalar(
                3.0 * self.offset * ((d2r / r3) + (dr * dr) / r5),
            ));
            d3ndir = d3ndir.added(&ndir.multiplied_scalar(
                self.offset * (6.0 * dr * dr / r5 + 6.0 * dr * d2r / r5 - 15.0 * dr * dr * dr / r7 - d3r),
            ));

            // V2 = P" (U) :
            d2ndir = d2ndir.multiplied_scalar(self.offset / r);
            d2ndir = d2ndir.subtracted(&dndir.multiplied_scalar(2.0 * self.offset * dr / r3));
            d2ndir = d2ndir.subtracted(&ndir.multiplied_scalar(
                self.offset * (((3.0 * dr * dr) / r5) - (d2r / r3)),
            ));

            // V1 = P' (U) :
            dndir = dndir.multiplied_scalar(self.offset / r);
            dndir = dndir.subtracted(&ndir.multiplied_scalar(self.offset * dr / r3));
        }

        let n = ndir.multiplied_scalar(self.offset / r);
        // P(u)
        *value = GpPnt2d::new(value.x() + n.x(), value.y() + n.y());
        // P'(u) :
        *d1 = d1.added(&dndir);
        // P"(u)
        *d2 = d2.added(&d2ndir);
        // P"'(u)
        if is_dir_change {
            *d3 = d3.reversed();
        }
        *d3 = d3.added(&d3ndir);
        true
    }
}

/// `Geom2d_OffsetCurveUtils::AdjustDerivative` (`Geom2d_OffsetCurveUtils.pxx:296-364`):
/// the 2D transcription of the 3D algorithm — first non-vanishing `EvalDN` order
/// in `{2, 3}` (`aMaxDerivOrder`), chord-sign selection with `aDelta =
/// max((u_sup - u_inf) * 1e-3, 1e-7)`, then `theD2..theD4` from the `EvalDN`
/// values one and two orders above, times the sign, and
/// `theIsDirectionChange = V.Dot(V1) < 0`.
///
/// OCCT returns `false` only when an `EvalDN` throws; the `Curve2d` trait has no
/// failure channel (see the header), so this always returns `true`.
fn adjust_derivative(
    curve: &dyn Curve2d,
    max_derivative: i32,
    u: f64,
    d1: &mut GpVec2d,
    d2: &mut GpVec2d,
    d3: &mut GpVec2d,
    d4: &mut GpVec2d,
    is_direction_change: &mut bool,
) -> bool {
    /// `gp::Resolution()` (`pxx:305`).
    const A_TOL: f64 = GP_RESOLUTION;
    /// `aMinStep` (`pxx:306`).
    const MIN_STEP: f64 = 1e-7;
    /// `aMaxDerivOrder` (`pxx:307`).
    const MAX_DERIV_ORDER: i32 = 3;
    /// `DivisionFactor` (`pxx:313`).
    const DIVISION_FACTOR: f64 = 1.0e-3;

    *is_direction_change = false;
    let u_inf = curve.first_parameter();
    let u_sup = curve.last_parameter();
    // `RealLast()` / `RealFirst()` (`Standard_Real.hxx`): an unbounded range
    // gives `du = 0` (so `aDelta = aMinStep`).
    let du = if u_sup >= f64::MAX || u_inf <= f64::MIN { 0.0 } else { u_sup - u_inf };
    let delta = (du * DIVISION_FACTOR).max(MIN_STEP);

    // Derivative is approximated by Taylor-series (`pxx:326-333`).
    let mut index = 1;
    let v = loop {
        index += 1;
        let current = curve.eval_dn(u, index);
        if current.square_magnitude() > A_TOL || index >= MAX_DERIV_ORDER {
            break current;
        }
    };

    let u_shift = if u - u_inf < delta { u + delta } else { u - delta };
    let p1 = curve.d0(u.min(u_shift));
    let p2 = curve.d0(u.max(u_shift));
    let v1 = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
    *is_direction_change = v.dot(&v1) < 0.0;
    let sign = if *is_direction_change { -1.0 } else { 1.0 };

    *d1 = v.multiplied_scalar(sign);
    let derivs = [d2, d3, d4];
    for i in 1..max_derivative {
        let dn = curve.eval_dn(u, index + i);
        *derivs[(i - 1) as usize] = dn.multiplied_scalar(sign);
    }
    true
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
        let mut a_d3 = d3;
        let mut is_direction_change = false;
        if a_d1.square_magnitude() <= GP_RESOLUTION {
            // `cxx:258-280`: the basis D1 is singular, so `D1..D3` and
            // `isDirectionChange` come from `AdjustDerivative(..., 3, ...)`.
            let mut a_dummy_d4 = GpVec2d::zero();
            // A `false` return throws `Geom2d_UndefinedDerivative` (`cxx:265-268`);
            // `adjust_derivative` has no failure channel and always succeeds.
            let _ = adjust_derivative(
                self.basis.as_ref(),
                3,
                u,
                &mut a_d1,
                &mut a_d2,
                &mut a_d3,
                &mut a_dummy_d4,
                &mut is_direction_change,
            );
        }
        if !self.calculate_d2(&mut value, &mut a_d1, &mut a_d2, &a_d3, is_direction_change) {
            // `cxx:282` throws `Geom2d_UndefinedDerivative`.
            return (p, d1, d2);
        }
        (value, a_d1, a_d2)
    }

    fn first_parameter(&self) -> f64 { self.basis.first_parameter() }
    fn last_parameter(&self) -> f64 { self.basis.last_parameter() }

    /// `Geom2d_OffsetCurve::EvalD3` (`Geom2d_OffsetCurve.cxx:289-327`) →
    /// `Geom2d_OffsetCurveUtils::EvaluateD3` (`pxx:460-484`) →
    /// [`Self::calculate_d3`]. The basis fourth derivative is `EvalDN(U, 4)`
    /// (`pxx:474`).
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        let (p, d1, d2, d3) = self.basis.d3(u);
        let mut a_d1 = d1;
        let mut a_d2 = d2;
        let mut a_d3 = d3;
        let mut a_d4 = self.basis.eval_dn(u, 4);
        let mut is_direction_change = false;
        if a_d1.square_magnitude() <= GP_RESOLUTION {
            // `cxx:302-311`: `AdjustDerivative(theBasisCurve, 4, ...)`.
            let _ = adjust_derivative(
                self.basis.as_ref(),
                4,
                u,
                &mut a_d1,
                &mut a_d2,
                &mut a_d3,
                &mut a_d4,
                &mut is_direction_change,
            );
        }
        let mut value = p;
        if !self.calculate_d3(
            &mut value,
            &mut a_d1,
            &mut a_d2,
            &mut a_d3,
            &a_d4,
            is_direction_change,
        ) {
            // `cxx:325` throws `Geom2d_UndefinedDerivative`.
            return (p, d1, d2, d3);
        }
        (value, a_d1, a_d2, a_d3)
    }

    /// `Geom2d_OffsetCurve::EvalDN` (`Geom2d_OffsetCurve.cxx:332-356`): orders
    /// 1..3 from `EvalD1`/`EvalD2`/`EvalD3`, **higher orders forwarded to the
    /// basis curve** (`cxx:355`). `N < 1` returns a zero vector instead of the
    /// OCCT throw (the `Curve2d::eval_dn` convention, `curve.rs`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        match n {
            i32::MIN..=0 => GpVec2d::zero(),
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            _ => self.basis.eval_dn(u, n),
        }
    }
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
    /// `Geom2d_OffsetCurve::TransformedParameter` (`Geom2d_OffsetCurve.cxx:423-426`).
    fn transformed_parameter(&self, u: f64, t: &GpTrsf2d) -> f64 {
        self.basis.transformed_parameter(u, t)
    }

    /// `Geom2d_OffsetCurve::BasisCurve()` (`Geom2d_OffsetCurve.cxx:174-177`).
    fn offset_basis(&self) -> Option<&dyn Curve2d> { Some(&*self.basis) }
}
