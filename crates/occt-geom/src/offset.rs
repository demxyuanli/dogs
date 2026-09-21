//! Offset 3D curve. Source: `Geom_OffsetCurve.hxx`
//!
//! Point and derivative formulas are a port of `Geom_OffsetCurveUtils.pxx`:
//! `CalculateD0` (`:47-63`), `CalculateD1` (`:74-116`), `CalculateD2` (`:129-201`),
//! `CalculateD3` (`:215-306`), `AdjustDerivative` (`:322-390`).
//! The offset point moves along the **local normal** `Ndir = D1 ^ Direction`
//! (not along `Direction`): `P(u) = p(u) + Offset * Ndir / ||Ndir||`.
//! `EvalD0/EvalD1/EvalD2/EvalD3/EvalDN` wrappers mirror
//! `Geom_OffsetCurve.cxx:261-410`.
//!
//! UNPORTED (documented deviations, no OCCT invention added):
//! - `AdjustDerivative` consumes basis `EvalDN` orders up to 5 (`EvalD3` feeds it
//!   `EvalDN(U, 4)`, `pxx:510`). `Curve::eval_dn` is faithful for B-spline
//!   (`BSplCLib::DN`, `bspline_curve.rs:159`) and for
//!   line/circle/ellipse/hyperbola/parabola (`ElCLib::*DN`, `clib.rs`); the
//!   remaining bases still fall back to the trait default (zero above order 3) —
//!   `Geom_BezierCurve::EvalDN` (`Geom_BezierCurve.cxx:601-617`) and
//!   `Geom_Curve`'s non-elementary subclasses are not ported to that depth.
//! - OCCT throws `Geom_UndefinedValue`/`Geom_UndefinedDerivative` when a
//!   `CalculateD*` call returns false (`cxx:271-274`, `:291-294`, `:324-336`).
//!   The `Curve` trait has no failure channel, so the port returns the basis
//!   value, matching the existing convention in `bspline_surface.rs:182-185`
//!   and `offset_surface.rs:160-165`. `AdjustDerivative`'s own `false` return
//!   is unreachable in the port for the same reason (`pxx:320`).
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

    /// `Geom_OffsetCurveUtils::CalculateD3` (`Geom_OffsetCurveUtils.pxx:215-306`).
    /// `d4` is the basis fourth derivative (`EvaluateD3` passes `EvalDN(U, 4)`,
    /// `pxx:510`). The `theIsDirChange` arm reverses `D3` before adding the
    /// normal term (`pxx:301-305`), mirroring `CalculateD2`'s `D2` handling.
    #[allow(clippy::too_many_arguments)]
    fn calculate_d3(
        &self,
        value: &mut GpPnt,
        d1: &mut GpVec,
        d2: &mut GpVec,
        d3: &mut GpVec,
        d4: &GpVec,
        is_dir_change: bool,
        tolerance: f64,
    ) -> bool {
        let mut ndir = d1.xyz().crossed(self.direction.xyz());
        let mut dndir = d2.xyz().crossed(self.direction.xyz());
        let mut d2ndir = d3.xyz().crossed(self.direction.xyz());
        let mut d3ndir = d4.xyz().crossed(self.direction.xyz());
        let r2 = ndir.square_modulus();
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r4 = r2 * r2;
        let r5 = r3 * r2;
        let r6 = r3 * r3;
        let r7 = r5 * r2;
        let dr = ndir.dot(&dndir);
        let d2r = ndir.dot(&d2ndir) + dndir.dot(&dndir);
        let d3r = ndir.dot(&d3ndir) + 3.0 * dndir.dot(&d2ndir);

        if r7 <= tolerance {
            if r6 <= tolerance {
                return false;
            }
            // We try another computation but the stability is not very good
            // dixit ISG.
            // V3 = P"' (U) :
            d3ndir = d3ndir.subtracted(&d2ndir.multiplied(3.0 * dr / r2));
            d3ndir = d3ndir.subtracted(&dndir.multiplied(3.0 * ((d2r / r2) + (dr * dr / r4))));
            d3ndir = d3ndir.added(&ndir.multiplied(
                6.0 * dr * dr / r4 + 6.0 * dr * d2r / r4 - 15.0 * dr * dr * dr / r6 - d3r,
            ));
            d3ndir = d3ndir.multiplied(self.offset / r);

            // V2 = P" (U) :
            d2ndir = d2ndir.subtracted(&dndir.multiplied(2.0 * dr / r2));
            d2ndir = d2ndir.subtracted(&ndir.multiplied((3.0 * dr * dr / r4) - (d2r / r2)));
            d2ndir = d2ndir.multiplied(self.offset / r);

            // V1 = P' (U) :
            dndir = dndir.multiplied(r);
            dndir = dndir.subtracted(&ndir.multiplied(dr / r));
            dndir = dndir.multiplied(self.offset / r2);
        } else {
            // Same computation as IICURV in EUCLID-IS because the stability is better.
            // V3 = P"' (U) :
            d3ndir = d3ndir.divided(r);
            d3ndir = d3ndir.subtracted(&d2ndir.multiplied(3.0 * dr / r3));
            d3ndir = d3ndir.subtracted(&dndir.multiplied(3.0 * ((d2r / r3) + (dr * dr) / r5)));
            d3ndir = d3ndir.added(&ndir.multiplied(
                6.0 * dr * dr / r5 + 6.0 * dr * d2r / r5 - 15.0 * dr * dr * dr / r7 - d3r,
            ));
            d3ndir = d3ndir.multiplied(self.offset);

            // V2 = P" (U) :
            d2ndir = d2ndir.divided(r);
            d2ndir = d2ndir.subtracted(&dndir.multiplied(2.0 * dr / r3));
            d2ndir = d2ndir.subtracted(&ndir.multiplied((3.0 * dr * dr / r5) - (d2r / r3)));
            d2ndir = d2ndir.multiplied(self.offset);

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
        *d2 = d2.added(&GpVec::from_xyz(&d2ndir));
        // P"'(u) :
        if is_dir_change {
            *d3 = d3.reversed();
        }
        *d3 = d3.added(&GpVec::from_xyz(&d3ndir));
        true
    }
}

/// `Geom_OffsetCurveUtils::AdjustDerivative` (`Geom_OffsetCurveUtils.pxx:322-390`):
/// at a singular parameter (basis `D1` magnitude `<= gp::Resolution()`) the
/// tangent is rebuilt from the first `EvalDN` order in `{2, 3}` (`aMaxDerivOrder`)
/// whose magnitude exceeds `gp::Resolution()`; its sign is chosen so that it
/// agrees with the chord `P(u - aDelta) P(u + aDelta)` (`aDelta =
/// max((u_sup - u_inf) * 1e-3, 1e-7)`), and `theD2..theD4` are the `EvalDN`
/// values one and two orders above it, times that same sign.
/// `theIsDirectionChange` is `V.Dot(V1) < 0`.
///
/// OCCT returns `false` only when an `EvalDN`/`EvalD0` throws; the `Curve` trait
/// has no failure channel (see the header), so this always returns `true`.
fn adjust_derivative(
    curve: &dyn Curve,
    max_derivative: i32,
    u: f64,
    d1: &mut GpVec,
    d2: &mut GpVec,
    d3: &mut GpVec,
    d4: &mut GpVec,
    is_direction_change: &mut bool,
) -> bool {
    /// `gp::Resolution()` (`pxx:331`).
    const A_TOL: f64 = GP_RESOLUTION;
    /// `aMinStep` (`pxx:332`).
    const MIN_STEP: f64 = 1e-7;
    /// `aMaxDerivOrder` (`pxx:333`).
    const MAX_DERIV_ORDER: i32 = 3;
    /// `DivisionFactor` (`pxx:339`).
    const DIVISION_FACTOR: f64 = 1.0e-3;

    *is_direction_change = false;
    let u_inf = curve.first_parameter();
    let u_sup = curve.last_parameter();
    // `RealLast()` / `RealFirst()` (`Standard_Real.hxx`): an unbounded range
    // gives `du = 0` (so `aDelta = aMinStep`).
    let du = if u_sup >= f64::MAX || u_inf <= f64::MIN { 0.0 } else { u_sup - u_inf };
    let delta = (du * DIVISION_FACTOR).max(MIN_STEP);

    // Derivative is approximated by Taylor-series (`pxx:352-359`).
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
    let v1 = GpVec::from_pnts(&p1, &p2);
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

impl Curve for GeomOffsetCurve {
    /// `Geom_OffsetCurve::BasisCurve()` / `Geom_OffsetCurve::Offset()`
    /// (`GeomAdaptor_Curve::GetType() == GeomAbs_OffsetCurve` companion).
    fn offset_curve(&self) -> Option<(std::sync::Arc<dyn Curve>, f64)> {
        Some((self.basis.clone(), self.offset))
    }

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
        let mut a_d3 = d3;
        let mut is_direction_change = false;
        if a_d1.square_magnitude() <= GP_RESOLUTION {
            // `cxx:311-330`: the basis D1 is singular, so `D1..D3` and
            // `isDirectionChange` come from `AdjustDerivative(..., 3, ...)`.
            let mut a_dummy_d4 = GpVec::zero();
            // A `false` return throws `Geom_UndefinedDerivative` (`cxx:317-320`);
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
        if !self.calculate_d2(
            &mut value,
            &mut a_d1,
            &mut a_d2,
            &a_d3,
            is_direction_change,
            GP_RESOLUTION,
        ) {
            // `cxx:328` throws `Geom_UndefinedDerivative`.
            return (p, d1, d2);
        }
        (value, a_d1, a_d2)
    }

    fn first_parameter(&self) -> f64 { self.basis.first_parameter() }
    fn last_parameter(&self) -> f64 { self.basis.last_parameter() }

    /// `Geom_OffsetCurve::EvalD3` (`Geom_OffsetCurve.cxx:342-380`) →
    /// `Geom_OffsetCurveUtils::EvaluateD3` (`pxx:495-529`) → [`Self::calculate_d3`].
    /// The basis fourth derivative is `EvalDN(U, 4)` (`pxx:510`) and feeds the
    /// singular arm's `AdjustDerivative(..., 4, ...)`.
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (p, d1, d2, d3) = self.basis.d3(u);
        let mut a_d1 = d1;
        let mut a_d2 = d2;
        let mut a_d3 = d3;
        let mut a_d4 = self.basis.eval_dn(u, 4);
        let mut is_direction_change = false;
        if a_d1.square_magnitude() <= GP_RESOLUTION {
            // `cxx:357-366`: `AdjustDerivative(theBasisCurve, 4, ...)`.
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
            GP_RESOLUTION,
        ) {
            // `cxx:379` throws `Geom_UndefinedDerivative`.
            return (p, d1, d2, d3);
        }
        (value, a_d1, a_d2, a_d3)
    }

    /// `Geom_OffsetCurve::EvalDN` (`Geom_OffsetCurve.cxx:386-410`): orders 1..3
    /// come from `EvalD1`/`EvalD2`/`EvalD3`, **every higher order is forwarded to
    /// the basis curve** (`cxx:409`). `N < 1` returns a zero vector instead of
    /// the OCCT throw (the `Curve::eval_dn` convention, `curve.rs:12-25`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        match n {
            i32::MIN..=0 => GpVec::zero(),
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            _ => self.basis.eval_dn(u, n),
        }
    }
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
