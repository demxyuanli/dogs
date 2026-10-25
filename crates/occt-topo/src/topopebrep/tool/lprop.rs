//! `LProp_CurveUtils` — the curve-local-property helpers `BRepLProp_CLProps`
//! is built from.
//!
//! Sources:
//! - `TKGeomBase/LProp/LProp_CurveUtils.hxx:175-220` (`ComputeTangent`)
//! - `:226-242` (`ComputeCurvature`)
//! - `:246-255` (`ComputeNormal`)
//! - `:296-328` (`SetParameter` / `EnsureDeriv` / `EvalDerivatives`)
//! - `:332-378` (`IsTangentDefined`)
//! - `:390-405` (`Tangent`), `:410-427` (`Curvature`), `:433-442` (`Normal`)
//!
//! `BRepLProp_CLProps` is `GeomLProp_CLPropsBase<gp_Pnt, gp_Vec, gp_Dir,
//! BRepAdaptor_Curve>` (`BRepLProp_CLProps.hxx:21`), i.e. the same template
//! instantiated on the edge adaptor, so only one implementation is needed.
//!
//! `LProp_NotDefined` is raised by OCCT when the tangent is undefined
//! (`LProp_CurveUtils.hxx:399-400`, `:422`) and by `Normal` when the curvature
//! is null or infinite (`:439-440`); those are `Err` here.

use occt_core::gp::{GpDir, GpPnt, GpVec};

use super::adaptor::CurveAdaptor;

/// `LProp_Status` (`LProp_Status.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LPropStatus {
    /// `LProp_Undefined`
    Undefined = -1,
    /// `LProp_Undecided`
    Undecided = 0,
    /// `LProp_Defined`
    Defined = 1,
    /// `LProp_Computed`
    Computed = 2,
}

/// `GeomLProp_CLPropsBase` state: `myPnt` + `myDerivArr[3]` + `myTangentStatus`
/// + `mySignificantFirstDerivativeOrder`.
#[derive(Debug, Clone)]
pub struct CLProps {
    /// `myPnt`
    pub pnt: GpPnt,
    /// `myDerivArr[0..3]` (D1, D2, D3); a zero vector means "not evaluated".
    pub deriv: [GpVec; 3],
    /// `myLinTol` = the `Resolution` ctor argument.
    pub lin_tol: f64,
    /// `myCN` = `Adaptor3d_Curve::Continuity()`.
    pub cn: i32,
    /// `myTangentStatus`
    pub tangent_status: LPropStatus,
    /// `mySignificantFirstDerivativeOrder`
    pub sig_order: i32,
    /// `myCurvature`
    pub curvature: f64,
}

impl CLProps {
    /// `GeomLProp_CLPropsBase(C, N, Resolution)` followed by
    /// `SetParameter(U)` (`LProp_CurveUtils.hxx:286-299`).
    ///
    /// `N` is the highest derivative order OCCT is asked for
    /// (`Standard_OutOfRange_Raise_if(N < 0 || N > 3)`); an out-of-range `N` is
    /// `Err`. OCCT defers the evaluation to the first `D1/D2/D3` call
    /// (`EnsureDeriv`, `hxx:308-323`); the port evaluates `N` orders here,
    /// which is the same set of values.
    pub fn new(c: &CurveAdaptor, u: f64, n: i32, resolution: f64) -> Result<Self, String> {
        if !(0..=3).contains(&n) {
            return Err("GeomLProp_CLProps: N < 0 || N > 3".into());
        }
        let mut p = Self {
            pnt: GpPnt::new(0.0, 0.0, 0.0),
            deriv: [
                GpVec::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
            ],
            lin_tol: resolution,
            cn: c.continuity() as i32,
            tangent_status: LPropStatus::Undecided,
            sig_order: 0,
            curvature: 0.0,
        };
        p.set_parameter(c, u, n)?;
        Ok(p)
    }

    /// `SetParameter(U)` (`LProp_CurveUtils.hxx:286-299`): evaluates the
    /// derivatives up to `der_order` and resets the tangent status to
    /// `LProp_Undecided`.
    pub fn set_parameter(&mut self, c: &CurveAdaptor, u: f64, der_order: i32) -> Result<(), String> {
        self.eval_derivatives(c, u, der_order)?;
        self.tangent_status = LPropStatus::Undecided;
        Ok(())
    }

    /// `EvalDerivatives` (`LProp_CurveUtils.hxx:277-291`): `D0` for order 0,
    /// otherwise `D1`/`D2`/`D3`.
    fn eval_derivatives(&mut self, c: &CurveAdaptor, u: f64, der_order: i32) -> Result<(), String> {
        match der_order {
            0 => {
                self.pnt = c.d0(u).ok_or_else(|| "CLProps: no curve".to_string())?;
            }
            1 => {
                let (p, d1) = c.d1(u).ok_or_else(|| "CLProps: no curve".to_string())?;
                self.pnt = p;
                self.deriv[0] = d1;
            }
            2 => {
                let (p, d1, d2) = c.d2(u).ok_or_else(|| "CLProps: no curve".to_string())?;
                self.pnt = p;
                self.deriv[0] = d1;
                self.deriv[1] = d2;
            }
            _ => {
                let (p, d1, d2, d3) = c.d3(u).ok_or_else(|| "CLProps: no curve".to_string())?;
                self.pnt = p;
                self.deriv[0] = d1;
                self.deriv[1] = d2;
                self.deriv[2] = d3;
            }
        }
        Ok(())
    }

    /// `EnsureDeriv` (`LProp_CurveUtils.hxx:308-323`): upgrade the derivative
    /// set when a higher order is required.
    fn ensure_deriv(&mut self, c: &CurveAdaptor, u: f64, required: i32, der_order: &mut i32) -> GpVec {
        if *der_order < required {
            *der_order = required;
            let _ = self.eval_derivatives(c, u, required);
        }
        self.deriv[(required - 1) as usize].clone()
    }

    /// `IsTangentDefined()` (`LProp_CurveUtils.hxx:332-378`): the first
    /// non-null derivative among D1..D3, bounded by the continuity order.
    pub fn is_tangent_defined(&mut self, c: &CurveAdaptor, u: f64, der_order: &mut i32) -> bool {
        if self.tangent_status == LPropStatus::Undefined {
            return false;
        }
        if self.tangent_status >= LPropStatus::Defined {
            return true;
        }
        let tol_sq = self.lin_tol * self.lin_tol;
        let mut an_order = 0;
        while an_order < 4 {
            an_order += 1;
            if self.cn >= an_order {
                if an_order > 3 {
                    // `default: theTanStatus = LProp_Undefined; return false;`
                    self.tangent_status = LPropStatus::Undefined;
                    return false;
                }
                let v = self.ensure_deriv(c, u, an_order, der_order);
                if v.square_magnitude() > tol_sq {
                    self.sig_order = an_order;
                    self.tangent_status = LPropStatus::Defined;
                    return true;
                }
            } else {
                self.tangent_status = LPropStatus::Undefined;
                return false;
            }
        }
        false
    }

    /// `Tangent(D)` (`LProp_CurveUtils.hxx:390-405`): `Err` is
    /// `LProp_NotDefined`.
    pub fn tangent(&mut self, c: &CurveAdaptor, u: f64, der_order: &mut i32) -> Result<GpDir, String> {
        if !self.is_tangent_defined(c, u, der_order) {
            return Err("LProp_NotDefined".into());
        }
        compute_tangent(c, u, &self.deriv, self.sig_order)
    }

    /// `Curvature()` (`LProp_CurveUtils.hxx:410-427`): `RealLast()` when the
    /// first significant derivative is above order 1.
    pub fn curvature_of(&mut self, c: &CurveAdaptor, u: f64, der_order: &mut i32) -> Result<f64, String> {
        if !self.is_tangent_defined(c, u, der_order) {
            return Err("LProp_NotDefined".into());
        }
        if self.sig_order > 1 {
            return Ok(f64::MAX);
        }
        self.curvature = compute_curvature(&self.deriv[0], &self.deriv[1], self.lin_tol * self.lin_tol);
        Ok(self.curvature)
    }

    /// `Normal(N)` (`LProp_CurveUtils.hxx:433-442`): `Err` is
    /// `LProp_NotDefined("CLProps::Normal(): Curvature is null or infinity")`.
    pub fn normal(&mut self, c: &CurveAdaptor, u: f64, der_order: &mut i32) -> Result<GpDir, String> {
        let curvature = self.curvature_of(c, u, der_order)?;
        if curvature == f64::MAX || curvature.abs() <= self.lin_tol {
            return Err("CLProps::Normal(): Curvature is null or infinity".into());
        }
        compute_normal(&self.deriv[0], &self.deriv[1])
    }
}

/// `ComputeTangent` (`LProp_CurveUtils.hxx:175-220`).
///
/// `sig_order == 1` is `Dir(D1)`; otherwise the significant derivative is
/// oriented against a chord of length `max(span * 1e-3, 1e-7)` taken from `u`.
fn compute_tangent(c: &CurveAdaptor, u: f64, deriv: &[GpVec; 3], sig_order: i32) -> Result<GpDir, String> {
    if sig_order == 1 {
        return GpDir::from_vec(&deriv[0]).map_err(|e| e.to_string());
    }
    const THE_DIVISION_FACTOR: f64 = 1.0e-3;
    const THE_MIN_STEP: f64 = 1.0e-7;

    let u_supremum = c.last_parameter();
    let u_infimum = c.first_parameter();
    let a_du = if (u_supremum >= f64::MAX) || (u_infimum <= -f64::MAX) {
        0.0
    } else {
        u_supremum - u_infimum
    };
    let a_delta = (a_du * THE_DIVISION_FACTOR).max(THE_MIN_STEP);

    let mut a_v = deriv[(sig_order - 1) as usize].clone();
    let other_u = if u - u_infimum < a_delta {
        u + a_delta
    } else {
        u - a_delta
    };
    let p1 = c.d0(u.min(other_u)).ok_or_else(|| "ComputeTangent: no curve".to_string())?;
    let p2 = c.d0(u.max(other_u)).ok_or_else(|| "ComputeTangent: no curve".to_string())?;
    let a_chord = GpVec::from_pnts(&p1, &p2);
    if a_v.dot(&a_chord) < 0.0 {
        a_v = a_v.multiplied_scalar(-1.0);
    }
    GpDir::from_vec(&a_v).map_err(|e| e.to_string())
}

/// `ComputeCurvature` (`LProp_CurveUtils.hxx:226-242`):
/// `|D1 x D2| / |D1|^3`, with 0 when `D2` is null or the derivatives are
/// collinear.
fn compute_curvature(d1: &GpVec, d2: &GpVec, tol_sq: f64) -> f64 {
    let dd1 = d1.square_magnitude();
    let dd2 = d2.square_magnitude();
    if dd2 <= tol_sq {
        return 0.0;
    }
    let n = d1.crossed(d2).square_magnitude();
    let t = n / dd1 / dd2;
    if t <= tol_sq {
        return 0.0;
    }
    n.sqrt() / dd1 / dd1.sqrt()
}

/// `ComputeNormal` (`LProp_CurveUtils.hxx:246-255`):
/// `D2 * (D1 . D1) - D1 * (D1 . D2)`.
fn compute_normal(d1: &GpVec, d2: &GpVec) -> Result<GpDir, String> {
    let d1d1 = d1.dot(d1);
    let d1d2 = d1.dot(d2);
    let a_norm = d2
        .multiplied_scalar(d1d1)
        .subtracted(&d1.multiplied_scalar(d1d2));
    GpDir::from_vec(&a_norm).map_err(|e| e.to_string())
}

