//! `GeomLProp_CLProps` — local properties of a 3D curve at one parameter.
//!
//! Source: `ModelingData/TKGeomBase/GeomLProp/GeomLProp_CLProps.hxx:47-219`
//! (the `GeomLProp_CLPropsBase` template instantiated as
//! `GeomLProp_CLPropsBase<gp_Pnt, gp_Vec, gp_Dir, occ::handle<Geom_Curve>>`)
//! and `ModelingData/TKGeomBase/LProp/LProp_CurveUtils.hxx`
//! (`SetParameter` `:284-298`, `EnsureDeriv` `:305-321`, `IsTangentDefined`
//! `:332-380`, `Tangent` `:390-402`, `ComputeTangent` `:171-218`, `Curvature`
//! `:412-427`, `ComputeCurvature` `:228-241`, `Normal` `:435-443`,
//! `ComputeNormal` `:250-253`) — OCCT 8.0.0.
//!
//! Translation notes:
//! - The OCCT `LProp_NotDefined` throws of `Tangent` / `Curvature` / `Normal`
//!   become `Err(String)`.
//! - `CentreOfCurvature` (`hxx:196-204`, `LProp_CurveUtils.hxx:451-463`) is not
//!   ported: nothing in this tree calls it.

use occt_core::gp::{GpDir, GpPnt, GpVec};

use crate::curve::Curve;

/// `LProp_Status` (`LProp_Status.hxx:22-28`). The discriminants matter: the
/// OCCT comparisons are `theStatus >= LProp_Defined` / `== LProp_Undefined`,
/// so `Undefined` is 1 and `Defined` is 2.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum LPropStatus {
    /// `LProp_Undecided`.
    Undecided,
    /// `LProp_Undefined`.
    Undefined,
    /// `LProp_Defined`.
    Defined,
    /// `LProp_Computed`.
    Computed,
}

/// `RealLast()` (`Standard_Real.hxx`).
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx`).
const REAL_FIRST: f64 = -f64::MAX;

/// `GeomLProp_CLProps`: point, derivatives up to order 3, tangent, curvature
/// and normal of a curve at a parameter (`GeomLProp_CLProps.hxx:47-219`).
pub struct ClProps<'a> {
    curve: &'a dyn Curve,
    u: f64,
    der_order: i32,
    cn: i32,
    lin_tol: f64,
    pnt: GpPnt,
    deriv: [GpVec; 3],
    curvature: f64,
    tangent_status: LPropStatus,
    sig_order: i32,
}

impl<'a> ClProps<'a> {
    /// `GeomLProp_CLPropsBase(C, U, N, Resolution)` (`hxx:79-88`): stores the
    /// parameter and evaluates the derivatives up to `n` (`SetParameter`,
    /// `LProp_CurveUtils.hxx:284-298`).
    pub fn new(curve: &'a dyn Curve, u: f64, n: i32, resolution: f64) -> Result<Self, String> {
        // `Standard_OutOfRange_Raise_if(N < 0 || N > 3, ...)` (`hxx:86`).
        if !(0..=3).contains(&n) {
            return Err("GeomLProp_CLProps: invalid derivative order".into());
        }
        let mut s = Self {
            curve,
            u,
            der_order: n,
            cn: 4,
            lin_tol: resolution,
            pnt: GpPnt::new(0.0, 0.0, 0.0),
            deriv: [GpVec::zero(), GpVec::zero(), GpVec::zero()],
            curvature: 0.0,
            tangent_status: LPropStatus::Undecided,
            sig_order: 0,
        };
        s.eval_derivatives(n);
        Ok(s)
    }

    /// `Value()` (`hxx:129`).
    pub fn value(&self) -> GpPnt {
        self.pnt
    }

    /// `D1()` (`hxx:134-137`, `EnsureDeriv` for order 1).
    pub fn d1(&mut self) -> GpVec {
        self.ensure_deriv(1)
    }

    /// `D2()` (`hxx:141-144`, `EnsureDeriv` for order 2).
    pub fn d2(&mut self) -> GpVec {
        self.ensure_deriv(2)
    }

    /// `D3()` (`hxx:148-151`, `EnsureDeriv` for order 3).
    pub fn d3(&mut self) -> GpVec {
        self.ensure_deriv(3)
    }

    /// `IsTangentDefined()` (`hxx:156-164`, `LProp_CurveUtils.hxx:332-380`):
    /// the first derivative of order 1..3 whose squared magnitude exceeds
    /// `lin_tol^2` gives the significant order.
    pub fn is_tangent_defined(&mut self) -> bool {
        match self.tangent_status {
            LPropStatus::Undefined => return false,
            LPropStatus::Defined | LPropStatus::Computed => return true,
            LPropStatus::Undecided => {}
        }
        let tol_sq = self.lin_tol * self.lin_tol;
        for order in 1..=4 {
            if self.cn < order {
                self.tangent_status = LPropStatus::Undefined;
                return false;
            }
            let av = match order {
                1 => self.d1(),
                2 => self.d2(),
                3 => self.d3(),
                _ => {
                    self.tangent_status = LPropStatus::Undefined;
                    return false;
                }
            };
            if av.square_magnitude() > tol_sq {
                self.sig_order = order;
                self.tangent_status = LPropStatus::Defined;
                return true;
            }
        }
        false
    }

    /// `Tangent(D)` (`hxx:167-176`, `LProp_CurveUtils.hxx:390-402` and
    /// `ComputeTangent` `:171-218`). `Err` where OCCT throws
    /// `LProp_NotDefined` or the `gp_Dir` constructor rejects a null vector.
    pub fn tangent(&mut self) -> Result<GpDir, String> {
        if !self.is_tangent_defined() {
            return Err("LProp_NotDefined: GeomLProp_CLProps::Tangent".into());
        }
        let sig = self.sig_order;
        if sig == 1 {
            return GpDir::from_vec(&self.deriv[0])
                .map_err(|_| "GeomLProp_CLProps::Tangent: null derivative".to_string());
        }

        // `ComputeTangent` (`LProp_CurveUtils.hxx:186-218`).
        let u = self.u;
        let u_sup = self.curve.last_parameter();
        let u_inf = self.curve.first_parameter();
        let du = if u_sup >= REAL_LAST || u_inf <= REAL_FIRST {
            0.0
        } else {
            u_sup - u_inf
        };
        let delta = (du * 1.0e-3).max(1.0e-7);
        let mut av = self.deriv[(sig - 1) as usize];
        let other_u = if u - u_inf < delta { u + delta } else { u - delta };
        let p1 = self.curve.d0(u.min(other_u));
        let p2 = self.curve.d0(u.max(other_u));
        let chord = GpVec::from_pnts(&p1, &p2);
        if av.dot(&chord) < 0.0 {
            av = av.multiply_scalar(-1.0);
        }
        GpDir::from_vec(&av).map_err(|_| "GeomLProp_CLProps::Tangent: null derivative".to_string())
    }

    /// `Curvature()` (`hxx:179-186`, `LProp_CurveUtils.hxx:412-427`): returns
    /// `RealLast` when the significant derivative order is above 1.
    pub fn curvature(&mut self) -> Result<f64, String> {
        if !self.is_tangent_defined() {
            return Err("LProp_NotDefined: GeomLProp_CLProps::Curvature".into());
        }
        if self.sig_order > 1 {
            return Ok(REAL_LAST);
        }
        self.curvature = compute_curvature(
            &self.deriv[0],
            &self.deriv[1],
            self.lin_tol * self.lin_tol,
        );
        Ok(self.curvature)
    }

    /// `Normal(N)` (`hxx:189-193`, `LProp_CurveUtils.hxx:435-443` and
    /// `ComputeNormal` `:250-253`).
    pub fn normal(&mut self) -> Result<GpDir, String> {
        let curvature = self.curvature()?;
        if curvature == REAL_LAST || curvature.abs() <= self.lin_tol {
            return Err(
                "LProp_NotDefined: GeomLProp_CLProps::Normal: null or infinite curvature".into(),
            );
        }
        let d1 = self.deriv[0];
        let d2 = self.deriv[1];
        let a = d1.dot(&d1);
        let b = d1.dot(&d2);
        let norm = d2.multiply_scalar(a).added(&d1.multiply_scalar(-b));
        GpDir::from_vec(&norm).map_err(|_| "GeomLProp_CLProps::Normal: null normal".to_string())
    }

    /// `CentreOfCurvature(P)` (`hxx:196-204`, `LProp_CurveUtils.hxx:451-463`
    /// and `ComputeCentreOfCurvature` `:264-274`).
    ///
    /// `|Curvature()| <= myLinTol` throws `LProp_NotDefined`; the centre is
    /// `P + normalize(D2*|D1|^2 - D1*(D1 . D2)) / curvature`.
    pub fn centre_of_curvature(&mut self) -> Result<GpPnt, String> {
        let curvature = self.curvature()?;
        if curvature.abs() <= self.lin_tol {
            return Err("LProp_NotDefined: GeomLProp_CLProps::CentreOfCurvature".into());
        }
        let d1 = self.deriv[0];
        let d2 = self.deriv[1];
        let a = d1.dot(&d1);
        let b = d1.dot(&d2);
        let mut norm = d2.multiply_scalar(a).subtracted(&d1.multiply_scalar(b));
        norm.normalize();
        Ok(self.pnt.translated_vec(&norm.divided(curvature)))
    }

    /// `SetParameter(U)` (`hxx:103-113`, `LProp_CurveUtils.hxx:287-298`):
    /// store the parameter, re-evaluate up to `myDerOrder`, reset the tangent
    /// status.
    pub fn set_parameter(&mut self, u: f64) {
        self.u = u;
        self.eval_derivatives(self.der_order);
        self.tangent_status = LPropStatus::Undecided;
    }

    /// `EvalDerivatives` (`LProp_CurveUtils.hxx:145-162`).
    fn eval_derivatives(&mut self, order: i32) {
        let u = self.u;
        match order {
            0 => self.pnt = self.curve.d0(u),
            1 => {
                let (p, v1) = self.curve.d1(u);
                self.pnt = p;
                self.deriv[0] = v1;
            }
            2 => {
                let (p, v1, v2) = self.curve.d2(u);
                self.pnt = p;
                self.deriv[0] = v1;
                self.deriv[1] = v2;
            }
            3 => {
                let (p, v1, v2, v3) = self.curve.d3(u);
                self.pnt = p;
                self.deriv[0] = v1;
                self.deriv[1] = v2;
                self.deriv[2] = v3;
            }
            _ => {}
        }
    }

    /// `EnsureDeriv` (`LProp_CurveUtils.hxx:305-321`).
    fn ensure_deriv(&mut self, required: i32) -> GpVec {
        if self.der_order < required {
            self.der_order = required;
            self.eval_derivatives(self.der_order);
        }
        self.deriv[(required - 1) as usize]
    }
}

/// `ComputeCurvature` (`LProp_CurveUtils.hxx:228-241`): `|D1 x D2| / |D1|^3`,
/// or 0 when `D2` or the turn are below the squared tolerance.
fn compute_curvature(d1: &GpVec, d2: &GpVec, tol_sq: f64) -> f64 {
    let dd1 = d1.square_magnitude();
    let dd2 = d2.square_magnitude();
    if dd2 <= tol_sq {
        return 0.0;
    }
    let cross = d1.crossed(d2);
    let n = cross.square_magnitude();
    let t = n / dd1 / dd2;
    if t <= tol_sq {
        return 0.0;
    }
    n.sqrt() / dd1 / dd1.sqrt()
}
