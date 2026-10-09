//! Port of `GeomLProp_CLProps2d` (`GeomLProp_CLProps.hxx:222-223`): the
//! `GeomLProp_CLPropsBase` template instantiated on 2D points, vectors and
//! directions, with the `LProp_CurveUtils::DirectAccess` policy.
//!
//! Sources: `LProp_CurveUtils.hxx` (`SetParameter`, `EnsureDeriv`,
//! `IsTangentDefined`, `Tangent`, `Curvature`, `Normal`, `ComputeTangent`,
//! `ComputeCurvature`, `ComputeNormal`).
//!
//! Not ported: `CentreOfCurvature` (no caller in the BRepClass chain).
//! The third derivative relies on `Curve2d::d3`, whose trait default returns
//! zero for curve types that do not override it.

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};

use crate::curve::Curve2d;

/// `LProp_Status` (`LProp_Status.hxx`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LPropStatus {
    Undecided,
    Defined,
    Undefined,
}

/// `RealLast()` (`Standard_Real.hxx`).
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx`).
const REAL_FIRST: f64 = -f64::MAX;

/// `GeomLProp_CLProps2d`: local properties of a 2D curve at a parameter.
pub struct GeomLPropCLProps2d<'a> {
    curve: &'a dyn Curve2d,
    u: f64,
    der_order: i32,
    cn: i32,
    lin_tol: f64,
    pnt: GpPnt2d,
    deriv: [GpVec2d; 3],
    curvature: f64,
    tangent_status: LPropStatus,
    sig_order: i32,
}

impl<'a> GeomLPropCLProps2d<'a> {
    /// `GeomLProp_CLPropsBase(C, U, N, Resolution)`.
    pub fn new(curve: &'a dyn Curve2d, u: f64, n: i32, resolution: f64) -> Self {
        assert!(
            (0..=3).contains(&n),
            "GeomLProp_CLProps::GeomLProp_CLProps(): invalid derivative order"
        );
        let mut s = Self {
            curve,
            u,
            der_order: n,
            cn: 4,
            lin_tol: resolution,
            pnt: GpPnt2d::new(0.0, 0.0),
            deriv: [GpVec2d::zero(), GpVec2d::zero(), GpVec2d::zero()],
            curvature: 0.0,
            tangent_status: LPropStatus::Undecided,
            sig_order: 0,
        };
        // `SetParameter(U)`: evaluate up to the requested order.
        s.eval_derivatives(n);
        s
    }

    /// `Value()`.
    pub fn value(&self) -> GpPnt2d {
        self.pnt
    }

    /// `D1()` (`EnsureDeriv` with order 1).
    pub fn d1(&mut self) -> GpVec2d {
        self.ensure_deriv(1)
    }

    /// `D2()` (`EnsureDeriv` with order 2).
    pub fn d2(&mut self) -> GpVec2d {
        self.ensure_deriv(2)
    }

    /// `D3()` (`EnsureDeriv` with order 3).
    pub fn d3(&mut self) -> GpVec2d {
        self.ensure_deriv(3)
    }

    /// `IsTangentDefined()` (`LProp_CurveUtils.hxx:333-381`): the first
    /// derivative of order 1..3 whose squared magnitude exceeds
    /// `lin_tol^2` gives the significant order.
    pub fn is_tangent_defined(&mut self) -> bool {
        match self.tangent_status {
            LPropStatus::Undefined => return false,
            LPropStatus::Defined => return true,
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

    /// `Tangent(D)` (`LProp_CurveUtils.hxx:391-401`, `ComputeTangent`
    /// `:179-217`). Panics where OCCT throws `LProp_NotDefined`.
    pub fn tangent(&mut self) -> GpDir2d {
        if !self.is_tangent_defined() {
            panic!("LProp_NotDefined: GeomLProp_CLProps2d::Tangent");
        }
        let sig = self.sig_order;
        if sig == 1 {
            return GpDir2d::from_vec2d(&self.deriv[0])
                .expect("GeomLProp_CLProps2d::Tangent: null first derivative");
        }
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
        let other_u = if u - u_inf < delta {
            u + delta
        } else {
            u - delta
        };
        let p1 = self.curve.d0(u.min(other_u));
        let p2 = self.curve.d0(u.max(other_u));
        let chord = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
        if av.dot(&chord) < 0.0 {
            av.multiply_scalar(-1.0);
        }
        GpDir2d::from_vec2d(&av).expect("GeomLProp_CLProps2d::Tangent: null derivative")
    }

    /// `Curvature()` (`LProp_CurveUtils.hxx:415-428`). Returns `RealLast`
    /// when the significant order is above 1. Panics where OCCT throws.
    pub fn curvature(&mut self) -> f64 {
        if !self.is_tangent_defined() {
            panic!("LProp_NotDefined: GeomLProp_CLProps2d::Curvature");
        }
        if self.sig_order > 1 {
            return REAL_LAST;
        }
        self.curvature = compute_curvature(
            &self.deriv[0],
            &self.deriv[1],
            self.lin_tol * self.lin_tol,
        );
        self.curvature
    }

    /// `Normal(N)` (`LProp_CurveUtils.hxx:436-443`, `ComputeNormal`
    /// `:249-253`). Panics where OCCT throws.
    pub fn normal(&mut self) -> GpDir2d {
        let curvature = self.curvature();
        if curvature == REAL_LAST || curvature.abs() <= self.lin_tol {
            panic!("LProp_NotDefined: GeomLProp_CLProps2d::Normal: null or infinite curvature");
        }
        let d1 = self.deriv[0];
        let d2 = self.deriv[1];
        let a = d1.dot(&d1);
        let b = d1.dot(&d2);
        let norm = GpVec2d::new(d2.x() * a - d1.x() * b, d2.y() * a - d1.y() * b);
        GpDir2d::from_vec2d(&norm).expect("GeomLProp_CLProps2d::Normal: null normal")
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

    /// `EnsureDeriv` (`LProp_CurveUtils.hxx:306-321`).
    fn ensure_deriv(&mut self, required: i32) -> GpVec2d {
        if self.der_order < required {
            self.der_order = required;
            self.eval_derivatives(self.der_order);
        }
        self.deriv[(required - 1) as usize]
    }
}

/// `ComputeCurvature` (`LProp_CurveUtils.hxx:229-241`): `|D1 x D2| / |D1|^3`,
/// or 0 when `D2` or the turn is below the squared tolerance.
fn compute_curvature(d1: &GpVec2d, d2: &GpVec2d, tol_sq: f64) -> f64 {
    let dd1 = d1.square_magnitude();
    let dd2 = d2.square_magnitude();
    if dd2 <= tol_sq {
        return 0.0;
    }
    let cross = d1.crossed(d2);
    let n = cross * cross;
    let t = n / dd1 / dd2;
    if t <= tol_sq {
        return 0.0;
    }
    n.sqrt() / dd1 / dd1.sqrt()
}
