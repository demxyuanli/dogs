//! `GeomFill_FunctionDraft` — the residual `G(w(t)) - S(u(t),v(t))` of the
//! draft-location "generatrix against stop surface" solve.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_FunctionDraft.cxx:36-212` and
//! `GeomFill_FunctionDraft.hxx:33-95` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Adaptor3d_Curve` / `Adaptor3d_Surface` are [`Curve`] / [`Surface`].
//! - `math_Vector` / `math_Matrix` are [`MathVector`] / [`MathMatrix`]; the
//!   coordinate loop over `1..=3` keeps OCCT's `i` index.
//! - `Value` / `Derivatives` / `Values` (`cxx:64-121`) are the
//!   `math_FunctionSetWithDerivatives` overrides, so they are the
//!   [`MathFunctionSetWithDerivatives`] methods here.
//! - `DerivT` (`cxx:126-146`), `Deriv2T` (`cxx:151-167`), `DerivTX`
//!   (`cxx:175-187`) and `Deriv2X` (`cxx:193-212`) keep the OCCT signatures
//!   (curve passed in explicitly, because the caller may solve on the trimmed
//!   curve rather than on `TheCurve`).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use crate::curve::Curve;
use crate::geom_fill::tensor::Tensor;
use crate::surface::Surface;

/// `gp_Vec::Coord(i)` with OCCT's 1-based index.
fn comp(v: &GpVec, i: usize) -> f64 {
    match i {
        1 => v.coord.x,
        2 => v.coord.y,
        _ => v.coord.z,
    }
}

/// `gp_Pnt::Coord(i)` with OCCT's 1-based index.
fn comp_p(p: &GpPnt, i: usize) -> f64 {
    match i {
        1 => p.coord.x,
        2 => p.coord.y,
        _ => p.coord.z,
    }
}

/// `GeomFill_FunctionDraft` (`GeomFill_FunctionDraft.hxx:33-95`).
pub struct FunctionDraft {
    /// `TheCurve` (`hxx`): the generatrix.
    curve: Arc<dyn Curve>,
    /// `TheSurface` (`hxx`): the stop surface.
    surface: Arc<dyn Surface>,
}

impl FunctionDraft {
    /// `GeomFill_FunctionDraft(S, C)` (`cxx:36-42`).
    pub fn new(surface: Arc<dyn Surface>, curve: Arc<dyn Curve>) -> Self {
        Self { curve, surface }
    }

    /// `NbVariables` (`cxx:46-51`): `(t, u, v)`.
    pub fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations` (`cxx:56-61`).
    pub fn nb_equations(&self) -> usize {
        3
    }

    fn point_at(&self, x: &MathVector) -> GpPnt {
        self.curve.d0(x.value(1))
    }

    fn surface_point(&self, x: &MathVector) -> GpPnt {
        self.surface.d0(x.value(2), x.value(3))
    }

    /// `Value(X, F)` (`cxx:64-77`).
    pub fn value(&self, x: &MathVector, f: &mut MathVector) -> bool {
        let p = self.point_at(x);
        let p1 = self.surface_point(x);
        f.set_value(1, comp_p(&p, 1) - comp_p(&p1, 1));
        f.set_value(2, comp_p(&p, 2) - comp_p(&p1, 2));
        f.set_value(3, comp_p(&p, 3) - comp_p(&p1, 3));
        true
    }

    /// `Derivatives(X, D)` (`cxx:82-97`).
    pub fn derivatives(&self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (_, dp) = self.curve.d1(x.value(1));
        let (_, dp1u, dp1v) = self.surface.d1(x.value(2), x.value(3));
        for i in 1..=3 {
            d.set_value(i, 1, comp(&dp, i));
            d.set_value(i, 2, -comp(&dp1u, i));
            d.set_value(i, 3, -comp(&dp1v, i));
        }
        true
    }

    /// `Values(X, F, D)` (`cxx:104-121`).
    pub fn values(&self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (p, dp) = self.curve.d1(x.value(1));
        let (p1, dp1u, dp1v) = self.surface.d1(x.value(2), x.value(3));
        for i in 1..=3 {
            f.set_value(i, comp_p(&p, i) - comp_p(&p1, i));
            d.set_value(i, 1, comp(&dp, i));
            d.set_value(i, 2, -comp(&dp1u, i));
            d.set_value(i, 3, -comp(&dp1v, i));
        }
        true
    }

    /// `DerivT(C, Param, W, dN, teta, F)` (`cxx:126-146`):
    /// `F = C'(Param) + W*sin(teta)*dN`.
    pub fn deriv_t(
        &self,
        c: &dyn Curve,
        param: f64,
        w: f64,
        dn: &GpVec,
        teta: f64,
        f: &mut MathVector,
    ) -> bool {
        let (_, dp) = c.d1(param);
        let s = w * teta.sin();
        for i in 1..=3 {
            f.set_value(i, comp(&dp, i) + s * comp(dn, i));
        }
        true
    }

    /// `Deriv2T(C, Param, W, d2N, teta, F)` (`cxx:151-167`):
    /// `F = C''(Param) + W*sin(teta)*d2N`.
    pub fn deriv2_t(
        &self,
        c: &dyn Curve,
        param: f64,
        w: f64,
        d2n: &GpVec,
        teta: f64,
        f: &mut MathVector,
    ) -> bool {
        let (_, _, d2p) = c.d2(param);
        let s = w * teta.sin();
        for i in 1..=3 {
            f.set_value(i, comp(&d2p, i) + s * comp(d2n, i));
        }
        true
    }

    /// `DerivTX(dN, teta, D)` (`cxx:175-187`): only the `W` column is filled.
    pub fn deriv_tx(&self, dn: &GpVec, teta: f64, d: &mut MathMatrix) -> bool {
        let s = teta.sin();
        for i in 1..=3 {
            d.set_value(i, 1, comp(dn, i) * s);
            d.set_value(i, 2, 0.0);
            d.set_value(i, 3, 0.0);
        }
        true
    }

    /// `Deriv2X(X, T)` (`cxx:193-212`): the surface second derivatives enter
    /// the `(u, u)`, `(u, v)` and `(v, v)` planes of the tensor.
    pub fn deriv2_x(&self, x: &MathVector, t: &mut Tensor) -> bool {
        let (_, _, _, d2pu, d2pv, d2puv) = self.surface.d2(x.value(2), x.value(3));
        t.init(0.0);
        for i in 1..=3 {
            t.set_value(i, 2, 2, -comp(&d2pu, i));
            t.set_value(i, 3, 2, -comp(&d2puv, i));
            t.set_value(i, 2, 3, -comp(&d2puv, i));
            t.set_value(i, 3, 3, -comp(&d2pv, i));
        }
        true
    }
}

impl MathFunctionSetWithDerivatives for FunctionDraft {
    fn nb_variables(&self) -> usize {
        FunctionDraft::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        FunctionDraft::nb_equations(self)
    }

    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        FunctionDraft::value(self, x, f)
    }

    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        FunctionDraft::derivatives(self, x, d)
    }

    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        FunctionDraft::values(self, x, f, d)
    }
}
