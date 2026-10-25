//! `GeomFill_PlanFunc` — the scalar function whose root is the parameter on the
//! guide curve of the point lying in the plane normal to the path.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_PlanFunc.hxx:32-76` and
//! `GeomFill_PlanFunc.cxx:22-91` (OCCT 8.0.0): constructor `cxx:22-29`,
//! `Value` `cxx:31-37`, `Derivative` `cxx:39-46`, `Values` `cxx:48-57`,
//! `D2` `cxx:63`, `DEDT` `cxx:65-70`, `D2E` `cxx:72-91`.
//!
//! Translation notes:
//! - `myCurve->D0` / `D1` are [`Curve::d0`] / [`Curve::d1`].
//! - `myPnt` / `myVec` are the plane point and normal (`hxx:67-68`); `V` and
//!   `G` (`hxx:69-71`) are scratch members of OCCT. They stay fields here so
//!   the `DEDT` / `D2E` bodies keep OCCT's evaluation order.
//! - `D2` (`cxx:63`) has an empty body in OCCT — the commented-out block at
//!   `cxx:59-62` is dead code — and is ported as a no-op.
//! - `math_FunctionWithDerivative` is [`FunctionWithDerivative`], the trait
//!   `occt_math::MathFunctionRoot` drives.
//! - `D2E` (`cxx:72-91`) evaluates `D1` into `G` but never reads the tangent,
//!   and builds `DVDT` from `G - DP` rather than from that tangent. Both are
//!   reproduced verbatim.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_math::FunctionWithDerivative;

use crate::curve::Curve;

/// `GeomFill_PlanFunc` (`GeomFill_PlanFunc.hxx:32-76`).
pub struct PlanFunc {
    /// `myPnt` (`hxx:67`): the point of the plane.
    my_pnt: GpPnt,
    /// `myVec` (`hxx:68`): the normal of the plane.
    my_vec: GpVec,
    /// `V` (`hxx:69`): scratch, `G - myPnt`.
    v: GpVec,
    /// `G` (`hxx:70`): scratch, the last point evaluated on the curve.
    g: GpPnt,
    /// `myCurve` (`hxx:71`).
    my_curve: Arc<dyn Curve>,
}

impl PlanFunc {
    /// `GeomFill_PlanFunc(theP, theV, theC)` (`cxx:22-29`): stores the point,
    /// the normal and the curve; `V` / `G` stay unset as in OCCT.
    pub fn new(the_p: GpPnt, the_v: GpVec, the_c: Arc<dyn Curve>) -> Self {
        Self {
            my_pnt: the_p,
            my_vec: the_v,
            v: GpVec::zero(),
            g: GpPnt::new(0.0, 0.0, 0.0),
            my_curve: the_c,
        }
    }

    /// `GeomFill_PlanFunc::D2` (`cxx:63`): empty in OCCT.
    pub fn d2(&mut self, _x: f64, _f: &mut f64, _d1: &mut f64, _d2: &mut f64) {}

    /// `GeomFill_PlanFunc::DEDT` (`cxx:65-70`): the partial derivative of `F`
    /// with respect to the plane itself, the plane moving by `DPnt` with
    /// normal velocity `DVec`.
    pub fn dedt(&mut self, x: f64, dpnt: &GpVec, dvec: &GpVec, dfdt: &mut f64) {
        self.g = self.my_curve.d0(x);
        self.v = GpVec::from_pnts(&self.my_pnt, &self.g);
        *dfdt = dvec.dot(&self.v) - self.my_vec.dot(dpnt);
    }

    /// `GeomFill_PlanFunc::D2E` (`cxx:72-91`). The two `D2` arguments `D2P` /
    /// `D2V` and the two outputs `D2FDT2` / `D2FDTDX` are unused by OCCT.
    #[allow(clippy::too_many_arguments)]
    pub fn d2e(
        &mut self,
        x: f64,
        dp: &GpVec,
        _d2p: &GpVec,
        dv: &GpVec,
        _d2v: &GpVec,
        dfdt: &mut f64,
        _d2fdt2: &mut f64,
        _d2fdt2dx: &mut f64,
    ) {
        let (g, _dg) = self.my_curve.d1(x);
        self.g = g;
        self.v = GpVec::from_pnts(&self.my_pnt, &self.g);
        let dvdt = GpVec::from_xyz(self.g.xyz()).subtracted(dp);
        *dfdt = dv.dot(&self.v) + self.my_vec.dot(&dvdt);
    }
}

impl FunctionWithDerivative for PlanFunc {
    /// `GeomFill_PlanFunc::Value` (`cxx:31-37`): `F = myVec . (G - myPnt)`.
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        self.g = self.my_curve.d0(x);
        self.v = GpVec::from_pnts(&self.my_pnt, &self.g);
        *f = self.my_vec.dot(&self.v);
        true
    }

    /// `GeomFill_PlanFunc::Derivative` (`cxx:39-46`): `D = myVec . dG/dX`.
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (g, dg) = self.my_curve.d1(x);
        self.g = g;
        *d = self.my_vec.dot(&dg);
        true
    }

    /// `GeomFill_PlanFunc::Values` (`cxx:48-57`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (g, dg) = self.my_curve.d1(x);
        self.g = g;
        self.v = GpVec::from_pnts(&self.my_pnt, &self.g);
        *f = self.my_vec.dot(&self.v);
        *d = self.my_vec.dot(&dg);
        true
    }
}
