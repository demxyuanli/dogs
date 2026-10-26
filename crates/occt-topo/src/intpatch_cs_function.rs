//! `IntPatch_CSFunction` - the residual of a curve-on-surface / surface
//! intersection.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_CSFunction.hxx` and
//! `.cxx`.
//!
//! The unknown `X = (u, v, w)` is the parameter on `S1`, `S2` and the 2d curve
//! respectively. `F1..F3` are the three coordinates of
//! `S1(u, v) - S2(C(w))`, so `Root()` is the squared 3d distance of the two
//! evaluation points, and `Point()` is their midpoint (`cxx:56-66`).
//!
//! Ports `NbVariables` / `NbEquations` / `Value` / `Derivatives` / `Values` /
//! `Point` / `Root` / `AuxillarSurface` / `AuxillarCurve` (`cxx:44-142`).
//!
//! OCCT copies the two surface handles and the curve handle into `void*` slots
//! (`cxx:35-41`) and the analytic `IntPatch_CurvIntSurf` copies the whole
//! function by value (`IntImp_IntCS.gxx:27-39`). The port borrows them instead,
//! which is the same data without the handle indirection.

use occt_core::gp::GpPnt;
use occt_geom::Surface;
use occt_geom2d::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::intpatch_h_curve2d_tool as curve2d_tool;

/// `IntPatch_CSFunction` (`IntPatch_CSFunction.hxx:27-65`).
///
/// `Copy` mirrors OCCT's by-value copy of the function in
/// `IntImp_IntCS.gxx:35` (`myFunction(F)`); it only duplicates three borrowed
/// handles plus the last residual and midpoint.
#[derive(Clone, Copy)]
pub struct IntPatchCSFunction<'a> {
    /// `surface1`: the surface on which the intersection is searched
    /// (`IntPatch_CSFunction.cxx:17-19`).
    surface1: &'a dyn Surface,
    /// `surface2`: the surface `C` lives on.
    surface2: &'a dyn Surface,
    /// `curve`: `C`, a curve on `surface2`.
    curve: &'a dyn Curve2d,
    /// `p`: midpoint of the last evaluation (`cxx:123`).
    p: GpPnt,
    /// `f`: squared residual of the last evaluation (`cxx:121`).
    f: f64,
}

impl<'a> IntPatchCSFunction<'a> {
    /// `IntPatch_CSFunction(S1, C, S2)` (`cxx:33-42`). `f = 0.` there.
    pub fn new(surface1: &'a dyn Surface, curve: &'a dyn Curve2d, surface2: &'a dyn Surface) -> Self {
        Self {
            surface1,
            surface2,
            curve,
            p: GpPnt::zero(),
            f: 0.0,
        }
    }

    /// `Point()` (`cxx:125-128`).
    pub fn point(&self) -> &GpPnt {
        &self.p
    }

    /// `Root()` (`cxx:130-133`).
    pub fn root(&self) -> f64 {
        self.f
    }

    /// `AuxillarSurface()` (`cxx:135-138`).
    pub fn auxillar_surface(&self) -> &'a dyn Surface {
        self.surface1
    }

    /// `AuxillarCurve()` (`cxx:140-143`).
    pub fn auxillar_curve(&self) -> &'a dyn Curve2d {
        self.curve
    }

    /// The tail shared by `Value` (`cxx:56-66`) and `Values` (`cxx:114-122`):
    /// evaluate both surfaces, store the residual `f` and the midpoint `p`.
    fn set_state(&mut self, psurf: &GpPnt, pcurv: &GpPnt, f1: f64, f2: f64, f3: f64) {
        self.f = f1 * f1 + f2 * f2 + f3 * f3;
        self.p = GpPnt::from_xyz(&psurf.xyz().added(pcurv.xyz()).multiplied(0.5));
    }
}

impl MathFunctionSetWithDerivatives for IntPatchCSFunction<'_> {
    /// `NbVariables()` (`cxx:44-47`): 3.
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:49-52`): 3.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:54-67`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let psurf = self.surface1.d0(x.value(1), x.value(2));
        let p2d = curve2d_tool::value(self.curve, x.value(3));
        let pcurv = self.surface2.d0(p2d.x(), p2d.y());

        let f1 = psurf.x() - pcurv.x();
        let f2 = psurf.y() - pcurv.y();
        let f3 = psurf.z() - pcurv.z();
        f.set_value(1, f1);
        f.set_value(2, f2);
        f.set_value(3, f3);
        self.set_state(&psurf, &pcurv, f1, f2, f3);
        true
    }

    /// `Derivatives(X, D)` (`cxx:69-92`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (_, d1u, d1v) = self.surface1.d1(x.value(1), x.value(2));
        let (p2d, d2d) = curve2d_tool::d1(self.curve, x.value(3));
        let (_, du, dv) = self.surface2.d1(p2d.x(), p2d.y());
        let d1w = du.multiplied_scalar(d2d.x()).added(&dv.multiplied_scalar(d2d.y()));

        d.set_value(1, 1, d1u.x());
        d.set_value(1, 2, d1v.x());
        d.set_value(1, 3, -d1w.x());
        d.set_value(2, 1, d1u.y());
        d.set_value(2, 2, d1v.y());
        d.set_value(2, 3, -d1w.y());
        d.set_value(3, 1, d1u.z());
        d.set_value(3, 2, d1v.z());
        d.set_value(3, 3, -d1w.z());
        true
    }

    /// `Values(X, F, D)` (`cxx:94-124`): one `D1` pair of evaluations instead
    /// of separate `Derivatives` + `Value` passes.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (psurf, d1u, d1v) = self.surface1.d1(x.value(1), x.value(2));
        let (p2d, d2d) = curve2d_tool::d1(self.curve, x.value(3));
        let (pcurv, du, dv) = self.surface2.d1(p2d.x(), p2d.y());
        let d1w = du.multiplied_scalar(d2d.x()).added(&dv.multiplied_scalar(d2d.y()));

        d.set_value(1, 1, d1u.x());
        d.set_value(1, 2, d1v.x());
        d.set_value(1, 3, -d1w.x());
        d.set_value(2, 1, d1u.y());
        d.set_value(2, 2, d1v.y());
        d.set_value(2, 3, -d1w.y());
        d.set_value(3, 1, d1u.z());
        d.set_value(3, 2, d1v.z());
        d.set_value(3, 3, -d1w.z());
        let f1 = psurf.x() - pcurv.x();
        let f2 = psurf.y() - pcurv.y();
        let f3 = psurf.z() - pcurv.z();
        f.set_value(1, f1);
        f.set_value(2, f2);
        f.set_value(3, f3);
        self.set_state(&psurf, &pcurv, f1, f2, f3);
        true
    }
}
