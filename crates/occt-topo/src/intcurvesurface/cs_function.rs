//! `IntCurveSurface_TheCSFunctionOfHInter`.
//!
//! Source: `IntImp_ZerCSParFunc.gxx` as instantiated by
//! `IntCurveSurface_TheCSFunctionOfHInter_0.cxx` (`ThePSurfaceTool =
//! Adaptor3d_HSurfaceTool`, `TheCurveTool = IntCurveSurface_TheHCurveTool`).
//! The function set is the 3-equation / 3-variable system `S(u,v) - C(w)`
//! used by `IntCurveSurface_TheExactHInter` (instantiation of
//! `IntImp_IntCS.gxx`) to refine a start point coming out of the polygon /
//! polyhedron interference.

use occt_core::gp::GpPnt;
use occt_geom::{Curve, Surface};

use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

/// `IntImp_ZerCSParFunc` / `IntCurveSurface_TheCSFunctionOfHInter`
/// (`...hxx:30-61`).
pub struct TheCSFunctionOfHInter<'a> {
    surface: &'a dyn Surface,
    curve: &'a dyn Curve,
    p: GpPnt,
    f: f64,
}

impl<'a> TheCSFunctionOfHInter<'a> {
    /// `IntImp_ZerCSParFunc(S, C)` (`...gxx:24-30`).
    pub fn new(surface: &'a dyn Surface, curve: &'a dyn Curve) -> Self {
        Self {
            surface,
            curve,
            p: GpPnt::new(0.0, 0.0, 0.0),
            f: 0.0,
        }
    }

    /// `IntImp_ZerCSParFunc::Point` (`...gxx:97-100`).
    pub fn point(&self) -> &GpPnt {
        &self.p
    }

    /// `IntImp_ZerCSParFunc::Root` (`...gxx:102-105`): the squared distance
    /// between the surface point and the curve point.
    pub fn root(&self) -> f64 {
        self.f
    }

    /// `IntImp_ZerCSParFunc::AuxillarSurface` (`...gxx:107-110`).
    pub fn auxillar_surface(&self) -> &dyn Surface {
        self.surface
    }

    /// `IntImp_ZerCSParFunc::AuxillarCurve` (`...gxx:112-115`).
    pub fn auxillar_curve(&self) -> &dyn Curve {
        self.curve
    }
}

impl MathFunctionSetWithDerivatives for TheCSFunctionOfHInter<'_> {
    /// `IntImp_ZerCSParFunc::NbVariables` (`...gxx:31-34`).
    fn nb_variables(&self) -> usize {
        3
    }

    /// `IntImp_ZerCSParFunc::NbEquations` (`...gxx:36-39`).
    fn nb_equations(&self) -> usize {
        3
    }

    /// `IntImp_ZerCSParFunc::Value` (`...gxx:41-54`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let psurf = self.surface.d0(x.value(1), x.value(2));
        let pcurv = self.curve.d0(x.value(3));
        let f1 = psurf.x() - pcurv.x();
        let f2 = psurf.y() - pcurv.y();
        let f3 = psurf.z() - pcurv.z();
        f.set_value(1, f1);
        f.set_value(2, f2);
        f.set_value(3, f3);
        self.f = f1 * f1 + f2 * f2 + f3 * f3;
        self.p = GpPnt::new(
            0.5 * (psurf.x() + pcurv.x()),
            0.5 * (psurf.y() + pcurv.y()),
            0.5 * (psurf.z() + pcurv.z()),
        );
        true
    }

    /// `IntImp_ZerCSParFunc::Derivatives` (`...gxx:56-71`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (_, d1u, d1v) = self.surface.d1(x.value(1), x.value(2));
        let (_, d1w) = self.curve.d1(x.value(3));
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

    /// `IntImp_ZerCSParFunc::Values` (`...gxx:73-95`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (psurf, d1u, d1v) = self.surface.d1(x.value(1), x.value(2));
        let (pcurv, d1w) = self.curve.d1(x.value(3));
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
        self.f = f1 * f1 + f2 * f2 + f3 * f3;
        self.p = GpPnt::new(
            0.5 * (psurf.x() + pcurv.x()),
            0.5 * (psurf.y() + pcurv.y()),
            0.5 * (psurf.z() + pcurv.z()),
        );
        true
    }
}
