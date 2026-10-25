//! `BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox` — the scalar
//! function "distance from the parametric surface point to the quadric" whose
//! zero is a point of the intersection line.
//!
//! Sources:
//! - `ModelingAlgorithms/TKGeomAlgo/IntImp/IntImp_ZerImpFunc.gxx:25-…` (the
//!   body) and `IntImp_ZerImpFunc.lxx:19-84` (the accessors)
//! - `BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox_0.cxx:16-30`: the
//!   `#define ThePSurface BRepAdaptor_Surface / ThePSurfaceTool
//!   BRepApprox_SurfaceTool / TheISurface IntSurf_Quadric / TheISurfaceTool
//!   IntSurf_QuadricTool` binding
//! - `BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox.hxx:28-78`
//!
//! Translation notes:
//! - OCCT stores the two surfaces as `void*` (`SURF` / `FUNC`) and reads them
//!   back as references (`lxx:62-84`). The port stores owned copies
//!   ([`SurfaceAdaptor`] is `Clone`, [`IntSurfQuadric`] is `Clone`): the only
//!   construction site is `ApproxInt_ImpPrmSvSurfaces`'s constructors, which
//!   pass its own surfaces and never mutate them afterwards, so a copy is
//!   observationally the same and avoids a self-referential struct.
//! - `EpsAng` (`gxx:16`) is unused by the body; `EpsAng2` (`:17`) and
//!   `Tolpetit` (`:18`) are the two constants used by `IsTangent`.
//! - `Value` does **not** refresh `gradient`, `Derivatives` / `Values` do
//!   (`gxx:66-108`); `IsTangent` reads whatever the last call left, as in OCCT.
//! - `Direction3d` / `Direction2d` (`lxx:48-60`) throw
//!   `StdFail_UndefinedDerivative` when `IsTangent()`; the port returns
//!   `Option` (`None` on that branch) because `IsTangent` mutates the cached
//!   state.

use occt_core::gp::{GpDir2d, GpPnt, GpVec};
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::brep_approx::SurfaceTool;
use crate::int_surf_quadric::IntSurfQuadric;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `#define EpsAng2 1.e-16` (`gxx:17`).
const EPS_ANG2: f64 = 1.0e-16;
/// `#define Tolpetit 1.e-16` (`gxx:18`).
const TOL_PETIT: f64 = 1.0e-16;

/// `BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox` (`hxx:28`).
#[derive(Clone)]
pub struct ZerImpFunc {
    /// `surf` (`gxx:32`): `ThePSurface`, the parametric half.
    surf: SurfaceAdaptor,
    /// `func` (`gxx:33`): `TheISurface`, the implicit half.
    func: IntSurfQuadric,
    /// `u`.
    u: f64,
    /// `v`.
    v: f64,
    /// `tol`.
    tol: f64,
    /// `pntsol`.
    pntsol: GpPnt,
    /// `valf`.
    valf: f64,
    /// `computed`.
    computed: bool,
    /// `tangent`.
    tangent: bool,
    /// `tgdu`.
    tgdu: f64,
    /// `tgdv`.
    tgdv: f64,
    /// `derived`.
    derived: bool,
    /// `gradient`.
    gradient: GpVec,
    /// `d1u`.
    d1u: GpVec,
    /// `d1v`.
    d1v: GpVec,
    /// `d3d`.
    d3d: GpVec,
    /// `d2d`.
    d2d: GpDir2d,
}

impl ZerImpFunc {
    /// `IntImp_ZerImpFunc(const ThePSurface& PS, const TheISurface& IS)`
    /// (`gxx:44-60`).
    pub fn new(ps: SurfaceAdaptor, is: IntSurfQuadric) -> Self {
        Self {
            surf: ps,
            func: is,
            u: 0.0,
            v: 0.0,
            tol: 0.0,
            pntsol: GpPnt::zero(),
            valf: 0.0,
            computed: false,
            tangent: false,
            tgdu: 0.0,
            tgdv: 0.0,
            derived: false,
            gradient: GpVec::zero(),
            d1u: GpVec::zero(),
            d1v: GpVec::zero(),
            d3d: GpVec::zero(),
            d2d: GpDir2d::default(),
        }
    }

    /// `Root()` (`lxx:40-43`).
    pub fn root(&self) -> f64 {
        self.valf
    }

    /// `Tolerance()` (`lxx:45-48`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `Point()` (`lxx:50-53`).
    pub fn point(&self) -> GpPnt {
        self.pntsol
    }

    /// `Direction3d()` (`lxx:55-60`): `StdFail_UndefinedDerivative` when the
    /// point is tangent.
    pub fn direction_3d(&mut self) -> Option<GpVec> {
        if self.is_tangent() {
            return None;
        }
        Some(self.d3d)
    }

    /// `Direction2d()` (`lxx:56-61`), the 2D mirror of
    /// [`ZerImpFunc::direction_3d`].
    pub fn direction_2d(&mut self) -> Option<GpDir2d> {
        if self.is_tangent() {
            return None;
        }
        Some(self.d2d)
    }

    /// `PSurface()` (`lxx:63-66`).
    pub fn psurface(&self) -> &SurfaceAdaptor {
        &self.surf
    }

    /// `ISurface()` (`lxx:68-71`).
    pub fn isurface(&self) -> &IntSurfQuadric {
        &self.func
    }

    /// `IsTangent()` (`gxx:110-141`).
    pub fn is_tangent(&mut self) -> bool {
        if !self.computed {
            self.computed = true;
            if !self.derived {
                SurfaceTool::d1(
                    &self.surf,
                    self.u,
                    self.v,
                    &mut self.pntsol,
                    &mut self.d1u,
                    &mut self.d1v,
                );
                self.derived = true;
            }

            self.tgdu = self.gradient.dot(&self.d1v);
            self.tgdv = -self.gradient.dot(&self.d1u);
            let n2grad = self.gradient.square_magnitude();
            let n2grad_eps_ang2 = n2grad * EPS_ANG2;
            let n2d1u = self.d1u.square_magnitude();
            let n2d1v = self.d1v.square_magnitude();
            self.tangent = (self.tgdu * self.tgdu <= n2grad_eps_ang2 * n2d1v)
                && (self.tgdv * self.tgdv <= n2grad_eps_ang2 * n2d1u);
            if !self.tangent {
                self.d3d = self
                    .d3d
                    .set_linear_form_2(self.tgdu, &self.d1u, self.tgdv, &self.d1v);
                self.d2d = GpDir2d::new(self.tgdu, self.tgdv)
                    .expect("gp_Dir2d: IntImp_ZerImpFunc.gxx:133");
                if self.d3d.magnitude() <= TOL_PETIT {
                    // jag
                    self.tangent = true;
                }
            }
        }
        self.tangent
    }
}

impl MathFunctionSetWithDerivatives for ZerImpFunc {
    /// `NbVariables()` (`gxx:84-87`).
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`gxx:89-92`).
    fn nb_equations(&self) -> usize {
        1
    }

    /// `Value(X, F)` (`gxx:94-104`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.u = x.value(1);
        self.v = x.value(2);
        self.pntsol = SurfaceTool::value(&self.surf, self.u, self.v);
        self.valf = self.func.distance(&self.pntsol);
        f.set_value(1, self.valf);
        self.computed = false;
        self.derived = false;
        true
    }

    /// `Derivatives(X, D)` (`gxx:106-118`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.u = x.value(1);
        self.v = x.value(2);
        SurfaceTool::d1(
            &self.surf,
            self.u,
            self.v,
            &mut self.pntsol,
            &mut self.d1u,
            &mut self.d1v,
        );
        self.gradient = self.func.gradient(&self.pntsol);
        d.set_value(1, 1, self.d1u.dot(&self.gradient));
        d.set_value(1, 2, self.d1v.dot(&self.gradient));
        self.computed = false;
        self.derived = true;
        true
    }

    /// `Values(X, F, D)` (`gxx:120-134`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.u = x.value(1);
        self.v = x.value(2);
        SurfaceTool::d1(
            &self.surf,
            self.u,
            self.v,
            &mut self.pntsol,
            &mut self.d1u,
            &mut self.d1v,
        );
        let (val, grad) = self.func.val_and_grad(&self.pntsol);
        self.valf = val;
        self.gradient = grad;
        f.set_value(1, self.valf);
        d.set_value(1, 1, self.d1u.dot(&self.gradient));
        d.set_value(1, 2, self.d1v.dot(&self.gradient));
        self.computed = false;
        self.derived = true;
        true
    }
}
