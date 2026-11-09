//! `Contap_SurfFunction` - the 2-variable contour function `F(u, v) = 0` on a
//! parametric surface: the derivative of the surface normal along the
//! projection direction, normalised by the mean normal magnitude.
//!
//! Sources: `Contap_SurfFunction.hxx:30-118`, `.cxx:29-285`, `.lxx:19-107`.
//!
//! Deviation (handle -> value): `mySurf` is an `occ::handle` in OCCT; the port
//! stores the adaptor by value and answers `false` from `Value` /
//! `Derivatives` / `Values` when it was never set, where OCCT dereferences a
//! null handle.

use occt_core::gp::{GpDir, GpDir2d, GpPnt, GpVec};
use occt_core::precision::RESOLUTION;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::sample::SampleGrid;
use super::surf_props;
use super::types::ContapTFunction;

/// `Contap_SurfFunction` (`hxx:30-118`).
#[derive(Clone)]
pub struct ContapSurfFunction {
    /// `mySurf` (`hxx:99`).
    surf: Option<SurfaceAdaptor>,
    /// `myMean` (`hxx:100`).
    my_mean: f64,
    /// `myType` (`hxx:101`).
    my_type: ContapTFunction,
    /// `myDir` (`hxx:102`).
    my_dir: GpDir,
    /// `myEye` (`hxx:103`).
    my_eye: GpPnt,
    /// `myAng` (`hxx:105`).
    my_ang: f64,
    /// `myCosAng` (`hxx:106`).
    my_cos_ang: f64,
    /// `tol` (`hxx:107`).
    tol: f64,
    /// `solpt` (`hxx:108`).
    solpt: GpPnt,
    /// `valf` (`hxx:109`).
    valf: f64,
    /// `Usol` (`hxx:110`).
    u_sol: f64,
    /// `Vsol` (`hxx:111`).
    v_sol: f64,
    /// `Fpu` (`hxx:112`).
    fpu: f64,
    /// `Fpv` (`hxx:113`).
    fpv: f64,
    /// `d2d` (`hxx:114`).
    d2d: Option<GpDir2d>,
    /// `d3d` (`hxx:115`).
    d3d: GpVec,
    /// `tangent` (`hxx:116`).
    tangent: bool,
    /// `computed` (`hxx:117`).
    computed: bool,
    /// `derived` (`hxx:118`).
    derived: bool,
}

impl Default for ContapSurfFunction {
    /// `Contap_SurfFunction()` (`cxx:29-45`).
    fn default() -> Self {
        Self::new()
    }
}

impl ContapSurfFunction {
    /// `Contap_SurfFunction()` (`cxx:29-45`): `myMean = 1`, `myType =
    /// Contap_ContourStd`, `myDir = +Z`, `myAng = myCosAng = 0`, `tol = 1e-6`.
    pub fn new() -> Self {
        Self {
            surf: None,
            my_mean: 1.0,
            my_type: ContapTFunction::ContourStd,
            my_dir: GpDir::new(0.0, 0.0, 1.0).expect("+Z"),
            my_eye: GpPnt::default(),
            my_ang: 0.0,
            my_cos_ang: 0.0,
            tol: 1.0e-6,
            solpt: GpPnt::default(),
            valf: 0.0,
            u_sol: 0.0,
            v_sol: 0.0,
            fpu: 0.0,
            fpv: 0.0,
            d2d: None,
            d3d: GpVec::new(0.0, 0.0, 0.0),
            tangent: false,
            computed: false,
            derived: false,
        }
    }

    /// `Set(S)` (`cxx:47-69`): stores the surface and takes the mean normal
    /// magnitude over the sample points.
    pub fn set_surface(&mut self, s: SurfaceAdaptor) {
        let (grid, nbs) = SampleGrid::new(&s);
        if nbs > 0 {
            self.my_mean = 0.0;
            for i in 1..=nbs {
                let (u, v) = grid.sample_point(&s, i);
                if let Some((p, norm)) = surf_props::normale(&s, u, v) {
                    self.solpt = p;
                    self.my_mean += norm.magnitude();
                }
            }
            self.my_mean /= f64::from(nbs);
        }
        self.computed = false;
        self.derived = false;
        self.surf = Some(s);
    }

    /// `Set(Eye)` (`lxx:19-23`): `Contap_ContourPrs`.
    pub fn set_eye(&mut self, eye: GpPnt) {
        self.my_type = ContapTFunction::ContourPrs;
        self.my_eye = eye;
        self.my_ang = 0.0;
    }

    /// `Set(Direction)` (`lxx:25-30`): `Contap_ContourStd`.
    pub fn set_direction(&mut self, direction: GpDir) {
        self.my_type = ContapTFunction::ContourStd;
        self.my_dir = direction;
        self.my_ang = 0.0;
    }

    /// `Set(Direction, Angle)` (`lxx:32-39`): `Contap_DraftStd`.
    pub fn set_direction_angle(&mut self, direction: GpDir, angle: f64) {
        self.my_type = ContapTFunction::DraftStd;
        self.my_dir = direction;
        self.my_ang = angle;
        self.my_cos_ang = (std::f64::consts::FRAC_PI_2 + angle).cos();
    }

    /// `Set(Eye, Angle)` (`lxx:41-48`): `Contap_DraftPrs`.
    pub fn set_eye_angle(&mut self, eye: GpPnt, angle: f64) {
        self.my_type = ContapTFunction::DraftPrs;
        self.my_eye = eye;
        self.my_ang = angle;
        self.my_cos_ang = (std::f64::consts::FRAC_PI_2 + angle).cos();
    }

    /// `Set(Tolerance)` (`lxx:50-53`): `tol = max(Tolerance, 1e-12)`.
    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.tol = tolerance.max(1.0e-12);
    }

    /// `Root()` (`lxx:58-61`).
    pub fn root(&self) -> f64 {
        self.valf
    }

    /// `Tolerance()` (`lxx:63-66`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `Point()` (`lxx:55-58`).
    pub fn point(&self) -> &GpPnt {
        &self.solpt
    }

    /// `FunctionType()` (`lxx:101-104`).
    pub fn function_type(&self) -> ContapTFunction {
        self.my_type
    }

    /// `Eye()` (`lxx:88-91`).
    pub fn eye(&self) -> &GpPnt {
        &self.my_eye
    }

    /// `Direction()` (`lxx:93-96`).
    pub fn direction(&self) -> &GpDir {
        &self.my_dir
    }

    /// `Angle()` (`lxx:98-101`).
    pub fn angle(&self) -> f64 {
        self.my_ang
    }

    /// `Surface()` (`lxx:81-84`).
    pub fn surface(&self) -> Option<&SurfaceAdaptor> {
        self.surf.as_ref()
    }

    /// `Direction3d()` (`lxx:68-73`): `StdFail_UndefinedDerivative` when
    /// `IsTangent()`.
    pub fn direction3d(&mut self) -> Result<&GpVec, &'static str> {
        if self.is_tangent() {
            return Err("Contap_SurfFunction::Direction3d: StdFail_UndefinedDerivative");
        }
        Ok(&self.d3d)
    }

    /// `Direction2d()` (`lxx:75-80`).
    pub fn direction2d(&mut self) -> Result<GpDir2d, &'static str> {
        if self.is_tangent() {
            return Err("Contap_SurfFunction::Direction2d: StdFail_UndefinedDerivative");
        }
        self.d2d
            .ok_or("Contap_SurfFunction::Direction2d: StdFail_UndefinedDerivative")
    }

    /// `IsTangent()` (`cxx:214-284`): the cached tangency test, computing the
    /// first derivatives on the way when `Derivatives` was not called.
    pub fn is_tangent(&mut self) -> bool {
        if self.computed {
            return self.tangent;
        }
        self.computed = true;

        if !self.derived {
            if let Some(s) = self.surf.as_ref() {
                if let Some((_, norm, dnu, dnv)) =
                    surf_props::norm_and_dn(s, self.u_sol, self.v_sol)
                {
                    let dir = GpVec::from_xyz(self.my_dir.xyz());
                    match self.my_type {
                        ContapTFunction::ContourStd => {
                            self.fpu = dnu.dot(&dir) / self.my_mean;
                            self.fpv = dnv.dot(&dir) / self.my_mean;
                        }
                        ContapTFunction::ContourPrs => {
                            let ep = GpVec::from_pnts(&self.my_eye, &self.solpt);
                            self.fpu = dnu.dot(&ep) / self.my_mean;
                            self.fpv = dnv.dot(&ep) / self.my_mean;
                        }
                        ContapTFunction::DraftStd => {
                            let mut n = norm;
                            n.normalize();
                            self.fpu =
                                (dnu.dot(&dir) - self.my_cos_ang * dnu.dot(&n)) / self.my_mean;
                            self.fpv =
                                (dnv.dot(&dir) - self.my_cos_ang * dnv.dot(&n)) / self.my_mean;
                        }
                        // `cxx:257-259`: nothing is written.
                        ContapTFunction::DraftPrs => {}
                    }
                }
            }
            self.derived = true;
        }

        self.tangent = false;
        let d = (self.fpu * self.fpu + self.fpv * self.fpv).sqrt();
        if d <= RESOLUTION {
            self.tangent = true;
        } else {
            self.d2d = GpDir2d::new(-self.fpv, self.fpu).ok();
            // `Adaptor3d_HSurfaceTool::D1(mySurf, Usol, Vsol, solpt, d1u, d1v)`
            // (`cxx:271`).
            let (d1u, d1v) = match self.surf.as_ref().and_then(|s| s.d1(self.u_sol, self.v_sol)) {
                Some((_, d1u, d1v)) => (d1u, d1v),
                None => {
                    self.tangent = true;
                    return self.tangent;
                }
            };
            self.d3d = d1v
                .multiplied_scalar(self.fpu)
                .added(&d1u.multiplied_scalar(-self.fpv));
            if self.d3d.magnitude() <= self.tol {
                self.tangent = true;
            }
        }
        self.tangent
    }

    /// The value arm of `Value` / `Values` (`cxx:81-109`, `:169-204`).
    fn eval_value(&mut self, u: f64, v: f64) -> Option<f64> {
        self.u_sol = u;
        self.v_sol = v;
        let s = self.surf.as_ref()?;
        let (p, norm) = surf_props::normale(s, self.u_sol, self.v_sol)?;
        self.solpt = p;
        let dir = GpVec::from_xyz(self.my_dir.xyz());
        let f = match self.my_type {
            ContapTFunction::ContourStd => norm.dot(&dir),
            ContapTFunction::ContourPrs => {
                norm.dot(&GpVec::from_pnts(&self.my_eye, &self.solpt))
            }
            ContapTFunction::DraftStd => norm.dot(&dir) - self.my_cos_ang * norm.magnitude(),
            // `cxx:102-105`: the value is left untouched.
            ContapTFunction::DraftPrs => 0.0,
        };
        Some(f / self.my_mean)
    }

    /// The gradient arm of `Derivatives` / `Values` (`cxx:112-156`,
    /// `:158-212`).
    fn eval_derivatives(&mut self, u: f64, v: f64) -> Option<(f64, f64, f64)> {
        self.u_sol = u;
        self.v_sol = v;
        let s = self.surf.as_ref()?;
        let (p, norm, dnu, dnv) = surf_props::norm_and_dn(s, self.u_sol, self.v_sol)?;
        self.solpt = p;
        let dir = GpVec::from_xyz(self.my_dir.xyz());
        let (f, dfu, dfv) = match self.my_type {
            ContapTFunction::ContourStd => (norm.dot(&dir), dnu.dot(&dir), dnv.dot(&dir)),
            ContapTFunction::ContourPrs => {
                let ep = GpVec::from_pnts(&self.my_eye, &self.solpt);
                (norm.dot(&ep), dnu.dot(&ep), dnv.dot(&ep))
            }
            ContapTFunction::DraftStd => {
                let f = norm.dot(&dir) - self.my_cos_ang * norm.magnitude();
                let mut n = norm;
                n.normalize();
                (
                    f,
                    dnu.dot(&dir) - self.my_cos_ang * dnu.dot(&n),
                    dnv.dot(&dir) - self.my_cos_ang * dnv.dot(&n),
                )
            }
            // `cxx:149-151`, `:201-204`: nothing is written.
            ContapTFunction::DraftPrs => (0.0, 0.0, 0.0),
        };
        Some((f / self.my_mean, dfu / self.my_mean, dfv / self.my_mean))
    }
}

impl MathFunctionSetWithDerivatives for ContapSurfFunction {
    /// `NbVariables()` (`cxx:71-74`): 2.
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`cxx:76-79`): 1.
    fn nb_equations(&self) -> usize {
        1
    }

    /// `Value(X, F)` (`cxx:81-109`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let v = match self.eval_value(x.value(1), x.value(2)) {
            Some(v) => v,
            None => return false,
        };
        self.valf = v;
        f.set_value(1, v);
        self.computed = false;
        self.derived = false;
        true
    }

    /// `Derivatives(X, Grad)` (`cxx:111-156`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (_, dfu, dfv) = match self.eval_derivatives(x.value(1), x.value(2)) {
            Some(v) => v,
            None => return false,
        };
        d.set_value(1, 1, dfu);
        d.set_value(1, 2, dfv);
        self.fpu = dfu;
        self.fpv = dfv;
        self.computed = false;
        self.derived = true;
        true
    }

    /// `Values(X, F, Grad)` (`cxx:158-212`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (v, dfu, dfv) = match self.eval_derivatives(x.value(1), x.value(2)) {
            Some(v) => v,
            None => return false,
        };
        f.set_value(1, v);
        d.set_value(1, 1, dfu);
        d.set_value(1, 2, dfv);
        self.valf = v;
        self.fpu = dfu;
        self.fpv = dfv;
        self.computed = false;
        self.derived = true;
        true
    }
}

