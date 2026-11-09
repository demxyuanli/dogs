//! `Contap_ArcFunction` - the 1-variable contour function on a restriction arc
//! of the surface: `F(u) = 0` with `(u, v) = Arc(u)`.
//!
//! Sources: `Contap_ArcFunction.hxx:28-82`, `.cxx:22-205`, `.lxx:16-64`.
//!
//! Deviation (handle -> value): `mySurf` / `myArc` are `occ::handle`s in OCCT;
//! the port stores the adaptors by value and answers `false` from the function
//! calls when either was never set, where OCCT dereferences a null handle.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_math::function_set_root::FunctionWithDerivative;

use occt_topo::brep_adaptor_curve2d::Curve2dAdaptor;
use occt_topo::int_surf_quadric::IntSurfQuadric;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::sample::{self, SampleGrid};
use super::surf_props;
use super::types::ContapTFunction;

/// `Contap_ArcFunction` (`hxx:28-82`).
pub struct ContapArcFunction {
    /// `myArc` (`hxx:70`).
    arc: Option<Arc<Curve2dAdaptor>>,
    /// `mySurf` (`hxx:71`).
    surf: Option<SurfaceAdaptor>,
    /// `myMean` (`hxx:72`).
    my_mean: f64,
    /// `myType` (`hxx:73`).
    my_type: ContapTFunction,
    /// `myDir` (`hxx:74`).
    my_dir: GpDir,
    /// `myCosAng` (`hxx:75`).
    my_cos_ang: f64,
    /// `myEye` (`hxx:76`).
    my_eye: GpPnt,
    /// `solpt` (`hxx:77`): the surface point of the last evaluation.
    solpt: GpPnt,
    /// `seqpt` (`hxx:78`): the points collected through `GetStateNumber`.
    seqpt: Vec<GpPnt>,
    /// `myQuad` (`hxx:79`): never set by `Contap_ArcFunction` itself, so it
    /// stays the default-constructed `IntSurf_Quadric` (`cxx:37-46`).
    my_quad: IntSurfQuadric,
}

impl Default for ContapArcFunction {
    /// `Contap_ArcFunction()` (`cxx:22-30`).
    fn default() -> Self {
        Self::new()
    }
}

impl ContapArcFunction {
    /// `Contap_ArcFunction()` (`cxx:22-30`): `myMean = 1`, `myDir = +Z`,
    /// `myCosAng = 0`.
    pub fn new() -> Self {
        Self {
            arc: None,
            surf: None,
            my_mean: 1.0,
            my_type: ContapTFunction::ContourStd,
            my_dir: GpDir::new(0.0, 0.0, 1.0).expect("+Z"),
            my_cos_ang: 0.0,
            my_eye: GpPnt::default(),
            solpt: GpPnt::default(),
            seqpt: Vec::new(),
            my_quad: IntSurfQuadric::default(),
        }
    }

    /// `Set(S)` (`cxx:32-55`): stores the surface and takes the mean normal
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
        self.surf = Some(s);
    }

    /// `Set(Direction)` (`lxx:32-36`): `Contap_ContourStd`.
    pub fn set_direction(&mut self, direction: GpDir) {
        self.my_type = ContapTFunction::ContourStd;
        self.my_dir = direction;
    }

    /// `Set(Direction, Angle)` (`lxx:17-23`): `Contap_DraftStd`.
    pub fn set_direction_angle(&mut self, direction: GpDir, angle: f64) {
        self.my_type = ContapTFunction::DraftStd;
        self.my_dir = direction;
        self.my_cos_ang = (std::f64::consts::FRAC_PI_2 + angle).cos();
    }

    /// `Set(Eye)` (`lxx:38-42`): `Contap_ContourPrs`.
    pub fn set_eye(&mut self, eye: GpPnt) {
        self.my_type = ContapTFunction::ContourPrs;
        self.my_eye = eye;
    }

    /// `Set(Eye, Angle)` (`lxx:25-30`): `Contap_DraftPrs`.
    pub fn set_eye_angle(&mut self, eye: GpPnt, angle: f64) {
        self.my_type = ContapTFunction::DraftPrs;
        self.my_eye = eye;
        self.my_cos_ang = (std::f64::consts::FRAC_PI_2 + angle).cos();
    }

    /// `Set(A)` (`lxx:44-48`).
    pub fn set_arc(&mut self, a: Arc<Curve2dAdaptor>) {
        self.arc = Some(a);
        self.seqpt.clear();
    }

    /// `NbSamples()` (`cxx:192-198`):
    /// `max(max(NbSamplesU, NbSamplesV), NbSamplesOnArc)`.
    pub fn nb_samples(&self) -> i32 {
        let (nu, nv) = match self.surf.as_ref() {
            Some(s) => (sample::nb_samples_u(s), sample::nb_samples_v(s)),
            None => (0, 0),
        };
        let na = self
            .arc
            .as_ref()
            .and_then(|a| a.curve())
            .map(|c| sample::nb_samples_on_arc(c.as_ref()))
            .unwrap_or(0);
        nu.max(nv).max(na)
    }

    /// `Valpoint(Index)` (`lxx:50-53`), 1-based.
    pub fn valpoint(&self, index: usize) -> Option<&GpPnt> {
        index.checked_sub(1).and_then(|i| self.seqpt.get(i))
    }

    /// `Quadric()` (`cxx:204-208`).
    pub fn quadric(&self) -> &IntSurfQuadric {
        &self.my_quad
    }

    /// `Surface()` (`lxx:55-58`).
    pub fn surface(&self) -> Option<&SurfaceAdaptor> {
        self.surf.as_ref()
    }

    /// `LastComputedPoint()` (`lxx:60-63`).
    pub fn last_computed_point(&self) -> &GpPnt {
        &self.solpt
    }

    /// The `(u, v)` point of the arc at `u` (`Contap_HCurve2dTool::Value`,
    /// `cxx:58`).
    fn arc_point(&self, u: f64) -> Option<occt_core::gp::GpPnt2d> {
        self.arc.as_ref().and_then(|a| sample::arc_value(a, u))
    }

    /// The value of `F` and the surface point it was evaluated at.
    fn eval_value(&mut self, u: f64) -> Option<f64> {
        let pt2d = self.arc_point(u)?;
        let s = self.surf.as_ref()?;
        let (p, norm) = surf_props::normale(s, pt2d.x(), pt2d.y())?;
        self.solpt = p;
        let f = match self.my_type {
            ContapTFunction::ContourStd => norm.dot(&GpVec::from_xyz(self.my_dir.xyz())),
            ContapTFunction::ContourPrs => {
                norm.dot(&GpVec::from_pnts(&self.my_eye, &self.solpt))
            }
            ContapTFunction::DraftStd => {
                norm.dot(&GpVec::from_xyz(self.my_dir.xyz()))
                    - self.my_cos_ang * norm.magnitude()
            }
            // `cxx:79-82`: the value is left untouched.
            ContapTFunction::DraftPrs => return Some(0.0),
        };
        Some(f / self.my_mean)
    }

    /// The `F` / `dF/du` pair of `Value` / `Derivative` (`cxx:57-133`).
    fn eval_values(&mut self, u: f64) -> Option<(f64, f64)> {
        let (pt2d, d2d) = self.arc.as_ref().and_then(|a| sample::arc_d1(a, u))?;
        let s = self.surf.as_ref()?;
        let (p, norm, dnu, dnv) = surf_props::norm_and_dn(s, pt2d.x(), pt2d.y())?;
        self.solpt = p;
        let dir = GpVec::from_xyz(self.my_dir.xyz());
        let (f, dfu, dfv) = match self.my_type {
            ContapTFunction::ContourStd => (
                norm.dot(&dir),
                dnu.dot(&dir),
                dnv.dot(&dir),
            ),
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
            // `cxx:79-82`, `:177-180`: nothing is written.
            ContapTFunction::DraftPrs => return Some((0.0, 0.0)),
        };
        let d = d2d.x() * dfu / self.my_mean + d2d.y() * dfv / self.my_mean;
        Some((f / self.my_mean, d))
    }
}

impl FunctionWithDerivative for ContapArcFunction {
    /// `Value(U, F)` (`cxx:57-84`).
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        match self.eval_value(x) {
            Some(v) => {
                *f = v;
                true
            }
            None => false,
        }
    }

    /// `Derivative(U, D)` (`cxx:86-133`).
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        match self.eval_values(x) {
            Some((_, dv)) => {
                *d = dv;
                true
            }
            None => false,
        }
    }

    /// `Values(U, F, D)` (`cxx:135-185`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        match self.eval_values(x) {
            Some((v, dv)) => {
                *f = v;
                *d = dv;
                true
            }
            None => false,
        }
    }

    /// `GetStateNumber()` (`cxx:187-191`): appends the last computed point and
    /// returns the new sequence length.
    fn state_number(&mut self) -> i32 {
        self.seqpt.push(self.solpt);
        self.seqpt.len() as i32
    }
}

impl ContapArcFunction {
    /// `Value(U, F)` for callers that hold the concrete type.
    pub fn eval(&mut self, x: f64, f: &mut f64) -> bool {
        <Self as FunctionWithDerivative>::value(self, x, f)
    }

    /// `Values(U, F, D)` for callers that hold the concrete type.
    pub fn eval_f_and_d(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        <Self as FunctionWithDerivative>::values(self, x, f, d)
    }

    /// `GetStateNumber()` for callers that hold the concrete type (it appends
    /// the last computed point and answers the new sequence length, 1-based).
    pub fn get_state_number(&mut self) -> i32 {
        <Self as FunctionWithDerivative>::state_number(self)
    }
}

impl occt_core::math_fn::MathFunction for ContapArcFunction {
    /// `Value(X, F)` of `math_Function` (`Contap_ArcFunction.cxx:57-84`).
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        <Self as FunctionWithDerivative>::value(self, x, f)
    }

    /// `GetStateNumber()` (`Contap_ArcFunction.cxx:187-191`).
    fn state_number(&mut self) -> i32 {
        <Self as FunctionWithDerivative>::state_number(self)
    }
}

impl occt_core::math_fn::MathFunctionWithDerivative for ContapArcFunction {
    /// `Derivative(X, D)` (`Contap_ArcFunction.cxx:86-133`).
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        <Self as FunctionWithDerivative>::derivative(self, x, d)
    }

    /// `Values(X, F, D)` (`Contap_ArcFunction.cxx:135-185`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        <Self as FunctionWithDerivative>::values(self, x, f, d)
    }
}

