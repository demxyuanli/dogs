//! `BlendFunc_CSCircular` — the circular blending function between a curve and
//! a surface, driven by a guide line and a parameter law.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSCircular.hxx:36-237`,
//! `BlendFunc_CSCircular.cxx:42-1029` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(U, V)`, the
//! parameters of the point on the surface. With `nplan` the unit tangent of the
//! guide line at the current parameter, `prmc` / `dprmc` the value and the
//! derivative of the law at that same parameter, `ns` the surface normal
//! projected into the section plane, `ray` the signed radius and
//! `vref = ray * ns + (pts - ptc)`, the two equations are
//! 1. `nplan . (pts - ptc)` — the two points share a section plane normal to
//!    the guide line,
//! 2. `|vref|^2 - ray^2` — the curve point is at the distance `|ray|` from the
//!    centre of the arc.
//!
//! `BlendFunc_CSConstRad` is the same system with the curve parameter as the
//! third unknown and a constant radius; here the curve parameter is imposed by
//! the law, which is why `Set(Param)` (and not `Value`) evaluates it.
//!
//! Translation notes:
//! - `Eps` (`cxx:39`) guards the projected normal in `Value` / `Derivatives` /
//!   `Values` / `IsSolution` (`if (norm < Eps) norm = 1.;`, `cxx:263-266`,
//!   `:297-300`, `:353-356`, `:178-181`) but **not** in `GetSection`
//!   (`cxx:593-595`), `Section(Param, ...)` (`cxx:511-513`),
//!   `Section(P, Poles, ...)` (`cxx:822-824`) nor `Tangent` (`cxx:475`); the
//!   difference is kept.
//! - `set_radius` keeps the sign rule of `Set(Radius, Choix)` (`cxx:77-90`):
//!   quadrants 3 and 4 give a positive radius, every other choice a negative
//!   one.
//! - `Set(First, Last)` throws `Standard_NotImplemented` (`cxx:114-117`); the
//!   port panics with the same message.
//! - `Section(P, ..., D2Poles, ...)` (`cxx:533-554`) forwards to
//!   `Blend_CSFunction::Section` (`Blend_CSFunction.cxx:35-47`), which answers
//!   `false`; the port does the same instead of panicking.
//! - `TangentOnS` / `TangentOnC` / `Tangent2d` throw `Standard_DomainError`
//!   when `istangent` (`cxx:436-464`); the trait returns a value and every OCCT
//!   caller guards the call with `!IsTangencyPoint()`, so the stored value is
//!   returned.
//! - `GetSection` (`cxx:557-707`) and the analytic `Section(Param, ...)`
//!   (`cxx:490-531`) are kept as inherent methods.
//! - `myTConv` (`hxx:236`) is not in the constructor init list (`cxx:44-59`),
//!   so OCCT leaves it indeterminate until the first `BlendFunc::GetShape`
//!   writes it through its `TConv&` parameter (`BlendFunc.cxx:41-87`). The port
//!   starts it at `TgtThetaOver2` to stay deterministic; `get_shape` overwrites
//!   it before `Knots` / `Mults` / `GetTolerance` read it.
//! - `GeomFill::GetCircle` answers `()` here (OCCT answers `true`); the section
//!   routines therefore do not forward it.

pub(super) mod constants;
mod equations;
mod section;
mod solution;

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::law::LawFunction;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendCsFunction, BlendFuncSectionShape, BlendPoint};

/// `class BlendFunc_CSCircular : public Blend_CSFunction` (`hxx:36-237`).
pub struct CsCircular {
    /// `surf` (`hxx:213`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:214`).
    curv: Arc<dyn Curve>,
    /// `guide` (`hxx:215`).
    guide: Arc<dyn Curve>,
    /// `law` (`hxx:216`): the change of parameter between `curv` and `guide`.
    law: Arc<dyn LawFunction>,
    /// `pts` (`hxx:217`): the point on the surface.
    pts: GpPnt,
    /// `ptc` (`hxx:218`): the point on the curve.
    ptc: GpPnt,
    /// `pt2d` (`hxx:219`).
    pt2d: GpPnt2d,
    /// `prmc` (`hxx:220`): the parameter on the curve, imposed by the law.
    prmc: f64,
    /// `dprmc` (`hxx:221`): the derivative of the law at the guide parameter.
    dprmc: f64,
    /// `istangent` (`hxx:222`).
    istangent: bool,
    /// `tgs` (`hxx:223`).
    tgs: GpVec,
    /// `tg2d` (`hxx:224`).
    tg2d: GpVec2d,
    /// `tgc` (`hxx:225`).
    tgc: GpVec,
    /// `ray` (`hxx:226`): the signed radius of the section.
    ray: f64,
    /// `choix` (`hxx:227`): the fillet quadrant.
    choix: i32,
    /// `d1gui` (`hxx:228`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:229`).
    d2gui: GpVec,
    /// `nplan` (`hxx:230`): the unit tangent of the guide line.
    nplan: GpVec,
    /// `normtg` (`hxx:231`).
    normtg: f64,
    /// `maxang` (`hxx:232`): the widest section angle seen by `IsSolution`.
    maxang: f64,
    /// `minang` (`hxx:233`): the narrowest section angle.
    minang: f64,
    /// `mySShape` (`hxx:234`).
    my_sshape: BlendFuncSectionShape,
    /// `myTConv` (`hxx:235`); see the module note on its initialisation.
    my_tconv: ParameterisationType,
    /// Scratch storage for `Value` / `Derivatives` / `Values`: `hxx` declares
    /// no `E` / `DEDX` for this class (OCCT uses local `gp_Vec`s, and
    /// `IsSolution` its own `valsol` / `gradsol`, `cxx:154-155`), so this is
    /// only the port's output buffer, never read across calls.
    e: MathVector,
    /// Output buffer for the jacobian; see `e`.
    dedx: MathMatrix,
}

impl CsCircular {
    /// `BlendFunc_CSCircular(S, C, CG, L)` (`cxx:44-59`). The constructor also
    /// assigns `law = L` a second time in its body; the port keeps the single
    /// assignment.
    pub fn new(
        surf: Arc<dyn Surface>,
        curv: Arc<dyn Curve>,
        guide: Arc<dyn Curve>,
        law: Arc<dyn LawFunction>,
    ) -> Self {
        Self {
            surf,
            curv,
            guide,
            law,
            pts: GpPnt::zero(),
            ptc: GpPnt::zero(),
            pt2d: GpPnt2d::new(0.0, 0.0),
            prmc: 0.0,
            dprmc: 0.0,
            istangent: true,
            tgs: GpVec::zero(),
            tg2d: GpVec2d::new(0.0, 0.0),
            tgc: GpVec::zero(),
            ray: 0.0,
            choix: 0,
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            maxang: f64::MIN,
            minang: f64::MAX,
            my_sshape: BlendFuncSectionShape::Rational,
            my_tconv: ParameterisationType::TgtThetaOver2,
            e: MathVector::new(1, 2),
            dedx: MathMatrix::new(1, 2, 1, 2),
        }
    }

    /// `Set(Radius, Choix)` (`cxx:77-90`): the radius is positive for the
    /// quadrants 3 and 4, negative otherwise.
    pub fn set_radius(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        self.ray = match choix {
            3 | 4 => radius.abs(),
            _ => -radius.abs(),
        };
    }

    /// `Set(TypeSection)` (`cxx:95-98`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `Set(Param)` (`cxx:102-109`): the section plane and the curve parameter
    /// imposed by the law. `ptgui` is local in OCCT and stays local here.
    pub fn set_param(&mut self, param: f64) {
        let (_ptgui, d1gui, d2gui) = self.guide.d2(param);
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        let (prmc, dprmc) = self.law.d1(param);
        self.prmc = prmc;
        self.dprmc = dprmc;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
    }

    /// `ray` (`hxx:226`).
    pub fn ray(&self) -> f64 {
        self.ray
    }

    /// `choix` (`hxx:227`).
    pub fn choix(&self) -> i32 {
        self.choix
    }

    /// `maxang` (`hxx:232`).
    pub fn max_ang(&self) -> f64 {
        self.maxang
    }

    /// `minang` (`hxx:233`).
    pub fn min_ang(&self) -> f64 {
        self.minang
    }
}

impl MathFunctionSetWithDerivatives for CsCircular {
    /// `NbVariables()` (`cxx:64-67`): `2`, the surface parameters; the curve
    /// parameter comes from the law.
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`cxx:72-75`): `2`.
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`cxx:255-277`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let ok = self.evaluate_value(x);
        f.initialized(&self.e);
        ok
    }

    /// `Derivatives(X, D)` (`cxx:282-337`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let ok = self.evaluate_derivatives(x);
        for i in 1..=2 {
            for j in 1..=2 {
                d.set_value(i, j, self.dedx.value(i, j));
            }
        }
        ok
    }

    /// `Values(X, F, D)` (`cxx:339-397`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let ok = self.evaluate_values(x);
        f.initialized(&self.e);
        for i in 1..=2 {
            for j in 1..=2 {
                d.set_value(i, j, self.dedx.value(i, j));
            }
        }
        ok
    }
}

impl BlendAppFunction for CsCircular {
    /// `Set(Param)` (`cxx:102-109`).
    fn set(&mut self, param: f64) {
        self.set_param(param);
    }

    /// `Set(First, Last)` (`cxx:114-117`) throws `Standard_NotImplemented`.
    fn set_interval(&mut self, _first: f64, _last: f64) {
        panic!("Standard_NotImplemented: BlendFunc_CSCircular::Set");
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:120-125`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:130-148`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:152-248`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()`: `Blend_CSFunction::GetMinimalDistance` throws
    /// `Standard_NotImplemented` (`Blend_CSFunction.cxx:49-52`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Err("Standard_NotImplemented: Blend_CSFunction::GetMinimalDistance".into())
    }

    /// `Pnt1()`: `Blend_CSFunction::Pnt1()` = `PointOnC()`
    /// (`Blend_CSFunction.cxx:25-28`).
    fn pnt1(&self) -> &GpPnt {
        &self.ptc
    }

    /// `Pnt2()`: `Blend_CSFunction::Pnt2()` = `PointOnS()`
    /// (`Blend_CSFunction.cxx:30-33`).
    fn pnt2(&self) -> &GpPnt {
        &self.pts
    }

    /// `IsRational()` (`cxx:713-716`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:719-722`): `maxang * |ray|`.
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weigths)` (`cxx:727-731`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:736-739`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.curve_nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:743-746`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.curve_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:749-753`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        self.section_shape(nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)`
    /// (`cxx:757-771`).
    fn get_tolerance_3d(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut MathVector,
        tol1d: &mut MathVector,
    ) {
        self.fill_section_tolerance(bound_tol, surf_tol, angle_tol, tol3d, tol1d);
    }

    /// `Knots(TKnots)` (`cxx:774-777`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:782-785`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`cxx:845-1026`).
    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        self.section_with_derivatives(p, poles, dpoles, poles2d, dpoles2d, weights, dweights)
    }

    /// `Section(P, Poles, Poles2d, Weigths)` (`cxx:788-842`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:533-554`) forwards to
    /// `Blend_CSFunction::Section`, which answers `false`.
    #[allow(clippy::too_many_arguments)]
    fn section_with_derivatives2(
        &mut self,
        _p: &BlendPoint,
        _poles: &mut [GpPnt],
        _dpoles: &mut [GpVec],
        _d2poles: &mut [GpVec],
        _poles2d: &mut [GpPnt2d],
        _dpoles2d: &mut [GpVec2d],
        _d2poles2d: &mut [GpVec2d],
        _weights: &mut [f64],
        _dweights: &mut [f64],
        _d2weights: &mut [f64],
    ) -> bool {
        false
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1024-1029`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendCsFunction for CsCircular {
    /// `PointOnS()` (`cxx:402-405`).
    fn point_on_s(&self) -> &GpPnt {
        &self.pts
    }

    /// `PointOnC()` (`cxx:409-412`).
    fn point_on_c(&self) -> &GpPnt {
        &self.ptc
    }

    /// `Pnt2d()` (`cxx:416-419`).
    fn pnt2d(&self) -> &GpPnt2d {
        &self.pt2d
    }

    /// `ParameterOnC()` (`cxx:423-426`).
    fn parameter_on_c(&self) -> f64 {
        self.prmc
    }

    /// `IsTangencyPoint()` (`cxx:430-433`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS()` (`cxx:436-442`).
    fn tangent_on_s(&self) -> &GpVec {
        &self.tgs
    }

    /// `Tangent2d()` (`cxx:458-464`).
    fn tangent2d(&self) -> &GpVec2d {
        &self.tg2d
    }

    /// `TangentOnC()` (`cxx:446-452`).
    fn tangent_on_c(&self) -> &GpVec {
        &self.tgc
    }

    /// `Tangent(U, V, TgS, NmS)` (`cxx:468-487`).
    fn tangent(&self, u: f64, v: f64, tg_s: &mut GpVec, norm_s: &mut GpVec) {
        self.compute_tangent(u, v, tg_s, norm_s);
    }
}
