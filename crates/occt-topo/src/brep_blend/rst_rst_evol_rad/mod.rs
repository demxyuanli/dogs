//! `BRepBlend_RstRstEvolRad` — the evolutive-radius blending function between
//! two restrictions (pcurves on two surfaces), driven by a guide line and a
//! radius law.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstEvolRad.hxx:41-268`,
//! `.cxx:117-1092` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(u, v)`: `u` on
//! the restriction `rst1` over `surf1`, `v` on `rst2` over `surf2`. With
//! `nplan` the unit tangent of the guide at the current parameter and
//! `theD = -(nplan . ptgui)`, the two equations are
//! `nplan . P1 + theD` and `nplan . P2 + theD`: both restriction points lie in
//! the section plane. The circle through `P1` and `P2` of radius `|ray|` has its
//! centre on the perpendicular bisector of `P1P2`, inside the section plane
//! (`CenterCircleRst1Rst2`), and the section normals are the two radii
//! `Center -> P1` / `Center -> P2`.
//!
//! Differences from [`super::rst_rst_const_rad`]
//! (`BRepBlend_RstRstConstRad`):
//! - `ray` / `dray` come from `tevol->Value` / `tevol->D1` and keep the sign the
//!   law gives; there is no `Set(Radius, Choix)`, only `Set(Choix)`
//!   (`cxx:554-557`) which stores the quadrant index.
//! - `Section(P, ..., DPoles, ...)` normalises the two section radii *after* the
//!   derivative block (`cxx:921-974`), so the unnormalised vectors enter
//!   `d1n1` / `d1n2`, and passes `dray` through to `GeomFill::GetCircle`
//!   (`cxx:1033`).
//! - `NbIntervals` / `Intervals` merge the guide's intervals with the law's
//!   (`FusionneIntervalles`, `cxx:54-113`).
//!
//! Translation notes:
//! - `Set(First, Last)` (`cxx:228-232`) narrows both `tguide` and `tevol` with
//!   `Trim(First, Last, 1.e-12)`. Every `Law_*` implementation in OCCT keeps the
//!   parameterisation on the kept interval, so `tguide` aliases `guide` as in
//!   [`super::surf_rst_const_rad`], while `tevol` is replaced by
//!   `fevol.trim(..)` so the port still narrows the law the way OCCT does.
//! - `myTConv` is left uninitialised by OCCT's constructor (`cxx:117-139` sets
//!   only `mySShape`); `GetShape` / `Set(TypeSection)` write it before any use,
//!   and the port starts from `Convert_TgtThetaOver2` (the enumerator `0`).
//! - `ray`, `dray` and `choix` are uninitialised in OCCT until the first
//!   evaluation / `Set(Choix)`; the port starts them at `0` and documents it
//!   rather than inventing a quadrant.
//! - `TangentOnRst1` / `Tangent2dOnRst1` / `TangentOnRst2` / `Tangent2dOnRst2`
//!   raise `Standard_DomainError` when `istangent` (`cxx:428-470`); the trait
//!   returns a reference, and every caller in OCCT guards the call with
//!   `!IsTangencyPoint()`, so the stored value is returned.

pub mod constants;
pub mod equations;
pub mod queries;
pub mod section;
pub mod solution;

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::law::LawFunction;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{
    BlendAppFunction, BlendDecrochStatus, BlendFuncSectionShape, BlendPoint, BlendRstRstFunction,
};

pub use constants::{
    DECROCH_DOT_EPS, DIST_EPS, EPS, FUSIONNE_EPS_FACTOR, GAUSS_DEFAULT_EPS, GAUSS_SECTION_EPS,
    SVD_EPS,
};

/// `class BRepBlend_RstRstEvolRad : public Blend_RstRstFunction`
/// (`hxx:41-268`).
pub struct SurfRstRstEvolRad {
    /// `surf1` (`hxx:229`): the surface `rst1` lives on.
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:230`): the surface `rst2` lives on.
    surf2: Arc<dyn Surface>,
    /// `rst1` (`hxx:231`): the first restriction pcurve.
    rst1: Arc<dyn Curve2d>,
    /// `rst2` (`hxx:232`): the second restriction pcurve.
    rst2: Arc<dyn Curve2d>,
    /// `cons1` (`hxx:233`): `rst1` seen as a 3D curve on `surf1`.
    cons1: CurveOnSurface,
    /// `cons2` (`hxx:234`): `rst2` seen as a 3D curve on `surf2`.
    cons2: CurveOnSurface,
    /// `guide` (`hxx:235`): the guide line.
    guide: Arc<dyn Curve>,
    /// `tguide` (`hxx:236`): `guide->Trim(First, Last, 1e-12)`; see the module
    /// note -- evaluation is identical, so it aliases `guide`.
    tguide: Arc<dyn Curve>,
    /// `ptrst1` (`hxx:237`).
    ptrst1: GpPnt,
    /// `ptrst2` (`hxx:238`).
    ptrst2: GpPnt,
    /// `pt2drst1` (`hxx:239`).
    pt2drst1: GpPnt2d,
    /// `pt2drst2` (`hxx:240`).
    pt2drst2: GpPnt2d,
    /// `prmrst1` (`hxx:241`).
    prmrst1: f64,
    /// `prmrst2` (`hxx:242`).
    prmrst2: f64,
    /// `istangent` (`hxx:243`).
    istangent: bool,
    /// `tgrst1` (`hxx:244`).
    tgrst1: GpVec,
    /// `tg2drst1` (`hxx:245`).
    tg2drst1: GpVec2d,
    /// `tgrst2` (`hxx:246`).
    tgrst2: GpVec,
    /// `tg2drst2` (`hxx:247`).
    tg2drst2: GpVec2d,
    /// `ray` (`hxx:248`): the radius the law gave at the current parameter.
    ray: f64,
    /// `dray` (`hxx:249`): the derivative of `ray` with respect to the guide
    /// parameter.
    dray: f64,
    /// `choix` (`hxx:250`).
    choix: i32,
    /// `ptgui` (`hxx:251`).
    ptgui: GpPnt,
    /// `d1gui` (`hxx:252`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:253`).
    d2gui: GpVec,
    /// `nplan` (`hxx:254`).
    nplan: GpVec,
    /// `normtg` (`hxx:255`).
    normtg: f64,
    /// `theD` (`hxx:256`).
    the_d: f64,
    /// `surfref1` (`hxx:257`), set by `Set(SurfRef1, RstRef1, ..)`.
    surfref1: Option<Arc<dyn Surface>>,
    /// `rstref1` (`hxx:258`), set by `Set(SurfRef1, RstRef1, ..)`.
    rstref1: Option<Arc<dyn Curve2d>>,
    /// `surfref2` (`hxx:259`), set by `Set(.., SurfRef2, RstRef2)`.
    surfref2: Option<Arc<dyn Surface>>,
    /// `rstref2` (`hxx:260`), set by `Set(.., SurfRef2, RstRef2)`.
    rstref2: Option<Arc<dyn Curve2d>>,
    /// `maxang` (`hxx:261`).
    maxang: f64,
    /// `minang` (`hxx:262`).
    minang: f64,
    /// `distmin` (`hxx:263`).
    distmin: f64,
    /// `mySShape` (`hxx:264`).
    my_sshape: BlendFuncSectionShape,
    /// `myTConv` (`hxx:265`); see the module note.
    my_tconv: ParameterisationType,
    /// `tevol` (`hxx:266`): the radius law, trimmed to the guide interval.
    tevol: Arc<dyn LawFunction>,
    /// `fevol` (`hxx:267`): the untrimmed radius law `Set(First, Last)` trims.
    fevol: Arc<dyn LawFunction>,
}

impl SurfRstRstEvolRad {
    /// `BRepBlend_RstRstEvolRad(Surf1, Rst1, Surf2, Rst2, CGuide, Evol)`
    /// (`cxx:117-139`).
    pub fn new(
        surf1: Arc<dyn Surface>,
        rst1: Arc<dyn Curve2d>,
        surf2: Arc<dyn Surface>,
        rst2: Arc<dyn Curve2d>,
        cguide: Arc<dyn Curve>,
        evol: Arc<dyn LawFunction>,
    ) -> Self {
        Self {
            cons1: CurveOnSurface::new(rst1.clone(), surf1.clone()),
            cons2: CurveOnSurface::new(rst2.clone(), surf2.clone()),
            surf1,
            surf2,
            rst1,
            rst2,
            tguide: cguide.clone(),
            guide: cguide,
            ptrst1: GpPnt::zero(),
            ptrst2: GpPnt::zero(),
            pt2drst1: GpPnt2d::zero(),
            pt2drst2: GpPnt2d::zero(),
            prmrst1: 0.0,
            prmrst2: 0.0,
            istangent: true,
            tgrst1: GpVec::zero(),
            tg2drst1: GpVec2d::zero(),
            tgrst2: GpVec::zero(),
            tg2drst2: GpVec2d::zero(),
            // OCCT leaves `ray` / `dray` uninitialised until the first
            // `Set` / `Section`; see the module note.
            ray: 0.0,
            dray: 0.0,
            choix: 0,
            ptgui: GpPnt::zero(),
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            the_d: 0.0,
            surfref1: None,
            rstref1: None,
            surfref2: None,
            rstref2: None,
            // `RealFirst()` / `RealLast()` (`hxx:261-263`).
            maxang: f64::MIN,
            minang: f64::MAX,
            distmin: f64::MAX,
            my_sshape: BlendFuncSectionShape::Rational,
            // OCCT leaves `myTConv` uninitialised; see the module note.
            my_tconv: ParameterisationType::TgtThetaOver2,
            tevol: evol.clone(),
            fevol: evol,
        }
    }

    /// `Set(SurfRef1, RstRef1, SurfRef2, RstRef2)` (`cxx:198-208`): the faces
    /// and pcurves `Decroch` measures the restriction normals on.
    pub fn set_refs(
        &mut self,
        surf_ref1: Arc<dyn Surface>,
        rst_ref1: Arc<dyn Curve2d>,
        surf_ref2: Arc<dyn Surface>,
        rst_ref2: Arc<dyn Curve2d>,
    ) {
        self.surfref1 = Some(surf_ref1);
        self.surfref2 = Some(surf_ref2);
        self.rstref1 = Some(rst_ref1);
        self.rstref2 = Some(rst_ref2);
    }

    /// `Set(Choix)` (`cxx:554-557`): only the quadrant index is stored; the
    /// radius keeps the sign the law gives it.
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
    }

    /// `Set(TypeSection)` (`cxx:561-564`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `ray` (`hxx:248`): the section radius the law gave last.
    pub fn ray(&self) -> f64 {
        self.ray
    }

    /// `dray` (`hxx:249`): the radius derivative the law gave last.
    pub fn dray(&self) -> f64 {
        self.dray
    }

    /// `choix` (`hxx:250`).
    pub fn choix(&self) -> i32 {
        self.choix
    }

    /// `maxang` (`hxx:261`), the widest section angle seen by `IsSolution`.
    pub fn max_ang(&self) -> f64 {
        self.maxang
    }

    /// `minang` (`hxx:262`), the narrowest section angle seen by `IsSolution`.
    pub fn min_ang(&self) -> f64 {
        self.minang
    }

    /// `surf1` (`hxx:229`).
    pub fn surf1(&self) -> &Arc<dyn Surface> {
        &self.surf1
    }

    /// `surf2` (`hxx:230`).
    pub fn surf2(&self) -> &Arc<dyn Surface> {
        &self.surf2
    }

    /// `tevol` (`hxx:266`), the radius law currently in force.
    pub fn law(&self) -> &Arc<dyn LawFunction> {
        &self.tevol
    }
}

impl MathFunctionSetWithDerivatives for SurfRstRstEvolRad {
    /// `NbVariables()` (`cxx:143-146`): `2`.
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`cxx:150-153`): `2`.
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`cxx:157-167`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:170-186`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d)
    }

    /// `Values(X, F, D)` (`cxx:188-195`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, f, d)
    }
}

impl BlendAppFunction for SurfRstRstEvolRad {
    /// `Set(Param)` (`cxx:211-224`): the section plane and the radius law.
    fn set(&mut self, param: f64) {
        self.d1gui = GpVec::zero();
        self.nplan = GpVec::zero();
        let (ptgui, d1gui, d2gui) = self.tguide.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        let the_d = self.nplan.dot(&equations::pnt_vec(&self.ptgui));
        self.the_d = -the_d;
        let (ray, dray) = self.tevol.d1(param);
        self.ray = ray;
        self.dray = dray;
    }

    /// `Set(First, Last)` (`cxx:228-232`): `tguide = guide->Trim(First, Last,
    /// 1.e-12)` (see the module note) and `tevol = fevol->Trim(First, Last,
    /// 1.e-12)`.
    fn set_interval(&mut self, first: f64, last: f64) {
        match self.fevol.trim(first, last, 1.0e-12) {
            Ok(trimmed) => self.tevol = trimmed,
            Err(_) => {
                // OCCT's `Law_Function::Trim` has no failure arm; the port keeps
                // the untrimmed law rather than inventing one.
            }
        }
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:236-240`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:244-251`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:254-368`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:372-375`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.distmin)
    }

    /// `Pnt1()`: `PointOnRst1()` (`Blend_RstRstFunction.cxx:22-25`).
    fn pnt1(&self) -> &GpPnt {
        &self.ptrst1
    }

    /// `Pnt2()`: `PointOnRst2()` (`Blend_RstRstFunction.cxx:27-30`).
    fn pnt2(&self) -> &GpPnt {
        &self.ptrst2
    }

    /// `IsRational()` (`cxx:658-661`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:665-668`).
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:672-676`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:680-698`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.nb_intervals_impl(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:703-726`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.intervals_impl(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:731-735`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        self.section_shape(nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1D)`
    /// (`cxx:742-756`).
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

    /// `Knots(TKnots)` (`cxx:759-762`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:766-769`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`cxx:839-1059`).
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

    /// `Section(P, Poles, Poles2d, Weigths)` (`cxx:773-836`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:1063-1074`): `false`.
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
        self.no_second_derivative_section()
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1077-1092`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendRstRstFunction for SurfRstRstEvolRad {
    /// `PointOnRst1()` (`cxx:379-382`).
    fn point_on_rst1(&self) -> &GpPnt {
        &self.ptrst1
    }

    /// `PointOnRst2()` (`cxx:386-389`).
    fn point_on_rst2(&self) -> &GpPnt {
        &self.ptrst2
    }

    /// `Pnt2dOnRst1()` (`cxx:393-396`).
    fn pnt2d_on_rst1(&self) -> &GpPnt2d {
        &self.pt2drst1
    }

    /// `Pnt2dOnRst2()` (`cxx:400-403`).
    fn pnt2d_on_rst2(&self) -> &GpPnt2d {
        &self.pt2drst2
    }

    /// `ParameterOnRst1()` (`cxx:407-410`).
    fn parameter_on_rst1(&self) -> f64 {
        self.prmrst1
    }

    /// `ParameterOnRst2()` (`cxx:414-417`).
    fn parameter_on_rst2(&self) -> f64 {
        self.prmrst2
    }

    /// `IsTangencyPoint()` (`cxx:421-424`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnRst1()` (`cxx:428-435`); OCCT raises `Standard_DomainError`
    /// when `istangent` (see the module note).
    fn tangent_on_rst1(&self) -> &GpVec {
        &self.tgrst1
    }

    /// `Tangent2dOnRst1()` (`cxx:439-446`).
    fn tangent2d_on_rst1(&self) -> &GpVec2d {
        &self.tg2drst1
    }

    /// `TangentOnRst2()` (`cxx:450-457`).
    fn tangent_on_rst2(&self) -> &GpVec {
        &self.tgrst2
    }

    /// `Tangent2dOnRst2()` (`cxx:461-468`).
    fn tangent2d_on_rst2(&self) -> &GpVec2d {
        &self.tg2drst2
    }

    /// `Decroch(Sol, NRst1, TgRst1, NRst2, TgRst2)` (`cxx:472-552`).
    fn decroch(
        &self,
        sol: &MathVector,
        n_rst1: &mut GpVec,
        tg_rst1: &mut GpVec,
        n_rst2: &mut GpVec,
        tg_rst2: &mut GpVec,
    ) -> BlendDecrochStatus {
        self.compute_decroch(sol, n_rst1, tg_rst1, n_rst2, tg_rst2)
    }
}
