//! `BRepBlend_RstRstConstRad` — the constant-radius blending function between
//! two restrictions (pcurves on two surfaces), driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstConstRad.hxx:41-265`,
//! `.cxx:53-981` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(u, v)`: `u` on
//! the restriction `rst1` over `surf1`, `v` on `rst2` over `surf2`. With
//! `nplan` the unit tangent of the guide at the current parameter and
//! `theD = -(nplan . ptgui)`, the two equations are
//! `nplan . P1 + theD` and `nplan . P2 + theD`: both restriction points lie in
//! the section plane.
//!
//! Unlike [`super::surf_rst_const_rad`], the section is not anchored on a surface
//! normal: the circle through `P1` and `P2` of radius `ray` has its centre on the
//! perpendicular bisector of `P1P2`, inside the section plane
//! (`CenterCircleRst1Rst2`), and the section normals are the two radii
//! `Center -> P1` / `Center -> P2`.
//!
//! Translation notes:
//! - `Set(Radius, Choix)` (`cxx:485-490`) stores `ray = |Radius|` — the sign is
//!   *not* applied here, unlike the surface/restriction sibling. `choix` only
//!   drives the `choix % 2` / `choix > 2` switches inside the section and
//!   decroch code.
//! - `Set(First, Last)` (`cxx:160-163`) only narrows `tguide` through
//!   `Adaptor3d_Curve::Trim`; see [`super::surf_rst_const_rad`] for why the port
//!   aliases `guide`.
//! - `CenterCircleRst1Rst2` returns `false` when the two points are further
//!   apart than the circle diameter (`cxx:522-525`), with a `1.E-07` dead band
//!   on `ray^2 - |P1P2|^2/4`; the port returns `None` there and the two callers
//!   that ignore the flag (`Section`) fall back to the origin point, exactly
//!   like OCCT's untouched `gp_Pnt`.
//! - `myTConv` is left uninitialised by OCCT's constructor (`cxx:53-79` sets
//!   only `mySShape`); `GetShape` / `Set(TypeSection)` write it before any use,
//!   and the port starts from `Convert_TgtThetaOver2` (the enumerator `0`).
//! - `TangentOnRst1` / `Tangent2dOnRst1` / `TangentOnRst2` / `Tangent2dOnRst2`
//!   raise `Standard_DomainError` when `istangent` (`cxx:358-400`); the trait
//!   returns a reference, and every caller in OCCT guards the call with
//!   `!IsTangencyPoint()`, so the stored value is returned.

pub mod constants;
pub mod equations;
pub mod section;
pub mod solution;

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{
    BlendAppFunction, BlendDecrochStatus, BlendFuncSectionShape, BlendPoint, BlendRstRstFunction,
};

pub use constants::{
    DECROCH_DOT_EPS, DIST_EPS, EPS, GAUSS_DEFAULT_EPS, GAUSS_SECTION_EPS, SVD_EPS,
};

/// `class BRepBlend_RstRstConstRad : public Blend_RstRstFunction`
/// (`hxx:41-265`).
pub struct SurfRstRstConstRad {
    /// `surf1` (`hxx:227`): the surface `rst1` lives on.
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:228`): the surface `rst2` lives on.
    surf2: Arc<dyn Surface>,
    /// `rst1` (`hxx:229`): the first restriction pcurve.
    rst1: Arc<dyn Curve2d>,
    /// `rst2` (`hxx:230`): the second restriction pcurve.
    rst2: Arc<dyn Curve2d>,
    /// `cons1` (`hxx:231`): `rst1` seen as a 3D curve on `surf1`.
    cons1: CurveOnSurface,
    /// `cons2` (`hxx:232`): `rst2` seen as a 3D curve on `surf2`.
    cons2: CurveOnSurface,
    /// `guide` (`hxx:233`): the guide line.
    guide: Arc<dyn Curve>,
    /// `tguide` (`hxx:234`): `guide->Trim(First, Last, 1e-12)`; see the module
    /// note -- evaluation is identical, so it aliases `guide`.
    tguide: Arc<dyn Curve>,
    /// `ptrst1` (`hxx:235`).
    ptrst1: GpPnt,
    /// `ptrst2` (`hxx:236`).
    ptrst2: GpPnt,
    /// `pt2drst1` (`hxx:237`).
    pt2drst1: GpPnt2d,
    /// `pt2drst2` (`hxx:238`).
    pt2drst2: GpPnt2d,
    /// `prmrst1` (`hxx:239`).
    prmrst1: f64,
    /// `prmrst2` (`hxx:240`).
    prmrst2: f64,
    /// `istangent` (`hxx:241`).
    istangent: bool,
    /// `tgrst1` (`hxx:242`).
    tgrst1: GpVec,
    /// `tg2drst1` (`hxx:243`).
    tg2drst1: GpVec2d,
    /// `tgrst2` (`hxx:244`).
    tgrst2: GpVec,
    /// `tg2drst2` (`hxx:245`).
    tg2drst2: GpVec2d,
    /// `ray` (`hxx:246`): `|Radius|`, set by `Set(Radius, Choix)`.
    ray: f64,
    /// `choix` (`hxx:247`).
    choix: i32,
    /// `ptgui` (`hxx:248`).
    ptgui: GpPnt,
    /// `d1gui` (`hxx:249`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:250`).
    d2gui: GpVec,
    /// `nplan` (`hxx:251`).
    nplan: GpVec,
    /// `normtg` (`hxx:252`).
    normtg: f64,
    /// `theD` (`hxx:253`).
    the_d: f64,
    /// `surfref1` (`hxx:255`), set by `Set(SurfRef1, RstRef1, ..)`.
    surfref1: Option<Arc<dyn Surface>>,
    /// `rstref1` (`hxx:256`), set by `Set(SurfRef1, RstRef1, ..)`.
    rstref1: Option<Arc<dyn Curve2d>>,
    /// `surfref2` (`hxx:257`), set by `Set(.., SurfRef2, RstRef2)`.
    surfref2: Option<Arc<dyn Surface>>,
    /// `rstref2` (`hxx:258`), set by `Set(.., SurfRef2, RstRef2)`.
    rstref2: Option<Arc<dyn Curve2d>>,
    /// `maxang` (`hxx:258`).
    maxang: f64,
    /// `minang` (`hxx:259`).
    minang: f64,
    /// `distmin` (`hxx:260`).
    distmin: f64,
    /// `mySShape` (`hxx:261`).
    my_sshape: BlendFuncSectionShape,
    /// `myTConv` (`hxx:262`); see the module note.
    my_tconv: ParameterisationType,
}

impl SurfRstRstConstRad {
    /// `BRepBlend_RstRstConstRad(Surf1, Rst1, Surf2, Rst2, CGuide)`
    /// (`cxx:53-79`).
    pub fn new(
        surf1: Arc<dyn Surface>,
        rst1: Arc<dyn Curve2d>,
        surf2: Arc<dyn Surface>,
        rst2: Arc<dyn Curve2d>,
        cguide: Arc<dyn Curve>,
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
            ray: 0.0,
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
            // `RealFirst()` / `RealLast()` (`hxx:258-260`).
            maxang: f64::MIN,
            minang: f64::MAX,
            distmin: f64::MAX,
            my_sshape: BlendFuncSectionShape::Rational,
            // OCCT leaves `myTConv` uninitialised; see the module note.
            my_tconv: ParameterisationType::TgtThetaOver2,
        }
    }

    /// `Set(SurfRef1, RstRef1, SurfRef2, RstRef2)` (`cxx:138-148`): the faces
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

    /// `Set(Radius, Choix)` (`cxx:485-490`): `ray = |Radius|`, `choix` as is.
    pub fn set_radius(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        self.ray = radius.abs();
    }

    /// `Set(TypeSection)` (`cxx:494-497`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `ray` (`hxx:246`): the section radius, always non-negative.
    pub fn ray(&self) -> f64 {
        self.ray
    }

    /// `choix` (`hxx:247`).
    pub fn choix(&self) -> i32 {
        self.choix
    }

    /// `maxang` (`hxx:258`), the widest section angle seen by `IsSolution`.
    pub fn max_ang(&self) -> f64 {
        self.maxang
    }

    /// `minang` (`hxx:259`), the narrowest section angle seen by `IsSolution`.
    pub fn min_ang(&self) -> f64 {
        self.minang
    }

    /// `surf1` (`hxx:227`).
    pub fn surf1(&self) -> &Arc<dyn Surface> {
        &self.surf1
    }

    /// `surf2` (`hxx:228`).
    pub fn surf2(&self) -> &Arc<dyn Surface> {
        &self.surf2
    }
}

impl MathFunctionSetWithDerivatives for SurfRstRstConstRad {
    /// `NbVariables()` (`cxx:83-86`): `2`.
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`cxx:90-93`): `2`.
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`cxx:95-107`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:111-127`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d)
    }

    /// `Values(X, F, D)` (`cxx:131-135`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, f, d)
    }
}

impl BlendAppFunction for SurfRstRstConstRad {
    /// `Set(Param)` (`cxx:152-157`): the section plane.
    fn set(&mut self, param: f64) {
        let (ptgui, d1gui, d2gui) = self.tguide.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        let the_d = self.nplan.dot(&equations::pnt_vec(&self.ptgui));
        self.the_d = -the_d;
    }

    /// `Set(First, Last)` (`cxx:160-163`): `tguide = guide->Trim(First, Last,
    /// 1.e-12)`. See the module note: `Trim` does not change the evaluated
    /// values, so `tguide` keeps aliasing `guide`.
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// `GetTolerance(Tolerance, Tol)` (`cxx:167-171`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:175-183`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:187-297`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:305-308`).
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

    /// `IsRational()` (`cxx:590-593`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:597-600`).
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:604-608`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:612-615`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.guide_nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:619-622`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.guide_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:626-630`).
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
    /// (`cxx:637-650`).
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

    /// `Knots(TKnots)` (`cxx:654-657`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:661-664`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`cxx:732-947`).
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

    /// `Section(P, Poles, Poles2d, Weigths)` (`cxx:666-729`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:951-965`): `false`.
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

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:967-981`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendRstRstFunction for SurfRstRstConstRad {
    /// `PointOnRst1()` (`cxx:310-313`).
    fn point_on_rst1(&self) -> &GpPnt {
        &self.ptrst1
    }

    /// `PointOnRst2()` (`cxx:317-320`).
    fn point_on_rst2(&self) -> &GpPnt {
        &self.ptrst2
    }

    /// `Pnt2dOnRst1()` (`cxx:324-327`).
    fn pnt2d_on_rst1(&self) -> &GpPnt2d {
        &self.pt2drst1
    }

    /// `Pnt2dOnRst2()` (`cxx:331-334`).
    fn pnt2d_on_rst2(&self) -> &GpPnt2d {
        &self.pt2drst2
    }

    /// `ParameterOnRst1()` (`cxx:338-341`).
    fn parameter_on_rst1(&self) -> f64 {
        self.prmrst1
    }

    /// `ParameterOnRst2()` (`cxx:344-347`).
    fn parameter_on_rst2(&self) -> f64 {
        self.prmrst2
    }

    /// `IsTangencyPoint()` (`cxx:352-355`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnRst1()` (`cxx:358-365`); OCCT raises `Standard_DomainError`
    /// when `istangent` (see the module note).
    fn tangent_on_rst1(&self) -> &GpVec {
        &self.tgrst1
    }

    /// `Tangent2dOnRst1()` (`cxx:369-376`).
    fn tangent2d_on_rst1(&self) -> &GpVec2d {
        &self.tg2drst1
    }

    /// `TangentOnRst2()` (`cxx:381-388`).
    fn tangent_on_rst2(&self) -> &GpVec {
        &self.tgrst2
    }

    /// `Tangent2dOnRst2()` (`cxx:392-399`).
    fn tangent2d_on_rst2(&self) -> &GpVec2d {
        &self.tg2drst2
    }

    /// `Decroch(Sol, NRst1, TgRst1, NRst2, TgRst2)` (`cxx:405-483`).
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
