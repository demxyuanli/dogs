//! `BRepBlend_SurfRstConstRad` — the constant-radius blending function between a
//! surface and a pcurve on another surface, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstConstRad.hxx:41-244`,
//! `.cxx:51-1045` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(u, v, w)`: `u, v`
//! on `surf`, `w` on the restriction `rst` (a pcurve on `surfrst`). The three
//! equations, with `nplan` the unit tangent of the guide line at the current
//! parameter and `theD = -(nplan . ptgui)`, are
//! 1. `nplan . pts + theD` — the surface point lies in the section plane,
//! 2. `nplan . ptrst + theD` — the restriction point lies in the same plane,
//! 3. `(|ray * ns + (pts - ptrst)| / ray)^2 - 1` — the two points are `|ray|`
//!    apart, `ns` being the surface normal projected into the section plane.
//!
//! Translation notes:
//! - `Set(First, Last)` (`cxx:270-273`) only narrows `tguide` through
//!   `Adaptor3d_Curve::Trim`. `GeomAdaptor_Curve::Trim` keeps the basis
//!   parameterisation (`myFirst` / `myLast` only bound the interval), and this
//!   class never asks `tguide` for its bounds, so the port evaluates the guide
//!   itself and keeps the interval as documentation. See
//!   [`occt_geom::curve_on_surface`] for the same reasoning on `cons`.
//! - `myTConv` is left uninitialised by OCCT's constructor (`cxx:51-73` sets
//!   only `mySShape`); `GetShape` / `Set(TypeSection)` write it before any use,
//!   and the port starts from `Convert_TgtThetaOver2` (the enumerator `0`).
//! - `TangentOnS` / `Tangent2dOnS` / `TangentOnRst` / `Tangent2dOnRst` raise
//!   `Standard_DomainError` when `istangent` (`cxx:478-518`); the trait returns
//!   a reference, and every caller in OCCT guards the call with
//!   `!IsTangencyPoint()`, so the stored value is returned.

pub mod constants;
pub mod equations;
pub mod section;
pub mod solution;

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::convert::ParameterisationType;
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom2d::curve::Curve2d;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{
    BlendAppFunction, BlendFuncSectionShape, BlendPoint, BlendSurfRstFunction,
};

pub use constants::{
    DECROCH_DOT_EPS, DECROCH_NT_EPS, EPS, GAUSS_DEFAULT_EPS, GAUSS_SECTION_EPS, SVD_EPS,
};

/// `class BRepBlend_SurfRstConstRad : public Blend_SurfRstFunction`
/// (`hxx:41-244`).
pub struct SurfRstConstRad {
    /// `surf` (`hxx:212`): the surface the section starts from.
    surf: Arc<dyn Surface>,
    /// `surfrst` (`hxx:213`): the surface carrying the restriction.
    surfrst: Arc<dyn Surface>,
    /// `rst` (`hxx:214`): the restriction pcurve on `surfrst`.
    rst: Arc<dyn Curve2d>,
    /// `cons` (`hxx:215`): `rst` seen as a 3D curve on `surfrst`.
    cons: CurveOnSurface,
    /// `guide` (`hxx:216`): the guide line.
    guide: Arc<dyn Curve>,
    /// `tguide` (`hxx:217`): `guide->Trim(First, Last, 1e-12)`; see the module
    /// note -- evaluation is identical, so it aliases `guide`.
    tguide: Arc<dyn Curve>,
    /// `pts` (`hxx:218`).
    pts: GpPnt,
    /// `ptrst` (`hxx:219`).
    ptrst: GpPnt,
    /// `pt2ds` (`hxx:220`).
    pt2ds: GpPnt2d,
    /// `pt2drst` (`hxx:221`).
    pt2drst: GpPnt2d,
    /// `prmrst` (`hxx:222`).
    prmrst: f64,
    /// `istangent` (`hxx:223`).
    istangent: bool,
    /// `tgs` (`hxx:224`).
    tgs: GpVec,
    /// `tg2ds` (`hxx:225`).
    tg2ds: GpVec2d,
    /// `tgrst` (`hxx:226`).
    tgrst: GpVec,
    /// `tg2drst` (`hxx:227`).
    tg2drst: GpVec2d,
    /// `ray` (`hxx:228`).
    ray: f64,
    /// `choix` (`hxx:229`).
    choix: i32,
    /// `ptgui` (`hxx:230`).
    ptgui: GpPnt,
    /// `d1gui` (`hxx:231`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:232`).
    d2gui: GpVec,
    /// `nplan` (`hxx:233`).
    nplan: GpVec,
    /// `normtg` (`hxx:234`).
    normtg: f64,
    /// `theD` (`hxx:235`).
    the_d: f64,
    /// `surfref` (`hxx:236`), set by `Set(SurfRef, RstRef)`.
    surfref: Option<Arc<dyn Surface>>,
    /// `rstref` (`hxx:237`), set by `Set(SurfRef, RstRef)`.
    rstref: Option<Arc<dyn Curve2d>>,
    /// `maxang` (`hxx:238`).
    maxang: f64,
    /// `minang` (`hxx:239`).
    minang: f64,
    /// `distmin` (`hxx:240`).
    distmin: f64,
    /// `mySShape` (`hxx:241`).
    my_sshape: BlendFuncSectionShape,
    /// `myTConv` (`hxx:242`); see the module note.
    my_tconv: ParameterisationType,
}

impl SurfRstConstRad {
    /// `BRepBlend_SurfRstConstRad(Surf, SurfRst, Rst, CGuide)` (`cxx:51-73`).
    pub fn new(
        surf: Arc<dyn Surface>,
        surf_rst: Arc<dyn Surface>,
        rst: Arc<dyn Curve2d>,
        cguide: Arc<dyn Curve>,
    ) -> Self {
        Self {
            cons: CurveOnSurface::new(rst.clone(), surf_rst.clone()),
            surf,
            surfrst: surf_rst,
            rst,
            tguide: cguide.clone(),
            guide: cguide,
            pts: GpPnt::zero(),
            ptrst: GpPnt::zero(),
            pt2ds: GpPnt2d::zero(),
            pt2drst: GpPnt2d::zero(),
            prmrst: 0.0,
            istangent: true,
            tgs: GpVec::zero(),
            tg2ds: GpVec2d::zero(),
            tgrst: GpVec::zero(),
            tg2drst: GpVec2d::zero(),
            ray: 0.0,
            choix: 0,
            ptgui: GpPnt::zero(),
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            the_d: 0.0,
            surfref: None,
            rstref: None,
            // `RealFirst()` / `RealLast()` (`hxx:238-240`).
            maxang: f64::MIN,
            minang: f64::MAX,
            distmin: f64::MAX,
            my_sshape: BlendFuncSectionShape::Rational,
            // OCCT leaves `myTConv` uninitialised; see the module note.
            my_tconv: ParameterisationType::TgtThetaOver2,
        }
    }

    /// `Set(SurfRef, RstRef)` (`cxx:246-251`): the face and the pcurve
    /// `Decroch` measures the restriction normal on.
    pub fn set_refs(&mut self, surf_ref: Arc<dyn Surface>, rst_ref: Arc<dyn Curve2d>) {
        self.surfref = Some(surf_ref);
        self.rstref = Some(rst_ref);
    }

    /// `Set(Radius, Choix)` (`cxx:579-596`): the sign of `ray` follows the
    /// fillet quadrant. The `default:` arm repeats the `1`/`2` arm.
    pub fn set_radius(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => self.ray = -radius.abs(),
            3 | 4 => self.ray = radius.abs(),
            _ => self.ray = -radius.abs(),
        }
    }

    /// `Set(TypeSection)` (`cxx:600-603`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `ray` (`hxx:228`): the signed section radius.
    pub fn ray(&self) -> f64 {
        self.ray
    }

    /// `choix` (`hxx:229`).
    pub fn choix(&self) -> i32 {
        self.choix
    }

    /// `maxang` (`hxx:238`), the widest section angle seen by `IsSolution`.
    pub fn max_ang(&self) -> f64 {
        self.maxang
    }

    /// `minang` (`hxx:239`), the narrowest section angle seen by `IsSolution`.
    pub fn min_ang(&self) -> f64 {
        self.minang
    }

    /// `surf` (`hxx:212`).
    pub fn surf(&self) -> &Arc<dyn Surface> {
        &self.surf
    }

    /// `surfrst` (`hxx:213`).
    pub fn surf_rst(&self) -> &Arc<dyn Surface> {
        &self.surfrst
    }
}

impl MathFunctionSetWithDerivatives for SurfRstConstRad {
    /// `NbVariables()` (`cxx:76-79`): `3`.
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:82-86`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:90-108`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:113-170`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d)
    }

    /// `Values(X, F, D)` (`cxx:175-241`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, f, d)
    }
}

impl BlendAppFunction for SurfRstConstRad {
    /// `Set(Param)` (`cxx:255-267`): the section plane.
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
    }

    /// `Set(First, Last)` (`cxx:270-273`): `tguide = guide->Trim(First, Last,
    /// 1.e-12)`. See the module note: `Trim` does not change the evaluated
    /// values, so `tguide` keeps aliasing `guide`.
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// `GetTolerance(Tolerance, Tol)` (`cxx:277-282`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:286-308`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:311-427`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:430-433`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.distmin)
    }

    /// `Pnt1()`: `PointOnS()` (`Blend_SurfRstFunction.cxx:22-25`).
    fn pnt1(&self) -> &GpPnt {
        &self.pts
    }

    /// `Pnt2()`: `PointOnRst()` (`Blend_SurfRstFunction.cxx:27-30`).
    fn pnt2(&self) -> &GpPnt {
        &self.ptrst
    }

    /// `IsRational()` (`cxx:661-664`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:668-671`).
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:675-680`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:684-687`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.guide_nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:691-695`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.guide_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:697-702`).
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
    /// (`cxx:710-723`).
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

    /// `Knots(TKnots)` (`cxx:727-730`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:734-737`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`cxx:802-1014`).
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
        self.section_with_derivatives(
            p, poles, dpoles, poles2d, dpoles2d, weights, dweights,
        )
    }

    /// `Section(P, Poles, Poles2d, Weigths)` (`cxx:742-798`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:1017-1031`): `false`.
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

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1033-1045`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendSurfRstFunction for SurfRstConstRad {
    /// `PointOnS()` (`cxx:437-441`).
    fn point_on_s(&self) -> &GpPnt {
        &self.pts
    }

    /// `PointOnRst()` (`cxx:444-448`).
    fn point_on_rst(&self) -> &GpPnt {
        &self.ptrst
    }

    /// `Pnt2dOnS()` (`cxx:450-455`).
    fn pnt2d_on_s(&self) -> &GpPnt2d {
        &self.pt2ds
    }

    /// `Pnt2dOnRst()` (`cxx:457-462`).
    fn pnt2d_on_rst(&self) -> &GpPnt2d {
        &self.pt2drst
    }

    /// `ParameterOnRst()` (`cxx:464-468`).
    fn parameter_on_rst(&self) -> f64 {
        self.prmrst
    }

    /// `IsTangencyPoint()` (`cxx:470-474`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS()` (`cxx:476-484`); OCCT raises `Standard_DomainError` when
    /// `istangent` (see the module note).
    fn tangent_on_s(&self) -> &GpVec {
        &self.tgs
    }

    /// `Tangent2dOnS()` (`cxx:486-496`).
    fn tangent2d_on_s(&self) -> &GpVec2d {
        &self.tg2ds
    }

    /// `TangentOnRst()` (`cxx:498-508`).
    fn tangent_on_rst(&self) -> &GpVec {
        &self.tgrst
    }

    /// `Tangent2dOnRst()` (`cxx:510-520`).
    fn tangent2d_on_rst(&self) -> &GpVec2d {
        &self.tg2drst
    }

    /// `Decroch(Sol, NS, TgS)` (`cxx:523-575`).
    fn decroch(&self, sol: &MathVector, ns: &mut GpVec, tg_s: &mut GpVec) -> bool {
        self.compute_decroch(sol, ns, tg_s)
    }
}
