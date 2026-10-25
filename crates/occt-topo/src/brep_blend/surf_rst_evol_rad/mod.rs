//! `BRepBlend_SurfRstEvolRad` — the evolutive-radius blending function between
//! a surface and a pcurve on another surface, driven by a guide line and a
//! radius law.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstEvolRad.hxx:41-248`,
//! `.cxx:113-1166` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(u, v, w)`: `u, v`
//! on `surf`, `w` on the restriction `rst` (a pcurve on `surfrst`). The three
//! equations, with `nplan` the unit tangent of the guide at the current
//! parameter, `theD = -(nplan . ptgui)` and `ray` the signed radius the law
//! `tevol` gives at that parameter, are
//! 1. `nplan . pts + theD` — the surface point lies in the section plane,
//! 2. `nplan . ptrst + theD` — the restriction point lies in the same plane,
//! 3. `|ray * ns + (pts - ptrst)|^2 - ray^2` — the two points are `|ray|` apart,
//!    `ns` being the surface normal projected into the section plane.
//!
//! Differences from [`super::surf_rst_const_rad`] (`BRepBlend_SurfRstConstRad`):
//! - `ray` / `dray` come from `tevol->Value` / `tevol->D1` and are multiplied by
//!   `sg1`, which only `Set(Choix)` (`cxx:645-660`) writes; the constant-radius
//!   sibling stores the signed radius once in `Set(Radius, Choix)`.
//! - the third equation is the plain `vref.SquareMagnitude() - ray * ray`
//!   (`cxx:168`, `:265`), not the `(vref / ray)^2 - 1` form scaled by `ray^2`
//!   that `SurfRstConstRad` uses (`cxx:208`) — the two agree algebraically but
//!   not in rounding.
//! - `NbIntervals` / `Intervals` merge the guide's intervals with the law's
//!   (`FusionneIntervalles`, `cxx:52-110`).
//!
//! Translation notes:
//! - `Set(First, Last)` (`cxx:330-334`) narrows both `tguide` and `tevol` with
//!   `Trim(First, Last, 1.e-12)`. Every `Law_*` implementation in OCCT keeps the
//!   parameterisation on the kept interval (`Law_BSpline::Trim` segments the
//!   basis B-spline, `Law_Linear` / `Law_Constant` build the same values), so
//!   `tguide` aliases `guide` as in [`super::surf_rst_const_rad`], while
//!   `tevol` is replaced by `fevol.trim(..)` so the port still narrows the law
//!   the way OCCT does.
//! - `myTConv` is left uninitialised by OCCT's constructor (`cxx:113-131` sets
//!   only `mySShape`); `GetShape` / `Set(TypeSection)` write it before any use,
//!   and the port starts from `Convert_TgtThetaOver2` (the enumerator `0`).
//! - `sg1`, `ray`, `dray` and `choix` are uninitialised in OCCT until
//!   `Set(Choix)` runs; the port starts them at `0` / `0.` and documents it
//!   rather than inventing a quadrant.
//! - `TangentOnS` / `Tangent2dOnS` / `TangentOnRst` / `Tangent2dOnRst` raise
//!   `Standard_DomainError` when `istangent` (`cxx:545-585`); the trait returns
//!   a reference, and every caller in OCCT guards the call with
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

use crate::blend::{BlendAppFunction, BlendFuncSectionShape, BlendPoint, BlendSurfRstFunction};

pub use constants::{DECROCH_DOT_EPS, DECROCH_NT_EPS, EPS, GAUSS_DEFAULT_EPS, SVD_EPS};

/// `class BRepBlend_SurfRstEvolRad : public Blend_SurfRstFunction`
/// (`hxx:41-248`).
pub struct SurfRstEvolRad {
    /// `surf` (`hxx:211`): the surface the section starts from.
    surf: Arc<dyn Surface>,
    /// `surfrst` (`hxx:212`): the surface carrying the restriction.
    surfrst: Arc<dyn Surface>,
    /// `rst` (`hxx:213`): the restriction pcurve on `surfrst`.
    rst: Arc<dyn Curve2d>,
    /// `cons` (`hxx:214`): `rst` seen as a 3D curve on `surfrst`.
    cons: CurveOnSurface,
    /// `guide` (`hxx:215`): the guide line.
    guide: Arc<dyn Curve>,
    /// `tguide` (`hxx:216`): `guide->Trim(First, Last, 1e-12)`; see the module
    /// note -- evaluation is identical, so it aliases `guide`.
    tguide: Arc<dyn Curve>,
    /// `pts` (`hxx:217`).
    pts: GpPnt,
    /// `ptrst` (`hxx:218`).
    ptrst: GpPnt,
    /// `pt2ds` (`hxx:219`).
    pt2ds: GpPnt2d,
    /// `pt2drst` (`hxx:220`).
    pt2drst: GpPnt2d,
    /// `prmrst` (`hxx:221`).
    prmrst: f64,
    /// `istangent` (`hxx:222`).
    istangent: bool,
    /// `tgs` (`hxx:223`).
    tgs: GpVec,
    /// `tg2ds` (`hxx:224`).
    tg2ds: GpVec2d,
    /// `tgrst` (`hxx:225`).
    tgrst: GpVec,
    /// `tg2drst` (`hxx:226`).
    tg2drst: GpVec2d,
    /// `ray` (`hxx:227`): the signed radius at the current parameter.
    ray: f64,
    /// `dray` (`hxx:228`): the derivative of `ray` with respect to the guide
    /// parameter, also signed by `sg1`.
    dray: f64,
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
    /// `tevol` (`hxx:243`): the radius law, trimmed to the guide interval.
    tevol: Arc<dyn LawFunction>,
    /// `fevol` (`hxx:244`): the untrimmed radius law `Set(First, Last)` trims.
    fevol: Arc<dyn LawFunction>,
    /// `sg1` (`hxx:245`): `-1` for the quadrants `1`/`2`, `+1` for `3`/`4`.
    sg1: f64,
}

impl SurfRstEvolRad {
    /// `BRepBlend_SurfRstEvolRad(Surf, SurfRst, Rst, CGuide, Evol)`
    /// (`cxx:113-131`).
    pub fn new(
        surf: Arc<dyn Surface>,
        surf_rst: Arc<dyn Surface>,
        rst: Arc<dyn Curve2d>,
        cguide: Arc<dyn Curve>,
        evol: Arc<dyn LawFunction>,
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
            dray: 0.0,
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
            tevol: evol.clone(),
            fevol: evol,
            // OCCT leaves `sg1` uninitialised until `Set(Choix)`.
            sg1: 0.0,
        }
    }

    /// `Set(SurfRef, RstRef)` (`cxx:302-306`): the face and the pcurve
    /// `Decroch` measures the restriction normal on.
    pub fn set_refs(&mut self, surf_ref: Arc<dyn Surface>, rst_ref: Arc<dyn Curve2d>) {
        self.surfref = Some(surf_ref);
        self.rstref = Some(rst_ref);
    }

    /// `Set(Choix)` (`cxx:645-660`): the sign of the radius follows the
    /// fillet quadrant. The `default:` arm repeats the `1`/`2` arm.
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => self.sg1 = -1.0,
            3 | 4 => self.sg1 = 1.0,
            _ => self.sg1 = -1.0,
        }
    }

    /// `Set(TypeSection)` (`cxx:664-667`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `ray` (`hxx:227`): the signed section radius.
    pub fn ray(&self) -> f64 {
        self.ray
    }

    /// `dray` (`hxx:228`): the signed radius derivative.
    pub fn dray(&self) -> f64 {
        self.dray
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

    /// `surf` (`hxx:211`).
    pub fn surf(&self) -> &Arc<dyn Surface> {
        &self.surf
    }

    /// `surfrst` (`hxx:212`).
    pub fn surf_rst(&self) -> &Arc<dyn Surface> {
        &self.surfrst
    }

    /// `tevol` (`hxx:243`), the radius law currently in force.
    pub fn law(&self) -> &Arc<dyn LawFunction> {
        &self.tevol
    }
}

impl MathFunctionSetWithDerivatives for SurfRstEvolRad {
    /// `NbVariables()` (`cxx:138-141`): `3`.
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:145-148`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:153-171`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:175-231`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d)
    }

    /// `Values(X, F, D)` (`cxx:235-296`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, f, d)
    }
}

impl BlendAppFunction for SurfRstEvolRad {
    /// `Set(Param)` (`cxx:308-325`): the section plane and the radius law.
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
        self.ray = self.sg1 * ray;
        self.dray = self.sg1 * dray;
    }

    /// `Set(First, Last)` (`cxx:330-334`): `tguide = guide->Trim(First, Last,
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

    /// `GetTolerance(Tolerance, Tol)` (`cxx:337-343`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:347-368`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:371-491`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:495-497`).
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

    /// `IsRational()` (`cxx:728-731`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:735-738`).
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:742-747`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:751-769`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.nb_intervals_impl(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:773-797`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.intervals_impl(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:801-805`).
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
    /// (`cxx:807-821`).
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

    /// `Knots(TKnots)` (`cxx:825-827`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:831-833`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`cxx:839-1070`).
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

    /// `Section(P, Poles, Poles2d, Weigths)` (`cxx:1087-1149`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:1076-1085`): `false`.
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

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1151-1166`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendSurfRstFunction for SurfRstEvolRad {
    /// `PointOnS()` (`cxx:501-504`).
    fn point_on_s(&self) -> &GpPnt {
        &self.pts
    }

    /// `PointOnRst()` (`cxx:508-511`).
    fn point_on_rst(&self) -> &GpPnt {
        &self.ptrst
    }

    /// `Pnt2dOnS()` (`cxx:514-517`).
    fn pnt2d_on_s(&self) -> &GpPnt2d {
        &self.pt2ds
    }

    /// `Pnt2dOnRst()` (`cxx:521-524`).
    fn pnt2d_on_rst(&self) -> &GpPnt2d {
        &self.pt2drst
    }

    /// `ParameterOnRst()` (`cxx:529-532`).
    fn parameter_on_rst(&self) -> f64 {
        self.prmrst
    }

    /// `IsTangencyPoint()` (`cxx:535-538`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS()` (`cxx:541-550`); OCCT raises `Standard_DomainError` when
    /// `istangent` (see the module note).
    fn tangent_on_s(&self) -> &GpVec {
        &self.tgs
    }

    /// `Tangent2dOnS()` (`cxx:553-564`).
    fn tangent2d_on_s(&self) -> &GpVec2d {
        &self.tg2ds
    }

    /// `TangentOnRst()` (`cxx:567-576`).
    fn tangent_on_rst(&self) -> &GpVec {
        &self.tgrst
    }

    /// `Tangent2dOnRst()` (`cxx:579-588`).
    fn tangent2d_on_rst(&self) -> &GpVec2d {
        &self.tg2drst
    }

    /// `Decroch(Sol, NS, TgS)` (`cxx:588-640`).
    fn decroch(&self, sol: &MathVector, ns: &mut GpVec, tg_s: &mut GpVec) -> bool {
        self.compute_decroch(sol, ns, tg_s)
    }
}
