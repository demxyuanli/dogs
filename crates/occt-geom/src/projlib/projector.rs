//! `ProjLib_Projector` — the base of the projector family.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_Projector.hxx:31-186`,
//! `.cxx:34-251` (OCCT 8.0.0).
//!
//! The C++ class stores the projection result of one curve in the members
//! `myType` / `myLin` / `myCirc` / `myElips` / `myHypr` / `myParab` /
//! `myBezier` / `myBSpline` / `myIsPeriodic` / `isDone`, and lets every
//! subclass write them directly. The port keeps that state in
//! [`ProjectorState`] and exposes it through the [`Projector`] trait, whose
//! default `Project` overloads are the bodies of the C++ virtuals that set
//! `myType = GeomAbs_OtherCurve` (`cxx:145-176`).
//!
//! `myCirc` and friends are read back through `Line()` / `Circle()` / ...,
//! which in C++ raise `Standard_NoSuchObject` when `myType` does not match
//! (`cxx:88-141`); the port answers `Option<&...>` in that case.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{
    GpCirc, GpCirc2d, GpDir, GpDir2d, GpElips, GpElips2d, GpHypr, GpHypr2d, GpLin, GpLin2d,
    GpParab, GpParab2d, GpPnt2d, GpVec, GpVec2d,
};
use occt_core::kernel::geomabs::CurveType;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

/// The data members of `ProjLib_Projector` (`hxx:157-171`).
#[derive(Clone)]
pub struct ProjectorState {
    /// `myType` (`hxx:160`), `GeomAbs_BSplineCurve` at construction
    /// (`cxx:38`).
    pub curve_type: CurveType,
    /// `myIsPeriodic` (`hxx:159`), `false` at construction (`cxx:35`).
    pub periodic: bool,
    /// `isDone` (`hxx:158`), `false` at construction (`cxx:37`).
    pub done: bool,
    /// `myLin` (`hxx:161`).
    pub lin: GpLin2d,
    /// `myCirc` (`hxx:162`).
    pub circ: GpCirc2d,
    /// `myElips` (`hxx:163`).
    pub elips: GpElips2d,
    /// `myHypr` (`hxx:164`).
    pub hypr: GpHypr2d,
    /// `myParab` (`hxx:165`).
    pub parab: GpParab2d,
    /// `myBezier` (`hxx:167`).
    pub bezier: Option<Arc<Geom2dBezierCurve>>,
    /// `myBSpline` (`hxx:168`).
    pub bspline: Option<Arc<Geom2dBSplineCurve>>,
}

impl ProjectorState {
    /// `ProjLib_Projector()` (`cxx:34-39`).
    pub fn new() -> Self {
        Self {
            curve_type: CurveType::BSplineCurve,
            periodic: false,
            done: false,
            lin: GpLin2d::default(),
            circ: GpCirc2d::default(),
            elips: GpElips2d::default(),
            hypr: GpHypr2d::default(),
            parab: GpParab2d::default(),
            bezier: None,
            bspline: None,
        }
    }

    /// `ProjLib_Projector::UFrame(CFirst, CLast, UFirst, Period)`
    /// (`cxx:179-199`): only the line case moves the result, by the multiple
    /// of `Period` that brings its first point into
    /// `[UFirst, UFirst + Period)`.
    pub fn u_frame(&mut self, c_first: f64, _c_last: f64, u_first: f64, period: f64) {
        if self.curve_type != CurveType::Line {
            return;
        }
        let p_first = clib::line2d_value(&self.lin, c_first);
        let u = p_first.x();
        let new_u = clib::in_period(u, u_first, u_first + period);
        self.lin.translate_vec(&GpVec2d::new(new_u - u, 0.0));
    }

    /// `ProjLib_Projector::VFrame(CFirst, CLast, VFirst, Period)`
    /// (`cxx:204-224`).
    pub fn v_frame(&mut self, c_first: f64, _c_last: f64, v_first: f64, period: f64) {
        if self.curve_type != CurveType::Line {
            return;
        }
        let p_first = clib::line2d_value(&self.lin, c_first);
        let v = p_first.y();
        let new_v = clib::in_period(v, v_first, v_first + period);
        self.lin.translate_vec(&GpVec2d::new(0.0, new_v - v));
    }
}

impl Default for ProjectorState {
    fn default() -> Self {
        Self::new()
    }
}

/// `class ProjLib_Projector` (`hxx:31-186`): the state plus the virtual
/// `Project` family.
///
/// A projector is created empty (`Init` marks it `GeomAbs_OtherCurve`), then
/// `Project` is called with the 3D curve; `IsDone()` says whether an analytic
/// form was found. `ProjLib_ProjectedCurve` sends any projector that answers
/// `IsDone() == false` to its approximation branch.
pub trait Projector {
    /// The `ProjLib_Projector` data members.
    fn state(&self) -> &ProjectorState;
    /// Mutable access to the same members (the C++ subclasses write them).
    fn state_mut(&mut self) -> &mut ProjectorState;

    /// `IsDone()` (`cxx:47-50`).
    fn is_done(&self) -> bool {
        self.state().done
    }

    /// `Done()` (`cxx:54-57`).
    fn done(&mut self) {
        self.state_mut().done = true;
    }

    /// `GetType()` (`cxx:61-64`).
    fn get_type(&self) -> CurveType {
        self.state().curve_type
    }

    /// `SetType(Type)` (`cxx:68-71`).
    fn set_type(&mut self, curve_type: CurveType) {
        self.state_mut().curve_type = curve_type;
    }

    /// `IsPeriodic()` (`cxx:75-78`).
    fn is_periodic(&self) -> bool {
        self.state().periodic
    }

    /// `SetPeriodic()` (`cxx:82-85`).
    fn set_periodic(&mut self) {
        self.state_mut().periodic = true;
    }

    /// `Line()` (`cxx:88-95`): `None` where OCCT raises `NoSuchObject`.
    fn line(&self) -> Option<&GpLin2d> {
        let s = self.state();
        (s.curve_type == CurveType::Line).then_some(&s.lin)
    }

    /// `Circle()` (`cxx:99-106`).
    fn circle(&self) -> Option<&GpCirc2d> {
        let s = self.state();
        (s.curve_type == CurveType::Circle).then_some(&s.circ)
    }

    /// `Ellipse()` (`cxx:110-117`).
    fn ellipse(&self) -> Option<&GpElips2d> {
        let s = self.state();
        (s.curve_type == CurveType::Ellipse).then_some(&s.elips)
    }

    /// `Hyperbola()` (`cxx:121-128`).
    fn hyperbola(&self) -> Option<&GpHypr2d> {
        let s = self.state();
        (s.curve_type == CurveType::Hyperbola).then_some(&s.hypr)
    }

    /// `Parabola()` (`cxx:132-139`).
    fn parabola(&self) -> Option<&GpParab2d> {
        let s = self.state();
        (s.curve_type == CurveType::Parabola).then_some(&s.parab)
    }

    /// `Project(const gp_Lin&)` (`cxx:145-148`), the inherited body.
    fn project_lin(&mut self, _l: &GpLin) {
        self.state_mut().curve_type = CurveType::OtherCurve;
    }

    /// `Project(const gp_Circ&)` (`cxx:152-155`).
    fn project_circ(&mut self, _c: &GpCirc) {
        self.state_mut().curve_type = CurveType::OtherCurve;
    }

    /// `Project(const gp_Elips&)` (`cxx:158-161`).
    fn project_elips(&mut self, _e: &GpElips) {
        self.state_mut().curve_type = CurveType::OtherCurve;
    }

    /// `Project(const gp_Parab&)` (`cxx:164-167`).
    fn project_parab(&mut self, _p: &GpParab) {
        self.state_mut().curve_type = CurveType::OtherCurve;
    }

    /// `Project(const gp_Hypr&)` (`cxx:170-173`).
    fn project_hypr(&mut self, _h: &GpHypr) {
        self.state_mut().curve_type = CurveType::OtherCurve;
    }

    /// `UFrame(CFirst, CLast, UFirst, Period)` (`cxx:179-199`): only the line
    /// case moves the result, by the multiple of `Period` that brings its
    /// first point into `[UFirst, UFirst + Period)`.
    fn u_frame(&mut self, c_first: f64, c_last: f64, u_first: f64, period: f64) {
        self.state_mut().u_frame(c_first, c_last, u_first, period);
    }

    /// `VFrame(CFirst, CLast, VFirst, Period)` (`cxx:204-224`).
    fn v_frame(&mut self, c_first: f64, c_last: f64, v_first: f64, period: f64) {
        self.state_mut().v_frame(c_first, c_last, v_first, period);
    }

    /// `SetBezier(C)` (`cxx:228-231`).
    fn set_bezier(&mut self, c: Arc<Geom2dBezierCurve>) {
        self.state_mut().bezier = Some(c);
    }

    /// `Bezier()` (`cxx:234-237`).
    fn bezier(&self) -> Option<Arc<Geom2dBezierCurve>> {
        self.state().bezier.clone()
    }

    /// `SetBSpline(C)` (`cxx:241-244`).
    fn set_bspline(&mut self, c: Arc<Geom2dBSplineCurve>) {
        self.state_mut().bspline = Some(c);
    }

    /// `BSpline()` (`cxx:247-250`).
    fn bspline(&self) -> Option<Arc<Geom2dBSplineCurve>> {
        self.state().bspline.clone()
    }
}

/// `gp_Dir2d(gp_Vec2d(P1, P2))` (`gp_Vec2d.cxx:29-32`, `gp_Dir2d.hxx:304-320`):
/// the unit direction from `p1` to `p2`. `None` where `gp_Dir2d` raises
/// `Standard_ConstructionError` (coincident points).
///
/// Shared by the sphere / torus projectors, which build the 2D line direction
/// out of two evaluated points.
pub(crate) fn dir2d_from_pnts(p1: &GpPnt2d, p2: &GpPnt2d) -> Option<GpDir2d> {
    GpDir2d::new(p2.x() - p1.x(), p2.y() - p1.y()).ok()
}

/// `gp_Dir(gp_Vec)` — the normalizing direction constructor
/// (`gp_Dir.hxx:364`). `None` on a zero vector, where OCCT raises
/// `Standard_ConstructionError`.
pub(crate) fn dir_from_vec(v: &GpVec) -> Option<GpDir> {
    GpDir::from_xyz(v.xyz()).ok()
}
