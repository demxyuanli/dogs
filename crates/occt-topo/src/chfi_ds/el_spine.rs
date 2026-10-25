//! `ChFiDS_ElSpine` — the elementary spine of one edge of the guideline.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_ElSpine.hxx:45-163`,
//! `ChFiDS_ElSpine.cxx:30-384` (OCCT 8.0.0).
//!
//! OCCT makes this class an `Adaptor3d_Curve` over an inner
//! `GeomAdaptor_Curve`, so every `Adaptor3d_Curve` query below is an inline
//! delegation (`cxx:105-199`, `cxx:336-384`) and the stored `[pfirst, plast]`
//! only overrides `FirstParameter` / `LastParameter` (`cxx:76-88`).
//!
//! UNPORTED, with their OCCT location and the reason:
//! - `Trim` (`cxx:126-131`): the port's [`Curve`] trait carries no
//!   `GeomAdaptor_Curve::Trim` (the overload reuses the same adaptor when the
//!   range is unchanged), and no ported `TKFillet` code calls
//!   `ElSpine::Trim`,
//! - `Bezier` / `BSpline` (`cxx:374-384`): they `down_cast` the inner curve and
//!   build an adaptor of the converted curve; no ported `TKFillet` caller.
//!   `SetOrigin` (`cxx:230-244`) **is** ported: [`ElSpine::set_origin`] uses
//!   [`Curve::as_bspline_curve`] (`occt-geom/src/curve.rs`) plus
//!   `GeomBSplineCurve::set_origin_u`
//!   (`occt-geom/src/bspline_curve/periodic.rs:162`), the same re-knotting as
//!   `Geom_BSplineCurve::SetOrigin(U, Tol)` (`Geom_BSplineCurve.cxx:913-970`),
//!
//! Divergences that are not heuristics:
//! - `occ::handle<T>` becomes `Option<Arc<T>>` (null handle) or `Arc<T>`,
//! - the accessors that dereference the null `GeomAdaptor_Curve` handle in
//!   OCCT (`Standard_NullObject`) answer `Option` instead; the [`Curve`] impl,
//!   which cannot, unwraps with the same message.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCirc, GpElips, GpHypr, GpLin, GpParab, GpPnt, GpVec};
use occt_core::precision::INFINITE;
use occt_geom::curve::Curve;

use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};

use super::surf_data::SurfData;

/// Message OCCT produces by dereferencing the null `myCurve` handle.
const NO_CURVE: &str = "ChFiDS_ElSpine: no curve loaded (Standard_NullObject)";

/// Source: `ChFiDS_ElSpine.hxx:46-163`.
#[derive(Clone)]
pub struct ElSpine {
    /// `GeomAdaptor_Curve curve` — default-constructed (null) until `SetCurve`.
    curve: Option<Arc<dyn Curve>>,
    /// `ptfirst`
    pt_first: GpPnt,
    /// `ptlast`
    pt_last: GpPnt,
    /// `tgfirst`
    tg_first: GpVec,
    /// `tglast`
    tg_last: GpVec,
    /// `VerticesWithTangents` — `NCollection_Sequence<gp_Ax1>`, 1-based in OCCT.
    vertices_with_tangents: Vec<GpAx1>,
    /// `previous`
    previous: Option<Arc<SurfData>>,
    /// `next`
    next: Option<Arc<SurfData>>,
    /// `pfirst`
    p_first: f64,
    /// `plast`
    p_last: f64,
    /// `period`
    period: f64,
    /// `periodic`
    periodic: bool,
    /// `pfirstsav`
    p_first_sav: f64,
    /// `plastsav`
    p_last_sav: f64,
}

impl std::fmt::Debug for ElSpine {
    /// `dyn Curve` carries no `Debug`; the curve prints as loaded / not loaded.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ElSpine")
            .field("has_curve", &self.curve.is_some())
            .field("pt_first", &self.pt_first)
            .field("pt_last", &self.pt_last)
            .field("nb_vertices_with_tangents", &self.vertices_with_tangents.len())
            .field("has_previous", &self.previous.is_some())
            .field("has_next", &self.next.is_some())
            .field("p_first", &self.p_first)
            .field("p_last", &self.p_last)
            .field("period", &self.period)
            .field("periodic", &self.periodic)
            .field("p_first_sav", &self.p_first_sav)
            .field("p_last_sav", &self.p_last_sav)
            .finish()
    }
}

impl Default for ElSpine {
    fn default() -> Self {
        Self::new()
    }
}

impl ElSpine {
    /// `ChFiDS_ElSpine()` (`cxx:35-45`): the two saved parameters start at
    /// `Precision::Infinite()`, everything else at `0.` / `false`.
    pub fn new() -> Self {
        Self {
            curve: None,
            pt_first: GpPnt::zero(),
            pt_last: GpPnt::zero(),
            tg_first: GpVec::zero(),
            tg_last: GpVec::zero(),
            vertices_with_tangents: Vec::new(),
            previous: None,
            next: None,
            p_first: 0.0,
            p_last: 0.0,
            period: 0.0,
            periodic: false,
            p_first_sav: INFINITE,
            p_last_sav: INFINITE,
        }
    }

    /// `ShallowCopy()` (`cxx:48-72`): a copy that shares the inner handle and
    /// the two `SurfData` handles. OCCT returns it as an `Adaptor3d_Curve`
    /// handle; the port returns the concrete copy.
    pub fn shallow_copy(&self) -> Self {
        self.clone()
    }

    /// `FirstParameter()` (`cxx:76-80`): the stored `pfirst`, not the inner
    /// curve's.
    pub fn first_parameter(&self) -> f64 {
        self.p_first
    }

    /// `LastParameter()` (`cxx:83-87`): the stored `plast`.
    pub fn last_parameter(&self) -> f64 {
        self.p_last
    }

    /// `GetSavedFirstParameter()` (`cxx:90-93`).
    pub fn get_saved_first_parameter(&self) -> f64 {
        self.p_first_sav
    }

    /// `GetSavedLastParameter()` (`cxx:96-99`).
    pub fn get_saved_last_parameter(&self) -> f64 {
        self.p_last_sav
    }

    /// `Continuity()` (`cxx:105-108`) -> `GeomAdaptor_Curve::Continuity()`.
    pub fn continuity(&self) -> Option<u8> {
        self.curve.as_ref().map(|c| c.continuity())
    }

    /// `NbIntervals(S)` (`cxx:112-115`).
    pub fn nb_intervals(&self, s: u8) -> Option<i32> {
        self.curve.as_ref().map(|c| c.nb_intervals(s))
    }

    /// `Intervals(T, S)` (`cxx:118-122`): the `S`-continuity knot sequence,
    /// including the range ends (`NbIntervals + 1` values in OCCT).
    pub fn intervals(&self, s: u8) -> Option<Vec<f64>> {
        self.curve.as_ref().map(|c| c.parameter_intervals(s))
    }

    /// `Resolution(R3d)` (`cxx:134-137`).
    pub fn resolution(&self, r3d: f64) -> Option<f64> {
        self.curve.as_ref().map(|c| c.resolution(r3d))
    }

    /// `GetType()` (`cxx:141-144`) -> `GeomAdaptor_Curve::GetType()`, i.e. the
    /// `myTypeCurve` the adaptor computed in `Load` (`GeomAdaptor_Curve.cxx:238-313`).
    pub fn get_type(&self) -> Option<GeomAbsCurveType> {
        self.curve.as_ref().map(|c| curve_kind(c.as_ref()))
    }

    /// `IsPeriodic()` (`cxx:149-152`): the stored flag, not the inner curve's.
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// `SetPeriodic(I)` (`cxx:155-159`): sets the flag **and** derives
    /// `period` from the stored range.
    pub fn set_periodic(&mut self, i: bool) {
        self.periodic = i;
        self.period = self.p_last - self.p_first;
    }

    /// `Period()` (`cxx:162-171`): OCCT throws `Standard_Failure` when the flag
    /// is false; the port answers `None`.
    pub fn period(&self) -> Option<f64> {
        if !self.periodic {
            return None;
        }
        Some(self.period)
    }

    /// `SetOrigin(O)` (`cxx:230-244`): move the origin of a periodic spine to
    /// `O`.
    ///
    /// OCCT throws `Standard_Failure("Elspine non periodique")` when the flag
    /// is false, and does nothing when the inner curve does not `down_cast` to
    /// a `Geom_BSplineCurve` (`cxx:236-241`). The port answers the throw as
    /// `Err` and keeps the silent non-B-spline case silent.
    ///
    /// `Geom_BSplineCurve::SetOrigin(U, Tol)` (`Geom_BSplineCurve.cxx:913-970`)
    /// re-knots the curve in place, so the stored handle is replaced by the
    /// re-knotted clone, exactly like `curve.Load(bs)` (`cxx:240`) reloads the
    /// adaptor.
    pub fn set_origin(&mut self, o: f64) -> Result<(), String> {
        if !self.periodic {
            return Err("Standard_Failure: Elspine non periodique".into());
        }
        let Some(curve) = self.curve.as_ref() else {
            return Ok(());
        };
        let Some(bs) = curve.as_bspline_curve() else {
            return Ok(());
        };
        let mut bs = bs.clone();
        bs.set_origin_u(o, occt_core::precision::PCONFUSION)
            .map_err(str::to_string)?;
        self.curve = Some(Arc::new(bs));
        Ok(())
    }

    /// `EvalD0(theAbsC)` (`cxx:174-177`).
    pub fn eval_d0(&self, u: f64) -> Option<GpPnt> {
        self.curve.as_ref().map(|c| c.d0(u))
    }

    /// `EvalD1(theAbsC)` (`cxx:181-184`).
    pub fn eval_d1(&self, u: f64) -> Option<(GpPnt, GpVec)> {
        self.curve.as_ref().map(|c| c.d1(u))
    }

    /// `EvalD2(theAbsC)` (`cxx:188-191`).
    pub fn eval_d2(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        self.curve.as_ref().map(|c| c.d2(u))
    }

    /// `EvalD3(theAbsC)` (`cxx:195-198`).
    pub fn eval_d3(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
        self.curve.as_ref().map(|c| c.d3(u))
    }

    /// `FirstParameter(P)` (`cxx:205-208`): rewrites the stored `pfirst`
    /// without touching the inner curve.
    pub fn set_first_parameter(&mut self, p: f64) {
        self.p_first = p;
    }

    /// `LastParameter(P)` (`cxx:210-213`).
    pub fn set_last_parameter(&mut self, p: f64) {
        self.p_last = p;
    }

    /// `SaveFirstParameter()` (`cxx:216-219`).
    pub fn save_first_parameter(&mut self) {
        self.p_first_sav = self.p_first;
    }

    /// `SaveLastParameter()` (`cxx:221-224`).
    pub fn save_last_parameter(&mut self) {
        self.p_last_sav = self.p_last;
    }

    /// `SetFirstPointAndTgt(P, T)` (`cxx:246-250`).
    pub fn set_first_point_and_tgt(&mut self, p: GpPnt, t: GpVec) {
        self.pt_first = p;
        self.tg_first = t;
    }

    /// `SetLastPointAndTgt(P, T)` (`cxx:254-258`).
    pub fn set_last_point_and_tgt(&mut self, p: GpPnt, t: GpVec) {
        self.pt_last = p;
        self.tg_last = t;
    }

    /// `AddVertexWithTangent(anAx1)` (`cxx:262-265`).
    pub fn add_vertex_with_tangent(&mut self, an_ax1: GpAx1) {
        self.vertices_with_tangents.push(an_ax1);
    }

    /// `FirstPointAndTgt(P, T)` (`cxx:268-272`).
    pub fn first_point_and_tgt(&self) -> (GpPnt, GpVec) {
        (self.pt_first, self.tg_first)
    }

    /// `LastPointAndTgt(P, T)` (`cxx:276-281`).
    pub fn last_point_and_tgt(&self) -> (GpPnt, GpVec) {
        (self.pt_last, self.tg_last)
    }

    /// `NbVertices()` (`cxx:285-288`).
    pub fn nb_vertices(&self) -> usize {
        self.vertices_with_tangents.len()
    }

    /// `VertexWithTangent(Index)` (`cxx:292-295`), `Index` 1-based. OCCT
    /// raises `Standard_OutOfRange`; the port answers `None`.
    pub fn vertex_with_tangent(&self, index: usize) -> Option<&GpAx1> {
        if index == 0 {
            return None;
        }
        self.vertices_with_tangents.get(index - 1)
    }

    /// `SetCurve(C)` (`cxx:299-302`) -> `curve.Load(C)`.
    pub fn set_curve(&mut self, c: Arc<dyn Curve>) {
        self.curve = Some(c);
    }

    /// `Previous()` (`cxx:306-309`).
    pub fn previous(&self) -> Option<&Arc<SurfData>> {
        self.previous.as_ref()
    }

    /// `ChangePrevious()` (`cxx:313-316`).
    pub fn previous_mut(&mut self) -> &mut Option<Arc<SurfData>> {
        &mut self.previous
    }

    /// `Next()` (`cxx:320-323`).
    pub fn next(&self) -> Option<&Arc<SurfData>> {
        self.next.as_ref()
    }

    /// `ChangeNext()` (`cxx:327-330`).
    pub fn next_mut(&mut self) -> &mut Option<Arc<SurfData>> {
        &mut self.next
    }

    /// `Line()` (`cxx:339-342`): OCCT's `GeomAdaptor_Curve::Line` down-casts
    /// and throws `Standard_DomainError` on a type mismatch.
    pub fn line(&self) -> Option<GpLin> {
        self.curve.as_ref().and_then(|c| c.gp_line())
    }

    /// `Circle()` (`cxx:346-349`).
    pub fn circle(&self) -> Option<GpCirc> {
        self.curve.as_ref().and_then(|c| c.gp_circ())
    }

    /// `Ellipse()` (`cxx:353-356`).
    pub fn ellipse(&self) -> Option<GpElips> {
        self.curve.as_ref().and_then(|c| c.gp_ellipse())
    }

    /// `Hyperbola()` (`cxx:360-363`).
    pub fn hyperbola(&self) -> Option<GpHypr> {
        self.curve.as_ref().and_then(|c| c.gp_hyperbola())
    }

    /// `Parabola()` (`cxx:367-370`).
    pub fn parabola(&self) -> Option<GpParab> {
        self.curve.as_ref().and_then(|c| c.gp_parabola())
    }
}

/// `ChFiDS_ElSpine` *is* an `Adaptor3d_Curve` (`hxx:46`); this is that view for
/// the algorithms that take one (`BRepBlend`, `Approx`).
///
/// `transform` / `reverse` are no-ops: `Adaptor3d_Curve` has no such operations
/// (only `Geom_Curve` does), exactly like the other adaptor views in the port
/// (`occt-topo/src/shhealing/transfer_params.rs:169-171`).
impl Curve for ElSpine {
    fn d0(&self, u: f64) -> GpPnt {
        self.eval_d0(u).expect(NO_CURVE)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.eval_d1(u).expect(NO_CURVE)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.eval_d2(u).expect(NO_CURVE)
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        self.eval_d3(u).expect(NO_CURVE)
    }

    fn first_parameter(&self) -> f64 {
        ElSpine::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        ElSpine::last_parameter(self)
    }

    fn continuity(&self) -> u8 {
        ElSpine::continuity(self).expect(NO_CURVE)
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        ElSpine::nb_intervals(self, continuity).expect(NO_CURVE)
    }

    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        ElSpine::intervals(self, continuity).expect(NO_CURVE)
    }

    fn resolution(&self, r3d: f64) -> f64 {
        ElSpine::resolution(self, r3d).expect(NO_CURVE)
    }

    fn is_periodic(&self) -> bool {
        ElSpine::is_periodic(self)
    }

    fn period(&self) -> f64 {
        // OCCT raises `Standard_Failure` here when the flag is false
        // (`cxx:164-168`); a `Curve` accessor cannot, so the stored value is
        // returned and `is_periodic` stays the guard.
        self.period
    }

    fn circle_radius(&self) -> Option<f64> {
        self.curve.as_ref().and_then(|c| c.circle_radius())
    }

    fn gp_line(&self) -> Option<GpLin> {
        self.line()
    }

    fn gp_circ(&self) -> Option<GpCirc> {
        self.circle()
    }

    fn gp_ellipse(&self) -> Option<GpElips> {
        self.ellipse()
    }

    fn gp_hyperbola(&self) -> Option<GpHypr> {
        self.hyperbola()
    }

    fn gp_parabola(&self) -> Option<GpParab> {
        self.parabola()
    }

    fn is_line(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_line())
    }

    fn nurbs_degree(&self) -> Option<usize> {
        self.curve.as_ref().and_then(|c| c.nurbs_degree())
    }

    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bspline_poles())
    }

    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bspline_weights())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bezier_poles())
    }

    fn bezier_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bezier_weights())
    }

    fn transform(&mut self, _t: &occt_core::gp::GpTrsf) {}

    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
