//! Port of `ShapeAnalysis_TransferParameters` and
//! `ShapeAnalysis_TransferParametersProj` (TKShHealing, `ShapeAnalysis`
//! package: `ShapeAnalysis_TransferParameters.cxx` / `.hxx`,
//! `ShapeAnalysis_TransferParametersProj.cxx` / `.hxx`; the family is exactly
//! these two classes, `ShapeAnalysis/FILES.cmake:38-41`).
//!
//! The tool transfers parameters between the 3D curve of an edge and its
//! pcurve. The base class uses the linear map
//! `T2d = myShift + myScale * T3d` (`hxx:34-40`); the `Proj` subclass projects
//! instead whenever the edge is not SameParameter or its tolerance is not below
//! `MaxTolerance()` (`Proj.cxx:112-116`). Call sites ported next to this module:
//! `ShapeFix_Wire::FixNotchedEdges` (`ShapeFix_Wire.cxx:4030-4065`) and, later,
//! `FixSelfIntersectingEdge` / `FixIntersectingEdges` (`cxx:2449`, `:2597`).
//!
//! `ShapeBuild_Edge::CopyRanges` (`ShapeBuild_Edge.cxx:206-334`) is ported here
//! as [`copy_ranges`], together with the `IsPeriodic` unwrap
//! (`ShapeBuild_Edge.cxx:164-204`) and `AdjustByPeriod`
//! (`ShapeBuild_Edge.cxx:148-162`); both `TransferRange` overloads write
//! through it.
//!
//! UNPORTED in this module:
//! * `ShapeAnalysis_TransferParametersProj::CopyNMVertex` (`Proj.cxx:562-710`
//!   and `:715-803`): both overloads walk the vertex's
//!   `BRep_PointRepresentation` list (`BRep_TVertex::ChangePoints`,
//!   `BRep_PointOnCurve` / `BRep_PointOnSurface` / `BRep_PointOnCurveOnSurface`)
//!   and rebuild it with `BRep_Builder::UpdateVertex`. This port stores vertex
//!   geometry as a single point in the `GeometryRegistry` (no point
//!   representations), so there is nothing to copy; none of the three wire
//!   passes in scope calls it (`ShapeUpgrade_WireDivide.cxx:852`,
//!   `ShapeFix_Wire_1.cxx:813`, `ShapeFix_Wireframe.cxx:1051`,
//!   `ShapeFix_ComposeShell.cxx:3241` do, none of which is ported).
//! * `CorrectParameter`'s `Geom2d_BSplineCurve` knot snap
//!   (`Proj.cxx:268-279`) - see [`correct_parameter`].
//! * `myLocation` (`Proj.cxx:97`, `:199`, `:209`, `:318`, `:322`, `:472`): the
//!   port's healing path keeps locations identity, exactly like the rest of
//!   `shhealing` (`wire_fix.rs` reads edges through `BRepTool::edge_curve` /
//!   `curve_on_surface_oriented`, never the `_world` variants), so every
//!   `Transformed(myLocation)` / `Inverted()` is a no-op here.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{epsilon, Precision, CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

use super::adjust_by_period;
use super::shape_analysis_curve::{next_project, project_adaptor, project_range};

/// `GeomAdaptor_Curve(C, first, last)` (`GeomAdaptor_Curve.cxx:239-245`,
/// `cxx:679-691`): the same parameter space as `curve`, restricted to
/// `[first, last]`. `Value(u)` is `curve.d0(u)`, never a rescaled one, which is
/// what `NextProject` / `Project` expect for the `GeomAdaptor_Curve GAC(C3d,
/// first, last)` built at `Proj.cxx:400` (and for the COS adaptor built from a
/// representation range at `Proj.cxx:462-464`).
struct CurveRange {
    curve: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl CurveRange {
    fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { curve, first, last }
    }
}

impl Curve for CurveRange {
    fn d0(&self, u: f64) -> GpPnt {
        self.curve.d0(u)
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.curve.d1(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.curve.d2(u)
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }
    fn period(&self) -> f64 {
        self.curve.period()
    }
    fn continuity(&self) -> u8 {
        self.curve.continuity()
    }
    fn circle_radius(&self) -> Option<f64> {
        self.curve.circle_radius()
    }
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> {
        self.curve.gp_circ()
    }
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> {
        self.curve.gp_ellipse()
    }
    fn is_line(&self) -> bool {
        self.curve.is_line()
    }
    fn nurbs_degree(&self) -> Option<usize> {
        self.curve.nurbs_degree()
    }
    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curve.bspline_poles()
    }
    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curve.bspline_weights()
    }
    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curve.bezier_poles()
    }
    fn resolution(&self, r3d: f64) -> f64 {
        self.curve.resolution(r3d)
    }
    fn transform(&mut self, _t: &occt_core::gp::GpTrsf) {}
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(CurveRange::new(self.curve.clone(), self.first, self.last))
    }
}

/// `ShapeBuild_Edge::IsPeriodic(handle(Geom_Curve))`
/// (`ShapeBuild_Edge.cxx:164-183`): unwrap `Geom_OffsetCurve` /
/// `Geom_TrimmedCurve` down to the basis, then `IsPeriodic`. Returns the basis
/// period with the basis `FirstParameter()` / `LastParameter()` the caller
/// reads at `cxx:301-303`.
///
/// The `Geom_OffsetCurve` unwrap (`cxx:170-177`) has no `Curve` counterpart in
/// the port (no offset query), so only the trim unwrap runs
/// (`Curve::untrimmed_basis`, `GeomAdaptor_Curve.cxx:252-254`).
fn periodic_basis(curve: &Arc<dyn Curve>) -> Option<(f64, f64, f64)> {
    let (basis, first, last) = match curve.untrimmed_basis() {
        Some((basis, f, l)) => (basis, f, l),
        None => (curve.clone(), curve.first_parameter(), curve.last_parameter()),
    };
    if basis.is_periodic() {
        Some((basis.period(), first, last))
    } else {
        None
    }
}

/// `ShapeBuild_Edge::IsPeriodic(handle(Geom2d_Curve))`
/// (`ShapeBuild_Edge.cxx:185-204`) plus the range read at `cxx:311-317`.
fn periodic_basis_2d(curve: &Arc<dyn Curve2d>) -> Option<(f64, f64, f64)> {
    let mut basis: &dyn Curve2d = curve.as_ref();
    // `Geom2d_OffsetCurve` / `Geom2d_TrimmedCurve` unwrap (`cxx:191-201`).
    while let Some(inner) = basis.trimmed_basis().or_else(|| basis.offset_basis()) {
        basis = inner;
    }
    if basis.is_periodic() {
        Some((basis.period(), basis.first_parameter(), basis.last_parameter()))
    } else {
        None
    }
}

/// `ShapeBuild_Edge::CopyRanges(toedge, fromedge, alpha, beta)`
/// (`ShapeBuild_Edge.cxx:206-334`).
///
/// Walks the `fromedge` representations (`cxx:226-232`) and writes the matching
/// representation of `toedge` (`cxx:264-283`): the `BRep_Curve3D` range
/// (`cxx:285-289`, `cxx:326-331`) and every `BRep_CurveOnSurface` range
/// (`cxx:310-318`, `cxx:527-532`). A periodic basis whose shifted range leaves
/// its own range is re-shifted by a whole period, and `SameRange` /
/// `SameParameter` are cleared on `toedge` (`cxx:291-328`).
///
/// The port matches representations by the face key of the pcurve map, which IS
/// the surface identity `surface != toGC->Surface() || L != toGC->Location()`
/// tested at `cxx:281`.
pub(crate) fn copy_ranges(to: &Edge, from: &Edge, alpha: f64, beta: f64) {
    let reg = GeometryRegistry::global();
    let Some(from_geom) = reg.edge_geom(&from.0) else {
        return;
    };
    let Some(mut to_geom) = reg.edge_geom(&to.0) else {
        return;
    };

    // `cxx:234-245`: the `BRep_Curve3D` representation (`isC3d`).
    {
        let first = from_geom.first; // `cxx:285-286`
        let last = from_geom.last;
        let len = last - first;
        let mut new_first = first + alpha * len; // `cxx:289-290`
        let mut new_last = first + beta * len;
        if let Some((period, basis_first, basis_last)) = periodic_basis(&to_geom.curve) {
            // `cxx:321-328`.
            if ((new_first - basis_first).abs() > PCONFUSION && new_first < basis_first)
                || new_first >= basis_last
            {
                let shift = adjust_by_period(new_first, 0.5 * (basis_first + basis_last), period);
                new_first += shift;
                new_last += shift;
                reg.set_same_range(&to.0, false);
                reg.set_same_parameter(&to.0, false);
            }
        }
        to_geom.first = new_first;
        to_geom.last = new_last;
    }

    // `cxx:248-282`: every `BRep_CurveOnSurface` representation.
    let from_keys: Vec<usize> = from_geom.pcurves.keys().copied().collect();
    for key in from_keys {
        if !to_geom.pcurves.contains_key(&key) {
            continue; // `cxx:281`
        }
        let (first, last) = from_geom
            .pcurve_ranges
            .get(&key)
            .copied()
            .unwrap_or((from_geom.first, from_geom.last));
        let len = last - first;
        let mut new_first = first + alpha * len;
        let mut new_last = first + beta * len;
        let pcurve = to_geom.pcurves.get(&key).and_then(|v| v.first().cloned());
        if let Some((period, basis_first, basis_last)) = pcurve.and_then(|pc| periodic_basis_2d(&pc))
        {
            if ((new_first - basis_first).abs() > PCONFUSION && new_first < basis_first)
                || new_first >= basis_last
            {
                let shift = adjust_by_period(new_first, 0.5 * (basis_first + basis_last), period);
                new_first += shift;
                new_last += shift;
                reg.set_same_range(&to.0, false);
                reg.set_same_parameter(&to.0, false);
            }
        }
        to_geom.pcurve_ranges.insert(key, (new_first, new_last));
    }

    reg.set_edge(&to.0, to_geom);
}

/// `CorrectParameter` (`ShapeAnalysis_TransferParametersProj.cxx:255-281`).
///
/// UNPORTED: the `Geom2d_BSplineCurve` arm (`cxx:268-279`) snaps `param` onto a
/// knot within `Precision::PConfusion()`. `Curve2d` exposes no knot sequence
/// (only `bspline_degree` / `bspline_poles2d`), so `param` is returned
/// unchanged; the two unwrap arms (`cxx:258-266`) alone change nothing.
fn correct_parameter(_c2d: &Arc<dyn Curve2d>, param: f64) -> f64 {
    param
}

/// `Precision::IsInfinite(p.X()) || IsInfinite(p.Y()) || IsInfinite(p.Z())`
/// (`Proj.cxx:319-321`, `:325-327`, `:342-343`, `:348-351`).
fn point_is_infinite(p: &GpPnt) -> bool {
    Precision::is_infinite(p.x()) || Precision::is_infinite(p.y()) || Precision::is_infinite(p.z())
}

/// `ShapeAnalysis_TransferParameters` (`ShapeAnalysis_TransferParameters.cxx`).
#[derive(Clone)]
pub struct TransferParameters {
    /// `myFirst` / `myLast`: the 3D range (`hxx:85-86`).
    pub(crate) first: f64,
    pub(crate) last: f64,
    /// `myEdge` (`hxx:87`).
    pub(crate) edge: Edge,
    /// `myMaxTolerance` (`hxx:88`).
    pub(crate) max_tolerance: f64,
    /// `myFirst2d` / `myLast2d` (`hxx:92-93`).
    pub(crate) first2d: f64,
    pub(crate) last2d: f64,
    /// `myFace` (`hxx:94`).
    face: Option<Face>,
    /// `myShift` - private in OCCT (`hxx:91`).
    shift: f64,
    /// `myScale` - private in OCCT (`hxx:90`).
    scale: f64,
}

impl Default for TransferParameters {
    fn default() -> Self {
        Self::new()
    }
}

impl TransferParameters {
    /// `ShapeAnalysis_TransferParameters()` (`cxx:29-33`): `myShift = 0`,
    /// `myScale = 1`. OCCT leaves `myFirst` / `myLast` / `myFirst2d` /
    /// `myLast2d` / `myMaxTolerance` uninitialised here (they are only set by
    /// `Init`, `cxx:45-74`, or `SetMaxTolerance`, `cxx:79-82`); the port uses
    /// the `Proj` class's own fallbacks (`Proj.cxx:52` sets `myMaxTolerance = 1`)
    /// plus `0` / `1` for the ranges, which every call site overwrites.
    pub fn new() -> Self {
        Self {
            first: 0.0,
            last: 1.0,
            edge: Edge::new(),
            max_tolerance: 1.0,
            first2d: 0.0,
            last2d: 0.0,
            face: None,
            shift: 0.0,
            scale: 1.0,
        }
    }

    /// `ShapeAnalysis_TransferParameters(E, F)` (`cxx:37-41`).
    pub fn with_edge_face(e: &Edge, f: Option<&Face>) -> Self {
        let mut t = Self::new();
        t.init(e, f);
        t
    }

    /// `Init` (`cxx:45-74`).
    ///
    /// OCCT leaves `myFirst` / `myLast` uninitialised when the edge carries no
    /// 3D curve (`cxx:48` declares `double l, f;` and `cxx:55` only writes them
    /// on success) and `myFirst2d` / `myLast2d` at `0` when `F.IsNull()`
    /// (`cxx:49`); the port keeps the constructor's `0` / `1` / `0` / `0` in
    /// those cases.
    pub fn init(&mut self, e: &Edge, f: Option<&Face>) {
        self.scale = 1.0; // `cxx:46-47`
        self.shift = 0.0;
        let mut first2d = 0.0;
        let mut last2d = 0.0;
        self.edge = e.clone(); // `cxx:52`

        // `sae.Curve3d(E, curve3d, f, l, false)` (`cxx:54-56`): `orient = false`,
        // so the range is the stored one, not swapped for a REVERSED edge
        // (`ShapeAnalysis_Edge.cxx:100-125`).
        let mut curve3d: Option<Arc<dyn Curve>> = None;
        if let Some(curve) = BRepTool::edge_curve(e) {
            let (a, b) = BRepTool::edge_parameters(e);
            self.first = a;
            self.last = b;
            curve3d = Some(curve);
        }

        // `sae.PCurve(E, F, curve2d, f2d, l2d, false)` (`cxx:58-62`), skipped for
        // a null face (`cxx:59`: "process free edges").
        let mut curve2d: Option<Arc<dyn Curve2d>> = None;
        if let Some(f) = f {
            if let Some((pc, a, b)) = curve_on_surface_oriented(e, f, false) {
                curve2d = Some(pc);
                first2d = a;
                last2d = b;
            }
        }
        self.first2d = first2d;
        self.last2d = last2d;
        self.face = f.cloned();

        if curve3d.is_none() || curve2d.is_none() {
            return; // `cxx:66-69`
        }
        let ln2d = last2d - first2d; // `cxx:71`
        let ln3d = self.last - self.first; // `cxx:72`
        self.scale = if ln3d <= REAL_SMALL {
            1.0 // `cxx:73`: `ln3d <= gp::Resolution()`
        } else {
            ln2d / ln3d
        };
        self.shift = first2d - self.first * self.scale; // `cxx:74`
    }

    /// `SetMaxTolerance` (`cxx:79-82`).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.max_tolerance = maxtol;
    }

    /// `Perform(Params, To2d)` (`cxx:86-96`).
    pub fn perform(&self, params: &[f64], to2d: bool) -> Vec<f64> {
        params.iter().map(|&p| self.perform_param(p, to2d)).collect()
    }

    /// `Perform(Param, To2d)` (`cxx:100-112`).
    pub fn perform_param(&self, param: f64, to2d: bool) -> f64 {
        if to2d {
            self.shift + param * self.scale // `cxx:105`
        } else {
            -self.shift / self.scale + param / self.scale // `cxx:109`
        }
    }

    /// `TransferRange(newEdge, prevPar, currPar, Is2d)` (`cxx:116-146`).
    pub fn transfer_range(&self, new_edge: &mut Edge, prev_par: f64, curr_par: f64, is2d: bool) {
        if is2d {
            let span2d = self.last2d - self.first2d; // `cxx:124`
            let (tmp1, tmp2) = if prev_par > curr_par {
                (curr_par, prev_par)
            } else {
                (prev_par, curr_par)
            }; // `cxx:126-132`
            let alpha = (tmp1 - self.first2d) / span2d; // `cxx:133`
            let beta = (tmp2 - self.first2d) / span2d; // `cxx:134`
            copy_ranges(new_edge, &self.edge, alpha, beta); // `cxx:135`
        } else {
            let alpha = (prev_par - self.first) / (self.last - self.first); // `cxx:142`
            let beta = (curr_par - self.first) / (self.last - self.first); // `cxx:143`
            copy_ranges(new_edge, &self.edge, alpha, beta); // `cxx:144`
        }
    }

    /// `IsSameRange` (`cxx:150-153`).
    pub fn is_same_range(&self) -> bool {
        self.shift == 0.0 && self.scale == 1.0
    }
}

/// `ShapeAnalysis_TransferParametersProj`
/// (`ShapeAnalysis_TransferParametersProj.cxx`).
#[derive(Clone)]
pub struct TransferParametersProj {
    base: TransferParameters,
    /// `myCurve` (`hxx:97`).
    curve: Option<Arc<dyn Curve>>,
    /// `myCurve2d` (`hxx:98`).
    curve2d: Option<Arc<dyn Curve2d>>,
    /// The surface behind `myAC3d` (`Proj.cxx:97`).
    surface: Option<Arc<dyn Surface>>,
    /// `myAC3d`'s window, `(f2d, l2d)` (`Proj.cxx:94`, `:99`).
    ac_first: f64,
    ac_last: f64,
    /// `myPrecision` (`hxx:102`) - `BRep_Tool::Tolerance(E)` (`Proj.cxx:74`).
    precision: f64,
    /// `myForceProj` (`hxx:104`).
    force_proj: bool,
    /// `myInitOK` (`hxx:105`).
    init_ok: bool,
}

impl Default for TransferParametersProj {
    fn default() -> Self {
        Self::new()
    }
}

impl TransferParametersProj {
    /// `ShapeAnalysis_TransferParametersProj()` (`Proj.cxx:49-54`):
    /// `myPrecision = 0`, `myMaxTolerance = 1`, `myForceProj = false`,
    /// `myInitOK = false`.
    pub fn new() -> Self {
        Self {
            base: TransferParameters::new(),
            curve: None,
            curve2d: None,
            surface: None,
            ac_first: 0.0,
            ac_last: 0.0,
            precision: 0.0,
            force_proj: false,
            init_ok: false,
        }
    }

    /// `ShapeAnalysis_TransferParametersProj(E, F)` (`Proj.cxx:59-66`).
    pub fn with_edge_face(e: &Edge, f: &Face) -> Self {
        let mut t = Self::new();
        t.init(e, f);
        t
    }

    /// `Init` (`Proj.cxx:69-103`).
    pub fn init(&mut self, e: &Edge, f: &Face) {
        self.init_ok = false; // `Proj.cxx:71`
        self.base.init(e, Some(f)); // `Proj.cxx:72`
        // `myEdge = E` is already done by the base `Init` (`cxx:52`).
        self.precision = BRepTool::edge_tolerance(e); // `Proj.cxx:74`

        // `myCurve = BRep_Tool::Curve(E, myFirst, myLast)` (`Proj.cxx:77`); the
        // base already read the same range at `cxx:54-56`.
        self.curve = BRepTool::edge_curve(e);
        if self.curve.is_none() {
            self.base.first = 0.0; // `Proj.cxx:79-82`
            self.base.last = 1.0;
            return;
        }
        if self.base.face.is_none() {
            return; // `Proj.cxx:84-87` (`F.IsNull()`)
        }
        // `sae.PCurve(E, F, myCurve2d, f2d, l2d, false)` (`Proj.cxx:91-92`).
        if let Some((c2d, f2d, l2d)) = curve_on_surface_oriented(e, f, false) {
            self.surface = BRepTool::face_surface(f); // `Proj.cxx:97`
            self.curve2d = Some(c2d);
            self.ac_first = f2d;
            self.ac_last = l2d;
            self.init_ok = true; // `Proj.cxx:102`
        }
    }

    /// `SetMaxTolerance` (`cxx:79-82`).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.base.set_max_tolerance(maxtol);
    }

    /// `ForceProjection` (`Proj.cxx:555-558`).
    pub fn force_projection(&mut self) -> &mut bool {
        &mut self.force_proj
    }

    /// The `!myInitOK || (!myForceProj && myPrecision < myMaxTolerance &&
    /// BRep_Tool::SameParameter(myEdge))` gate repeated at `Proj.cxx:112-116`,
    /// `:186-189`, `:226-228`, `:292-294`, `:542-545`.
    fn use_linear(&self) -> bool {
        !self.init_ok
            || (!self.force_proj
                && self.precision < self.base.max_tolerance
                && BRepTool::same_parameter(&self.base.edge))
    }

    /// `myAC3d` (`Proj.cxx:95-101`): `surface(pcurve(t))` over `[f2d, l2d]`,
    /// `f2d` / `l2d` coming from the COS representation.
    fn ac3d(&self) -> Option<CurveOnSurface> {
        Some(CurveOnSurface::new(
            self.curve2d.clone()?,
            self.surface.clone()?,
            self.ac_first,
            self.ac_last,
        ))
    }

    /// `Adaptor3d_CurveOnSurface Ad1(AC2d(myCurve2d, First, Last), AdS)`
    /// (`Proj.cxx:200-201`, `:462-464`).
    fn ac3d_range(&self, first: f64, last: f64) -> Option<CurveOnSurface> {
        Some(CurveOnSurface::new(
            self.curve2d.clone()?,
            self.surface.clone()?,
            first,
            last,
        ))
    }

    /// `(FirstParameter(), LastParameter())` of `myAC3d` (`Proj.cxx:122-123`,
    /// `:232`, `:241-242`, `:332-336`).
    fn ac_range(&self) -> (f64, f64) {
        (self.ac_first, self.ac_last)
    }

    /// `myCurve->IsClosed()` (`Proj.cxx:145`). The port's substitute is the
    /// coincident-end test (see `shape_analysis_curve`), so a trimmed arc of a
    /// periodic basis is not treated as closed.
    fn curve_is_closed(&self) -> bool {
        let Some(c) = self.curve.as_ref() else {
            return false;
        };
        let (f, l) = (c.first_parameter(), c.last_parameter());
        f.is_finite() && l.is_finite() && c.d0(f).distance(&c.d0(l)) <= CONFUSION
    }

    /// `Perform(Knots, To2d)` (`Proj.cxx:107-176`).
    pub fn perform(&self, knots: &[f64], to2d: bool) -> Vec<f64> {
        if self.use_linear() {
            return self.base.perform(knots, to2d); // `Proj.cxx:112-116`
        }
        let len = knots.len();
        let preci = 2.0 * PCONFUSION; // `Proj.cxx:121`
        let (first, last) = if to2d {
            self.ac_range()
        } else {
            (self.base.first, self.base.last)
        }; // `Proj.cxx:122-123`
        let mut max_par = first; // `Proj.cxx:124`
        let last_par = last;
        let mut prev_par = max_par;
        let mut res_knots = Vec::with_capacity(len);
        for &knot in knots {
            // `Proj.cxx:130-142`
            let par = self.preform_segment(knot, to2d, prev_par, last_par);
            prev_par = par;
            if prev_par > last_par {
                prev_par -= preci; // `Proj.cxx:135-138`
            }
            res_knots.push(par);
            if par > max_par {
                max_par = par;
            }
        }
        // `Proj.cxx:145-160`: pdn correcting on periodic.
        if self.curve_is_closed() {
            for i in (0..len).rev() {
                if res_knots[i] < max_par {
                    let base_last = if to2d {
                        self.ac_last
                    } else {
                        self.curve.as_ref().map(|c| c.last_parameter()).unwrap_or(last)
                    };
                    res_knots[i] = base_last - (len - i - 1) as f64 * preci;
                } else {
                    break;
                }
            }
        }
        // `Proj.cxx:162-174`: pdn correction on range.
        for k in res_knots.iter_mut() {
            *k = k.max(first).min(last);
        }
        res_knots
    }

    /// `Perform(Knot, To2d)` (`Proj.cxx:222-252`).
    pub fn perform_param(&self, knot: f64, to2d: bool) -> f64 {
        if self.use_linear() {
            return self.base.perform_param(knot, to2d); // `Proj.cxx:226-228`
        }
        let (first, last) = if to2d {
            self.ac_range() // `Proj.cxx:232`
        } else {
            (self.base.first, self.base.last) // `Proj.cxx:236`
        };
        let res = self.preform_segment(knot, to2d, first, last);
        // `Proj.cxx:240-251`: pdn correction on range.
        res.max(first).min(last)
    }

    /// `PreformSegment` (`Proj.cxx:180-217`).
    fn preform_segment(&self, param: f64, to2d: bool, first: f64, last: f64) -> f64 {
        let lin_par = self.base.perform_param(param, to2d); // `Proj.cxx:185`
        if self.use_linear() {
            return lin_par; // `Proj.cxx:186-190`
        }
        let (proj_param, lin_dev, proj_dev);
        if to2d {
            // `Proj.cxx:198-203`.
            let (Some(curve), Some(ad1)) = (self.curve.as_ref(), self.ac3d_range(first, last))
            else {
                return lin_par;
            };
            let p1 = curve.d0(param); // `myCurve->Value(Param)`
            let projected = project_adaptor(&ad1, &p1, self.precision, true);
            proj_param = projected.param;
            proj_dev = projected.distance;
            lin_dev = p1.distance(&ad1.d0(lin_par));
        } else {
            // `Proj.cxx:205-210`.
            let (Some(curve), Some(ad1)) =
                (self.curve.as_ref(), self.ac3d_range(self.ac_first, self.ac_last))
            else {
                return lin_par;
            };
            let p1 = ad1.d0(param); // `myAC3d.Value(Param)`
            let projected = project_range(curve.as_ref(), &p1, self.precision, first, last, false);
            proj_param = projected.param;
            proj_dev = projected.distance;
            lin_dev = p1.distance(&curve.d0(lin_par));
        }
        // `Proj.cxx:212-215`.
        if lin_dev <= proj_dev || (lin_dev < self.precision && lin_dev <= 2.0 * proj_dev) {
            return lin_par;
        }
        proj_param
    }

    /// `TransferRange(newEdge, prevPar, currPar, Is2d)` (`Proj.cxx:285-536`).
    pub fn transfer_range(&self, new_edge: &mut Edge, prev_par: f64, curr_par: f64, is2d: bool) {
        if self.use_linear() {
            // `Proj.cxx:291-295`
            self.base.transfer_range(new_edge, prev_par, curr_par, is2d);
            return;
        }
        let reg = GeometryRegistry::global();
        let mut samerange = true; // `Proj.cxx:297`
        copy_ranges(new_edge, &self.base.edge, 0.0, 1.0); // `Proj.cxx:300`
        let preci = PCONFUSION; // `Proj.cxx:305`
        let (first_par, last_par) = if prev_par < curr_par {
            (prev_par, curr_par)
        } else {
            (curr_par, prev_par)
        }; // `Proj.cxx:307-314`

        let (p1, p2, alpha, beta);
        if is2d {
            // `Proj.cxx:317-337`.
            let Some(ad) = self.ac3d() else {
                return;
            };
            p1 = ad.d0(first_par); // `myAC3d.Value(firstPar)`
            if point_is_infinite(&p1) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            p2 = ad.d0(last_par);
            if point_is_infinite(&p2) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            let fact = self.ac_last - self.ac_first; // `Proj.cxx:332`
            let (mut a, mut b) = (0.0, 1.0);
            if fact > epsilon(self.ac_last) {
                a = (first_par - self.ac_first) / fact; // `Proj.cxx:335`
                b = (last_par - self.ac_first) / fact; // `Proj.cxx:336`
            }
            alpha = a;
            beta = b;
        } else {
            // `Proj.cxx:339-361`.
            let Some(curve) = self.curve.as_ref() else {
                return;
            };
            p1 = curve.d0(first_par); // `myCurve->Value(firstPar)`
            if point_is_infinite(&p1) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            p2 = curve.d0(last_par);
            if point_is_infinite(&p2) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            let fact = self.base.last - self.base.first; // `Proj.cxx:356`
            let (mut a, mut b) = (0.0, 1.0);
            if fact > epsilon(self.base.last) {
                a = (first_par - self.base.first) / fact; // `Proj.cxx:359`
                b = (last_par - self.base.first) / fact; // `Proj.cxx:360`
            }
            alpha = a;
            beta = b;
        }
        let use_linear_first = alpha < preci; // `Proj.cxx:362`
        let use_linear_last = 1.0 - beta < preci; // `Proj.cxx:363`

        // `Proj.cxx:367-376`: walk every representation of `newEdge`.
        let Some(mut to_geom) = reg.edge_geom(&new_edge.0) else {
            return;
        };

        // The `BRep_Curve3D` representation (`Proj.cxx:381-451`).
        {
            let (mut ppar1, mut ppar2) = if !is2d {
                (first_par, last_par) // `Proj.cxx:383-387`
            } else {
                let c3d = to_geom.curve.clone(); // `toGC->Curve3D()` (`Proj.cxx:390`)
                let first = to_geom.first; // `Proj.cxx:393-395`
                let last = to_geom.last;
                let len = last - first;
                let lin_first = first + alpha * len; // `Proj.cxx:405`
                let lin_last = first + beta * len; // `Proj.cxx:406`
                let gac = CurveRange::new(c3d.clone(), first, last); // `Proj.cxx:400`
                let proj1 = next_project(&gac, lin_first, &p1, self.precision); // `Proj.cxx:407`
                let proj2 = next_project(&gac, lin_last, &p2, self.precision); // `Proj.cxx:408`
                let use_linear = (proj1.param - proj2.param).abs() < preci; // `Proj.cxx:404`
                let pos1 = c3d.d0(lin_first); // `Proj.cxx:409-410`
                let pos2 = c3d.d0(lin_last);
                let d01 = pos1.distance(&p1); // `Proj.cxx:411-412`
                let d02 = pos2.distance(&p2);
                let (mut a, mut b) = (proj1.param, proj2.param);
                // `Proj.cxx:413-421`
                if use_linear_first
                    || use_linear
                    || d01 <= proj1.distance
                    || (d01 < self.precision && d01 <= 2.0 * proj1.distance)
                {
                    a = lin_first;
                }
                if use_linear_last
                    || use_linear
                    || d02 <= proj2.distance
                    || (d02 < self.precision && d02 <= 2.0 * proj2.distance)
                {
                    b = lin_last;
                }
                (a, b)
            };
            // `Proj.cxx:423-443`: order and nudge.
            if ppar1 > ppar2 {
                std::mem::swap(&mut ppar1, &mut ppar2);
            }
            if ppar2 - ppar1 < preci {
                if ppar1 - to_geom.first < preci {
                    ppar2 += 2.0 * preci;
                } else if to_geom.last - ppar2 < preci {
                    ppar1 -= 2.0 * preci;
                } else {
                    ppar1 -= preci;
                    ppar2 += preci;
                }
            }
            to_geom.first = ppar1; // `Proj.cxx:444` (`toGC->SetRange`)
            to_geom.last = ppar2;
            if ppar1 != first_par || ppar2 != last_par {
                samerange = false; // `Proj.cxx:447-450`
            }
        }

        // The `BRep_CurveOnSurface` representations (`Proj.cxx:452-531`).
        let keys: Vec<usize> = to_geom.pcurves.keys().copied().collect();
        for key in keys {
            let Some(c2d) = to_geom.pcurves.get(&key).and_then(|v| v.first().cloned()) else {
                continue;
            };
            let Some(surface) = self.surface.clone() else {
                continue;
            };
            let (first, last) = to_geom
                .pcurve_ranges
                .get(&key)
                .copied()
                .unwrap_or((to_geom.first, to_geom.last)); // `Proj.cxx:457-459`
            let len = last - first;
            let ad1 = CurveOnSurface::new(c2d.clone(), surface, first, last); // `Proj.cxx:462-464`
            let lin_first = first + alpha * len; // `Proj.cxx:471`
            let lin_last = first + beta * len; // `Proj.cxx:472`
            let proj1 = next_project(&ad1, lin_first, &p1, self.precision); // `Proj.cxx:473`
            let proj2 = next_project(&ad1, lin_last, &p2, self.precision); // `Proj.cxx:474`
            // `Proj.cxx:476-478`
            let is_first_on_end = (proj1.param - first) / len < PCONFUSION;
            let is_last_on_end = (last - proj2.param) / len < PCONFUSION;
            let use_linear = (proj1.param - proj2.param).abs() < PCONFUSION;
            let mut local_linear_first = use_linear_first;
            let mut local_linear_last = use_linear_last;
            if is_first_on_end && !local_linear_first {
                local_linear_first = true; // `Proj.cxx:479-482`
            }
            if is_last_on_end && !local_linear_last {
                local_linear_last = true; // `Proj.cxx:483-486`
            }
            let pos1 = ad1.d0(lin_first); // `Proj.cxx:488-489`
            let pos2 = ad1.d0(lin_last);
            let d01 = pos1.distance(&p1); // `Proj.cxx:490-491`
            let d02 = pos2.distance(&p2);
            let mut ppar1 = proj1.param;
            let mut ppar2 = proj2.param;
            // `Proj.cxx:492-499`
            if local_linear_first
                || use_linear
                || d01 <= proj1.distance
                || (d01 < self.precision && d01 <= 2.0 * proj1.distance)
            {
                ppar1 = lin_first;
            }
            if local_linear_last
                || use_linear
                || d02 <= proj2.distance
                || (d02 < self.precision && d02 <= 2.0 * proj2.distance)
            {
                ppar2 = lin_last;
            }
            if ppar1 > ppar2 {
                std::mem::swap(&mut ppar1, &mut ppar2); // `Proj.cxx:501-508`
            }
            ppar1 = correct_parameter(&c2d, ppar1); // `Proj.cxx:509-510`
            ppar2 = correct_parameter(&c2d, ppar2);
            if ppar2 - ppar1 < preci {
                // `Proj.cxx:511-524`
                if ppar1 - first < preci {
                    ppar2 += 2.0 * preci;
                } else if last - ppar2 < preci {
                    ppar1 -= 2.0 * preci;
                } else {
                    ppar1 -= preci;
                    ppar2 += preci;
                }
            }
            to_geom.pcurve_ranges.insert(key, (ppar1, ppar2)); // `Proj.cxx:527`
            if ppar1 != first_par || ppar2 != last_par {
                samerange = false; // `Proj.cxx:528-531`
            }
        }
        reg.set_edge(&new_edge.0, to_geom);
        reg.set_same_range(&new_edge.0, samerange); // `Proj.cxx:534`
    }

    /// `IsSameRange` (`Proj.cxx:539-551`).
    pub fn is_same_range(&self) -> bool {
        if self.use_linear() {
            self.base.is_same_range()
        } else {
            false // `Proj.cxx:549`
        }
    }
}
