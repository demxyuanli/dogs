//! `HLRBRep_InterCSurf` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_InterCSurf.hxx:49-208`,
//! `.cxx:59-577`).
//!
//! "Compute the intersection between the curve and the surface"
//! (`hxx:56-57`). The class is nothing but an instantiation of the
//! `IntCurveSurface_InterImpl` template
//! (`HLRBRep_InterCSurf.cxx:39-56`):
//!
//! ```text
//! TheCurve       = gp_Lin                  TheCurveTool   = HLRBRep_LineTool
//! TheSurface     = HLRBRep_Surface*        TheSurfaceTool = HLRBRep_SurfaceTool
//! ThePolygon     = HLRBRep_ThePolygonOfInterCSurf
//! ThePolyhedron  = HLRBRep_ThePolyhedronOfInterCSurf
//! TheInterference/TheCSFunction/TheExactInter/TheQuadCurvExactInter
//! ```
//!
//! so every `HLRBRep_InterCSurf::X` body in the `.cxx` is a call to
//! `IntCurveSurface_InterImpl::X` with those parameters and nothing else
//! (`cxx:106-577`). The port therefore reuses the shared
//! [`occt_topo::intcurvesurface`] engine (the same template already ported for
//! `IntCurveSurface_HInter`) and adds the two pieces that make it an HLR
//! instantiation:
//!
//! * [`HlrBRepLineCurve`] - `HLRBRep_LineTool` seen as an
//!   [`occt_geom::Curve`] (`HLRBRep_LineTool.hxx:47-179`),
//! * [`HlrBRepSurfaceView`] - `HLRBRep_SurfaceTool` seen as an
//!   [`occt_geom::Surface`] (`hlrbrep_surface/adaptor.rs`).
//!
//! The sampling counts (`HLRBRep_SurfaceTool::NbSamplesU/V`, `.cxx:22-155`)
//! are the HLR ones, not `Adaptor3d_HSurfaceTool`'s, so the polyhedron grid
//! of [perform_polygon] / [perform_polygon_polyhedron] goes through
//! [`crate::hlrbrep_surface_tool`].
//!
//! ## Ported
//!
//! * `Perform(theCurve)` / `Perform(theCurve, theSurface)`
//!   (`cxx:106-113` = `InterImpl::Perform`, `Inter.pxx:63-87`),
//! * `Perform(theCurve, theSurface, U1, V1, U2, V2)` (`cxx:122-152` =
//!   `PerformBounds`, `Inter.pxx:105-182`),
//! * `Perform(theCurve, thePolygon, theSurface)` (`cxx:156-168` =
//!   `PerformPolygon`, `Inter.pxx:193-217`),
//! * `Perform(theCurve, theSurface, thePolyhedron)` (`cxx:174-191` =
//!   `PerformPolyhedron`, `Inter.pxx:229-242`),
//! * `Perform(theCurve, thePolygon, theSurface, thePolyhedron)`
//!   (`cxx:195-220` = `PerformPolygonPolyhedron`, `Inter.pxx:253-268`),
//! * `Perform(..., theBndBSB)` (`cxx:224-253` = `PerformPolygonPolyhedronBSB`,
//!   `Inter.pxx:279-295`),
//! * `InternalPerform` both overloads (`cxx:256-317` = `InternalPerformBSB`
//!   `Inter.pxx:308-355` and `InternalPerform` `Inter.pxx:367-417`),
//! * `PerformConicSurf(gp_Lin, ...)` (`cxx:369-405` =
//!   `PerformConicSurfLine`, `Inter.pxx:523-708`),
//! * `AppendIntAna` (`cxx:540-550` = `InterImpl::AppendIntAna`,
//!   `Inter.pxx:968-1003` + `ProcessIntAna`, `InterUtils.pxx:1188-1227`),
//! * `AppendSegment` (`cxx:573-577`): OCCT's body is `// Not implemented`,
//!   so the port keeps an empty [`HlrBRepInterCSurf::append_segment`] and no
//!   caller can add a segment,
//! * `IsDone` / `IsParallel` / `NbPoints` / `Point` / `NbSegments` / `Segment`
//!   and `ResetFields` (`IntCurveSurface_Intersection.cxx:26-133`).
//!
//! ## UNPORTED
//!
//! * `InternalPerformCurveQuadric` (`cxx:321-330`) and
//!   `InternalPerform(theCurve, thePolygon, theSurface, U1, V1, U2, V2)`
//!   (`cxx:334-366` = `InternalPerformPolygonBounds`, `Inter.pxx:451-520`).
//!   Both are only reached from the `gp_Circ` / `gp_Elips` / `gp_Parab` /
//!   `gp_Hypr` arms of `PerformBounds` (`Inter.pxx:126-136`), and
//!   `HLRBRep_LineTool::GetType` answers `GeomAbs_Line` unconditionally
//!   (`HLRBRep_LineTool.lxx:127-130`), so in this instantiation they are dead
//!   code. Their B-spline arm (`Inter.pxx:487-510`) would additionally need
//!   `Adaptor3d_TopolTool::SamplePnts`, which `HLRBRep_SurfaceTool` does not
//!   provide.
//! * `EstLimForInfExtr` / `EstLimForInfRevl` / `EstLimForInfOffs`
//!   (`IntCurveSurface_InterUtils.pxx:148-700`) inside the non-analytic arm of
//!   `PerformConicSurfLine` (`Inter.pxx:625-669`). They need
//!   `ProjLib_Plane` + `IntAna2d_AnaIntersection` and, for the extrusion /
//!   revolution / offset arms, `DoSurface` / `DoNewBounds`
//!   (`InterUtils.pxx:930-1114`) with `HLRBRep_SurfaceTool::BasisSurface`
//!   (which `HLRBRep_SurfaceTool.hxx` does not declare). The branch is
//!   unreachable for an HLR surface: `Precision::IsInfinite` is true only past
//!   `1e100` and every `HLRBRep_SurfaceTool::First/LastU/VParameter` is a
//!   `BRepAdaptor_Surface` face bound (`HLRBRep_Surface.lxx:37-61`), i.e. a
//!   finite number. The port keeps the branch, ports the generic
//!   `EstLimForInfSurf` fallback (`InterUtils.pxx:137-143`) and returns early
//!   with this note for the three special cases, so nothing is invented.
//! * `DoSurface` / `DoNewBounds` (`cxx:62-102`): only the unported estimators
//!   call them.

use occt_core::bnd::BndBox;
use occt_core::elib::{clib, slib};
use occt_core::gp::{GpLin, GpPnt, GpTrsf, GpVec};
use occt_core::intf::{IntfPolyhGrid, IntfTool};
use occt_core::precision;
use occt_geom::curve::Curve;
use occt_geom::intana::{line_torus_intersect, IntAnaIntConicQuad, IntAnaQuadric};
use occt_geom::surface::Surface;
use occt_math::MathFunctionSetRoot;

use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::hlrbrep_line_tool;
use crate::hlrbrep_surface::{HlrBRepSurface, HlrBRepSurfaceView};
use crate::hlrbrep_surface_tool;
use occt_topo::intcurvesurface::{
    clamp_uv_parameters, collect_interference_points, compute_append_point,
    decompose_surface_intervals, process_sorted_points, sort_start_points, surface_params,
    IntersectionPoint, IntersectionSegment, SortedStartPoints, TheCSFunctionOfHInter,
    TheExactHInter, TheInterferenceOfHInter, ThePolygonOfHInter, ThePolyhedronOfHInter,
    THE_TOLERANCE_ANGULAIRE, THE_TOLTANGENCY,
};

/// `HLRBRep_LineTool` (`HLRBRep_LineTool.hxx:47-179`) as an
/// [`occt_geom::Curve`], i.e. the `TheCurveTool` half of the instantiation
/// (`HLRBRep_InterCSurf.cxx:44-45`).
///
/// Every method forwards to [`crate::hlrbrep_line_tool`], so the answers are the
/// `.lxx` ones: `FirstParameter` = `RealFirst()`, `LastParameter` =
/// `RealLast()` (`lxx:29-35`) - not `+/-inf` as `Geom_Line` answers - and
/// `Period()` is `0.` (`lxx:80-83`) although the range is `2 * RealLast()`.
#[derive(Debug, Clone)]
pub struct HlrBRepLineCurve {
    lin: GpLin,
}

impl HlrBRepLineCurve {
    /// The `gp_Lin` the tool is instantiated on.
    pub fn new(lin: GpLin) -> Self {
        Self { lin }
    }

    /// `HLRBRep_LineTool::Line(C)` (`lxx:132-135`): the curve itself.
    pub fn lin(&self) -> &GpLin {
        &self.lin
    }
}

impl Curve for HlrBRepLineCurve {
    fn d0(&self, u: f64) -> GpPnt {
        hlrbrep_line_tool::d0(&self.lin, u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        hlrbrep_line_tool::d1(&self.lin, u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        hlrbrep_line_tool::d2(&self.lin, u)
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        hlrbrep_line_tool::d3(&self.lin, u)
    }

    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        hlrbrep_line_tool::dn(&self.lin, u, n)
    }

    fn first_parameter(&self) -> f64 {
        hlrbrep_line_tool::first_parameter(&self.lin)
    }

    fn last_parameter(&self) -> f64 {
        hlrbrep_line_tool::last_parameter(&self.lin)
    }

    fn is_periodic(&self) -> bool {
        hlrbrep_line_tool::is_periodic(&self.lin)
    }

    fn period(&self) -> f64 {
        hlrbrep_line_tool::period(&self.lin)
    }

    fn is_closed(&self) -> bool {
        hlrbrep_line_tool::is_closed(&self.lin)
    }

    fn continuity(&self) -> u8 {
        hlrbrep_line_tool::continuity(&self.lin)
    }

    fn nb_intervals(&self, s: u8) -> i32 {
        hlrbrep_line_tool::nb_intervals(&self.lin, s) as i32
    }

    fn resolution(&self, r3d: f64) -> f64 {
        hlrbrep_line_tool::resolution(&self.lin, r3d)
    }

    fn is_line(&self) -> bool {
        hlrbrep_line_tool::get_type(&self.lin) == hlrbrep_line_tool::CurveType::Line
    }

    fn gp_line(&self) -> Option<GpLin> {
        Some(self.lin)
    }

    fn transform(&mut self, t: &GpTrsf) {
        self.lin.transform(t);
    }

    fn reverse(&mut self) {
        self.lin.pos.vdir.reverse();
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}

/// `HLRBRep_InterCSurf` (`HLRBRep_InterCSurf.hxx:49-208`).
///
/// The state is the `IntCurveSurface_Intersection` base
/// (`IntCurveSurface_Intersection.hxx:33-136`): `done`, `myIsParallel` and the
/// point / segment sequences.
#[derive(Debug, Clone, Default)]
pub struct HlrBRepInterCSurf {
    /// `done` (`IntCurveSurface_Intersection.hxx:127`).
    done: bool,
    /// `myIsParallel` (`hxx:128`): "Curve is parallel surface. This case is
    /// recognized only for some pairs of analytical curves and surfaces
    /// (plane - line, ...)" (`hxx:128-131`).
    is_parallel: bool,
    /// `lpnt` (`hxx:134`).
    points: Vec<IntersectionPoint>,
    /// `lseg` (`hxx:135`).
    segments: Vec<IntersectionSegment>,
}

impl HlrBRepInterCSurf {
    /// `HLRBRep_InterCSurf()` (`cxx:59` = `IntCurveSurface_Intersection()`,
    /// `IntCurveSurface_Intersection.cxx:26-30`): `done = false`,
    /// `myIsParallel = false`, empty sequences.
    pub fn new() -> Self {
        Self::default()
    }

    /// `ResetFields()` (`IntCurveSurface_Intersection.cxx:119-128`): clear the
    /// sequences and both flags, but only when the last run finished.
    fn reset_fields(&mut self) {
        if self.done {
            self.segments.clear();
            self.points.clear();
            self.done = false;
            self.is_parallel = false;
        }
    }

    /// `Append(const IntCurveSurface_IntersectionPoint&)`
    /// (`IntCurveSurface_Intersection.cxx:88-113`): skip a point whose
    /// `(u, v, w)` and transition all match an existing one (`PARAMEQUAL` is
    /// `1e-8`).
    fn append(&mut self, pt: IntersectionPoint) {
        /// `#define PARAMEQUAL(a, b) (std::abs((a) - (b)) < (1.e-8))`
        /// (`IntCurveSurface_Intersection.cxx:21`).
        const PARAMEQUAL: f64 = 1.0e-8;
        for q in &self.points {
            if (pt.u - q.u).abs() < PARAMEQUAL
                && (pt.v - q.v).abs() < PARAMEQUAL
                && (pt.param - q.param).abs() < PARAMEQUAL
                && q.state == pt.state
            {
                return;
            }
        }
        self.points.push(pt);
    }

    /// `Append(const IntCurveSurface_IntersectionSegment&)`
    /// (`IntCurveSurface_Intersection.cxx:116-119`): a plain push, no merge.
    ///
    /// `HLRBRep_InterCSurf::AppendSegment` (`cxx:573-577`) is `// Not
    /// implemented`, and no conic arm of `PerformConicSurfLine` appends a
    /// segment, so nothing in this class reaches this method; it stays for the
    /// base-class contract.
    #[allow(dead_code)]
    fn append_segment(&mut self, seg: IntersectionSegment) {
        self.segments.push(seg);
    }

    /// `IsDone()` (`IntCurveSurface_Intersection.cxx:32-35`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`IntCurveSurface_Intersection.cxx:37-40`).
    pub fn is_parallel(&self) -> bool {
        self.is_parallel
    }

    /// `NbPoints()` (`IntCurveSurface_Intersection.cxx:42-49`): raises
    /// `StdFail_NotDone` when `!IsDone()`, a panic here.
    pub fn nb_points(&self) -> usize {
        assert!(self.done, "IntCurveSurface_Intersection::NbPoints: NotDone");
        self.points.len()
    }

    /// `Point(N)` (`IntCurveSurface_Intersection.cxx:51-58`). OCCT is 1-based;
    /// the port is 0-based, as every other sequence accessor in the crate.
    pub fn point(&self, n: usize) -> &IntersectionPoint {
        assert!(self.done, "IntCurveSurface_Intersection::Point: NotDone");
        &self.points[n]
    }

    /// `NbSegments()` (`IntCurveSurface_Intersection.cxx:60-67`).
    pub fn nb_segments(&self) -> usize {
        assert!(
            self.done,
            "IntCurveSurface_Intersection::NbSegments: NotDone"
        );
        self.segments.len()
    }

    /// `Segment(N)` (`IntCurveSurface_Intersection.cxx:69-76`), 0-based.
    pub fn segment(&self, n: usize) -> &IntersectionSegment {
        assert!(
            self.done,
            "IntCurveSurface_Intersection::Segment: NotDone"
        );
        &self.segments[n]
    }

    /// `Perform(theCurve, theSurface)` (`cxx:106-113` = `InterImpl::Perform`,
    /// `Inter.pxx:63-87`): reset, then `PerformBounds` once per C2 interval of
    /// the surface.
    ///
    /// `None` stands for the state where `HLRBRep_SurfaceTool::D0` would throw
    /// `Standard_NullObject` because the HLR surface carries no underlying
    /// `Geom_Surface` (`BRepAdaptor_Surface::Surface` on an uninitialized
    /// adaptor); the port answers "no intersection" instead of throwing.
    pub fn perform(&mut self, curve: &GpLin, surface: &HlrBRepSurface) {
        self.reset_fields();
        self.done = true;

        let Some(view) = HlrBRepSurfaceView::new(surface) else {
            return;
        };
        let mut intervals = Vec::new();
        decompose_surface_intervals(&view, &mut intervals);
        for iv in intervals {
            self.perform_bounds(curve, &view, iv.u0, iv.v0, iv.u1, iv.v1);
        }
    }

    /// `Perform(theCurve, theSurface, U1, V1, U2, V2)` (`cxx:122-152` =
    /// `PerformBounds`, `Inter.pxx:105-182`).
    ///
    /// `HLRBRep_LineTool::GetType` always answers `GeomAbs_Line`, so the only
    /// reachable arm is the `GeomAbs_Line` one (`Inter.pxx:122-124`); the
    /// `default` arm of the switch (`Inter.pxx:138-180`) is dead here for the
    /// reason given in the module notes.
    fn perform_bounds(
        &mut self,
        curve: &GpLin,
        surface: &HlrBRepSurfaceView<'_>,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        // `Protection from double type overflow (bug26525)` (`Inter.pxx:116-118`).
        let (mut uu1, mut uu2, mut vv1, mut vv2) = (u1, u2, v1, v2);
        clamp_uv_parameters(&mut uu1, &mut uu2, &mut vv1, &mut vv2);

        self.perform_conic_surf_line(curve, surface, uu1, vv1, uu2, vv2);
    }

    /// `Perform(theCurve, thePolygon, theSurface)` (`cxx:156-168` =
    /// `PerformPolygon`, `Inter.pxx:193-217`): build the polyhedron from the
    /// HLR surface tool's sample counts and hand both to
    /// `PerformPolygonPolyhedron`.
    pub fn perform_with_polygon(
        &mut self,
        curve: &GpLin,
        polygon: &ThePolygonOfHInter,
        surface: &HlrBRepSurface,
    ) {
        self.reset_fields();
        self.done = true;

        let Some(view) = HlrBRepSurfaceView::new(surface) else {
            return;
        };
        let u1 = hlrbrep_surface_tool::first_u_parameter(surface);
        let v1 = hlrbrep_surface_tool::first_v_parameter(surface);
        let u2 = hlrbrep_surface_tool::last_u_parameter(surface);
        let v2 = hlrbrep_surface_tool::last_v_parameter(surface);
        let mut nbsu = hlrbrep_surface_tool::nb_samples_u_range(surface, u1, u2);
        let mut nbsv = hlrbrep_surface_tool::nb_samples_v_range(surface, v1, v2);
        if nbsu > 40 {
            nbsu = 40;
        }
        if nbsv > 40 {
            nbsv = 40;
        }
        let polyhedron = ThePolyhedronOfHInter::new(
            &view,
            nbsu.max(1) as usize,
            nbsv.max(1) as usize,
            u1,
            v1,
            u2,
            v2,
        );
        self.perform_polygon_polyhedron(curve, polygon, surface, &polyhedron);
    }

    /// `Perform(theCurve, theSurface, thePolyhedron)` (`cxx:174-191` =
    /// `PerformPolyhedron`, `Inter.pxx:229-242`): sample the polygon over the
    /// curve tool's own range (`HLRBRep_LineTool::FirstParameter` =
    /// `RealFirst()`, `.lxx:29-32`) with
    /// `HLRBRep_LineTool::NbSamples` = `3` (`.lxx:201-204`).
    pub fn perform_with_polyhedron(
        &mut self,
        curve: &GpLin,
        surface: &HlrBRepSurface,
        polyhedron: &ThePolyhedronOfHInter,
    ) {
        self.reset_fields();
        self.done = true;

        let u1 = hlrbrep_line_tool::first_parameter(curve);
        let u2 = hlrbrep_line_tool::last_parameter(curve);
        let line = HlrBRepLineCurve::new(*curve);
        let polygon =
            ThePolygonOfHInter::new(&line, hlrbrep_line_tool::nb_samples(curve, u1, u2));
        self.perform_polygon_polyhedron(curve, &polygon, surface, polyhedron);
    }

    /// `Perform(theCurve, thePolygon, theSurface, thePolyhedron)`
    /// (`cxx:195-220` = `PerformPolygonPolyhedron`, `Inter.pxx:253-268`): the
    /// UV window for `InternalPerform` is the surface tool's full range, not
    /// the polyhedron's.
    pub fn perform_polygon_polyhedron(
        &mut self,
        curve: &GpLin,
        polygon: &ThePolygonOfHInter,
        surface: &HlrBRepSurface,
        polyhedron: &ThePolyhedronOfHInter,
    ) {
        self.reset_fields();
        self.done = true;

        let Some(view) = HlrBRepSurfaceView::new(surface) else {
            return;
        };
        let u1 = hlrbrep_surface_tool::first_u_parameter(surface);
        let v1 = hlrbrep_surface_tool::first_v_parameter(surface);
        let u2 = hlrbrep_surface_tool::last_u_parameter(surface);
        let v2 = hlrbrep_surface_tool::last_v_parameter(surface);
        let line = HlrBRepLineCurve::new(*curve);
        self.internal_perform(&line, polygon, &view, polyhedron, u1, v1, u2, v2);
    }

    /// `Perform(theCurve, thePolygon, theSurface, thePolyhedron, theBndBSB)`
    /// (`cxx:224-253` = `PerformPolygonPolyhedronBSB`, `Inter.pxx:279-295`).
    pub fn perform_polygon_polyhedron_bsb(
        &mut self,
        curve: &GpLin,
        polygon: &ThePolygonOfHInter,
        surface: &HlrBRepSurface,
        polyhedron: &ThePolyhedronOfHInter,
        bsb: &IntfPolyhGrid,
    ) {
        self.reset_fields();
        self.done = true;

        let Some(view) = HlrBRepSurfaceView::new(surface) else {
            return;
        };
        let u1 = hlrbrep_surface_tool::first_u_parameter(surface);
        let v1 = hlrbrep_surface_tool::first_v_parameter(surface);
        let u2 = hlrbrep_surface_tool::last_u_parameter(surface);
        let v2 = hlrbrep_surface_tool::last_v_parameter(surface);
        let line = HlrBRepLineCurve::new(*curve);
        self.internal_perform_bsb(&line, polygon, &view, polyhedron, u1, v1, u2, v2, bsb);
    }

    /// `InternalPerform(..., theBSB)` (`cxx:256-286` = `InternalPerformBSB`,
    /// `Inter.pxx:308-355`): the only difference to the BSB-less overload is
    /// the interference constructed from the bounding sort box.
    #[allow(clippy::too_many_arguments)]
    fn internal_perform_bsb(
        &mut self,
        curve: &dyn Curve,
        polygon: &ThePolygonOfHInter,
        surface: &dyn Surface,
        polyhedron: &ThePolyhedronOfHInter,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
        bsb: &IntfPolyhGrid,
    ) {
        let interference =
            TheInterferenceOfHInter::of_polygon_polyhedron_grid(polygon, polyhedron, bsb);
        let func = TheCSFunctionOfHInter::new(surface, curve);
        let mut exact_inter = TheExactHInter::new(func, THE_TOLTANGENCY);
        let mut rsnld = MathFunctionSetRoot::with_iterations(exact_inter.function(), 100);

        let winf = polygon.inf_parameter();
        let wsup = polygon.sup_parameter();

        let mut start_points = SortedStartPoints::new();
        collect_interference_points(&interference, polyhedron, polygon, &mut start_points);
        sort_start_points(&mut start_points);

        let mut pts: Vec<IntersectionPoint> = Vec::new();
        process_sorted_points(
            &mut exact_inter,
            &mut rsnld,
            &start_points,
            u0,
            u1,
            v0,
            v1,
            winf,
            wsup,
            curve,
            surface,
            &mut pts,
        );
        for p in pts {
            self.append(p);
        }
    }

    /// `InternalPerform(...)` (`cxx:289-317` = `InternalPerform`,
    /// `Inter.pxx:367-417`): polygon/polyhedron interference for the start
    /// points, then `IntCurveSurface_TheExactHInter` Newton refinement inside
    /// the `[u0, u1] x [v0, v1]` window.
    #[allow(clippy::too_many_arguments)]
    fn internal_perform(
        &mut self,
        curve: &dyn Curve,
        polygon: &ThePolygonOfHInter,
        surface: &dyn Surface,
        polyhedron: &ThePolyhedronOfHInter,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
    ) {
        let interference = TheInterferenceOfHInter::of_polygon_polyhedron(polygon, polyhedron);
        let func = TheCSFunctionOfHInter::new(surface, curve);
        let mut exact_inter = TheExactHInter::new(func, THE_TOLTANGENCY);
        let mut rsnld = MathFunctionSetRoot::with_iterations(exact_inter.function(), 100);

        let winf = polygon.inf_parameter();
        let wsup = polygon.sup_parameter();

        let mut start_points = SortedStartPoints::new();
        collect_interference_points(&interference, polyhedron, polygon, &mut start_points);
        sort_start_points(&mut start_points);

        let mut pts: Vec<IntersectionPoint> = Vec::new();
        process_sorted_points(
            &mut exact_inter,
            &mut rsnld,
            &start_points,
            u0,
            u1,
            v0,
            v1,
            winf,
            wsup,
            curve,
            surface,
            &mut pts,
        );
        for p in pts {
            self.append(p);
        }
    }

    /// `PerformConicSurf(gp_Lin, ...)` (`cxx:369-405` =
    /// `PerformConicSurfLine`, `Inter.pxx:523-708`).
    ///
    /// The analytic arms hand the line and the surface's own quadric to
    /// `IntAna_IntConicQuad` (`Inter.pxx:547-597`); the line/torus arm runs
    /// `IntAna_IntLinTorus` through `ProcessLinTorus`
    /// (`InterUtils.pxx:1282-1317`) and falls back to the sampled path when the
    /// solver is not done (`Inter.pxx:576-582`). Everything else sets
    /// `isAnaProcessed = false` and takes the polyhedron arm
    /// (`Inter.pxx:601-707`).
    fn perform_conic_surf_line(
        &mut self,
        lin: &GpLin,
        surface: &HlrBRepSurfaceView<'_>,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        let hlr = surface.hlr();
        let line = HlrBRepLineCurve::new(*lin);
        let mut is_ana_processed = true;

        match hlrbrep_surface_tool::get_type(hlr) {
            // `case GeomAbs_Plane:` (`Inter.pxx:547-554`):
            // `IntAna_IntConicQuad(Lin, Pln, THE_TOLERANCE_ANGULAIRE)`.
            GeomAbsSurfaceType::Plane => match hlrbrep_surface_tool::plane(hlr) {
                Some(pln) => {
                    let ana = IntAnaIntConicQuad::line_plane(
                        lin,
                        &pln,
                        THE_TOLERANCE_ANGULAIRE,
                        0.0,
                        0.0,
                    );
                    self.append_int_ana(&line, surface, &ana);
                }
                None => is_ana_processed = false,
            },
            // `case GeomAbs_Cylinder:` (`Inter.pxx:555-558`) - the
            // `IntAna_IntConicQuad(const gp_Lin&, const gp_Cylinder&)` ctor is
            // `myQuadric.SetCylinder(...)` + `Perform(Lin)`
            // (`IntAna_IntConicQuad.cxx:212-219`).
            GeomAbsSurfaceType::Cylinder => match hlrbrep_surface_tool::cylinder(hlr) {
                Some(cyl) => {
                    let ana =
                        IntAnaIntConicQuad::line_quadric(lin, &IntAnaQuadric::from_cylinder(&cyl));
                    self.append_int_ana(&line, surface, &ana);
                }
                None => is_ana_processed = false,
            },
            // `case GeomAbs_Sphere:` (`Inter.pxx:559-563`).
            GeomAbsSurfaceType::Sphere => match hlrbrep_surface_tool::sphere(hlr) {
                Some(sph) => {
                    let ana =
                        IntAnaIntConicQuad::line_quadric(lin, &IntAnaQuadric::from_sphere(&sph));
                    self.append_int_ana(&line, surface, &ana);
                }
                None => is_ana_processed = false,
            },
            // `case GeomAbs_Torus:` (`Inter.pxx:564-583`): `ProcessLinTorus`.
            GeomAbsSurfaceType::Torus => {
                let processed = match hlrbrep_surface_tool::torus(hlr) {
                    Some(tor) => {
                        let inter = line_torus_intersect(lin, &tor);
                        if !inter.is_done() {
                            false
                        } else {
                            for p in inter.points() {
                                let w = clib::line_parameter(lin, p);
                                // `ParamOnTorus(i, Fi, Theta)` (`...lxx:31-46`):
                                // `ElSLib::TorusParameters` on the same torus.
                                let (fi, theta) = slib::torus_parameters(
                                    tor.position(),
                                    tor.major_radius(),
                                    tor.minor_radius(),
                                    p,
                                );
                                if let Some(pt) =
                                    compute_append_point(&line, w, surface, fi, theta)
                                {
                                    self.append(pt);
                                }
                            }
                            true
                        }
                    }
                    None => false,
                };
                if !processed {
                    is_ana_processed = false;
                }
            }
            // `case GeomAbs_Cone:` (`Inter.pxx:584-596`): the analytic solve
            // only applies while the half angle stays off `PI/2`.
            GeomAbsSurfaceType::Cone => {
                let correction = 1.0e5 * precision::ANGULAR;
                match hlrbrep_surface_tool::cone(hlr) {
                    Some(cn) if cn.semi_angle().abs() < std::f64::consts::FRAC_PI_2 - correction => {
                        let ana =
                            IntAnaIntConicQuad::line_quadric(lin, &IntAnaQuadric::from_cone(&cn));
                        self.append_int_ana(&line, surface, &ana);
                    }
                    _ => is_ana_processed = false,
                }
            }
            // `default: isAnaProcessed = false;` (`Inter.pxx:598-599`).
            _ => is_ana_processed = false,
        }

        if !is_ana_processed {
            self.perform_polyhedron_fallback(lin, surface, u1, v1, u2, v2);
        }
    }

    /// The `!isAnaProcessed` arm of `PerformConicSurfLine`
    /// (`Inter.pxx:601-707`): sample the surface on the full `[u1, u2] x
    /// [v1, v2]` window with `HLRBRep_SurfaceTool`'s counts (at least `20` each),
    /// clip the line against that polyhedron's box (`Intf_Tool::LinBox`) and run
    /// `InternalPerform` on every resulting line segment.
    fn perform_polyhedron_fallback(
        &mut self,
        lin: &GpLin,
        surface: &HlrBRepSurfaceView<'_>,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        let hlr = surface.hlr();
        let line = HlrBRepLineCurve::new(*lin);

        let mut nbsu = hlrbrep_surface_tool::nb_samples_u_range(hlr, u1, u2);
        let mut nbsv = hlrbrep_surface_tool::nb_samples_v_range(hlr, v1, v2);
        if nbsu < 1 {
            nbsu = 1;
        }
        if nbsv < 1 {
            nbsv = 1;
        }

        let u1inf = precision::Precision::is_infinite(u1);
        let u2inf = precision::Precision::is_infinite(u2);
        let v1inf = precision::Precision::is_infinite(v1);
        let v2inf = precision::Precision::is_infinite(v2);

        let (mut u1new, mut u2new, mut v1new, mut v2new) = (u1, u2, v1, v2);

        if u1inf || u2inf || v1inf || v2inf {
            match hlrbrep_surface_tool::get_type(hlr) {
                // UNPORTED: `EstLimForInfExtr` (`InterUtils.pxx:148-374`). See
                // the module notes; unreachable because an HLR face's UV bounds
                // are finite (`HLRBRep_Surface.lxx:37-61`).
                GeomAbsSurfaceType::SurfaceOfExtrusion => return,
                // UNPORTED: `EstLimForInfRevl` (`InterUtils.pxx:376-497`).
                GeomAbsSurfaceType::SurfaceOfRevolution => return,
                // UNPORTED: `EstLimForInfOffs` (`InterUtils.pxx:499-700`).
                GeomAbsSurfaceType::OffsetSurface => return,
                // `EstLimForInfSurf` (`InterUtils.pxx:137-143`).
                _ => {
                    u1new = u1new.max(-1.0e10);
                    u2new = u2new.min(1.0e10);
                    v1new = v1new.max(-1.0e10);
                    v2new = v2new.min(1.0e10);
                }
            }
        }
        // `if (NoIntersection) return;` (`Inter.pxx:674-677`): only the three
        // unported estimators above can raise it, so it stays false here.

        if nbsu < 20 {
            nbsu = 20;
        }
        if nbsv < 20 {
            nbsv = 20;
        }

        let polyhedron = ThePolyhedronOfHInter::new(
            surface,
            nbsu as usize,
            nbsv as usize,
            u1new,
            v1new,
            u2new,
            v2new,
        );
        let mut bnd_tool = IntfTool::new();
        let mut box_line = BndBox::new();
        bnd_tool.lin_box(lin, polyhedron.bounding(), &mut box_line);

        for nbseg in 1..=bnd_tool.nb_segments() {
            let mut pinf = bnd_tool.begin_param(nbseg);
            let mut psup = bnd_tool.end_param(nbseg);
            if (psup - pinf) < 1.0e-10 {
                pinf -= 1.0e-10;
                psup += 1.0e-10;
            }
            let polygon = ThePolygonOfHInter::with_range(&line, pinf, psup, 2);
            self.internal_perform(
                &line, &polygon, surface, &polyhedron, u1new, v1new, u2new, v2new,
            );
        }
    }

    /// `AppendIntAna` (`cxx:540-550` = `InterImpl::AppendIntAna`,
    /// `Inter.pxx:968-1003` + `ProcessIntAna`, `InterUtils.pxx:1188-1227`).
    ///
    /// `ProcessIntAna` reports `IsInQuadric` / `IsParallel` through its
    /// `theIsParallel` out-parameter instead of points, and `AppendIntAna` then
    /// only sets `myIsParallel` - it never clears it and appends nothing.
    fn append_int_ana(
        &mut self,
        curve: &dyn Curve,
        surface: &HlrBRepSurfaceView<'_>,
        ana: &IntAnaIntConicQuad,
    ) {
        if !ana.is_done() {
            return;
        }
        if ana.is_in_quadric() || ana.is_parallel() {
            self.is_parallel = true;
            return;
        }
        for i in 1..=ana.nb_points() {
            let w = ana.param_on_conic(i);
            let p = ana.point(i);
            // `ComputeParamsOnQuadric` (`InterUtils.pxx:900-928`).
            let Some((u, v)) = surface_params(surface, p) else {
                continue;
            };
            // `ComputeAppendPoint` (`InterUtils.pxx:1117-1176`).
            if let Some(pt) = compute_append_point(curve, w, surface, u, v) {
                self.append(pt);
            }
        }
    }
}
