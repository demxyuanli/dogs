//! Face–face surface intersection — thin orchestrator over `crate::intpatch`.
//!
//! Ports the orchestration of `IntTools_FaceFace` (TKBO) without the IntPatch
//! walker and `GeomInt_WLApprox`: the two faces' underlying surfaces are
//! classified, supported analytic pairs dispatch to the exact closed forms in
//! `crate::intpatch` (plane∩plane, plane∩sphere, sphere∩sphere, plane∩cylinder)
//! and `occt_geom::intana` — the `IntAna_QuadQuadGeo` conics (plane∩cone,
//! plane∩torus, cylinder×cylinder, cylinder×sphere, sphere×cone, cone×cone,
//! cylinder×cone). Everything else falls back to the intpatch sampling tracer.
//!
//! The output mirrors `IntTools_Curve`: a 3D curve plus the per-face 2D
//! pcurves (`pcurve_full::make_pcurve_full`) over a valid parameter range.
//!
//! ponytail: the cone×cone common-generatrix branch and torus×* (non-plane)
//! pairs still fall through to the grid tracer — OCCT uses the numeric IntPatch
//! walker there, so the boundary is the same; the plane×cylinder generatrix
//! pair collapses to a single line curve. Add per-pair closed forms when a
//! caller needs them.

use std::cmp::Ordering;
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::intpatch::{self, IntersectionCurve, SurfaceIntersection};
use crate::inttools_data::{CurveKind, IntRange};
use crate::pcurve_full::classify_surface_kind;
use crate::shape::Face;

/// Default intersection tolerance used when the caller does not call
/// [`FaceFace::set_tolerance`].
const DEFAULT_TOL: f64 = 1e-6;

// ---------------------------------------------------------------------------
// Data classes
// ---------------------------------------------------------------------------

/// One face/face intersection curve (mirrors `IntTools_Curve`).
///
/// Carries the 3D intersection curve, the parameter range over which it is
/// valid on both faces, and (when computable) the pcurve of the 3D curve on
/// each face — the same payload `IntTools_Curve::SetCurves` stores.
#[derive(Clone)]
pub struct FaceFaceCurve {
    /// The analytic family of the 3D curve.
    pub kind: CurveKind,
    /// The 3D intersection curve.
    pub curve: Arc<dyn Curve>,
    /// Parameter range `[first, last]` over which the curve is valid.
    pub range: IntRange,
    /// Index of the first face (0 — the sorted `FaceFace::face1`).
    pub face1_idx: usize,
    /// Index of the second face (1 — the sorted `FaceFace::face2`).
    pub face2_idx: usize,
    /// Pcurve of the 3D curve on face1 (if computable).
    pub pcurve1: Option<Arc<dyn Curve2d>>,
    /// Pcurve of the 3D curve on face2 (if computable).
    pub pcurve2: Option<Arc<dyn Curve2d>>,
    /// 3D tolerance reached (`IntTools_Curve::Tolerance`).
    pub tolerance: f64,
    /// Tangential tolerance (`IntTools_Curve::TangentialTolerance`).
    pub tangential_tolerance: f64,
}

impl std::fmt::Debug for FaceFaceCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut d = f.debug_struct("FaceFaceCurve");
        d.field("kind", &self.kind).field("range", &self.range);
        d.field("face1_idx", &self.face1_idx).field("face2_idx", &self.face2_idx);
        // Sample a few points so the debug output stays human-readable.
        let (a, b) = (self.range.first, self.range.last);
        let pts: Vec<GpPnt> = (0..3)
            .map(|i| self.curve.d0(a + (b - a) * i as f64 / 2.0))
            .collect();
        d.field("points", &pts)
            .field("has_pcurve1", &self.pcurve1.is_some())
            .field("has_pcurve2", &self.pcurve2.is_some())
            .field("tolerance", &self.tolerance)
            .field("tangential_tolerance", &self.tangential_tolerance)
            .finish()
    }
}

/// The sequence of intersection curves of two faces (`IntTools_FaceFace::Lines`).
#[derive(Debug, Clone, Default)]
pub struct FaceFaceResult {
    curves: Vec<FaceFaceCurve>,
}

impl FaceFaceResult {
    /// Empty result.
    pub fn new() -> Self {
        Self { curves: Vec::new() }
    }

    /// Number of intersection curves.
    pub fn nb_curves(&self) -> usize {
        self.curves.len()
    }

    /// Whether no intersection curves were found.
    pub fn is_empty(&self) -> bool {
        self.curves.is_empty()
    }

    /// The `i`-th intersection curve (panics when `i >= nb_curves`).
    pub fn curve(&self, i: usize) -> &FaceFaceCurve {
        &self.curves[i]
    }

    /// The full slice of intersection curves.
    pub fn curves(&self) -> &[FaceFaceCurve] {
        &self.curves
    }
}

/// Placeholder for `IntTools_Context`.
///
/// The OCCT context caches surface adaptors and topology tools for the IntPatch
/// walker; the thin FaceFace port only needs the two faces and a tolerance, so
/// the context is accepted (for API compatibility) and otherwise unused.
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceFaceContext;

#[path = "int_face_face_helpers.rs"]
mod helpers;
pub(crate) use helpers::{
    cone_from_surface, curve_range, cylinder_from_surface, cylinder_params, lin2d_through,
    line_in_uv_rect, mid, pcurve_of_curve, plane_cylinder_kind, plane_torus_circles, same_curve,
    sphere_from_surface, surface_sort_index,
};

#[path = "int_face_face_bounds.rs"]
pub(crate) mod bounds;
#[path = "int_face_face_make_curve.rs"]
mod make_curve;
#[path = "int_face_face_analytic.rs"]
mod analytic;

// ---------------------------------------------------------------------------
// FaceFace
// ---------------------------------------------------------------------------

/// Intersects the underlying surfaces of two faces.
///
/// Mirrors `IntTools_FaceFace`: `SetParameters` stores the 3D/2D approximation
/// flags (defaults all-on, `1e-7`) for MakeCurve / WLine consumers;
/// `Perform` classifies the surfaces, dispatches to the analytic
/// `intpatch` closed forms when possible and the sampling tracer otherwise,
/// and stores the resulting curves in a [`FaceFaceResult`].
pub struct FaceFace {
    face1: Option<Face>,
    face2: Option<Face>,
    tol: f64,
    fuzzy_value: f64,
    approx: bool,
    approx1: bool,
    approx2: bool,
    tol_approx: f64,
    done: bool,
    tangent_faces: bool,
    result: FaceFaceResult,
    list_of_pnts: Vec<(f64, f64, f64, f64)>,
}

impl FaceFace {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            face1: None,
            face2: None,
            tol: DEFAULT_TOL,
            fuzzy_value: CONFUSION,
            approx: true,
            approx1: true,
            approx2: true,
            tol_approx: 1.0e-7,
            done: false,
            tangent_faces: false,
            result: FaceFaceResult::new(),
            list_of_pnts: Vec::new(),
        }
    }

    /// Sets the first face.
    pub fn set_face1(&mut self, face: Face) {
        self.face1 = Some(face);
    }

    /// Sets the second face.
    pub fn set_face2(&mut self, face: Face) {
        self.face2 = Some(face);
    }

    /// Sets the intersection tolerance.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tol = tol;
    }

    /// `IntTools_FaceFace::SetParameters`.
    pub fn set_parameters(
        &mut self,
        to_approx_c3d: bool,
        to_approx_c2d_on_s1: bool,
        to_approx_c2d_on_s2: bool,
        approximation_tolerance: f64,
    ) {
        self.approx = to_approx_c3d;
        self.approx1 = to_approx_c2d_on_s1;
        self.approx2 = to_approx_c2d_on_s2;
        self.tol_approx = approximation_tolerance;
    }

    /// `IntTools_FaceFace::SetFuzzyValue`.
    pub fn set_fuzzy_value(&mut self, fuzz: f64) {
        self.fuzzy_value = fuzz.max(CONFUSION);
    }

    /// Approximation flags stored by [`Self::set_parameters`].
    pub fn parameters(&self) -> (bool, bool, bool, f64) {
        (self.approx, self.approx1, self.approx2, self.tol_approx)
    }

    /// `IntTools_FaceFace::FuzzyValue`.
    pub fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value
    }

    /// `IntTools_FaceFace::SetList` — EF vertices used as IntPatch start points.
    pub fn set_list(&mut self, list: Vec<(f64, f64, f64, f64)>) {
        self.list_of_pnts = list;
    }

    /// Sets the intersection context (accepted for API compatibility; the thin
    /// port does not use it).
    pub fn set_context(&mut self, _ctx: FaceFaceContext) {
        // no-op: the context caches IntPatch walker state we do not have
    }

    /// Intersects the underlying surfaces of the two faces.
    ///
    /// The surfaces are ordered internally so the "higher" analytic type becomes
    /// `Face1` (`IntTools_FaceFace.cxx:351-357`, `SortTypes`), and the
    /// intersectors express their pcurves in *that* order. OCCT swaps the
    /// pcurves back when `bReverse` is set — `cxx:420-436` for the plane/plane
    /// early return, `cxx:550-563` for the general path (the intersection points
    /// are re-bound the same way, `cxx:595-604`) — so on output
    /// `pcurve1` is the pcurve on the face passed as `face1` and `pcurve2` the
    /// one on `face2`, in the *caller's* order. That is the order
    /// `BOPAlgo_PaveFiller::MakeBlocks` consumes them in: it binds
    /// `aIC.FirstCurve2d()` to `myDS->Shape(nF1)` and `SecondCurve2d()` to
    /// `Shape(nF2)` (`BOPAlgo_PaveFiller_6.cxx:747-748`, `:914`).
    /// Output curves are sorted by curve parameter and deduplicated.
    pub fn perform(&mut self) -> Result<(), String> {
        self.done = false;
        let fa = self.face1.clone().ok_or("FaceFace::perform: face1 not set")?;
        let fb = self.face2.clone().ok_or("FaceFace::perform: face2 not set")?;
        let sa = BRepTool::face_surface(&fa).ok_or("FaceFace::perform: face1 has no surface")?;
        let sb = BRepTool::face_surface(&fb).ok_or("FaceFace::perform: face2 has no surface")?;
        let ka = classify_surface_kind(sa.as_ref());
        let kb = classify_surface_kind(sb.as_ref());

        // OCCT orders the faces so the higher analytic type becomes Face1.
        let reverse = surface_sort_index(ka) < surface_sort_index(kb);
        let (fa, fb, sa, sb, ka, kb) = if reverse {
            let tmp = self.approx1;
            self.approx1 = self.approx2;
            self.approx2 = tmp;
            (fb, fa, sb, sa, kb, ka)
        } else {
            (fa, fb, sa, sb, ka, kb)
        };

        let tol_fa = BRepTool::face_tolerance(&fa);
        let tol_fb = BRepTool::face_tolerance(&fb);
        let tol = self.tol.max(tol_fa + tol_fb);

        self.tangent_faces = false;
        let mut curves =
            self.intersect_surfaces(ka, kb, sa.as_ref(), sb.as_ref(), &fa, &fb, tol)?;

        // `IntTools_FaceFace.cxx:420-436` / `:550-563`: undo the sort for the
        // pcurves, so they come out bound to the caller's Face1/Face2.
        if reverse {
            for c in &mut curves {
                std::mem::swap(&mut c.pcurve1, &mut c.pcurve2);
            }
        }

        if !self.list_of_pnts.is_empty() {
            for &(u1, v1, u2, v2) in &self.list_of_pnts {
                let p1 = sa.d0(u1, v1);
                let p2 = sb.d0(u2, v2);
                if p1.distance(&p2) <= tol {
                    self.tangent_faces = true;
                    break;
                }
            }
        }

        // Output sorted by curve parameter, then deduplicated.
        curves.sort_by(|x, y| x.range.first.partial_cmp(&y.range.first).unwrap_or(Ordering::Equal));
        curves.dedup_by(|a, b| same_curve(a, b));

        self.result = FaceFaceResult { curves };
        if !self.tangent_faces {
            self.compute_tol_reached_3d(&fa, &fb);
        }
        self.done = true;
        Ok(())
    }

    /// Whether the intersection was computed successfully.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// The intersection curves (empty when not done or no intersection).
    pub fn result(&self) -> FaceFaceResult {
        self.result.clone()
    }

    /// The sorted first face.
    pub fn face1(&self) -> Option<&Face> {
        self.face1.as_ref()
    }

    /// The sorted second face.
    pub fn face2(&self) -> Option<&Face> {
        self.face2.as_ref()
    }

    /// Whether the faces are tangent (coincident over a region).
    pub fn tangent_faces(&self) -> bool {
        self.tangent_faces
    }

    // ------------------------------------------------------------------
    // Internal dispatch
    // ------------------------------------------------------------------

    fn intersect_surfaces(
        &mut self,
        kind_a: SurfaceKind,
        kind_b: SurfaceKind,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        match (kind_a, kind_b) {
            (SurfaceKind::Plane, SurfaceKind::Plane) => self.plane_plane(sa, sb, fa, fb, tol),
            (SurfaceKind::Plane, SurfaceKind::Sphere)
            | (SurfaceKind::Sphere, SurfaceKind::Plane) => {
                self.plane_sphere(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
                self.sphere_sphere(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Plane, SurfaceKind::Cylinder)
            | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
                self.plane_cylinder(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Plane, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Plane) => {
                self.plane_cone(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Plane, SurfaceKind::Torus) | (SurfaceKind::Torus, SurfaceKind::Plane) => {
                self.plane_torus(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Cylinder, SurfaceKind::Cylinder) => {
                self.cylinder_cylinder(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Sphere, SurfaceKind::Cylinder) => {
                self.cylinder_sphere(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Sphere, SurfaceKind::Cone) => self.sphere_cone(sa, sb, fa, fb, tol),
            (SurfaceKind::Cone, SurfaceKind::Cylinder) => {
                self.cylinder_cone(sa, sb, fa, fb, tol)
            }
            (SurfaceKind::Cone, SurfaceKind::Cone) => self.cone_cone(sa, sb, fa, fb, tol),
            _ => self.general(sa, sb, fa, fb, tol),
        }
    }

    /// General fallback: `GeomInt_IntSS` over `CorrectSurfaceBoundaries` UV
    /// boxes (domains from the faces). Walking / analytic lines go through
    /// `MakeCurve` (degree-1 B-spline when WLApprox is not done). Marks
    /// `tangent_faces` on coincidence.
    fn general(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let sa_arc = BRepTool::face_surface(fa).unwrap_or_else(|| Arc::from(sa.clone_dyn()));
        let sb_arc = BRepTool::face_surface(fb).unwrap_or_else(|| Arc::from(sb.clone_dyn()));
        let ka = classify_surface(sa_arc.as_ref());
        let kb = classify_surface(sb_arc.as_ref());
        let (u1min, u1max, v1min, v1max) = bounds::corrected_uv_box(fa, kb, ka, self.tol);
        let (u2min, u2max, v2min, v2max) = bounds::corrected_uv_box(fb, ka, kb, self.tol);
        let sa_b: Arc<dyn Surface> =
            Arc::new(bounds::BoundedSurface::new(sa_arc, u1min, u1max, v1min, v1max));
        let sb_b: Arc<dyn Surface> =
            Arc::new(bounds::BoundedSurface::new(sb_arc, u2min, u2max, v2min, v2max));

        let mut iss = crate::geom_int::IntSS::new();
        iss.load(
            crate::geom_int::TopolTool::from_face(fa, sa_b.as_ref()),
            crate::geom_int::TopolTool::from_face(fb, sb_b.as_ref()),
            sa_b.clone(),
            sb_b.clone(),
        );
        iss.perform_loaded(tol, self.approx, self.approx1, self.approx2);
        if iss.tangent_faces() {
            self.tangent_faces = true;
            return Ok(Vec::new());
        }
        if !iss.lines().is_empty() {
            let mut out = Vec::new();
            for l in iss.lines() {
                out.push(FaceFaceCurve {
                    kind: CurveKind::BSpline,
                    curve: l.curve.clone(),
                    range: curve_range(l.curve.as_ref()),
                    face1_idx: 0,
                    face2_idx: 1,
                    pcurve1: l.pcurve1.clone(),
                    pcurve2: l.pcurve2.clone(),
                    tolerance: iss.tol_reached_3d(),
                    tangential_tolerance: 0.0,
                });
            }
            let lim = (tol * 4.0).max(0.2);
            if curves_lie_on_both(&out, sa_b.as_ref(), sb_b.as_ref(), lim) {
                return Ok(out);
            }
        }

        let mut out = Vec::new();
        match intpatch::surface_surface_intersection(sa_b.as_ref(), sb_b.as_ref(), tol) {
            SurfaceIntersection::Curves(ics) => {
                for ic in ics {
                    let wl = crate::int_tools_wline::WLine::from_intersection_curve(&ic);
                    if wl.nb_pnts() >= 2 {
                        let pieces =
                            self.make_curve_walking(&wl, sa_b.as_ref(), sb_b.as_ref(), fa, fb);
                        if pieces.is_empty() {
                            out.push(self.curve_from_ic(ic, CurveKind::BSpline, fa, fb));
                        } else {
                            out.extend(pieces);
                        }
                    } else {
                        out.push(self.curve_from_ic(ic, CurveKind::BSpline, fa, fb));
                    }
                }
            }
            SurfaceIntersection::Coincident => {
                self.tangent_faces = true;
            }
            SurfaceIntersection::None => {}
        }
        Ok(out)
    }

    /// Wrap an `intpatch` intersection curve into a `FaceFaceCurve`, building
    /// the per-face pcurves over the curve's parameter range.
    fn curve_from_ic(
        &self,
        ic: IntersectionCurve,
        kind: CurveKind,
        fa: &Face,
        fb: &Face,
    ) -> FaceFaceCurve {
        let curve = ic.curve;
        let range = curve_range(curve.as_ref());
        self.finish_curve(curve, kind, range, fa, fb)
    }

    /// Assemble a `FaceFaceCurve` from a curve, kind, range and the two faces.
    fn finish_curve(
        &self,
        curve: Arc<dyn Curve>,
        kind: CurveKind,
        range: IntRange,
        fa: &Face,
        fb: &Face,
    ) -> FaceFaceCurve {
        let pcurve1 = pcurve_of_curve(&curve, range, fa);
        let pcurve2 = pcurve_of_curve(&curve, range, fb);
        FaceFaceCurve {
            kind,
            curve,
            range,
            face1_idx: 0,
            face2_idx: 1,
            pcurve1,
            pcurve2,
            tolerance: 0.0,
            tangential_tolerance: 0.0,
        }
    }
}

impl Default for FaceFace {
    fn default() -> Self {
        Self::new()
    }
}

/// Reject IntSS 3D curves that leave either surface; the general tracer runs
/// instead (`GeomInt` `RejetLigne` analogue).
fn curves_lie_on_both(
    curves: &[FaceFaceCurve],
    sa: &dyn Surface,
    sb: &dyn Surface,
    lim: f64,
) -> bool {
    if curves.is_empty() {
        return false;
    }
    for c in curves {
        let (a, b) = (c.range.first, c.range.last);
        if !a.is_finite() || !b.is_finite() || b <= a {
            return false;
        }
        for i in 0..=12 {
            let t = a + (b - a) * i as f64 / 12.0;
            let p = c.curve.d0(t);
            if intpatch::distance_to_surface(&p, sa) > lim
                || intpatch::distance_to_surface(&p, sb) > lim
            {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
#[path = "int_face_face_tests.rs"]
mod tests;
