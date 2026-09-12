//! Edge/face full intersection — Phase 17b.
//!
//! Port of `IntTools_EdgeFace` (TKBO): computes the common parts between an
//! edge and a face in 3D space.
//!
//! A common part is either:
//!
//! * an **edge** — a (sub-)range of the edge that lies on the face surface
//!   (coincidence), encoded as [`CommonPartType::Edge`];
//! * a **point** — a single parameter at which the edge touches/crosses the
//!   face, encoded as [`CommonPartType::Vertex`] with
//!   [`CommonPrt::vertex_parameter1`] (the range is kept, matching
//!   `IntTools_EdgeFace::MakeType`).
//!
//! Flow (mirrors `IntTools_EdgeFace::Perform`):
//!
//! 1. `check_data` — edge geometry accessibility (degenerated / non-geometric).
//! 2. tolerance preparation (`tol_e + tol_f`, with the B-spline special case);
//! 3. quick coincidence check ([`EdgeFace::is_coincident`]) when requested;
//! 4. [`BeanFaceIntersector`] over the bean range → candidate ranges;
//! 5. `IsProjectable` filter (distance + face 2D restriction);
//! 6. `MakeType` + `CheckTouch` classification per range (point vs edge);
//! 7. line/cylinder and circle/plane touch refinement.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::{Curve, Surface};

use crate::bean_face::BeanFaceIntersector;
use crate::brep_surface::{is_planar, surface_closest_params, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::edge_face_kind::*;
use crate::inttools_data::{CommonPartType, CommonPrt, IntRange};
use crate::shape::{Edge, Face};

/// Edge/face intersection algorithm. Port of `IntTools_EdgeFace`.
///
/// Not `Debug`-derived: the cached `Arc<dyn Curve>`/`Arc<dyn Surface>` handles
/// are not `Debug`. `Clone` is provided (Arc handles clone cheaply).
#[derive(Clone)]
pub struct EdgeFace {
    pub(crate) edge: Edge,
    pub(crate) face: Face,
    pub(crate) range: IntRange,
    pub(crate) fuzzy_value: f64,
    pub(crate) quick_coincidence_check: bool,
    pub(crate) curve: Option<Arc<dyn Curve>>,
    pub(crate) surface: Option<Arc<dyn Surface>>,
    pub(crate) criteria: f64,
    pub(crate) is_done: bool,
    pub(crate) error_status: i32,
    pub(crate) common_parts: Vec<CommonPrt>,
    pub(crate) face_classifier: Option<FClass2d>,
    pub(crate) min_distance: f64,
}

impl Default for EdgeFace {
    fn default() -> Self {
        Self::new()
    }
}

impl EdgeFace {
    /// Empty constructor (`IntTools_EdgeFace()`).
    pub fn new() -> Self {
        Self {
            edge: Edge::new(),
            face: Face::new(),
            range: IntRange::new_unchecked(f64::NEG_INFINITY, f64::INFINITY),
            fuzzy_value: CONFUSION,
            quick_coincidence_check: false,
            curve: None,
            surface: None,
            criteria: CONFUSION,
            is_done: false,
            error_status: 1,
            common_parts: Vec::new(),
            face_classifier: None,
            min_distance: f64::MAX,
        }
    }

    // ---- setters / getters --------------------------------------------------

    /// Sets the edge for intersection.
    pub fn set_edge(&mut self, edge: Edge) {
        self.edge = edge;
    }

    /// Returns the edge.
    pub fn edge(&self) -> &Edge {
        &self.edge
    }

    /// Sets the face for intersection.
    pub fn set_face(&mut self, face: Face) {
        self.face = face;
    }

    /// Returns the face.
    pub fn face(&self) -> &Face {
        &self.face
    }

    /// Sets the boundaries of the edge to process.
    pub fn set_range(&mut self, first: f64, last: f64) {
        self.range = IntRange::new_unchecked(first, last);
    }

    /// Returns the processing range.
    pub fn range(&self) -> IntRange {
        self.range
    }

    /// Sets the fuzzy value (clamped to `Precision::Confusion`).
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.fuzzy_value = v.max(CONFUSION);
    }

    /// Returns the fuzzy value.
    pub fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value
    }

    /// Sets the quick coincidence check flag.
    pub fn set_quick_coincidence_check(&mut self, b: bool) {
        self.quick_coincidence_check = b;
    }

    /// Returns the quick coincidence check flag.
    pub fn is_coincidence_checked_quickly(&self) -> bool {
        self.quick_coincidence_check
    }

    // ---- performing ----------------------------------------------------------

    /// Launches the intersection.
    ///
    /// `Err` is returned only for hard setup failures (missing geometry, invalid
    /// range). Algorithmic failures set [`error_status`](Self::error_status)
    /// (`2`/`3`/`4` per OCCT) and still return `Ok`.
    pub fn perform(&mut self) -> Result<(), String> {
        self.common_parts.clear();
        self.error_status = 0;
        self.min_distance = f64::MAX;
        self.check_data();
        if self.error_status != 0 {
            return Ok(());
        }

        let curve = BRepTool::edge_curve_world(&self.edge)
            .ok_or_else(|| "EdgeFace::perform: edge has no curve".to_string())?;
        let surface = BRepTool::face_surface_world(&self.face)
            .ok_or_else(|| "EdgeFace::perform: face has no surface".to_string())?;
        let (ef, el) = BRepTool::edge_parameters(&self.edge);
        if !ef.is_finite() || !el.is_finite() || el - ef <= 1e-15 {
            return Err("EdgeFace::perform: edge has an empty or unbounded parameter range".into());
        }
        // Default the processing range to the whole edge when not explicitly set.
        if !self.range.is_valid() || !self.range.first.is_finite() || !self.range.last.is_finite() {
            self.range = IntRange::new_unchecked(ef, el);
        }

        self.is_done = false;
        self.curve = Some(curve.clone());
        self.surface = Some(surface.clone());

        let c_kind = curve_kind(curve.as_ref());
        let s_kind = surface_kind(surface.as_ref());

        // Prepare myCriteria.
        let fuzz = self.fuzzy_value * 0.5;
        let tol_f = BRepTool::face_tolerance(&self.face) + fuzz;
        let tol_e = BRepTool::edge_tolerance(&self.edge) + fuzz;
        self.criteria = match c_kind {
            CurveKind::BSpline => {
                let diff1 = tol_e / tol_f.max(1e-300);
                let diff2 = tol_f / tol_e.max(1e-300);
                if diff1 > 100.0 || diff2 > 100.0 {
                    tol_e.max(tol_f)
                } else {
                    1.5 * tol_e + tol_f
                }
            }
            _ => tol_e + tol_f,
        };

        // 2D classifier for the face's UV restriction (used by coincidence and
        // projectability checks).
        let cl_tol = BRepTool::face_tolerance(&self.face).max(PCONFUSION);
        self.face_classifier = Some(FClass2d::new(&self.face, cl_tol)?);

        if self.quick_coincidence_check && self.is_coincident() {
            let mut cp = CommonPrt::new();
            cp.part_type = CommonPartType::Edge;
            cp.range = self.range;
            cp.face = Some(self.face.0.clone());
            let p1 = curve.d0(self.range.first);
            let p2 = curve.d0(self.range.last);
            cp.set_bounding_points(p1, p2);
            self.common_parts.push(cp);
            self.is_done = true;
            return Ok(());
        }

        let mut intersector = BeanFaceIntersector::new();
        intersector.initialize(curve.clone(), surface.clone(), tol_e, tol_f);
        intersector.set_bean_parameters(self.range.first, self.range.last);
        // `IntTools_EdgeFace::Perform` does not call `SetSurfaceParameters`;
        // the adaptor/surface ranges come from `BeanFaceIntersector::Init`.
        intersector.perform()?;
        if !intersector.is_done() {
            return Ok(());
        }
        self.min_distance = intersector.minimal_square_distance().sqrt();

        for r in intersector.result() {
            let mid = 0.5 * (r.first + r.last);
            if self.is_projectable(mid) {
                let mut cp = CommonPrt::new();
                cp.range = r;
                cp.face = Some(self.face.0.clone());
                let p1 = curve.d0(r.first);
                let p2 = curve.d0(r.last);
                cp.set_bounding_points(p1, p2);
                self.common_parts.push(cp);
            }
        }

        let nb = self.common_parts.len();
        for i in 0..nb {
            let mut cp = self.common_parts[i].clone();
            self.make_type(&mut cp);
            self.common_parts[i] = cp;
        }

        // Line/Cylinder and Circle/Plane common-part refinement.
        let special = (c_kind == CurveKind::Line && s_kind == SurfaceKind::Cylinder)
            || (c_kind == CurveKind::Circle
                && s_kind == SurfaceKind::Plane
                && !is_coplanar(curve.as_ref(), surface.as_ref())
                && !is_radius(curve.as_ref(), surface.as_ref(), self.criteria));
        if special {
            self.refine_touch_parts();
        }

        self.is_done = true;
        Ok(())
    }

    /// Whether the computation was successful.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// Completion code: `0` success; `1` not started; `2`/`3` invalid input;
    /// `4` projection failed.
    pub fn error_status(&self) -> i32 {
        self.error_status
    }

    /// The resulting common parts.
    pub fn common_parts(&self) -> &[CommonPrt] {
        &self.common_parts
    }

    /// `IntTools_EdgeFace::MinimalDistance`.
    pub fn minimal_distance(&self) -> f64 {
        self.min_distance
    }

    /// For each common part, `Some(t)` when the part is a `TopAbs_VERTEX`
    /// at edge parameter `t`, `None` when it is a coincident edge sub-range.
    pub fn point_parameters(&self) -> Vec<Option<f64>> {
        self.common_parts
            .iter()
            .map(|cp| {
                if cp.part_type == CommonPartType::Vertex {
                    cp.vertex_parameter1
                        .or_else(|| Some(0.5 * (cp.range.first + cp.range.last)))
                } else {
                    None
                }
            })
            .collect()
    }

    // ---- internals -----------------------------------------------------------

    /// `CheckData`: sets `error_status` to `2` (degenerated edge) or `3`
    /// (non-geometric edge).
    fn check_data(&mut self) {
        if BRepTool::is_degenerated(&self.edge) {
            self.error_status = 2;
        }
        if BRepTool::edge_curve(&self.edge).is_none() {
            self.error_status = 3;
        }
    }

    /// Whether the edge is entirely on the face, sampled over the range.
    /// Port of `IntTools_EdgeFace::IsCoincident`.
    fn is_coincident(&self) -> bool {
        let curve = self.curve.clone().expect("curve set");
        let surface = self.surface.clone().expect("surface set");
        let a_nb_seg = if curve_kind(curve.as_ref()) == CurveKind::Line
            && surface_kind(surface.as_ref()) == SurfaceKind::Plane
        {
            2
        } else {
            23
        };
        let a_tresh = 0.5;
        let a_tresh_idx_f = ((a_nb_seg + 1) as f64 * 0.25) as i32;
        let a_tresh_idx_l = ((a_nb_seg + 1) as f64 * 0.75) as i32;

        let (mut a_t1, mut a_t2) = (self.range.first, self.range.last);
        if a_t2 - a_t1 <= 1e-12 {
            return false;
        }
        let a_bnd_shift = 0.01 * (a_t2 - a_t1);
        a_t1 += a_bnd_shift;
        a_t2 -= a_bnd_shift;
        if a_t2 <= a_t1 {
            return false;
        }
        let d_t = (a_t2 - a_t1) / a_nb_seg as f64;

        let mut is_classified = false;
        let mut i_cnt = 0;
        for i in 0..=a_nb_seg {
            let a_t = a_t1 + i as f64 * d_t;
            let a_p = curve.d0(a_t);
            let (u, v, a_d) = self.project_point(&a_p);
            if a_d > self.criteria {
                if a_d > 100.0 * self.criteria {
                    return false;
                }
                continue;
            }
            i_cnt += 1;
            if ((0 < i) && (i < a_tresh_idx_f)) || ((a_tresh_idx_l < i) && (i < a_nb_seg)) {
                continue;
            }
            if is_classified && (i != a_nb_seg) {
                continue;
            }
            let state = self
                .face_classifier
                .as_ref()
                .map(|c| c.perform(GpPnt2d::new(u, v)))
                .unwrap_or(FaceState::In);
            if state == FaceState::Out {
                return false;
            }
            if i != 0 {
                is_classified = true;
            }
        }
        let a_coeff = i_cnt as f64 / (a_nb_seg + 1) as f64;
        a_coeff > a_tresh
    }

    /// Whether the curve point at parameter `t` is on the face within
    /// `myCriteria` and inside its 2D restriction.
    /// Port of `IntTools_EdgeFace::IsProjectable`.
    fn is_projectable(&self, t: f64) -> bool {
        let curve = self.curve.clone().expect("curve set");
        let p = curve.d0(t);
        let (u, v, dist) = self.project_point(&p);
        if dist > self.criteria {
            return false;
        }
        match &self.face_classifier {
            Some(cl) => cl.perform(GpPnt2d::new(u, v)) != FaceState::Out,
            None => true,
        }
    }

    /// Signed distance from the curve point at `t` to the surface, minus
    /// `myCriteria`. Port of `IntTools_EdgeFace::DistanceFunction`.
    pub(crate) fn distance_function(&self, t: f64) -> f64 {
        let curve = self.curve.clone().expect("curve set");
        let surface = self.surface.clone().expect("surface set");
        let p = curve.d0(t);
        if let Some(d) = is_eq_distance(&p, surface.as_ref()) {
            return d - self.criteria;
        }
        let (_, _, dist) = self.project_point(&p);
        dist - self.criteria
    }

    /// Project a point onto the face surface: `(u, v, distance)`. Planar faces
    /// are projected analytically (exact distance); others use the grid+refine
    /// projector.
    fn project_point(&self, p: &GpPnt) -> (f64, f64, f64) {
        let surface = self.surface.clone().expect("surface set");
        if is_planar(surface.as_ref(), 6, 6, 1e-6) {
            plane_projection(surface.as_ref(), p)
        } else {
            let (u, v) = surface_closest_params(surface.as_ref(), p, 24, 24);
            let q = surface.d0(u, v);
            (u, v, p.distance(&q))
        }
    }

    /// Distance from the curve point at `t` to the face surface (unsigned).
    fn surface_distance_at(&self, t: f64) -> f64 {
        let curve = self.curve.clone().expect("curve set");
        let p = curve.d0(t);
        let (_, _, d) = self.project_point(&p);
        d
    }

    /// Sample the curve–surface distance over `[t0, t1]`. Returns
    /// `(min, max, param_at_min)` with golden-section refinement around the min.
    pub(crate) fn distance_profile(&self, t0: f64, t1: f64, n: usize) -> (f64, f64, f64) {
        let curve = self.curve.clone().expect("curve set");
        if t1 <= t0 {
            let p = curve.d0(t0);
            let (_, _, d) = self.project_point(&p);
            return (d, d, t0);
        }
        let mut min_d = f64::INFINITY;
        let mut max_d = 0.0;
        let mut min_t = t0;
        for i in 0..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let p = curve.d0(t);
            let (_, _, d) = self.project_point(&p);
            if d < min_d {
                min_d = d;
                min_t = t;
            }
            if d > max_d {
                max_d = d;
            }
        }
        let h = (t1 - t0) / n as f64;
        let lo = (min_t - h).max(t0);
        let hi = (min_t + h).min(t1);
        if hi > lo {
            let (rt, rd) = golden_1d(&|t: f64| self.surface_distance_at(t), lo, hi, 1e-10);
            if rd < min_d {
                min_d = rd;
                min_t = rt;
            }
        }
        (min_d, max_d, min_t)
    }

}
/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    const GOLD: f64 = 0.6180339887498949;
    let mut a = lo;
    let mut b = hi;
    let mut c = b - GOLD * (b - a);
    let mut d = a + GOLD * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    while (b - a) > eps {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - GOLD * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GOLD * (b - a);
            fd = f(d);
        }
    }
    let x = 0.5 * (a + b);
    (x, f(x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::inttools::edge_face_intersections;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn box_face(idx: usize) -> (crate::brep_extrema::test_box::UnitBox, Face) {
        let ub = unit_box();
        let face = ub.faces[idx].clone();
        (ub, face)
    }

    #[test]
    fn default_state_is_not_started() {
        let ef = EdgeFace::new();
        assert!(!ef.is_done());
        assert_eq!(ef.error_status(), 1);
        assert!(ef.common_parts().is_empty());
        assert_eq!(ef.fuzzy_value(), CONFUSION);
        assert!(!ef.is_coincidence_checked_quickly());
    }

    #[test]
    fn line_in_face_plane_is_edge_part() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Edge lying in the face plane, crossing the face interior.
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.0), &GpPnt::new(1.0, 0.5, 0.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);

        let cps = ef.common_parts();
        assert_eq!(cps.len(), 1, "common parts: {cps:?}");
        assert_eq!(cps[0].part_type, CommonPartType::Edge);
        assert!(cps[0].range.length() > 0.5, "coincident range: {cps:?}");
        let params = ef.point_parameters();
        assert_eq!(params.len(), 1);
        assert!(params[0].is_none(), "in-plane edge is not a point: {params:?}");

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn line_piercing_face_is_point_and_matches_approximate() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 2.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);

        let cps = ef.common_parts();
        assert_eq!(cps.len(), 1, "common parts: {cps:?}");
        let params = ef.point_parameters();
        assert_eq!(params.len(), 1);
        let t = params[0].expect("piercing line is a point");
        let curve = BRepTool::edge_curve(&e).expect("edge curve");
        let p = curve.d0(t);
        assert!(
            p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-4,
            "intersection point {p:?} at t={t}"
        );

        // Cross-check with the existing approximate solver.
        let hits = edge_face_intersections(&e, &face, 1e-9);
        assert_eq!(hits.len(), 1, "approximate hits: {hits:?}");
        assert!(p.distance(&hits[0].1) < 1e-4, "{p:?} vs {:?}", hits[0].1);

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn line_outside_face_is_empty() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Pierces the plane but outside the face's UV domain.
        let e = b.make_edge_segment(&GpPnt::new(2.0, 2.0, -1.0), &GpPnt::new(2.0, 2.0, 1.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 2.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);
        assert!(
            ef.common_parts().is_empty(),
            "common parts: {:?}",
            ef.common_parts()
        );

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn edge_endpoint_on_face_is_point() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Starts exactly on the face, then leaves along +Z.
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, 0.0), &GpPnt::new(0.5, 0.5, 1.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);

        let params = ef.point_parameters();
        assert_eq!(params.len(), 1, "parts: {:?}", ef.common_parts());
        let t = params[0].expect("endpoint touch is a point");
        assert!(t.abs() < 1e-4, "touch parameter {t}, want ~0");

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn quick_coincidence_detects_in_plane_edge() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.0), &GpPnt::new(1.0, 0.5, 0.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.set_quick_coincidence_check(true);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        let cps = ef.common_parts();
        assert_eq!(cps.len(), 1, "common parts: {cps:?}");
        assert_eq!(cps[0].part_type, CommonPartType::Edge);
        assert!(cps[0].range.length() > 0.5);

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn quick_coincidence_rejects_offset_line() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Parallel to the plane, 0.5 above it — not coincident.
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.5), &GpPnt::new(1.0, 0.5, 0.5));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.set_quick_coincidence_check(true);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        // Offset line never touches the surface: no common parts.
        assert!(
            ef.common_parts().is_empty(),
            "common parts: {:?}",
            ef.common_parts()
        );

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn non_geometric_edge_sets_error_status() {
        let (ub, face) = box_face(0);
        // An edge with no registered curve is "non-geometric" → error status 3.
        let e = Edge::new();
        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.perform().expect("perform does not hard-fail");
        assert!(!ef.is_done());
        assert_eq!(ef.error_status(), 3);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn fuzzy_value_clamped_to_confusion() {
        let mut ef = EdgeFace::new();
        ef.set_fuzzy_value(0.0);
        assert_eq!(ef.fuzzy_value(), CONFUSION);
        ef.set_fuzzy_value(1e-3);
        assert_eq!(ef.fuzzy_value(), 1e-3);
    }
}
