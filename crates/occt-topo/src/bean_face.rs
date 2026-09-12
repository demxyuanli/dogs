//! Edge bean / face intersection ranges — the numerical heart of the precise
//! boolean.
//!
//! Port of `IntTools_BeanFaceIntersector` (TKBO). Given a curve *bean* (the
//! part of an edge over a parameter window) and a surface, compute the
//! parameter ranges of the curve that lie **on** the surface within the
//! combined edge/face tolerance. The real face boundary is *not* taken into
//! account (as in OCCT): the results are ranges of the curve whose points are
//! within `tol_e + tol_f` of the face surface, inside the surface parameter
//! window `[umin,umax]×[vmin,vmax]`.
//!
//! The `Perform()` flow mirrors the OCCT dispatch:
//!
//! 1. **Line/plane** → `ComputeLinePlane` (substitute the line into the plane
//!    equation → one root, or the whole range when the line lies in the plane).
//! 2. **`FastComputeAnalytic`** — analytic coincidence / no-intersection
//!    verdicts for conic curves against quadrics (plane/circle/ellipse,
//!    cylinder/line/circle, sphere/line). Returns `true` when it can conclude
//!    decisively; otherwise computation continues.
//! 3. **Range manager** (`IntTools_MarkedRangeSet` semantics) + coincidence
//!    scan (`TestComputeCoinside`).
//! 4. **Localized path** for high-degree NURBS surfaces: box culling of
//!    subdivided curve cells, then exact curve–surface intersection per
//!    surviving cell.
//! 5. **Exact path**: `IntCurveSurface_HInter` (`crate::intcurvesurface`) for
//!    discrete points / coincidence segments, then `Extrema_ExtCS`
//!    (`occt_geom::extrema_surf`) for near-surface spans, then boundary
//!    completion (`ComputeNearRangeBoundaries`).
//!
//! `ponytail:` `dyn Surface` cannot be downcast, so the analytic fast paths
//! classify curves/surfaces by sampling geometric invariants (the established
//! pattern in `crate::intcurvesurface`), and the local point-to-surface solve
//! of OCCT's `Extrema_GenLocateExtPS` is replaced by the grid+refine projector
//! `crate::brep_surface::surface_closest_params`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::{Curve, Surface};

use crate::bean_face_kind::{
    classify_fast_curve, plane_geometry, sphere_geometry, FastCurveKind,
};
use crate::bean_face_range::MarkedRangeSet;
use crate::brep_surface::{classify_surface, surface_closest_params, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::inttools_data::IntRange;
use crate::inttools_range::IntContext;
use crate::meshing::range_splitter::{classify_surface as classify_surface_mesh, SurfaceType};

const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// BeanFaceIntersector
// ---------------------------------------------------------------------------

/// Computes the parameter ranges of a curve bean that lie on a surface within
/// the combined edge/face tolerance. Port of `IntTools_BeanFaceIntersector`.
pub struct BeanFaceIntersector {
    pub(crate) curve: Option<Arc<dyn Curve>>,
    pub(crate) surface: Option<Arc<dyn Surface>>,
    pub(crate) first_parameter: f64,
    pub(crate) last_parameter: f64,
    pub(crate) umin: f64,
    pub(crate) umax: f64,
    pub(crate) vmin: f64,
    pub(crate) vmax: f64,
    pub(crate) bean_tolerance: f64,
    pub(crate) face_tolerance: f64,
    pub(crate) curve_resolution: f64,
    pub(crate) criteria: f64,
    pub(crate) range_manager: MarkedRangeSet,
    pub(crate) context: Option<IntContext>,
    pub(crate) results: Vec<IntRange>,
    pub(crate) is_done: bool,
    pub(crate) min_sq_distance: f64,
}

impl BeanFaceIntersector {
    /// Create an empty intersector; call [`initialize`](Self::initialize)
    /// before [`perform`](Self::perform).
    pub fn new() -> Self {
        Self {
            curve: None,
            surface: None,
            first_parameter: 0.0,
            last_parameter: 0.0,
            umin: 0.0,
            umax: 0.0,
            vmin: 0.0,
            vmax: 0.0,
            bean_tolerance: 0.0,
            face_tolerance: 0.0,
            curve_resolution: PCONFUSION,
            criteria: CONFUSION,
            range_manager: MarkedRangeSet::new(),
            context: None,
            results: Vec::new(),
            is_done: false,
            min_sq_distance: f64::MAX,
        }
    }

    /// Initialise with a curve bean and a face surface plus the edge/face
    /// tolerances. The surface parameter window defaults to the surface's own
    /// ranges (override with [`set_surface_parameters`](Self::set_surface_parameters)).
    pub fn initialize(&mut self, curve: Arc<dyn Curve>, surface: Arc<dyn Surface>, tol_e: f64, tol_f: f64) {
        self.curve = Some(curve);
        self.surface = Some(surface);
        self.bean_tolerance = tol_e;
        self.face_tolerance = tol_f;
        self.criteria = tol_e + tol_f;
        let a = self.curve().first_parameter();
        let b = self.curve().last_parameter();
        let (first, last) = if a.is_finite() && b.is_finite() && b >= a { (a, b) } else { (0.0, 0.0) };
        self.first_parameter = first;
        self.last_parameter = last;
        self.curve_resolution = self.resolution(self.criteria);
        let (u0, u1) = self.surface().u_range();
        let (v0, v1) = self.surface().v_range();
        self.set_surface_parameters(u0, u1, v0, v1);
        self.results.clear();
        self.is_done = false;
        self.min_sq_distance = f64::MAX;
    }

    /// Initialise from an edge and a face (`IntTools_BeanFaceIntersector::Init`
    /// of `TopoDS_Edge` / `TopoDS_Face`). Criteria include `Precision::Confusion`.
    pub fn initialize_edge_face(&mut self, edge: &crate::shape::Edge, face: &crate::shape::Face) {
        let Some(curve) = BRepTool::edge_curve(edge) else {
            return;
        };
        let Some(surface) = BRepTool::face_surface(face) else {
            return;
        };
        let tol_e = BRepTool::edge_tolerance(edge);
        let tol_f = BRepTool::face_tolerance(face);
        self.initialize(curve, surface, tol_e, tol_f);
        self.criteria = tol_e + tol_f + CONFUSION;
        self.curve_resolution = self.resolution(self.criteria);
        let (ef, el) = BRepTool::edge_parameters(edge);
        if ef.is_finite() && el.is_finite() && el >= ef {
            self.set_bean_parameters(ef, el);
        }
    }

    /// Restrict the curve bean to `[first, last]`.
    pub fn set_bean_parameters(&mut self, first: f64, last: f64) {
        self.first_parameter = first;
        self.last_parameter = last;
    }

    /// Restrict the surface parameter window to `[umin,umax]×[vmin,vmax]`.
    pub fn set_surface_parameters(&mut self, umin: f64, umax: f64, vmin: f64, vmax: f64) {
        self.umin = umin;
        self.umax = umax;
        self.vmin = vmin;
        self.vmax = vmax;
    }

    /// Attach an intersection context (cached, for interface compatibility; the
    /// projector used here operates directly on the supplied surface).
    pub fn set_context(&mut self, ctx: IntContext) {
        self.context = Some(ctx);
    }

    /// Whether the last [`perform`](Self::perform) completed.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// The parameter ranges of the curve lying on the surface.
    pub fn result(&self) -> Vec<IntRange> {
        self.results.clone()
    }

    /// The minimal (squared) distance found between the curve and the surface.
    pub fn minimal_square_distance(&self) -> f64 {
        self.min_sq_distance
    }

    // ---- accessors -------------------------------------------------------

    pub(crate) fn curve(&self) -> &dyn Curve {
        self.curve.as_ref().expect("BeanFaceIntersector: curve not set").as_ref()
    }

    pub(crate) fn surface(&self) -> &dyn Surface {
        self.surface.as_ref().expect("BeanFaceIntersector: surface not set").as_ref()
    }

    pub(crate) fn curve_d0(&self, u: f64) -> GpPnt {
        self.curve().d0(u)
    }

    pub(crate) fn surface_d0(&self, u: f64, v: f64) -> GpPnt {
        self.surface().d0(u, v)
    }

    /// Approximate `BRepAdaptor_Curve::Resolution(tol)`: the parameter
    /// increment over which the curve deviates from a straight chord by ~`tol`.
    /// Uses the current bean window (unbounded curves must have it set via
    /// [`set_bean_parameters`](Self::set_bean_parameters)).
    pub(crate) fn resolution(&self, tol: f64) -> f64 {
        let c = self.curve();
        let a = self.first_parameter;
        let b = self.last_parameter;
        let span = if a.is_finite() && b.is_finite() && b > a { b - a } else { 0.0 };
        if span.abs() <= 1e-30 || tol <= 0.0 {
            return PCONFUSION;
        }
        let n = 8;
        let h = span / n as f64 * 0.5;
        let mut max_speed = 0.0f64;
        for i in 0..=n {
            let u = a + span * i as f64 / n as f64;
            let speed = c.d0(u - h).distance(&c.d0(u + h)) / (2.0 * h);
            if speed > max_speed {
                max_speed = speed;
            }
        }
        if max_speed <= 1e-30 {
            return tol.max(PCONFUSION);
        }
        (tol / max_speed).max(PCONFUSION)
    }

    // ---- Distance ----------------------------------------------------------

    /// Closest surface parameters `(u, v)` of `p` and the surface distance.
    ///
    /// Planes and spheres are projected **exactly** (the reconstructed analytic
    /// frame), because the grid+refine projector `surface_closest_params` is
    /// only ~1e-3 accurate on analytic quadrics — far too coarse for the
    /// tolerance-driven range walk (`criteria` ~ 1e-7). All other surfaces use
    /// the grid+refine projector.
    pub(crate) fn closest_params_dist(&self, p: &GpPnt) -> (f64, f64, f64) {
        match classify_surface(self.surface()) {
            SurfaceKind::Plane => {
                if let Some((ploc, px, py, pn)) = plane_geometry(self.surface()) {
                    let d = GpVec::from_pnts(&ploc, p);
                    let u = d.dot(&GpVec::from_xyz(px.xyz()));
                    let v = d.dot(&GpVec::from_xyz(py.xyz()));
                    let dist = d.dot(&GpVec::from_xyz(pn.xyz())).abs();
                    return (u, v, dist);
                }
            }
            SurfaceKind::Sphere => {
                if let Some((sc, sr)) = sphere_geometry(self.surface()) {
                    let d = GpVec::from_pnts(&sc, p);
                    let m = d.magnitude();
                    if m > 1e-30 {
                        let u = d.y().atan2(d.x()).rem_euclid(2.0 * PI);
                        let v = (d.z() / m).clamp(-1.0, 1.0).asin();
                        return (u, v, (m - sr).abs());
                    }
                }
            }
            _ => {}
        }
        let (u, v) = surface_closest_params(self.surface(), p, 24, 24);
        let dist = self.surface_d0(u, v).distance(p);
        (u, v, dist)
    }

    // ---- Range expansion ----------------------------------------------------

    /// Expand a result range from a known on-surface point, increasing or
    /// decreasing the parameter. Port of `ComputeRangeFromStartPoint(bool, …)`.
    pub(crate) fn compute_range_from_start_point(&mut self, to_increase: bool, parameter: f64, u: f64, v: f64) {
        let found = self.range_manager.get_index(parameter, to_increase);
        if found < 0 {
            return;
        }
        self.compute_range_from_start_point_idx(to_increase, parameter, u, v, found as usize);
    }

    /// Expand a result range from a known on-surface point, walking the
    /// parameter by `curve_resolution` and bisecting on the distance until the
    /// curve leaves the surface. Port of
    /// `ComputeRangeFromStartPoint(bool, double, double, double, int)`.
    pub(crate) fn compute_range_from_start_point_idx(
        &mut self,
        to_increase: bool,
        parameter: f64,
        u_param: f64,
        v_param: f64,
        mut valid_index: usize,
    ) {
        if self.range_manager.flag(valid_index) > 0 {
            return;
        }
        let mut min_delta = self.curve_resolution * 0.5;
        let mut delta_restrictor = 0.1 * (self.last_parameter - self.first_parameter);
        if min_delta > delta_restrictor {
            min_delta = delta_restrictor * 0.5;
        }
        let ten_of_min_delta = min_delta * 10.0;
        let mut delta = self.curve_resolution;
        let mut cur_par = if to_increase { parameter + delta } else { parameter - delta };
        let mut prev_par = parameter;
        let mut current_range = self.range_manager.range(valid_index);
        let mut boundary_condition = if to_increase {
            cur_par > current_range.last
        } else {
            cur_par < current_range.first
        };
        if boundary_condition {
            cur_par = if to_increase { current_range.last } else { current_range.first };
            boundary_condition = false;
        }
        let mut loop_counter = 0;
        let _ = (u_param, v_param); // initial guess consumed by the local projector
        let mut another_solution_found = false;
        let mut is_boundary_index = false;
        let mut is_valid_index = true;

        while delta >= min_delta && loop_counter <= 10 {
            let a_point = self.curve_d0(cur_par);
            let (su, sv, sd) = self.closest_params_dist(&a_point);
            let point_found = if sd < self.criteria {
                true
            } else {
                self.distance(cur_par) < self.criteria
            };
            let _ = (su, sv);

            if point_found {
                prev_par = cur_par;
                another_solution_found = true;
                if boundary_condition && (is_boundary_index || !is_valid_index) {
                    break;
                }
            } else {
                delta_restrictor = delta;
            }

            delta = if point_found { delta * 2.0 } else { delta * 0.5 };
            delta = if delta < delta_restrictor { delta } else { delta_restrictor };
            cur_par = if to_increase { prev_par + delta } else { prev_par - delta };

            if cur_par == prev_par {
                break;
            }

            boundary_condition = if to_increase {
                cur_par > current_range.last
            } else {
                cur_par < current_range.first
            };
            is_boundary_index = false;
            is_valid_index = true;

            if boundary_condition {
                is_boundary_index = (!to_increase && valid_index == 0)
                    || (to_increase && valid_index + 1 == self.range_manager.len());
                if !is_boundary_index {
                    if point_found {
                        let adj_flag = if to_increase {
                            self.range_manager.flag(valid_index + 1)
                        } else {
                            self.range_manager.flag(valid_index - 1)
                        };
                        if adj_flag == 0 {
                            valid_index = if to_increase { valid_index + 1 } else { valid_index - 1 };
                            current_range = self.range_manager.range(valid_index);
                            if (to_increase && cur_par > current_range.last)
                                || (!to_increase && cur_par < current_range.first)
                            {
                                cur_par = 0.5 * (current_range.first + current_range.last);
                                delta *= 0.5;
                            }
                        } else {
                            is_valid_index = false;
                            cur_par = if to_increase { current_range.last } else { current_range.first };
                        }
                    }
                } else {
                    cur_par = if to_increase { current_range.last } else { current_range.first };
                }
                if delta < ten_of_min_delta {
                    loop_counter += 1;
                } else {
                    loop_counter = 0;
                }
            }
        }

        if another_solution_found {
            if to_increase {
                self.range_manager.insert_range(parameter, prev_par, 2);
            } else {
                self.range_manager.insert_range(prev_par, parameter, 2);
            }
        }
    }

    /// Insert a degenerate (point) result range when a crossing point added no
    /// interval — a tangency point. Port of the static `SetEmptyResultRange`.
    pub(crate) fn set_empty_result_range(&mut self, parameter: f64) {
        let indices = self.range_manager.get_indices(parameter);
        let mut add = !indices.is_empty();
        for &k in &indices {
            if self.range_manager.flag(k) == 2 {
                add = false;
                break;
            }
        }
        if add {
            self.range_manager.insert_range(parameter, parameter, 2);
        }
    }

    // ---- Main algorithm ------------------------------------------------------

    /// Launch the intersection. Sets `is_done` and fills `result()`.
    pub fn perform(&mut self) -> Result<(), String> {
        if self.curve.is_none() || self.surface.is_none() {
            return Err("BeanFaceIntersector::perform: curve and surface must be initialized".into());
        }
        self.is_done = false;
        self.results.clear();
        self.min_sq_distance = f64::MAX;

        // Fast computation of the Line/Plane case.
        if classify_fast_curve(self.curve()) == FastCurveKind::Line
            && classify_surface(self.surface()) == SurfaceKind::Plane
        {
            self.compute_line_plane();
            return Ok(());
        }

        // Fast analytic coincidence / no-intersection verdict.
        if self.fast_compute_analytic() {
            self.is_done = true;
            return Ok(());
        }

        // Range manager over the bean window.
        self.range_manager.set_boundaries(self.first_parameter, self.last_parameter, 0);

        // Coincidence scan over the whole bean.
        if self.test_compute_coinside() {
            self.results.push(IntRange::new_unchecked(self.first_parameter, self.last_parameter));
            self.is_done = true;
            return Ok(());
        }

        // Localized path for high-degree NURBS-like surfaces.
        let b_localize = self.surface_is_localizable();
        let is_localized = b_localize && self.compute_localized();

        if !is_localized {
            self.compute_around_exact_intersection();
            self.compute_using_extremum();
            self.compute_near_range_boundaries();
        }

        self.is_done = true;

        // Collect the flagged (== 2) ranges, merging adjacent ones.
        for i in 0..self.range_manager.len() {
            if self.range_manager.flag(i) != 2 {
                continue;
            }
            let r = self.range_manager.range(i);
            if let Some(last) = self.results.last_mut() {
                if (r.first - last.last).abs() > PCONFUSION {
                    self.results.push(r);
                } else {
                    last.last = last.last.max(r.last);
                }
            } else {
                self.results.push(r);
            }
        }
        Ok(())
    }

    /// Whether the localized path should be attempted for this surface type.
    fn surface_is_localizable(&self) -> bool {
        let finite = self.umin.is_finite()
            && self.umax.is_finite()
            && self.vmin.is_finite()
            && self.vmax.is_finite();
        if !finite {
            return false;
        }
        // OCCT: Bezier / Other / (BSpline with degree>2 and knots>2).
        // Without knot access a BSpline is treated as Other and still localizes.
        match classify_surface_mesh(self.surface()) {
            SurfaceType::BezierSurface | SurfaceType::OtherSurface | SurfaceType::BSplineSurface => {
                true
            }
            _ => false,
        }
    }
}

impl Default for BeanFaceIntersector {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Static helpers
// ---------------------------------------------------------------------------

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
pub(crate) fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    const GOLD: f64 = 0.618_033_988_749_894_9;
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
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpLin, GpPln, GpSphere as GpSphereT};
    use occt_geom::{GeomBSplineCurve, GeomLine, GeomPlane, GeomSphere};

    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_tool::BRepTool;

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    fn sphere_at(origin: GpPnt, r: f64) -> GeomSphere {
        let ax3 = GpAx3::new(origin, dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        GeomSphere::new(GpSphereT::new(ax3, r).unwrap())
    }

    fn unit_sphere() -> Arc<dyn Surface> {
        Arc::new(sphere_at(GpPnt::zero(), 1.0))
    }

    /// A line edge through the box bottom face along +Z at (0.5, 0.5).
    fn z_line(px: f64, py: f64, z0: f64, z1: f64) -> Arc<dyn Curve> {
        let dir = if z1 >= z0 { dir(0.0, 0.0, 1.0) } else { dir(0.0, 0.0, -1.0) };
        Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(px, py, z0), dir)))
    }

    /// The bottom face (z = 0) surface of the unit box.
    fn box_bottom_surface() -> Arc<dyn Surface> {
        let b = unit_box();
        let face = &b.faces[0];
        BRepTool::face_surface(face).expect("box bottom face surface")
    }

    fn run(curve: Arc<dyn Curve>, surface: Arc<dyn Surface>, u: (f64, f64, f64, f64)) -> Vec<IntRange> {
        let mut bfi = BeanFaceIntersector::new();
        let (cf, cl) = (curve.first_parameter(), curve.last_parameter());
        let (cf, cl) = if cf.is_finite() && cl.is_finite() { (cf, cl) } else { (-10.0, 10.0) };
        bfi.initialize(curve, surface, 1e-7, 1e-7);
        bfi.set_bean_parameters(cf, cl);
        bfi.set_surface_parameters(u.0, u.1, u.2, u.3);
        bfi.perform().expect("perform succeeds");
        assert!(bfi.is_done());
        bfi.result()
    }

    #[test]
    fn line_through_box_face_single_range() {
        let surf = box_bottom_surface();
        let curve = z_line(0.5, 0.5, -1.0, 1.0);
        let ranges = run(curve, surf, (0.0, 1.0, 0.0, 1.0));
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        // Line from (0.5,0.5,-1) along +Z: the plane z = 0 is met at t = 1.
        let r = ranges[0];
        assert!(r.first <= 1.0 && r.last >= 1.0, "range {r:?} must bracket the crossing");
        assert!((r.first - 1.0).abs() < 1e-3, "first {r:?}");
        assert!((r.last - 1.0).abs() < 1e-3, "last {r:?}");
    }

    #[test]
    fn line_missing_box_face_empty() {
        let surf = box_bottom_surface();
        // Line outside the face's UV domain (x beyond the box).
        let curve = z_line(5.0, 5.0, -1.0, 1.0);
        let ranges = run(curve, surf, (0.0, 1.0, 0.0, 1.0));
        assert!(ranges.is_empty(), "ranges: {ranges:?}");
    }

    #[test]
    fn line_in_box_plane_covers_bean() {
        let surf = box_bottom_surface();
        // Line lying in z = 0. `ComputeLinePlane` in-plane returns the whole
        // bean range; the face 2D restriction is applied by EdgeFace.
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 0.5, 0.0),
            dir(1.0, 0.0, 0.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surf, 1e-7, 1e-7);
        bfi.set_bean_parameters(-1.0, 2.0);
        bfi.set_surface_parameters(0.0, 1.0, 0.0, 1.0);
        bfi.perform().expect("perform");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        let r = ranges[0];
        assert!((r.first + 1.0).abs() < 1e-9, "first {r:?}");
        assert!((r.last - 2.0).abs() < 1e-9, "last {r:?}");
    }

    #[test]
    fn bspline_crosses_sphere_two_ranges() {
        let surface = unit_sphere();
        // A B-spline through (-2,0,0) and (2,0,0) bowing toward z = 0.5; it
        // enters/exits the unit sphere near x = ±1.
        let poles = vec![
            GpPnt::new(-2.0, 0.0, 0.0),
            GpPnt::new(-0.5, 0.0, 0.6),
            GpPnt::new(0.5, 0.0, 0.6),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let curve = Arc::new(
            GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        let ranges = run(curve, surface, (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0));
        assert!(!ranges.is_empty(), "ranges: {ranges:?}");
        // Every returned parameter must be within tolerance of the sphere.
        let mut bfi = BeanFaceIntersector::new();
        let surface2 = unit_sphere();
        let poles2 = vec![
            GpPnt::new(-2.0, 0.0, 0.0),
            GpPnt::new(-0.5, 0.0, 0.6),
            GpPnt::new(0.5, 0.0, 0.6),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let curve2 = Arc::new(
            GeomBSplineCurve::new(poles2, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        bfi.initialize(curve2, surface2, 1e-7, 1e-7);
        bfi.set_bean_parameters(0.0, 1.0);
        bfi.set_surface_parameters(0.0, 2.0 * PI, -PI / 2.0, PI / 2.0);
        bfi.perform().expect("perform");
        for r in bfi.result() {
            let mut u = 0.0;
            let mut v = 0.0;
            let d = bfi.distance_with_uv(r.first, &mut u, &mut v);
            assert!(d < 1e-4, "range start {r:?} distance {d}");
            let d2 = bfi.distance_with_uv(r.last, &mut u, &mut v);
            assert!(d2 < 1e-4, "range end {r:?} distance {d2}");
        }
        // The minimal squared distance should be (essentially) zero.
        assert!(bfi.minimal_square_distance() < 1e-10, "min_sq {}", bfi.minimal_square_distance());
    }

    #[test]
    fn sphere_missed_by_bspline_is_empty() {
        let surface = unit_sphere();
        // A straight-ish B-spline passing well above the sphere.
        let poles = vec![
            GpPnt::new(-2.0, 0.0, 3.0),
            GpPnt::new(-0.5, 0.0, 3.0),
            GpPnt::new(0.5, 0.0, 3.0),
            GpPnt::new(2.0, 0.0, 3.0),
        ];
        let curve = Arc::new(
            GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        let ranges = run(curve, surface, (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0));
        assert!(ranges.is_empty(), "ranges: {ranges:?}");
    }

    #[test]
    fn minimal_square_distance_fast_reject_stays_large() {
        // A line at y = 3 against the unit sphere: `FastComputeAnalytic`
        // conclusively rejects the pair (no intersection) and returns early,
        // leaving `minimal_square_distance()` at its sentinel — exactly as OCCT
        // leaves `myMinSqDistance` at `RealLast()` on the fast path.
        let surface = unit_sphere();
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 3.0, 0.0),
            dir(1.0, 0.0, 0.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surface, 1e-7, 1e-7);
        bfi.set_bean_parameters(-3.0, 3.0);
        bfi.set_surface_parameters(0.0, 2.0 * PI, -PI / 2.0, PI / 2.0);
        bfi.perform().expect("perform");
        assert!(bfi.result().is_empty(), "line at y=3 misses the unit sphere");
        assert!(bfi.minimal_square_distance() > 1e100, "fast-reject leaves the sentinel");
    }

    #[test]
    fn circle_in_plane_full_range() {
        let ax3 = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        let plane = Arc::new(GeomPlane::new(GpPln::new(ax3))) as Arc<dyn Surface>;
        let circle = Arc::new(occt_geom::GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0)))
            as Arc<dyn Curve>;
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(circle, plane, 1e-7, 1e-7);
        bfi.set_bean_parameters(0.0, 2.0 * PI);
        bfi.set_surface_parameters(-5.0, 5.0, -5.0, 5.0);
        bfi.perform().expect("perform");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        assert!((ranges[0].first - 0.0).abs() < 1e-9);
        assert!((ranges[0].last - 2.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn line_crosses_bspline_surface_localized_path() {
        // A curved BSpline patch (poles at (u, v, u²) → z ≈ u²) over [0,1]²,
        // classified as a NURBS surface so the localized path runs. A vertical
        // line at (0.5, 0.5) crosses it somewhere over the patch; the result
        // must be one range whose midpoint lies on the surface.
        use occt_geom::bspline_surface::{bspline_surface_uniform_knots, GeomBSplineSurface};
        let (ku, kv) = bspline_surface_uniform_knots(3, 3, 2, 2);
        let poles: Vec<Vec<GpPnt>> = (0..3)
            .map(|i| {
                let u = i as f64 * 0.5;
                (0..3)
                    .map(|j| {
                        let v = j as f64 * 0.5;
                        GpPnt::new(u, v, u * u)
                    })
                    .collect()
            })
            .collect();
        let surface = Arc::new(GeomBSplineSurface::new(poles, ku, kv, 2, 2).unwrap()) as Arc<dyn Surface>;
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.5, 0.5, -1.0),
            dir(0.0, 0.0, 1.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surface, 1e-7, 1e-7);
        bfi.set_bean_parameters(-1.0, 3.0);
        bfi.set_surface_parameters(0.0, 1.0, 0.0, 1.0);
        bfi.perform().expect("perform succeeds on BSpline surface");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        let r = ranges[0];
        // The crossing point must genuinely lie on the surface.
        let mid = 0.5 * (r.first + r.last);
        let p = GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.5, 0.5, -1.0), dir(0.0, 0.0, 1.0)))
            .d0(mid);
        let (_, _, d) = bfi.closest_params_dist(&p);
        assert!(d < 1e-4, "crossing at {mid} is {d} off the BSpline patch");
    }

    #[test]
    fn edge_face_cross_check_matches_approximate_solver() {
        // Line through the box bottom face: compare with the existing
        // approximate `inttools::edge_face_intersections`.
        let bx = unit_box();
        let face = &bx.faces[0];
        let b = crate::builder::TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        let approx = crate::inttools::edge_face_intersections(&edge, face, 1e-9);
        assert_eq!(approx.len(), 1, "approx: {approx:?}");
        assert!(approx[0].1.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6);

        let surf = BRepTool::face_surface(face).expect("face surface");
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.5, 0.5, -1.0),
            dir(0.0, 0.0, 1.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surf, 1e-7, 1e-7);
        bfi.set_bean_parameters(-10.0, 10.0);
        bfi.set_surface_parameters(0.0, 1.0, 0.0, 1.0);
        bfi.perform().expect("perform");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        // Approximate solver hits at edge parameter 1.0 (z = 0); the bean's
        // crossing is at curve parameter 1.0 too. Compare the 3D points.
        let mid = 0.5 * (ranges[0].first + ranges[0].last);
        let p = GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.5, 0.5, -1.0), dir(0.0, 0.0, 1.0)))
            .d0(mid);
        assert!(p.distance(&approx[0].1) < 1e-3, "bean midpoint {p:?} vs approx {:?}", approx[0].1);
    }
}
