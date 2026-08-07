//! Face–face surface intersection — thin orchestrator over `crate::intpatch`.
//!
//! Ports the orchestration of `IntTools_FaceFace` (TKBO) without the IntPatch
//! walker and `GeomInt_WLApprox`: the two faces' underlying surfaces are
//! classified, supported analytic pairs dispatch to the exact closed forms in
//! `crate::intpatch` (plane∩plane, plane∩sphere, sphere∩sphere, plane∩cylinder)
//! and `occt_geom::intana` (plane∩cone, plane∩torus), and everything else falls
//! back to the intpatch sampling tracer.
//!
//! The output mirrors `IntTools_Curve`: a 3D curve plus the per-face 2D
//! pcurves (`pcurve_full::make_pcurve_full`) over a valid parameter range.
//!
//! ponytail: cone×cone / cone×sphere / torus×* (non-plane) pairs still fall
//! through to the grid tracer — OCCT uses the numeric IntPatch walker there, so
//! the boundary is the same; the plane×cylinder generatrix pair collapses to a
//! single line curve. Add per-pair closed forms when a caller needs them.

use std::cmp::Ordering;
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpDir, GpDir2d, GpLin, GpLin2d, GpPln, GpPnt, GpPnt2d,
    GpVec, GpVec2d,
};
use occt_geom::intana::{quadric_quadric_plane_cone, QuadricIntersection};
use occt_geom::{
    Curve, GeomCircle, GeomEllipse, GeomHyperbola, GeomLine, GeomParabola, Surface,
};
use occt_geom2d::curve::Curve2d;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::intpatch::{self, IntersectionCurve, SurfaceIntersection};
use crate::inttools_data::{CurveKind, IntRange};
use crate::pcurve_full::{classify_surface_kind, cone_params, make_pcurve_full, torus_params};
use crate::shape::Face;
use crate::tgeometry::GeometryRegistry;

/// Default half-extent used when an intersection curve reports an unbounded
/// parameter range, so it can still be represented in a `FaceFaceCurve`.
const DEFAULT_WINDOW: f64 = 4.0;

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

// ---------------------------------------------------------------------------
// Surface / curve helpers
// ---------------------------------------------------------------------------

/// OCCT `IntTools_FaceFace::IndexType` — a total ordering used to sort the two
/// faces so the "higher" analytic type becomes `Face1`.
fn surface_sort_index(kind: SurfaceKind) -> usize {
    match kind {
        SurfaceKind::Plane => 0,
        SurfaceKind::Cylinder => 1,
        SurfaceKind::Cone => 2,
        SurfaceKind::Sphere => 3,
        SurfaceKind::Torus => 4,
        SurfaceKind::Other => 10,
    }
}

/// Cylinder parameters `(center, unit axis, radius)` recovered from a surface
/// by sampling its invariants (mirrors `pcurve_full::cylinder_params`).
fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axis.divided(m);
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let r = q0.distance(&q1) / 2.0;
    if r <= 1e-9 {
        return None;
    }
    let center = mid(&q0, &q1);
    let (nu, nv) = (8, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = (j as f64 - nv as f64 / 2.0) / 2.0;
            let p = s.d0(u, v);
            if (dist_to_axis(&p, &center, &ax) - r).abs() > 1e-3 * r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, r))
}

fn mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

/// The analytic family of a plane×cylinder intersection, recomputed from the
/// plane normal and the cylinder axis (matches `intersect_plane_cylinder`).
fn plane_cylinder_kind(pln: &GpPln, ax: &GpAx1) -> CurveKind {
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = GpVec::from_xyz(ax.direction().xyz()).normalized();
    let n_par = n.dot(&az);
    let n_perp = n.subtracted(&az.multiplied_scalar(n_par));
    if n_perp.magnitude() <= 1e-9 {
        CurveKind::Circle
    } else if n_par.abs() <= 1e-9 {
        CurveKind::Line
    } else {
        CurveKind::Ellipse
    }
}

// ---------------------------------------------------------------------------
// Cone / Torus helpers
// ---------------------------------------------------------------------------

/// A unit direction perpendicular to `z` (axis for a cone/torus frame).
fn perp_x_dir(z: &GpDir) -> GpDir {
    let base = if z.x().abs() < 0.9 {
        GpDir::new(1.0, 0.0, 0.0).expect("x")
    } else {
        GpDir::new(0.0, 1.0, 0.0).expect("y")
    };
    base.cross(z).unwrap_or(base)
}

/// Recover a `GpCone` from a cone surface by sampling its invariants
/// (apex + axis + semi-angle from `pcurve_full`). The cone is anchored at its
/// apex with radius 0 — `IntAna_QuadQuadGeo` only reads apex / axis / semi-angle
/// geometry (plus `radius` in one axis-flip heuristic, where 0 is safe).
fn cone_from_surface(s: &dyn Surface) -> Option<GpCone> {
    let (apex, ax, alpha) = cone_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(apex, z, &x).ok()?;
    GpCone::new(ax3, 0.0, alpha).ok()
}

/// Plane coefficients `A·x + B·y + C·z + D = 0` with unit normal (OCCT form).
fn plane_coeffs(p: &GpPln) -> (f64, f64, f64, f64) {
    let n = GpVec::from_xyz(p.axis().direction().xyz()).normalized();
    let loc = p.location();
    let d = -(n.x() * loc.x() + n.y() * loc.y() + n.z() * loc.z());
    (n.x(), n.y(), n.z(), d)
}

/// A `GpCirc` in the plane through `center` with normal `normal` and `radius`.
fn circle_gp(center: GpPnt, normal: GpDir, radius: f64) -> Option<GpCirc> {
    let x_dir = perp_x_dir(&normal);
    let ax2 = GpAx2::new(center, normal, x_dir).ok()?;
    Some(GpCirc::new(ax2, radius))
}

/// Plane ∩ torus circles. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln,
/// gp_Torus)` (IntAna_QuadQuadGeo.cxx): up to two circles when the torus axis
/// is parallel to the plane normal (perpendicular cut → radii `major ± dt`,
/// `dt = √(minor² − dist²)` from the center) or perpendicular to it (axis in
/// the plane through the center → two `minor` circles at `±major`). Returns
/// `None` for any other orientation — the general (non-planar) torus section,
/// which OCCT routes to the numeric walker.
fn plane_torus_circles(pln: &GpPln, center: GpPnt, ax: &GpVec, major: f64, minor: f64, tol: f64) -> Option<Vec<GpCirc>> {
    if minor >= major {
        return None; // degenerate torus → IntAna_NoGeometricSolution
    }
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = ax.normalized();
    let n_par = n.dot(&az);
    let (a, b, c, d) = plane_coeffs(pln);
    let dist = a * center.x() + b * center.y() + c * center.z() + d;

    if (n_par.abs() - 1.0).abs() <= 1e-12 {
        // Axis ∥ plane normal → perpendicular cut.
        let a_dr = dist.abs() - minor;
        if a_dr > 1e-13 {
            return None; // plane misses the tube → IntAna_Empty
        }
        let dist = if a_dr.abs() < 1e-13 {
            if dist < 0.0 { -minor } else { minor }
        } else {
            dist
        };
        let a_dt = (minor * minor - dist * dist).max(0.0).sqrt();
        let center_on_plane = center.translated_vec(&n.multiplied_scalar(-dist));
        let normal = GpDir::from_vec(&n).ok()?;
        let mut out = vec![circle_gp(center_on_plane, normal, major + a_dt)?];
        if a_dr < -1e-13 && a_dt > tol {
            out.push(circle_gp(center_on_plane, normal, (major - a_dt).max(0.0))?);
        }
        return Some(out);
    }

    if n_par.abs() > 1e-12 {
        return None; // oblique → IntAna_NoGeometricSolution (numeric)
    }
    // Axis ⊥ normal → plane must contain the torus axis through the center.
    if dist.abs() > 1e-14 {
        return None;
    }
    let a_dir = GpDir::from_vec(&az).ok()?;
    let normal = GpDir::from_vec(&n).ok()?;
    let e = a_dir.cross(&normal).ok()?;
    let c1 = center.translated_vec(&GpVec::from_xyz(e.xyz()).multiplied_scalar(major));
    let c2 = center.translated_vec(&GpVec::from_xyz(e.xyz()).multiplied_scalar(-major));
    Some(vec![
        circle_gp(c1, normal, minor)?,
        circle_gp(c2, normal, minor)?,
    ])
}

/// Parameter interval of the (unit-speed) 2D line that lies inside the UV
/// rectangle — a Liang–Barsky clip. Returns `None` when the line misses the
/// rectangle entirely. Unbounded rectangle dimensions impose no constraint.
fn line_in_uv_rect(lin: &GpLin2d, bounds: (f64, f64, f64, f64)) -> Option<(f64, f64)> {
    let (umin, umax, vmin, vmax) = bounds;
    let loc = lin.pos.loc;
    let dir = lin.pos.vdir;
    let (x0, y0) = (loc.x(), loc.y());
    let (dx, dy) = (dir.x, dir.y);
    let mut t0 = f64::NEG_INFINITY;
    let mut t1 = f64::INFINITY;
    for (lo, hi, p0, dp) in [(umin, umax, x0, dx), (vmin, vmax, y0, dy)] {
        if !lo.is_finite() || !hi.is_finite() {
            continue;
        }
        if dp.abs() < 1e-15 {
            if p0 < lo - 1e-9 || p0 > hi + 1e-9 {
                return None;
            }
            continue;
        }
        let ta = (lo - p0) / dp;
        let tb = (hi - p0) / dp;
        let (e, l) = if ta < tb { (ta, tb) } else { (tb, ta) };
        t0 = t0.max(e);
        t1 = t1.min(l);
        if t0 > t1 + 1e-9 {
            return None;
        }
    }
    if !t0.is_finite() || !t1.is_finite() || t1 - t0 <= 1e-9 {
        return None;
    }
    Some((t0, t1))
}

/// Unit-speed 2D line through two distinct points (`None` when coincident).
fn lin2d_through(q0: GpPnt2d, q1: GpPnt2d) -> Option<GpLin2d> {
    let v = GpVec2d::new(q1.x() - q0.x(), q1.y() - q0.y());
    let m = v.magnitude();
    if m < 1e-12 {
        return None;
    }
    let d = GpDir2d::from_vec2d(&v).ok()?;
    Some(GpLin2d::from_pnt_dir(q0, d))
}

/// The finite range of a curve, or a default window for unbounded curves.
fn curve_range(curve: &dyn Curve) -> IntRange {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        IntRange::new_unchecked(a, b)
    } else {
        IntRange::new_unchecked(-DEFAULT_WINDOW, DEFAULT_WINDOW)
    }
}

/// Pcurve of `curve` on `face` over `range`, via a temporary edge registered in
/// the geometry side-table (cleaned up immediately after). Best-effort: returns
/// `None` when the pcurve cannot be computed.
fn pcurve_of_curve(curve: &Arc<dyn Curve>, range: IntRange, face: &Face) -> Option<Arc<dyn Curve2d>> {
    if !range.is_valid() || range.length() < 1e-12 {
        return None;
    }
    let builder = TopoBuilder::new();
    let edge = builder.make_edge(curve.clone(), range.first, range.last);
    let pc = make_pcurve_full(&edge, face).ok();
    GeometryRegistry::global().clear_shape(&edge.0);
    pc
}

/// Whether two curves represent the same intersection (same kind, same range
/// and a coincident mid-point) — used to deduplicate the result.
fn same_curve(a: &FaceFaceCurve, b: &FaceFaceCurve) -> bool {
    if a.kind != b.kind {
        return false;
    }
    if (a.range.first - b.range.first).abs() > 1e-7 || (a.range.last - b.range.last).abs() > 1e-7 {
        return false;
    }
    let ta = 0.5 * (a.range.first + a.range.last);
    let tb = 0.5 * (b.range.first + b.range.last);
    a.curve.d0(ta).distance(&b.curve.d0(tb)) < 1e-6
}

// ---------------------------------------------------------------------------
// FaceFace
// ---------------------------------------------------------------------------

/// Intersects the underlying surfaces of two faces.
///
/// Mirrors `IntTools_FaceFace`: `SetParameters` is dropped (approximation is
/// always on), `Perform` classifies the surfaces, dispatches to the analytic
/// `intpatch` closed forms when possible and the sampling tracer otherwise,
/// and stores the resulting curves in a [`FaceFaceResult`].
pub struct FaceFace {
    face1: Option<Face>,
    face2: Option<Face>,
    tol: f64,
    done: bool,
    tangent_faces: bool,
    result: FaceFaceResult,
}

impl FaceFace {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            face1: None,
            face2: None,
            tol: DEFAULT_TOL,
            done: false,
            tangent_faces: false,
            result: FaceFaceResult::new(),
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

    /// Sets the intersection context (accepted for API compatibility; the thin
    /// port does not use it).
    pub fn set_context(&mut self, _ctx: FaceFaceContext) {
        // no-op: the context caches IntPatch walker state we do not have
    }

    /// Intersects the underlying surfaces of the two faces.
    ///
    /// The faces are ordered internally so the "higher" analytic type becomes
    /// `face1` (OCCT `SortTypes`); the result curves always carry `face1_idx =
    /// 0` and `face2_idx = 1` in that sorted order. Output curves are sorted by
    /// curve parameter and deduplicated.
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

        // Output sorted by curve parameter, then deduplicated.
        curves.sort_by(|x, y| x.range.first.partial_cmp(&y.range.first).unwrap_or(Ordering::Equal));
        curves.dedup_by(|a, b| same_curve(a, b));

        self.result = FaceFaceResult { curves };
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
            _ => self.general(sa, sb, fa, fb, tol),
        }
    }

    /// Plane ∩ plane — an exact intersection line trimmed to the overlap of
    /// the two faces' UV domains (mirrors `PerformPlanes`).
    fn plane_plane(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let pa = intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction (a)")?;
        let pb = intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction (b)")?;
        let n1 = GpVec::from_xyz(pa.axis().direction().xyz());
        let n2 = GpVec::from_xyz(pb.axis().direction().xyz());
        if n1.cross_magnitude(&n2) <= tol {
            // Parallel planes: coincident → tangent, otherwise no intersection.
            let dd = n1
                .normalized()
                .dot(&GpVec::from_pnts(&pa.location(), &pb.location()))
                .abs();
            if dd <= tol {
                self.tangent_faces = true;
            }
            return Ok(Vec::new());
        }
        let (origin, dir) = crate::face_face::plane_plane_intersection(&pa, &pb)
            .ok_or("FaceFace: plane-plane intersection failed")?;
        match self.trim_plane_plane_line(origin, dir, sa, sb, fa, fb, tol)? {
            Some(c) => Ok(vec![c]),
            None => Ok(Vec::new()),
        }
    }

    /// Trim the infinite plane/plane intersection line to the segment lying on
    /// both bounded faces, by projecting the line into each face's UV domain
    /// and clipping against the face UV bounds.
    fn trim_plane_plane_line(
        &self,
        origin: GpPnt,
        dir: GpVec,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Option<FaceFaceCurve>, String> {
        let dir_u = dir.normalized();
        let p0 = origin;
        let p1 = origin.translated_vec(&dir_u);
        let proj_a0 = intpatch::project_params(sa, &p0);
        let proj_a1 = intpatch::project_params(sa, &p1);
        let proj_b0 = intpatch::project_params(sb, &p0);
        let proj_b1 = intpatch::project_params(sb, &p1);
        let a0 = GpPnt2d::new(proj_a0.0, proj_a0.1);
        let a1 = GpPnt2d::new(proj_a1.0, proj_a1.1);
        let b0 = GpPnt2d::new(proj_b0.0, proj_b0.1);
        let b1 = GpPnt2d::new(proj_b1.0, proj_b1.1);
        let ba = crate::wireframe::face_uv_bounds(fa, sa);
        let bb = crate::wireframe::face_uv_bounds(fb, sb);

        // For each face, clip the projected line against the face UV rectangle.
        // The 3D line is p(t) = origin + t·dir (unit dir); the projected 2D
        // line maps p(t) → q0 + t·(q1 − q0), so a 2D parameter interval
        // [s0, s1] corresponds to the 3D interval [s0/k, s1/k], k = |q1 − q0|.
        let mut t_lo = f64::NEG_INFINITY;
        let mut t_hi = f64::INFINITY;
        for (q0, q1, bounds) in [(a0, a1, ba), (b0, b1, bb)] {
            let Some(lin) = lin2d_through(q0, q1) else { continue };
            let k = q0.distance(&q1);
            if k < 1e-12 {
                continue;
            }
            match line_in_uv_rect(&lin, bounds) {
                Some((s0, s1)) => {
                    t_lo = t_lo.max(s0 / k);
                    t_hi = t_hi.min(s1 / k);
                }
                None => return Ok(None),
            }
        }
        if !t_lo.is_finite() || !t_hi.is_finite() || t_hi - t_lo <= tol {
            return Ok(None);
        }

        let lin = GpLin::from_pnt_dir(origin, GpDir::from_vec(&dir_u).map_err(|e| e.to_string())?);
        let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(lin));
        let range = IntRange::new_unchecked(t_lo, t_hi);
        Ok(Some(self.finish_curve(curve, CurveKind::Line, range, fa, fb)))
    }

    /// Plane ∩ sphere — an exact circle (or no curve when the plane misses the
    /// sphere), falling back to the tracer when the closed form is empty.
    fn plane_sphere(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, (c, r)) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                intpatch::sphere_params(sb).ok_or("FaceFace: sphere extraction")?,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                intpatch::sphere_params(sa).ok_or("FaceFace: sphere extraction")?,
            )
        };
        let mut out = Vec::new();
        if let Some(ic) = intpatch::intersect_plane_sphere(&pln, c, r) {
            out.push(self.curve_from_ic(ic, CurveKind::Circle, fa, fb));
        }
        if out.is_empty() {
            out = self.general(sa, sb, fa, fb, tol)?;
        }
        Ok(out)
    }

    /// Sphere ∩ sphere — an exact circle (or no curve when disjoint/nested),
    /// falling back to the tracer when the closed form is empty.
    fn sphere_sphere(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let (c1, r1) = intpatch::sphere_params(sa).ok_or("FaceFace: sphere extraction (a)")?;
        let (c2, r2) = intpatch::sphere_params(sb).ok_or("FaceFace: sphere extraction (b)")?;
        let mut out = Vec::new();
        if let Some(ic) = intpatch::intersect_sphere_sphere(c1, r1, c2, r2) {
            out.push(self.curve_from_ic(ic, CurveKind::Circle, fa, fb));
        }
        if out.is_empty() {
            out = self.general(sa, sb, fa, fb, tol)?;
        }
        Ok(out)
    }

    /// Plane ∩ cylinder — an exact circle / generatrix line / ellipse, falling
    /// back to the tracer when the closed form is empty.
    fn plane_cylinder(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, cyl) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                sb,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                sa,
            )
        };
        let (center, ax, rad) = cylinder_params(cyl).ok_or("FaceFace: cylinder extraction")?;
        let ax1 = GpAx1::new(center, GpDir::from_vec(&ax).map_err(|e| e.to_string())?);
        let kind = plane_cylinder_kind(&pln, &ax1);
        let mut out = Vec::new();
        if let Some(ic) = intpatch::intersect_plane_cylinder(&pln, &ax1, rad) {
            out.push(self.curve_from_ic(ic, kind, fa, fb));
        }
        if out.is_empty() {
            out = self.general(sa, sb, fa, fb, tol)?;
        }
        Ok(out)
    }

    /// Plane ∩ cone — the exact conic section (circle / ellipse / parabola /
    /// hyperbola / two generatrix lines) via `IntAna_QuadQuadGeo::Perform
    /// (gp_Pln, gp_Cone)` (intana::quadric_quadric_plane_cone), falling back to
    /// the tracer when the closed form yields nothing (tangent apex point,
    /// degenerate or disjoint).
    fn plane_cone(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, cone) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                cone_from_surface(sb).ok_or("FaceFace: cone extraction")?,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                cone_from_surface(sa).ok_or("FaceFace: cone extraction")?,
            )
        };
        let qi = quadric_quadric_plane_cone(&pln, &cone, 1e-12, 1e-7);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Plane ∩ torus — the exact circles of `IntAna_QuadQuadGeo::Perform
    /// (gp_Pln, gp_Torus)` (up to two: a perpendicular cut at `major ± dt`, or
    /// the two `minor` circles when the axis lies in the plane), falling back
    /// to the tracer for an oblique section.
    fn plane_torus(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, (center, ax, major, minor)) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                torus_params(sb).ok_or("FaceFace: torus extraction")?,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                torus_params(sa).ok_or("FaceFace: torus extraction")?,
            )
        };
        let mut out = Vec::new();
        if let Some(circs) = plane_torus_circles(&pln, center, &ax, major, minor, tol) {
            for c in circs {
                out.push(self.conic_curve(Arc::new(GeomCircle::new(c)), CurveKind::Circle, fa, fb));
            }
        }
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Wrap an `intana` conic (analytic section) into a `FaceFaceCurve`.
    fn conic_curve(
        &self,
        curve: Arc<dyn Curve>,
        kind: CurveKind,
        fa: &Face,
        fb: &Face,
    ) -> FaceFaceCurve {
        let range = curve_range(curve.as_ref());
        self.finish_curve(curve, kind, range, fa, fb)
    }

    /// Convert a `QuadricIntersection` (exact conic section) into section
    /// curves. A tangent `Point`, `Same` and `None` produce no curve — the
    /// caller falls back to the tracer.
    fn conics_to_curves(&self, qi: QuadricIntersection, fa: &Face, fb: &Face) -> Vec<FaceFaceCurve> {
        use QuadricIntersection::*;
        match qi {
            Line(l) => vec![self.conic_curve(Arc::new(GeomLine::new(l)), CurveKind::Line, fa, fb)],
            TwoLines(l1, l2) => vec![
                self.conic_curve(Arc::new(GeomLine::new(l1)), CurveKind::Line, fa, fb),
                self.conic_curve(Arc::new(GeomLine::new(l2)), CurveKind::Line, fa, fb),
            ],
            Circle(c) => vec![self.conic_curve(Arc::new(GeomCircle::new(c)), CurveKind::Circle, fa, fb)],
            Ellipse(e) => vec![self.conic_curve(Arc::new(GeomEllipse::new(e)), CurveKind::Ellipse, fa, fb)],
            Parabola(p) => vec![self.conic_curve(Arc::new(GeomParabola::new(p)), CurveKind::Parabola, fa, fb)],
            Hyperbola(h) => vec![self.conic_curve(Arc::new(GeomHyperbola::new(h)), CurveKind::Hyperbola, fa, fb)],
            Point(_) | Same | None => Vec::new(),
        }
    }

    /// General fallback: the `intpatch` sampling tracer (analytic dispatcher +
    /// marching-squares grid tracer). Marks `tangent_faces` on coincidence.
    fn general(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let mut out = Vec::new();
        match intpatch::surface_surface_intersection(sa, sb, tol) {
            SurfaceIntersection::Curves(ics) => {
                for ic in ics {
                    out.push(self.curve_from_ic(ic, CurveKind::BSpline, fa, fb));
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
        }
    }
}

impl Default for FaceFace {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpAx3, GpCylinder, GpSphere, GpTorus};
    use occt_geom::{GeomCone, GeomCylinder, GeomSphere, GeomTorus};

    const TOL: f64 = 1e-6;

    fn plane_z(z: f64) -> GpPln {
        GpPln::new(
            GpAx3::new(
                GpPnt::new(0.0, 0.0, z),
                GpDir::new(0.0, 0.0, 1.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
        )
    }

    fn plane_face(pln: &GpPln) -> Face {
        TopoBuilder::new().make_face_plane(pln)
    }

    fn sphere_face(center: GpPnt, r: f64) -> Face {
        let ax3 = GpAx3::new(
            center,
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        TopoBuilder::new().make_face(Arc::new(GeomSphere::new(GpSphere::new(ax3, r).unwrap())), &[])
    }

    fn cylinder_face(radius: f64) -> Face {
        let ax3 = GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        TopoBuilder::new().make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax3, radius).unwrap())), &[])
    }

    /// A cone face, apex at the origin, axis +Z, with the given semi-angle
    /// (location ring radius 1 — the geometric tip is the location).
    fn cone_face(semi_angle: f64) -> Face {
        let cone = GpCone::new(GpAx3::standard(), 1.0, semi_angle).unwrap();
        TopoBuilder::new().make_face(Arc::new(GeomCone::new(cone)), &[])
    }

    /// A torus face centered at the origin, axis +Z.
    fn torus_face(major: f64, minor: f64) -> Face {
        let torus = GpTorus::new(GpAx3::standard(), major, minor).unwrap();
        TopoBuilder::new().make_face(Arc::new(GeomTorus::new(torus)), &[])
    }

    /// A plane face through `origin` with the given unit normal (x-dir is
    /// `+X`, valid for normals with zero x-component).
    fn plane_face_normal(origin: GpPnt, normal: GpDir) -> Face {
        let ax3 = GpAx3::new(origin, normal, &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        TopoBuilder::new().make_face_plane(&GpPln::new(ax3))
    }

    /// Sample a FaceFaceCurve and assert every sample lies on both faces'
    /// surfaces (within `eps`).
    fn assert_points_on_both(curve: &FaceFaceCurve, sa: &dyn Surface, sb: &dyn Surface, eps: f64) {
        let (a, b) = (curve.range.first, curve.range.last);
        assert!(a.is_finite() && b.is_finite(), "range {a}..{b}");
        for i in 0..=16 {
            let t = a + (b - a) * i as f64 / 16.0;
            let p = curve.curve.d0(t);
            assert!(intpatch::distance_to_surface(&p, sa) < eps, "on surface a: {p:?}");
            assert!(intpatch::distance_to_surface(&p, sb) < eps, "on surface b: {p:?}");
        }
    }

    #[test]
    fn not_done_before_perform() {
        let ff = FaceFace::new();
        assert!(!ff.is_done());
        assert!(ff.result().is_empty());
    }

    #[test]
    fn perform_requires_faces() {
        let mut ff = FaceFace::new();
        ff.set_face1(plane_face(&plane_z(0.0)));
        assert!(ff.perform().is_err(), "missing face2 must fail");
    }

    #[test]
    fn plane_sphere_intersects_in_circle() {
        let pln = plane_face(&plane_z(0.5));
        let sph = sphere_face(GpPnt::zero(), 1.0);
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(sph);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        assert!(ff.is_done());
        let res = ff.result();
        assert_eq!(res.nb_curves(), 1, "one circle of intersection");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Circle, "plane∩sphere is a circle");
        let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(plane_z(0.5)));
        let s_sph: Arc<dyn Surface> = Arc::new(GeomSphere::new(
            GpSphere::new(
                GpAx3::new(
                    GpPnt::zero(),
                    GpDir::new(0.0, 0.0, 1.0).unwrap(),
                    &GpDir::new(1.0, 0.0, 0.0).unwrap(),
                )
                .unwrap(),
                1.0,
            )
            .unwrap(),
        ));
        for i in 0..=16 {
            let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
            let p = c.curve.d0(t);
            assert!((p.distance(&GpPnt::zero()) - 1.0).abs() < 1e-6, "on sphere: {p:?}");
            assert!((p.z() - 0.5).abs() < 1e-6, "in plane: {p:?}");
        }
        assert_points_on_both(c, s_pln.as_ref(), s_sph.as_ref(), 1e-4);
    }

    #[test]
    fn sphere_sphere_intersects_in_circle() {
        let s1 = sphere_face(GpPnt::zero(), 1.0);
        let s2 = sphere_face(GpPnt::new(1.5, 0.0, 0.0), 1.0);
        let mut ff = FaceFace::new();
        ff.set_face1(s1);
        ff.set_face2(s2);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 1, "one circle of intersection");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Circle);
        let expected_r = (1.0f64 - 0.75 * 0.75).sqrt();
        for i in 0..=16 {
            let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
            let p = c.curve.d0(t);
            assert!((p.x() - 0.75).abs() < 1e-6, "plane x: {p:?}");
            let r = GpPnt::new(0.0, p.y(), p.z()).distance(&GpPnt::zero());
            assert!((r - expected_r).abs() < 1e-6, "radius {r} (expected {expected_r})");
        }
    }

    #[test]
    fn two_box_faces_intersect_in_line_segment() {
        let b = crate::brep_extrema::test_box::unit_box();
        let bottom = b.faces[0].clone(); // z = 0
        let front = b.faces[2].clone(); // y = 0
        let mut ff = FaceFace::new();
        ff.set_face1(bottom);
        ff.set_face2(front);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 1, "adjacent box faces meet in one line");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Line);
        // The line is the box edge (0,0,0)-(1,0,0).
        let (a, b) = (c.range.first, c.range.last);
        let p0 = c.curve.d0(a);
        let p1 = c.curve.d0(b);
        assert!((p0.y()).abs() < 1e-6 && (p0.z()).abs() < 1e-6, "start {p0:?}");
        assert!((p1.y()).abs() < 1e-6 && (p1.z()).abs() < 1e-6, "end {p1:?}");
        assert!((p1.x() - p0.x()).abs() > 0.5, "segment spans the box edge");
        assert!(p0.x() >= -1e-6 && p1.x() <= 1.0 + 1e-6, "within x∈[0,1]");
    }

    #[test]
    fn plane_cylinder_intersects_in_circle() {
        let pln = plane_face(&plane_z(1.0));
        let cyl = cylinder_face(2.0);
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(cyl);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 1, "plane ⊥ cylinder axis meets in one circle");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Circle);
        for i in 0..=16 {
            let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
            let p = c.curve.d0(t);
            assert!((p.z() - 1.0).abs() < 1e-6, "in plane: {p:?}");
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - 2.0).abs() < 1e-6, "radius {r}");
        }
    }

    #[test]
    fn plane_cylinder_parallel_generatrices() {
        // Plane x = 0.5 parallel to the Z-axis cylinder radius 1 → generatrix
        // lines at y = ±√0.75. The FaceFace curve must agree with the intpatch
        // closed form it delegates to (对拍).
        let pln = GpPln::new(
            GpAx3::new(
                GpPnt::new(0.5, 0.0, 0.0),
                GpDir::new(1.0, 0.0, 0.0).unwrap(),
                &GpDir::new(0.0, 1.0, 0.0).unwrap(),
            )
            .unwrap(),
        );
        let ax1 = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ref_ic = intpatch::intersect_plane_cylinder(&pln, &ax1, 1.0).expect("reference lines");
        let ref_mid = ref_ic.curve.d0(0.0);

        let mut ff = FaceFace::new();
        ff.set_face1(plane_face(&pln));
        ff.set_face2(cylinder_face(1.0));
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert!(res.nb_curves() >= 1, "at least one generatrix line");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Line);
        let (a, b) = (c.range.first, c.range.last);
        let mid = c.curve.d0(0.5 * (a + b));
        assert!(mid.distance(&ref_mid) < 1e-6, "FaceFace line {mid:?} != intpatch {ref_mid:?}");
    }

    #[test]
    fn plane_cone_intersects_in_circle() {
        // Cone with semi-angle atan(0.5), geometric tip at the origin; plane
        // z = 1 perpendicular to the axis cuts a circle of radius tan(atan0.5)
        // = 0.5 at (0,0,1).
        let pln = plane_face(&plane_z(1.0));
        let cone = cone_face(0.5f64.atan());
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(cone);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 1, "plane ⊥ cone axis meets in one circle");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Circle);
        for i in 0..=16 {
            let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
            let p = c.curve.d0(t);
            assert!((p.z() - 1.0).abs() < 1e-6, "in plane: {p:?}");
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - 0.5).abs() < 1e-6, "radius {r}");
        }
        let s_cone: Arc<dyn Surface> = Arc::new(GeomCone::new(GpCone::new(GpAx3::standard(), 1.0, 0.5f64.atan()).unwrap()));
        let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(plane_z(1.0)));
        assert_points_on_both(c, s_pln.as_ref(), s_cone.as_ref(), 1e-4);
    }

    #[test]
    fn plane_cone_intersects_in_ellipse() {
        // Oblique plane (not through the apex, not perpendicular to the axis,
        // not parallel to a generatrix) → exact ellipse, sampled on both faces.
        let pln = plane_face_normal(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.3, 0.954).unwrap());
        let cone = cone_face(0.5f64.atan());
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(cone);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 1, "one conic section");
        let c = res.curve(0);
        assert_eq!(c.kind, CurveKind::Ellipse, "oblique cut is an ellipse");
        let s_cone: Arc<dyn Surface> = Arc::new(GeomCone::new(GpCone::new(GpAx3::standard(), 1.0, 0.5f64.atan()).unwrap()));
        let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(GpPln::new(
            GpAx3::new(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.3, 0.954).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(),
        )));
        assert_points_on_both(c, s_pln.as_ref(), s_cone.as_ref(), 1e-4);
    }

    #[test]
    fn plane_torus_axis_in_plane_two_circles() {
        // Torus R=3 r=1, plane y = 0 contains the axis → two circles of radius
        // 1 centered at (±3, 0, 0).
        let pln = plane_face_normal(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let tor = torus_face(3.0, 1.0);
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(tor);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 2, "axis-in-plane cut gives two minor circles");
        for c in res.curves() {
            assert_eq!(c.kind, CurveKind::Circle);
            // The circle's center (midpoint of two opposite points) is on the
            // major ring: (±3, 0, 0).
            let (a, b) = (c.range.first, c.range.last);
            let p0 = c.curve.d0(a);
            let p1 = c.curve.d0(0.5 * (a + b)); // opposite point (circle is 2π-periodic)
            let center = mid(&p0, &p1);
            assert!(center.y().abs() < 1e-6, "center in plane: {center:?}");
            assert!((GpPnt::new(center.x(), 0.0, center.z()).distance(&GpPnt::zero()) - 3.0).abs() < 1e-6,
                "center on major ring {center:?}");
            for i in 0..=16 {
                let t = a + (b - a) * i as f64 / 16.0;
                let p = c.curve.d0(t);
                assert!(p.y().abs() < 1e-6, "in plane: {p:?}");
            }
        }
        let s_tor: Arc<dyn Surface> = Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap()));
        let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(GpPln::new(
            GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(),
        )));
        for c in res.curves() {
            assert_points_on_both(c, s_pln.as_ref(), s_tor.as_ref(), 1e-4);
        }
    }

    #[test]
    fn plane_torus_perpendicular_cut_two_circles() {
        // Torus R=3 r=1, plane z = 0 perpendicular to the axis through the
        // center → two circles of radius 4 and 2.
        let pln = plane_face(&plane_z(0.0));
        let tor = torus_face(3.0, 1.0);
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(tor);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        let res = ff.result();
        assert_eq!(res.nb_curves(), 2, "perpendicular cut through the center gives two circles");
        let mut radii: Vec<f64> = res
            .curves()
            .iter()
            .map(|c| c.curve.d0(0.0).distance(&GpPnt::zero()))
            .collect();
        radii.sort_by(f64::total_cmp);
        assert!((radii[0] - 2.0).abs() < 1e-6, "minor circle radius {}", radii[0]);
        assert!((radii[1] - 4.0).abs() < 1e-6, "major circle radius {}", radii[1]);
        let s_tor: Arc<dyn Surface> = Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap()));
        let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(plane_z(0.0)));
        for c in res.curves() {
            assert_points_on_both(c, s_pln.as_ref(), s_tor.as_ref(), 1e-4);
        }
    }

    #[test]
    fn general_bspline_sphere_uses_fallback() {
        // A genuinely curved B-spline patch (classifies Other, not Plane) × a
        // sphere goes through the sampling tracer. The traced polyline points
        // must lie near both surfaces.
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                row.push(GpPnt::new(u, v, u * u + v * v));
            }
            grid.push(row);
        }
        let patch_surf = intpatch::make_bspline_surface_from_grid(&grid, 1, 1).expect("fit");
        let patch_surf_check = patch_surf.clone();
        let patch_face = TopoBuilder::new().make_face(patch_surf, &[]);
        let sph = sphere_face(GpPnt::new(0.5, 0.5, 0.5), 1.0);

        let mut ff = FaceFace::new();
        ff.set_face1(patch_face);
        ff.set_face2(sph);
        ff.set_tolerance(0.05);
        ff.perform().expect("perform");
        let res = ff.result();
        assert!(res.nb_curves() >= 1, "tracer found an intersection");
        let s_sph: Arc<dyn Surface> = Arc::new(GeomSphere::new(
            GpSphere::new(
                GpAx3::new(
                    GpPnt::new(0.5, 0.5, 0.5),
                    GpDir::new(0.0, 0.0, 1.0).unwrap(),
                    &GpDir::new(1.0, 0.0, 0.0).unwrap(),
                )
                .unwrap(),
                1.0,
            )
            .unwrap(),
        ));
        for c in res.curves() {
            let (a, b) = (c.range.first, c.range.last);
            for i in 0..=12 {
                let t = a + (b - a) * i as f64 / 12.0;
                let p = c.curve.d0(t);
                assert!(intpatch::distance_to_surface(&p, patch_surf_check.as_ref()) < 0.2, "near patch: {p:?}");
                assert!(intpatch::distance_to_surface(&p, s_sph.as_ref()) < 0.2, "near sphere: {p:?}");
            }
        }
    }

    #[test]
    fn disjoint_sphere_plane_yields_no_curves() {
        let pln = plane_face(&plane_z(2.0)); // above the unit sphere
        let sph = sphere_face(GpPnt::zero(), 1.0);
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(sph);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        assert!(ff.result().is_empty(), "no intersection when the plane misses the sphere");
    }

    #[test]
    fn coplanar_planes_are_tangent() {
        let p1 = plane_face(&plane_z(0.0));
        let p2 = plane_face(&plane_z(0.0));
        let mut ff = FaceFace::new();
        ff.set_face1(p1);
        ff.set_face2(p2);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        assert!(ff.tangent_faces(), "coincident planes are tangent faces");
        assert!(ff.result().is_empty());
    }

    #[test]
    fn parallel_planes_do_not_intersect() {
        let p1 = plane_face(&plane_z(0.0));
        let p2 = plane_face(&plane_z(1.0));
        let mut ff = FaceFace::new();
        ff.set_face1(p1);
        ff.set_face2(p2);
        ff.set_tolerance(TOL);
        ff.perform().expect("perform");
        assert!(!ff.tangent_faces());
        assert!(ff.result().is_empty());
    }

    #[test]
    fn result_is_sorted_and_deduplicated() {
        let pln = plane_face(&plane_z(0.5));
        let sph = sphere_face(GpPnt::zero(), 1.0);
        let mut ff = FaceFace::new();
        ff.set_face1(pln);
        ff.set_face2(sph);
        ff.perform().expect("perform");
        let res = ff.result();
        let ranges: Vec<f64> = res.curves().iter().map(|c| c.range.first).collect();
        let mut sorted = ranges.clone();
        sorted.sort_by(f64::total_cmp);
        assert_eq!(ranges, sorted, "curves sorted by first parameter");
    }

    #[test]
    fn context_is_optional_noop() {
        let mut ff = FaceFace::new();
        ff.set_context(FaceFaceContext);
        assert!(ff.tangent_faces() == false);
        assert!(!ff.is_done());
    }
}
