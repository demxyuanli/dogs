//! Walking-line helpers. Source: `IntTools_WLineTool.cxx`.
//!
//! `NotUseSurfacesForApprox` reports whether a walking-line segment starts or
//! ends in a degenerated zone of either face, in which case surface normals
//! must not drive the 3D approximation.
//!
//! `DecompositionOfWLine` splits a walking line at periodic-boundary crossings
//! so each piece stays in one period of each surface.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec2d};
use occt_core::precision::{ANGULAR, PCONFUSION, RESOLUTION};
use occt_geom::Surface;

use crate::brep_tool::BRepTool;
use crate::intpatch::IntersectionCurve;
use crate::shape::Face;

#[path = "int_tools_wline_decomp.rs"]
mod decomp;
pub use decomp::decomposition_of_wline;

/// One point of a walking line (`IntSurf_PntOn2S`).
#[derive(Debug, Clone, Copy)]
pub struct PntOn2S {
    /// 3D intersection point.
    pub p: GpPnt,
    /// Parameters on surface 1.
    pub u1: f64,
    pub v1: f64,
    /// Parameters on surface 2.
    pub u2: f64,
    pub v2: f64,
}

impl PntOn2S {
    /// 3D value (`IntSurf_PntOn2S::Value`).
    pub fn value(&self) -> GpPnt {
        self.p
    }

    /// `(u1, v1, u2, v2)` (`Parameters`).
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.v1, self.u2, self.v2)
    }

    /// Parameters on surface 1 (`ParametersOnS1`).
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        (self.u1, self.v1)
    }

    /// Parameters on surface 2 (`ParametersOnS2`).
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        (self.u2, self.v2)
    }

    /// `SetValue(onFirst, u, v)` — replace one surface's UV.
    pub fn set_uv_on(&mut self, on_first: bool, u: f64, v: f64) {
        if on_first {
            self.u1 = u;
            self.v1 = v;
        } else {
            self.u2 = u;
            self.v2 = v;
        }
    }

    /// `SetValue(P, u1, v1, u2, v2)`.
    pub fn set_value(&mut self, p: GpPnt, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.p = p;
        self.u1 = u1;
        self.v1 = v1;
        self.u2 = u2;
        self.v2 = v2;
    }

    /// `IntSurf_PntOn2S::IsSame`. `tol2d < 0` compares 3D only.
    pub fn is_same(&self, other: &Self, tol3d: f64, tol2d: f64) -> bool {
        if self.p.square_distance(&other.p) > tol3d * tol3d {
            return false;
        }
        if tol2d < 0.0 {
            return true;
        }
        let d1 = (self.u1 - other.u1).hypot(self.v1 - other.v1);
        if d1 > tol2d {
            return false;
        }
        let d2 = (self.u2 - other.u2).hypot(self.v2 - other.v2);
        d2 <= tol2d
    }
}

/// `IntSurf_TypeTrans`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TransType {
    In,
    Out,
    Touch,
    #[default]
    Undecided,
}

/// How an [`WLine`] was created (`IntPatch_WLine::IntPatch_WL*`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WLineWay {
    ImpImp,
    ImpPrm,
    #[default]
    PrmPrm,
}

/// `IntPatch_Point` subset used by GeomInt LineConstructor / LineTool.
#[derive(Debug, Clone, Copy)]
pub struct PatchPoint {
    /// 3D point.
    pub p: GpPnt,
    /// Parameter on the intersection line (`ParameterOnLine`).
    pub param_on_line: f64,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    pub on_dom_s1: bool,
    pub on_dom_s2: bool,
    pub trans1: TransType,
    pub trans2: TransType,
    pub is_multiple: bool,
    pub tolerance: f64,
}

impl PatchPoint {
    /// Point with UV on both surfaces and no domain flags.
    pub fn new(p: GpPnt, param_on_line: f64, u1: f64, v1: f64, u2: f64, v2: f64) -> Self {
        Self {
            p,
            param_on_line,
            u1,
            v1,
            u2,
            v2,
            on_dom_s1: false,
            on_dom_s2: false,
            trans1: TransType::Undecided,
            trans2: TransType::Undecided,
            is_multiple: false,
            tolerance: 0.0,
        }
    }

    /// `ParameterOnLine`.
    pub fn parameter_on_line(&self) -> f64 {
        self.param_on_line
    }

    /// `ParametersOnS1`.
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        (self.u1, self.v1)
    }

    /// `ParametersOnS2`.
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        (self.u2, self.v2)
    }
}

/// Walking line (`IntPatch_WLine` / `IntSurf_LineOn2S`).
#[derive(Debug, Clone, Default)]
pub struct WLine {
    /// Sampled points, 0-based. OCCT `Point(i)` is 1-based.
    pub points: Vec<PntOn2S>,
    /// Vertices on the walking line (`NbVertex` / `Vertex`).
    pub vertices: Vec<PatchPoint>,
    /// `HasFirstPoint` / `HasLastPoint`.
    pub has_first_point: bool,
    pub has_last_point: bool,
    /// Creating intersector (`GetCreatingWay`).
    pub creating_way: WLineWay,
}

impl WLine {
    /// Empty walking line.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            creating_way: WLineWay::PrmPrm,
        }
    }

    /// `NbVertex`.
    pub fn nb_vertex(&self) -> i32 {
        self.vertices.len() as i32
    }

    /// 1-based `Vertex(i)`.
    pub fn vertex(&self, i1: i32) -> &PatchPoint {
        &self.vertices[(i1 as usize).saturating_sub(1)]
    }

    /// `GetCreatingWay`.
    pub fn creating_way(&self) -> WLineWay {
        self.creating_way
    }

    /// `SetCreatingWayInfo`.
    pub fn set_creating_way(&mut self, way: WLineWay) {
        self.creating_way = way;
    }

    /// Number of points (`NbPnts`).
    pub fn nb_pnts(&self) -> i32 {
        self.points.len() as i32
    }

    /// 1-based point access (`Point(i)`).
    pub fn point(&self, i1: i32) -> &PntOn2S {
        &self.points[(i1 as usize).saturating_sub(1)]
    }

    /// 1-based mutable point access.
    pub fn point_mut(&mut self, i1: i32) -> &mut PntOn2S {
        &mut self.points[(i1 as usize).saturating_sub(1)]
    }

    /// Append a point (`IntSurf_LineOn2S::Add`).
    pub fn add(&mut self, p: PntOn2S) {
        self.points.push(p);
    }

    /// 1-based `IntSurf_LineOn2S::InsertBefore`.
    pub fn insert_before(&mut self, i1: i32, p: PntOn2S) {
        let i = (i1 as usize).saturating_sub(1).min(self.points.len());
        self.points.insert(i, p);
    }

    /// 1-based `IntSurf_LineOn2S::RemovePoint`.
    pub fn remove_point(&mut self, i1: i32) {
        let i = (i1 as usize).saturating_sub(1);
        if i < self.points.len() {
            self.points.remove(i);
        }
    }

    /// First/last vertices when none were attached (`ComputeVertexParameters` empty case).
    pub fn ensure_end_vertices(&mut self) {
        if !self.vertices.is_empty() || self.points.is_empty() {
            return;
        }
        let a = self.points[0];
        self.vertices
            .push(PatchPoint::new(a.p, 1.0, a.u1, a.v1, a.u2, a.v2));
        self.has_first_point = true;
        if self.points.len() > 1 {
            let n = self.points.len();
            let b = self.points[n - 1];
            self.vertices
                .push(PatchPoint::new(b.p, n as f64, b.u1, b.v1, b.u2, b.v2));
            self.has_last_point = true;
        }
    }

    /// Build from an `intpatch` intersection curve (3D samples + UV on both faces).
    pub fn from_intersection_curve(ic: &IntersectionCurve) -> Self {
        Self::from_intersection_curve_way(ic, WLineWay::PrmPrm)
    }

    /// `from_intersection_curve` with `IntPatch_WLine::SetCreatingWay`.
    pub fn from_intersection_curve_way(ic: &IntersectionCurve, way: WLineWay) -> Self {
        let n = ic.points.len().min(ic.on_a.len()).min(ic.on_b.len());
        let mut points = Vec::with_capacity(n);
        for i in 0..n {
            points.push(PntOn2S {
                p: ic.points[i],
                u1: ic.on_a[i].0,
                v1: ic.on_a[i].1,
                u2: ic.on_b[i].0,
                v2: ic.on_b[i].1,
            });
        }
        let mut vertices = Vec::new();
        let has_ends = n >= 1;
        if has_ends {
            let a = points[0];
            let b = points[n - 1];
            vertices.push(PatchPoint::new(a.p, 1.0, a.u1, a.v1, a.u2, a.v2));
            if n > 1 {
                vertices.push(PatchPoint::new(
                    b.p,
                    n as f64,
                    b.u1,
                    b.v1,
                    b.u2,
                    b.v2,
                ));
            }
        }
        Self {
            points,
            vertices,
            has_first_point: has_ends,
            has_last_point: has_ends && n > 1,
            creating_way: way,
        }
    }
}

/// Whether the walking-line ends at `ifprm` / `ilprm` sit in a degenerated
/// zone of `f1` or `f2`. Port of `IntTools_WLineTool::NotUseSurfacesForApprox`.
///
/// `first` / `last` are the `(u1, v1, u2, v2)` parameters of the walking-line
/// points at those indices (`IntSurf_PntOn2S::Parameters`).
pub fn not_use_surfaces_for_approx(
    f1: &Face,
    f2: &Face,
    first: (f64, f64, f64, f64),
    last: (f64, f64, f64, f64),
) -> bool {
    if is_point_in_degenerated_zone(first, f1, f2) {
        return true;
    }
    is_point_in_degenerated_zone(last, f1, f2)
}

/// 1-based walking-line overload used by `MakeCurve`.
pub fn not_use_surfaces_for_approx_wline(
    f1: &Face,
    f2: &Face,
    wl: &WLine,
    ifprm: i32,
    ilprm: i32,
) -> bool {
    if wl.nb_pnts() < 1 || ifprm < 1 || ilprm < 1 {
        return false;
    }
    let n = wl.nb_pnts();
    let i0 = ifprm.min(n);
    let i1 = ilprm.min(n);
    not_use_surfaces_for_approx(f1, f2, wl.point(i0).parameters(), wl.point(i1).parameters())
}

fn is_degenerated_zone(p2d: &GpPnt2d, s: &dyn Surface, i_dir: i32) -> bool {
    let b_flag = true;
    let (us1, us2) = s.u_range();
    let (vs1, vs2) = s.v_range();
    let xm = p2d.x();
    let ym = p2d.y();
    let pm = s.d0(xm, ym);
    let d_x = 1.0e-5;
    let d_y = 1.0e-5;
    let d_d = 1.0e-12;
    if i_dir == 1 {
        let xb = xm;
        let xe = xm;
        let mut yb = ym - d_y;
        if yb < vs1 {
            yb = vs1;
        }
        let mut ye = ym + d_y;
        if ye > vs2 {
            ye = vs2;
        }
        let pb = s.d0(xb, yb);
        let pe = s.d0(xe, ye);
        let d1 = pm.distance(&pb);
        let d2 = pm.distance(&pe);
        if d1 < d_d && d2 < d_d {
            return b_flag;
        }
        return !b_flag;
    } else if i_dir == 2 {
        let yb = ym;
        let ye = ym;
        let mut xb = xm - d_x;
        if xb < us1 {
            xb = us1;
        }
        let mut xe = xm + d_x;
        if xe > us2 {
            xe = us2;
        }
        let pb = s.d0(xb, yb);
        let pe = s.d0(xe, ye);
        let d1 = pm.distance(&pb);
        let d2 = pm.distance(&pe);
        if d1 < d_d && d2 < d_d {
            return b_flag;
        }
        return !b_flag;
    }
    !b_flag
}

pub(crate) fn u_resolution(s: &dyn Surface, tol3d: f64) -> f64 {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let u = if u0.is_finite() && u1.is_finite() {
        0.5 * (u0 + u1)
    } else {
        0.0
    };
    let v = if v0.is_finite() && v1.is_finite() {
        0.5 * (v0 + v1)
    } else {
        0.0
    };
    let (_, du, _) = s.d1(u, v);
    let mag = du.magnitude();
    if mag > 1e-30 {
        (tol3d / mag).max(1e-12)
    } else {
        tol3d
    }
}

pub(crate) fn v_resolution(s: &dyn Surface, tol3d: f64) -> f64 {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let u = if u0.is_finite() && u1.is_finite() {
        0.5 * (u0 + u1)
    } else {
        0.0
    };
    let v = if v0.is_finite() && v1.is_finite() {
        0.5 * (v0 + v1)
    } else {
        0.0
    };
    let (_, _, dv) = s.d1(u, v);
    let mag = dv.magnitude();
    if mag > 1e-30 {
        (tol3d / mag).max(1e-12)
    } else {
        tol3d
    }
}

fn is_point_in_degenerated_zone(
    p2s: (f64, f64, f64, f64),
    f1: &Face,
    f2: &Face,
) -> bool {
    let mut b_flag = true;
    let Some(s1) = BRepTool::face_surface(f1) else {
        return false;
    };
    let Some(s2) = BRepTool::face_surface(f2) else {
        return false;
    };
    let (us11, us12) = s1.u_range();
    let (vs11, vs12) = s1.v_range();
    // OCCT: `aS1->Bounds(US21, ...)` after fetching S2 (copy as written).
    let (us21, us22) = s1.u_range();
    let (vs21, vs22) = s1.v_range();
    let (u1, v1, u2, v2) = p2s;
    let a_delta = 1.0e-7;
    let mut a_d = u_resolution(s1.as_ref(), a_delta);
    let mut a_p2d = GpPnt2d::new(u1, v1);
    if (u1 - us11).abs() < a_d {
        b_flag = is_degenerated_zone(&a_p2d, s1.as_ref(), 1);
        if b_flag {
            return b_flag;
        }
    }
    if (u1 - us12).abs() < a_d {
        b_flag = is_degenerated_zone(&a_p2d, s1.as_ref(), 1);
        if b_flag {
            return b_flag;
        }
    }
    a_d = v_resolution(s1.as_ref(), a_delta);
    if (v1 - vs11).abs() < a_delta {
        b_flag = is_degenerated_zone(&a_p2d, s1.as_ref(), 2);
        if b_flag {
            return b_flag;
        }
    }
    if (v1 - vs12).abs() < a_delta {
        b_flag = is_degenerated_zone(&a_p2d, s1.as_ref(), 2);
        if b_flag {
            return b_flag;
        }
    }
    a_d = u_resolution(s2.as_ref(), a_delta);
    a_p2d = GpPnt2d::new(u2, v2);
    if (u2 - us21).abs() < a_delta {
        b_flag = is_degenerated_zone(&a_p2d, s2.as_ref(), 1);
        if b_flag {
            return b_flag;
        }
    }
    if (u2 - us22).abs() < a_delta {
        b_flag = is_degenerated_zone(&a_p2d, s2.as_ref(), 1);
        if b_flag {
            return b_flag;
        }
    }
    a_d = v_resolution(s2.as_ref(), a_delta);
    if (v2 - vs21).abs() < a_delta {
        b_flag = is_degenerated_zone(&a_p2d, s2.as_ref(), 2);
        if b_flag {
            return b_flag;
        }
    }
    if (v2 - vs22).abs() < a_delta {
        b_flag = is_degenerated_zone(&a_p2d, s2.as_ref(), 2);
        if b_flag {
            return b_flag;
        }
    }
    let _ = a_d;
    !b_flag
}

/// `GeomInt::AdjustPeriodic` (`GeomInt.cxx:21`). Default `theEps` is 0.
pub fn adjust_periodic(
    the_par: f64,
    the_par_min: f64,
    the_par_max: f64,
    the_period: f64,
    the_eps: f64,
) -> (f64, f64) {
    let mut the_offset = 0.0;
    let mut the_new_par = the_par;
    let b_min = the_par_min - the_par > the_eps;
    let b_max = the_par - the_par_max > the_eps;
    if b_min || b_max {
        let dp = if b_min {
            the_par_max - the_par
        } else {
            the_par_min - the_par
        };
        let a_nb_per = (dp / the_period).trunc();
        the_offset = a_nb_per * the_period;
        the_new_par += the_offset;
    }
    (the_new_par, the_offset)
}

/// Geometric U period, or `None` when the surface is not U-periodic.
pub fn u_period(s: &dyn Surface) -> Option<f64> {
    if !s.is_u_periodic() {
        return None;
    }
    let (a, b) = s.u_range();
    if a.is_finite() && b.is_finite() && (b - a).abs() > 1e-30 {
        Some((b - a).abs())
    } else {
        Some(2.0 * std::f64::consts::PI)
    }
}

/// Geometric V period, or `None` when the surface is not V-periodic.
pub fn v_period(s: &dyn Surface) -> Option<f64> {
    if !s.is_v_periodic() {
        return None;
    }
    let (a, b) = s.v_range();
    if a.is_finite() && b.is_finite() && (b - a).abs() > 1e-30 {
        Some((b - a).abs())
    } else {
        Some(2.0 * std::f64::consts::PI)
    }
}

/// `IsPointOnBoundary` — static subfunction of `DecompositionOfWLine`.
pub fn is_point_on_boundary(
    the_parameter: f64,
    the_first_boundary: f64,
    the_second_boundary: f64,
    the_resolution: f64,
    is_on_first_boundary: &mut bool,
) -> bool {
    let b_ret = false;
    for i in 0..2 {
        *is_on_first_boundary = i == 0;
        let adist = if *is_on_first_boundary {
            (the_parameter - the_first_boundary).abs()
        } else {
            (the_parameter - the_second_boundary).abs()
        };
        if adist < the_resolution {
            return !b_ret;
        }
    }
    b_ret
}

/// `AdjustByNeighbour` — pick the periodic image of `original` nearest `neighbour`.
pub fn adjust_by_neighbour(
    the_neighbour: &GpPnt2d,
    the_original: &GpPnt2d,
    s: &dyn Surface,
) -> GpPnt2d {
    let ap1 = *the_neighbour;
    let mut ap2 = *the_original;
    if let Some(a_period) = u_period(s) {
        let mut a_sq_dist_min = 1.0e100;
        for p_it in -1..=1 {
            let a_p_test = GpPnt2d::new(the_original.x() + a_period * p_it as f64, the_original.y());
            let dd = ap1.square_distance(&a_p_test);
            if dd < a_sq_dist_min {
                ap2 = a_p_test;
                a_sq_dist_min = dd;
            }
        }
    }
    if let Some(a_period) = v_period(s) {
        let mut a_sq_dist_min = 1.0e100;
        for p_it in -1..=1 {
            let a_p_test = GpPnt2d::new(ap2.x(), the_original.y() + a_period * p_it as f64);
            let dd = ap1.square_distance(&a_p_test);
            if dd < a_sq_dist_min {
                ap2 = a_p_test;
                a_sq_dist_min = dd;
            }
        }
    }
    ap2
}

/// `RefineVector` — snap a 2D vector whose one component is nearly ±1.
pub fn refine_vector(a_v2d: &mut GpVec2d) {
    let a_eps = f64::EPSILON;
    let a_r1 = 1.0 - a_eps;
    let a_r2 = 1.0 + a_eps;
    let mut a_c = [a_v2d.x(), a_v2d.y()];
    for k in 0..2 {
        let m = (k + 1) % 2;
        let a_num = a_c[k].abs();
        if a_num > a_r1 && a_num < a_r2 {
            a_c[k] = if a_c[k] < 0.0 { -1.0 } else { 1.0 };
            a_c[m] = 0.0;
            break;
        }
    }
    a_v2d.set_coord(a_c[0], a_c[1]);
}

/// `FindPoint` — intersect the last→first walk with the UV rectangle.
pub fn find_point(
    the_first_point: &GpPnt2d,
    the_last_point: &GpPnt2d,
    the_umin: f64,
    the_umax: f64,
    the_vmin: f64,
    the_vmax: f64,
) -> Option<GpPnt2d> {
    let a_vec = GpVec2d::new(
        the_last_point.x() - the_first_point.x(),
        the_last_point.y() - the_first_point.y(),
    );
    for i in 0..4 {
        let (an_other_vec, an_other_vec_normal, mut aprojpoint) = if i % 2 == 0 {
            let mut p = *the_last_point;
            if i < 2 {
                p.set_x(the_umin);
            } else {
                p.set_x(the_umax);
            }
            (GpVec2d::new(0.0, 1.0), GpVec2d::new(1.0, 0.0), p)
        } else {
            let mut p = *the_last_point;
            if i < 2 {
                p.set_y(the_vmin);
            } else {
                p.set_y(the_vmax);
            }
            (GpVec2d::new(1.0, 0.0), GpVec2d::new(0.0, 1.0), p)
        };
        let Ok(mut anormvec) = a_vec.normalized() else {
            continue;
        };
        refine_vector(&mut anormvec);
        let adot1 = anormvec.dot(&an_other_vec_normal);
        if adot1.abs() < ANGULAR {
            continue;
        }
        let (adist, b_is_out) = if i % 2 == 0 {
            if i < 2 {
                (
                    (the_last_point.x() - the_umin).abs(),
                    the_last_point.x() < the_umin,
                )
            } else {
                (
                    (the_last_point.x() - the_umax).abs(),
                    the_last_point.x() > the_umax,
                )
            }
        } else if i < 2 {
            (
                (the_last_point.y() - the_vmin).abs(),
                the_last_point.y() < the_vmin,
            )
        } else {
            (
                (the_last_point.y() - the_vmax).abs(),
                the_last_point.y() > the_vmax,
            )
        };
        let mut anoffset = adist * an_other_vec.dot(&anormvec) / adot1;
        for j in 0..2 {
            anoffset = if j == 0 { anoffset } else { -anoffset };
            let shift = an_other_vec.multiplied_scalar(anoffset);
            let acurpoint = GpPnt2d::new(aprojpoint.x() + shift.x(), aprojpoint.y() + shift.y());
            let mut acurvec = GpVec2d::new(
                acurpoint.x() - the_last_point.x(),
                acurpoint.y() - the_last_point.y(),
            );
            if b_is_out {
                acurvec.reverse();
            }
            let a_dot_x = a_vec.dot(&acurvec);
            let an_angle_x = a_vec.angle(&acurvec);
            if a_dot_x > 0.0 && an_angle_x.abs() < PCONFUSION {
                if i % 2 == 0 {
                    if acurpoint.y() >= the_vmin && acurpoint.y() <= the_vmax {
                        return Some(acurpoint);
                    }
                } else if acurpoint.x() >= the_umin && acurpoint.x() <= the_umax {
                    return Some(acurpoint);
                }
            }
        }
        let _ = aprojpoint;
    }
    let _ = RESOLUTION;
    None
}
