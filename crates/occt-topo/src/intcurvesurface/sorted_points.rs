//! `IntCurveSurface_InterUtils`: start-point collection, sorting, exact
//! refinement and the parameter/transition helpers.
//!
//! Source: `IntCurveSurface_InterUtils.pxx:1346-1519` (`CollectInterferencePoints`,
//! `SortStartPoints`, `ProcessSortedPoints`, `SortedStartPoints`),
//! `:1116-1176` (`ComputeAppendPoint`), `:855-894` (`ComputeTransitions`),
//! `:1522-1609` (`UVBounds`, `DecomposeSurfaceIntervals`) and `:1611-1633`
//! (`ClampUVParameters`).

use occt_core::intf::{IntfInterference, IntfTangentZone};
use occt_core::precision;

use occt_math::MathFunctionSetRoot;

use super::exact_inter::TheExactHInter;
use super::interference::TheInterferenceOfHInter;
use super::polygon::ThePolygonOfHInter;
use super::polyhedron::ThePolyhedronOfHInter;
use super::prelude::*;
use super::section_point_params::section_point_to_parameters;
use super::types::{classify_curve, CurveKind, IntersectionPoint, State};

/// `THE_TOLERANCE_ANGULAIRE` (`IntCurveSurface_InterUtils.pxx:18`).
pub const THE_TOLERANCE_ANGULAIRE: f64 = occt_core::precision::ANGULAR;

/// `IntCurveSurface_InterUtils::SortedStartPoints`
/// (`IntCurveSurface_InterUtils.pxx:19-108`). 1-based arrays so the OCCT index
/// arithmetic of `SortStartPoints` is preserved literally.
#[derive(Debug, Default, Clone)]
pub struct SortedStartPoints {
    pub tab_u: Vec<f64>,
    pub tab_v: Vec<f64>,
    pub tab_w: Vec<f64>,
}

impl SortedStartPoints {
    /// `Clear` / default constructor.
    pub fn new() -> Self {
        Self {
            // Index 0 is unused (OCCT arrays are 1-based).
            tab_u: vec![0.0, 0.0],
            tab_v: vec![0.0, 0.0],
            tab_w: vec![0.0, 0.0],
        }
    }

    /// `Clear` (`...pxx:28-34`).
    pub fn clear(&mut self) {
        self.tab_u.clear();
        self.tab_v.clear();
        self.tab_w.clear();
    }

    /// `Append` (`...pxx:36-50`).
    pub fn append(&mut self, u: f64, v: f64, w: f64) {
        if self.tab_u.is_empty() {
            self.tab_u.push(0.0);
            self.tab_v.push(0.0);
            self.tab_w.push(0.0);
        }
        self.tab_u.push(u);
        self.tab_v.push(v);
        self.tab_w.push(w);
    }

    /// `Size` (`...pxx:52-60`).
    pub fn size(&self) -> usize {
        self.tab_u.len().saturating_sub(1)
    }

    /// `TabU(i)` (1-based).
    pub fn tab_u(&self, i: usize) -> f64 {
        self.tab_u[i]
    }
    /// `TabV(i)` (1-based).
    pub fn tab_v(&self, i: usize) -> f64 {
        self.tab_v[i]
    }
    /// `TabW(i)` (1-based).
    pub fn tab_w(&self, i: usize) -> f64 {
        self.tab_w[i]
    }
}

/// `CollectInterferencePoints` (`IntCurveSurface_InterUtils.pxx:1346-1375`).
pub fn collect_interference_points(
    interference: &TheInterferenceOfHInter,
    polyhedron: &ThePolyhedronOfHInter,
    polygon: &ThePolygonOfHInter,
    points: &mut SortedStartPoints,
) {
    points.clear();

    let base: &IntfInterference = interference.base();
    let nb_section_points = base.nb_section_points();
    let nb_tangent_zones = base.nb_tangent_zones();

    // `Intf_Interference::PntValue` / `ZoneValue` in this port store their
    // vectors 0-based (see `interference.rs:52,72`), while
    // `Intf_TangentZone::GetPoint` keeps OCCT's 1-based convention.
    for i in 0..nb_section_points {
        let sp = base.pnt_value(i);
        let (u, v, w) = section_point_to_parameters(sp, polyhedron, polygon);
        points.append(u, v, w);
    }

    for i in 0..nb_tangent_zones {
        let tz: &IntfTangentZone = base.zone_value(i);
        let nbpnts = tz.number_of_points();
        for j in 1..=nbpnts {
            let sp = tz.get_point(j);
            let (u, v, w) = section_point_to_parameters(sp, polyhedron, polygon);
            points.append(u, v, w);
        }
    }
}

/// `SortStartPoints` (`IntCurveSurface_InterUtils.pxx:1380-1446`).
pub fn sort_start_points(points: &mut SortedStartPoints) {
    let nb_start_points = points.size();
    if nb_start_points == 0 {
        return;
    }

    let ptol = 10.0 * precision::PCONFUSION;

    // Sort by W.
    let mut triok;
    loop {
        triok = true;
        for i in 2..=nb_start_points {
            let im1 = i - 1;
            if points.tab_w(i) < points.tab_w(im1) {
                points.tab_w.swap(i, im1);
                points.tab_u.swap(i, im1);
                points.tab_v.swap(i, im1);
                triok = false;
            }
        }
        if triok {
            break;
        }
    }

    // Sort by U for same W.
    loop {
        triok = true;
        for i in 2..=nb_start_points {
            let im1 = i - 1;
            if (points.tab_w(i) - points.tab_w(im1)) < ptol {
                points.tab_w[i] = points.tab_w(im1);
                if points.tab_u(i) < points.tab_u(im1) {
                    points.tab_u.swap(i, im1);
                    points.tab_v.swap(i, im1);
                    triok = false;
                }
            }
        }
        if triok {
            break;
        }
    }

    // Sort by V for same W and U.
    loop {
        triok = true;
        for i in 2..=nb_start_points {
            let im1 = i - 1;
            if ((points.tab_w(i) - points.tab_w(im1)) < ptol)
                && ((points.tab_u(i) - points.tab_u(im1)) < ptol)
            {
                points.tab_u[i] = points.tab_u(im1);
                if points.tab_v(i) < points.tab_v(im1) {
                    points.tab_v.swap(i, im1);
                    triok = false;
                }
            }
        }
        if triok {
            break;
        }
    }
}

/// `ProcessSortedPoints` (`IntCurveSurface_InterUtils.pxx:1452-1519`).
#[allow(clippy::too_many_arguments)]
pub fn process_sorted_points<'a>(
    exact_inter: &mut TheExactHInter<'a>,
    rsnld: &mut MathFunctionSetRoot,
    points: &SortedStartPoints,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    winf: f64,
    wsup: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
    result: &mut Vec<IntersectionPoint>,
) {
    result.clear();

    let nb_start_points = points.size();
    if nb_start_points == 0 {
        return;
    }

    let ptol = 10.0 * precision::PCONFUSION;
    let mut su = 0.0;
    let mut sv = 0.0;
    let mut sw = 0.0;

    for i in 1..=nb_start_points {
        let mut u = points.tab_u(i);
        let mut v = points.tab_v(i);
        let mut w = points.tab_w(i);

        if i == 1 {
            su = u - 1.0;
        }

        if (u - su).abs() > ptol || (v - sv).abs() > ptol || (w - sw).abs() > ptol {
            exact_inter.perform(u, v, w, rsnld, u0, u1, v0, v1, winf, wsup);
            if exact_inter.is_done() && !exact_inter.is_empty() {
                w = exact_inter.parameter_on_curve();
                let (uu, vv) = exact_inter.parameter_on_surface();
                u = uu;
                v = vv;

                if let Some(pt) = compute_append_point(curve, w, surface, u, v) {
                    result.push(pt);
                }
            }
        }
        su = points.tab_u(i);
        sv = points.tab_v(i);
        sw = points.tab_w(i);
    }
}

/// `ComputeAppendPoint` (`IntCurveSurface_InterUtils.pxx:1117-1176`). Returns
/// `None` where OCCT returns `false` (out of the curve / surface ranges).
pub fn compute_append_point(
    curve: &dyn Curve,
    lw: f64,
    surface: &dyn Surface,
    su: f64,
    sv: f64,
) -> Option<IntersectionPoint> {
    let w0 = curve.first_parameter();
    let w1 = curve.last_parameter();
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();

    let mut w = lw;
    let mut u = su;
    let mut v = sv;

    let a_ctype = classify_curve(curve);
    if curve.is_periodic() || a_ctype == CurveKind::Circle || a_ctype == CurveKind::Ellipse {
        w = in_period(w, w0, w0 + curve.period());
    }

    if (w0 - w) >= super::exact_inter::THE_TOLTANGENCY
        || (w - w1) >= super::exact_inter::THE_TOLTANGENCY
    {
        return None;
    }

    let a_kind = crate::brep_surface::classify_surface(surface);
    if surface.is_u_periodic()
        || a_kind == crate::brep_surface::SurfaceKind::Cylinder
        || a_kind == crate::brep_surface::SurfaceKind::Cone
        || a_kind == crate::brep_surface::SurfaceKind::Sphere
    {
        u = in_period(u, u0, u0 + surface.u_period());
    }

    if surface.is_v_periodic() {
        v = in_period(v, v0, v0 + surface.v_period());
    }

    if (u0 - u) >= super::exact_inter::THE_TOLTANGENCY
        || (u - u1) >= super::exact_inter::THE_TOLTANGENCY
    {
        return None;
    }
    if (v0 - v) >= super::exact_inter::THE_TOLTANGENCY
        || (v - v1) >= super::exact_inter::THE_TOLTANGENCY
    {
        return None;
    }

    let trans_on_curve = compute_transitions(curve, w, surface, u, v);
    let p = curve.d0(w);
    Some(IntersectionPoint::new(w, u, v, p, trans_on_curve))
}

/// `ComputeTransitions` (`IntCurveSurface_InterUtils.pxx:856-894`).
pub fn compute_transitions(
    curve: &dyn Curve,
    w: f64,
    surface: &dyn Surface,
    u: f64,
    v: f64,
) -> State {
    let (_, d1u, d1v) = surface.d1(u, v);
    let n_surf = d1u.crossed(&d1v);
    let (_, _d1w) = curve.d1(w);

    let norm = n_surf.magnitude();
    if norm > THE_TOLERANCE_ANGULAIRE && d1u.square_magnitude() > THE_TOLERANCE_ANGULAIRE {
        // `D1U.Normalize(); CosDir = NSurf.Dot(D1U); CosDir /= Norm;`
        let cos_dir = n_surf.normalized().dot(&d1u.normalized());
        if -cos_dir > THE_TOLERANCE_ANGULAIRE {
            //  --Curve-->    <----Surface----
            State::In
        } else if cos_dir > THE_TOLERANCE_ANGULAIRE {
            //  --Curve-->  ----Surface-->
            State::Out
        } else {
            State::On
        }
    } else {
        State::On
    }
}

/// `UVBounds` (`IntCurveSurface_InterUtils.pxx:1522-1544`).
#[derive(Debug, Clone, Copy)]
pub struct UVBounds {
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

impl UVBounds {
    pub fn new(u0: f64, u1: f64, v0: f64, v1: f64) -> Self {
        Self { u0, u1, v0, v1 }
    }
}

/// `DecomposeSurfaceIntervals` (`IntCurveSurface_InterUtils.pxx:1549-1609`).
pub fn decompose_surface_intervals(surface: &dyn Surface, intervals: &mut Vec<UVBounds>) {
    intervals.clear();

    let nb_u_on_s = surface.nb_u_intervals(2);
    let nb_v_on_s = surface.nb_v_intervals(2);

    if nb_u_on_s > 1 {
        let tab_u = surface.u_intervals(2);
        for iu in 0..nb_u_on_s as usize {
            let u0 = tab_u[iu];
            let u1 = tab_u[iu + 1];
            if nb_v_on_s > 1 {
                let tab_v = surface.v_intervals(2);
                for iv in 0..nb_v_on_s as usize {
                    intervals.push(UVBounds::new(u0, u1, tab_v[iv], tab_v[iv + 1]));
                }
            } else {
                let (v0, v1) = surface.v_range();
                intervals.push(UVBounds::new(u0, u1, v0, v1));
            }
        }
    } else if nb_v_on_s > 1 {
        let (u0, u1) = surface.u_range();
        let tab_v = surface.v_intervals(2);
        for iv in 0..nb_v_on_s as usize {
            intervals.push(UVBounds::new(u0, u1, tab_v[iv], tab_v[iv + 1]));
        }
    } else {
        let (u0, u1) = surface.u_range();
        let (v0, v1) = surface.v_range();
        intervals.push(UVBounds::new(u0, u1, v0, v1));
    }
}

/// `ClampUVParameters` (`IntCurveSurface_InterUtils.pxx:1613-1633`).
/// Protection from double type overflow in square magnitude computation.
pub fn clamp_uv_parameters(u1: &mut f64, u2: &mut f64, v1: &mut f64, v2: &mut f64) {
    const THE_PARAM_LIMIT: f64 = 1.0e50;
    if *u1 < -THE_PARAM_LIMIT {
        *u1 = -THE_PARAM_LIMIT;
    }
    if *u2 > THE_PARAM_LIMIT {
        *u2 = THE_PARAM_LIMIT;
    }
    if *v1 < -THE_PARAM_LIMIT {
        *v1 = -THE_PARAM_LIMIT;
    }
    if *v2 > THE_PARAM_LIMIT {
        *v2 = THE_PARAM_LIMIT;
    }
}
