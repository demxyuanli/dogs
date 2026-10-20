//! `ShapeAnalysis_Curve::Project` (`ShapeAnalysis_Curve.cxx`) port.

use super::super::prelude::*;
use super::*;

/// `ShapeAnalysis_Curve::ProjectOnSegments` (`ShapeAnalysis_Curve.cxx:81-122`).
fn project_on_segments(
    curve: &dyn Curve,
    point: &GpPnt,
    segment_count: i32,
    start: &mut f64,
    end: &mut f64,
    proj_dist: &mut f64,
    proj_param: &mut f64,
) {
    if segment_count <= 0 {
        return;
    }
    let step = (*end - *start) / segment_count as f64;
    let mut min_sq = *proj_dist * *proj_dist;
    let mut changed = false;
    for i in 0..=segment_count {
        let u = *start + step * i as f64;
        let q = curve.d0(u);
        let sq = point.square_distance(&q);
        if sq < min_sq {
            min_sq = sq;
            *proj_param = u;
            changed = true;
        }
    }
    if changed {
        *proj_dist = min_sq.sqrt();
    }
    *end = (*end).min(*proj_param + step);
    *start = (*start).max(*proj_param - step);
}

/// `ElCLib::LineParameter` (`ElCLib.cxx:1192-1195`).
fn elclib_line_parameter(curve: &dyn Curve, p: &GpPnt) -> f64 {
    let origin = curve.d0(0.0);
    let dir = curve.d1(0.0).1;
    p.coord
        .subtracted(&origin.coord)
        .dot(dir.xyz())
}

/// `ShapeAnalysis_Curve::ProjectAct` (`cxx:265-497`).
///
/// OCCT's `ProjectAct` returns a parameter on every path (`ShapeAnalysis_Curve.cxx`
/// has no failure channel), so this port has no `Option` either.
fn shape_analysis_project_act(
    curve: &dyn Curve,
    point: &GpPnt,
    preci: f64,
    u_min: f64,
    u_max: f64,
) -> (f64, f64) {
    let mut ok = false;
    let mut proj_param = 0.0;
    let mut computed_param = 0.0;
    let mut computed_dist = f64::INFINITY;
    let mut have_old = false;
    let mut old_param = 0.0;
    let mut old_dist = f64::INFINITY;
    let mut mod_min = f64::INFINITY;
    if let Some((t, d)) = crate::int_tools_vertex_line::extrema_project_in_range(
        curve, point, u_min, u_max,
    ) {
        have_old = true;
        old_param = t;
        old_dist = d;
        computed_param = t;
        computed_dist = d;
        proj_param = t;
        mod_min = d;
        ok = d <= preci;
    }

    let mut closed = false;
    let mut period = 0.0;
    if ok {
        // `ShapeAnalysis_Curve.cxx:340-344`: `theCurve.IsClosed()` only.
        // Do not treat `is_periodic` alone as closed (invent vs cxx) — that
        // forced AdjustByPeriod on Extrema hits and could shrink / wrap
        // edge ranges after always-Project (`TranslateEdge.cxx:442-444`).
        let lo = curve.d0(u_min);
        let hi = curve.d0(u_max);
        if lo.distance(&hi) <= occt_core::precision::CONFUSION {
            closed = true;
            period = u_max - u_min;
        }
    }

    if !ok {
        // `ShapeAnalysis_Curve.cxx:355-477`: `switch (theCurve.GetType())` with a
        // precise `ElCLib::Parameter` arm for every analytic type, and the
        // segmented `ProjectOnSegments` search only for the `default:` case
        // (B-spline / Bezier / other).
        if let Some(c) = curve.gp_circ() {
            let loc = c.position().location();
            if c.radius() <= occt_core::precision::RESOLUTION
                || point.square_distance(&loc) <= occt_core::precision::RESOLUTION
            {
                proj_param = u_min;
            } else {
                proj_param = occt_core::elib::clib::circle_parameter(&c.position(), point);
            }
            closed = true;
            period = 2.0 * std::f64::consts::PI;
        } else if let Some(h) = curve.gp_hyperbola() {
            proj_param = occt_core::elib::clib::parameter_hypr(&h, point);
        } else if let Some(prb) = curve.gp_parabola() {
            proj_param = occt_core::elib::clib::parameter_parab(&prb, point);
        } else if curve.is_line() {
            proj_param = elclib_line_parameter(curve, point);
        } else if let Some(el) = curve.gp_ellipse() {
            proj_param = occt_core::elib::clib::parameter_elips(&el, point);
            closed = true;
            period = 2.0 * std::f64::consts::PI;
        } else {
            let mut seg_lo = u_min;
            let mut seg_hi = u_max;
            let mut dist = f64::INFINITY;
            project_on_segments(
                curve,
                point,
                25,
                &mut seg_lo,
                &mut seg_hi,
                &mut dist,
                &mut proj_param,
            );
            if dist <= preci {
                return (proj_param, dist);
            }
            if let Some((t, q)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(
                curve, point, proj_param, u_min, u_max,
            ) {
                let d_newton = point.distance(&q);
                if d_newton < mod_min {
                    return (t, d_newton);
                }
            }
            for n in [40, 20, 25, 40] {
                project_on_segments(
                    curve,
                    point,
                    n,
                    &mut seg_lo,
                    &mut seg_hi,
                    &mut dist,
                    &mut proj_param,
                );
                if dist <= preci {
                    return (proj_param, dist);
                }
            }
            if dist > mod_min && have_old {
                return (computed_param, computed_dist);
            }
            return (proj_param, dist);
        }
    }

    if closed && (proj_param < u_min || proj_param > u_max) && period.abs() > 1e-100 {
        // `ShapeAnalysis::AdjustByPeriod` (`ShapeAnalysis.cxx:48-62`).
        let to_val = 0.5 * (u_min + u_max);
        let diff = proj_param - to_val;
        let d = diff.abs();
        let p = period.abs();
        if d > 0.5 * p {
            let shift = if diff > 0.0 { -p } else { p } * (d / p + 0.5).floor();
            proj_param += shift;
        }
    }

    let q = curve.d0(proj_param);
    let new_dist = point.distance(&q);
    if have_old && old_dist * old_dist < new_dist * new_dist {
        return (old_param, old_dist);
    }
    (proj_param, new_dist)
}

/// `ShapeAnalysis_Curve::Project` (`cxx:147-201`) with `AdjustToEnds=false`
/// as `MakeFromCurve3D` (`TranslateEdge.cxx:443-444`).
pub(in crate::step) fn shape_analysis_project(curve: &dyn Curve, p: &GpPnt, preci: f64) -> (f64, f64) {
    let (mut u_min, mut u_max) = (curve.first_parameter(), curve.last_parameter());
    if u_min > u_max {
        std::mem::swap(&mut u_min, &mut u_max);
    }
    if u_min.is_finite() && u_max.is_finite() {
        let low = curve.d0(u_min);
        let high = curve.d0(u_max);
        let dl = low.distance(p);
        if dl <= occt_core::precision::CONFUSION {
            return (u_min, dl);
        }
        let dh = high.distance(p);
        if dh <= occt_core::precision::CONFUSION {
            return (u_max, dh);
        }
        let closed = low.distance(&high) <= occt_core::precision::CONFUSION;
        if !closed {
            let delta = crate::int_tools_vertex_line::adaptor_resolution(curve, preci)
                .min((u_max - u_min) * 0.1);
            u_min -= delta;
            u_max += delta;
        }
    }
    shape_analysis_project_act(curve, p, preci, u_min, u_max)
}
