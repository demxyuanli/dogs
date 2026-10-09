//! Port of the `ShapeAnalysis_Curve` projection helpers
//! (`ShapeAnalysis_Curve.cxx`, TKShHealing).
//!
//! `ShapeAnalysis_TransferParametersProj` (`ShapeAnalysis_TransferParametersProj.cxx:200`,
//! `:209`, `:404`, `:405`, `:473`, `:474`) and the wire analysers ported next to
//! it (`ShapeAnalysis_Wire::CheckNotchedEdges` `cxx:1993`, `ProjectInside`
//! `cxx:1842`) project a 3D point onto a 3D curve or onto a curve-on-surface.
//! Only the overloads those call sites use are ported:
//!
//! * `Project(const Adaptor3d_Curve&, ...)` (`cxx:205-261`) - [`project_adaptor`]
//! * `Project(const Handle(Geom_Curve)&, ..., cf, cl, AdjustToEnds)`
//!   (`cxx:147-201`) - [`project_range`]
//! * `ProjectAct` (`cxx:265-497`) - [`project_act`]
//! * `ProjectOnSegments` (`cxx:81-122`) - [`project_on_segments`]
//! * `NextProject(paramPrev, const Adaptor3d_Curve&, ...)` (`cxx:561-579`) -
//!   [`next_project`]
//!
//! The `Geom_Curve` `NextProject` overload (`cxx:504-557`) and the `Geom_Curve`
//! `Project` overload without an explicit range (`cxx:126-142`, which only
//! forwards to `cxx:147` with the curve's own range) are already ported in
//! `pcurve_full/common.rs:1423` (`next_project_on_curve_range`), so they are not
//! duplicated here.
//!
//! The `!OK` switch of [`project_act`] is **fully ported** (batch 88, task
//! T-17): `GeomAbs_Circle` (`cxx:355-374`), `GeomAbs_Hyperbola` (`cxx:376-380`),
//! `GeomAbs_Parabola` (`cxx:382-386`), `GeomAbs_Line` (`cxx:388-392`) and
//! `GeomAbs_Ellipse` (`cxx:394-399`, including its `anIsClosedCurve` /
//! `aCurvePeriod = 2π` writes), each through the matching
//! `ElCLib::Parameter`/`Value` pair (`clib::{circle,hyperbola,parabola,ellipse}_parameter`
//! and `clib::{circle,hyperbola,parabola,ellipse,line}_value`;
//! `ElCLib::LineParameter` `ElCLib.cxx:1192-1195`). The three conic arms were
//! parked until the `Curve` trait gained `gp_hyperbola`/`gp_parabola`/`gp_ellipse`
//! (A20/A29) and `clib` gained the matching `*_parameter` helpers (T-56).
//!
//! `theCurve.IsClosed()` (`cxx:340`, `cxx:187`) has no `Curve` counterpart; a
//! curve whose two range ends coincide stands in for it. A trimmed arc of a
//! periodic basis must NOT count as closed: the `AdjustByPeriod` block
//! (`cxx:479-484`) would then wrap its range (`shhealing/wire_fix.rs`-style deviation,
//! see `step/read_geometry.rs:336-344`).

use occt_core::elib::clib;
use occt_core::gp::GpPnt;
use occt_core::precision::{Precision, COMPUTATIONAL, CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::extrema_pc::extrema_ext_pc_min_in_range;
use occt_geom::Curve;

use std::sync::Arc;

use super::adjust_by_period;

/// A projection result: OCCT's `dist`, `proj` and `param` out-parameters of
/// `ShapeAnalysis_Curve::Project` / `NextProject`.
#[derive(Clone, Copy, Debug)]
pub struct Projection {
    pub distance: f64,
    pub point: GpPnt,
    pub param: f64,
}

impl Projection {
    fn new(distance: f64, point: GpPnt, param: f64) -> Self {
        Self { distance, point, param }
    }
}

/// `ShapeAnalysis_Curve::ProjectOnSegments` (`ShapeAnalysis_Curve.cxx:81-122`).
///
/// `start` / `end` are narrowed around the best parameter found, so a repeated
/// call probes a smaller window (`cxx:118-121`). The projected point is not
/// carried: it is `curve.d0(proj_param)` (`cxx:107`).
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
    let step = (*end - *start) / f64::from(segment_count);
    let mut min_sq = *proj_dist * *proj_dist;
    let mut changed = false;
    for i in 0..=segment_count {
        let u = *start + step * f64::from(i);
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

/// `ElCLib::LineParameter` (`ElCLib.cxx:1192-1195`) without a `gp_Lin`
/// accessor: `(P - Origin) . Direction` on the curve's own parameterisation.
fn line_parameter(curve: &dyn Curve, p: &GpPnt) -> f64 {
    let origin = curve.d0(0.0);
    let dir = curve.d1(0.0).1;
    p.coord.subtracted(&origin.coord).dot(dir.xyz())
}

/// `C3D->IsKind(STANDARD_TYPE(Geom_BoundedCurve))` (`ShapeAnalysis_Curve.cxx:161`).
/// `Geom_BoundedCurve` is the parent of `Geom_BezierCurve`, `Geom_BSplineCurve`
/// and `Geom_TrimmedCurve` only - `Geom_Conic` (`Geom_Circle`, `Geom_Ellipse`,
/// `Geom_Hyperbola`, `Geom_Parabola`) is not one of them (`Geom_Curve.hxx:118`,
/// `Geom_BoundedCurve.hxx:24`, `Geom_BSplineCurve.hxx:42`,
/// `Geom_BezierCurve.hxx:41`, `Geom_TrimmedCurve.hxx:45`), which is exactly the
/// `is_geom_trimmed() || nurbs_degree().is_some()` pair here.
fn is_bounded_curve(c: &dyn Curve) -> bool {
    c.is_geom_trimmed() || c.nurbs_degree().is_some()
}

/// `ShapeAnalysis_Curve::ProjectAct` (`ShapeAnalysis_Curve.cxx:265-497`).
fn project_act(
    curve: &dyn Curve,
    point: &GpPnt,
    tolerance: f64,
    u_inf: f64,
    u_sup: f64,
) -> Projection {
    // `cxx:275-303`: `Extrema_ExtPC` on the range, closest `IsMin` solution.
    let mut proj_param = 0.0;
    let mut proj_point = GpPnt::zero();
    let mut proj_distance = f64::INFINITY;
    let mut mod_min = f64::INFINITY;
    let mut ok = false;
    if let Some((u, q, d)) = extrema_ext_pc_min_in_range(curve, point, u_inf, u_sup) {
        proj_param = u;
        proj_point = q;
        proj_distance = d;
        mod_min = d;
        // `cxx:332-338`: a solution farther than `theTolerance` is not accepted.
        ok = proj_distance <= tolerance;
    }

    // `cxx:322-328`: "remember the computed values ... used in case the
    // projection is not successful" (a default `gp_Pnt` when there is no
    // solution at all).
    let computed_param = proj_param;
    let computed_point = proj_point;
    let have_old = ok;
    let old_param = proj_param;
    let old_point = proj_point;

    let mut is_closed = false;
    let mut period = 0.0;
    if ok {
        // `cxx:340-344`: `theCurve.IsClosed()`.
        if curve.d0(u_inf).distance(&curve.d0(u_sup)) <= CONFUSION {
            is_closed = true;
            period = u_sup - u_inf; // `cxx:343`
        }
    }

    if !ok {
        // `cxx:351`: "Generally speaking, we try to ALWAYS return a result
        // that's NOT EVEN GOOD."
        proj_param = 0.0;
        if let Some(circ) = curve.gp_circ() {
            // `cxx:355-374`: `GeomAbs_Circle`.
            let location = circ.position().location();
            proj_point = location;
            if circ.radius() <= REAL_SMALL || point.square_distance(&location) <= REAL_SMALL {
                proj_param = u_inf;
                proj_point = GpPnt::from_xyz(
                    &location
                        .coord
                        .added(&circ.position().x_direction().xyz().multiplied(circ.radius())),
                );
            } else {
                proj_param = clib::circle_parameter(&circ.position(), point);
                proj_point = clib::circle_value(&circ, proj_param);
            }
            is_closed = true;
            period = 2.0 * std::f64::consts::PI;
        } else if let Some(hypr) = curve.gp_hyperbola() {
            // `cxx:376-380` (`GeomAbs_Hyperbola`).
            proj_param =
                clib::hyperbola_parameter(&hypr.pos, hypr.major_radius, hypr.minor_radius, point);
            proj_point = clib::hyperbola_value(&hypr, proj_param);
        } else if let Some(parab) = curve.gp_parabola() {
            // `cxx:382-386` (`GeomAbs_Parabola`).
            proj_param = clib::parabola_parameter(&parab.pos, point);
            proj_point = clib::parabola_value(&parab, proj_param);
        } else if let Some(lin) = curve.gp_line() {
            // `cxx:388-392` (`GeomAbs_Line`). The parameter is taken on the port
            // curve's own parameterisation (`ElCLib::Parameter`), the point from
            // the `gp_Lin` (`ElCLib::Value`).
            proj_param = line_parameter(curve, point);
            proj_point = clib::line_value(&lin, proj_param);
        } else if let Some(elips) = curve.gp_ellipse() {
            // `cxx:394-399` (`GeomAbs_Ellipse`).
            proj_param = clib::ellipse_parameter(
                &elips.pos,
                elips.major_radius,
                elips.minor_radius,
                point,
            );
            proj_point = clib::ellipse_value(&elips, proj_param);
            is_closed = true;
            period = 2.0 * std::f64::consts::PI;
        } else {
            // `cxx:401-477`: the `default` arm - `ProjectOnSegments(25)`, one
            // `Extrema_LocateExtPC`, then the 40/20/25/40 probe loop. Each arm
            // returns directly, skipping `cxx:479-496`.
            let mut seg_lo = u_inf;
            let mut seg_hi = u_sup;
            let mut seg_dist = f64::INFINITY; // `cxx:406`
            project_on_segments(
                curve,
                point,
                25,
                &mut seg_lo,
                &mut seg_hi,
                &mut seg_dist,
                &mut proj_param,
            );
            if seg_dist <= tolerance {
                // `cxx:417-420`
                return Projection::new(seg_dist, curve.d0(proj_param), proj_param);
            }
            // `cxx:422-437`: `Extrema_LocateExtPC(P, C, theProjParam, uMin, uMax, theTolerance)`.
            if let Some((t, q)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(
                curve, point, proj_param, u_inf, u_sup,
            ) {
                let newton_dist = point.distance(&q);
                if newton_dist < mod_min {
                    return Projection::new(newton_dist, q, t);
                }
            }
            for segment_count in [40, 20, 25, 40] {
                // `cxx:449-465`
                project_on_segments(
                    curve,
                    point,
                    segment_count,
                    &mut seg_lo,
                    &mut seg_hi,
                    &mut seg_dist,
                    &mut proj_param,
                );
                if seg_dist <= tolerance {
                    return Projection::new(seg_dist, curve.d0(proj_param), proj_param);
                }
            }
            // `cxx:469-477`: "we return the closest point found so far".
            if seg_dist > mod_min {
                return Projection::new(mod_min, computed_point, computed_param);
            }
            return Projection::new(seg_dist, curve.d0(proj_param), proj_param);
        }
    }

    // `cxx:479-484`: correcting on a periodic curve.
    if is_closed && (proj_param < u_inf || proj_param > u_sup) {
        proj_param += adjust_by_period(proj_param, 0.5 * (u_inf + u_sup), period);
    }

    // `cxx:486-495`: keep the better of the two solutions.
    if have_old {
        let old_sq = old_point.square_distance(point);
        let new_sq = proj_point.square_distance(point);
        if old_sq < new_sq {
            proj_point = old_point;
            proj_param = old_param;
        }
    }
    Projection::new(proj_point.distance(point), proj_point, proj_param)
}

/// `ShapeAnalysis_Curve::Project(const Adaptor3d_Curve& C3D, P3D, preci, proj,
/// param, AdjustToEnds)` (`ShapeAnalysis_Curve.cxx:205-261`).
pub fn project_adaptor(
    curve: &dyn Curve,
    point: &GpPnt,
    preci: f64,
    adjust_to_ends: bool,
) -> Projection {
    let u_min = curve.first_parameter();
    let u_max = curve.last_parameter();
    if Precision::is_infinite(u_min) && Precision::is_infinite(u_max) {
        return project_act(curve, point, preci, u_min, u_max); // `cxx:216-220`
    }
    let low = curve.d0(u_min); // `cxx:224-229`
    let high = curve.d0(u_max);
    let dist_low = low.distance(point);
    let dist_high = high.distance(point);
    // `cxx:221-222`: `prec = AdjustToEnds ? preci : Precision::Confusion()`.
    let prec = if adjust_to_ends { preci } else { CONFUSION };
    if dist_low <= prec {
        // `cxx:231-235`
        return Projection::new(dist_low, low, u_min);
    }
    if dist_high <= prec {
        // `cxx:237-242`
        return Projection::new(dist_high, high, u_max);
    }
    let proj = project_act(curve, point, preci, u_min, u_max); // `cxx:245`
    if proj.distance < dist_low + CONFUSION && proj.distance < dist_high + CONFUSION {
        return proj; // `cxx:246-250`
    }
    if dist_low < dist_high {
        Projection::new(dist_low, low, u_min)
    } else {
        Projection::new(dist_high, high, u_max)
    }
}

/// `ShapeAnalysis_Curve::Project(const Handle(Geom_Curve)& C3D, P3D, preci,
/// proj, param, cf, cl, AdjustToEnds)` (`ShapeAnalysis_Curve.cxx:147-201`).
pub fn project_range(
    curve: &dyn Curve,
    point: &GpPnt,
    preci: f64,
    cf: f64,
    cl: f64,
    adjust_to_ends: bool,
) -> Projection {
    let (u_min, u_max) = if cf < cl { (cf, cl) } else { (cl, cf) }; // `cxx:158-159`
    if is_bounded_curve(curve) {
        // `cxx:161-182`
        let prec = if adjust_to_ends { preci } else { CONFUSION };
        let low = curve.d0(u_min);
        let dist_low = low.distance(point);
        if dist_low <= prec {
            return Projection::new(dist_low, low, u_min);
        }
        let high = curve.d0(u_max);
        let dist_high = high.distance(point);
        if dist_high <= prec {
            return Projection::new(dist_high, high, u_max);
        }
    }
    let (mut lo, mut hi) = (u_min, u_max);
    if curve.d0(u_min).distance(&curve.d0(u_max)) > CONFUSION {
        // `cxx:184-198`: `!C3D->IsClosed()` widens the range before `ProjectAct`.
        let delta = curve.resolution(preci).min((u_max - u_min) * 0.1);
        lo -= delta;
        hi += delta;
    }
    project_act(curve, point, preci, lo, hi)
}

/// `ShapeAnalysis_Curve::NextProject(paramPrev, const Adaptor3d_Curve& C3D,
/// P3D, preci, proj, param)` (`ShapeAnalysis_Curve.cxx:561-579`).
pub fn next_project(curve: &dyn Curve, prev: f64, point: &GpPnt, preci: f64) -> Projection {
    let u_min = curve.first_parameter(); // `cxx:566-571`
    let u_max = curve.last_parameter();
    if let Some((t, q)) =
        crate::int_tools_vertex_line::extrema_locate_ext_pc(curve, point, prev, u_min, u_max)
    {
        return Projection::new(point.distance(&q), q, t);
    }
    // `cxx:578`: `Project(C3D, P3D, preci, proj, param, false)`.
    project_adaptor(curve, point, preci, false)
}

/// `Geom_Curve::IsClosed()` for the `Geom_BoundedCurve` members
/// `ShapeAnalysis_Curve::ValidateRange` branches on.
///
/// * `Geom_TrimmedCurve::IsClosed` (`Geom_TrimmedCurve.cxx:155-168`): a trim
///   covering an exact multiple of the basis period is closed, otherwise the
///   two range ends decide (`StartPoint`/`EndPoint` are the *basis* values at
///   `uTrim1`/`uTrim2`, `Geom_TrimmedCurve.cxx:271-287`).
/// * `Geom_BSplineCurve::IsClosed` (`Geom_BSplineCurve_1.cxx:146-149`) and
///   `Geom_BezierCurve::IsClosed` (`Geom_BezierCurve.cxx`): the two range ends
///   coincide within `Precision::Computational()`.
///
/// `Geom_Circle` / `Geom_Ellipse` report closed unconditionally, but
/// `ValidateRange`'s first branch only tests `IsClosed` inside
/// `IsKind(Geom_BoundedCurve)`, which excludes the conics.
fn geom_curve_is_closed(c: &dyn Curve) -> bool {
    if c.is_geom_trimmed() {
        if let Some((basis, _, _)) = c.untrimmed_basis() {
            let period = basis.period();
            let length = c.last_parameter() - c.first_parameter();
            if basis.is_periodic()
                && period > 0.0
                && length > PCONFUSION
                && (length - period * (length / period).round()).abs() <= PCONFUSION
            {
                return true;
            }
        }
    }
    let (f, l) = (c.first_parameter(), c.last_parameter());
    if !f.is_finite() || !l.is_finite() {
        return false;
    }
    c.d0(f).square_distance(&c.d0(l)) <= occt_core::precision::COMPUTATIONAL
}

/// `ShapeAnalysis_Curve::IsPeriodic(const Handle(Geom_Curve)&)`
/// (`ShapeAnalysis_Curve.cxx:1450-1469`): unwrap `Geom_OffsetCurve` and
/// `Geom_TrimmedCurve` down to the basis, then ask `IsPeriodic`.
fn is_periodic_curve(c: &Arc<dyn Curve>) -> bool {
    let mut cur = c.clone();
    // The `while` in cxx:1456-1466 unwraps one level per iteration; the bound
    // only guards against a self-referential offset/trim cycle.
    for _ in 0..32 {
        if let Some((basis, _)) = cur.offset_curve() {
            cur = basis;
            continue;
        }
        if cur.is_geom_trimmed() {
            if let Some((basis, _, _)) = cur.untrimmed_basis() {
                cur = basis;
                continue;
            }
        }
        break;
    }
    cur.is_periodic()
}

/// `ShapeAnalysis_Curve::ValidateRange(const Handle(Geom_Curve)&, First, Last,
/// preci)` (`ShapeAnalysis_Curve.cxx:586-732`).
///
/// `First` / `Last` are the `double&` out-parameters; returns the `bool` result.
/// The two `theCurve->Reverse()` calls (`cxx:704`, `cxx:721`) mutate the curve
/// object in place through the OCCT handle. The port's `Curve` is an immutable
/// `Arc`, and both arms require `First > Last`, which [`validate_range`]'s only
/// caller (`ShapeFix_SplitTool::CutEdge`, `ShapeFix_SplitTool.cxx:283`) cannot
/// produce: it passes `na = min(pend, cut)` / `nb = max(pend, cut)` after
/// rejecting `|cut - pend| < 10 * PConfusion`, so `na < nb` strictly. The
/// numeric adjustment of `First` / `Last` is reproduced; the in-place curve
/// reversal is marked UNPORTED for those two unreachable arms.
pub fn validate_range(curve: &dyn Curve, first: &mut f64, last: &mut f64, preci: f64) -> bool {
    let cf = curve.first_parameter(); // `cxx:593-594`
    let cl = curve.last_parameter();

    if is_bounded_curve(curve) && !geom_curve_is_closed(curve) {
        // `cxx:599-615`: clamp `First` / `Last` onto `[cf, cl]`.
        if *first < cf {
            *first = cf;
        } else if *first > cl {
            *first = cl;
        }
        if *last < cf {
            *last = cf;
        } else if *last > cl {
            *last = cl;
        }
    }

    // `cxx:618-621`: `IsPeriodic` on the unwrapped basis, `PConfusion` as the
    // precision (the `preci` parameter is deliberately not used here).
    if is_periodic_curve(&Arc::<dyn Curve>::from(curve.clone_dyn())) {
        // `ElCLib::AdjustPeriodic` (`ElCLib.cxx:115-148`); the port keeps this
        // shared `ElCLib` routine in `clib2d`.
        occt_core::elib::clib2d::adjust_periodic(cf, cl, PCONFUSION, first, last);
        return true;
    }

    if *first < *last {
        // `cxx:623-626`: nothing to fix.
        return true;
    }

    if geom_curve_is_closed(curve) {
        // `cxx:628-666`: one of the projected points sits on the 3D curve's
        // parameterisation origin, so the algorithm returned `cl +- preci`
        // instead of `cf` or vice versa.
        if (*last - cf).abs() < PCONFUSION {
            *last = cl;
        } else if (*first - cl).abs() < PCONFUSION {
            *first = cf;
        } else {
            // `cxx:650-659`: S4136 - check in 3D.
            if curve.d0(*first).distance(&curve.d0(cf)) < preci {
                *first = cf;
            }
            if curve.d0(*last).distance(&curve.d0(cl)) < preci {
                *last = cl;
            }
        }
        if *first > *last {
            std::mem::swap(first, last);
        }
        return true;
    }

    // `cxx:668-711`: the curve is closed within the 3D tolerance.
    if curve.nurbs_degree().is_some() {
        // `cxx:670-672`: `aBSpline->StartPoint().Distance(aBSpline->EndPoint()) <= preci`.
        let closed = curve
            .d0(cf)
            .distance(&curve.d0(cl))
            <= preci;
        if closed {
            let (f, l) = (*first, *last);
            if (*last - cf).abs() < PCONFUSION {
                *last = cl;
            } else if (*first - cl).abs() < PCONFUSION {
                *first = cf;
            } else {
                *first = l;
                *last = f;
            }
        } else {
            // `cxx:704-708`: `First = theCurve->ReversedParameter(First); ...;
            // theCurve->Reverse();`. UNPORTED - `Curve` has no
            // `ReversedParameter` and the in-place `Reverse()` cannot be
            // expressed on an `Arc`. Unreachable from `CutEdge` (see the doc
            // comment), so `First` / `Last` are left as they are.
        }
        // `cxx:709-713`: PRO7656 - a degenerate range is reset to the full range.
        if *first == *last {
            *first = cf;
            *last = cl;
            return false;
        }
        return true;
    }

    // `cxx:715-731`: `else` arm - not a B-spline and not closed in 3D.
    if *first > *last {
        // `cxx:719-723`: `First = theCurve->ReversedParameter(First); ...;
        // theCurve->Reverse();`. UNPORTED for the same reason as `cxx:704-708`
        // above; unreachable from `CutEdge`.
    }
    // `cxx:723-728`: PTV OCC966 - widen a degenerate range by `PConfusion`.
    if *first == *last {
        *first -= PCONFUSION;
        *last += PCONFUSION;
    }
    false
}
