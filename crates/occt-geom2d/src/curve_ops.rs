//! 2D curve operations: arc length, closest point, intersections.
//!
//! **Faithful**: [`curve2d_closest_point`] / [`curve2d_distance_to_point`] go
//! through `crate::extrema2d::point_curve_extrema2d` (`Extrema_ExtPC2d`, the
//! engine behind `Geom2dAPI_ProjectPointOnCurve`; task R2-18a, batch 85).
//! `lin2d_intersection` / `segment_intersection` / `project_point_on_segment`
//! are the analytic point/segment solves they claim to be.
//!
//! **UNPORTED (audit A16, task R2-18b)**: [`curve2d_intersections`] is still the
//! port's 256×256 sampler with alternating 1-D minimization. The faithful route
//! is `IntAna2d_AnaIntersection` for the analytic pairs (`perform_lin_lin`,
//! `perform_lin_circ`, `perform_circ_circ`, `perform_{lin,circ,elips,parab,hypr}_conic`
//! — ported in `occt_core::intana2d::ana_intersection`) and `Extrema_ExtCC2d`
//! (`crate::extrema2d::curve_curve_extrema2d_all`, whose general-curve seeding is
//! itself the A7-family substitute) for the rest. [`curve2d_length`] /
//! [`curve2d_length_tol`] are the faithful `GCPnts_AbscissaPoint::Length`
//! (`GCPnts_AbscissaPoint.cxx:305-423`) → `CPnts_AbscissaPoint::Length`
//! (`CPnts_AbscissaPoint.cxx:148-207`) → `math_GaussSingleIntegration` path.

use occt_core::gp::{GpLin2d, GpPnt2d, GpXY};
use crate::curve::Curve2d;

/// Arc length of a bounded 2D curve — the faithful `GCPnts_AbscissaPoint` path.
///
/// `GCPnts_AbscissaPoint::Length(theC)` (`GCPnts_AbscissaPoint.cxx:312-315`) is
/// `Length(theC, FirstParameter, LastParameter)` (`:342-347`), i.e. [`length_2d`]
/// with no tolerance. Returns `NaN` for unbounded curves (a port guard; OCCT
/// would integrate over an infinite range).
pub fn curve2d_length(c: &dyn Curve2d) -> f64 {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        return f64::NAN;
    }
    length_2d(c, a, b, None)
}

/// `GCPnts_AbscissaPoint::Length(theC, theU1, theU2, theTol)`
/// (`GCPnts_AbscissaPoint.cxx:361-367`, body `:371-423`): the tolerance overload
/// of [`curve2d_length`]. On the `Parametrized` arm of [`compute_type_2d`] it is
/// exactly `CPnts_AbscissaPoint::Length(C, U1, U2, Tol)`
/// (`CPnts_AbscissaPoint.cxx:191-207`); the length-parametrised arm ignores
/// `tol`, and the composite arm refines each `GeomAbs_CN` interval separately.
pub fn curve2d_length_tol(c: &dyn Curve2d, u1: f64, u2: f64, tol: f64) -> f64 {
    length_2d(c, u1, u2, Some(tol))
}

/// `GCPnts_AbscissaType` (`GCPnts_AbscissaType.hxx`) plus the length ratio that
/// `computeType` (`GCPnts_AbscissaPoint.cxx:26-65`) fills in for the
/// length-parametrised case.
enum Abs2d {
    LengthParametrized(f64),
    Parametrized,
    Composite,
}

/// `computeType` (`GCPnts_AbscissaPoint.cxx:26-65`): more than one
/// `GeomAbs_CN` interval (`GeomAbs_CN = 6`, `occt-core/src/kernel/geomabs.rs`)
/// ⇒ `GCPnts_AbsComposite`; otherwise a line, a circle, or a two-pole
/// **non-rational** Bezier / B-spline is length-parametrised
/// (`aBz->NbPoles() == 2 && !aBz->IsRational()`, `:45`;
/// `aBs->NbPoles() == 2 && !aBs->IsRational()`, `:54`); everything else is
/// integrated.
fn compute_type_2d(c: &dyn Curve2d) -> Abs2d {
    if nb_intervals_cn(c) > 1 {
        return Abs2d::Composite;
    }
    let base = adaptor2d(c);
    if base.is_line() {
        return Abs2d::LengthParametrized(1.0);
    }
    if let Some(circ) = base.gp_circ2d() {
        return Abs2d::LengthParametrized(circ.radius());
    }
    if let Some(n) = base.bezier_nb_poles() {
        if n == 2 && !base.is_rational() {
            // `aBz->DN(0, 1)`: the derivative of the Bezier's own `[0, 1]`
            // parametrisation at 0.
            return Abs2d::LengthParametrized(base.d1(0.0).1.magnitude());
        }
    }
    if let Some((xs, _ys)) = base.bspline_poles2d() {
        if xs.len() == 2 && !base.is_rational() {
            // `aBs->DN(aBs->FirstParameter(), 1)`.
            return Abs2d::LengthParametrized(base.d1(base.first_parameter()).1.magnitude());
        }
    }
    Abs2d::Parametrized
}

/// `GCPnts_AbscissaPoint::length` (`GCPnts_AbscissaPoint.cxx:371-423`). The
/// length-parametrised arm returns `|U2 − U1| × ratio` (`:381-383`); the
/// parametrized arm forwards to `CPnts_AbscissaPoint::Length` (`:384-387`); the
/// composite arm walks the `GeomAbs_CN` break points and adds the length of the
/// part of each interval that overlaps `[min(U1, U2), max(U1, U2)]` (`:388-420`).
fn length_2d(c: &dyn Curve2d, u1: f64, u2: f64, tol: Option<f64>) -> f64 {
    match compute_type_2d(c) {
        Abs2d::LengthParametrized(ratio) => (u2 - u1).abs() * ratio,
        Abs2d::Parametrized => match tol {
            Some(t) => gauss_length_tol_2d(c, u1, u2, t),
            None => gauss_length_2d(c, u1, u2),
        },
        Abs2d::Composite => {
            let ti = parameter_intervals_cn(c);
            let uu1 = u1.min(u2);
            let uu2 = u1.max(u2);
            let mut total = 0.0;
            for w in ti.windows(2) {
                if w[0] > uu2 {
                    break;
                }
                if w[1] < uu1 {
                    continue;
                }
                let lo = w[0].max(uu1);
                let hi = w[1].min(uu2);
                total += match tol {
                    Some(t) => gauss_length_tol_2d(c, lo, hi, t),
                    None => gauss_length_2d(c, lo, hi),
                };
            }
            total
        }
    }
}

/// `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:285-288`) stores the
/// **basis** of a `Geom2d_TrimmedCurve` (recursively, `:287`): every
/// `GetType()`-driven query (`Line()`, `Circle()`, `Bezier()`, `BSpline()`,
/// `IsRational()`, `Degree()`, `NbPoles()`) resolves on the basis, while
/// `FirstParameter`/`LastParameter` stay the trimmed range.
fn adaptor2d(c: &dyn Curve2d) -> &dyn Curve2d {
    let mut cur = c;
    while let Some(basis) = cur.trimmed_basis() {
        cur = basis;
    }
    cur
}

/// `Geom2dAdaptor_Curve::NbIntervals(GeomAbs_CN)`
/// (`Geom2dAdaptor_Curve.cxx:409-486`), resolved as the stored adaptor curve: a
/// trimmed curve unwraps to its basis (`:285-288`), and an offset curve asks a
/// basis adaptor with `default: BaseS = GeomAbs_CN` (`:472-474`).
fn nb_intervals_cn(c: &dyn Curve2d) -> i32 {
    parameter_intervals_cn(c).len().saturating_sub(1).max(1) as i32
}

/// `Geom2dAdaptor_Curve::Intervals(GeomAbs_CN)` (`:490-573`), resolved the same
/// way as [`nb_intervals_cn`]. OCCT builds an offset's `GeomAbs_CN` break points
/// from a basis adaptor over its own `[First, Last]` (`:559-563`) and then
/// forces the two end points to the offset's range (`:564-565`); the basis break
/// points are a superset of that, and [`length_2d`] clips them, so the sum is the
/// same.
fn parameter_intervals_cn(c: &dyn Curve2d) -> Vec<f64> {
    if let Some(basis) = c.offset_basis() {
        return parameter_intervals_cn(basis);
    }
    if let Some(basis) = c.trimmed_basis() {
        return parameter_intervals_cn(basis);
    }
    c.parameter_intervals(6)
}

/// `order(const Adaptor2d_Curve2d&)` (`CPnts_AbscissaPoint.cxx:79-100`):
/// `Line` 2, `Parabola` 5, `BezierCurve` `min(24, 2·Degree)`,
/// `BSplineCurve` `min(24, 2·NbPoles − 1)`, everything else 10. The type queries
/// resolve on the adaptor's stored curve ([`adaptor2d`]).
fn gauss_order_2d(c: &dyn Curve2d) -> usize {
    let base = adaptor2d(c);
    if base.is_line() {
        return 2;
    }
    if base.gp_parab2d().is_some() {
        return 5;
    }
    if let Some(n) = base.bezier_nb_poles() {
        return (2 * n.saturating_sub(1)).min(24).max(1);
    }
    if let Some((xs, _ys)) = base.bspline_poles2d() {
        return (2 * xs.len()).saturating_sub(1).min(24).max(1);
    }
    10
}

/// `math::GaussPointsMax()` (`math.cxx:25-28`, returns 61), the clamp
/// `math_GaussSingleIntegration` puts on the requested order
/// (`math_GaussSingleIntegration.cxx:60`, `:70`).
const GAUSS_POINTS_MAX: usize = 61;

/// `CPnts_AbscissaPoint::Length(C, U1, U2)` (`CPnts_AbscissaPoint.cxx:148-161`):
/// `|math_GaussSingleIntegration(|C'|, U1, U2, order(C))|`, the order clamped to
/// `min(math::GaussPointsMax(), order(C))`.
fn gauss_length_2d(c: &dyn Curve2d, u1: f64, u2: f64) -> f64 {
    let n = gauss_order_2d(c).min(GAUSS_POINTS_MAX);
    let f = |u: f64| c.d1(u).1.magnitude();
    occt_math::gauss::integrate(&f, u1, u2, n).abs()
}

/// `CPnts_AbscissaPoint::Length(C, U1, U2, Tol)`
/// (`CPnts_AbscissaPoint.cxx:191-207`) → `math_GaussSingleIntegration(FG, U1, U2,
/// order(C), Tol)` (`math_GaussSingleIntegration.cxx:64-98`): repeat the rule on
/// `2^k` equal sub-intervals (`IterMax = 13`) until two successive totals differ
/// by at most `tol`.
///
/// Measured: the no-tolerance overload is only ~1.3e-3 accurate on a full ellipse
/// (order 10 ⇒ 14.551758936716; the true perimeter is 14.532672330892 by a
/// Jacobi-series reference). This overload refines to 14.532672330821.
fn gauss_length_tol_2d(c: &dyn Curve2d, u1: f64, u2: f64, tol: f64) -> f64 {
    const ITER_MAX: usize = 13;
    let n = gauss_order_2d(c).min(GAUSS_POINTS_MAX);
    let f = |u: f64| c.d1(u).1.magnitude();
    let mut len = occt_math::gauss::integrate(&f, u1, u2, n);
    let mut nb_interval = 1usize;
    for _ in 0..ITER_MAX {
        let old_len = len;
        len = 0.0;
        nb_interval *= 2;
        let du = (u2 - u1) / nb_interval as f64;
        for i in 0..nb_interval {
            len += occt_math::gauss::integrate(
                &f,
                u1 + i as f64 * du,
                u1 + (i + 1) as f64 * du,
                n,
            );
        }
        if (old_len - len).abs() <= tol {
            break;
        }
    }
    len.abs()
}

/// Point on the curve at parameter `u`.
pub fn curve2d_point(c: &dyn Curve2d, u: f64) -> GpPnt2d {
    c.d0(u)
}

/// Closest point on the curve to `p`. Returns `(parameter, point)`.
///
/// Faithful route: `Geom2dAPI_ProjectPointOnCurve::Perform`
/// (`Geom2dAPI_ProjectPointOnCurve.cxx:78-92` in 8.0.0 is
/// `myExtPC.Perform(P)`, taking the smallest-distance solution) →
/// `Extrema_ExtPC2d`, which [`crate::extrema2d::point_curve_extrema2d`] ports:
/// lines/circles through `Extrema_ExtPElC2d`'s analytic arms, everything else
/// through the Extrema engine (whose general-curve seeding is still the
/// A7-family substitute, see `extrema2d::curve_curve`). The previous body was a
/// 64-sample scan plus golden section with an expanding-window hack for
/// unbounded curves (audit A16, task R2-18a).
pub fn curve2d_closest_point(c: &dyn Curve2d, p: &GpPnt2d, _tol: f64) -> Option<(f64, GpPnt2d)> {
    let e = crate::extrema2d::point_curve_extrema2d(c, p);
    if e.u1.is_finite() && e.p2.x().is_finite() && e.p2.y().is_finite() {
        Some((e.u1, e.p2))
    } else {
        None
    }
}

/// Distance from `p` to the nearest point on the curve.
pub fn curve2d_distance_to_point(c: &dyn Curve2d, p: &GpPnt2d, tol: f64) -> f64 {
    curve2d_closest_point(c, p, tol).map_or(f64::INFINITY, |(_, q)| q.distance(p))
}

/// Approximate intersection points of two curves.
/// Returns `(u, v, point)` where `u`/`v` are the parameters on `a`/`b`.
///
/// **Analytic pairs go through `IntAna2d_AnaIntersection`** (the faithful port in
/// [`occt_core::intana2d`]): line/line (`perform_lin_lin`), line/circle
/// (`perform_lin_circ`), circle/circle (`perform_circ_circ`), line/conic
/// (`perform_lin_conic`), circle/conic (`perform_circ_conic`) and
/// ellipse/parabola/hyperbola against a conic (`perform_{elips,parab,hypr}_conic`).
/// Pairs whose first curve is not the specialized operand are evaluated with the
/// arguments swapped and their parameters swapped back, so the result is always
/// expressed on `a`/`b` as given.
///
/// **UNPORTED (audit A16, task R2-18b)**: every other pair — B-spline, Bezier,
/// offset, trimmed-of-those — still runs the port's 256×256 sampler with
/// alternating 1-D minimization, which is not root-exact. The faithful route for
/// them is `Extrema_ExtCC2d`
/// ([`crate::extrema2d::curve_curve_extrema2d_all`], whose general-curve seeding
/// is itself the A7-family substitute) — see the module header.
pub fn curve2d_intersections(a: &dyn Curve2d, b: &dyn Curve2d, tol: f64) -> Vec<(f64, f64, GpPnt2d)> {
    if let Some(points) = analytic_intersections2d(a, b, tol) {
        return points;
    }
    let tol = tol.max(1e-12);
    let na = 256;
    let nb = 256;
    let sa = sample_curve(a, b, na);
    let sb = sample_curve(b, a, nb);
    if sa.is_empty() || sb.is_empty() {
        return Vec::new();
    }
    let ua0 = sa[0].0;
    let ua1 = sa[sa.len() - 1].0;
    let vb0 = sb[0].0;
    let vb1 = sb[sb.len() - 1].0;
    let du = ((ua1 - ua0) / na as f64).max(1e-12);
    let dv = ((vb1 - vb0) / nb as f64).max(1e-12);

    let coarse = tol * 100.0 + 1e-9;
    let mut found: Vec<(f64, f64)> = Vec::new();
    for (ui, ai) in &sa {
        for (vj, bj) in &sb {
            if ai.square_distance(bj) < coarse * coarse {
                // Refine the candidate pair by alternating 1-D minimization.
                let mut u = *ui;
                let mut v = *vj;
                for _ in 0..6 {
                    let (nu, _) = minimize_1d(
                        &|t| a.d0(t).square_distance(&b.d0(v)),
                        (u - du).max(ua0),
                        (u + du).min(ua1),
                    );
                    u = nu;
                    let (nv, _) = minimize_1d(
                        &|t| a.d0(u).square_distance(&b.d0(t)),
                        (v - dv).max(vb0),
                        (v + dv).min(vb1),
                    );
                    v = nv;
                }
                if a.d0(u).distance(&b.d0(v)) <= tol {
                    found.push((u, v));
                }
            }
        }
    }

    // Dedupe candidates that converge to the same geometric point.
    let mut out: Vec<(f64, f64, GpPnt2d)> = Vec::new();
    for (u, v) in found {
        let pa = a.d0(u);
        let pb = b.d0(v);
        let mid = GpPnt2d::new((pa.x() + pb.x()) * 0.5, (pa.y() + pb.y()) * 0.5);
        let dup = out.iter().any(|(_, _, p)| p.distance(&mid) <= tol.max(1e-6));
        if !dup {
            out.push((u, v, mid));
        }
    }
    out
}

/// The analytic family of a 2D curve, in `Geom2dAdaptor_Curve::GetType` terms
/// (`Geom2dAdaptor_Curve.cxx:96-120`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Analytic2d {
    Lin,
    Circ,
    Elips,
    Parab,
    Hypr,
}

fn analytic_kind2d(c: &dyn Curve2d) -> Option<Analytic2d> {
    if c.gp_lin2d().is_some() {
        Some(Analytic2d::Lin)
    } else if c.gp_circ2d().is_some() {
        Some(Analytic2d::Circ)
    } else if c.gp_elips2d().is_some() {
        Some(Analytic2d::Elips)
    } else if c.gp_parab2d().is_some() {
        Some(Analytic2d::Parab)
    } else if c.gp_hypr2d().is_some() {
        Some(Analytic2d::Hypr)
    } else {
        None
    }
}

/// The `IntAna2d_Conic` behind a 2D curve (`IntAna2d_Conic::SetLin2d` and
/// friends, `IntAna2d_Conic.cxx`).
fn conic_of2d(c: &dyn Curve2d) -> Option<occt_core::intana2d::IntAna2dConic> {
    use occt_core::intana2d::IntAna2dConic;
    if let Some(l) = c.gp_lin2d() {
        return Some(IntAna2dConic::from_lin2d(&l));
    }
    if let Some(ci) = c.gp_circ2d() {
        return Some(IntAna2dConic::from_circ2d(&ci));
    }
    if let Some(e) = c.gp_elips2d() {
        return Some(IntAna2dConic::from_elips2d(&e));
    }
    if let Some(p) = c.gp_parab2d() {
        return Some(IntAna2dConic::from_parab2d(&p));
    }
    if let Some(h) = c.gp_hypr2d() {
        return Some(IntAna2dConic::from_hypr2d(&h));
    }
    None
}

/// Is `u` a parameter of the bounded curve `c`?
///
/// `IntAna2d_AnaIntersection` works on the *unbounded* conics, while callers of
/// [`curve2d_intersections`] hand in bounded curves (and `Geom2dAPI_InterCurveCurve`
/// intersects the bounded adaptors), so a conic root outside a curve's own range
/// is not an intersection of the two curves. Periodic curves are tested through
/// one period starting at their first parameter.
fn param_on_curve2d(c: &dyn Curve2d, u: f64, tol: f64) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return true;
    }
    let u = if c.is_periodic() {
        let p = c.period();
        if p > 0.0 {
            a + (u - a).rem_euclid(p)
        } else {
            u
        }
    } else {
        u
    };
    u >= a - tol && u <= b + tol
}

/// `IntAna2d_AnaIntersection` for the analytic pairs; `None` when the pair is not
/// analytic (the caller falls back to the sampler).
///
/// The **first** curve supplies the specialized operand of the OCCT overload —
/// exactly the eight `Perform` methods the port exposes (`Lin,Lin`),
/// (`Lin,Circ`), (`Circ,Circ`), (`Lin,Conic`), (`Circ,Conic`), (`Elips,Conic`),
/// (`Parab,Conic`), (`Hypr,Conic`) — and the second is taken as an
/// `IntAna2d_Conic`. The one exception is `(Circ, Lin)`, which OCCT writes as
/// `Perform(Lin, Circ)`; the port calls it with the arguments swapped and swaps
/// the resulting parameters back, so the answer is always on `a`/`b` as given.
fn analytic_intersections2d(
    a: &dyn Curve2d,
    b: &dyn Curve2d,
    tol: f64,
) -> Option<Vec<(f64, f64, GpPnt2d)>> {
    use occt_core::intana2d::IntAna2dAnaIntersection;
    let (ka, kb) = (analytic_kind2d(a)?, analytic_kind2d(b)?);
    let cb = conic_of2d(b)?;

    let mut ana = IntAna2dAnaIntersection::new();
    let swap;
    match (ka, kb) {
        (Analytic2d::Lin, Analytic2d::Lin) => {
            ana.perform_lin_lin(&a.gp_lin2d()?, &b.gp_lin2d()?);
            swap = false;
        }
        (Analytic2d::Lin, Analytic2d::Circ) => {
            ana.perform_lin_circ(&a.gp_lin2d()?, &b.gp_circ2d()?);
            swap = false;
        }
        (Analytic2d::Circ, Analytic2d::Lin) => {
            ana.perform_lin_circ(&b.gp_lin2d()?, &a.gp_circ2d()?);
            swap = true;
        }
        (Analytic2d::Circ, Analytic2d::Circ) => {
            ana.perform_circ_circ(&a.gp_circ2d()?, &b.gp_circ2d()?);
            swap = false;
        }
        (Analytic2d::Lin, _) => {
            ana.perform_lin_conic(&a.gp_lin2d()?, &cb);
            swap = false;
        }
        (Analytic2d::Circ, _) => {
            ana.perform_circ_conic(&a.gp_circ2d()?, &cb);
            swap = false;
        }
        (Analytic2d::Elips, _) => {
            ana.perform_elips_conic(&a.gp_elips2d()?, &cb);
            swap = false;
        }
        (Analytic2d::Parab, _) => {
            ana.perform_parab_conic(&a.gp_parab2d()?, &cb);
            swap = false;
        }
        (Analytic2d::Hypr, _) => {
            ana.perform_hypr_conic(&a.gp_hypr2d()?, &cb);
            swap = false;
        }
    }

    if !ana.is_done() {
        return Some(Vec::new());
    }
    let mut out: Vec<(f64, f64, GpPnt2d)> = Vec::new();
    for i in 1..=ana.nb_points() {
        let p = ana.point(i);
        let (u1, u2) = (p.param_on_first(), p.param_on_second());
        let (u, v) = if swap { (u2, u1) } else { (u1, u2) };
        if param_on_curve2d(a, u, tol) && param_on_curve2d(b, v, tol) {
            out.push((u, v, *p.value()));
        }
    }
    Some(out)
}

/// Exact intersection of two infinite lines (cross-product formula).
/// `None` when the lines are parallel.
pub fn lin2d_intersection(l1: &GpLin2d, l2: &GpLin2d) -> Option<GpPnt2d> {
    let d1 = l1.pos.vdir;
    let d2 = l2.pos.vdir;
    let denom = d1.crossed(&d2);
    if denom.abs() < 1e-15 {
        return None;
    }
    let w = l2.pos.loc.xy().subtracted(&l1.pos.loc.xy());
    let d2v = GpXY::new(d2.x, d2.y);
    let t = w.crossed(&d2v) / denom;
    Some(GpPnt2d::new(l1.pos.loc.x() + t * d1.x, l1.pos.loc.y() + t * d1.y))
}

/// Exact intersection of segments `p1-p2` and `p3-p4`.
/// `None` when parallel or the intersection falls outside either segment.
pub fn segment_intersection(p1: &GpPnt2d, p2: &GpPnt2d, p3: &GpPnt2d, p4: &GpPnt2d) -> Option<GpPnt2d> {
    let d1 = GpXY::new(p2.x() - p1.x(), p2.y() - p1.y());
    let d2 = GpXY::new(p4.x() - p3.x(), p4.y() - p3.y());
    let denom = d1.crossed(&d2);
    if denom.abs() < 1e-15 {
        return None;
    }
    let w = GpXY::new(p3.x() - p1.x(), p3.y() - p1.y());
    let t = w.crossed(&d2) / denom;
    let s = w.crossed(&d1) / denom;
    let eps = 1e-12;
    if t >= -eps && t <= 1.0 + eps && s >= -eps && s <= 1.0 + eps {
        Some(GpPnt2d::new(p1.x() + t * d1.x, p1.y() + t * d1.y))
    } else {
        None
    }
}

/// Total length of an open polyline.
pub fn polyline2d_length(pts: &[GpPnt2d]) -> f64 {
    pts.windows(2).map(|w| w[0].distance(&w[1])).sum()
}

/// Project `p` onto segment `a-b`. Returns the closest point and its
/// parameter `t` in `[0, 1]`.
pub fn project_point_on_segment(p: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d) -> (GpPnt2d, f64) {
    let abx = b.x() - a.x();
    let aby = b.y() - a.y();
    let len2 = abx * abx + aby * aby;
    let t = if len2 < 1e-30 {
        0.0
    } else {
        (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby) / len2).clamp(0.0, 1.0)
    };
    (GpPnt2d::new(a.x() + t * abx, a.y() + t * aby), t)
}

// --- internals -------------------------------------------------------------

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn minimize_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64) -> (f64, f64) {
    const GOLD: f64 = 0.6180339887498949;
    let mut a = lo;
    let mut b = hi;
    let mut c = b - GOLD * (b - a);
    let mut d = a + GOLD * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    while (b - a) > 1e-12 {
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

/// Approximate bounding box of a bounded curve by sampling: `(xmin, xmax, ymin, ymax)`.
fn curve_bbox(c: &dyn Curve2d, n: usize) -> Option<(f64, f64, f64, f64)> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let mut minx = f64::INFINITY;
    let mut maxx = f64::NEG_INFINITY;
    let mut miny = f64::INFINITY;
    let mut maxy = f64::NEG_INFINITY;
    for i in 0..=n {
        let q = c.d0(a + (b - a) * i as f64 / n as f64);
        minx = minx.min(q.x());
        maxx = maxx.max(q.x());
        miny = miny.min(q.y());
        maxy = maxy.max(q.y());
    }
    Some((minx, maxx, miny, maxy))
}

/// Sample a curve over a finite parameter window. Unbounded curves get a
/// window that comfortably covers the other curve's bounding box.
fn sample_curve(c: &dyn Curve2d, other: &dyn Curve2d, n: usize) -> Vec<(f64, GpPnt2d)> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    let (lo, hi) = if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        let span = match curve_bbox(other, 32) {
            Some((x0, x1, y0, y1)) => {
                let size = (x1 - x0).max(y1 - y0).max(1e-6);
                let speed = c.d1(0.0).1.magnitude().max(1e-30);
                (4.0 * size / speed).max(1.0)
            }
            None => 100.0,
        };
        (-span, span)
    };
    (0..=n)
        .map(|i| {
            let u = lo + (hi - lo) * i as f64 / n as f64;
            (u, c.d0(u))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::Geom2dCircle;
    use crate::line::Geom2dLine;
    use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d};

    #[test]
    fn circle_length() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let len = curve2d_length(&circle);
        assert!((len - 2.0 * std::f64::consts::PI).abs() < 1e-6, "len={len}");
    }

    #[test]
    fn circle_line_two_intersections() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let hits = curve2d_intersections(&circle as &dyn Curve2d, &line as &dyn Curve2d, 1e-6);
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        let mut xs: Vec<f64> = hits.iter().map(|(_, _, p)| p.x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] + 1.0).abs() < 1e-4, "x0={}", xs[0]);
        assert!((xs[1] - 1.0).abs() < 1e-4, "x1={}", xs[1]);
    }

    #[test]
    fn crossing_segments() {
        let a = GpPnt2d::new(0.0, 0.0);
        let b = GpPnt2d::new(2.0, 2.0);
        let c = GpPnt2d::new(0.0, 2.0);
        let d = GpPnt2d::new(2.0, 0.0);
        let p = segment_intersection(&a, &b, &c, &d).unwrap();
        assert!((p.x() - 1.0).abs() < 1e-12 && (p.y() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn parallel_lines_none() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 1.0), GpDir2d::new(1.0, 0.0).unwrap());
        assert!(lin2d_intersection(&l1, &l2).is_none());
    }

    #[test]
    fn line_intersection_point() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 1.0), GpDir2d::new(0.0, 1.0).unwrap());
        let p = lin2d_intersection(&l1, &l2).unwrap();
        assert!((p.x() - 0.0).abs() < 1e-12 && (p.y() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn closest_point_on_line_is_foot() {
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let p = GpPnt2d::new(3.0, 4.0);
        let (u, foot) = curve2d_closest_point(&line as &dyn Curve2d, &p, 1e-6).unwrap();
        assert!((u - 3.0).abs() < 1e-6, "u={u}");
        assert!((foot.x() - 3.0).abs() < 1e-6 && (foot.y() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn project_on_segment_clamps() {
        let p = GpPnt2d::new(3.0, 4.0);
        let a = GpPnt2d::new(0.0, 0.0);
        let b = GpPnt2d::new(10.0, 0.0);
        let (foot, t) = project_point_on_segment(&p, &a, &b);
        assert!((t - 0.3).abs() < 1e-12, "t={t}");
        assert!((foot.x() - 3.0).abs() < 1e-12 && (foot.y() - 0.0).abs() < 1e-12);
        // Point beyond the segment clamps to the end.
        let (foot2, t2) = project_point_on_segment(&GpPnt2d::new(20.0, 0.0), &a, &b);
        assert!((t2 - 1.0).abs() < 1e-12);
        assert!(foot2.distance(&b) < 1e-12);
    }

    #[test]
    fn polyline_length() {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(3.0, 4.0),
            GpPnt2d::new(3.0, 6.0),
        ];
        assert!((polyline2d_length(&pts) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn distance_to_point_on_circle() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let d = curve2d_distance_to_point(&circle as &dyn Curve2d, &GpPnt2d::new(3.0, 0.0), 1e-6);
        assert!((d - 2.0).abs() < 1e-6, "d={d}");
    }

    #[test]
    fn closest_point_on_circle() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let (_, q) = curve2d_closest_point(&circle as &dyn Curve2d, &GpPnt2d::new(3.0, 0.0), 1e-6).unwrap();
        assert!((q.x() - 1.0).abs() < 1e-6 && (q.y() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn polyline_length_empty() {
        assert_eq!(polyline2d_length(&[]), 0.0);
        assert_eq!(polyline2d_length(&[GpPnt2d::zero()]), 0.0);
    }
}
