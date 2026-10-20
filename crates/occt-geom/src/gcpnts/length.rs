use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::curve::Curve;

use super::abscissa::CpntsAbscissaPoint;

/// `math::GaussPointsMax()` (`math.cxx:24-27`).
const GAUSS_POINTS_MAX: usize = 61;

/// `GeomAbs_CN` (`GeomAbs_Shape.hxx:29-45`), the continuity `computeType`
/// and `Intervals` use (`GCPnts_AbscissaPoint.cxx:28`, `:97`, `:99`).
const GEOM_ABS_CN: u8 = 6;

/// `GCPnts_AbscissaType` (`GCPnts_AbscissaType.hxx:22-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum AbscissaType {
    LengthParametrized,
    Parametrized,
    AbsComposite,
}

/// `computeType` (`GCPnts_AbscissaPoint.cxx:25-65`): the curve's abscissa
/// type plus the length ratio `theRatio` the `GCPnts_LengthParametrized`
/// arm divides by (`:36` line, `:40` circle radius, `:47`/`:56` the first
/// derivative magnitude of a two-pole curve).
pub(super) fn compute_type(c: &dyn Curve) -> (AbscissaType, f64) {
    if c.nb_intervals(GEOM_ABS_CN) > 1 {
        return (AbscissaType::AbsComposite, 1.0);
    }
    if c.is_line() {
        // `cxx:36`: `aRatio = 1.0` on a `Geom_Line`, whose parameter *is* the
        // arc length. Translation boundary: the port's trimmed views
        // (`GeomTrimmedCurve`) reparametrise onto `[0, 1]` with the chain rule
        // already baked into `d1`, whereas OCCT's `GeomAdaptor_Curve` keeps
        // the basis parameter (`GeomAdaptor_Curve.cxx:252-254`) so its `|d1|`
        // is exactly `1.0`. Taking the view's constant speed `|d1|` reproduces
        // `1.0` on the adaptor view and stays correct under the `[0, 1]`
        // reparametrisation.
        return (
            AbscissaType::LengthParametrized,
            c.d1(c.first_parameter()).1.magnitude(),
        );
    }
    if c.circle_radius().is_some() {
        // `cxx:40`: `aRatio = R` on the adaptor view (angle parameter) — the
        // same view-vs-adaptor equivalence as the line arm above: the circle's
        // speed is constant, so `|d1(FirstParameter)|` is `R` there and
        // `R * span` on a `[0, 1]`-reparametrised view.
        return (
            AbscissaType::LengthParametrized,
            c.d1(c.first_parameter()).1.magnitude(),
        );
    }
    // `GeomAbs_BezierCurve` (`cxx:43-51`): a two-pole non-rational Bezier is
    // a straight segment. This port's `GeomBezierCurve` carries no weights
    // (see its module note), so `IsRational()` is always false and only the
    // pole count decides.
    if let Some(p) = c.bezier_poles() {
        return if p.len() == 2 {
            (AbscissaType::LengthParametrized, c.d1(c.first_parameter()).1.magnitude())
        } else {
            (AbscissaType::Parametrized, 1.0)
        };
    }
    // `GeomAbs_BSplineCurve` (`cxx:52-60`).
    if let Some(p) = c.bspline_poles() {
        return if p.len() == 2 && c.bspline_weights().is_none() {
            (AbscissaType::LengthParametrized, c.d1(c.first_parameter()).1.magnitude())
        } else {
            (AbscissaType::Parametrized, 1.0)
        };
    }
    (AbscissaType::Parametrized, 1.0)
}

/// `GCPnts_AbsComposite` (`GCPnts_AbscissaPoint.cxx:96-158`): walk the
/// `GeomAbs_CN` intervals accumulating `CPnts_AbscissaPoint::Length`, then
/// solve inside the interval that contains the target. The final
/// "push a little bit outside the limits" bracket is `:153-156`.
pub(super) fn compute_abs_composite(
    c: &dyn Curve,
    abscissa: f64,
    u0: f64,
    ui: f64,
    resolution: f64,
) -> Option<f64> {
    let nb = c.nb_intervals(GEOM_ABS_CN);
    let ti = c.parameter_intervals(GEOM_ABS_CN);
    if nb < 1 || ti.len() < 2 {
        return None;
    }
    let mut abscis = abscissa;
    let mut u_start = u0;
    let mut ui = ui;
    let mut sign = 1.0;
    // `BSplCLib::Hunt(aTI, theU0, anIndex)` (`cxx:102`) yields a 1-based
    // interval index; `knots::hunt` is the 0-based form of the same search.
    let mut index = occt_core::bspl::knots::hunt(&ti, u0) as i32 + 1;
    let mut direction = 1i32;
    if abscis < 0.0 {
        direction = 0;
        abscis = -abscis;
        sign = -1.0;
    }
    while index >= 1 && index <= nb {
        let l = cpnts_length(c, u_start, ti[(index + direction) as usize - 1]);
        if (l - abscis).abs() <= CONFUSION {
            return Some(ti[(index + direction) as usize - 1]);
        }
        if l > abscis {
            if ui < ti[(index - 1) as usize] || ui > ti[index as usize] {
                let du = (abscis / l) * (ti[index as usize] - u_start);
                ui = if direction != 0 { u_start + du } else { u_start - du };
            }
            let mut computer = CpntsAbscissaPoint::new(c);
            computer.init_range(c, ti[(index - 1) as usize], ti[index as usize]);
            computer.perform_with_guess(sign * abscis, u_start, ui, resolution);
            return if computer.is_done() { Some(computer.parameter()) } else { None };
        }
        u_start = ti[(index + direction) as usize - 1];
        abscis -= l;
        index += if direction != 0 { 1 } else { -1 };
    }
    // `cxx:153-156`.
    ui = u_start + 0.1;
    let mut computer = CpntsAbscissaPoint::new(c);
    computer.init_range(c, u_start, u_start + 0.2);
    computer.perform_with_guess(sign * abscis, u_start, ui, resolution);
    if computer.is_done() { Some(computer.parameter()) } else { None }
}

pub(super) fn speed(c: &dyn Curve, u: f64) -> f64 {
    let (_, d1) = c.d1(u);
    d1.magnitude()
}

/// `CPnts_AbscissaPoint.cxx:57-77` (`order`): `f3d` integrated with a Gauss rule
/// whose order depends on the curve type — `Line` 2, `Parabola` 5,
/// `BezierCurve` `min(24, 2*Degree)`, `BSplineCurve` `min(24, 2*NbPoles - 1)`,
/// everything else 10.
pub(super) fn gauss_order(c: &dyn Curve) -> usize {
    if c.is_line() {
        2
    } else if c.gp_parabola().is_some() {
        5
    } else if let Some(p) = c.bezier_poles() {
        (2 * p.len().saturating_sub(1)).min(24)
    } else if let Some(p) = c.bspline_poles() {
        (2 * p.len()).saturating_sub(1).min(24)
    } else {
        10
    }
}

/// `math_GaussSingleIntegration::Perform` (`math_GaussSingleIntegration.cxx:100-150`):
/// scale the `[-1, 1]` rule onto `[Lower, Upper]` and sum `Weight * F(Point)`.
/// OCCT sums the symmetric pairs explicitly and scales by `xr` at the end; the
/// port uses the already-scaled pairs from `gauss_legendre`, which differ only in
/// the floating-point accumulation order.
pub(super) fn gauss_single(f: &dyn Fn(f64) -> f64, lower: f64, upper: f64, order: usize) -> f64 {
    let order = order.clamp(1, GAUSS_POINTS_MAX);
    let (points, weights) = occt_math::gauss::gauss_legendre(lower, upper, order);
    let mut val = 0.0;
    for (p, w) in points.iter().zip(weights.iter()) {
        val += w * f(*p);
    }
    val
}

/// `math_GaussSingleIntegration` with a tolerance
/// (`math_GaussSingleIntegration.cxx:64-98`): repeat the rule on `2^k` equal
/// sub-intervals (`IterMax = 13`) until two successive totals differ by at most
/// `tol`.
pub(super) fn gauss_single_tol(
    f: &dyn Fn(f64) -> f64,
    lower: f64,
    upper: f64,
    order: usize,
    tol: f64,
) -> f64 {
    const ITER_MAX: usize = 13;
    let mut len = gauss_single(f, lower, upper, order);
    let mut nb_interval = 1usize;
    for _ in 1..ITER_MAX {
        let old_len = len;
        len = 0.0;
        nb_interval *= 2;
        let du = (upper - lower) / nb_interval as f64;
        for i in 0..nb_interval {
            len += gauss_single(f, lower + i as f64 * du, lower + (i + 1) as f64 * du, order);
        }
        if (old_len - len).abs() <= tol {
            break;
        }
    }
    len
}

/// Arc length of `c` over `[a, b]` (`CPnts_AbscissaPoint::Length(C, U1, U2,
/// Tol)`, `CPnts_AbscissaPoint.cxx:168-184`): `|math_GaussSingleIntegration(…)|
/// ` on the speed, with the per-type Gauss order.
pub fn curve_length_range(c: &dyn Curve, a: f64, b: f64, tol: f64) -> f64 {
    if !(a.is_finite() && b.is_finite()) {
        return f64::INFINITY;
    }
    if b <= a {
        return 0.0;
    }
    let order = gauss_order(c);
    let speed_fn = |u: f64| speed(c, u);
    gauss_single_tol(&speed_fn, a, b, order, tol).abs()
}

/// `GCPnts_AbscissaPoint::length` (`GCPnts_AbscissaPoint.cxx:371-423`), the
/// type dispatch all six public `Length` overloads (`:305-366`) forward to.
///
/// * `GCPnts_LengthParametrized` (`:381-383`): the parameter already *is* an
///   arc length up to the constant `aRatio` `computeType` returned (line,
///   circle, two-pole Bezier/BSpline), so the answer is `|U2 - U1| * aRatio`.
/// * `GCPnts_Parametrized` (`:384-387`): `theTol != nullptr` picks the
///   tolerance-aware `CPnts_AbscissaPoint::Length(C, U1, U2, Tol)`, otherwise
///   the fixed-order `CPnts_AbscissaPoint::Length(C, U1, U2)`.
/// * `GCPnts_AbsComposite` (`:388-420`): sum `CPnts_AbscissaPoint::Length`
///   over the `GeomAbs_CN` intervals, using only the part of each interval
///   that overlaps `[min(U1, U2), max(U1, U2)]`.
///
/// `tol == None` is OCCT's `nullptr` (`:337`, `:345`), i.e. the three-argument
/// `Length` overloads; `Some(t)` is the four-argument form.
pub fn curve_length_between(c: &dyn Curve, u1: f64, u2: f64, tol: Option<f64>) -> f64 {
    match compute_type(c) {
        (AbscissaType::LengthParametrized, ratio) => (u2 - u1).abs() * ratio,
        (AbscissaType::Parametrized, _) => match tol {
            Some(t) => curve_length_range(c, u1, u2, t),
            None => cpnts_length(c, u1, u2),
        },
        (AbscissaType::AbsComposite, _) => {
            let nb = c.nb_intervals(GEOM_ABS_CN);
            let ti = c.parameter_intervals(GEOM_ABS_CN);
            if nb < 1 || ti.len() < nb as usize + 1 {
                return 0.0;
            }
            let uu1 = u1.min(u2);
            let uu2 = u1.max(u2);
            let mut total = 0.0;
            for index in 0..nb as usize {
                // `aTI(anIndex) > aUU2` / `aTI(anIndex + 1) < aUU1`
                // (`cxx:398`, `:402`), 1-based there, 0-based here.
                if ti[index] > uu2 {
                    break;
                }
                if ti[index + 1] < uu1 {
                    continue;
                }
                let lo = ti[index].max(uu1);
                let hi = ti[index + 1].min(uu2);
                total += match tol {
                    Some(t) => curve_length_range(c, lo, hi, t),
                    None => cpnts_length(c, lo, hi),
                };
            }
            total
        }
    }
}

/// `GCPnts_AbscissaPoint::Length(C)` (`GCPnts_AbscissaPoint.cxx:305-308`):
/// the dispatch over the curve's own first/last parameter, no tolerance.
pub fn curve_length(c: &dyn Curve) -> f64 {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() {
        curve_length_between(c, a, b, None)
    } else {
        f64::INFINITY
    }
}

/// `CPnts_AbscissaPoint::Length(C, U1, U2)` (`CPnts_AbscissaPoint.cxx:148-161`):
/// `|math_GaussSingleIntegration(f3d, U1, U2, order(C))|` — the no-tolerance
/// overload, which is the one `CPnts_AbscissaPoint::Init` uses for `myL`
/// (`:307`).
pub(super) fn cpnts_length(c: &dyn Curve, u1: f64, u2: f64) -> f64 {
    let f = |u: f64| speed(c, u);
    gauss_single(&f, u1, u2, gauss_order(c)).abs()
}

/// `CPnts_AbscissaPoint::Length(C, U1, U2, Tol)` (`CPnts_AbscissaPoint.cxx:168-184`):
/// the same integral with the adaptive
/// `math_GaussSingleIntegration(..., Tol)` rule. `AdvCompute` uses this
/// overload for every step of the `GCPnts_AbsComposite` walk
/// (`GCPnts_AbscissaPoint.cxx:206`, `:232`).
pub(super) fn cpnts_length_tol(c: &dyn Curve, u1: f64, u2: f64, tol: f64) -> f64 {
    let f = |u: f64| speed(c, u);
    gauss_single_tol(&f, u1, u2, gauss_order(c), tol).abs()
}

/// `GCPnts_AbsComposite` on the `AdvCompute` path
/// (`GCPnts_AbscissaPoint.cxx:187-295`), the tolerance-aware twin of
/// [`compute_abs_composite`]. Differences from the plain `Compute` arm
/// (`:96-158`), all per the `.cxx`:
///
/// * every `CPnts_AbscissaPoint::Length` step takes `theEPSILON` (`:206`, `:232`);
/// * the `anIndex == 0 && aDirection > 0` preamble (`:202-229`) exists only here;
/// * the root search uses `Init(theC, U1, U2, theEPSILON)` + `AdvPerform`
///   (`:219-220`, `:254-255`, `:291-292`);
/// * the zero-length test inside the loop is `Precision::PConfusion()`
///   (`:232`) while the preamble's is `theEPSILON` (`:206`);
/// * the "push a little bit outside the limits" tail (`:273-293`) offsets both
///   the guess and the bracket by `aSign` and clamps them to the first/last
///   parameter for a non-periodic curve.
pub(super) fn compute_abs_composite_adv(
    c: &dyn Curve,
    abscissa: f64,
    u0: f64,
    ui: f64,
    tol: f64,
) -> Option<f64> {
    let nb = c.nb_intervals(GEOM_ABS_CN);
    let ti = c.parameter_intervals(GEOM_ABS_CN);
    if nb < 1 || ti.len() < 2 {
        return None;
    }
    // `NCollection_Array1<double> aTI(1, aNbIntervals + 1)` is 1-based.
    let at = |i: i32| ti[(i - 1) as usize];
    let mut abscis = abscissa;
    let mut u_start = u0;
    let mut ui = ui;
    let mut sign = 1.0;
    let mut index = occt_core::bspl::knots::hunt(&ti, u0) as i32 + 1;
    let mut direction = 1i32;
    if abscis < 0.0 {
        direction = 0;
        abscis = -abscis;
        sign = -1.0;
    }

    // `cxx:202-229`: `Hunt` returned 0 (the parameter is below the first
    // knot), so the first step is bracketed by `aTI(1)`.
    if index == 0 && direction > 0 {
        let l = cpnts_length_tol(c, u_start, at(index + direction), tol);
        if (l - abscis).abs() <= tol {
            return Some(at(index + direction));
        }
        if l > abscis {
            if ui > at(index + 1) {
                ui = u_start + (abscis / l) * (at(index + 1) - u_start);
            }
            let mut computer = CpntsAbscissaPoint::new(c);
            computer.init_range_tol(c, u_start, at(index + 1), tol);
            computer.adv_perform(sign * abscis, u_start, ui, tol);
            return if computer.is_done() { Some(computer.parameter()) } else { None };
        }
        u_start = at(index + direction);
        abscis -= l;
        index += 1;
    }

    while index >= 1 && index <= nb {
        let l = cpnts_length_tol(c, u_start, at(index + direction), tol);
        if (l - abscis).abs() <= PCONFUSION {
            return Some(at(index + direction));
        }
        if l > abscis {
            if ui < at(index) || ui > at(index + 1) {
                let du = (abscis / l) * (at(index + 1) - u_start);
                ui = if direction != 0 { u_start + du } else { u_start - du };
            }
            let mut computer = CpntsAbscissaPoint::new(c);
            computer.init_range_tol(c, at(index), at(index + 1), tol);
            computer.adv_perform(sign * abscis, u_start, ui, tol);
            return if computer.is_done() { Some(computer.parameter()) } else { None };
        }
        u_start = at(index + direction);
        abscis -= l;
        index += if direction != 0 { 1 } else { -1 };
    }

    // `cxx:273-293`: "Push a little bit outside the limits (hairy !!!)".
    let mut ui = u_start + sign * 0.1;
    let mut u1 = u_start + sign * 0.2;
    if !c.is_periodic() {
        let (first, last) = (c.first_parameter(), c.last_parameter());
        if sign > 0.0 {
            ui = ui.min(last);
            u1 = u1.min(last);
        } else {
            ui = ui.max(first);
            u1 = u1.max(first);
        }
    }
    let mut computer = CpntsAbscissaPoint::new(c);
    computer.init_range_tol(c, u_start, u1, tol);
    computer.adv_perform(sign * abscis, u_start, ui, tol);
    if computer.is_done() { Some(computer.parameter()) } else { None }
}
