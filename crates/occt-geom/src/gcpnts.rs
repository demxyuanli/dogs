//! Curve discretization points. Port of the `GCPnts`/`CPnts` packages
//! (TKGeomBase): `CPnts_AbscissaPoint`, `GCPnts_AbscissaPoint`,
//! `GCPnts_UniformAbscissa`, `GCPnts_QuasiUniformAbscissa`, plus a parameter
//! accessor for the tangential-deflection sampler.
//!
//! Arc length is the faithful `CPnts_AbscissaPoint::Length` path: the speed
//! `|C'(u)|` is integrated with `math_GaussSingleIntegration` at OCCT's
//! per-type order (`CPnts_AbscissaPoint.cxx:57-77`) and, for the tolerance
//! overload, with the 13-iteration interval-doubling loop
//! (`math_GaussSingleIntegration.cxx:64-98`). OCCT tabulates the Gauss nodes and
//! weights (`math.cxx` `Point[]`/`Weight[]`, `GaussPointsMax() = 61`); the port
//! computes the same values with [`occt_math::gauss::gauss_legendre`] and clamps
//! the order to 61 the same way.
//!
//! The abscissa *inversion* (`GCPnts_AbscissaPoint::Parameter`) is the faithful
//! `CPnts_AbscissaPoint::Init`/`Perform` path
//! (`CPnts_AbscissaPoint.cxx:270-313`, `:374-432`) driven by
//! `math_FunctionRoot` on `CPnts_MyRootFunction` / `CPnts_MyGaussFunction`
//! (`CPnts_MyRootFunction.cxx:19-90`, `CPnts_MyGaussFunction.cxx:17-27`);
//! `math_FunctionRoot` itself runs `math_FunctionSetRoot`
//! (`math_FunctionRoot.cxx:73-119`).
//!
//! **UNPORTED (A15/T-51 remainder)**: `GCPnts_AbscissaPoint::Compute`'s type
//! dispatch (`GCPnts_AbscissaPoint.cxx:26-65`, `:67-161`) is not reproduced —
//! only its `GCPnts_Parametrized` arm (`:91-95`, reached from `compute`
//! `:428-442`) is. The `GCPnts_LengthParametrized` (`:87-90`) and
//! `GCPnts_AbsComposite` (`:96-158`) arms read `GeomAdaptor_Curve`'s
//! `GetType`/`NbIntervals`/`Intervals`; the port's `GeomTrimmedCurve` view
//! remaps a trim onto `[0, 1]` instead of unwrapping to the basis parameter
//! range the adaptor keeps (`GeomAdaptor_Curve.cxx:239-254`), so those arms
//! cannot be expressed faithfully here. `uniform_abscissa` /
//! `quasi_uniform_abscissa` keep their own outer shape — `uniform_abscissa`
//! is "n intervals" (`n + 1` points), not `GCPnts_UniformAbscissa`'s
//! `NbPoints` point count.
//!
//! `GCPnts_UniformDeflection` / `GCPnts_QuasiUniformDeflection` are **UNPORTED**
//! (see `occt-core/src/gcpnts.rs` module docs) — this module no longer exposes a
//! fake "uniform deflection" accessor.

use occt_core::gcpnts::{perform_tangential_curve, CurveSecondDeriv};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_math::function_set_root::{MathFunctionRoot, MathFunctionWithDerivative};

use crate::curve::Curve;

/// `math::GaussPointsMax()` (`math.cxx:24-27`).
const GAUSS_POINTS_MAX: usize = 61;

/// `GeomAbs_CN` (`GeomAbs_Shape.hxx:29-45`), the continuity `computeType`
/// and `Intervals` use (`GCPnts_AbscissaPoint.cxx:28`, `:97`, `:99`).
const GEOM_ABS_CN: u8 = 6;

/// `GCPnts_AbscissaType` (`GCPnts_AbscissaType.hxx:22-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AbscissaType {
    LengthParametrized,
    Parametrized,
    AbsComposite,
}

/// `computeType` (`GCPnts_AbscissaPoint.cxx:25-65`): the curve's abscissa
/// type plus the length ratio `theRatio` the `GCPnts_LengthParametrized`
/// arm divides by (`:36` line, `:40` circle radius, `:47`/`:56` the first
/// derivative magnitude of a two-pole curve).
fn compute_type(c: &dyn Curve) -> (AbscissaType, f64) {
    if c.nb_intervals(GEOM_ABS_CN) > 1 {
        return (AbscissaType::AbsComposite, 1.0);
    }
    if c.is_line() {
        return (AbscissaType::LengthParametrized, 1.0);
    }
    if let Some(r) = c.circle_radius() {
        return (AbscissaType::LengthParametrized, r);
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
fn compute_abs_composite(
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

fn speed(c: &dyn Curve, u: f64) -> f64 {
    let (_, d1) = c.d1(u);
    d1.magnitude()
}

/// `CPnts_AbscissaPoint.cxx:57-77` (`order`): `f3d` integrated with a Gauss rule
/// whose order depends on the curve type — `Line` 2, `Parabola` 5,
/// `BezierCurve` `min(24, 2*Degree)`, `BSplineCurve` `min(24, 2*NbPoles - 1)`,
/// everything else 10.
fn gauss_order(c: &dyn Curve) -> usize {
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
fn gauss_single(f: &dyn Fn(f64) -> f64, lower: f64, upper: f64, order: usize) -> f64 {
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
fn gauss_single_tol(f: &dyn Fn(f64) -> f64, lower: f64, upper: f64, order: usize, tol: f64) -> f64 {
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

/// Total arc length of `c` over its parameter range.
pub fn curve_length(c: &dyn Curve) -> f64 {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() {
        curve_length_range(c, a, b, CONFUSION * 0.1)
    } else {
        f64::INFINITY
    }
}

/// `CPnts_AbscissaPoint::Length(C, U1, U2)` (`CPnts_AbscissaPoint.cxx:148-161`):
/// `|math_GaussSingleIntegration(f3d, U1, U2, order(C))|` — the no-tolerance
/// overload, which is the one `CPnts_AbscissaPoint::Init` uses for `myL`
/// (`:307`).
fn cpnts_length(c: &dyn Curve, u1: f64, u2: f64) -> f64 {
    let f = |u: f64| speed(c, u);
    gauss_single(&f, u1, u2, gauss_order(c)).abs()
}

/// `CPnts_MyRootFunction` (`CPnts_MyRootFunction.hxx:32-64`,
/// `CPnts_MyRootFunction.cxx:19-90`): `Value(X) = Integral(X0, X, |C'|) - L`,
/// `Derivative(X) = |C'(X)|`, the `math_FunctionWithDerivative` that
/// `math_FunctionRoot` drives. Its `myFunction` is `CPnts_MyGaussFunction`
/// wrapping the `f3d` speed (`CPnts_AbscissaPoint.cxx:41-47`); `Derivative` is
/// `myFunction.Value` (`CPnts_MyRootFunction.cxx:63-66`).
struct CpntsMyRootFunction<'a> {
    c: &'a dyn Curve,
    /// `myOrder` set by `Init(F, D, Order)` (`:19-23`).
    order: usize,
    /// `myX0` / `myL` set by `Init(X0, L)` (`:25-30`).
    x0: f64,
    l: f64,
}

impl<'a> CpntsMyRootFunction<'a> {
    fn new(c: &'a dyn Curve, order: usize) -> Self {
        Self {
            c,
            order,
            x0: 0.0,
            l: 0.0,
        }
    }

    /// `Init(const double X0, const double L)` (`:25-30`): `myTol = -1` "to
    /// suppress the tolerance", so `Integral` uses the no-tolerance
    /// `math_GaussSingleIntegration`.
    ///
    /// UNPORTED: the tolerance overload `Init(X0, L, Tol)` (`:32-37`) and the
    /// `AdvPerform` path that uses it (`CPnts_AbscissaPoint.cxx:436-474`) are
    /// not exposed by this port's public API.
    fn init(&mut self, x0: f64, l: f64) {
        self.x0 = x0;
        self.l = l;
    }

    /// `math_GaussSingleIntegration(myFunction, myX0, X, myOrder)`
    /// (`:41-45`, `math_GaussSingleIntegration.cxx:55-62`).
    fn integral(&self, x: f64) -> f64 {
        let f = |u: f64| speed(self.c, u);
        gauss_single(&f, self.x0, x, self.order)
    }
}

impl MathFunctionWithDerivative for CpntsMyRootFunction<'_> {
    /// `CPnts_MyRootFunction::Value` (`:39-61`). `math_GaussSingleIntegration`
    /// fails only when the integrand's `Value` returns false
    /// (`math_GaussSingleIntegration.cxx:128-147`), and
    /// `CPnts_MyGaussFunction::Value` always returns true
    /// (`CPnts_MyGaussFunction.cxx:23-27`), so the port always succeeds here.
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.integral(x) - self.l;
        true
    }

    /// `CPnts_MyRootFunction::Derivative` (`:63-66`) = `myFunction.Value(X)`.
    fn derivative(&mut self, x: f64, df: &mut f64) -> bool {
        *df = speed(self.c, x);
        true
    }

    /// `CPnts_MyRootFunction::Values` (`:68-90`).
    fn values(&mut self, x: f64, f: &mut f64, df: &mut f64) -> bool {
        *f = self.integral(x) - self.l;
        *df = speed(self.c, x);
        true
    }
}

/// `CPnts_AbscissaPoint` (`CPnts_AbscissaPoint.hxx:35-198`): the arc-length
/// inversion solving `Integral(U0, X, |C'|) = Abscissa` with
/// `math_FunctionRoot` (`CPnts_AbscissaPoint.cxx:395-432`).
struct CpntsAbscissaPoint<'a> {
    done: bool,
    l: f64,
    param: f64,
    umin: f64,
    umax: f64,
    f: CpntsMyRootFunction<'a>,
}

impl<'a> CpntsAbscissaPoint<'a> {
    /// `CPnts_AbscissaPoint()` (`:211-218`): everything zero except
    /// `myDone = false`. `Init` replaces `myF` and fills the rest.
    fn new(c: &'a dyn Curve) -> Self {
        Self {
            done: false,
            l: 0.0,
            param: 0.0,
            umin: 0.0,
            umax: 0.0,
            f: CpntsMyRootFunction::new(c, gauss_order(c)),
        }
    }

    /// `Init(const Adaptor3d_Curve& C)` (`:270-273`).
    fn init(&mut self, c: &'a dyn Curve) {
        self.init_range(c, c.first_parameter(), c.last_parameter());
    }

    /// `Init(C, U1, U2)` (`:301-313`): `myF.Init(f3d, &C, order(C))`,
    /// `myL = Length(C, U1, U2)`, and the widened bracket
    /// `myUMin = min - DU`, `myUMax = max + DU` with `DU = max - min`
    /// (`:308-312`).
    fn init_range(&mut self, c: &'a dyn Curve, u1: f64, u2: f64) {
        self.f = CpntsMyRootFunction::new(c, gauss_order(c));
        self.l = cpnts_length(c, u1, u2);
        let mut umin = u1.min(u2);
        let mut umax = u1.max(u2);
        let du = umax - umin;
        umin -= du;
        umax += du;
        self.umin = umin;
        self.umax = umax;
    }

    /// `Perform(Abscissa, U0, Resolution)` (`:374-391`): the guess
    /// `Ui = U0 + (Abscissa / myL) * (myUMax - myUMin) / 3` ("exercise : why
    /// 3 ?", `:387-388` — the `/3` undoes the `DU` widening on both sides),
    /// then the 4-argument `Perform`.
    fn perform(&mut self, abscissa: f64, u0: f64, resolution: f64) {
        if self.l < CONFUSION {
            self.done = true;
            self.param = u0;
        } else {
            let ui = u0 + (abscissa / self.l) * (self.umax - self.umin) / 3.0;
            self.perform_with_guess(abscissa, u0, ui, resolution);
        }
    }

    /// `Perform(Abscissa, U0, Ui, Resolution)` (`:395-432`). The validity test
    /// on `Solution.Value()` / `Derivative` is commented out in OCCT
    /// (`:415-425`); only `Solution.IsDone()` is consulted (`:426-430`), so the
    /// port does the same.
    fn perform_with_guess(&mut self, abscissa: f64, u0: f64, ui: f64, resolution: f64) {
        if self.l < CONFUSION {
            self.done = true;
            self.param = u0;
            return;
        }
        self.done = false;
        self.f.init(u0, abscissa);
        let (umin, umax) = (self.umin, self.umax);
        let solution =
            MathFunctionRoot::new_with_bounds(&mut self.f, ui, resolution, umin, umax, 100);
        if solution.is_done() {
            self.done = true;
            self.param = solution.root();
        }
    }

    /// `IsDone()` (`CPnts_AbscissaPoint.lxx:19-22`).
    fn is_done(&self) -> bool {
        self.done
    }

    /// `Parameter()` (`CPnts_AbscissaPoint.lxx:26-30`): OCCT raises
    /// `StdFail_NotDone` when `!myDone`; the caller checks `is_done()` first
    /// (the port's `Parameter` does not hard-crash on misuse).
    fn parameter(&self) -> f64 {
        self.param
    }
}

/// `GCPnts_AbscissaPoint::Compute` on the `GCPnts_Parametrized` arm
/// (`GCPnts_AbscissaPoint.cxx:77-95`): the `Precision::Confusion()` shortcut
/// (`:77-81`) then `Init(theC)` + `Perform(theAbscis, theU0, theUi,
/// theEPSILON)` (`:91-95`). `None` stands for OCCT's `!IsDone()`, whose
/// `Parameter()` raises `StdFail_NotDone`.
///
/// UNPORTED: `GCPnts_LengthParametrized` (`:87-90`) and
/// `GCPnts_AbsComposite` (`:96-158`); see the module header.
fn compute_with_guess(
    c: &dyn Curve,
    abscissa: f64,
    u0: f64,
    ui: f64,
    resolution: f64,
) -> Option<f64> {
    if abscissa.abs() <= CONFUSION {
        return Some(u0);
    }
    // `GCPnts_AbscissaPoint.cxx:83-158`: `computeType` picks the arm.
    let (ty, ratio) = compute_type(c);
    match ty {
        // `:87-90`: the abscissa *is* the parameter up to the ratio.
        AbscissaType::LengthParametrized => return Some(u0 + abscissa / ratio),
        // `:96-158`.
        AbscissaType::AbsComposite => {
            return compute_abs_composite(c, abscissa, u0, ui, resolution);
        }
        // `:91-95`: the iterative root search below.
        AbscissaType::Parametrized => {}
    }
    let mut computer = CpntsAbscissaPoint::new(c);
    computer.init(c);
    computer.perform_with_guess(abscissa, u0, ui, resolution);
    if computer.is_done() {
        Some(computer.parameter())
    } else {
        None
    }
}

/// Parameter `u` at arc length `abscissa` from the point of parameter `from` —
/// the `GCPnts_AbscissaPoint::Parameter()` entry
/// (`GCPnts_AbscissaPoint.hxx:153`).
///
/// Positive `abscissa` walks forward, negative backward. Mirrors
/// `GCPnts_AbscissaPoint(theC, theAbscissa, theU0)`
/// (`GCPnts_AbscissaPoint.cxx:446-451` → `compute` `:428-442` →
/// `Compute` `:67-161`). `Err` covers OCCT's `Standard_ConstructionError`
/// for a zero-length curve (`:433-436`) and `StdFail_NotDone` when the root
/// search does not converge (`CPnts_AbscissaPoint.lxx:26-30`). As OCCT documents
/// (`CPnts_AbscissaPoint.hxx:78`), the result may lie outside the curve's
/// parameter bounds: the `Perform` bracket is widened by `DU` on each side
/// (`CPnts_AbscissaPoint.cxx:310-312`).
pub fn abscissa_point(c: &dyn Curve, abscissa: f64, from: f64) -> Result<f64, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        // Port guard: OCCT would integrate over an infinite range.
        return Err("abscissa_point: unbounded curve".to_string());
    }
    // `aL = Length(theC)` and `Standard_ConstructionError` below
    // `Precision::Confusion()` (`GCPnts_AbscissaPoint.cxx:432-436`).
    let al = cpnts_length(c, a, b);
    if al < CONFUSION {
        return Err("abscissa_point: zero-length curve".to_string());
    }
    // `aUUi = theU0 + (anAbscis / aL) * (Last - First)` (`:440`),
    // `theC.Resolution(Precision::Confusion())` (`:441`).
    let ui = from + (abscissa / al) * (b - a);
    let resolution = c.resolution(CONFUSION);
    compute_with_guess(c, abscissa, from, ui, resolution)
        .ok_or_else(|| "abscissa_point: math_FunctionRoot did not converge".to_string())
}

/// `n + 1` parameters at equal arc-length spacing across the curve
/// (`GCPnts_UniformAbscissa`). A zero-length curve returns `[a; n + 1]`.
pub fn uniform_abscissa(c: &dyn Curve, n: usize) -> Result<Vec<f64>, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("uniform_abscissa: unbounded curve".to_string());
    }
    let n = n.max(1);
    let total = curve_length_range(c, a, b, CONFUSION * 0.1);
    if total <= CONFUSION {
        // Degenerate: every parameter maps to the same point.
        return Ok(vec![a; n + 1]);
    }
    let step = total / n as f64;
    let resolution = c.resolution(CONFUSION);
    let mut params = Vec::with_capacity(n + 1);
    params.push(a);
    let mut prev = a;
    for _ in 1..n {
        // Find u with arc length `step` from `prev`, inside `[prev, b]`,
        // through the faithful `CPnts_AbscissaPoint::Init(C, U1, U2)` +
        // `Perform(Abscissa, U0, Resolution)` (`CPnts_AbscissaPoint.cxx:301-313`,
        // `:374-391`), i.e. `math_FunctionRoot` on the integral of the speed.
        let mut computer = CpntsAbscissaPoint::new(c);
        computer.init_range(c, prev, b);
        computer.perform(step, prev, resolution);
        if !computer.is_done() {
            return Err(
                "uniform_abscissa: math_FunctionRoot did not converge".to_string(),
            );
        }
        let u = computer.parameter().min(b);
        params.push(u);
        prev = u;
    }
    params.push(b);
    Ok(params)
}

/// `n` parameters distributed at equal *chord* intervals over a polyline
/// approximation of the curve (port of `GCPnts_QuasiUniformAbscissa`, which
/// uses a 2n-sample chord-length table).
pub fn quasi_uniform_abscissa(c: &dyn Curve, n: usize) -> Result<Vec<f64>, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("quasi_uniform_abscissa: unbounded curve".to_string());
    }
    let n = n.max(2);
    let total = curve_length_range(c, a, b, CONFUSION * 0.1);
    if total <= CONFUSION {
        return Ok(vec![a; n]);
    }
    // Chord-length vs parameter table over 2n samples (last lands on b).
    let samples = 2 * n;
    let du = (b - a) / (samples - 1) as f64;
    let mut cum = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    cum.push(0.0);
    params.push(a);
    let mut prev = c.d0(a);
    for k in 1..samples {
        let u = a + k as f64 * du;
        let p = c.d0(u);
        let l = cum[k - 1] + prev.distance(&p);
        cum.push(l);
        params.push(u);
        prev = p;
    }
    let total = *cum.last().unwrap();
    let dcorde = total / (n - 1) as f64;
    let mut out = Vec::with_capacity(n);
    out.push(a);
    let mut idx = 1usize;
    for i in 1..n - 1 {
        let target = dcorde * i as f64;
        while idx < samples - 1 && cum[idx] < target {
            idx += 1;
        }
        let denom = cum[idx] - cum[idx - 1];
        let alpha = if denom.abs() > 1e-30 {
            (target - cum[idx - 1]) / denom
        } else {
            0.0
        };
        out.push(params[idx - 1] + alpha * (params[idx] - params[idx - 1]));
    }
    out.push(b);
    Ok(out)
}

/// Adapter from `&dyn Curve` to the `GCPnts_TangentialDeflection` engine.
struct Adapter<'a>(&'a dyn Curve);

impl CurveSecondDeriv for Adapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

/// Parameters refined by both chord deviation `tol` and tangent-angle change
/// `angle_tol`, via the faithful `GCPnts_TangentialDeflection` engine
/// (`GCPnts_TangentialDeflection.cxx:522-916`, ported in
/// `occt-core/src/gcpnts_perform.rs`).
///
/// The CN breakpoint set and the BSpline/Bezier minimum-point bump follow the
/// same derivation as `meshing::edge_discret::CurveTessellator::initialize`
/// (`GCPnts_TangentialDeflection::initialize`, `cxx:415-453`); `u_tol`/`min_len`
/// are the OCCT `Initialize` defaults (`cxx:302-322`, `1.0e-9` / `CONFUSION`).
pub fn tangential_deflection(c: &dyn Curve, tol: f64, angle_tol: f64) -> Vec<f64> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return vec![a, b];
    }
    let ang = if angle_tol.is_finite() && angle_tol > 0.0 {
        angle_tol
    } else {
        std::f64::consts::PI
    };
    let mut intervals = c.parameter_intervals(6);
    if intervals.len() < 2 {
        intervals = vec![a, b];
    }
    let degree_min_nb = c.nurbs_degree().map(|d| (d + 1).max(2)).unwrap_or(2);
    let (params, _) = perform_tangential_curve(
        &Adapter(c),
        a,
        b,
        ang,
        tol,
        2,
        PCONFUSION,
        CONFUSION,
        &intervals,
        degree_min_nb,
    );
    if params.len() < 2 {
        vec![a, b]
    } else {
        params
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeomCircle, GeomLine, GeomTrimmedCurve};
    use occt_core::elib::clib;
    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt};
    use std::sync::Arc;

    const PI: f64 = std::f64::consts::PI;

    fn unit_circle() -> GeomCircle {
        GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))
    }

    fn unit_circle_arc() -> GeomTrimmedCurve {
        GeomTrimmedCurve::new(Arc::new(unit_circle()), 0.0, PI / 2.0)
    }

    #[test]
    fn length_unit_circle() {
        let c = unit_circle();
        let l = curve_length(&c);
        assert!((l - 2.0 * PI).abs() < 1e-7, "length {l}");
    }

    #[test]
    fn length_quarter_arc() {
        let c = unit_circle_arc();
        let l = curve_length(&c);
        assert!((l - PI / 2.0).abs() < 1e-7, "length {l}");
    }

    #[test]
    fn abscissa_point_half_circle() {
        let c = unit_circle();
        // Arc length π from u=0 lands at u=π.
        let u = abscissa_point(&c, PI, 0.0).unwrap();
        assert!((u - PI).abs() < 1e-6, "u {u}");
        let p = c.d0(u);
        assert!(p.distance(&GpPnt::new(-1.0, 0.0, 0.0)) < 1e-6, "point {p:?}");
    }

    #[test]
    fn abscissa_point_total_length_returns_last_parameter() {
        let c = unit_circle();
        let u = abscissa_point(&c, 2.0 * PI, 0.0).unwrap();
        assert!((u - 2.0 * PI).abs() < 1e-6, "u {u}");
        // `CPnts_AbscissaPoint.hxx:78`: "The computed point can be outside
        // of the curve 's bounds". `CPnts_AbscissaPoint::Perform` widens the
        // root bracket to `[U1 - DU, U2 + DU]` (`CPnts_AbscissaPoint.cxx:310-312`),
        // so an abscissa past the total length is not rejected: on the periodic
        // circle the integral `Length(0, X)` still matches it at `X = 3*pi`.
        let u = abscissa_point(&c, 3.0 * PI, 0.0).unwrap();
        assert!((u - 3.0 * PI).abs() < 1e-6, "u {u}");
    }

    #[test]
    fn abscissa_point_negative_goes_backward() {
        let c = unit_circle();
        // Backward by π/2 from u=π lands at u=π/2 (the point (0,1,0)).
        let u = abscissa_point(&c, -PI / 2.0, PI).unwrap();
        assert!((u - PI / 2.0).abs() < 1e-6, "u {u}");
        let p = c.d0(u);
        assert!(p.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6, "point {p:?} at u {u}");
    }

    #[test]
    fn uniform_abscissa_unit_circle_quarter_steps() {
        let c = unit_circle();
        let params = uniform_abscissa(&c, 4).unwrap();
        assert_eq!(params.len(), 5);
        let expected = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0, 2.0 * PI];
        for (g, e) in params.iter().zip(expected.iter()) {
            assert!((g - e).abs() < 1e-6, "got {g} expected {e}: {params:?}");
        }
        // Points are equally spaced on the circle.
        let pts: Vec<GpPnt> = params.iter().map(|&u| c.d0(u)).collect();
        for w in pts.windows(2) {
            let chord = w[0].distance(&w[1]);
            assert!((chord - (2.0f64).sqrt()).abs() < 1e-5, "chord {chord}");
        }
    }

    #[test]
    fn uniform_abscissa_line_is_uniform() {
        // A 10-unit trimmed line has parameter range [0, 1]; equal arc-length
        // spacing of 10/4 = 2.5 lands at params 0, 0.25, 0.5, 0.75, 1.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 10.0);
        let params = uniform_abscissa(&c, 4).unwrap();
        for (i, u) in params.iter().enumerate() {
            assert!((u - 0.25 * i as f64).abs() < 1e-9, "params {params:?}");
        }
    }

    #[test]
    fn uniform_abscissa_degenerate_zero_length() {
        let line = GeomLine::from_pnt_dir(GpPnt::new(1., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 0.0);
        let params = uniform_abscissa(&c, 3).unwrap();
        assert_eq!(params, vec![0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn quasi_uniform_abscissa_line() {
        // 10-unit trimmed line over [0, 1]: 5 points at equal chord intervals
        // 2.5 → params 0, 0.25, 0.5, 0.75, 1.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 10.0);
        let params = quasi_uniform_abscissa(&c, 5).unwrap();
        assert_eq!(params.len(), 5);
        assert!((params[0] - 0.0).abs() < 1e-9);
        assert!((params[4] - 1.0).abs() < 1e-9);
        assert!((params[1] - 0.25).abs() < 1e-6, "params {params:?}");
        assert!((params[2] - 0.5).abs() < 1e-6, "params {params:?}");
        assert!((params[3] - 0.75).abs() < 1e-6, "params {params:?}");
    }

    #[test]
    fn sampled_circle_points_consistent() {
        let c = unit_circle();
        let params = uniform_abscissa(&c, 8).unwrap();
        for &u in &params {
            let p = c.d0(u);
            let q = clib::circle_value(c.circ(), u);
            assert!(p.distance(&q) < 1e-12);
        }
    }
}
