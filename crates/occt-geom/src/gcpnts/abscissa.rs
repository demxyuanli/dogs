use occt_core::precision::CONFUSION;
use occt_math::function_set_root::{FunctionWithDerivative, MathFunctionRoot};

use crate::curve::Curve;

use super::length::{
    compute_abs_composite, compute_abs_composite_adv, compute_type, cpnts_length,
    curve_length_range, gauss_order, gauss_single, gauss_single_tol, speed, AbscissaType,
};

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
    /// `myTol` (`CPnts_MyRootFunction.hxx`): `-1` from `Init(X0, L)` "to
    /// suppress the tolerance" (`:29`); `Tol` from `Init(X0, L, Tol)`
    /// (`:32-37`). `Value`/`Values` pick the tol-aware integration only when
    /// `myTol > 0` (`:43-50`, `:72-79`).
    tol: f64,
}

impl<'a> CpntsMyRootFunction<'a> {
    fn new(c: &'a dyn Curve, order: usize) -> Self {
        Self {
            c,
            order,
            x0: 0.0,
            l: 0.0,
            tol: -1.0,
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
        self.tol = -1.0;
    }

    /// `Init(const double X0, const double L, const double Tol)`
    /// (`CPnts_MyRootFunction.cxx:32-37`) -- the tolerance overload
    /// `AdvPerform` uses with `Resolution / 10` (`CPnts_AbscissaPoint.cxx:453`).
    fn init_tol(&mut self, x0: f64, l: f64, tol: f64) {
        self.x0 = x0;
        self.l = l;
        self.tol = tol;
    }

    /// `math_GaussSingleIntegration(myFunction, myX0, X, myOrder)`
    /// (`:41-45`, `math_GaussSingleIntegration.cxx:55-62`).
    fn integral(&self, x: f64) -> f64 {
        let f = |u: f64| speed(self.c, u);
        // `CPnts_MyRootFunction::Value` (`:43-50`): `myTol <= 0` keeps
        // `math_GaussSingleIntegration`; a positive `myTol` uses the
        // adaptive form `math_GaussSingleIntegration(..., myTol)`.
        if self.tol <= 0.0 {
            gauss_single(&f, self.x0, x, self.order)
        } else {
            gauss_single_tol(&f, self.x0, x, self.order, self.tol)
        }
    }
}

impl FunctionWithDerivative for CpntsMyRootFunction<'_> {
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
pub(super) struct CpntsAbscissaPoint<'a> {
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
    pub(super) fn new(c: &'a dyn Curve) -> Self {
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
    pub(super) fn init_range(&mut self, c: &'a dyn Curve, u1: f64, u2: f64) {
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

    /// `Init(C, U1, U2, Tol)` (`CPnts_AbscissaPoint.cxx:336-351`): the same
    /// as [`Self::init_range`] except `myL = Length(C, U1, U2, Tol)` -- the
    /// tol-aware integration (`CPnts_AbscissaPoint.cxx:345`).
    pub(super) fn init_range_tol(&mut self, c: &'a dyn Curve, u1: f64, u2: f64, tol: f64) {
        self.f = CpntsMyRootFunction::new(c, gauss_order(c));
        self.l = curve_length_range(c, u1, u2, tol);
        let mut umin = u1.min(u2);
        let mut umax = u1.max(u2);
        let du = umax - umin;
        umin -= du;
        umax += du;
        self.umin = umin;
        self.umax = umax;
    }

    /// `AdvPerform(Abscissa, U0, Ui, Resolution)` (`:436-...`): identical to
    /// [`Self::perform_with_guess`] except that the root function is given
    /// the tolerance `Resolution / 10` (rbv's modification, `:453`).
    pub(super) fn adv_perform(&mut self, abscissa: f64, u0: f64, ui: f64, resolution: f64) {
        if self.l < CONFUSION {
            self.done = true;
            self.param = u0;
            return;
        }
        self.done = false;
        self.f.init_tol(u0, abscissa, resolution / 10.0);
        let (umin, umax) = (self.umin, self.umax);
        let solution =
            MathFunctionRoot::new_with_bounds(&mut self.f, ui, resolution, umin, umax, 100);
        if solution.is_done() {
            self.done = true;
            self.param = solution.root();
        }
    }

    /// `Perform(Abscissa, U0, Resolution)` (`:374-391`): the guess
    /// `Ui = U0 + (Abscissa / myL) * (myUMax - myUMin) / 3` ("exercise : why
    /// 3 ?", `:387-388` — the `/3` undoes the `DU` widening on both sides),
    /// then the 4-argument `Perform`.
    pub(super) fn perform(&mut self, abscissa: f64, u0: f64, resolution: f64) {
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
    pub(super) fn perform_with_guess(
        &mut self,
        abscissa: f64,
        u0: f64,
        ui: f64,
        resolution: f64,
    ) {
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
    pub(super) fn is_done(&self) -> bool {
        self.done
    }

    /// `Parameter()` (`CPnts_AbscissaPoint.lxx:26-30`): OCCT raises
    /// `StdFail_NotDone` when `!myDone`; the caller checks `is_done()` first
    /// (the port's `Parameter` does not hard-crash on misuse).
    pub(super) fn parameter(&self) -> f64 {
        self.param
    }
}

/// `GCPnts_AbscissaPoint::Compute` on the `GCPnts_Parametrized` arm
/// (`GCPnts_AbscissaPoint.cxx:77-95`): the `Precision::Confusion()` shortcut
/// (`:77-81`) then `Init(theC)` + `Perform(theAbscis, theU0, theUi,
/// theEPSILON)` (`:91-95`). `None` stands for OCCT's `!IsDone()`, whose
/// `Parameter()` raises `StdFail_NotDone`.
///
/// All three `GCPnts_AbscissaType` arms are implemented below: the
/// `GCPnts_LengthParametrized` shortcut (`:87-90`), the `GCPnts_AbsComposite`
/// walk (`:96-158`) and the `GCPnts_Parametrized` root search (`:91-95`); the
/// dispatch is `compute_type` (`:26-65`).
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

/// `GCPnts_AbscissaPoint(C, Abscissa, U0, Ui, Tol)`
/// (`GCPnts_AbscissaPoint.cxx:533-541` → `AdvCompute`). Used by
/// `Approx_CurvlinFunc::GetUParameter` (`Approx_CurvlinFunc.cxx:560`).
pub fn abscissa_point_from_guess(
    c: &dyn Curve,
    abscissa: f64,
    u0: f64,
    ui: f64,
    tol: f64,
) -> Option<f64> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    if abscissa.abs() <= CONFUSION {
        return Some(u0);
    }
    let (ty, ratio) = compute_type(c);
    match ty {
        AbscissaType::LengthParametrized => return Some(u0 + abscissa / ratio),
        // `AdvCompute` `:187-295`: the tolerance-aware composite walk.
        AbscissaType::AbsComposite => {
            return compute_abs_composite_adv(c, abscissa, u0, ui, tol);
        }
        AbscissaType::Parametrized => {}
    }
    let mut computer = CpntsAbscissaPoint::new(c);
    computer.init_range_tol(c, a, b, tol);
    // `AdvCompute` passes `theTol` as the `theEPSILON` argument
    // (`GCPnts_AbscissaPoint.cxx:533-541`), so the root search and the
    // `myF.Init(U0, Abscissa, Resolution / 10)` tolerance are `tol`, not the
    // curve's `Resolution(Precision::Confusion())`.
    computer.adv_perform(abscissa, u0, ui, tol);
    if computer.is_done() {
        Some(computer.parameter())
    } else {
        None
    }
}

/// Same as [`abscissa_point`] but with the caller's tolerance -- the
/// `GCPnts_AbscissaPoint(theTol, theC, theAbscissa, theU0)` entry, which
/// reaches `advCompute` (`GCPnts_AbscissaPoint.cxx:465-487`, `:495-497`) and
/// therefore `CPnts_AbscissaPoint::Init(C, Tol)` (`:287-290`) +
/// `AdvPerform` (`CPnts_AbscissaPoint.cxx:436-474`). The tol-aware
/// integration replaces the fixed-order Gauss rule in the root function
/// (`CPnts_MyRootFunction.cxx:43-50`), and every branch of
/// `AdvCompute` (`GCPnts_AbscissaPoint.cxx:166-295`) is ported, including the
/// `GCPnts_AbsComposite` arm with its `anIndex == 0` preamble.
pub fn abscissa_point_with_tolerance(
    c: &dyn Curve,
    abscissa: f64,
    from: f64,
    tol: f64,
) -> Result<f64, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("abscissa_point_with_tolerance: unbounded curve".to_string());
    }
    // `advCompute` (`:470-481`): `aL = Length(theC, theTol)` (the
    // `Standard_ConstructionError` throw is commented out there).
    let al = curve_length_range(c, a, b, tol);
    let ui = if al >= CONFUSION {
        from + (abscissa / al) * (b - a)
    } else {
        0.0
    };
    if abscissa.abs() <= CONFUSION {
        return Ok(from);
    }
    let (ty, ratio) = compute_type(c);
    match ty {
        // `AdvCompute` `:177-180`.
        AbscissaType::LengthParametrized => return Ok(from + abscissa / ratio),
        // `AdvCompute` `:187-295`.
        AbscissaType::AbsComposite => {
            return compute_abs_composite_adv(c, abscissa, from, ui, tol)
                .ok_or_else(|| {
                    "abscissa_point_with_tolerance: math_FunctionRoot did not converge"
                        .to_string()
                });
        }
        // `:183-184`: `Init(theC, theEPSILON)` + `AdvPerform`.
        AbscissaType::Parametrized => {}
    }
    let mut computer = CpntsAbscissaPoint::new(c);
    computer.init_range_tol(c, a, b, tol);
    // `AdvCompute` passes `theTol` as the `theEPSILON` argument
    // (`GCPnts_AbscissaPoint.cxx:486`), so the root search tolerance is `tol`.
    computer.adv_perform(abscissa, from, ui, tol);
    if computer.is_done() {
        Ok(computer.parameter())
    } else {
        Err("abscissa_point_with_tolerance: math_FunctionRoot did not converge".to_string())
    }
}
