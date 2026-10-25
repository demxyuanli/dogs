//! `Law_Function` (the root trait) with `Law_Constant`, `Law_Linear` and
//! `Law_S`.
//!
//! Source: `Law_Function.hxx`/`Law_Function.cxx`, `Law_Constant.hxx`/
//! `Law_Constant.cxx`, `Law_Linear.hxx`/`Law_Linear.cxx`, `Law_S.hxx`/
//! `Law_S.cxx`.
//!
//! Translation notes:
//! - OCCT's `occ::handle<Law_Function>` becomes `Arc<dyn LawFunction>`; the
//!   `Trim` result of OCCT is a new handle, so the port returns a new
//!   `Arc` as well.
//! - `GeomAbs_Shape` is [`occt_core::kernel::geomabs::Shape`]; its
//!   `PartialOrd` reproduces OCCT's enum order
//!   (`C0 < G1 < C1 < G2 < C2 < C3 < CN`), which `Law_BSpFunc` relies on
//!   (`S > Continuity()`).
//! - `NCollection_Array1<double>& T` becomes a `Vec<f64>` the callee fills;
//!   OCCT writes 1-based into `T` with `T.Lower() == 1`, so the port's
//!   vector is 0-based and the OCCT index `I` is the vector index `I - 1`.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;

/// `Law_Function` (`Law_Function.hxx:30-72`): the root class of the evolution
/// laws.
///
/// A law evaluates a scalar $\lambda(t)$ with its first and second
/// derivatives, reports the parametric bounds of its domain, and can be
/// trimmed to a sub-interval.
pub trait LawFunction: Send + Sync {
    /// `Law_Function::Continuity` — the global continuity of the law.
    fn continuity(&self) -> Shape;

    /// `Law_Function::NbIntervals(S)` — the number of intervals on which the
    /// law reaches the continuity `S`. May be one when
    /// `continuity() >= S`.
    fn nb_intervals(&self, s: Shape) -> Result<i32, String>;

    /// `Law_Function::Intervals(T, S)` — the parameters bounding the
    /// intervals of continuity `S`.
    ///
    /// OCCT expects the array to hold `nb_intervals(S) + 1` entries; the port
    /// returns exactly that many.
    fn intervals(&self, s: Shape) -> Result<Vec<f64>, String>;

    /// `Law_Function::Value(X)`.
    fn value(&self, x: f64) -> f64;

    /// `Law_Function::D1(X, F, D)` — the value `F` and the first derivative
    /// `D` at `X`.
    fn d1(&self, x: f64) -> (f64, f64);

    /// `Law_Function::D2(X, F, D, D2)` — the value `F`, the first derivative
    /// `D` and the second derivative `D2` at `X`.
    fn d2(&self, x: f64) -> (f64, f64, f64);

    /// `Law_Function::Trim(PFirst, PLast, Tol)` — a law equivalent to `self`
    /// between `PFirst` and `PLast`.
    fn trim(&self, p_first: f64, p_last: f64, tol: f64) -> Result<Arc<dyn LawFunction>, String>;

    /// `Law_Function::Bounds(PFirst, PLast)` — the parametric bounds.
    fn bounds(&self) -> (f64, f64);

    /// `occ::down_cast<Law_Linear>` support: the port has no RTTI, so
    /// `Law_Linear` identifies itself here and every other law takes the
    /// default. Used by `GeomFill_ConstrainedFilling::mklin`
    /// (`GeomFill_ConstrainedFilling.cxx:91-101`).
    fn as_linear(&self) -> Option<&LawLinear> {
        None
    }
}

/// `Law_Constant` (`Law_Constant.hxx`/`Law_Constant.cxx:24-110`).
///
/// A constant evolution law: the value is `radius` over the whole domain.
#[derive(Debug, Clone, Copy, Default)]
pub struct LawConstant {
    radius: f64,
    first: f64,
    last: f64,
}

impl LawConstant {
    /// `Law_Constant::Law_Constant` (`cxx:27-32`): an empty law.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Law_Constant::Set(Radius, PFirst, PLast)` (`cxx:36-41`).
    pub fn set(&mut self, radius: f64, p_first: f64, p_last: f64) {
        self.radius = radius;
        self.first = p_first;
        self.last = p_last;
    }

    /// The constant value.
    pub fn radius(&self) -> f64 {
        self.radius
    }
}

impl LawFunction for LawConstant {
    /// `Law_Constant::Continuity` (`cxx:46-49`): `GeomAbs_CN`.
    fn continuity(&self) -> Shape {
        Shape::CN
    }

    /// `Law_Constant::NbIntervals` (`cxx:54-57`): 1.
    fn nb_intervals(&self, _s: Shape) -> Result<i32, String> {
        Ok(1)
    }

    /// `Law_Constant::Intervals` (`cxx:62-68`): the two bounds.
    fn intervals(&self, _s: Shape) -> Result<Vec<f64>, String> {
        Ok(vec![self.first, self.last])
    }

    /// `Law_Constant::Value` (`cxx:73-76`).
    fn value(&self, _x: f64) -> f64 {
        self.radius
    }

    /// `Law_Constant::D1` (`cxx:81-85`).
    fn d1(&self, _x: f64) -> (f64, f64) {
        (self.radius, 0.0)
    }

    /// `Law_Constant::D2` (`cxx:90-94`).
    fn d2(&self, _x: f64) -> (f64, f64, f64) {
        (self.radius, 0.0, 0.0)
    }

    /// `Law_Constant::Trim` (`cxx:99-104`): the same radius over the new
    /// bounds.
    fn trim(&self, p_first: f64, p_last: f64, _tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        let mut l = LawConstant::new();
        l.set(self.radius, p_first, p_last);
        Ok(Arc::new(l))
    }

    /// `Law_Constant::Bounds` (`cxx:109-113`).
    fn bounds(&self) -> (f64, f64) {
        (self.first, self.last)
    }
}

/// `Law_Linear` (`Law_Linear.hxx`/`Law_Linear.cxx:26-103`).
///
/// A linear evolution law between the values `valdeb` at `pdeb` and `valfin`
/// at `pfin`.
#[derive(Debug, Clone, Copy, Default)]
pub struct LawLinear {
    valdeb: f64,
    valfin: f64,
    pdeb: f64,
    pfin: f64,
}

impl LawLinear {
    /// `Law_Linear::Law_Linear` (`cxx:28-33`): an empty law.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Law_Linear::Set(Pdeb, Valdeb, Pfin, Valfin)` (`cxx:35-40`).
    pub fn set(&mut self, pdeb: f64, valdeb: f64, pfin: f64, valfin: f64) {
        self.pdeb = pdeb;
        self.pfin = pfin;
        self.valdeb = valdeb;
        self.valfin = valfin;
    }

    /// The value at the first bound.
    pub fn value_at_start(&self) -> f64 {
        self.valdeb
    }

    /// The value at the last bound.
    pub fn value_at_end(&self) -> f64 {
        self.valfin
    }
}

impl LawFunction for LawLinear {
    /// `Law_Linear::Continuity` (`cxx:45-48`): `GeomAbs_CN`.
    fn continuity(&self) -> Shape {
        Shape::CN
    }

    /// `Law_Linear::NbIntervals` (`cxx:53-56`): 1.
    fn nb_intervals(&self, _s: Shape) -> Result<i32, String> {
        Ok(1)
    }

    /// `Law_Linear::Intervals` (`cxx:61-66`): the two bounds.
    fn intervals(&self, _s: Shape) -> Result<Vec<f64>, String> {
        Ok(vec![self.pdeb, self.pfin])
    }

    /// `Law_Linear::Value` (`cxx:69-72`).
    fn value(&self, x: f64) -> f64 {
        ((x - self.pdeb) * self.valfin + (self.pfin - x) * self.valdeb) / (self.pfin - self.pdeb)
    }

    /// `Law_Linear::D1` (`cxx:74-78`).
    fn d1(&self, x: f64) -> (f64, f64) {
        let f = ((x - self.pdeb) * self.valfin + (self.pfin - x) * self.valdeb)
            / (self.pfin - self.pdeb);
        let d = (self.valfin - self.valdeb) / (self.pfin - self.pdeb);
        (f, d)
    }

    /// `Law_Linear::D2` (`cxx:80-84`).
    fn d2(&self, x: f64) -> (f64, f64, f64) {
        let f = ((x - self.pdeb) * self.valfin + (self.pfin - x) * self.valdeb)
            / (self.pfin - self.pdeb);
        let d = (self.valfin - self.valdeb) / (self.pfin - self.pdeb);
        (f, d, 0.0)
    }

    /// `Law_Linear::Trim` (`cxx:89-98`): the same line restricted to the new
    /// bounds.
    fn trim(&self, p_first: f64, p_last: f64, _tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        let mut l = LawLinear::new();
        let vdeb = ((p_first - self.pdeb) * self.valfin + (self.pfin - p_first) * self.valdeb)
            / (self.pfin - self.pdeb);
        let vfin =
            ((p_last - self.pdeb) * self.valfin + (self.pfin - p_last) * self.valdeb)
                / (self.pfin - self.pdeb);
        l.set(p_first, vdeb, p_last, vfin);
        Ok(Arc::new(l))
    }

    /// `Law_Linear::Bounds` (`cxx:102-106`).
    fn bounds(&self) -> (f64, f64) {
        (self.pdeb, self.pfin)
    }

    /// `occ::down_cast<Law_Linear>` (see [`LawFunction::as_linear`]).
    fn as_linear(&self) -> Option<&LawLinear> {
        Some(self)
    }
}

/// `Law_S` (`Law_S.hxx`/`Law_S.cxx:24-52`).
///
/// A cubic law with prescribed end values and end derivatives. OCCT builds it
/// as a `Law_BSpline` of degree 3 with two knots of multiplicity 4
/// (`cxx:36-50`), and `Law_S` derives from `Law_BSpFunc`, so the port holds
/// a [`LawBSpFunc`] and delegates to it.
#[derive(Debug, Clone, Default)]
pub struct LawS {
    func: crate::law::bspline::LawBSpFunc,
}

impl LawS {
    /// `Law_S::Law_S` (`cxx:28`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Law_S::Set(Pdeb, Valdeb, Pfin, Valfin)` (`cxx:30-32`): null
    /// derivatives at both ends.
    pub fn set(&mut self, pdeb: f64, valdeb: f64, pfin: f64, valfin: f64) -> Result<(), String> {
        self.set_with_derivatives(pdeb, valdeb, 0.0, pfin, valfin, 0.0)
    }

    /// `Law_S::Set(Pdeb, Valdeb, Ddeb, Pfin, Valfin, Dfin)` (`cxx:34-52`).
    pub fn set_with_derivatives(
        &mut self,
        pdeb: f64,
        valdeb: f64,
        ddeb: f64,
        pfin: f64,
        valfin: f64,
        dfin: f64,
    ) -> Result<(), String> {
        let mut poles = vec![0.0f64; 4];
        let mut knots = vec![0.0f64; 2];
        let mut mults = vec![0i32; 2];
        poles[0] = valdeb;
        poles[3] = valfin;
        let coe = (pfin - pdeb) / 3.0;
        poles[1] = valdeb + coe * ddeb;
        poles[2] = valfin - coe * dfin;
        knots[0] = pdeb;
        knots[1] = pfin;
        mults[0] = 4;
        mults[1] = 4;
        let curve = crate::law::bspline::LawBSpline::new(poles, knots, mults, 3, false)?;
        self.func.set_curve(curve);
        Ok(())
    }

    /// The underlying `Law_BSpFunc`.
    pub fn func(&self) -> &crate::law::bspline::LawBSpFunc {
        &self.func
    }
}

impl LawFunction for LawS {
    fn continuity(&self) -> Shape {
        self.func.continuity()
    }

    fn nb_intervals(&self, s: Shape) -> Result<i32, String> {
        self.func.nb_intervals(s)
    }

    fn intervals(&self, s: Shape) -> Result<Vec<f64>, String> {
        self.func.intervals(s)
    }

    fn value(&self, x: f64) -> f64 {
        self.func.value(x)
    }

    fn d1(&self, x: f64) -> (f64, f64) {
        self.func.d1(x)
    }

    fn d2(&self, x: f64) -> (f64, f64, f64) {
        self.func.d2(x)
    }

    fn trim(&self, p_first: f64, p_last: f64, tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        self.func.trim(p_first, p_last, tol)
    }

    fn bounds(&self) -> (f64, f64) {
        self.func.bounds()
    }
}
