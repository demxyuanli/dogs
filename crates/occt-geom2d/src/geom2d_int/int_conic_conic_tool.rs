//! Port of `IntCurve_IntConicConic_Tool`
//! (`IntCurve_IntConicConic_Tool.hxx:26-149`,
//! `IntCurve_IntConicConic_Tool.cxx:20-294`): the `Interval` /
//! `PeriodicInterval` interval algebra, `Determine_Transition_LC` and
//! `NormalizeOnCircleDomain`, all of which the `IntCurve_IntConicConic::Perform`
//! overloads of `IntCurve_IntConicConic_1.cxx` build on.

use std::f64::consts::PI;

use occt_core::gp::GpVec2d;
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dPosition, IntRes2dSituation, IntRes2dTransition, IntRes2dTypeTrans,
};
use occt_core::precision::COMPUTATIONAL;

/// `#define TOLERANCE_ANGULAIRE 0.00000001` (`Tool.cxx:20`). This is the
/// translation-unit-local value of `IntCurve_IntConicConic_Tool.cxx`;
/// `IntCurve_IntConicConic_1.cxx:41` redefines the same macro to `1.e-15`
/// inside its own unit, so the two constants must not be merged.
const TOLERANCE_ANGULAIRE: f64 = 0.00000001;

/// `static double PIpPI = M_PI + M_PI` (`Tool.hxx:26`).
pub(crate) fn pi_p_pi() -> f64 {
    PI + PI
}

/// `Determine_Transition_LC` (`Tool.cxx:26-93`). `Tan1` is normalized in the
/// TOUCH branch, as in OCCT, so it is taken by mutable reference.
#[allow(clippy::too_many_arguments)]
pub(crate) fn determine_transition_lc(
    pos1: IntRes2dPosition,
    tan1: &mut GpVec2d,
    norm1: &GpVec2d,
    t1: &mut IntRes2dTransition,
    pos2: IntRes2dPosition,
    tan2: &mut GpVec2d,
    norm2: &GpVec2d,
    t2: &mut IntRes2dTransition,
    _tol: f64,
) {
    let sgn = tan1.crossed(tan2);
    let norm = tan1.magnitude() * tan2.magnitude();

    if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
        // Transition TOUCH (######### in OCCT)
        let opos = tan1.dot(tan2) < 0.0;

        // Modified by Sergey KHROMOV - Thu Nov 2 17:57:15 2000
        let _ = tan1.normalize();
        let norm = GpVec2d::new(-tan1.y(), tan1.x());

        let val1 = norm.dot(norm1);
        let val2 = norm.dot(norm2);

        if (val1 - val2).abs() <= COMPUTATIONAL {
            t1.set_touch(true, pos1, IntRes2dSituation::Unknown, opos);
            t2.set_touch(true, pos2, IntRes2dSituation::Unknown, opos);
        } else if val2 > val1 {
            t2.set_touch(true, pos2, IntRes2dSituation::Inside, opos);
            if opos {
                t1.set_touch(true, pos1, IntRes2dSituation::Inside, opos);
            } else {
                t1.set_touch(true, pos1, IntRes2dSituation::Outside, opos);
            }
        } else {
            // Val1 > Val2
            t2.set_touch(true, pos2, IntRes2dSituation::Outside, opos);
            if opos {
                t1.set_touch(true, pos1, IntRes2dSituation::Outside, opos);
            } else {
                t1.set_touch(true, pos1, IntRes2dSituation::Inside, opos);
            }
        }
    } else if sgn < 0.0 {
        t1.set_in_out(false, pos1, IntRes2dTypeTrans::In);
        t2.set_in_out(false, pos2, IntRes2dTypeTrans::Out);
    } else {
        t1.set_in_out(false, pos1, IntRes2dTypeTrans::Out);
        t2.set_in_out(false, pos2, IntRes2dTypeTrans::In);
    }
}

/// `NormalizeOnCircleDomain` (`Tool.cxx:96-108`).
pub(crate) fn normalize_on_circle_domain(param: f64, the_domain: &IntRes2dDomain) -> f64 {
    let mut param = param;
    while param < the_domain.first_parameter() {
        param += pi_p_pi();
    }
    while param > the_domain.last_parameter() {
        param -= pi_p_pi();
    }
    param
}

/// `class Interval` (`Tool.hxx:47-61`, `Tool.cxx:190-294`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Interval {
    /// `Binf` (`Tool.hxx:50`)
    pub binf: f64,
    /// `Bsup` (`Tool.hxx:51`)
    pub bsup: f64,
    /// `HasFirstBound` (`Tool.hxx:52`)
    pub has_first_bound: bool,
    /// `HasLastBound` (`Tool.hxx:53`)
    pub has_last_bound: bool,
    /// `IsNull` (`Tool.hxx:54`)
    pub is_null: bool,
}

impl Default for Interval {
    /// `Interval()` (`Tool.cxx:190-198`).
    fn default() -> Self {
        Self { binf: 0.0, bsup: 0.0, has_first_bound: false, has_last_bound: false, is_null: true }
    }
}

impl Interval {
    /// `Interval(a, b)` (`Tool.cxx:199-214`).
    pub(crate) fn from_bounds(a: f64, b: f64) -> Self {
        let (binf, bsup) = if a < b { (a, b) } else { (b, a) };
        Self { binf, bsup, has_first_bound: true, has_last_bound: true, is_null: false }
    }

    /// `Interval(const IntRes2d_Domain&)` (`Tool.cxx:215-239`).
    pub(crate) fn from_domain(domain: &IntRes2dDomain) -> Self {
        let mut r = Self {
            binf: 0.0,
            bsup: 0.0,
            has_first_bound: false,
            has_last_bound: false,
            is_null: false,
        };
        if domain.has_first_point() {
            r.has_first_bound = true;
            r.binf = domain.first_parameter() - domain.first_tolerance();
        }
        if domain.has_last_point() {
            r.has_last_bound = true;
            r.bsup = domain.last_parameter() + domain.last_tolerance();
        }
        r
    }

    /// `Interval(a, hf, b, hl)` (`Tool.cxx:240-247`).
    pub(crate) fn from_bounds_flags(a: f64, hf: bool, b: f64, hl: bool) -> Self {
        Self { binf: a, bsup: b, has_first_bound: hf, has_last_bound: hl, is_null: false }
    }

    /// `Length()` (`Tool.cxx:249-252`).
    pub(crate) fn length(&self) -> f64 {
        if self.is_null {
            -1.0
        } else {
            (self.bsup - self.binf).abs()
        }
    }

    /// `IntersectionWithBounded` (`Tool.cxx:254-294`).
    pub(crate) fn intersection_with_bounded(&self, inter: &Interval) -> Interval {
        if self.is_null || inter.is_null {
            return Interval::default();
        }
        if !(self.has_first_bound || self.has_last_bound) {
            return Interval::from_bounds(inter.binf, inter.bsup);
        }

        let a;
        if self.has_first_bound {
            if inter.bsup < self.binf {
                return Interval::default();
            }
            a = if inter.binf < self.binf { self.binf } else { inter.binf };
        } else {
            a = inter.binf;
        }

        let b;
        if self.has_last_bound {
            if inter.binf > self.bsup {
                return Interval::default();
            }
            b = if inter.bsup > self.bsup { self.bsup } else { inter.bsup };
        } else {
            b = inter.bsup;
        }

        Interval::from_bounds(a, b)
    }
}

/// `class PeriodicInterval` (`Tool.hxx:67-148`, `Tool.cxx:111-187`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct PeriodicInterval {
    /// `Binf` (`Tool.hxx:70`)
    pub binf: f64,
    /// `Bsup` (`Tool.hxx:71`)
    pub bsup: f64,
    /// `isnull` (`Tool.hxx:72`)
    pub is_null: bool,
}

impl Default for PeriodicInterval {
    /// `PeriodicInterval()` (`Tool.hxx:108-112`).
    fn default() -> Self {
        Self { binf: 0.0, bsup: 0.0, is_null: true }
    }
}

impl PeriodicInterval {
    /// `PeriodicInterval(const IntRes2d_Domain&)` (`Tool.hxx:94-106`).
    pub(crate) fn from_domain(domain: &IntRes2dDomain) -> Self {
        let binf = if domain.has_first_point() { domain.first_parameter() } else { -1.0 };
        let bsup = if domain.has_last_point() { domain.last_parameter() } else { 20.0 };
        Self { binf, bsup, is_null: false }
    }

    /// `PeriodicInterval(a, b)` (`Tool.hxx:114-122`).
    pub(crate) fn from_bounds(a: f64, b: f64) -> Self {
        let mut r = Self { binf: a, bsup: b, is_null: false };
        if (b - a) < pi_p_pi() {
            r.normalize();
        }
        r
    }

    /// `SetValues(a, b)` (`Tool.hxx:124-130`).
    pub(crate) fn set_values(&mut self, a: f64, b: f64) {
        self.is_null = false;
        self.binf = a;
        self.bsup = b;
        if (b - a) < pi_p_pi() {
            self.normalize();
        }
    }

    /// `SetNull()` (`Tool.hxx:74`).
    pub(crate) fn set_null(&mut self) {
        self.is_null = true;
    }

    /// `IsNull()` (`Tool.hxx:76`).
    pub(crate) fn is_null(&self) -> bool {
        self.is_null
    }

    /// `Complement()` (`Tool.hxx:78-91`).
    pub(crate) fn complement(&mut self) {
        if !self.is_null {
            let t = self.binf;
            self.binf = self.bsup;
            self.bsup = t + pi_p_pi();
            if self.binf > pi_p_pi() {
                self.binf -= pi_p_pi();
                self.bsup -= pi_p_pi();
            }
        }
    }

    /// `Length()` (`Tool.hxx:93`).
    pub(crate) fn length(&self) -> f64 {
        if self.is_null {
            -100.0
        } else {
            (self.bsup - self.binf).abs()
        }
    }

    /// `Normalize()` (`Tool.hxx:132-145`).
    pub(crate) fn normalize(&mut self) {
        if !self.is_null {
            while self.binf > pi_p_pi() {
                self.binf -= pi_p_pi();
            }
            while self.binf < 0.0 {
                self.binf += pi_p_pi();
            }
            while self.bsup < self.binf {
                self.bsup += pi_p_pi();
            }
            while self.bsup >= (self.binf + pi_p_pi()) {
                self.bsup -= pi_p_pi();
            }
        }
    }

    /// `FirstIntersection` (`Tool.cxx:111-155`). OCCT takes `PInter` by
    /// non-const reference and shifts it in place, so it is mutated here too.
    pub(crate) fn first_intersection(&self, p_inter: &mut PeriodicInterval) -> PeriodicInterval {
        if p_inter.is_null || self.is_null {
            return PeriodicInterval::default();
        }
        if self.length() >= pi_p_pi() {
            return PeriodicInterval::from_bounds(p_inter.binf, p_inter.bsup);
        }
        if p_inter.length() >= pi_p_pi() {
            return PeriodicInterval::from_bounds(self.binf, self.bsup);
        }
        if p_inter.bsup <= self.binf {
            while p_inter.binf <= self.binf && p_inter.bsup <= self.binf {
                p_inter.binf += pi_p_pi();
                p_inter.bsup += pi_p_pi();
            }
        }
        if p_inter.binf >= self.bsup {
            while p_inter.binf >= self.bsup && p_inter.bsup >= self.bsup {
                p_inter.binf -= pi_p_pi();
                p_inter.bsup -= pi_p_pi();
            }
        }
        if p_inter.bsup < self.binf || p_inter.binf > self.bsup {
            return PeriodicInterval::default();
        }

        let a = if p_inter.binf > self.binf { p_inter.binf } else { self.binf };
        let b = if p_inter.bsup < self.bsup { p_inter.bsup } else { self.bsup };

        PeriodicInterval::from_bounds(a, b)
    }

    /// `SecondIntersection` (`Tool.cxx:159-187`).
    pub(crate) fn second_intersection(&self, p_inter: &mut PeriodicInterval) -> PeriodicInterval {
        if p_inter.is_null || self.is_null || self.length() >= pi_p_pi() || p_inter.length() >= pi_p_pi()
        {
            return PeriodicInterval::default();
        }

        let mut p_inter_inf = p_inter.binf + pi_p_pi();
        let mut p_inter_sup = p_inter.bsup + pi_p_pi();
        if p_inter_inf > self.bsup {
            p_inter_inf = p_inter.binf - pi_p_pi();
            p_inter_sup = p_inter.bsup - pi_p_pi();
        }
        if p_inter_sup < self.binf || p_inter_inf > self.bsup {
            return PeriodicInterval::default();
        }

        let a = if p_inter_inf > self.binf { p_inter_inf } else { self.binf };
        let b = if p_inter_sup < self.bsup { p_inter_sup } else { self.bsup };
        PeriodicInterval::from_bounds(a, b)
    }
}
