//! U1 domain of `cos(U2-FI2)=B*cos(U1-FI1)+C`, interval inscription,
//! critical U1, and U2(U1) monotonicity.
//! Source: `IntPatch_ImpImpIntersection.cxx` WorkWithBoundaries /
//! CriticalPointsComputing / CylCylMonotonicity / InscribeInterval.

use occt_core::precision::Precision;

use super::{inscribe_point, Coeffs, NUL_VALUE, PERIOD};

const REAL_SMALL: f64 = 1.0e-150;

/// `Bnd_Range` (void when `hi < lo`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct URange {
    lo: f64,
    hi: f64,
}

impl URange {
    pub(crate) fn void() -> Self {
        Self { lo: 0.0, hi: -1.0 }
    }

    pub(crate) fn is_void(&self) -> bool {
        self.hi < self.lo
    }

    pub(crate) fn add(&mut self, x: f64) {
        if self.is_void() {
            self.lo = x;
            self.hi = x;
        } else {
            self.lo = self.lo.min(x);
            self.hi = self.hi.max(x);
        }
    }

    pub(crate) fn delta(&self) -> f64 {
        self.hi - self.lo
    }

    pub(crate) fn get_bounds(&self) -> Option<(f64, f64)> {
        if self.is_void() {
            None
        } else {
            Some((self.lo, self.hi))
        }
    }

    pub(crate) fn get_min(&self) -> Option<f64> {
        if self.is_void() {
            None
        } else {
            Some(self.lo)
        }
    }

    pub(crate) fn get_max(&self) -> Option<f64> {
        if self.is_void() {
            None
        } else {
            Some(self.hi)
        }
    }

    pub(crate) fn set_void(&mut self) {
        self.lo = 0.0;
        self.hi = -1.0;
    }

    pub(crate) fn shift(&mut self, d: f64) {
        if !self.is_void() {
            self.lo += d;
            self.hi += d;
        }
    }

    pub(crate) fn common(&mut self, other: &Self) {
        if other.is_void() {
            self.set_void();
            return;
        }
        if self.is_void() {
            return;
        }
        self.lo = self.lo.max(other.lo);
        self.hi = self.hi.min(other.hi);
    }

    pub(crate) fn union_with(&mut self, other: &Self) -> bool {
        if self.is_void() || other.is_void() {
            return false;
        }
        if self.hi < other.lo || self.lo > other.hi {
            return false;
        }
        self.lo = self.lo.min(other.lo);
        self.hi = self.hi.max(other.hi);
        true
    }

    fn is_intersected_in(&self, val: f64, period: f64) -> bool {
        if self.is_void() {
            return false;
        }
        let period = period.abs();
        let df = self.lo - val;
        let dl = self.hi - val;
        if period <= REAL_SMALL {
            let delta = df * dl;
            if delta == 0.0 {
                return false;
            }
            return delta < 0.0;
        }
        let v1 = df / period;
        let v2 = dl / period;
        let p1 = v1.floor() as i32;
        let p2 = v2.floor() as i32;
        if p1 != p2 {
            return v2 != p2 as f64;
        }
        false
    }

    fn split(&self, val: f64, period: f64) -> Vec<URange> {
        let period = period.abs();
        if !self.is_intersected_in(val, period) {
            return vec![*self];
        }
        if period <= 0.0 {
            return vec![
                URange {
                    lo: self.lo,
                    hi: val,
                },
                URange {
                    lo: val,
                    hi: self.hi,
                },
            ];
        }
        let mut out = Vec::new();
        let mut prev = val + period * ((self.lo - val) / period).ceil();
        if prev > self.lo {
            out.push(URange {
                lo: self.lo,
                hi: prev,
            });
        }
        let mut cur = prev + period;
        while cur <= self.hi {
            out.push(URange { lo: prev, hi: cur });
            prev = cur;
            cur += period;
        }
        if prev < self.hi {
            out.push(URange {
                lo: prev,
                hi: self.hi,
            });
        }
        out
    }
}

fn clamp_acos_arg(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

/// `WorkWithBoundaries::BoundariesComputing`.
pub(crate) fn boundaries_computing(c: &Coeffs, period: f64) -> Option<[URange; 2]> {
    let mut r = [URange::void(), URange::void()];
    if c.b > 0.0 {
        if c.b + c.c.abs() < -1.0 {
            return None;
        } else if c.b + c.c.abs() <= 1.0 {
            r[0].add(c.fi1);
            r[0].add(period + c.fi1);
        } else if (1.0 + c.c <= c.b) && (c.b <= 1.0 - c.c) {
            let ang = clamp_acos_arg(-(c.c + 1.0) / c.b).acos();
            r[0].add(c.fi1);
            r[0].add(ang + c.fi1);
            r[1].add(period - ang + c.fi1);
            r[1].add(period + c.fi1);
        } else if (1.0 - c.c <= c.b) && (c.b <= 1.0 + c.c) {
            let ang = clamp_acos_arg((1.0 - c.c) / c.b).acos();
            r[0].add(ang + c.fi1);
            r[0].add(period - ang + c.fi1);
        } else if c.b - c.c.abs() >= 1.0 {
            let a1 = clamp_acos_arg((1.0 - c.c) / c.b).acos();
            let a2 = clamp_acos_arg(-(c.c + 1.0) / c.b).acos();
            r[0].add(a1 + c.fi1);
            r[0].add(a2 + c.fi1);
            r[1].add(period - a2 + c.fi1);
            r[1].add(period - a1 + c.fi1);
        } else {
            return None;
        }
    } else if c.b < 0.0 {
        if c.b + c.c.abs() > 1.0 {
            return None;
        } else if -c.b + c.c.abs() <= 1.0 {
            r[0].add(c.fi1);
            r[0].add(period + c.fi1);
        } else if (-c.c - 1.0 <= c.b) && (c.b <= c.c - 1.0) {
            let ang = clamp_acos_arg((1.0 - c.c) / c.b).acos();
            r[0].add(c.fi1);
            r[0].add(ang + c.fi1);
            r[1].add(period - ang + c.fi1);
            r[1].add(period + c.fi1);
        } else if (c.c - 1.0 <= c.b) && (c.b <= -c.b - 1.0) {
            let ang = clamp_acos_arg(-(c.c + 1.0) / c.b).acos();
            r[0].add(ang + c.fi1);
            r[0].add(period - ang + c.fi1);
        } else if -c.b - c.c.abs() >= 1.0 {
            let a1 = clamp_acos_arg(-(c.c + 1.0) / c.b).acos();
            let a2 = clamp_acos_arg((1.0 - c.c) / c.b).acos();
            r[0].add(a1 + c.fi1);
            r[0].add(a2 + c.fi1);
            r[1].add(period - a2 + c.fi1);
            r[1].add(period - a1 + c.fi1);
        } else {
            return None;
        }
    } else {
        return None;
    }
    Some(r)
}

/// `InscribeInterval`.
pub(crate) fn inscribe_interval(
    uf: f64,
    ul: f64,
    range: &mut URange,
    tol2d: f64,
    period: f64,
) -> bool {
    let Some(mut u) = range.get_min() else {
        return false;
    };
    let delta = range.delta();
    let force_min = (ul - u).abs() < tol2d;
    if inscribe_point(uf, ul, &mut u, tol2d, period, force_min) {
        range.set_void();
        range.add(u);
        range.add(u + delta);
        return true;
    }
    let Some(mut u) = range.get_max() else {
        return false;
    };
    let force_max = (uf - u).abs() < tol2d;
    if inscribe_point(uf, ul, &mut u, tol2d, period, force_max) {
        range.set_void();
        range.add(u);
        range.add(u - delta);
        return true;
    }
    false
}

fn exclude_near(arr: &mut [f64], u1f: f64, u1l: f64, tol: f64) -> bool {
    let mut changed = false;
    for i in 1..arr.len() {
        if Precision::is_infinite(arr[i]) {
            break;
        }
        if arr[i] - arr[i - 1] < tol {
            if arr[i - 1] != 0.0 && arr[i - 1] != u1f && arr[i - 1] != u1l {
                arr[i] = 0.5 * (arr[i] + arr[i - 1]);
            } else {
                arr[i] = arr[i - 1];
            }
            arr[i - 1] = Precision::INFINITE;
            changed = true;
        }
    }
    changed
}

fn crit_from_cos(c: &Coeffs, target: f64, sign: f64) -> f64 {
    if c.b.abs() <= NUL_VALUE {
        return Precision::INFINITE;
    }
    let a = (target - c.c) / c.b;
    if a.abs() < 1.0 {
        sign * a.acos() + c.fi1
    } else {
        Precision::INFINITE
    }
}

/// `CriticalPointsComputing`. Returns sorted finite U1 values in `[0, period)`.
pub(crate) fn critical_points(
    c: &Coeffs,
    u1f: f64,
    u1l: f64,
    u2f: f64,
    u2l: f64,
    period: f64,
    tol2d: f64,
) -> Vec<f64> {
    let mut u = [Precision::INFINITE; 12];
    u[0] = 0.0;
    u[1] = period;
    u[2] = u1f;
    u[3] = u1l;
    let cos_fi2 = c.fi2.cos();
    let bsb = c.b.abs();
    if (c.c - bsb <= cos_fi2) && (cos_fi2 <= c.c + bsb) && c.b.abs() > NUL_VALUE {
        let arg = ((cos_fi2 - c.c) / c.b).clamp(-1.0, 1.0);
        u[4] = -arg.acos() + c.fi1;
        u[5] = arg.acos() + c.fi1;
    }
    let mut sf = (u2f - c.fi2).cos();
    let mut sl = (u2l - c.fi2).cos();
    if sf > sl {
        std::mem::swap(&mut sf, &mut sl);
    }
    u[6] = crit_from_cos(c, sl, -1.0);
    u[7] = crit_from_cos(c, sf, -1.0);
    u[8] = crit_from_cos(c, sf, 1.0);
    u[9] = crit_from_cos(c, sl, 1.0);
    u[10] = c.fi1;
    u[11] = std::f64::consts::PI + c.fi1;
    for x in &mut u {
        if Precision::is_infinite(*x) {
            continue;
        }
        *x %= period;
        if *x < 0.0 {
            *x += period;
        }
    }
    let mut n = 12;
    loop {
        u[..n].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if !exclude_near(&mut u[..n], u1f, u1l, tol2d) {
            break;
        }
    }
    while n > 0 && Precision::is_infinite(u[n - 1]) {
        n -= 1;
    }
    if n > 1 {
        let a = u[0];
        let b = u[n - 1];
        if (b - a - period).abs() < tol2d {
            u[0] = (a + b - period) / 2.0;
            n -= 1;
        }
    }
    u[..n].to_vec()
}

/// `ComputationMethods::CylCylMonotonicity`.
pub(crate) fn cyl_cyl_monotonicity(u1: f64, wl: i32, c: &Coeffs, period: f64) -> Option<bool> {
    let is_plus = match wl {
        0 => true,
        1 => false,
        _ => return None,
    };
    let mut u_tmp = u1 - c.fi1;
    if !inscribe_point(0.0, period, &mut u_tmp, 0.0, period, false) {
        return None;
    }
    let mut increasing = true;
    if (std::f64::consts::PI - u_tmp) < REAL_SMALL && u_tmp < period {
        increasing = false;
    }
    if c.b < 0.0 {
        increasing = !increasing;
    }
    if !is_plus {
        increasing = !increasing;
    }
    Some(increasing)
}

/// Common inscribed length of analytic U1 ranges vs cylinder U bounds (`IntCyCy` aSumRange).
pub(crate) fn sum_inscribed_u(ranges: [URange; 2], uf: f64, ul: f64, tol2d: f64) -> f64 {
    let mut list = vec![ranges[0], ranges[1]];
    let splits = [uf, ul, 0.0];
    for s in splits {
        let tmp = list;
        list = Vec::new();
        for r in tmp {
            list.extend(r.split(s, PERIOD));
        }
    }
    let mut sum = 0.0;
    for mut cur in list {
        let mut bound = URange::void();
        bound.add(uf);
        bound.add(ul);
        if !inscribe_interval(uf, ul, &mut cur, tol2d, PERIOD) {
            let Some((f, l)) = cur.get_bounds() else {
                continue;
            };
            if l < uf {
                cur.shift(PERIOD);
            } else if f > ul {
                cur.shift(-PERIOD);
            }
        }
        bound.common(&cur);
        let d = bound.delta();
        if d > 0.0 {
            sum += d;
        }
    }
    sum
}

pub(crate) fn merge_ranges(ranges: &mut [URange; 2]) {
    let other = ranges[1];
    if ranges[0].union_with(&other) {
        ranges[1].set_void();
    }
}
