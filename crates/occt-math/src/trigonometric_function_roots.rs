//! `math_TrigonometricFunctionRoots` — roots of
//! `A*cos(x)^2 + 2*B*cos(x)*sin(x) + C*cos(x) + D*sin(x) + E = 0` on an interval.
//!
//! Source: `FoundationClasses/TKMath/math/math_TrigonometricFunctionRoots.cxx:32-572`
//! and `math_TrigonometricEquationFunction.hxx` (OCCT 8.0.0):
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `math_TrigonometricEquationFunction` | `hxx:25-92` | [`TrigonometricEquationFunction`] |
//! | `math_TrigonometricFunctionRoots(D, E, lo, hi)` | `cxx:32-47` | [`TrigonometricFunctionRoots::linear`] |
//! | `math_TrigonometricFunctionRoots(C, D, E, lo, hi)` | `cxx:49-59` | [`TrigonometricFunctionRoots::quadratic`] |
//! | `math_TrigonometricFunctionRoots(A..E, lo, hi)` | `cxx:60-75` | [`TrigonometricFunctionRoots::new`] |
//! | `Perform` | `cxx:77-550` | [`TrigonometricFunctionRoots::perform`] |
//!
//! Translation notes:
//! - `NCollection_Array1<double> ko(1, 5)` / `Zer(1, 4)` / `Sol(1, 4)` are kept
//!   as 1-based arrays with a wasted slot 0 and a little slack
//!   (`[f64; 9]`), so every index of the C++ body transcribes unchanged,
//! - `RealFirst()` / `RealLast()` are `-DBL_MAX` / `DBL_MAX`,
//! - `math_DirectPolynomialRoots` is
//!   [`crate::direct_polynomial_roots::DirectPolynomialRoots`]; a `!IsDone()`
//!   answer is the `Done = false` return of `cxx:218-221`, `:369-372`,
//! - `math_NewtonFunctionRoot` is
//!   [`crate::newton_function_root::NewtonFunctionRoot`]; its `IsDone()`
//!   wording is the `Resol.IsDone()` of `cxx:465`,
//! - `Epsilon(x)` is [`occt_core::precision::epsilon`],
//! - `Dump` (`cxx:552-572`) is `UNPORTED` (no `Standard_OStream` in the port).

use occt_core::precision::{epsilon, PCONFUSION};

use crate::direct_polynomial_roots::DirectPolynomialRoots;
use crate::function_set_root::FunctionWithDerivative;
use crate::newton_function_root::NewtonFunctionRoot;

/// 1-based access into the fixed arrays of the C++ body.
const N: usize = 9;

/// `math_TrigonometricEquationFunction` (`hxx:25-92`): the function and its
/// first derivative, as `math_NewtonFunctionRoot` needs them.
#[derive(Debug, Clone, Copy)]
pub struct TrigonometricEquationFunction {
    aa: f64,
    bb: f64,
    cc: f64,
    dd: f64,
    ee: f64,
}

impl TrigonometricEquationFunction {
    /// `math_TrigonometricEquationFunction(A, B, C, D, E)` (`hxx:36-44`).
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        Self { aa: a, bb: b, cc: c, dd: d, ee: e }
    }

    /// `Value(X, F)` (`hxx:46-53`).
    pub fn evaluate(&self, x: f64) -> f64 {
        let cn = x.cos();
        let sn = x.sin();
        // `F = CN * (AA*CN + (BB+BB)*SN + CC) + DD*SN + EE;`
        cn * (self.aa * cn + (self.bb + self.bb) * sn + self.cc) + self.dd * sn + self.ee
    }

    /// `Derivative(X, D)` (`hxx:55-63`).
    pub fn evaluate_derivative(&self, x: f64) -> f64 {
        let cn = x.cos();
        let sn = x.sin();
        // `D = -AA*CN*SN + BB*(CN*CN - SN*SN); D += D; D += -CC*SN + DD*CN;`
        let mut d = -self.aa * cn * sn + self.bb * (cn * cn - sn * sn);
        d += d;
        d += -self.cc * sn + self.dd * cn;
        d
    }
}

impl FunctionWithDerivative for TrigonometricEquationFunction {
    /// `Value(X, F)` (`hxx:46-53`).
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.evaluate(x);
        true
    }

    /// `Derivative(X, D)` (`hxx:55-63`).
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        *d = self.evaluate_derivative(x);
        true
    }

    /// `Values(X, F, D)` (`hxx:65-90`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let cn = x.cos();
        let sn = x.sin();
        let aa_cn = self.aa * cn;
        let bb_sn = self.bb * sn;

        *f = aa_cn * cn + bb_sn * (cn + cn) + self.cc * cn + self.dd * sn + self.ee;
        *d = -aa_cn * sn + self.bb * (cn * cn - sn * sn);
        *d += *d;
        *d += -self.cc * sn + self.dd * cn;
        true
    }
}

/// `math_TrigonometricFunctionRoots` (`math_TrigonometricFunctionRoots.hxx`).
#[derive(Debug, Clone)]
pub struct TrigonometricFunctionRoots {
    /// `int NbSol` — `-1` before `Done` (`cxx:36`).
    nb_sol: i32,
    /// `NCollection_Array1<double> Sol(1, 4)` (`cxx:37`), 1-based.
    sol: [f64; N],
    /// `bool InfiniteStatus` (`cxx:38`).
    infinite_status: bool,
    /// `bool Done` (`cxx:39`).
    done: bool,
}

impl Default for TrigonometricFunctionRoots {
    fn default() -> Self {
        Self { nb_sol: -1, sol: [0.0; N], infinite_status: false, done: false }
    }
}

impl TrigonometricFunctionRoots {
    /// `math_TrigonometricFunctionRoots(A, B, C, D, E, InfBound, SupBound)`
    /// (`cxx:60-75`).
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        let mut r = Self::default();
        r.perform(a, b, c, d, e, inf_bound, sup_bound);
        r
    }

    /// `math_TrigonometricFunctionRoots(D, E, InfBound, SupBound)`
    /// (`cxx:32-47`): `A = B = C = 0`.
    pub fn linear(d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        Self::new(0.0, 0.0, 0.0, d, e, inf_bound, sup_bound)
    }

    /// `math_TrigonometricFunctionRoots(C, D, E, InfBound, SupBound)`
    /// (`cxx:49-59`): `A = B = 0`.
    pub fn quadratic(c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        Self::new(0.0, 0.0, c, d, e, inf_bound, sup_bound)
    }

    /// `IsDone()` (`math_TrigonometricFunctionRoots.hxx`).
    pub fn is_done(&self) -> bool { self.done }

    /// `InfiniteRoots()`.
    pub fn infinite_roots(&self) -> bool { self.infinite_status }

    /// `NbSolutions()`.
    pub fn nb_solutions(&self) -> i32 { self.nb_sol }

    /// `Value(theIndex)` (`Sol(theIndex)`), 1-based as in OCCT.
    pub fn value(&self, index: i32) -> f64 {
        if index < 1 || index > self.nb_sol {
            return f64::NAN;
        }
        self.sol[index as usize]
    }

    /// The roots `1..=NbSolutions`, in OCCT order.
    pub fn roots(&self) -> &[f64] {
        let n = if self.nb_sol > 0 { self.nb_sol as usize } else { 0 };
        &self.sol[1..=n.min(N - 1)]
    }

    /// `Perform(A, B, C, D, E, InfBound, SupBound)` (`cxx:77-550`).
    pub fn perform(
        &mut self,
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        e: f64,
        inf_bound: f64,
        sup_bound: f64,
    ) {
        let mut nzer: usize;
        let nit = 10;
        // `NCollection_Array1<double> ko(1, 5), Zer(1, 4);` (`cxx:88-89`).
        let mut ko = [0.0f64; N];
        let mut zer = [0.0f64; N];
        let tol1 = 1.0e-15;

        self.infinite_status = false;
        self.done = true;

        // `cxx:93`
        let eps = 1.5e-12;

        // `cxx:95-123`: the interval and its period offset.
        let depi = std::f64::consts::PI + std::f64::consts::PI;
        let my_borne_inf;
        let mut delta;
        let modu;
        if inf_bound <= f64::MIN && sup_bound >= f64::MAX {
            my_borne_inf = 0.0;
            delta = depi;
            modu = 0.0;
        } else if sup_bound >= f64::MAX {
            my_borne_inf = inf_bound;
            delta = depi;
            modu = my_borne_inf / depi;
        } else if inf_bound <= f64::MIN {
            my_borne_inf = sup_bound - depi;
            delta = depi;
            modu = my_borne_inf / depi;
        } else {
            my_borne_inf = inf_bound;
            delta = sup_bound - inf_bound;
            modu = inf_bound / depi;
            if (sup_bound - inf_bound) > depi {
                delta = depi;
            }
        }

        // `cxx:125-237`: the cosine-sine coefficient arms.
        if a.abs() <= eps && b.abs() <= eps {
            if c.abs() <= eps {
                if d.abs() <= eps {
                    if e.abs() <= eps {
                        // `cxx:131-134`: infinite solutions.
                        self.infinite_status = true;
                        return;
                    } else {
                        // `cxx:136-140`
                        self.nb_sol = 0;
                        return;
                    }
                } else {
                    // `cxx:143-172`: d*sin(x) + e = 0.
                    self.nb_sol = 0;
                    let aa = -e / d;
                    if aa.abs() > 1.0 {
                        return;
                    }
                    zer[1] = aa.asin();
                    zer[2] = std::f64::consts::PI - zer[1];
                    nzer = 2;
                    for i in 1..=nzer {
                        if zer[i] <= -eps {
                            zer[i] = depi - zer[i].abs();
                        }
                        // `cxx:164-165`
                        zer[i] += modu.trunc() * depi;
                        let x = zer[i] - my_borne_inf;
                        if x > -epsilon(delta) && x < delta + epsilon(delta) {
                            self.nb_sol += 1;
                            self.sol[self.nb_sol as usize] = zer[i];
                        }
                    }
                    return;
                }
            } else if d.abs() <= eps {
                // `cxx:176-207`: c*cos(x) + e = 0.
                self.nb_sol = 0;
                let aa = -e / c;
                if aa.abs() > 1.0 {
                    return;
                }
                zer[1] = aa.acos();
                zer[2] = -zer[1];
                nzer = 2;
                for i in 1..=nzer {
                    if zer[i] <= -eps {
                        zer[i] = depi - zer[i].abs();
                    }
                    // `cxx:198`: `+= trunc(Mod) * 2 * PI`.
                    zer[i] += modu.trunc() * 2.0 * std::f64::consts::PI;
                    let x = zer[i] - my_borne_inf;
                    if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                        self.nb_sol += 1;
                        self.sol[self.nb_sol as usize] = zer[i];
                    }
                }
                return;
            } else {
                // `cxx:208-236`: the second-degree reduction AA*t^2 + BB*t + CC.
                let aa = e - c;
                let bb = 2.0 * d;
                let cc = e + c;
                let resol = DirectPolynomialRoots::quadratic(aa, bb, cc);
                if !resol.is_done() {
                    self.done = false;
                    return;
                } else if !resol.infinite_roots() {
                    nzer = resol.nb_solutions();
                    for i in 1..=nzer {
                        zer[i] = resol.value(i);
                    }
                } else {
                    // `cxx:231-235`
                    self.infinite_status = true;
                    return;
                }
            }
        } else {
            // `cxx:238-435`: the two additional analytical cases, then the
            // quartic.
            if a.abs() <= eps && e.abs() <= eps {
                if c.abs() <= eps {
                    // `cxx:243-296`: 2*B*sin*cos + D*sin = 0.
                    nzer = 2;
                    zer[1] = 0.0;
                    zer[2] = std::f64::consts::PI;
                    let aa = -d / (b * 2.0);
                    if aa.abs() <= 1.0 + PCONFUSION {
                        nzer = 4;
                        if aa >= 1.0 {
                            zer[3] = 0.0;
                            zer[4] = 0.0;
                        } else if aa <= -1.0 {
                            zer[3] = std::f64::consts::PI;
                            zer[4] = std::f64::consts::PI;
                        } else {
                            zer[3] = aa.acos();
                            zer[4] = depi - zer[3];
                        }
                    }

                    self.nb_sol = 0;
                    for i in 1..=nzer {
                        if zer[i] <= my_borne_inf - eps {
                            zer[i] += depi;
                        }
                        zer[i] += modu.trunc() * 2.0 * std::f64::consts::PI;
                        let x = zer[i] - my_borne_inf;
                        if x >= -PCONFUSION && x <= delta + PCONFUSION {
                            if zer[i] < inf_bound {
                                zer[i] = inf_bound;
                            }
                            if zer[i] > sup_bound {
                                zer[i] = sup_bound;
                            }
                            self.nb_sol += 1;
                            self.sol[self.nb_sol as usize] = zer[i];
                        }
                    }
                    return;
                }
                if d.abs() <= eps {
                    // `cxx:299-351`: 2*B*sin*cos + C*cos = 0.
                    nzer = 2;
                    zer[1] = std::f64::consts::PI / 2.0;
                    zer[2] = std::f64::consts::PI * 3.0 / 2.0;

                    let aa = -c / (b * 2.0);
                    if aa.abs() <= 1.0 + PCONFUSION {
                        nzer = 4;
                        if aa >= 1.0 {
                            zer[3] = std::f64::consts::PI / 2.0;
                            zer[4] = std::f64::consts::PI / 2.0;
                        } else if aa <= -1.0 {
                            zer[3] = std::f64::consts::PI * 3.0 / 2.0;
                            zer[4] = std::f64::consts::PI * 3.0 / 2.0;
                        } else {
                            zer[3] = aa.asin();
                            zer[4] = std::f64::consts::PI - zer[3];
                        }
                    }

                    self.nb_sol = 0;
                    for i in 1..=nzer {
                        if zer[i] <= my_borne_inf - eps {
                            zer[i] += depi;
                        }
                        zer[i] += modu.trunc() * 2.0 * std::f64::consts::PI;
                        let x = zer[i] - my_borne_inf;
                        if x >= -PCONFUSION && x <= delta + PCONFUSION {
                            if zer[i] < inf_bound {
                                zer[i] = inf_bound;
                            }
                            if zer[i] > sup_bound {
                                zer[i] = sup_bound;
                            }
                            self.nb_sol += 1;
                            self.sol[self.nb_sol as usize] = zer[i];
                        }
                    }
                    return;
                }
            }

            // `cxx:355-435`: the fourth-degree reduction in tan(x/2).
            ko[1] = a - c + e;
            ko[2] = 2.0 * d - 4.0 * b;
            ko[3] = 2.0 * e - 2.0 * a;
            ko[4] = 4.0 * b + 2.0 * d;
            ko[5] = a + c + e;
            let mut bko;
            loop {
                bko = false;
                let resol4 = DirectPolynomialRoots::quartic(ko[1], ko[2], ko[3], ko[4], ko[5]);
                if !resol4.is_done() {
                    self.done = false;
                    return;
                } else if !resol4.infinite_roots() {
                    nzer = resol4.nb_solutions();
                    for i in 1..=nzer {
                        zer[i] = resol4.value(i);
                    }
                } else {
                    // `cxx:381-385`
                    self.infinite_status = true;
                    return;
                }

                // `cxx:387-400`: sort ascending.
                let mut triok;
                loop {
                    triok = true;
                    for i in 1..nzer {
                        if zer[i] > zer[i + 1] {
                            zer.swap(i, i + 1);
                            triok = false;
                        }
                    }
                    if triok {
                        break;
                    }
                }

                // `cxx:402-424`: a near-double root whose residual is too large
                // means the leading coefficients are ill-conditioned; the body
                // shrinks them by 1e-4 and retries.
                for i in 1..nzer {
                    if (zer[i + 1] - zer[i]).abs() < eps {
                        let qw = zer[i + 1];
                        let va = ko[4] + qw * (2.0 * ko[3] + qw * (3.0 * ko[2] + qw * (4.0 * ko[1])));
                        if va.abs() > eps {
                            bko = true;
                            break;
                        }
                    }
                }
                if bko {
                    for k in 1..=5 {
                        ko[k] *= 0.0001;
                    }
                } else {
                    break;
                }
            }
        }

        // `cxx:437-504`: verify the roots against the bounds, polish each with
        // Newton and insert it in ascending order.
        let supm_infs100 = (sup_bound - inf_bound) * 0.01;
        self.nb_sol = 0;
        for i in 1..=nzer {
            let mut teta = zer[i].atan();
            teta += teta;
            if zer[i] <= -eps {
                teta = depi - teta.abs();
            }
            teta += modu.trunc() * depi;
            if teta - my_borne_inf < 0.0 {
                teta += depi;
            }

            let x = teta - my_borne_inf;
            if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                let x = teta;

                // `cxx:462-466`: the Newton call.
                let mut f = TrigonometricEquationFunction::new(a, b, c, d, e);
                let resol = NewtonFunctionRoot::new(&mut f, x, tol1, eps, nit);
                let mut teta_newton = teta;
                if resol.is_done() {
                    teta_newton = resol.root();
                }
                // `cxx:469-478`: the polished value is only kept when it stayed
                // near the bracketed root.
                let delta_newton = teta_newton - teta;
                if delta_newton > supm_infs100 || delta_newton < -supm_infs100 {
                    // `cxx:473`: OCCT only traces here.
                } else {
                    teta = teta_newton;
                }

                // `cxx:480-502`: insertion sort.
                let mut flag4 = false;
                let mut k = 1;
                while k <= self.nb_sol {
                    if teta < self.sol[k as usize] {
                        let mut l = k;
                        while l <= self.nb_sol {
                            let j = self.nb_sol - l + k;
                            self.sol[(j + 1) as usize] = self.sol[j as usize];
                            l += 1;
                        }
                        self.sol[k as usize] = teta;
                        self.nb_sol += 1;
                        flag4 = true;
                        break;
                    }
                    k += 1;
                }
                if !flag4 {
                    self.nb_sol += 1;
                    self.sol[self.nb_sol as usize] = teta;
                }
            }
        }

        // `cxx:505-550`: the special case of PI.
        if self.nb_sol < 4 {
            let start_index = self.nb_sol + 1;
            for sol_it in start_index..=4 {
                let teta = std::f64::consts::PI + modu.trunc() * 2.0 * std::f64::consts::PI;
                let x = teta - my_borne_inf;
                if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                    if (a - c + e).abs() <= eps {
                        let mut flag4 = false;
                        let mut j = 0;
                        let mut k = 1;
                        while k <= self.nb_sol {
                            j = k;
                            if teta < self.sol[k as usize] {
                                flag4 = true;
                                break;
                            }
                            if sol_it == start_index
                                && (teta - self.sol[k as usize]).abs() <= eps
                            {
                                return;
                            }
                            k += 1;
                        }

                        if !flag4 {
                            self.nb_sol += 1;
                            self.sol[self.nb_sol as usize] = teta;
                        } else {
                            let mut k = j;
                            while k <= self.nb_sol {
                                let i = self.nb_sol - k + j;
                                self.sol[(i + 1) as usize] = self.sol[i as usize];
                                k += 1;
                            }
                            self.sol[j as usize] = teta;
                            self.nb_sol += 1;
                        }
                    }
                }
            }
        }
    }
}
