//! `math_TrigonometricFunctionRoots` - solutions of
//! `A*cos(x)^2 + 2*B*cos(x)*sin(x) + C*cos(x) + D*sin(x) + E = 0`.
//! Source: `src/FoundationClasses/TKMath/math/math_TrigonometricFunctionRoots.cxx`,
//! `math_TrigonometricFunctionRoots.hxx` and
//! `math_TrigonometricFunctionRoots.lxx`.
//!
//! Used by `IntAna2d_AnaIntersection::Perform(gp_Circ2d, IntAna2d_Conic)`
//! (`IntAna2d_AnaIntersection_5.cxx:55`) and
//! `Perform(gp_Elips2d, IntAna2d_Conic)` (`IntAna2d_AnaIntersection_6.cxx:57`).

use crate::math_direct_poly_roots::DirectPolynomialRoots;
use crate::math_newton_function_root::{NewtonFunctionRoot, TrigonometricEquationFunction};
use crate::precision::{epsilon, PCONFUSION};

/// `RealFirst()` / `RealLast()` (`Standard_Real.hxx:128-130`).
const REAL_FIRST: f64 = -f64::MAX;
const REAL_LAST: f64 = f64::MAX;

const TWO_PI: f64 = std::f64::consts::PI * 2.0;

/// `math_TrigonometricFunctionRoots` (`math_TrigonometricFunctionRoots.hxx:33-102`).
#[derive(Debug, Clone, Copy)]
pub struct TrigonometricFunctionRoots {
    nb_sol: i32,
    sol: [f64; 4],
    infinite_status: bool,
    done: bool,
}

impl TrigonometricFunctionRoots {
    /// Ctor `(D, E, InfBound, SupBound)` (`cxx:32-43`).
    pub fn new_de(d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        Self::new_abcde(0.0, 0.0, 0.0, d, e, inf_bound, sup_bound)
    }

    /// Ctor `(C, D, E, InfBound, SupBound)` (`cxx:45-58`).
    pub fn new_cde(c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        Self::new_abcde(0.0, 0.0, c, d, e, inf_bound, sup_bound)
    }

    /// Ctor `(A, B, C, D, E, InfBound, SupBound)` (`cxx:60-73`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_abcde(a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        let mut r = Self {
            nb_sol: -1,
            sol: [0.0; 4],
            infinite_status: false,
            done: false,
        };
        r.perform(a, b, c, d, e, inf_bound, sup_bound);
        r
    }

    /// `IsDone` (`lxx:23-27`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `InfiniteRoots` (`lxx:19-22`).
    pub fn infinite_roots(&self) -> bool {
        self.infinite_status
    }

    /// `Value` (`lxx:35-44`); 1-based, `StdFail_InfiniteSolutions`,
    /// `StdFail_NotDone` and `Standard_OutOfRange` all raise.
    pub fn value(&self, index: i32) -> f64 {
        assert!(
            !self.infinite_status,
            "StdFail_InfiniteSolutions in math_TrigonometricFunctionRoots::Value"
        );
        assert!(
            self.done,
            "StdFail_NotDone in math_TrigonometricFunctionRoots::Value"
        );
        assert!(
            index <= self.nb_sol,
            "Standard_OutOfRange in math_TrigonometricFunctionRoots::Value"
        );
        self.sol[(index - 1) as usize]
    }

    /// `NbSolutions` (`lxx:46-51`).
    pub fn nb_solutions(&self) -> i32 {
        assert!(
            !self.infinite_status,
            "StdFail_InfiniteSolutions in math_TrigonometricFunctionRoots::NbSolutions"
        );
        assert!(
            self.done,
            "StdFail_NotDone in math_TrigonometricFunctionRoots::NbSolutions"
        );
        self.nb_sol
    }

    /// `Perform` (`cxx:75-551`).
    #[allow(clippy::too_many_arguments)]
    fn perform(&mut self, a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) {
        let mut i;
        let mut j = 0usize;
        let mut n_zer: i32;
        let nit = 10;
        let my_borne_inf;
        let mut delta;
        let mod_;
        let mut teta;
        let mut x;
        let tol1 = 1.0e-15;
        let mut ko = [0.0f64; 5];
        let mut zer = [0.0f64; 4];
        let mut flag4;

        self.infinite_status = false;
        self.done = true;

        let eps = 1.5e-12;

        let depi = TWO_PI;
        if inf_bound <= REAL_FIRST && sup_bound >= REAL_LAST {
            my_borne_inf = 0.0;
            delta = depi;
            mod_ = 0.0;
        } else if sup_bound >= REAL_LAST {
            my_borne_inf = inf_bound;
            delta = depi;
            mod_ = my_borne_inf / depi;
        } else if inf_bound <= REAL_FIRST {
            my_borne_inf = sup_bound - depi;
            delta = depi;
            mod_ = my_borne_inf / depi;
        } else {
            my_borne_inf = inf_bound;
            delta = sup_bound - inf_bound;
            mod_ = inf_bound / depi;
            if (sup_bound - inf_bound) > depi {
                delta = depi;
            }
        }

        if a.abs() <= eps && b.abs() <= eps {
            if c.abs() <= eps {
                if d.abs() <= eps {
                    if e.abs() <= eps {
                        // Infinite number of solutions.
                        self.infinite_status = true;
                        return;
                    } else {
                        self.nb_sol = 0;
                        return;
                    }
                } else {
                    // Equation of the type d*sin(x) + e = 0
                    self.nb_sol = 0;
                    let aa = -e / d;
                    if aa.abs() > 1.0 {
                        return;
                    }

                    zer[0] = aa.asin();
                    zer[1] = std::f64::consts::PI - zer[0];
                    n_zer = 2;
                    for i in 1..=n_zer as usize {
                        if zer[i - 1] <= -eps {
                            zer[i - 1] = depi - zer[i - 1].abs();
                        }
                        // Bring the solutions between InfBound and SupBound.
                        zer[i - 1] += mod_.trunc() * depi;
                        x = zer[i - 1] - my_borne_inf;
                        if x > -epsilon(delta) && x < delta + epsilon(delta) {
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = zer[i - 1];
                        }
                    }
                    return;
                }
            } else if d.abs() <= eps {
                // First degree equation of the form c*cos(x) + e = 0
                self.nb_sol = 0;
                let aa = -e / c;
                if aa.abs() > 1.0 {
                    return;
                }
                zer[0] = aa.acos();
                zer[1] = -zer[0];
                n_zer = 2;

                for i in 1..=n_zer as usize {
                    if zer[i - 1] <= -eps {
                        zer[i - 1] = depi - zer[i - 1].abs();
                    }
                    // Bring the solutions between InfBound and SupBound.
                    zer[i - 1] += mod_.trunc() * TWO_PI;
                    x = zer[i - 1] - my_borne_inf;
                    if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                        self.nb_sol += 1;
                        self.sol[(self.nb_sol - 1) as usize] = zer[i - 1];
                    }
                }
                return;
            } else {
                // Second degree equation.
                let aa = e - c;
                let bb = 2.0 * d;
                let cc = e + c;

                let resol = DirectPolynomialRoots::new3(aa, bb, cc);
                if !resol.is_done() {
                    self.done = false;
                    return;
                } else if !resol.infinite_roots() {
                    n_zer = resol.nb_solutions();
                    for i in 1..=n_zer as usize {
                        zer[i - 1] = resol.value(i as i32);
                    }
                } else {
                    self.infinite_status = true;
                    return;
                }
            }
        } else {
            // Two additional analytical cases.
            if a.abs() <= eps && e.abs() <= eps {
                if c.abs() <= eps {
                    // 2 * B * sin * cos + D * sin = 0
                    n_zer = 2;
                    zer[0] = 0.0;
                    zer[1] = std::f64::consts::PI;

                    let aa = -d / (b * 2.0);
                    if aa.abs() <= 1.0 + PCONFUSION {
                        n_zer = 4;
                        if aa >= 1.0 {
                            zer[2] = 0.0;
                            zer[3] = 0.0;
                        } else if aa <= -1.0 {
                            zer[2] = std::f64::consts::PI;
                            zer[3] = std::f64::consts::PI;
                        } else {
                            zer[2] = aa.acos();
                            zer[3] = depi - zer[2];
                        }
                    }

                    self.nb_sol = 0;
                    for i in 1..=n_zer as usize {
                        if zer[i - 1] <= my_borne_inf - eps {
                            zer[i - 1] += depi;
                        }
                        // Bring the solutions between InfBound and SupBound.
                        zer[i - 1] += mod_.trunc() * TWO_PI;
                        x = zer[i - 1] - my_borne_inf;
                        if x >= -PCONFUSION && x <= delta + PCONFUSION {
                            if zer[i - 1] < inf_bound {
                                zer[i - 1] = inf_bound;
                            }
                            if zer[i - 1] > sup_bound {
                                zer[i - 1] = sup_bound;
                            }
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = zer[i - 1];
                        }
                    }
                    return;
                }
                if d.abs() <= eps {
                    // 2 * B * sin * cos + C * cos = 0
                    n_zer = 2;
                    zer[0] = std::f64::consts::PI / 2.0;
                    zer[1] = std::f64::consts::PI * 3.0 / 2.0;

                    let aa = -c / (b * 2.0);
                    if aa.abs() <= 1.0 + PCONFUSION {
                        n_zer = 4;
                        if aa >= 1.0 {
                            zer[2] = std::f64::consts::PI / 2.0;
                            zer[3] = std::f64::consts::PI / 2.0;
                        } else if aa <= -1.0 {
                            zer[2] = std::f64::consts::PI * 3.0 / 2.0;
                            zer[3] = std::f64::consts::PI * 3.0 / 2.0;
                        } else {
                            zer[2] = aa.asin();
                            zer[3] = std::f64::consts::PI - zer[2];
                        }
                    }

                    self.nb_sol = 0;
                    for i in 1..=n_zer as usize {
                        if zer[i - 1] <= my_borne_inf - eps {
                            zer[i - 1] += depi;
                        }
                        // Bring the solutions between InfBound and SupBound.
                        zer[i - 1] += mod_.trunc() * TWO_PI;
                        x = zer[i - 1] - my_borne_inf;
                        if x >= -PCONFUSION && x <= delta + PCONFUSION {
                            if zer[i - 1] < inf_bound {
                                zer[i - 1] = inf_bound;
                            }
                            if zer[i - 1] > sup_bound {
                                zer[i - 1] = sup_bound;
                            }
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = zer[i - 1];
                        }
                    }
                    return;
                }
            }

            // Fourth degree equation.
            ko[0] = a - c + e;
            ko[1] = 2.0 * d - 4.0 * b;
            ko[2] = 2.0 * e - 2.0 * a;
            ko[3] = 4.0 * b + 2.0 * d;
            ko[4] = a + c + e;
            let mut bko;
            loop {
                bko = false;
                let resol4 = DirectPolynomialRoots::new5(ko[0], ko[1], ko[2], ko[3], ko[4]);
                if !resol4.is_done() {
                    self.done = false;
                    return;
                } else if !resol4.infinite_roots() {
                    n_zer = resol4.nb_solutions();
                    for i in 1..=n_zer as usize {
                        zer[i - 1] = resol4.value(i as i32);
                    }
                } else {
                    self.infinite_status = true;
                    return;
                }

                let mut triok;
                loop {
                    triok = true;
                    for i in 1..n_zer as usize {
                        if zer[i - 1] > zer[i] {
                            zer.swap(i - 1, i);
                            triok = false;
                        }
                    }
                    if triok {
                        break;
                    }
                }

                for i in 1..n_zer as usize {
                    if (zer[i] - zer[i - 1]).abs() < eps {
                        // Is it a double root or a numerical error?
                        let qw = zer[i];
                        let va = ko[3] + qw * (2.0 * ko[2] + qw * (3.0 * ko[1] + qw * (4.0 * ko[0])));
                        if va.abs() > eps {
                            bko = true;
                            break;
                        }
                    }
                }
                if bko {
                    // If there is a small coefficient, divide.
                    for c in ko.iter_mut() {
                        *c *= 0.0001;
                    }
                }
                if !bko {
                    break;
                }
            }
        }

        // Verification of the solutions against the bounds.
        let supm_infs100 = (sup_bound - inf_bound) * 0.01;
        self.nb_sol = 0;
        for i in 1..=n_zer as usize {
            teta = zer[i - 1].atan();
            teta += teta;
            if zer[i - 1] <= -eps {
                teta = depi - teta.abs();
            }
            teta += mod_.trunc() * depi;
            if teta - my_borne_inf < 0.0 {
                teta += depi;
            }

            x = teta - my_borne_inf;
            if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                x = teta;

                // Newton call.
                let mut teta_newton = teta;
                let mut my_f = TrigonometricEquationFunction::new(a, b, c, d, e);
                let resol = NewtonFunctionRoot::new(&mut my_f, x, tol1, eps, nit);
                if resol.is_done() {
                    teta_newton = resol.root();
                }
                // lbr le 7 mars 97 (Newton converges far from the initial solution)
                let delta_newton = teta_newton - teta;
                if delta_newton <= supm_infs100 && delta_newton >= -supm_infs100 {
                    teta = teta_newton;
                }

                flag4 = false;

                for k in 1..=self.nb_sol as usize {
                    // Sort the values in increasing order.
                    if teta < self.sol[k - 1] {
                        for l in k..=self.nb_sol as usize {
                            j = self.nb_sol as usize - l + k;
                            self.sol[j] = self.sol[j - 1];
                        }
                        self.sol[k - 1] = teta;
                        self.nb_sol += 1;
                        flag4 = true;
                        break;
                    }
                }
                if !flag4 {
                    self.nb_sol += 1;
                    self.sol[(self.nb_sol - 1) as usize] = teta;
                }
            }
        }

        // Special case of PI.
        if self.nb_sol < 4 {
            let start_index = self.nb_sol + 1;
            for sol_it in start_index..=4 {
                teta = std::f64::consts::PI + mod_.trunc() * TWO_PI;
                x = teta - my_borne_inf;
                if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                    if (a - c + e).abs() <= eps {
                        flag4 = false;
                        for k in 1..=self.nb_sol as usize {
                            j = k;
                            if teta < self.sol[k - 1] {
                                flag4 = true;
                                break;
                            }
                            if sol_it == start_index && (teta - self.sol[k - 1]).abs() <= eps {
                                return;
                            }
                        }

                        if !flag4 {
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = teta;
                        } else {
                            for k in j..=self.nb_sol as usize {
                                i = self.nb_sol as usize - k + j;
                                self.sol[i] = self.sol[i - 1];
                            }
                            self.sol[j - 1] = teta;
                            self.nb_sol += 1;
                        }
                    }
                }
            }
        }
    }
}
