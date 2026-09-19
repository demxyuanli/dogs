//! Helpers shared by the `IntAna2d` intersections.
//! Source: `IntAna2d_Outils.hxx` + `IntAna2d_Outils.cxx`.
use crate::gp::GpAx2d;
use crate::math_direct_poly_roots::DirectPolynomialRoots;
use crate::precision::epsilon;

use super::int_point::IntAna2dIntPoint;

/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

/// `MyDirectPolynomialRoots` (`hxx:25-47`). The two OCCT constructors are
/// exposed as [`Self::new5`] and [`Self::new3`].
#[derive(Debug, Clone, Copy)]
pub struct MyDirectPolynomialRoots {
    sol: [f64; 16],
    val: [f64; 16],
    nbsol: i32,
    same: bool,
}

impl MyDirectPolynomialRoots {
    /// Ctor `(A4, A3, A2, A1, A0)` (`cxx:21-214`).
    #[allow(clippy::too_many_arguments)]
    pub fn new5(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Self {
        let mut r = Self {
            sol: [REAL_LAST; 16],
            val: [REAL_LAST; 16],
            nbsol: 0,
            same: false,
        };

        let a_a = [a0.abs(), a1.abs(), a2.abs(), a3.abs(), a4.abs()];
        if a_a[0] + a_a[1] + a_a[2] + a_a[3] + a_a[4] < epsilon(10000.0) {
            r.same = true;
            return r;
        }

        let mut pb_possible = false;
        let mut nbsol_poly_complet = 0i32;

        let math_a43210 = DirectPolynomialRoots::new5(a4, a3, a2, a1, a0);
        if math_a43210.is_done() {
            let nbp = math_a43210.nb_solutions();
            nbsol_poly_complet = nbp;
            let tol = epsilon(100.0);
            for i in 1..=nbp {
                let x = math_a43210.value(i);
                r.val[r.nbsol as usize] = a0 + x * (a1 + x * (a2 + x * (a3 + x * a4)));
                r.sol[r.nbsol as usize] = x;
                let v = r.val[r.nbsol as usize];
                if v > tol || v < -tol {
                    pb_possible = true;
                }
                r.nbsol += 1;
            }
            if nbp & 1 != 0 {
                pb_possible = true;
            }
        } else {
            pb_possible = true;
        }

        if pb_possible {
            let mut an_a_min = REAL_LAST;
            let mut an_a_max = -1.0f64;
            let an_eps0 = f64::EPSILON;
            for &c in a_a.iter() {
                an_a_min = an_a_min.min(c.max(an_eps0));
                an_a_max = an_a_max.max(c.max(an_eps0));
            }
            let an_eps = 1.0e-4f64.min(epsilon(1000.0 * an_a_max / an_a_min));

            let math_a4321 = DirectPolynomialRoots::new4(a4, a3, a2, a1);
            if math_a4321.is_done() {
                let nbp = math_a4321.nb_solutions();
                for i in 1..=nbp {
                    let x = math_a4321.value(i);
                    let mut add = true;
                    for j in 0..r.nbsol {
                        if (r.sol[j as usize] - x).abs() < an_eps {
                            add = false;
                        }
                    }
                    if add {
                        r.val[r.nbsol as usize] = a0 + x * (a1 + x * (a2 + x * (a3 + x * a4)));
                        r.sol[r.nbsol as usize] = x;
                        r.nbsol += 1;
                    }
                }
            }

            let math_a3210 = DirectPolynomialRoots::new4(a3, a2, a1, a0);
            if math_a3210.is_done() {
                let nbp = math_a3210.nb_solutions();
                for i in 1..=nbp {
                    let x = math_a3210.value(i);
                    let mut add = true;
                    for j in 0..r.nbsol {
                        if (r.sol[j as usize] - x).abs() < an_eps {
                            add = false;
                        }
                    }
                    if add {
                        r.val[r.nbsol as usize] = a0 + x * (a1 + x * (a2 + x * (a3 + x * a4)));
                        r.sol[r.nbsol as usize] = x;
                        r.nbsol += 1;
                    }
                }
            }

            let math_a210 = DirectPolynomialRoots::new3(a3, a2, a1);
            if math_a210.is_done() {
                let nbp = math_a210.nb_solutions();
                for i in 1..=nbp {
                    let x = math_a210.value(i);
                    let mut add = true;
                    for j in 0..r.nbsol {
                        if (r.sol[j as usize] - x).abs() < an_eps {
                            add = false;
                        }
                    }
                    if add {
                        r.val[r.nbsol as usize] = a0 + x * (a1 + x * (a2 + x * (a3 + x * a4)));
                        r.sol[r.nbsol as usize] = x;
                        r.nbsol += 1;
                    }
                }
            }

            // Sort by increasing |val| (bubble passes, as in `cxx:177-195`).
            loop {
                let mut tri_ok = true;
                for i in 1..r.nbsol {
                    let (i, i1) = (i as usize, (i - 1) as usize);
                    if r.val[i].abs() < r.val[i1].abs() {
                        r.val.swap(i, i1);
                        r.sol.swap(i, i1);
                        tri_ok = false;
                    }
                }
                if tri_ok {
                    break;
                }
            }

            // Keep the leading roots with a small residual, at least as many as
            // the complete polynomial produced (`cxx:200-203`). The extra bound
            // on the `sol`/`val` arrays avoids the out-of-range read the OCCT
            // loop would perform if every one of the (at most 12) collected
            // roots passed the residual test.
            r.nbsol = 0;
            while r.nbsol < nbsol_poly_complet
                || (r.nbsol < 16 && r.val[r.nbsol as usize].abs() < epsilon(10000.0))
            {
                r.nbsol += 1;
            }
        }

        if r.nbsol == 0 {
            r.nbsol = -1;
        }
        if r.nbsol > 4 {
            r.same = true;
            r.nbsol = 0;
        }
        r
    }

    /// Ctor `(A2, A1, A0)` (`cxx:216-248`).
    pub fn new3(a2: f64, a1: f64, a0: f64) -> Self {
        let mut r = Self {
            sol: [REAL_LAST; 16],
            val: [REAL_LAST; 16],
            nbsol: 0,
            same: false,
        };
        if (a2.abs() + a1.abs() + a0.abs()) < epsilon(10000.0) {
            r.same = true;
            return r;
        }
        let math_a210 = DirectPolynomialRoots::new3(a2, a1, a0);
        if math_a210.is_done() {
            for i in 1..=math_a210.nb_solutions() {
                let x = math_a210.value(i);
                r.val[r.nbsol as usize] = a0 + x * (a1 + x * a2);
                r.sol[r.nbsol as usize] = x;
                r.nbsol += 1;
            }
        } else {
            r.nbsol = -1;
        }
        r
    }

    /// `NbSolutions` (`hxx:34-36`).
    pub fn nb_solutions(&self) -> i32 {
        self.nbsol
    }

    /// `Value(i)` (`hxx:38`), 1-based.
    pub fn value(&self, i: i32) -> f64 {
        self.sol[(i - 1) as usize]
    }

    /// `IsDone` (`hxx:40-43`).
    pub fn is_done(&self) -> bool {
        self.nbsol > -1
    }

    /// `InfiniteRoots` (`hxx:45-47`).
    pub fn infinite_roots(&self) -> bool {
        self.same
    }
}

/// `Points_Confondus` (`cxx:250-260`).
pub fn points_confondus(x1: f64, y1: f64, x2: f64, y2: f64) -> bool {
    if (x1 - x2).abs() < epsilon(x1) {
        if (y1 - y2).abs() < epsilon(y1) {
            return true;
        }
    }
    false
}

/// `Traitement_Points_Confondus` (`cxx:266-297`): coincident points are removed
/// and `nb_pts` is updated.
pub fn traitement_points_confondus(nb_pts: &mut i32, pts: &mut [IntAna2dIntPoint; 4]) {
    let mut i = *nb_pts;
    while i > 1 {
        let mut non_egalite = true;
        let mut j = i - 1;
        while j > 0 && non_egalite {
            if points_confondus(
                pts[(i - 1) as usize].value().x(),
                pts[(i - 1) as usize].value().y(),
                pts[(j - 1) as usize].value().x(),
                pts[(j - 1) as usize].value().y(),
            ) {
                non_egalite = false;
                let mut k = i;
                while k < *nb_pts {
                    let xk = pts[k as usize].value().x();
                    let yk = pts[k as usize].value().y();
                    let uk = pts[k as usize].param_on_first();
                    pts[(k - 1) as usize].set_value3(xk, yk, uk);
                    k += 1;
                }
                *nb_pts -= 1;
            }
            j -= 1;
        }
        i -= 1;
    }
}

/// `Coord_Ancien_Repere` (`cxx:302-321`): converts coordinates expressed in the
/// frame `Dir1` back to the "absolute" frame.
pub fn coord_ancien_repere(x1: &mut f64, y1: &mut f64, dir1: &GpAx2d) {
    let t11 = dir1.vdir.x;
    let t21 = dir1.vdir.y;
    let t13 = dir1.loc.x();
    let t23 = dir1.loc.y();
    let t22 = t11;
    let t12 = -t21;

    let x0 = t11 * *x1 + t12 * *y1 + t13;
    let y0 = t21 * *x1 + t22 * *y1 + t23;

    *x1 = x0;
    *y1 = y0;
}
