//! `AppParCurves::BernsteinMatrix`, `AppParCurves::Bernstein`,
//! `AppParCurves::SecondDerivativeBernstein` — the Bezier basis matrices the
//! least-squares resolvers fill.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves.hxx:20-48` and
//! `AppParCurves.cxx:25-157` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The OCCT local `math_Vector B(1, NbPoles - 1)` is a `vec` of
//!   `nb_poles - 1` slots used 1-based through a local `b(i)` helper; the
//!   out-of-range case (`nb_poles < 2`, an empty OCCT array read) panics with
//!   `Standard_OutOfRange` instead of reading past the array.
//! - `math_Vector U` / `math_Matrix A` / `DA` keep their own bounds, exactly
//!   like OCCT: the loop walks `U.Lower()..=U.Upper()` and writes the same row
//!   index of `A` / `DA`.
//! - `SecondDerivativeBernstein` keeps the `deg == 1` / `deg == 2` special
//!   cases of `cxx:115-126`; the general arm is `cxx:128-156`.

use occt_math::{MathMatrix, MathVector};

/// `AppParCurves` static helpers (`AppParCurves.hxx:20-48`).
pub struct AppParCurves;

impl AppParCurves {
    /// 1-based read of the local `B` array, `Standard_OutOfRange` when out of
    /// bounds (`math_Vector` accesses behave the same way).
    fn b(b: &[f64], i: usize) -> f64 {
        match i.checked_sub(1).and_then(|k| b.get(k)) {
            Some(v) => *v,
            None => panic!("Standard_OutOfRange: AppParCurves::Bernstein"),
        }
    }

    /// 1-based write of the local `B` array.
    fn set_b(b: &mut [f64], i: usize, v: f64) {
        match i.checked_sub(1).and_then(|k| b.get_mut(k)) {
            Some(slot) => *slot = v,
            None => panic!("Standard_OutOfRange: AppParCurves::Bernstein"),
        }
    }

    /// `AppParCurves::BernsteinMatrix(NbPoles, U, A)` (`cxx:25-58`).
    pub fn bernstein_matrix(nb_poles: usize, u: &MathVector, a: &mut MathMatrix) {
        let first = u.lower();
        let last = u.upper();
        let mut b = vec![0.0; nb_poles.saturating_sub(1)];
        for i in first..=last {
            let u0 = u.value(i);
            let u1 = 1.0 - u0;
            Self::set_b(&mut b, 1, 1.0);
            let mut id = 2;
            while id <= nb_poles - 1 {
                let mut y0 = Self::b(&b, 1);
                let mut y1 = u0 * y0;
                Self::set_b(&mut b, 1, y0 - y1);
                let mut j = 2;
                while j <= id - 1 {
                    let xs = y1;
                    y0 = Self::b(&b, j);
                    y1 = u0 * y0;
                    Self::set_b(&mut b, j, y0 - y1 + xs);
                    j += 1;
                }
                Self::set_b(&mut b, id, y1);
                id += 1;
            }
            a.set_value(i, 1, u1 * Self::b(&b, 1));
            a.set_value(i, nb_poles, u0 * Self::b(&b, nb_poles - 1));
            for j in 2..=nb_poles - 1 {
                a.set_value(i, j, u1 * Self::b(&b, j) + u0 * Self::b(&b, j - 1));
            }
        }
    }

    /// `AppParCurves::Bernstein(NbPoles, U, A, DA)` (`cxx:60-103`).
    pub fn bernstein(nb_poles: usize, u: &MathVector, a: &mut MathMatrix, da: &mut MathMatrix) {
        let n_deg = (nb_poles - 1) as f64;
        let first = u.lower();
        let last = u.upper();
        let mut b = vec![0.0; nb_poles.saturating_sub(1)];
        for i in first..=last {
            let u0 = u.value(i);
            let u1 = 1.0 - u0;
            Self::set_b(&mut b, 1, 1.0);
            let mut id = 2;
            while id <= nb_poles - 1 {
                let mut y0 = Self::b(&b, 1);
                let mut y1 = u0 * y0;
                Self::set_b(&mut b, 1, y0 - y1);
                let mut j = 2;
                while j <= id - 1 {
                    let xs = y1;
                    y0 = Self::b(&b, j);
                    y1 = u0 * y0;
                    Self::set_b(&mut b, j, y0 - y1 + xs);
                    j += 1;
                }
                Self::set_b(&mut b, id, y1);
                id += 1;
            }
            da.set_value(i, 1, -n_deg * Self::b(&b, 1));
            da.set_value(i, nb_poles, n_deg * Self::b(&b, nb_poles - 1));
            a.set_value(i, 1, u1 * Self::b(&b, 1));
            a.set_value(i, nb_poles, u0 * Self::b(&b, nb_poles - 1));
            for j in 2..=nb_poles - 1 {
                let bj = Self::b(&b, j);
                let bj1 = Self::b(&b, j - 1);
                da.set_value(i, j, n_deg * (bj1 - bj));
                a.set_value(i, j, u1 * bj + u0 * bj1);
            }
        }
    }

    /// `AppParCurves::SecondDerivativeBernstein(U, DDA)` (`cxx:105-157`).
    ///
    /// `DDA` holds `deg + 1` values, `deg = DDA.Length() - 1`.
    pub fn second_derivative_bernstein(u: f64, dda: &mut MathVector) {
        let nb_poles = dda.len();
        let deg = nb_poles - 1;
        let n4 = (deg * (deg.saturating_sub(1))) as f64;
        let mut b = vec![0.0; deg.saturating_sub(1)];
        if deg >= 2 {
            Self::set_b(&mut b, 1, 1.0);
        }

        if deg == 1 {
            dda.set_value(1, 0.0);
            dda.set_value(2, 0.0);
        } else if deg == 2 {
            dda.set_value(1, 2.0);
            dda.set_value(2, -4.0);
            dda.set_value(3, 2.0);
        } else {
            let mut id = 2;
            while id <= deg - 1 {
                let mut y0 = Self::b(&b, 1);
                let mut y1 = u * y0;
                Self::set_b(&mut b, 1, y0 - y1);
                let mut j = 2;
                while j <= id - 1 {
                    let xs = y1;
                    y0 = Self::b(&b, j);
                    y1 = u * y0;
                    Self::set_b(&mut b, j, y0 - y1 + xs);
                    j += 1;
                }
                Self::set_b(&mut b, id, y1);
                id += 1;
            }

            dda.set_value(1, n4 * Self::b(&b, 1));
            dda.set_value(2, n4 * (-2.0 * Self::b(&b, 1) + Self::b(&b, 2)));
            dda.set_value(
                deg,
                n4 * (Self::b(&b, deg - 2) - 2.0 * Self::b(&b, deg - 1)),
            );
            dda.set_value(deg + 1, n4 * Self::b(&b, deg - 1));

            for j in 2..=deg - 2 {
                dda.set_value(
                    j + 1,
                    n4 * (Self::b(&b, j - 1) - 2.0 * Self::b(&b, j) + Self::b(&b, j + 1)),
                );
            }
        }
    }
}
