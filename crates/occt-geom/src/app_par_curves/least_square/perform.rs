//! The four `AppParCurves_LeastSquare::Perform` entry points
//! (`AppParCurves_LeastSquare.gxx:510-1055`).
//!
//! `Perform(Parameters)` keeps OCCT's two-stage shape: the unconstrained
//! contact block runs when **neither** constraint is exactly `TangencyPoint`
//! (`gxx:528`), and the tangent-augmented block always runs after it
//! (`gxx:626-735`), so a `CurvaturePoint` arm still performs the first stage
//! exactly as OCCT does.

use occt_math::{
    dactcl_decompose, dactcl_solve, Householder, MathIntVector, MathMatrix, MathStatus, MathVector,
};

use super::*;

impl LeastSquare {
    /// `Perform(Parameters)` (`gxx:510-735`).
    pub fn perform_parameters(&mut self, parameters: &MathVector) {
        self.done = false;
        if !self.isready {
            return;
        }
        let ninc1 = self.ninc - 1;
        self.iscalculated = false;

        // Calculation of matrix A and DA of approximation functions:
        self.compute_function(parameters);

        let first_exact_t = self.first_constraint == Constraint::TangencyPoint;
        let last_exact_t = self.last_constraint == Constraint::TangencyPoint;
        if !first_exact_t && !last_exact_t {
            let should_return = self.perform_no_tangency_arms();
            if should_return {
                return;
            }

            // resolution:
            let nincx = self.resfin - self.resinit + 1;
            if nincx < 1 {
                self.done = true;
                return;
            }
            let mut index = MathIntVector::new(1, nincx as usize);
            self.search_index(&mut index);
            let mut mytab = MathMatrix::new(
                self.resinit.max(0) as usize,
                self.resfin.max(0) as usize,
                1,
                self.b2.col_count(),
            );
            let mut the_aa = MathVector::new(1, index.value(nincx as usize) as usize);
            let mut my_tabb = MathVector::new(1, nincx as usize);

            self.make_taa_matrix(&mut the_aa, &mut mytab);
            dactcl_decompose(&mut the_aa, &index, 1.0e-20);

            for j in 1..=self.b2.col_count() as i32 {
                let mut kk2 = 1;
                for i in self.resinit..=self.resfin {
                    my_tabb.set_value(kk2, mytab.value(i as usize, j as usize));
                    kk2 += 1;
                }
                dactcl_solve(&the_aa, &mut my_tabb, &index, 1.0e-20);

                let mut i2 = 1;
                for k in self.resinit..=self.resfin {
                    self.mypoles.set_value(k as usize, j as usize, my_tabb.value(i2));
                    i2 += 1;
                }
            }
            self.done = true;
        }

        // ===========================================================
        // cas de tangence:
        // ===========================================================
        let nincx = self.resfin - self.resinit + 1;
        let nincx2 = 2 * nincx;

        let mut internal_index = MathIntVector::new(1, nincx as usize);
        self.search_index(&mut internal_index);
        let mut index = MathIntVector::new(1, self.ninc as usize);

        let mut l = 1;
        if self.resinit <= self.resfin {
            for j in 0..=self.na - 1 {
                let deport = j * internal_index.value(nincx as usize);
                for i in 1..=nincx {
                    index.set_value(l, internal_index.value(i as usize) + deport);
                    l += 1;
                }
            }
        }

        if self.resinit > self.resfin {
            index.set_value(1, 1);
        }
        let both_t = at_least_tangency(self.first_constraint) && at_least_tangency(self.last_constraint);
        if ninc1 > 1 && both_t {
            let v = index.value((ninc1 - 1) as usize) + ninc1;
            index.set_value(ninc1 as usize, v);
        }
        if at_least_tangency(self.first_constraint) || at_least_tangency(self.last_constraint) {
            let v = index.value((self.ninc - 1) as usize) + self.ninc;
            index.set_value(self.ninc as usize, v);
        }

        let mut the_a = MathVector::new(1, index.value(self.ninc as usize) as usize);
        let mut my_tab = MathVector::new(1, self.ninc as usize);

        self.make_taa_vector(&mut the_a, &mut my_tab);

        dactcl_decompose(&mut the_a, &index, 1.0e-20);
        let error = dactcl_solve(&the_a, &mut my_tab, &index, 1.0e-20);

        if error == MathStatus::Ok {
            self.done = true;
        }

        if both_t {
            self.lambda1 = my_tab.value(ninc1 as usize);
            self.lambda2 = my_tab.value(self.ninc as usize);
        } else if at_least_tangency(self.first_constraint) {
            self.lambda1 = my_tab.value(self.ninc as usize);
        } else if at_least_tangency(self.last_constraint) {
            self.lambda2 = my_tab.value(self.ninc as usize);
        }

        // The results are stored in mypoles.
        // =========================================
        let nbpoles = self.nbpoles;
        let nbpol1 = nbpoles - 1;
        let first_t = at_least_tangency(self.first_constraint);
        let last_t = at_least_tangency(self.last_constraint);
        let mut k = 1;
        let mut i2: i32 = 1;
        for _ci in 1..=self.nbp {
            let k1 = k + 1;
            let k2 = k + 2;
            for j in self.resinit..=self.resfin {
                self.mypoles.set_value(j as usize, k as usize, my_tab.value(i2 as usize));
                self.mypoles.set_value(
                    j as usize,
                    k1 as usize,
                    my_tab.value((i2 + nincx) as usize),
                );
                self.mypoles.set_value(
                    j as usize,
                    k2 as usize,
                    my_tab.value((i2 + nincx2) as usize),
                );
                i2 += 1;
            }

            if first_t {
                let p = |s: &Self, off: i32| {
                    s.mypoints.value(s.myfirstp as usize, off as usize) + s.lambda1 * s.vec1t.value(off as usize)
                };
                self.mypoles.set_value(2, k as usize, p(self, k));
                self.mypoles.set_value(2, k1 as usize, p(self, k1));
                self.mypoles.set_value(2, k2 as usize, p(self, k2));
            }
            if last_t {
                let p = |s: &Self, off: i32| {
                    s.mypoints.value(s.mylastp as usize, off as usize) - s.lambda2 * s.vec2t.value(off as usize)
                };
                self.mypoles.set_value(nbpol1 as usize, k as usize, p(self, k));
                self.mypoles.set_value(nbpol1 as usize, k1 as usize, p(self, k1));
                self.mypoles.set_value(nbpol1 as usize, k2 as usize, p(self, k2));
            }
            k += 3;
            i2 += nincx2;
        }

        for _ci in 1..=self.nbp2d {
            let k1 = k + 1;
            for j in self.resinit..=self.resfin {
                self.mypoles.set_value(j as usize, k as usize, my_tab.value(i2 as usize));
                self.mypoles.set_value(
                    j as usize,
                    k1 as usize,
                    my_tab.value((i2 + nincx) as usize),
                );
                i2 += 1;
            }
            if first_t {
                let p = |s: &Self, off: i32| {
                    s.mypoints.value(s.myfirstp as usize, off as usize) + s.lambda1 * s.vec1t.value(off as usize)
                };
                self.mypoles.set_value(2, k as usize, p(self, k));
                self.mypoles.set_value(2, k1 as usize, p(self, k1));
            }
            if last_t {
                let p = |s: &Self, off: i32| {
                    s.mypoints.value(s.mylastp as usize, off as usize) - s.lambda2 * s.vec2t.value(off as usize)
                };
                self.mypoles.set_value(nbpol1 as usize, k as usize, p(self, k));
                self.mypoles.set_value(nbpol1 as usize, k1 as usize, p(self, k1));
            }
            k += 2;
            i2 += nincx;
        }
    }

    /// The `B2` residual arms of `Perform(Parameters)` (`gxx:531-583`).
    ///
    /// Returns `true` when the `NoConstraint` / `NoConstraint` arm ran, which
    /// returns straight out of `Perform` in OCCT (`gxx:536-545`).
    fn perform_no_tangency_arms(&mut self) -> bool {
        let nbpoles = self.nbpoles;
        let first_no = self.first_constraint == Constraint::NoConstraint;
        let first_pass = self.first_constraint == Constraint::PassPoint;
        let last_no = self.last_constraint == Constraint::NoConstraint;
        let last_pass = self.last_constraint == Constraint::PassPoint;

        if first_no && last_no {
            let hh = Householder::new(&self.a, &self.mypoints);
            if !hh.is_done() {
                self.done = false;
                return true;
            }
            self.done = true;
            let sol = hh.all_values();
            for i in 1..=self.mypoles.row_count() {
                for j in 1..=self.mypoles.col_count() {
                    self.mypoles.set_value(i, j, sol.value(i, j));
                }
            }
            return true;
        }

        if first_no && last_pass {
            for j in self.first_p..=self.last_p {
                let ad1 = self.a.value(j as usize, nbpoles as usize);
                for i in 1..=self.b2.col_count() as i32 {
                    let v = self.mypoints.value(j as usize, i as usize)
                        - ad1 * self.mypoles.value(nbpoles as usize, i as usize);
                    self.b2.set_value(j as usize, i as usize, v);
                }
            }
        } else if first_pass && last_no {
            for j in self.first_p..=self.last_p {
                let a0 = self.a.value(j as usize, 1);
                for i in 1..=self.b2.col_count() as i32 {
                    let v = self.mypoints.value(j as usize, i as usize)
                        - a0 * self.mypoles.value(1, i as usize);
                    self.b2.set_value(j as usize, i as usize, v);
                }
            }
        } else if first_pass && last_pass {
            for j in self.first_p..=self.last_p {
                let a0 = self.a.value(j as usize, 1);
                let ad1 = self.a.value(j as usize, nbpoles as usize);
                for i in 1..=self.b2.col_count() as i32 {
                    let v = self.mypoints.value(j as usize, i as usize)
                        - a0 * self.mypoles.value(1, i as usize)
                        - ad1 * self.mypoles.value(nbpoles as usize, i as usize);
                    self.b2.set_value(j as usize, i as usize, v);
                }
            }
        }
        false
    }

    /// `Perform(Parameters, V1t, V2t, l1, l2)` (`gxx:737-760`).
    pub fn perform_tangents(
        &mut self,
        parameters: &MathVector,
        v1t: &MathVector,
        v2t: &MathVector,
        l1: f64,
        l2: f64,
    ) {
        self.done = false;
        if !self.isready {
            return;
        }
        self.resinit = 3;
        self.resfin = self.nbpoles - 2;
        let nincx = self.resfin - self.resinit + 1;
        self.ninc = self.na * nincx + 2;
        self.first_constraint = Constraint::TangencyPoint;
        self.last_constraint = Constraint::TangencyPoint;
        copy_shifted(&mut self.vec1t, v1t);
        copy_shifted(&mut self.vec2t, v2t);
        self.perform_parameters_lambda(parameters, l1, l2);
    }

    /// `Perform(Parameters, V1t, V2t, V1c, V2c, l1, l2)` (`gxx:762-793`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_tangents_curvature(
        &mut self,
        parameters: &MathVector,
        v1t: &MathVector,
        v2t: &MathVector,
        v1c: &MathVector,
        v2c: &MathVector,
        l1: f64,
        l2: f64,
    ) {
        self.done = false;
        if !self.isready {
            return;
        }
        self.resinit = 4;
        self.resfin = self.nbpoles - 3;
        let nincx = self.resfin - self.resinit + 1;
        self.ninc = self.na * nincx + 2;
        self.first_constraint = Constraint::CurvaturePoint;
        self.last_constraint = Constraint::CurvaturePoint;

        copy_shifted(&mut self.vec1t, v1t);
        copy_shifted(&mut self.vec2t, v2t);
        copy_shifted(&mut self.vec1c, v1c);
        copy_shifted(&mut self.vec2c, v2c);
        self.perform_parameters_lambda(parameters, l1, l2);
    }

    /// `Perform(Parameters, l1, l2)` (`gxx:795-1055`).
    pub fn perform_parameters_lambda(&mut self, parameters: &MathVector, l1: f64, l2: f64) {
        self.done = false;
        if !self.isready {
            return;
        }
        if !at_least_tangency(self.first_constraint) && !at_least_tangency(self.last_constraint) {
            self.perform_parameters(parameters);
            return;
        }
        self.iscalculated = false;

        self.lambda1 = l1;
        self.lambda2 = l2;
        let nbpoles = self.nbpoles;
        let deg = self.deg;
        let nbpol1 = nbpoles - 1;
        let l11 = deg as f64 * l1;
        let l22 = deg as f64 * l2;

        self.compute_function(parameters);

        let first_t = at_least_tangency(self.first_constraint);
        let last_t = at_least_tangency(self.last_constraint);

        if first_t {
            for i in 1..=self.mypoles.col_count() as i32 {
                let v = self.mypoints.value(self.myfirstp as usize, i as usize)
                    + l1 * self.vec1t.value(i as usize);
                self.mypoles.set_value(2, i as usize, v);
            }
        }

        if self.first_constraint == Constraint::CurvaturePoint {
            for i in 1..=self.mypoles.col_count() as i32 {
                let v = 2.0 * self.mypoles.value(2, i as usize) - self.mypoles.value(1, i as usize)
                    + l11 * l11 * self.vec1c.value(i as usize) / (deg * (deg - 1)) as f64;
                self.mypoles.set_value(3, i as usize, v);
            }
        }

        if last_t {
            for i in 1..=self.mypoles.col_count() as i32 {
                let v = self.mypoints.value(self.mylastp as usize, i as usize)
                    - l2 * self.vec2t.value(i as usize);
                self.mypoles.set_value(nbpol1 as usize, i as usize, v);
            }
        }

        if self.last_constraint == Constraint::CurvaturePoint {
            for i in 1..=self.mypoles.col_count() as i32 {
                let v = 2.0 * self.mypoles.value(nbpol1 as usize, i as usize)
                    - self.mypoles.value(nbpoles as usize, i as usize)
                    + l22 * l22 * self.vec2c.value(i as usize) / (deg * (deg - 1)) as f64;
                self.mypoles.set_value((nbpoles - 2) as usize, i as usize, v);
            }
        }

        if self.resinit > self.resfin {
            self.done = true;
            return;
        }

        // The fixed-pole subtractors, one arm per (First, Last) pair exactly
        // as `gxx:852-1020`.
        self.subtract_fixed_poles();

        let nincx = self.resfin - self.resinit + 1;

        let mut mytab = MathMatrix::new(
            self.resinit.max(0) as usize,
            self.resfin.max(0) as usize,
            1,
            self.b2.col_count(),
        );
        let mut index = MathIntVector::new(1, nincx as usize);
        self.search_index(&mut index);
        let mut aa = MathVector::new(1, index.value(nincx as usize) as usize);
        self.make_taa_matrix(&mut aa, &mut mytab);

        let mut my_tabb = MathVector::new(1, nincx as usize);

        dactcl_decompose(&mut aa, &index, 1.0e-20);

        for j in 1..=self.b2.col_count() as i32 {
            let mut kk2 = 1;
            for i in self.resinit..=self.resfin {
                my_tabb.set_value(kk2, mytab.value(i as usize, j as usize));
                kk2 += 1;
            }

            dactcl_solve(&aa, &mut my_tabb, &index, 1.0e-20);

            let mut i2 = 1;
            for k in self.resinit..=self.resfin {
                self.mypoles.set_value(k as usize, j as usize, my_tabb.value(i2));
                i2 += 1;
            }
        }

        self.done = true;
    }

    /// The `B2(j, i) = mypoints - sum(A(j, idx) * mypoles(idx, i))` arms of
    /// `Perform(Parameters, l1, l2)` (`gxx:852-1020`).
    fn subtract_fixed_poles(&mut self) {
        let nbpoles = self.nbpoles;
        let fc = self.first_constraint;
        let lc = self.last_constraint;
        for j in self.first_p..=self.last_p {
            let a = |k: i32| self.a.value(j as usize, k as usize);
            let terms: Vec<(i32, f64)> = match (fc, lc) {
                (Constraint::NoConstraint, Constraint::TangencyPoint) => {
                    vec![(nbpoles, a(nbpoles)), (nbpoles - 1, a(nbpoles - 1))]
                }
                (Constraint::NoConstraint, Constraint::CurvaturePoint) => vec![
                    (nbpoles, a(nbpoles)),
                    (nbpoles - 1, a(nbpoles - 1)),
                    (nbpoles - 2, a(nbpoles - 2)),
                ],
                (Constraint::PassPoint, Constraint::TangencyPoint) => {
                    vec![(1, a(1)), (nbpoles, a(nbpoles)), (nbpoles - 1, a(nbpoles - 1))]
                }
                (Constraint::PassPoint, Constraint::CurvaturePoint) => vec![
                    (1, a(1)),
                    (nbpoles, a(nbpoles)),
                    (nbpoles - 1, a(nbpoles - 1)),
                    (nbpoles - 2, a(nbpoles - 2)),
                ],
                (Constraint::TangencyPoint, Constraint::NoConstraint) => {
                    vec![(1, a(1)), (2, a(2))]
                }
                (Constraint::TangencyPoint, Constraint::PassPoint) => {
                    vec![(1, a(1)), (nbpoles, a(nbpoles)), (2, a(2))]
                }
                (Constraint::TangencyPoint, Constraint::TangencyPoint) => vec![
                    (1, a(1)),
                    (nbpoles, a(nbpoles)),
                    (2, a(2)),
                    (nbpoles - 1, a(nbpoles - 1)),
                ],
                (Constraint::CurvaturePoint, Constraint::NoConstraint) => {
                    vec![(1, a(1)), (2, a(2)), (3, a(3))]
                }
                (Constraint::CurvaturePoint, Constraint::PassPoint) => {
                    vec![(1, a(1)), (2, a(2)), (3, a(3)), (nbpoles, a(nbpoles))]
                }
                (Constraint::CurvaturePoint, Constraint::TangencyPoint) => vec![
                    (1, a(1)),
                    (2, a(2)),
                    (3, a(3)),
                    (nbpoles, a(nbpoles)),
                    (nbpoles - 1, a(nbpoles - 1)),
                ],
                (Constraint::CurvaturePoint, Constraint::CurvaturePoint) => vec![
                    (1, a(1)),
                    (2, a(2)),
                    (3, a(3)),
                    (nbpoles, a(nbpoles)),
                    (nbpoles - 1, a(nbpoles - 1)),
                    (nbpoles - 2, a(nbpoles - 2)),
                ],
                // Not reachable: `Perform(Parameters, l1, l2)` falls back to
                // `Perform(Parameters)` unless one end is at least a tangency
                // point (`gxx:800-804`).
                _ => Vec::new(),
            };
            for i in 1..=self.b2.col_count() as i32 {
                let mut v = self.mypoints.value(j as usize, i as usize);
                for (idx, coeff) in &terms {
                    v -= coeff * self.mypoles.value(*idx as usize, i as usize);
                }
                self.b2.set_value(j as usize, i as usize, v);
            }
        }
    }
}

