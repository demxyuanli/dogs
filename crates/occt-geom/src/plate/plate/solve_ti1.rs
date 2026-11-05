//! `Plate_Plate::SolveTI1` — only `PinpointConstraint`s are loaded.
//!
//! Source: `Plate_Plate.cxx:370-506`.
//!
//! Deviations (no control-flow change):
//! - `Message_ProgressScope` / `aScope.UserBreak()` (`cxx:443-466`) is
//!   UNPORTED: there is no progress range in the port.
//! - `algo_gauss.Solve(...)` raises `StdFail_NotDone` in OCCT when the
//!   decomposition failed; the port sets `ok = false` and returns.

use occt_core::gp::GpXyz;
use occt_math::{MathGauss, MathMatrix, MathVector};

use super::fill_xyz::{multiply_mat_vec, set_coord};
use super::Plate;

impl Plate {
    /// `SolveTI1(IterationNumber)` (`cxx:370-506`).
    pub(super) fn solve_ti1(&mut self, iteration_number: i32) {
        // computation of square matrix members
        let n_dim = self.n_el + (self.order * (self.order + 1) / 2) as usize;
        self.n_dim = n_dim;
        let mut mat = MathMatrix::with_init(0, n_dim - 1, 0, n_dim - 1, 0.0);

        self.points = vec![occt_core::gp::GpXY::zero(); self.n_el];
        for i in 0..self.n_el {
            self.points[i] = *self.my_constraints[i].pnt2d();
        }

        self.deru = vec![0i32; self.n_el];
        self.derv = vec![0i32; self.n_el];
        for i in 0..self.n_el {
            self.deru[i] = self.my_constraints[i].idu();
            self.derv[i] = self.my_constraints[i].idv();
        }

        for i in 0..self.n_el {
            for j in 0..i {
                let mut signe = 1.0;
                if (self.deru(j) + self.derv(j)) % 2 == 1 {
                    signe = -1.0;
                }
                let iu = self.deru(i) + self.deru(j);
                let iv = self.derv(i) + self.derv(j);
                mat.set_value(
                    i,
                    j,
                    signe * self.sol_em(self.point_at(i).subtracted(&self.point_at(j)), iu, iv),
                );
            }
        }

        let mut i = self.n_el;
        for iu in 0..self.order {
            let mut iv = 0;
            while iu + iv < self.order {
                for j in 0..self.n_el {
                    let idu = self.deru(j);
                    let idv = self.derv(j);
                    mat.set_value(i, j, self.polm(self.point_at(j), iu, iv, idu, idv));
                }
                i += 1;
                iv += 1;
            }
        }

        for i in 0..n_dim {
            for j in (i + 1)..n_dim {
                let v = mat.value(j, i);
                mat.set_value(i, j, v);
            }
        }

        // initialisation of the Gauss algorithm
        let mut pivot_max = 1.0e-12;
        self.ok = true;
        let mut algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);

        if !algo_gauss.is_done() {
            let nbm = (self.order * (self.order + 1) / 2) as usize;
            for i in self.n_el..self.n_el + nbm {
                mat.set_value(i, i, 1.0e-8);
            }
            pivot_max = 1.0e-18;
            algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);
            self.ok = algo_gauss.is_done();
        }

        if self.ok {
            // computation of the linear system solution for the X, Y and Z
            // coordinates
            let mut sec_member = MathVector::with_init(0, n_dim - 1, 0.0);
            let mut sol = MathVector::new(0, n_dim - 1);

            let mut new_solution = vec![GpXyz::zero(); n_dim];

            for icoor in 1..=3usize {
                for i in 0..self.n_el {
                    sec_member.set_value(i, self.my_constraints[i].value().coord(icoor - 1));
                }
                if algo_gauss.solve(&sec_member, &mut sol).is_err() {
                    self.ok = false;
                    return;
                }
                // iteration to refine the solution
                {
                    let mut sol1 = MathVector::new(0, n_dim - 1);
                    for _ in 1..=iteration_number {
                        let prod = multiply_mat_vec(&mat, &sol);
                        let sec_member1 = sec_member.subtracted(&prod);
                        if algo_gauss.solve(&sec_member1, &mut sol1).is_err() {
                            self.ok = false;
                            return;
                        }
                        sol.add(&sol1);
                    }
                }

                for i in 0..n_dim {
                    set_coord(&mut new_solution[i], icoor, sol.value(i));
                }
            }
            self.solution = new_solution;
        }
    }
}
