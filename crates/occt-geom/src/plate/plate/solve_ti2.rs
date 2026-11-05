//! `Plate_Plate::SolveTI2` — `LinearXYZConstraint`s are loaded, but no
//! `LinearScalarConstraint`.
//!
//! Source: `Plate_Plate.cxx:512-665`.
//!
//! Deviations (no control-flow change): `Message_ProgressScope`
//! (`cxx:561-586`) is UNPORTED, and a failed `algo_gauss.Solve` sets
//! `ok = false` instead of raising `StdFail_NotDone`.

use occt_core::gp::{GpXY, GpXyz};
use occt_math::{MathGauss, MathMatrix, MathVector};

use super::fill_xyz::{multiply_mat_vec, set_coord};
use super::Plate;

impl Plate {
    /// `SolveTI2(IterationNumber)` (`cxx:512-665`).
    pub(super) fn solve_ti2(&mut self, iteration_number: i32) {
        // computation of square matrix members
        let n_cc1 = self.my_constraints.len();
        let mut n_cc2 = 0usize;
        for c in &self.my_lxyz_constraints {
            n_cc2 += c.col_length();
        }

        let n_dimat = n_cc1 + n_cc2 + (self.order * (self.order + 1) / 2) as usize;

        self.points = vec![GpXY::zero(); self.n_el];
        self.deru = vec![0i32; self.n_el];
        self.derv = vec![0i32; self.n_el];

        for i in 0..n_cc1 {
            self.points[i] = *self.my_constraints[i].pnt2d();
            self.deru[i] = self.my_constraints[i].idu();
            self.derv[i] = self.my_constraints[i].idv();
        }

        let mut k = n_cc1;
        for c in &self.my_lxyz_constraints {
            for ppc in c.ppc() {
                self.points[k] = *ppc.pnt2d();
                self.deru[k] = ppc.idu();
                self.derv[k] = ppc.idv();
                k += 1;
            }
        }

        let mut mat = MathMatrix::with_init(0, n_dimat - 1, 0, n_dimat - 1, 0.0);

        self.fill_xyz_matrix(&mut mat, 0, 0, n_cc1, n_cc2);

        // initialisation of the Gauss algorithm
        let mut pivot_max = 1.0e-12;
        self.ok = true;
        let mut algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);

        if !algo_gauss.is_done() {
            for i in (n_cc1 + n_cc2)..n_dimat {
                mat.set_value(i, i, 1.0e-8);
            }
            pivot_max = 1.0e-18;
            algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);
            self.ok = algo_gauss.is_done();
        }

        if self.ok {
            // computation of the linear system solution for the X, Y and Z
            // coordinates
            let mut sec_member = MathVector::with_init(0, n_dimat - 1, 0.0);
            let mut sol = MathVector::new(0, n_dimat - 1);

            let n_dim = self.n_el + (self.order * (self.order + 1) / 2) as usize;
            self.n_dim = n_dim;
            let mut new_solution = vec![GpXyz::zero(); n_dim];

            for icoor in 1..=3usize {
                for i in 0..n_cc1 {
                    sec_member.set_value(i, self.my_constraints[i].value().coord(icoor - 1));
                }

                let mut k = n_cc1;
                for c in &self.my_lxyz_constraints {
                    let row_len = c.row_length();
                    let col_len = c.col_length();
                    for irow in 1..=col_len {
                        for icol in 1..=row_len {
                            let v = sec_member.value(k)
                                + c.coeff(irow, icol) * c.ppc()[icol - 1].value().coord(icoor - 1);
                            sec_member.set_value(k, v);
                        }
                        k += 1;
                    }
                }

                if algo_gauss.solve(&sec_member, &mut sol).is_err() {
                    self.ok = false;
                    return;
                }
                // iteration to refine the solution
                {
                    let mut sol1 = MathVector::new(0, n_dimat - 1);
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

                for i in 0..n_cc1 {
                    set_coord(&mut new_solution[i], icoor, sol.value(i));
                }

                let mut k_solution = n_cc1;
                let mut ksol = n_cc1;

                for c in &self.my_lxyz_constraints {
                    let row_len = c.row_length();
                    let col_len = c.col_length();
                    for icol in 1..=row_len {
                        let mut vsol = 0.0;
                        for irow in 1..=col_len {
                            vsol += c.coeff(irow, icol) * sol.value(ksol + irow - 1);
                        }
                        set_coord(&mut new_solution[k_solution], icoor, vsol);
                        k_solution += 1;
                    }
                    ksol += col_len;
                }

                let nbm = (self.order * (self.order + 1) / 2) as usize;
                for i in 0..nbm {
                    set_coord(&mut new_solution[self.n_el + i], icoor, sol.value(ksol + i));
                }
            }
            self.solution = new_solution;
        }
    }
}
