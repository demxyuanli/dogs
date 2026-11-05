//! `Plate_Plate::SolveTI3` — the general case: `PinpointConstraint`,
//! `LinearXYZConstraint` and `LinearScalarConstraint` all present.
//!
//! Source: `Plate_Plate.cxx:669-1056`.
//!
//! The system is `3 * n_dimsousmat + n_cc3` wide (`cxx:690-691`): three copies
//! of the `(n_cc1 + n_cc2 + nbm)` block for X / Y / Z, plus one row per
//! `LinearScalarConstraint` row that couples the three coordinates through the
//! `gp_XYZ` coefficients' dot products.
//!
//! Deviations (no control-flow change): `Message_ProgressScope`
//! (`cxx:916-944`) is UNPORTED, and a failed `algo_gauss.Solve` sets
//! `ok = false` instead of raising `StdFail_NotDone`.

use occt_core::gp::{GpXY, GpXyz};
use occt_math::{MathGauss, MathMatrix, MathVector};

use super::fill_xyz::{multiply_mat_vec, set_coord};
use super::Plate;

/// `mat(i, j) += value`.
fn add(mat: &mut MathMatrix, i: usize, j: usize, value: f64) {
    let old = mat.value(i, j);
    mat.set_value(i, j, old + value);
}

impl Plate {
    /// `SolveTI3(IterationNumber)` (`cxx:669-1056`).
    pub(super) fn solve_ti3(&mut self, iteration_number: i32) {
        // computation of square matrix members
        let n_cc1 = self.my_constraints.len();

        let mut n_cc2 = 0usize;
        for c in &self.my_lxyz_constraints {
            n_cc2 += c.col_length();
        }

        let mut n_cc3 = 0usize;
        for c in &self.my_lscalar_constraints {
            n_cc3 += c.col_length();
        }

        let nbm = (self.order * (self.order + 1) / 2) as usize;
        let n_dimsousmat = n_cc1 + n_cc2 + nbm;
        let n_dimat = 3 * n_dimsousmat + n_cc3;

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
        let n_ppc2 = k;
        for c in &self.my_lscalar_constraints {
            for ppc in c.ppc() {
                self.points[k] = *ppc.pnt2d();
                self.deru[k] = ppc.idu();
                self.derv[k] = ppc.idv();
                k += 1;
            }
        }

        let mut mat = MathMatrix::with_init(0, n_dimat - 1, 0, n_dimat - 1, 0.0);

        self.fill_xyz_matrix(&mut mat, 0, 0, n_cc1, n_cc2);
        self.fill_xyz_matrix(&mut mat, n_dimsousmat, n_dimsousmat, n_cc1, n_cc2);
        self.fill_xyz_matrix(&mut mat, 2 * n_dimsousmat, 2 * n_dimsousmat, n_cc1, n_cc2);

        let mut k = 3 * n_dimsousmat;
        let mut kppc = n_ppc2;
        for i in 0..self.my_lscalar_constraints.len() {
            let row_len = self.my_lscalar_constraints[i].row_length();
            let col_len = self.my_lscalar_constraints[i].col_length();

            // `cxx:741-769`: the LSC rows against every pinpoint constraint.
            for j in 0..n_cc1 {
                let mut vmat = MathVector::new(1, row_len);

                for ippc in 1..=row_len {
                    let mut signe = 1.0;
                    if (self.deru(j) + self.derv(j)) % 2 == 1 {
                        signe = -1.0;
                    }
                    let iu = self.deru(kppc + ippc - 1) + self.deru(j);
                    let iv = self.derv(kppc + ippc - 1) + self.derv(j);
                    vmat.set_value(
                        ippc,
                        signe
                            * self.sol_em(
                                self.point_at(kppc + ippc - 1).subtracted(&self.point_at(j)),
                                iu,
                                iv,
                            ),
                    );
                }

                for irow in 1..=col_len {
                    for icol in 1..=row_len {
                        let cf = self.my_lscalar_constraints[i].coeff(irow, icol);
                        let v = vmat.value(icol);
                        add(&mut mat, k + irow - 1, j, cf.x * v);
                        add(&mut mat, k + irow - 1, n_dimsousmat + j, cf.y * v);
                        add(&mut mat, k + irow - 1, 2 * n_dimsousmat + j, cf.z * v);
                    }
                }
            }

            // `cxx:771-822`: the LSC rows against every LinearXYZ constraint.
            let mut k2 = n_cc1;
            let mut kppc2 = n_cc1;
            for i2 in 0..self.my_lxyz_constraints.len() {
                let row_len2 = self.my_lxyz_constraints[i2].row_length();
                let col_len2 = self.my_lxyz_constraints[i2].col_length();

                let mut tmpmat = MathMatrix::new(1, row_len, 1, row_len2);

                for ippc in 1..=row_len {
                    for ippc2 in 1..=row_len2 {
                        let mut signe = 1.0;
                        if (self.deru(kppc2 + ippc2 - 1) + self.derv(kppc2 + ippc2 - 1)) % 2 == 1
                        {
                            signe = -1.0;
                        }
                        let iu = self.deru(kppc + ippc - 1) + self.deru(kppc2 + ippc2 - 1);
                        let iv = self.derv(kppc + ippc - 1) + self.derv(kppc2 + ippc2 - 1);
                        tmpmat.set_value(
                            ippc,
                            ippc2,
                            signe
                                * self.sol_em(
                                    self.point_at(kppc + ippc - 1)
                                        .subtracted(&self.point_at(kppc2 + ippc2 - 1)),
                                    iu,
                                    iv,
                                ),
                        );
                    }
                }

                for irow in 1..=col_len {
                    for irow2 in 1..=col_len2 {
                        for icol in 1..=row_len {
                            for icol2 in 1..=row_len2 {
                                let c1 = self.my_lscalar_constraints[i].coeff(irow, icol);
                                let c2 = self.my_lxyz_constraints[i2].coeff(irow2, icol2);
                                let t = tmpmat.value(icol, icol2);
                                add(&mut mat, k + irow - 1, k2 + irow2 - 1, c1.x * c2 * t);
                                add(
                                    &mut mat,
                                    k + irow - 1,
                                    n_dimsousmat + k2 + irow2 - 1,
                                    c1.y * c2 * t,
                                );
                                add(
                                    &mut mat,
                                    k + irow - 1,
                                    2 * n_dimsousmat + k2 + irow2 - 1,
                                    c1.z * c2 * t,
                                );
                            }
                        }
                    }
                }

                k2 += self.my_lxyz_constraints[i2].col_length();
                kppc2 += self.my_lxyz_constraints[i2].row_length();
            }

            // `cxx:824-852`: the LSC rows against the polynomial block.
            let mut j = n_cc1 + n_cc2;
            for iu in 0..self.order {
                let mut iv = 0;
                while iu + iv < self.order {
                    let mut vmat = MathVector::new(1, row_len);
                    for ippc in 1..=row_len {
                        let idu = self.deru(kppc + ippc - 1);
                        let idv = self.derv(kppc + ippc - 1);
                        vmat.set_value(
                            ippc,
                            self.polm(self.point_at(kppc + ippc - 1), iu, iv, idu, idv),
                        );
                    }

                    for irow in 1..=col_len {
                        for icol in 1..=row_len {
                            let cf = self.my_lscalar_constraints[i].coeff(irow, icol);
                            let v = vmat.value(icol);
                            add(&mut mat, k + irow - 1, j, cf.x * v);
                            add(&mut mat, k + irow - 1, n_dimsousmat + j, cf.y * v);
                            add(&mut mat, k + irow - 1, 2 * n_dimsousmat + j, cf.z * v);
                        }
                    }

                    j += 1;
                    iv += 1;
                }
            }

            // `cxx:854-898`: the LSC rows against every LSC constraint up to `i`.
            let mut k2 = 3 * n_dimsousmat;
            let mut kppc2 = n_ppc2;
            for i2 in 0..=i {
                let row_len2 = self.my_lscalar_constraints[i2].row_length();
                let col_len2 = self.my_lscalar_constraints[i2].col_length();

                let mut tmpmat = MathMatrix::new(1, row_len, 1, row_len2);

                for ippc in 1..=row_len {
                    for ippc2 in 1..=row_len2 {
                        let mut signe = 1.0;
                        if (self.deru(kppc2 + ippc2 - 1) + self.derv(kppc2 + ippc2 - 1)) % 2 == 1
                        {
                            signe = -1.0;
                        }
                        let a_iu = self.deru(kppc + ippc - 1) + self.deru(kppc2 + ippc2 - 1);
                        let iv = self.derv(kppc + ippc - 1) + self.derv(kppc2 + ippc2 - 1);
                        tmpmat.set_value(
                            ippc,
                            ippc2,
                            signe
                                * self.sol_em(
                                    self.point_at(kppc + ippc - 1)
                                        .subtracted(&self.point_at(kppc2 + ippc2 - 1)),
                                    a_iu,
                                    iv,
                                ),
                        );
                    }
                }

                for irow in 1..=col_len {
                    for irow2 in 1..=col_len2 {
                        for icol in 1..=row_len {
                            for icol2 in 1..=row_len2 {
                                let c1 = self.my_lscalar_constraints[i].coeff(irow, icol);
                                let c2 = self.my_lscalar_constraints[i2].coeff(irow2, icol2);
                                let t = tmpmat.value(icol, icol2);
                                add(&mut mat, k + irow - 1, k2 + irow2 - 1, c1.dot(&c2) * t);
                            }
                        }
                    }
                }

                k2 += self.my_lscalar_constraints[i2].col_length();
                kppc2 += self.my_lscalar_constraints[i2].row_length();
            }

            k += self.my_lscalar_constraints[i].col_length();
            kppc += self.my_lscalar_constraints[i].row_length();
        }

        for j in (3 * n_dimsousmat)..n_dimat {
            for i in 0..j {
                let v = mat.value(j, i);
                mat.set_value(i, j, v);
            }
        }

        // initialisation of the Gauss algorithm
        let mut pivot_max = 1.0e-12;
        self.ok = true;
        let mut algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);

        if !algo_gauss.is_done() {
            for i in (n_cc1 + n_cc2)..(n_cc1 + n_cc2 + nbm) {
                mat.set_value(i, i, 1.0e-8);
                mat.set_value(n_dimsousmat + i, n_dimsousmat + i, 1.0e-8);
                mat.set_value(2 * n_dimsousmat + i, 2 * n_dimsousmat + i, 1.0e-8);
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

            let n_dim = self.n_el + nbm;
            self.n_dim = n_dim;
            let mut new_solution = vec![GpXyz::zero(); n_dim];

            for icoor in 1..=3usize {
                for i in 0..n_cc1 {
                    sec_member
                        .set_value((icoor - 1) * n_dimsousmat + i, self.my_constraints[i].value().coord(icoor - 1));
                }

                let mut k = n_cc1;
                for c in &self.my_lxyz_constraints {
                    let row_len = c.row_length();
                    let col_len = c.col_length();
                    for irow in 1..=col_len {
                        for icol in 1..=row_len {
                            let old = sec_member.value((icoor - 1) * n_dimsousmat + k);
                            sec_member.set_value(
                                (icoor - 1) * n_dimsousmat + k,
                                old
                                    + c.coeff(irow, icol)
                                        * c.ppc()[icol - 1].value().coord(icoor - 1),
                            );
                        }
                        k += 1;
                    }
                }
            }
            let mut k = 3 * n_dimsousmat;
            for c in &self.my_lscalar_constraints {
                let row_len = c.row_length();
                let col_len = c.col_length();
                for irow in 1..=col_len {
                    for icol in 1..=row_len {
                        let old = sec_member.value(k);
                        sec_member.set_value(
                            k,
                            old + c.coeff(irow, icol).dot(c.ppc()[icol - 1].value()),
                        );
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

            for icoor in 1..=3usize {
                for i in 0..n_cc1 {
                    set_coord(
                        &mut new_solution[i],
                        icoor,
                        sol.value((icoor - 1) * n_dimsousmat + i),
                    );
                }

                let mut k_solution = n_cc1;
                let mut ksol = n_cc1;

                for c in &self.my_lxyz_constraints {
                    let row_len = c.row_length();
                    let col_len = c.col_length();
                    for icol in 1..=row_len {
                        let mut vsol = 0.0;
                        for irow in 1..=col_len {
                            vsol += c.coeff(irow, icol)
                                * sol.value((icoor - 1) * n_dimsousmat + ksol + irow - 1);
                        }
                        set_coord(&mut new_solution[k_solution], icoor, vsol);
                        k_solution += 1;
                    }
                    ksol += col_len;
                }

                ksol = n_cc1 + n_cc2;
                for i in 0..nbm {
                    set_coord(
                        &mut new_solution[self.n_el + i],
                        icoor,
                        sol.value((icoor - 1) * n_dimsousmat + ksol + i),
                    );
                }
            }

            let mut ksol = 3 * n_dimsousmat;
            let mut k_solution = n_ppc2;
            for c in &self.my_lscalar_constraints {
                let row_len = c.row_length();
                let col_len = c.col_length();
                for icol in 1..=row_len {
                    let mut vsol = GpXyz::zero();
                    for irow in 1..=col_len {
                        vsol = vsol.added(&c.coeff(irow, icol).multiplied(sol.value(ksol + irow - 1)));
                    }
                    new_solution[k_solution] = vsol;
                    k_solution += 1;
                }
                ksol += col_len;
            }

            self.solution = new_solution;
        }
    }
}
