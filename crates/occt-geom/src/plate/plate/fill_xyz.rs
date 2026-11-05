//! `Plate_Plate::fillXYZmatrix` — the shared assembly of the X / Y / Z block
//! of the linear system.
//!
//! Source: `Plate_Plate.cxx:1060-1216`.
//!
//! The routine fills the square block starting at `(i0, j0)`: the pinpoint part
//! (`cxx:1065-1080`), the `LinearXYZConstraint` rows and columns
//! (`cxx:1082-1161`), the polynomial rows and columns (`cxx:1163-1202`) and
//! finally mirrors the block across the diagonal (`cxx:1204-1214`).

use occt_core::gp::GpXyz;
use occt_math::{MathMatrix, MathVector};

use super::Plate;

/// `xyz.SetCoord(icoor, value)`, 1-based.
pub(super) fn set_coord(xyz: &mut GpXyz, icoor: usize, value: f64) {
    match icoor {
        1 => xyz.x = value,
        2 => xyz.y = value,
        _ => xyz.z = value,
    }
}

/// `math_Vector` times `math_Matrix` (`sec_member1 = mat * sol`,
/// `Plate_Plate.cxx:492`, `:628`, `:1000`).
pub(super) fn multiply_mat_vec(mat: &MathMatrix, vec: &MathVector) -> MathVector {
    let mut out = MathVector::new(mat.lower_row(), mat.upper_row());
    for r in mat.lower_row()..=mat.upper_row() {
        let mut sum = 0.0;
        for c in mat.lower_col()..=mat.upper_col() {
            sum += mat.value(r, c) * vec.value(c);
        }
        out.set_value(r, sum);
    }
    out
}

impl Plate {
    /// `fillXYZmatrix(mat, i0, j0, ncc1, ncc2)` (`cxx:1060-1216`).
    pub(super) fn fill_xyz_matrix(
        &self,
        mat: &mut MathMatrix,
        i0: usize,
        j0: usize,
        ncc1: usize,
        ncc2: usize,
    ) {
        for i in 0..ncc1 {
            for j in 0..i {
                let mut signe = 1.0;
                if (self.deru(j) + self.derv(j)) % 2 == 1 {
                    signe = -1.0;
                }
                let iu = self.deru(i) + self.deru(j);
                let iv = self.derv(i) + self.derv(j);
                mat.set_value(
                    i0 + i,
                    j0 + j,
                    signe * self.sol_em(self.point_at(i).subtracted(&self.point_at(j)), iu, iv),
                );
            }
        }

        let mut k = ncc1;
        let mut kppc = ncc1;
        for i in 0..self.my_lxyz_constraints.len() {
            let row_len = self.my_lxyz_constraints[i].row_length();
            let col_len = self.my_lxyz_constraints[i].col_length();

            for a_j in 0..ncc1 {
                let mut vmat = MathVector::new(1, row_len);

                for ippc in 1..=row_len {
                    let mut signe = 1.0;
                    if (self.deru(a_j) + self.derv(a_j)) % 2 == 1 {
                        signe = -1.0;
                    }
                    let iu = self.deru(kppc + ippc - 1) + self.deru(a_j);
                    let iv = self.derv(kppc + ippc - 1) + self.derv(a_j);
                    vmat.set_value(
                        ippc,
                        signe
                            * self.sol_em(
                                self.point_at(kppc + ippc - 1).subtracted(&self.point_at(a_j)),
                                iu,
                                iv,
                            ),
                    );
                }

                for irow in 1..=col_len {
                    for icol in 1..=row_len {
                        let old = mat.value(i0 + k + irow - 1, j0 + a_j);
                        mat.set_value(
                            i0 + k + irow - 1,
                            j0 + a_j,
                            old + self.my_lxyz_constraints[i].coeff(irow, icol) * vmat.value(icol),
                        );
                    }
                }
            }

            let mut k2 = ncc1;
            let mut kppc2 = ncc1;
            for i2 in 0..=i {
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
                                let old =
                                    mat.value(i0 + k + irow - 1, j0 + k2 + irow2 - 1);
                                mat.set_value(
                                    i0 + k + irow - 1,
                                    j0 + k2 + irow2 - 1,
                                    old
                                        + self.my_lxyz_constraints[i].coeff(irow, icol)
                                            * self.my_lxyz_constraints[i2].coeff(irow2, icol2)
                                            * tmpmat.value(icol, icol2),
                                );
                            }
                        }
                    }
                }

                k2 += self.my_lxyz_constraints[i2].col_length();
                kppc2 += self.my_lxyz_constraints[i2].row_length();
            }

            k += self.my_lxyz_constraints[i].col_length();
            kppc += self.my_lxyz_constraints[i].row_length();
        }

        let mut i = ncc1 + ncc2;
        for iu in 0..self.order {
            let mut iv = 0;
            while iu + iv < self.order {
                for a_j in 0..ncc1 {
                    let idu = self.deru(a_j);
                    let idv = self.derv(a_j);
                    mat.set_value(
                        i0 + i,
                        j0 + a_j,
                        self.polm(self.point_at(a_j), iu, iv, idu, idv),
                    );
                }

                let mut k2 = ncc1;
                let mut kppc2 = ncc1;
                for i2 in 0..self.my_lxyz_constraints.len() {
                    let row_len2 = self.my_lxyz_constraints[i2].row_length();
                    let col_len2 = self.my_lxyz_constraints[i2].col_length();

                    let mut vmat = MathVector::new(1, row_len2);
                    for ippc2 in 1..=row_len2 {
                        let idu = self.deru(kppc2 + ippc2 - 1);
                        let idv = self.derv(kppc2 + ippc2 - 1);
                        vmat.set_value(
                            ippc2,
                            self.polm(self.point_at(kppc2 + ippc2 - 1), iu, iv, idu, idv),
                        );
                    }

                    for irow2 in 1..=col_len2 {
                        for icol2 in 1..=row_len2 {
                            let old = mat.value(i0 + i, j0 + k2 + irow2 - 1);
                            mat.set_value(
                                i0 + i,
                                j0 + k2 + irow2 - 1,
                                old
                                    + self.my_lxyz_constraints[i2].coeff(irow2, icol2)
                                        * vmat.value(icol2),
                            );
                        }
                    }

                    k2 += self.my_lxyz_constraints[i2].col_length();
                    kppc2 += self.my_lxyz_constraints[i2].row_length();
                }

                i += 1;
                iv += 1;
            }
        }

        let n_dimat = ncc1 + ncc2 + (self.order * (self.order + 1) / 2) as usize;

        for i in 0..n_dimat {
            for a_j in (i + 1)..n_dimat {
                let v = mat.value(i0 + a_j, j0 + i);
                mat.set_value(i0 + i, j0 + a_j, v);
            }
        }
    }
}
