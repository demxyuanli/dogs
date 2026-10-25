//! `AppParCurves_ResolConstraint::ConstraintDerivative`
//! (`AppParCurves_ResolConstraint.gxx:573-946`): the derivative of the
//! constraint matrix with respect to the parameters, used by the gradient of
//! the approximation objective.
//!
//! Translation notes:
//! - `DeCont` is zeroed and filled with the pass-point and tangency rows
//!   (`gxx:606-716`), but the "same tangency ratio" rows write into the
//!   **member `Cont`**, not `DeCont` (`gxx:760-941`), and the function returns
//!   `DeCont`. That is OCCT's actual control flow and is kept verbatim.
//! - The 2d tangency vector of column `k` is read as `tabV2d(k)`
//!   (`gxx:697`), unlike the constructor's `tabV2d(k - nb3d)` (`gxx:247`);
//!   both spellings are kept, including the out-of-range read they imply
//!   when `k > mynb2d`.
//! - `Ibont` is filled here as well (`gxx:640-651`, `gxx:699-704`), so the
//!   ratio rows below use this call's own component choice.

use occt_core::gp::{GpVec, GpVec2d};
use occt_math::MathVector;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::bernstein::AppParCurves;
use crate::app_par_curves::resol_constraint::{ibont, set_block, set_ibont, zero_matrix, ResolConstraint};

use super::perform::{tangency_2d, tangency_3d};

/// 1-based `tabV(k)` (`Standard_OutOfRange` when OCCT would throw).
fn at3(v: &[GpVec], k: i32) -> GpVec {
    match v.get((k - 1).max(0) as usize) {
        Some(val) if k >= 1 => *val,
        _ => panic!("Standard_OutOfRange: AppParCurves_ResolConstraint::ConstraintDerivative"),
    }
}

/// 1-based `tabV2d(k)` (`Standard_OutOfRange` when OCCT would throw).
fn at2(v: &[GpVec2d], k: i32) -> GpVec2d {
    match v.get((k - 1).max(0) as usize) {
        Some(val) if k >= 1 => *val,
        _ => panic!("Standard_OutOfRange: AppParCurves_ResolConstraint::ConstraintDerivative"),
    }
}

/// `ConstraintDerivative(SSP, Parameters, Deg, DA)` (`gxx:573-946`).
pub(super) fn constraint_derivative<'a, L: LineTool>(
    this: &'a mut ResolConstraint,
    ssp: &L,
    parameters: &MathVector,
    deg: i32,
    da: &occt_math::MathMatrix,
) -> &'a occt_math::MathMatrix {
    let inc_pass = this.inc_pass;
    let inc_tan = this.inc_tan;
    let nb3d = ssp.nb_p3d();
    let nb2d = ssp.nb_p2d();
    let nb_cu = ssp.nb_p3d() + ssp.nb_p2d();
    let npol = deg + 1;
    let npol2 = 2 * npol;
    let ccol = nb3d * 3 + nb2d * 2;
    let mynb3d = if nb3d == 0 { 1 } else { nb3d };
    let mynb2d = if nb2d == 0 { 1 } else { nb2d };

    let mut ibonts = vec![0i32; (nb_cu * inc_tan as i32).max(0) as usize];
    let mut tab_v = vec![GpVec::zero(); mynb3d.max(0) as usize];
    let mut tab_v2d = vec![GpVec2d::zero(); mynb2d.max(0) as usize];
    let mut dda = MathVector::new(1, npol as usize);
    let mut de_c_init = zero_matrix(inc_pass, npol as usize);

    // `for (i...) for (j...) DeCont(i, j) = 0.0;` (`gxx:606-608`).
    for i in 1..=this.de_cont.row_count() {
        for j in 1..=this.de_cont.col_count() {
            this.de_cont.set_value(i, j, 0.0);
        }
    }

    // Remplissage de DK pour les points de passages (`gxx:610-621`).
    for i in 1..=inc_pass as i32 {
        let npt = this.ipas[(i - 1) as usize];
        for j in 1..=npol {
            de_c_init.set_value(i as usize, j as usize, da.value(npt as usize, j as usize));
        }
    }
    for i in 1..=ccol {
        set_block(
            &mut this.de_cont,
            inc_pass as i32 * (i - 1) + 1,
            inc_pass as i32 * i,
            npol * (i - 1) + 1,
            npol * i,
            &de_c_init,
        );
    }

    // Pour les points de tangence (`gxx:623-716`).
    let mut inc3: i32 = ccol * inc_pass as i32 + 1;
    let mut inc_col: i32 = 0;
    for k in 1..=nb_cu {
        if k <= nb3d {
            for i in 1..=inc_tan as i32 {
                let npt = this.itan[(i - 1) as usize];
                tangency_3d(ssp, npt, &mut tab_v);
                let v = at3(&tab_v, k);
                let (t1, t2, t3) = (v.x(), v.y(), v.z());
                let mut tmax = t1.abs();
                let mut bound = 1;
                if t2.abs() > tmax {
                    tmax = t2.abs();
                    bound = 2;
                }
                if t3.abs() > tmax {
                    tmax = t3.abs();
                    bound = 3;
                }
                set_ibont(&mut ibonts, inc_tan, k, i, bound);
                AppParCurves::second_derivative_bernstein(
                    parameters.value(npt as usize),
                    &mut dda,
                );
                if bound == 3 {
                    for j in 1..=npol {
                        let ddaij = dda.value(j as usize);
                        this.de_cont.set_value(
                            inc3 as usize,
                            (j + npol + inc_col) as usize,
                            ddaij * t3 / tmax,
                        );
                        this.de_cont.set_value(
                            inc3 as usize,
                            (j + npol2 + inc_col) as usize,
                            -ddaij * t2 / tmax,
                        );
                        this.de_cont.set_value(
                            (inc3 + 1) as usize,
                            (j + inc_col) as usize,
                            ddaij * t3 / tmax,
                        );
                        this.de_cont.set_value(
                            (inc3 + 1) as usize,
                            (j + npol2 + inc_col) as usize,
                            -ddaij * t1 / tmax,
                        );
                    }
                } else if bound == 1 {
                    for j in 1..=npol {
                        let ddaij = dda.value(j as usize);
                        this.de_cont.set_value(
                            inc3 as usize,
                            (j + inc_col) as usize,
                            ddaij * t3 / tmax,
                        );
                        this.de_cont.set_value(
                            inc3 as usize,
                            (j + npol2 + inc_col) as usize,
                            -ddaij * t1 / tmax,
                        );
                        this.de_cont.set_value(
                            (inc3 + 1) as usize,
                            (j + inc_col) as usize,
                            ddaij * t2 / tmax,
                        );
                        this.de_cont.set_value(
                            (inc3 + 1) as usize,
                            (j + npol + inc_col) as usize,
                            -ddaij * t1 / tmax,
                        );
                    }
                } else if bound == 2 {
                    for j in 1..=npol {
                        let ddaij = dda.value(j as usize);
                        this.de_cont.set_value(
                            inc3 as usize,
                            (j + npol + inc_col) as usize,
                            ddaij * t3 / tmax,
                        );
                        this.de_cont.set_value(
                            inc3 as usize,
                            (j + npol2 + inc_col) as usize,
                            -ddaij * t2 / tmax,
                        );
                        this.de_cont.set_value(
                            (inc3 + 1) as usize,
                            (j + inc_col) as usize,
                            ddaij * t2 / tmax,
                        );
                        this.de_cont.set_value(
                            (inc3 + 1) as usize,
                            (j + npol + inc_col) as usize,
                            -ddaij * t1 / tmax,
                        );
                    }
                }
                inc3 += 2;
            }
            inc_col += 3 * npol;
        } else {
            for i in 1..=inc_tan as i32 {
                let npt = this.itan[(i - 1) as usize];
                AppParCurves::second_derivative_bernstein(
                    parameters.value(npt as usize),
                    &mut dda,
                );
                tangency_2d(ssp, npt, &mut tab_v2d);
                let v2d = at2(&tab_v2d, k);
                let (t1, t2) = (v2d.x(), v2d.y());
                let mut tmax = t1.abs();
                let mut bound = 1;
                if t2.abs() > tmax {
                    tmax = t2.abs();
                    bound = 2;
                }
                set_ibont(&mut ibonts, inc_tan, k, i, bound);
                for j in 1..=npol {
                    let ddaij = dda.value(j as usize);
                    this.de_cont.set_value(
                        inc3 as usize,
                        (j + inc_col) as usize,
                        ddaij * t2,
                    );
                    this.de_cont.set_value(
                        inc3 as usize,
                        (j + npol + inc_col) as usize,
                        -ddaij * t1,
                    );
                }
                inc3 += 1;
            }
        }
    }

    // Equations exprimant le meme rapport de tangence sur chaque courbe
    // (`gxx:718-945`). These rows are written into the member `Cont`.
    let mut inc3 = inc3 - 1;
    for i in 1..=inc_tan as i32 {
        let mut inc_col = 0;
        let npt = this.itan[(i - 1) as usize];
        AppParCurves::second_derivative_bernstein(parameters.value(npt as usize), &mut dda);
        for k in 1..=nb_cu - 1 {
            inc3 += 1;
            let a_val = ibont(&ibonts, inc_tan, k, i);
            let (t1, ip, first_off) = match a_val {
                1 => {
                    if k <= nb3d {
                        tangency_3d(ssp, npt, &mut tab_v);
                        (at3(&tab_v, k).x(), 3 * npol, 0)
                    } else {
                        tangency_2d(ssp, npt, &mut tab_v2d);
                        (at2(&tab_v2d, k).x(), 2 * npol, 0)
                    }
                }
                2 => {
                    if k <= nb3d {
                        tangency_3d(ssp, npt, &mut tab_v);
                        (at3(&tab_v, k).y(), 3 * npol, npol)
                    } else {
                        tangency_2d(ssp, npt, &mut tab_v2d);
                        (at2(&tab_v2d, k).y(), 2 * npol, npol)
                    }
                }
                3 => {
                    tangency_3d(ssp, npt, &mut tab_v);
                    (at3(&tab_v, k).z(), 3 * npol, 2 * npol)
                }
                _ => panic!("AppParCurves_ResolConstraint: Ibont out of range"),
            };
            let a_next = ibont(&ibonts, inc_tan, k + 1, i);
            let (t2, add_index_2) = match a_next {
                1 => {
                    if k + 1 <= nb3d {
                        tangency_3d(ssp, npt, &mut tab_v);
                        (at3(&tab_v, k + 1).x(), 0)
                    } else {
                        tangency_2d(ssp, npt, &mut tab_v2d);
                        (at2(&tab_v2d, k + 1).x(), 0)
                    }
                }
                2 => {
                    if k + 1 <= nb3d {
                        tangency_3d(ssp, npt, &mut tab_v);
                        (at3(&tab_v, k + 1).y(), npol)
                    } else {
                        tangency_2d(ssp, npt, &mut tab_v2d);
                        (at2(&tab_v2d, k + 1).y(), npol)
                    }
                }
                3 => {
                    tangency_3d(ssp, npt, &mut tab_v);
                    (at3(&tab_v, k + 1).z(), 2 * npol)
                }
                _ => panic!("AppParCurves_ResolConstraint: Ibont out of range"),
            };
            for j in 1..=npol {
                let daij = dda.value(j as usize);
                this.cont.set_value(
                    inc3 as usize,
                    (j + inc_col + first_off) as usize,
                    daij * t2,
                );
                this.cont.set_value(
                    inc3 as usize,
                    (j + ip + inc_col + add_index_2) as usize,
                    -daij * t1,
                );
            }
            inc_col += ip;
        }
    }

    &this.de_cont
}




