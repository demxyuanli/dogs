//! The constructor body of `AppParCurves_ResolConstraint`
//! (`AppParCurves_ResolConstraint.gxx:118-503`): sizing of the local
//! Bernstein blocks, the tangency rows, the right-hand side, the starting
//! vector, the Uzawa call and the pole write-back.
//!
//! Translation notes:
//! - `ToolLine::Tangency(SSP, Npt, tabV)` only writes `tabV` when the point
//!   carries tangency data (`AppDef_MyLineTool.cxx:88-106`); the helpers below
//!   keep that "leave the previous contents" behaviour instead of clearing the
//!   array first.
//! - `ToolLine::Value` always fills the array, so the point loops read the
//!   returned `Vec` directly.
//! - `CTCinv(i, j) = CTCinv(j, i)` for `i <= j` (`gxx:462-469`) fills the
//!   upper triangle from the lower one; `math_Uzawa::InverseCont` only
//!   returns the lower triangle (see `crates/occt-math/src/uzawa.rs`).

use occt_core::gp::{GpVec, GpVec2d};
use occt_math::{MathMatrix, MathVector};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::app_par_curves::resol_constraint::{set_block, set_ibont, ResolConstraint};

use super::ratios;

/// `ToolLine::Tangency(SSP, Npt, tabV)` (`AppDef_MyLineTool.cxx:88-106`):
/// fills only when the point is a tangency point.
pub(super) fn tangency_3d<L: LineTool>(ssp: &L, npt: i32, tab_v: &mut [GpVec]) {
    if let Some(v) = ssp.tangency_3d(npt) {
        for (slot, val) in tab_v.iter_mut().zip(v.iter()) {
            *slot = *val;
        }
    }
}

/// `ToolLine::Tangency(SSP, Npt, tabV2d)` (`AppDef_MyLineTool.cxx:108-127`).
pub(super) fn tangency_2d<L: LineTool>(ssp: &L, npt: i32, tab_v2d: &mut [GpVec2d]) {
    if let Some(v) = ssp.tangency_2d(npt) {
        for (slot, val) in tab_v2d.iter_mut().zip(v.iter()) {
            *slot = *val;
        }
    }
}

/// `gxx:118-503`, run from [`ResolConstraint::new`] once the counts are known.
pub(super) fn perform<L: LineTool>(
    this: &mut ResolConstraint,
    ssp: &L,
    scurv: &mut MultiCurve,
    bern: &MathMatrix,
    derivative_bern: &MathMatrix,
    tolerance: f64,
    fc: Constraint,
    lc: Constraint,
) {
    let inc_pass = this.inc_pass;
    let inc_tan = this.inc_tan;
    let nb3d = ssp.nb_p3d();
    let nb2d = ssp.nb_p2d();
    let nb_cu = scurv.nb_curves() as i32;
    let def = scurv.nb_poles() as i32 - 1;
    let npol = def + 1;
    let npol2 = 2 * npol;
    let ccol = nb3d * 3 + nb2d * 2;
    let mynb3d = if nb3d == 0 { 1 } else { nb3d };
    let mynb2d = if nb2d == 0 { 1 } else { nb2d };

    // `math_Matrix ContInit(1, IncPass, 1, Npol); math_Vector Start(1,
    // CCol * Npol); NCollection_Array2<int> Ibont(1, NbCu, 1, IncTan);`
    // (`gxx:129-131`).
    let mut cont_init = MathMatrix::with_init(1, inc_pass, 1, npol as usize, 0.0);
    let mut start = MathVector::new(1, (ccol * npol) as usize);
    let mut ibonts = vec![0i32; (nb_cu * inc_tan as i32).max(0) as usize];

    // Filling Cont for the passing points (`gxx:135-146`).
    for i in 1..=inc_pass as i32 {
        let npt = this.ipas[(i - 1) as usize];
        for j in 1..=npol {
            cont_init.set_value(
                i as usize,
                j as usize,
                bern.value(npt as usize, j as usize),
            );
        }
    }
    for i in 1..=ccol {
        set_block(
            &mut this.cont,
            inc_pass as i32 * (i - 1) + 1,
            inc_pass as i32 * i,
            npol * (i - 1) + 1,
            npol * i,
            &cont_init,
        );
    }

    // `NCollection_Array1` working buffers (`gxx:152-155`).
    let mut tab_v = vec![GpVec::zero(); mynb3d.max(0) as usize];
    let mut tab_v2d = vec![GpVec2d::zero(); mynb2d.max(0) as usize];

    let mut inc3: i32 = ccol * inc_pass as i32 + 1;
    let mut inc_col: i32 = 0;
    let mut inc_sec: i32 = 0;

    // Tangency rows, right-hand side and starting vector, curve by curve
    // (`gxx:160-286`).
    for k in 1..=nb_cu {
        if k <= nb3d {
            for i in 1..=inc_tan as i32 {
                let npt = this.itan[(i - 1) as usize];
                // Choix du maximum de tangence pour exprimer la colinearite
                // (`gxx:168-182`).
                tangency_3d(ssp, npt, &mut tab_v);
                let v = tab_v[(k - 1) as usize];
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
                if bound == 3 {
                    for j in 1..=npol {
                        let daij = derivative_bern.value(npt as usize, j as usize);
                        this.cont.set_value(
                            inc3 as usize,
                            (j + npol + inc_col) as usize,
                            daij * t3 / tmax,
                        );
                        this.cont.set_value(
                            inc3 as usize,
                            (j + npol2 + inc_col) as usize,
                            -daij * t2 / tmax,
                        );
                        this.cont.set_value(
                            (inc3 + 1) as usize,
                            (j + inc_col) as usize,
                            daij * t3 / tmax,
                        );
                        this.cont.set_value(
                            (inc3 + 1) as usize,
                            (j + npol2 + inc_col) as usize,
                            -daij * t1 / tmax,
                        );
                    }
                } else if bound == 1 {
                    for j in 1..=npol {
                        let daij = derivative_bern.value(npt as usize, j as usize);
                        this.cont.set_value(
                            inc3 as usize,
                            (j + inc_col) as usize,
                            daij * t3 / tmax,
                        );
                        this.cont.set_value(
                            inc3 as usize,
                            (j + npol2 + inc_col) as usize,
                            -daij * t1 / tmax,
                        );
                        this.cont.set_value(
                            (inc3 + 1) as usize,
                            (j + inc_col) as usize,
                            daij * t2 / tmax,
                        );
                        this.cont.set_value(
                            (inc3 + 1) as usize,
                            (j + npol + inc_col) as usize,
                            -daij * t1 / tmax,
                        );
                    }
                } else if bound == 2 {
                    // `for (j = 1; j <= Def + 1; j++)` (`gxx:208`), i.e. Npol.
                    for j in 1..=def + 1 {
                        let daij = derivative_bern.value(npt as usize, j as usize);
                        this.cont.set_value(
                            inc3 as usize,
                            (j + npol + inc_col) as usize,
                            daij * t3 / tmax,
                        );
                        this.cont.set_value(
                            inc3 as usize,
                            (j + npol2 + inc_col) as usize,
                            -daij * t2 / tmax,
                        );
                        this.cont.set_value(
                            (inc3 + 1) as usize,
                            (j + inc_col) as usize,
                            daij * t2 / tmax,
                        );
                        this.cont.set_value(
                            (inc3 + 1) as usize,
                            (j + npol + inc_col) as usize,
                            -daij * t1 / tmax,
                        );
                    }
                }
                inc3 += 2;
            }

            // Remplissage du second membre (`gxx:220-228`).
            for i in 1..=inc_pass as i32 {
                let tab_p = ssp.value_3d(this.ipas[(i - 1) as usize]);
                let poi = tab_p[(k - 1) as usize];
                this.secont.set_value((i + inc_sec) as usize, poi.x());
                this.secont
                    .set_value((i + inc_pass as i32 + inc_sec) as usize, poi.y());
                this.secont
                    .set_value((i + 2 * inc_pass as i32 + inc_sec) as usize, poi.z());
            }
            inc_sec += 3 * inc_pass as i32;

            // Vecteur de depart (`gxx:231-238`).
            for j in 1..=npol {
                let poi = *scurv.value(j).point(k);
                start.set_value((j + inc_col) as usize, poi.x());
                start.set_value((j + npol + inc_col) as usize, poi.y());
                start.set_value((j + npol2 + inc_col) as usize, poi.z());
            }
            inc_col += 3 * npol;
        } else {
            for i in 1..=inc_tan as i32 {
                let npt = this.itan[(i - 1) as usize];
                tangency_2d(ssp, npt, &mut tab_v2d);
                let v2d = tab_v2d[(k - nb3d - 1) as usize];
                let (t1, t2) = (v2d.x(), v2d.y());
                let mut tmax = t1.abs();
                let mut bound = 1;
                if t2.abs() > tmax {
                    tmax = t2.abs();
                    bound = 2;
                }
                set_ibont(&mut ibonts, inc_tan, k, i, bound);
                for j in 1..=npol {
                    let daij = derivative_bern.value(npt as usize, j as usize);
                    this.cont
                        .set_value(inc3 as usize, (j + inc_col) as usize, daij * t2);
                    this.cont.set_value(
                        inc3 as usize,
                        (j + npol + inc_col) as usize,
                        -daij * t1,
                    );
                }
                inc3 += 1;
            }

            // Remplissage du second membre (`gxx:267-275`). OCCT reads
            // `tabP2d(i - nb3d)` while the loop index runs over the passing
            // points (`gxx:270`), which is one slot before the array start
            // when `i <= nb3d`; with `No_Standard_OutOfRange` (`gxx:21`) that
            // unchecked read is UB, so the port clamps it to the first slot.
            for i in 1..=inc_pass as i32 {
                let tab_p2d = ssp.value_2d(this.ipas[(i - 1) as usize]);
                let poi2d = tab_p2d[(i - nb3d - 1).max(0) as usize];
                this.secont.set_value((i + inc_sec) as usize, poi2d.x());
                this.secont
                    .set_value((i + inc_pass as i32 + inc_sec) as usize, poi2d.y());
            }
            inc_sec += 2 * inc_pass as i32;

            // Remplissage du vecteur de depart (`gxx:277-283`).
            for j in 1..=npol {
                let poi2d = *scurv.value(j).point2d(k);
                start.set_value((j + inc_col) as usize, poi2d.x());
                start.set_value((j + npol + inc_col) as usize, poi2d.y());
            }
            inc_col += npol2;
        }
    }

    // Equations exprimant le meme rapport de tangence sur chaque courbe
    // (`gxx:288-406`).
    ratios::fill_same_ratio(
        this,
        ssp,
        &mut tab_v,
        &mut tab_v2d,
        &ibonts,
        derivative_bern,
        nb_cu,
        nb3d,
        npol,
        inc3,
    );

    // Resolution par Uzawa (`gxx:453-457`). `EpsLix = Tolerance`, `EpsLic`
    // and `NbIterations` keep the `math_Uzawa.hxx:49-58` defaults.
    let uza = occt_math::Uzawa::new_with_tolerances(
        &this.cont,
        &this.secont,
        &start,
        tolerance,
        1.0e-06,
        500,
    );
    if !uza.is_done() {
        this.done = false;
        return;
    }

    // `CTCinv = UzaResol.InverseCont(); UzaResol.Duale(Vardua);` and the
    // mirroring of the lower triangle (`gxx:461-469`).
    this.ctc_inv = uza.inverse_cont().clone();
    this.vardua = uza.duale().clone();
    let nrows = this.ctc_inv.row_count();
    for i in 1..=nrows {
        for j in i..=nrows {
            let v = this.ctc_inv.value(j, i);
            this.ctc_inv.set_value(i, j, v);
        }
    }
    this.done = true;

    // Pole write-back (`gxx:471-502`).
    let mut polinit = 1;
    let mut polfin = npol;
    if super::at_least_pass(fc) {
        polinit = 2;
    }
    if super::at_least_pass(lc) {
        polfin = npol - 1;
    }

    let vec_poles = uza.value().clone();
    for i in polinit..=polfin {
        let mut inc_col = 0;
        let mut mpol = crate::app_par_curves::multi_point::MultiPoint::with_size(
            nb3d.max(0) as usize,
            nb2d.max(0) as usize,
        );
        for k in 1..=nb_cu {
            if k <= nb3d {
                let pol = occt_core::gp::GpPnt::new(
                    vec_poles.value((inc_col + i) as usize),
                    vec_poles.value((inc_col + npol + i) as usize),
                    vec_poles.value((inc_col + 2 * npol + i) as usize),
                );
                mpol.set_point(k, &pol);
                inc_col += 3 * npol;
            } else {
                let pol2d = occt_core::gp::GpPnt2d::new(
                    vec_poles.value((inc_col + i) as usize),
                    vec_poles.value((inc_col + npol + i) as usize),
                );
                mpol.set_point2d(k, &pol2d);
                inc_col += 2 * npol;
            }
        }
        scurv.set_value(i, &mpol);
    }
}




