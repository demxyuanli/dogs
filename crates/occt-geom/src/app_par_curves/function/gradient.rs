//! `AppParCurves_Function::Perform` (`AppParCurves_Function.gxx:295-646`), the
//! shared body of `Gradient` and `Values`.
//!
//! Translation notes:
//! - `K` is bound as a reference before `Resol.ConstraintDerivative` runs
//!   (`gxx:382-383`), and `ConstraintDerivative` writes its "same tangency
//!   ratio" rows into the very same `Cont` matrix. OCCT therefore reads the
//!   **mutated** `K` from `gxx:385` onwards; the port clones `K` after the
//!   derivative call to keep that order.
//! - `TK` (`gxx:385-386`) is assigned and never read again; the port keeps the
//!   computation as `_tk` to preserve the control flow and the dimension
//!   checks it performs.
//! - `DPTCO(...) = DPTCO(...) - RESTM(j, k) * DPTCO1(Pi, j + Inc)`
//!   (`gxx:595-601`) indexes `DPTCO1` with `j` where the accumulation
//!   variable is `k`; `RESTM(j, k)` is still multiplied. Kept verbatim.
//! - `Scal = ...` (assignment, not `+=`) inside the last `j` loop
//!   (`gxx:629-635`) keeps only the last index' contribution. Kept verbatim.
//! - `Cons` (`gxx:326`) is only reassigned when a constraint couple matches
//!   `Pi`; a point with no couple keeps the previous point's value. Kept.

use occt_math::{MathMatrix, MathVector};

use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::resol_constraint::ResolConstraint;

use super::matrix_ops::{col, dot, inverse, mat_mul, mat_vec, row, transposed};
use super::{at_least_pass, ParFunction};
use crate::approx::line_tool::LineTool;

/// `Perform(X)` (`gxx:295-646`).
///
/// The `unused_assignments` cases below are OCCT's own: `inc` is reset before
/// its pre-loop value is read (`gxx:568`) and `Erruzaz` is refreshed for every
/// 2d curve while only the 3d arm consumes it (`gxx:614-625`).
#[allow(unused_assignments)]
pub(super) fn perform<L: LineTool>(this: &mut ParFunction<L>, x: &MathVector) {
    this.my_parameters = x.clone();

    // Resolution moindres carres (`gxx:301-309`).
    this.my_least_square.perform_parameters(&this.my_parameters);
    if !this.my_least_square.is_done() {
        this.done = false;
        return;
    }

    for j in this.my_parameters.lower()..=this.my_parameters.upper() {
        this.val_grad_f.set_value(j, 0.0);
    }

    if !this.contraintes {
        // `MyLeastSquare.ErrorGradient(ValGrad_F, FVal, ERR3d, ERR2d)`
        // (`gxx:316-319`).
        let (f, err3d, err2d) = this.my_least_square.error_gradient(&mut this.val_grad_f);
        this.f_val = f;
        this.err3d = err3d;
        this.err2d = err2d;
        return;
    }

    // Resolution avec contraintes (`gxx:320-645`).
    let npol = this.degre + 1;
    this.err3d = 0.0;
    this.err2d = 0.0;

    let ncu = this.nb_cu.max(0) as usize;
    let mut ptcox = MathMatrix::new(1, npol.max(1) as usize, 1, ncu.max(1));
    let mut ptcoy = MathMatrix::new(1, npol.max(1) as usize, 1, ncu.max(1));
    let mut ptcoz = MathMatrix::new(1, npol.max(1) as usize, 1, ncu.max(1));
    let mut ptcx = MathMatrix::new(1, npol.max(1) as usize, 1, ncu.max(1));
    let mut ptcy = MathMatrix::new(1, npol.max(1) as usize, 1, ncu.max(1));
    let mut ptcz = MathMatrix::new(1, npol.max(1) as usize, 1, ncu.max(1));

    this.my_multi_curve = this.my_least_square.bezier_value();

    // `gxx:339-355`.
    for ci in 1..=this.nb_cu {
        let dimen = this.tabdim[(ci - 1) as usize];
        for j in 1..=npol {
            let mp = this.my_multi_curve.value(j);
            if dimen == 3 {
                let p = mp.point(ci);
                ptcox.set_value(j as usize, ci as usize, p.x());
                ptcoy.set_value(j as usize, ci as usize, p.y());
                ptcoz.set_value(j as usize, ci as usize, p.z());
            } else {
                let p = mp.point2d(ci);
                ptcox.set_value(j as usize, ci as usize, p.x());
                ptcoy.set_value(j as usize, ci as usize, p.y());
                ptcoz.set_value(j as usize, ci as usize, 0.0);
            }
        }
    }

    this.a = this.my_least_square.function_matrix().clone();
    this.da = this.my_least_square.derivative_function_matrix().clone();

    // Resolution avec contraintes (`gxx:358-367`).
    let mut resol = ResolConstraint::new(
        &this.my_multi_line,
        &mut this.my_multi_curve,
        this.first_p,
        this.last_p,
        &this.my_constraints,
        &this.a,
        &this.da,
        1.0e-10,
    );
    if !resol.is_done() {
        this.done = false;
        return;
    }

    // `gxx:371-380`.
    this.f_val = 0.0;
    for j in this.first_p..=this.last_p {
        this.val_grad_f.set_value(j as usize, 0.0);
    }

    // `const math_Matrix& K = Resol.ConstraintMatrix();` (`gxx:382`) is read
    // only after `ConstraintDerivative` (`gxx:383`) has rewritten it, so the
    // clone below comes after the call; `DK` is the returned derivative.
    let dk_owned = resol
        .constraint_derivative(&this.my_multi_line, x, this.degre, &this.da)
        .clone();
    let k = resol.constraint_matrix().clone();
    let vardua = resol.duale().clone();
    let kk = mat_mul(&transposed(&k), resol.inverse_matrix());
    let dtk = transposed(&dk_owned);
    let tra = transposed(&this.a);
    let trda = transposed(&this.da);
    // `RESTM = ((A.Transposed() * A).Inverse())` (`gxx:392`).
    let restm = inverse(&mat_mul(&transposed(&this.a), &this.a));

    // `gxx:395-398`.
    let mut dptco = MathVector::new(1, k.col_count());
    let mut dptco1 = MathMatrix::with_init(
        this.first_p.max(0) as usize,
        this.last_p.max(0) as usize,
        1,
        k.col_count(),
        0.0,
    );
    let mut dkptc = MathVector::new(1, k.row_count());

    // `gxx:399-415`: the moved poles, curve by curve.
    this.f_val = 0.0;
    for ci in 1..=this.nb_cu {
        let dimen = this.tabdim[(ci - 1) as usize];
        for j in 1..=npol {
            let mp = this.my_multi_curve.value(j);
            if dimen == 3 {
                let p = mp.point(ci);
                ptcx.set_value(j as usize, ci as usize, p.x());
                ptcy.set_value(j as usize, ci as usize, p.y());
                ptcz.set_value(j as usize, ci as usize, p.z());
            } else {
                let p = mp.point2d(ci);
                ptcx.set_value(j as usize, ci as usize, p.x());
                ptcy.set_value(j as usize, ci as usize, p.y());
                ptcz.set_value(j as usize, ci as usize, 0.0);
            }
        }
    }

    // Calcul du gradient sans contraintes (`gxx:417-471`).
    for ci in 1..=this.nb_cu {
        let dimen = this.tabdim[(ci - 1) as usize];
        for i in this.adeb..=this.afin {
            let mut aa = 0.0;
            let mut bb = 0.0;
            let mut cc = 0.0;
            let mut daa = 0.0;
            let mut dbb = 0.0;
            let mut dcc = 0.0;
            for j in 1..=npol {
                let aij = this.a.value(i as usize, j as usize);
                let daij = this.da.value(i as usize, j as usize);
                let px = ptcx.value(j as usize, ci as usize);
                let py = ptcy.value(j as usize, ci as usize);
                aa += aij * px;
                bb += aij * py;
                daa += daij * px;
                dbb += daij * py;
                if dimen == 3 {
                    let pz = ptcz.value(j as usize, ci as usize);
                    cc += aij * pz;
                    dcc += daij * pz;
                }
            }
            let fx = aa - this.ptlx.value(i as usize, ci as usize);
            let fy = bb - this.ptly.value(i as usize, ci as usize);
            let mut myf = fx * fx + fy * fy;
            let mut gr = 2.0 * (daa * fx + dbb * fy);
            if dimen == 3 {
                let fz = cc - this.ptlz.value(i as usize, ci as usize);
                myf += fz * fz;
                gr += 2.0 * dcc * fz;
                this.my_f.set_value(i as usize, ci as usize, myf);
                let fi = myf;
                if fi.sqrt() > this.err3d {
                    this.err3d = fi.sqrt();
                }
            } else {
                this.my_f.set_value(i as usize, ci as usize, myf);
                let fi = myf;
                if fi.sqrt() > this.err2d {
                    this.err2d = fi.sqrt();
                }
            }
            this.f_val += myf;
            let g = this.val_grad_f.value(i as usize) + gr;
            this.val_grad_f.set_value(i as usize, g);
        }
    }

    // Calcul de DK*PTC (`gxx:473-495`).
    for i in 1..=k.row_count() {
        let mut inc: i32 = 0;
        for ci in 1..=this.nb_cu {
            let dimen = this.tabdim[(ci - 1) as usize];
            let mut acc = 0.0;
            for j in 1..=npol {
                acc += dk_owned.value(i, (j + inc) as usize)
                    * ptcx.value(j as usize, ci as usize)
                    + dk_owned.value(i, (j + inc + npol) as usize)
                        * ptcy.value(j as usize, ci as usize);
                if dimen == 3 {
                    acc += dk_owned.value(i, (j + inc + 2 * npol) as usize)
                        * ptcz.value(j as usize, ci as usize);
                }
            }
            dkptc.set_value(i, acc);
            if dimen == 3 {
                inc += 3 * npol;
            } else {
                inc += 2 * npol;
            }
        }
    }

    // `DERR = (DTK)*Vardua - KK * ((DKPTC) + K * (DTK)*Vardua)`
    // (`gxx:497-498`).
    let dtk_vardua = mat_vec(&dtk, &vardua);
    let k_dtk_vardua = mat_vec(&k, &dtk_vardua);
    let mut inner = MathVector::new(dkptc.lower(), dkptc.upper());
    for n in 0..dkptc.len() {
        let idx = dkptc.lower() + n;
        inner.set_value(
            idx,
            dkptc.value(idx) + k_dtk_vardua.value(k_dtk_vardua.lower() + n),
        );
    }
    let kk_inner = mat_vec(&kk, &inner);
    let mut derr = MathVector::new(dtk_vardua.lower(), dtk_vardua.upper());
    for n in 0..dtk_vardua.len() {
        let idx = dtk_vardua.lower() + n;
        derr.set_value(idx, dtk_vardua.value(idx) - kk_inner.value(kk_inner.lower() + n));
    }

    // rajout du gradient avec contraintes (`gxx:500-563`).
    let mut inc: i32 = 0;
    let mut scalx: MathVector;
    let mut scaly: MathVector;
    let mut scalz: MathVector;
    let mut errx: MathVector;
    let mut erry: MathVector;
    let mut errz: MathVector;
    let mut erruzax: MathVector;
    let mut erruzay: MathVector;
    let mut erruzaz: MathVector;
    let mut trdapi: MathVector;
    let mut tr_api: MathVector;
    let mut ptcxci: MathVector;
    let mut ptcyci: MathVector;
    let mut ptczci: MathVector;

    for ci in 1..=this.nb_cu {
        let dimen = this.tabdim[(ci - 1) as usize];
        let ptcoxc_i = col(&ptcox, ci as usize);
        let ptcoyc_i = col(&ptcoy, ci as usize);
        let ptcozc_i = col(&ptcoz, ci as usize);
        ptcxci = col(&ptcx, ci as usize);
        ptcyci = col(&ptcy, ci as usize);
        ptczci = col(&ptcz, ci as usize);

        // `Errx = (A * PTCOXCI - PTLX.Col(Ci))` (`gxx:528-530`).
        let ptlx_col = col(&this.ptlx, ci as usize);
        let ptly_col = col(&this.ptly, ci as usize);
        let ptlz_col = col(&this.ptlz, ci as usize);
        errx = vec_sub(&mat_vec(&this.a, &ptcoxc_i), &ptlx_col);
        erry = vec_sub(&mat_vec(&this.a, &ptcoyc_i), &ptly_col);
        errz = vec_sub(&mat_vec(&this.a, &ptcozc_i), &ptlz_col);
        // `Scalx = (DA * PTCOXCI)` (`gxx:531-533`), over `DA`'s own bounds.
        scalx = mat_vec(&this.da, &ptcoxc_i);
        scaly = mat_vec(&this.da, &ptcoyc_i);
        scalz = mat_vec(&this.da, &ptcozc_i);
        erruzax = vec_sub(&ptcxci, &ptcoxc_i);
        erruzay = vec_sub(&ptcyci, &ptcoyc_i);
        erruzaz = vec_sub(&ptczci, &ptcozc_i);

        for pi in this.first_p..=this.last_p {
            trdapi = col(&trda, pi as usize);
            tr_api = col(&tra, pi as usize);
            let taa = dot(&tr_api, &row(&this.a, pi as usize));
            let mut scal = 0.0;
            for j in 1..=npol {
                let vx = scaled_add(&trdapi, errx.value(pi as usize), &tr_api, scalx.value(pi as usize));
                let vy = scaled_add(&trdapi, erry.value(pi as usize), &tr_api, scaly.value(pi as usize));
                dptco1.set_value(pi as usize, (j + inc) as usize, vx.value(j as usize));
                dptco1.set_value(pi as usize, (j + inc + npol) as usize, vy.value(j as usize));
                scal += dptco1.value(pi as usize, (j + inc) as usize)
                    * taa
                    * erruzax.value(j as usize)
                    + dptco1.value(pi as usize, (j + inc + npol) as usize)
                        * taa
                        * erruzay.value(j as usize);
                if dimen == 3 {
                    let vz = scaled_add(
                        &trdapi,
                        errz.value(pi as usize),
                        &tr_api,
                        scalz.value(pi as usize),
                    );
                    dptco1.set_value(pi as usize, (j + inc + 2 * npol) as usize, vz.value(j as usize));
                    scal += dptco1.value(pi as usize, (j + inc + 2 * npol) as usize)
                        * taa
                        * erruzaz.value(j as usize);
                }
            }
            let g = this.val_grad_f.value(pi as usize) - 2.0 * scal;
            this.val_grad_f.set_value(pi as usize, g);
        }
        if dimen == 3 {
            inc += 3 * npol;
        } else {
            inc += 2 * npol;
        }
    }

    // on calcule DPTCO = - RESTM * DPTCO1 (`gxx:565-644`).
    let mut cons = Constraint::NoConstraint;
    let mut inc: i32 = 0;
    for pi in this.first_p..=this.last_p {
        for couple in &this.my_constraints {
            if couple.index() == pi {
                cons = couple.constraint();
                break;
            }
        }
        if at_least_pass(cons) {
            inc = 0;
            for ci in 1..=this.nb_cu {
                let dimen = this.tabdim[(ci - 1) as usize];
                for j in 1..=npol {
                    let mut dx = 0.0;
                    let mut dy = 0.0;
                    let mut dz = 0.0;
                    for kk2 in 1..=npol {
                        dx -= restm.value(j as usize, kk2 as usize)
                            * dptco1.value(pi as usize, (j + inc) as usize);
                        dy -= restm.value(j as usize, kk2 as usize)
                            * dptco1.value(pi as usize, (j + inc + npol) as usize);
                        if dimen == 3 {
                            dz -= restm.value(j as usize, kk2 as usize)
                                * dptco1.value(pi as usize, (j + inc + 2 * npol) as usize);
                        }
                    }
                    dptco.set_value((j + inc) as usize, dx);
                    dptco.set_value((j + inc + npol) as usize, dy);
                    if dimen == 3 {
                        dptco.set_value((j + inc + 2 * npol) as usize, dz);
                    }
                }
                if dimen == 3 {
                    inc += 3 * npol;
                } else {
                    inc += 2 * npol;
                }
            }

            // `DERR = DERR - KK * K * DPTCO;` (`gxx:610`).
            let k_dptco = mat_vec(&k, &dptco);
            let kk_k_dptco = mat_vec(&kk, &k_dptco);
            derr = vec_sub(&derr, &kk_k_dptco);

            inc = 0;
            for ci in 1..=this.nb_cu {
                let dimen = this.tabdim[(ci - 1) as usize];
                let ptcoxc_i = col(&ptcox, ci as usize);
                let ptcoyc_i = col(&ptcoy, ci as usize);
                let ptcozc_i = col(&ptcoz, ci as usize);
                ptcxci = col(&ptcx, ci as usize);
                ptcyci = col(&ptcy, ci as usize);
                ptczci = col(&ptcz, ci as usize);
                erruzax = vec_sub(&ptcxci, &ptcoxc_i);
                erruzay = vec_sub(&ptcyci, &ptcoyc_i);
                erruzaz = vec_sub(&ptczci, &ptcozc_i);
                let mut scal = 0.0;

                for j in 1..=npol {
                    let apij = this.a.value(pi as usize, j as usize);
                    // `Scal = ...` — an assignment, not an accumulation
                    // (`gxx:629-630`).
                    scal = (apij * erruzax.value(j as usize))
                        * (apij * derr.value((j + inc) as usize))
                        + (apij * erruzay.value(j as usize))
                            * (apij * derr.value((j + inc + npol) as usize));
                    if dimen == 3 {
                        scal += (apij * erruzax.value(j as usize))
                            * (apij * derr.value((j + inc + 2 * npol) as usize));
                    }
                }

                let g = this.val_grad_f.value(pi as usize) + 2.0 * scal;
                this.val_grad_f.set_value(pi as usize, g);
                if dimen == 3 {
                    inc += 3 * npol;
                } else {
                    inc += 2 * npol;
                }
            }
        }
    }
}

/// `math_Vector` difference with OCCT's dimension check.
fn vec_sub(a: &MathVector, b: &MathVector) -> MathVector {
    if a.len() != b.len() {
        panic!("Standard_DimensionError: math_Vector::operator-");
    }
    let mut r = MathVector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set_value(i, a.value(i) - b.value(i));
    }
    r
}

/// `TrDAPI * Errx(Pi) + TrAPI * Scalx(Pi)` (`gxx:546-547`).
fn scaled_add(a: &MathVector, sa: f64, b: &MathVector, sb: f64) -> MathVector {
    if a.len() != b.len() {
        panic!("Standard_DimensionError: math_Vector::operator+");
    }
    let mut r = MathVector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set_value(i, a.value(i) * sa + b.value(i) * sb);
    }
    r
}


