//! `AppDef_Variational::TheMotor` and `Adjusting` — the optimization engine.
//!
//! Source: `AppDef_Variational.cxx:1160-1591` (`TheMotor`) and
//! `cxx:2816-2973` (`Adjusting`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `goto L8000` (`cxx:1347`) is the `break 'again` of the outer loop.
//! - `NCollection_Array1<double>` point parameters are [`MathVector`]s with
//!   lower bound 1; `myParameters` itself stays keyed by point number.
//! - The curve is shared through `Rc<RefCell<Curve>>`, which is how the
//!   criterion and this driver see one curve (`cxx:1210-1220`).

use std::cell::RefCell;
use std::rc::Rc;

use occt_math::fem_tool::{Assembly, Curve};
use occt_math::vector::MathVector;

use super::split::near_index;
use super::Variational;
use crate::app_def::linear_criteria::LinearCriteria;
use crate::app_def::smooth_criterion::SmoothCriterion;

/// `BigValue` / `SmallValue` / `SmallestValue` (`cxx:1169`).
const BIG_VALUE: f64 = 1.0e37;
const SMALL_VALUE: f64 = 1.0e-6;
const SMALLEST_VALUE: f64 = 1.0e-9;
/// `MaxNbEst` (`cxx:1177`).
const MAX_NB_EST: i32 = 5;

impl Variational {
    /// `TheMotor(J, WQuadratic, WQuality, TheCurve, Ecarts)` (`cxx:1160-1591`).
    /// The `WQuadratic` argument is unused in the OCCT body (`cxx:1162-1163`
    /// comments it out), so it is not taken here either.
    pub(super) fn the_motor(
        &mut self,
        w_quality: f64,
        the_curve: &mut Rc<RefCell<Curve>>,
        ecarts: &mut MathVector,
    ) {
        let nbr_pnt = self.last_point - self.first_point + 1;
        let nbr_constraint = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;

        // (0) Init (`cxx:1192-1220`).
        let (mut e1, mut e2, mut e3) = self.smooth_criterion.estimation();
        let j1min = 1.0e-8;
        let j2min = (e1 + 1.0e-8) * 1.0e-6;
        let j3min = j2min;
        if e1 < j1min {
            e1 = j1min;
        }
        if e2 < j2min {
            e2 = j2min;
        }
        if e3 < j3min {
            e3 = j3min;
        }
        self.smooth_criterion.set_estimation(e1, e2, e3);

        let mut ccurrent = the_curve.clone();
        let mut current_ti = MathVector::from_slice(self.parameters.as_slice());
        let mut old_ti = current_ti.clone();
        let mut cold = ccurrent.clone();
        let mut cblong = *self.smooth_criterion.est_length();
        let mut lnold = cblong;
        let dependence = self.smooth_criterion.dependence_table();

        self.smooth_criterion.set_curve(&ccurrent);
        let mut the_assembly = Assembly::new(&dependence, &self.smooth_criterion.assembly_table());

        let mut vocri = [BIG_VALUE; 3];
        let mut erold = BIG_VALUE;
        let mut valcri = [0.0f64; 3];
        let mut errmax = BIG_VALUE;
        let mut errmoy = 0.0f64;
        let mut errqua = 0.0f64;
        let mut isnear = false;
        let mut lconst;
        let mut eps_length = SMALL_VALUE;
        let mut eps_deg;

        // ============ Optimization (`cxx:1223`).
        'again: loop {
            // (1) Loop Optimization / Estimation (`cxx:1229-1376`).
            let mut lestim = true;
            lconst = true;
            let mut nb_est = 0i32;

            self.smooth_criterion.set_curve(&ccurrent);

            while lestim {
                // (1.1) Curve's Optimization.
                eps_length = SMALL_VALUE * cblong / nbr_pnt as f64;
                let cnew = {
                    let cu = ccurrent.borrow();
                    let mut c = Curve::new(
                        cu.dimension() as i32,
                        cu.nb_elements() as i32,
                        cu.base(),
                        eps_length,
                    );
                    c.knots_mut().copy_from_slice(cu.knots());
                    Rc::new(RefCell::new(c))
                };

                self.smooth_criterion.set_parameters(&current_ti);
                eps_deg = (w_quality * 0.1).min(cblong * 0.001);

                self.optimization(&mut the_assembly, lconst, eps_deg, &cnew, &current_ti);
                lconst = false;

                // (1.2) Quality criteria and estimation improvement.
                let mut v0 = valcri[0];
                let mut v1 = valcri[1];
                let mut v2 = valcri[2];
                let icdana =
                    self.smooth_criterion
                        .quality_values(j1min, j2min, j3min, &mut v0, &mut v1, &mut v2);
                valcri = [v0, v1, v2];
                // OCCT assigns lconst = true here and overwrites it at the top
                // of L1000 before any read (cxx:1230 / cxx:1358); kept
                // verbatim.
                #[allow(unused_assignments)]
                if icdana > 0 {
                    lconst = true;
                }

                let mut emax = 0.0;
                let mut equa = 0.0;
                let mut emoy = 0.0;
                self.smooth_criterion
                    .error_values(&mut emax, &mut equa, &mut emoy);
                errmax = emax;
                errqua = equa;
                errmoy = emoy;

                isnear = ((errqua / nbr_pnt as f64).sqrt() < 2.0 * w_quality)
                    && (self.nb_iterations > 1);

                // (1.3) Optimization of ti by orthogonal projection.
                let mut num_pnt = 0i32;
                let new_ti = if isnear {
                    let mut nti = MathVector::new(1, current_ti.len());
                    self.project(
                        &cnew,
                        &current_ti,
                        &mut nti,
                        ecarts,
                        &mut num_pnt,
                        &mut errmax,
                        &mut errqua,
                        &mut errmoy,
                        2,
                    );
                    nti
                } else {
                    current_ti.clone()
                };

                // (1.4) Progression's test.
                let mut iprog = 0i32;
                if erold > w_quality && errmax < 0.95 * erold {
                    iprog += 1;
                }
                if erold > w_quality && errmax < 0.8 * erold {
                    iprog += 1;
                }
                if erold > w_quality && errmax < w_quality {
                    iprog += 1;
                }
                if erold > w_quality && errmax < 0.99 * erold && errmax < 1.1 * w_quality {
                    iprog += 1;
                }
                if valcri[0] < 0.975 * vocri[0] {
                    iprog += 1;
                }
                if valcri[0] < 0.9 * vocri[0] {
                    iprog += 1;
                }
                if valcri[1] < 0.95 * vocri[1] {
                    iprog += 1;
                }
                if valcri[1] < 0.8 * vocri[1] {
                    iprog += 1;
                }
                if valcri[2] < 0.95 * vocri[2] {
                    iprog += 1;
                }
                if valcri[2] < 0.8 * vocri[2] {
                    iprog += 1;
                }
                if vocri[1] > SMALLEST_VALUE && vocri[2] > SMALLEST_VALUE {
                    if valcri[1] / vocri[1] + 2.0 * valcri[2] / vocri[2] < 2.8 {
                        iprog += 1;
                    }
                }

                if iprog < 2 && nb_est == 0 {
                    // (1.5) Invalidation of new knots, then `goto L8000`.
                    valcri = vocri;
                    errmax = erold;
                    cblong = lnold;
                    ccurrent = cold.clone();
                    current_ti = old_ti.clone();
                    break 'again;
                }

                vocri = valcri;
                lnold = cblong;
                erold = errmax;

                ccurrent = cnew;
                current_ti = new_ti;

                // (1.6) Test if the estimations seem OK, else repeat.
                nb_est += 1;
                lestim = nb_est < MAX_NB_EST && icdana == 2 && iprog > 0;

                if lestim && isnear {
                    // (1.7) Optimization of ti by ACR.
                    current_ti.as_mut_slice().sort_by(|a, b| a.partial_cmp(b).unwrap());

                    let decima = 4;
                    cblong = ccurrent.borrow_mut().length_of(0.0, 1.0);
                    *self.smooth_criterion.est_length() = cblong;

                    self.acr(&mut ccurrent, &mut current_ti, decima);
                    lconst = true;
                }
            }

            // (2) Loop of parametric / geometric optimization.
            let mut iter = 1;
            let mut to_optim = iter < self.nb_iterations && isnear;

            while to_optim {
                iter += 1;
                // (2.1) Save current results.
                vocri = valcri;
                erold = errmax;
                lnold = cblong;
                cold = ccurrent.clone();
                old_ti = current_ti.clone();

                // (2.2) Optimization of ti by ACR.
                current_ti.as_mut_slice().sort_by(|a, b| a.partial_cmp(b).unwrap());

                let decima = 4;
                cblong = ccurrent.borrow_mut().length_of(0.0, 1.0);
                *self.smooth_criterion.est_length() = cblong;

                self.acr(&mut ccurrent, &mut current_ti, decima);
                lconst = true;

                // (2.3) Optimization of curves.
                eps_length = SMALL_VALUE * cblong / nbr_pnt as f64;
                let cnew = {
                    let cu = ccurrent.borrow();
                    let mut c = Curve::new(
                        cu.dimension() as i32,
                        cu.nb_elements() as i32,
                        cu.base(),
                        eps_length,
                    );
                    c.knots_mut().copy_from_slice(cu.knots());
                    Rc::new(RefCell::new(c))
                };

                self.smooth_criterion.set_parameters(&current_ti);
                eps_deg = (w_quality * 0.1).min(cblong * 0.001);
                self.optimization(&mut the_assembly, lconst, eps_deg, &cnew, &current_ti);

                ccurrent = cnew;

                // (2.4) Quality criteria and estimation improvement.
                let mut v0 = valcri[0];
                let mut v1 = valcri[1];
                let mut v2 = valcri[2];
                let icdana =
                    self.smooth_criterion
                        .quality_values(j1min, j2min, j3min, &mut v0, &mut v1, &mut v2);
                valcri = [v0, v1, v2];
                // OCCT assigns lconst = true here and overwrites it at the top
                // of L1000 before any read (cxx:1230 / cxx:1358); kept
                // verbatim.
                #[allow(unused_assignments)]
                if icdana > 0 {
                    lconst = true;
                }

                let (_e1, _e2, _e3) = self.smooth_criterion.estimation();

                // (2.5) Optimization of ti by orthogonal projection.
                let mut num_pnt = 0i32;
                let mut new_ti = MathVector::new(1, current_ti.len());
                self.project(
                    &ccurrent,
                    &current_ti,
                    &mut new_ti,
                    ecarts,
                    &mut num_pnt,
                    &mut errmax,
                    &mut errqua,
                    &mut errmoy,
                    2,
                );

                // (2.6) Test de non regression.
                let mut iregre = 0i32;
                if nbr_constraint < nbr_pnt {
                    if errmax > w_quality && errmax > 1.05 * erold {
                        iregre += 1;
                    }
                    if errmax > w_quality && errmax > 2.0 * erold {
                        iregre += 1;
                    }
                    if erold > w_quality && errmax <= 0.5 * erold {
                        iregre -= 1;
                    }
                }
                let (ne1, ne2, ne3) = self.smooth_criterion.estimation();
                if valcri[0] > ne1 && valcri[0] > 1.1 * vocri[0] {
                    iregre += 1;
                }
                if valcri[1] > ne2 && valcri[1] > 1.1 * vocri[1] {
                    iregre += 1;
                }
                if valcri[2] > ne3 && valcri[2] > 1.1 * vocri[2] {
                    iregre += 1;
                }

                if iregre >= 2 {
                    // (2.7) Restore the previous iteration.
                    valcri = vocri;
                    errmax = erold;
                    cblong = lnold;
                    ccurrent = cold.clone();
                    current_ti = old_ti.clone();
                    to_optim = false;
                } else {
                    // Iteration is OK: `CCurrent = CNew; CurrentTi = NewTi;`.
                    current_ti = new_ti;
                }
                if iter >= self.nb_iterations {
                    to_optim = false;
                }
            }

            // (3) Optional splitting (`cxx:1501-1552`).
            if ccurrent.borrow().nb_elements() < self.max_segment as usize && self.with_cutting {
                // (3.1) Save the previous state.
                vocri = valcri;
                erold = errmax;
                cold = ccurrent.clone();
                old_ti = current_ti.clone();

                // (3.2) Arrange the ti: sort + rescale to (0, 1).
                current_ti.as_mut_slice().sort_by(|a, b| a.partial_cmp(b).unwrap());

                if current_ti.value(1) != 0.0 || current_ti.value(nbr_pnt as usize) != 1.0 {
                    let delat_t = 1.0 / (current_ti.value(nbr_pnt as usize) - current_ti.value(1));
                    let t1 = current_ti.value(1);
                    for ii in 2..nbr_pnt as usize {
                        let t = (current_ti.value(ii) - t1) * delat_t;
                        current_ti.set_value(ii, t);
                    }
                    current_ti.set_value(1, 0.0);
                    current_ti.set_value(nbr_pnt as usize, 1.0);
                }

                // (3.3) Insert new knots.
                let (cnew, iscut) =
                    self.split_curve(&ccurrent, &current_ti, eps_length);
                if !iscut {
                    break 'again;
                }
                ccurrent = cnew;
                // New knots => new assembly.
                self.smooth_criterion.set_curve(&ccurrent);
                the_assembly =
                    Assembly::new(&dependence, &self.smooth_criterion.assembly_table());
            } else {
                break 'again;
            }
        }

        // (4) Compute the best error (`cxx:1556-1568`).
        let mut num_pnt = 0i32;
        let mut new_ti = MathVector::new(1, current_ti.len());
        self.project(
            &ccurrent,
            &current_ti,
            &mut new_ti,
            ecarts,
            &mut num_pnt,
            &mut errmax,
            &mut errqua,
            &mut errmoy,
            10,
        );

        // (5) Field's update (`cxx:1570-1589`).
        *the_curve = ccurrent;
        *self.smooth_criterion.est_length() = cblong;
        for k in 0..new_ti.len() {
            self.parameters
                .set_value(self.first_point as usize + k, new_ti.value(k + 1));
        }
        self.criterium[0] = errqua;
        self.criterium[1] = valcri[0].sqrt();
        self.criterium[2] = valcri[1].sqrt();
        self.criterium[3] = valcri[2].sqrt();
        self.max_error = errmax;
        self.max_error_index = num_pnt;
        if nbr_pnt > nbr_constraint {
            self.average_error = errmoy / (nbr_pnt - nbr_constraint) as f64;
        } else {
            self.average_error = errmoy / nbr_constraint as f64;
        }
    }

    /// `Adjusting(J, WQuadratic, WQuality, TheCurve, Ecarts)` (`cxx:2816-2973`).
    pub(super) fn adjusting(
        &mut self,
        mut w_quadratic: f64,
        w_quality: f64,
        the_curve: &mut Rc<RefCell<Curve>>,
        ecarts: &mut MathVector,
    ) {
        const MXITER: i32 = 2;
        const EPS1: f64 = 1.0e-6;
        let nbr_pnt = self.last_point - self.first_point + 1;
        let nbr_constraint = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;
        let curv_tol = EPS1 * *self.smooth_criterion.est_length() / nbr_pnt as f64;

        let mut tbpoid = MathVector::with_init(
            self.first_point as usize,
            self.last_point as usize,
            1.0,
        );

        let mut loptim = true;
        let mut iter = 0;

        let mut vtest = w_quality * 0.9;
        let mut j1cibl =
            (self.criterium[0] / (nbr_pnt - nbr_constraint) as f64).sqrt();

        let mut ccurve = the_curve.clone();

        while loptim {
            iter += 1;

            // (1) Save the previous state.
            let vocri = self.criterium;
            let erold = self.max_error;
            let emold = self.average_error;

            // (2) Increase the least squares weight.
            if j1cibl > vtest {
                w_quadratic = j1cibl / vtest * w_quadratic;
            }

            // (3) Increase the weight associated with problematic points.
            let vseuil = w_quality * 0.88;
            for ipnt in self.first_point..=self.last_point {
                if ecarts.value(ipnt as usize) > vtest {
                    let ecart = (ecarts.value(ipnt as usize) - vseuil) / w_quality;
                    let w = (ecart * 3.0 + 1.0) * tbpoid.value(ipnt as usize);
                    tbpoid.set_value(ipnt as usize, w);
                }
            }

            // (4) Forced splitting.
            let cnew = if ccurve.borrow().nb_elements() < self.max_segment as usize
                && self.with_cutting
            {
                let (numint, _flag) = near_index(
                    self.param(self.max_error_index),
                    ccurve.borrow().knots(),
                    0.0,
                );
                let numint = numint as usize;
                let tpara = (ccurve.borrow().knots()[numint]
                    + ccurve.borrow().knots()[numint + 1]
                    + self.param(self.max_error_index) * 2.0)
                    / 4.0;

                let cu = ccurve.borrow();
                let mut c = Curve::new(
                    self.dimension,
                    cu.nb_elements() as i32 + 1,
                    cu.base(),
                    curv_tol,
                );
                let n = cu.knots().len();
                for i in 1..=numint {
                    c.knots_mut()[i] = cu.knots()[i];
                }
                for i in numint + 1..=n - 1 {
                    c.knots_mut()[i + 1] = cu.knots()[i];
                }
                c.knots_mut()[numint + 1] = tpara;
                Rc::new(RefCell::new(c))
            } else {
                let cu = ccurve.borrow();
                let mut c = Curve::new(
                    self.dimension,
                    cu.nb_elements() as i32,
                    cu.base(),
                    curv_tol,
                );
                c.knots_mut().copy_from_slice(cu.knots());
                Rc::new(RefCell::new(c))
            };

            let mut jnew: Box<dyn SmoothCriterion> =
                Box::new(LinearCriteria::new(&self.ssp, self.first_point, self.last_point));
            *jnew.est_length() = *self.smooth_criterion.est_length();
            let (e1, e2, e3) = self.smooth_criterion.estimation();
            jnew.set_estimation(e1, e2, e3);
            jnew.set_parameters(&self.parameters);
            jnew.set_weight(
                w_quadratic,
                w_quality,
                [self.percent[0], self.percent[1], self.percent[2]],
            );
            jnew.set_pnt_weight(&tbpoid);
            jnew.set_curve(&cnew);

            // (5) Re-smoothing (`cxx:2939`): `TheMotor(JNew, ...)` with
            // `JNew` as the active criterion.
            let saved = std::mem::replace(&mut self.smooth_criterion, jnew);
            self.the_motor(w_quality, &mut ccurve, ecarts);
            let jnew = std::mem::replace(&mut self.smooth_criterion, saved);

            // (6) Rejection tests.
            j1cibl = (self.criterium[0] / (nbr_pnt - nbr_constraint) as f64).sqrt();
            let vseuil2 = vocri[1].sqrt() + (erold - self.max_error) * 4.0;

            let lrejet = (self.max_error > w_quality && self.max_error > erold * 1.01)
                || (self.criterium[1].sqrt() > vseuil2 * 1.05);

            if lrejet {
                self.criterium = vocri;
                self.max_error = erold;
                self.average_error = emold;
                loptim = false;
            } else {
                self.smooth_criterion = jnew;
                self.smooth_criterion.set_curve(&ccurve);
            }

            // (7) Convergence test.
            if (iter >= MXITER && self.max_segment == ccurve.borrow().nb_elements() as i32)
                || self.max_error < w_quality
            {
                loptim = false;
            }
            let _ = &mut vtest;
        }

        *the_curve = ccurve;
    }
}
