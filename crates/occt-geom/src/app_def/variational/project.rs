//! `AppDef_Variational::Project` — orthogonal projection of the point
//! parameters onto the curve.
//!
//! Source: `AppDef_Variational.cxx:1691-1827` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `C->D0/D1/D2` take zero-based arrays here; OCCT passes
//!   `NCollection_Array1(1, myDimension)` and reads `ValOfC(i)`, so the
//!   element `i` of this port is `val_of_c[i - 1]`.
//! - `NumPoints` is not initialized by OCCT before the loop; the caller
//!   starts it at 0 (see [`super::motor::Variational::the_motor`]).

use std::cell::RefCell;
use std::rc::Rc;

use occt_math::fem_tool::Curve;
use occt_math::vector::MathVector;

use super::Variational;

impl Variational {
    /// `Project(C, Ti, ProjTi, Distance, NumPoints, MaxErr, QuaErr, AveErr,
    /// NbIterations)` (`cxx:1691-1827`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn project(
        &self,
        c: &Rc<RefCell<Curve>>,
        ti: &MathVector,
        proj_ti: &mut MathVector,
        distance: &mut MathVector,
        num_points: &mut i32,
        max_err: &mut f64,
        qua_err: &mut f64,
        ave_err: &mut f64,
        nb_iterations: i32,
    ) {
        const SEUIL: f64 = 1.0e-9;
        const EPS: f64 = 1.0e-12;

        *max_err = 0.0;
        *qua_err = 0.0;
        *ave_err = 0.0;

        let dim = self.dimension;
        let d0 = distance.lower() as i32 - 1;
        let mut val_of_c = vec![0.0f64; dim as usize];
        let mut first_der = vec![0.0f64; dim as usize];
        let mut secnd_der = vec![0.0f64; dim as usize];

        for ipnt in 1..=proj_ti.len() as i32 {
            let i0 = (ipnt - 1) * dim;

            let mut t_new = ti.value(ipnt as usize);

            let mut en_cour = true;
            let mut nit_cv = 0;
            let mut iter = 0;
            c.borrow_mut().d0(t_new, &mut val_of_c);

            let mut dist = 0.0;
            for i in 0..dim as usize {
                let aux = val_of_c[i] - self.tab_points[(i0 + i as i32 + 1) as usize];
                dist += aux * aux;
            }
            dist = dist.sqrt();

            // Newton's method for solving (C'(t), C(t) - P) = 0.
            while en_cour {
                iter += 1;
                let t0 = t_new;
                let dist0 = dist;

                {
                    let mut cu = c.borrow_mut();
                    cu.d2(t_new, &mut secnd_der);
                    cu.d1(t_new, &mut first_der);
                }

                let mut f1 = 0.0;
                let mut f2 = 0.0;
                for i in 0..dim as usize {
                    let aux = val_of_c[i] - self.tab_points[(i0 + i as i32 + 1) as usize];
                    let df = first_der[i];
                    f1 += aux * df;
                    f2 += df * df + aux * secnd_der[i];
                }

                if f2.abs() < EPS {
                    en_cour = false;
                } else {
                    t_new -= f1 / f2;
                    if t_new < 0.0 {
                        t_new = 0.0;
                    }
                    if t_new > 1.0 {
                        t_new = 1.0;
                    }

                    // Analysis of the result.
                    c.borrow_mut().d0(t_new, &mut val_of_c);

                    dist = 0.0;
                    for i in 0..dim as usize {
                        let aux = val_of_c[i] - self.tab_points[(i0 + i as i32 + 1) as usize];
                        dist += aux * aux;
                    }
                    dist = dist.sqrt();

                    let ecart = dist0 - dist;

                    if ecart <= -SEUIL {
                        // No improvement: stop here.
                        en_cour = false;
                        t_new = t0;
                        dist = dist0;
                    } else if ecart <= SEUIL {
                        // Convergence.
                        nit_cv += 1;
                    } else {
                        nit_cv = 0;
                    }

                    if nit_cv >= 2 || iter >= nb_iterations {
                        en_cour = false;
                    }
                }
            }

            proj_ti.set_value(ipnt as usize, t_new);
            distance.set_value((d0 + ipnt) as usize, dist);
            if dist > *max_err {
                *max_err = dist;
                *num_points = ipnt;
            }
            *qua_err += dist * dist;
            *ave_err += dist;
        }

        // `NumPoints` back to the interval [myFirstPoint, myLastPoint].
        *num_points = *num_points + self.first_point - 1;
    }
}
