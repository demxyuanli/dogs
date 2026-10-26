//! `AppDef_Variational::InitCutting` / `InitTthetaF`.
//!
//! Sources: `AppDef_Variational.cxx:2684-2813` (`InitCutting`),
//! `cxx:3315-3403` (`InitTthetaF`) and the file-local `NotParallel`
//! helper at `cxx:2976-2991`.

use occt_core::gp::GpVec;
use occt_math::fem_tool::Curve;
use occt_math::plib_hermit_jacobi::HermitJacobi;

use super::Variational;
use crate::app_par_curves::constraint::Constraint;

/// `static bool NotParallel(gp_Vec& T, gp_Vec& V)` (`cxx:2976-2991`).
fn not_parallel(t: &GpVec, v: &mut GpVec) -> bool {
    *v = *t;
    *v = GpVec::new(v.x() + 1.0, v.y(), v.z());
    if v.cross_magnitude(t) > 1.0e-12 {
        return true;
    }
    *v = GpVec::new(v.x(), v.y() + 1.0, v.z());
    if v.cross_magnitude(t) > 1.0e-12 {
        return true;
    }
    *v = GpVec::new(v.x(), v.y(), v.z() + 1.0);
    v.cross_magnitude(t) > 1.0e-12
}

impl Variational {
    /// `InitCutting(aBase, CurvTol, aCurve)` (`cxx:2684-2813`).
    pub(super) fn init_cutting(
        &self,
        a_base: &HermitJacobi,
        curv_tol: f64,
    ) -> Result<Curve, ()> {
        let mut orc_mx = -1i32;
        let mut n_cont = 0i32;
        let nb_constr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;

        for i in 1..=nb_constr {
            let kk = self.typ_constraints[(2 * i) as usize].abs() + 1;
            orc_mx = orc_mx.max(kk);
            n_cont += kk;
        }
        if orc_mx > self.max_degree - self.niv_cont {
            return Err(());
        }

        let mut n_libre = (self.max_degree - self.niv_cont - (self.max_degree + 1) / 4)
            .max(self.niv_cont + 1);
        let mut nb_elem = if n_cont % n_libre == 0 {
            n_cont / n_libre
        } else {
            n_cont / n_libre + 1
        };
        while nb_elem > self.max_segment && n_libre < self.max_degree - self.niv_cont {
            n_libre += 1;
            nb_elem = if n_cont % n_libre == 0 {
                n_cont / n_libre
            } else {
                n_cont / n_libre + 1
            };
        }
        if nb_elem > self.max_segment {
            return Err(());
        }

        let mut a_curve = Curve::new(self.dimension, nb_elem, a_base, curv_tol);
        let nc = nb_elem - 1;

        let mut nc_nt = (n_cont - 1) / nb_elem + 1;
        let mut n_plus = nb_elem - (nc_nt * nb_elem - n_cont);
        let _ = nc;

        let (mut ideb, mut ifin) = (0i32, nb_constr + 1);
        let (mut ndeb, mut nfin) = (0i32, 0i32);
        let lower = 1usize;
        let upper = a_curve.knots().len() - 1;
        let mut ind_el = lower;
        let mut nb_el = 0i32;

        a_curve.knots_mut()[ind_el] = self.param(self.first_point);
        a_curve.knots_mut()[upper] = self.param(self.last_point);

        while nb_elem - nb_el > 1 {
            ind_el += 1;
            nb_el += 1;
            if n_plus == 0 {
                nc_nt -= 1;
            }

            while ndeb < nc_nt && ideb < ifin {
                ideb += 1;
                ndeb += self.typ_constraints[(2 * ideb) as usize].abs() + 1;
            }

            if ndeb == nc_nt {
                ndeb = 0;
                if n_plus == 1
                    && self.param(self.typ_constraints[(2 * ideb - 1) as usize])
                        > a_curve.knots()[ind_el - 1]
                {
                    a_curve.knots_mut()[ind_el] =
                        self.param(self.typ_constraints[(2 * ideb - 1) as usize]);
                } else {
                    let v = (self.param(self.typ_constraints[(2 * ideb - 1) as usize])
                        + self.param(self.typ_constraints[(2 * ideb + 1) as usize]))
                        / 2.0;
                    a_curve.knots_mut()[ind_el] = v;
                }
            } else {
                ndeb -= nc_nt;
                a_curve.knots_mut()[ind_el] =
                    self.param(self.typ_constraints[(2 * ideb - 1) as usize]);
            }

            n_plus -= 1;
            if n_plus == 0 {
                nc_nt -= 1;
            }

            if nb_elem - nb_el == 1 {
                break;
            }

            nb_el += 1;

            while nfin < nc_nt && ideb < ifin {
                ifin -= 1;
                nfin += self.typ_constraints[(2 * ifin) as usize].abs() + 1;
            }

            if nfin == nc_nt {
                nfin = 0;
                let v = (self.param(self.typ_constraints[(2 * ifin - 1) as usize])
                    + self.param(self.typ_constraints[(2 * ifin - 3) as usize]))
                    / 2.0;
                a_curve.knots_mut()[upper + 1 - ind_el] = v;
            } else {
                nfin -= nc_nt;
                if self.param(self.typ_constraints[(2 * ifin - 1) as usize])
                    < a_curve.knots()[upper - ind_el + 1]
                {
                    a_curve.knots_mut()[upper + 1 - ind_el] =
                        self.param(self.typ_constraints[(2 * ifin - 1) as usize]);
                } else {
                    let v = (self.param(self.typ_constraints[(2 * ifin - 1) as usize])
                        + self.param(self.typ_constraints[(2 * ifin - 3) as usize]))
                        / 2.0;
                    a_curve.knots_mut()[upper + 1 - ind_el] = v;
                }
            }
        }
        Ok(a_curve)
    }

    /// `InitTthetaF(ndimen, typcon, begin, jndex)` (`cxx:3315-3403`).
    pub(super) fn init_ttheta_f(
        &mut self,
        ndimen: i32,
        typcon: Constraint,
        begin: i32,
        jndex: i32,
    ) -> bool {
        if !(2..=3).contains(&ndimen) {
            return false;
        }

        if typcon == Constraint::TangencyPoint || typcon == Constraint::CurvaturePoint {
            let tx = self.tab_constraints[jndex as usize];
            let ty = self.tab_constraints[(jndex + 1) as usize];
            let tz = if ndimen == 3 {
                self.tab_constraints[(jndex + 2) as usize]
            } else {
                0.0
            };
            let t = GpVec::new(tx, ty, tz);

            let mut v = GpVec::zero();
            if ndimen == 2 {
                v = GpVec::new(0.0, 0.0, 1.0);
            }
            if ndimen == 3 && !not_parallel(&t, &mut v) {
                return false;
            }

            let mut theta1 = v.crossed(&t);
            theta1.normalize();
            let mut theta2 = GpVec::zero();
            self.ttheta[begin as usize] = theta1.x();
            self.ttheta[(begin + 1) as usize] = theta1.y();
            if ndimen == 3 {
                theta2 = t.crossed(&theta1);
                theta2.normalize();
                self.ttheta[(begin + 2) as usize] = theta1.z();
                self.ttheta[(begin + 3) as usize] = theta2.x();
                self.ttheta[(begin + 4) as usize] = theta2.y();
                self.ttheta[(begin + 5) as usize] = theta2.z();
            }

            if typcon == Constraint::CurvaturePoint {
                let xx = tx * tx;
                let xy = tx * ty;
                let yy = ty * ty;
                if ndimen == 2 {
                    let fx = yy * theta1.x() - xy * theta1.y();
                    let fy = xx * theta1.y() - xy * theta1.x();
                    self.tfthet[begin as usize] = fx;
                    self.tfthet[(begin + 1) as usize] = fy;
                }
                if ndimen == 3 {
                    let xz = tx * tz;
                    let yz = ty * tz;
                    let zz = tz * tz;

                    let fx = (zz + yy) * theta1.x() - xy * theta1.y() - xz * theta1.z();
                    let fy = (xx + zz) * theta1.y() - xy * theta1.x() - yz * theta1.z();
                    let fz = (xx + yy) * theta1.z() - xz * theta1.x() - yz * theta1.y();
                    self.tfthet[begin as usize] = fx;
                    self.tfthet[(begin + 1) as usize] = fy;
                    self.tfthet[(begin + 2) as usize] = fz;

                    let fx = (zz + yy) * theta2.x() - xy * theta2.y() - xz * theta2.z();
                    let fy = (xx + zz) * theta2.y() - xy * theta2.x() - yz * theta2.z();
                    let fz = (xx + yy) * theta2.z() - xz * theta2.x() - yz * theta2.y();
                    self.tfthet[(begin + 3) as usize] = fx;
                    self.tfthet[(begin + 4) as usize] = fy;
                    self.tfthet[(begin + 5) as usize] = fz;
                }
            }
        }
        true
    }
}
