//! `AppDef_Variational::Init` and the initialization helpers.
//!
//! Sources: `AppDef_Variational.cxx:224-487` (`Init`), `cxx:2129-2204`
//! (`InitSmoothCriterion`), `cxx:2209-2253` (`InitParameters`),
//! `cxx:2256-2395` (`InitCriterionEstimations`), `cxx:2401-2581`
//! (`EstTangent`), `cxx:2585-2679` (`EstSecnd`), `cxx:2684-2813`
//! (`InitCutting`), `cxx:3315-3403` (`InitTthetaF`) and the file-local
//! `NotParallel` helper at `cxx:2976-2991`.

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::gp::{GpVec, GpVec2d};
use occt_math::fem_tool::Curve;
use occt_math::plib_hermit_jacobi::HermitJacobi;
use occt_math::vector::MathVector;

use super::{Variational, ANGULAR, CONFUSION};
use crate::app_def::my_line_tool::MyLineTool;
use crate::app_par_curves::constraint::Constraint;
impl Variational {
    /// `Init()` (`cxx:224-487`). Returns `Err(())` where OCCT throws
    /// `Standard_ConstructionError`.
    pub(super) fn init(&mut self) -> Result<(), ()> {
        let nb_constraints = self.constraints.len() as i32;
        if nb_constraints < 0 {
            return Err(()); // `cxx:239-242`.
        }
        self.nb_constraints = nb_constraints;

        // `cxx:243-249`.
        self.typ_constraints = vec![0i32; (2 * nb_constraints).max(1) as usize + 1];
        self.tab_constraints =
            vec![0.0f64; (2 * self.dimension * nb_constraints).max(1) as usize + 1];
        let n_theta = ((2 * self.nb_p2d + 6 * self.nb_p3d) * nb_constraints).max(1) as usize;
        self.ttheta = vec![0.0f64; n_theta + 1];
        self.tfthet = vec![0.0f64; n_theta + 1];

        let mut index = 1usize;
        let mut jndex = 1usize;
        let mut cur_multy_point: i32 = 1;
        self.nb_pass_points = 0;
        self.nb_tang_points = 0;
        self.nb_curv_points = 0;

        for iconstr in 0..self.constraints.len() {
            let couple = self.constraints[iconstr];
            let ipoint = couple.index();
            let valcontr = couple.constraint();
            let (nb_p2d, nb_p3d) = (self.nb_p2d, self.nb_p3d);
            match valcontr {
                Constraint::NoConstraint => {
                    // `cxx:270`: the placeholder entry does not advance the
                    // global position (undone after the match).
                    cur_multy_point -= nb_p3d * 6 + nb_p2d * 2;
                }
                Constraint::PassPoint => {
                    self.typ_constraints[index] = ipoint;
                    index += 1;
                    self.typ_constraints[index] = 0;
                    index += 1;
                    self.nb_pass_points += 1;
                    if nb_p2d != 0 {
                        jndex += (4 * nb_p2d) as usize;
                    }
                    if nb_p3d != 0 {
                        jndex += (6 * nb_p3d) as usize;
                    }
                }
                Constraint::TangencyPoint => {
                    self.typ_constraints[index] = ipoint;
                    index += 1;
                    self.typ_constraints[index] = 1;
                    index += 1;
                    self.nb_tang_points += 1;
                    if nb_p2d != 0 && nb_p3d == 0 {
                        // `cxx:289-305`.
                        let tab_v2d = match MyLineTool::tangency(&self.ssp, ipoint) {
                            Some((_, v2d)) => v2d,
                            None => return Err(()),
                        };
                        for jp2d in 0..nb_p2d as usize {
                            let vt2d = match tab_v2d[jp2d].normalized() {
                                Ok(v) => v,
                                Err(_) => return Err(()),
                            };
                            self.tab_constraints[jndex] = vt2d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt2d.y();
                            jndex += 1;
                            jndex += 2;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + jp2d as i32 * 2,
                                (jndex - 4) as i32,
                            );
                        }
                    }
                    if nb_p3d != 0 && nb_p2d == 0 {
                        // `cxx:306-324`.
                        let tab_v3d = match MyLineTool::tangency(&self.ssp, ipoint) {
                            Some((v3d, _)) => v3d,
                            None => return Err(()),
                        };
                        for jp3d in 0..nb_p3d as usize {
                            let vt3d = tab_v3d[jp3d].normalized();
                            self.tab_constraints[jndex] = vt3d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.z();
                            jndex += 1;
                            jndex += 3;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + jp3d as i32 * 6,
                                (jndex - 6) as i32,
                            );
                        }
                    }
                    if nb_p3d != 0 && nb_p2d != 0 {
                        // `cxx:325-350`.
                        let (tab_v3d, tab_v2d) =
                            match MyLineTool::tangency(&self.ssp, ipoint) {
                                Some(v) => v,
                                None => return Err(()),
                            };
                        for jp3d in 0..nb_p3d as usize {
                            let vt3d = tab_v3d[jp3d].normalized();
                            self.tab_constraints[jndex] = vt3d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.z();
                            jndex += 1;
                            jndex += 3;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + jp3d as i32 * 6,
                                (jndex - 6) as i32,
                            );
                        }
                        for jp2d in 0..nb_p2d as usize {
                            let vt2d = match tab_v2d[jp2d].normalized() {
                                Ok(v) => v,
                                Err(_) => return Err(()),
                            };
                            self.tab_constraints[jndex] = vt2d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt2d.y();
                            jndex += 1;
                            jndex += 2;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + jp2d as i32 * 2 + nb_p3d * 6,
                                (jndex - 4) as i32,
                            );
                        }
                    }
                }
                Constraint::CurvaturePoint => {
                    self.typ_constraints[index] = ipoint;
                    index += 1;
                    self.typ_constraints[index] = 2;
                    index += 1;
                    self.nb_curv_points += 1;
                    if nb_p2d != 0 && nb_p3d == 0 {
                        // `cxx:358-383`.
                        let tab_v2d = match MyLineTool::tangency(&self.ssp, ipoint) {
                            Some((_, v2d)) => v2d,
                            None => return Err(()),
                        };
                        let tab_v2d_curv = match MyLineTool::curvature_2d(&self.ssp, ipoint) {
                            Some(v2d) => v2d,
                            None => return Err(()),
                        };
                        for jp2d in 0..nb_p2d as usize {
                            let vt2d = match tab_v2d[jp2d].normalized() {
                                Ok(v) => v,
                                Err(_) => return Err(()),
                            };
                            let vc2d: GpVec2d = tab_v2d_curv[jp2d];
                            if (vc2d.angle(&vt2d).abs() - std::f64::consts::PI / 2.0).abs()
                                > ANGULAR
                            {
                                return Err(());
                            }
                            self.tab_constraints[jndex] = vt2d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt2d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc2d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc2d.y();
                            jndex += 1;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + jp2d as i32 * 2,
                                (jndex - 4) as i32,
                            );
                        }
                    }
                    if nb_p3d != 0 && nb_p2d == 0 {
                        // `cxx:385-411`.
                        let tab_v3d = match MyLineTool::tangency(&self.ssp, ipoint) {
                            Some((v3d, _)) => v3d,
                            None => return Err(()),
                        };
                        let tab_v3d_curv = match MyLineTool::curvature_3d(&self.ssp, ipoint) {
                            Some(v3d) => v3d,
                            None => return Err(()),
                        };
                        for jp3d in 0..nb_p3d as usize {
                            let vt3d = tab_v3d[jp3d].normalized();
                            let vc3d: GpVec = tab_v3d_curv[jp3d];
                            if !is_normal(&vc3d.normalized(), &vt3d) {
                                return Err(());
                            }
                            self.tab_constraints[jndex] = vt3d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.z();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc3d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc3d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc3d.z();
                            jndex += 1;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + jp3d as i32 * 6,
                                (jndex - 6) as i32,
                            );
                        }
                    }
                    if nb_p3d != 0 && nb_p2d != 0 {
                        // `cxx:412-455`.
                        let (tab_v3d, tab_v2d) =
                            match MyLineTool::tangency(&self.ssp, ipoint) {
                                Some(v) => v,
                                None => return Err(()),
                            };
                        let (tab_v3d_curv, tab_v2d_curv) =
                            match MyLineTool::curvature(&self.ssp, ipoint) {
                                Some(v) => v,
                                None => return Err(()),
                            };
                        for jp3d in 0..nb_p3d as usize {
                            let vt3d = tab_v3d[jp3d].normalized();
                            let vc3d: GpVec = tab_v3d_curv[jp3d];
                            if !is_normal(&vc3d.normalized(), &vt3d) {
                                return Err(());
                            }
                            self.tab_constraints[jndex] = vt3d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt3d.z();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc3d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc3d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc3d.z();
                            jndex += 1;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + jp3d as i32 * 6,
                                (jndex - 6) as i32,
                            );
                        }
                        for jp2d in 0..nb_p2d as usize {
                            let vt2d = match tab_v2d[jp2d].normalized() {
                                Ok(v) => v,
                                Err(_) => return Err(()),
                            };
                            let vc2d: GpVec2d = tab_v2d_curv[jp2d];
                            if (vc2d.angle(&vt2d).abs() - std::f64::consts::PI / 2.0).abs()
                                > ANGULAR
                            {
                                return Err(());
                            }
                            self.tab_constraints[jndex] = vt2d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vt2d.y();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc2d.x();
                            jndex += 1;
                            self.tab_constraints[jndex] = vc2d.y();
                            jndex += 1;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + jp2d as i32 * 2 + nb_p3d * 6,
                                (jndex - 4) as i32,
                            );
                        }
                    }
                }
            }
            cur_multy_point += self.nb_p3d * 6 + self.nb_p2d * 2;
        }

        // Over-constraint detection (`cxx:462-487`).
        let max_seg = if self.with_cutting { self.max_segment } else { 1 };
        let bound = (self.max_degree - self.niv_cont) * max_seg
            - self.nb_pass_points
            - 2 * self.nb_tang_points
            - 3 * self.nb_curv_points;
        if bound < 0 {
            self.is_over_constr = true;
            self.is_created = false;
        } else {
            self.init_smooth_criterion()?;
            self.is_created = true;
        }
        Ok(())
    }

    /// `InitSmoothCriterion()` (`cxx:2129-2205`).
    pub(super) fn init_smooth_criterion(&mut self) -> Result<(), ()> {
        const EPS2: f64 = 1.0e-6;
        const EPS3: f64 = 1.0e-9;

        let length = self.init_parameters()?;
        self.smooth_criterion.set_parameters(&self.parameters);

        let (e1, e2, e3) = self.init_criterion_estimations(length);

        *self.smooth_criterion.est_length() = length;
        self.smooth_criterion.set_estimation(e1, e2, e3);

        let w_quality = if !self.with_min_max && self.tolerance != 0.0 {
            self.tolerance
        } else if self.tolerance == 0.0 {
            1.0
        } else {
            self.tolerance.max(EPS2 * length)
        };

        let nb_constr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;
        let mut w_quadratic = ((self.nb_points - nb_constr) as f64).sqrt() * w_quality;
        if w_quadratic > EPS3 {
            w_quadratic = 1.0 / w_quadratic;
        }
        if w_quadratic == 0.0 {
            w_quadratic = e1.sqrt().max(1.0);
        }
        self.smooth_criterion.set_weight(
            w_quadratic,
            w_quality,
            [self.percent[0], self.percent[1], self.percent[2]],
        );

        let the_base =
            HermitJacobi::new(self.max_degree, self.continuity).map_err(|_| ())?;
        let curv_tol = EPS2 * length / self.nb_points as f64;

        let the_curve = if self.with_cutting && nb_constr != 0 {
            self.init_cutting(&the_base, curv_tol)?
        } else {
            let mut c = Curve::new(self.dimension, 1, &the_base, curv_tol);
            let n = c.knots().len();
            c.knots_mut()[1] = self.param(self.first_point);
            c.knots_mut()[n - 1] = self.param(self.last_point);
            c
        };

        self.smooth_criterion
            .set_curve(&Rc::new(RefCell::new(the_curve)));
        Ok(())
    }

    /// `InitParameters(Length)` (`cxx:2209-2253`).
    fn init_parameters(&mut self) -> Result<f64, ()> {
        let eps1 = CONFUSION * 0.01;

        let mut length = 0.0;
        self.set_param(self.first_point, length);

        let mut i1 = 0usize;
        for ipoint in (self.first_point + 1)..=self.last_point {
            let i0 = i1;
            i1 += self.dimension as usize;
            let mut dist = 0.0;
            for i in 1..=self.dimension as usize {
                let aux = self.tab_points[i1 + i] - self.tab_points[i0 + i];
                dist += aux * aux;
            }
            length += dist.sqrt();
            self.set_param(ipoint, length);
        }
        if length <= eps1 {
            return Err(());
        }

        for ipoint in (self.first_point + 1)..self.last_point {
            self.set_param(ipoint, self.param(ipoint) / length);
        }
        self.set_param(self.last_point, 1.0);

        if self.nb_points < 10 {
            length *= 1.0 + 0.1 / (self.nb_points - 1) as f64;
        }
        Ok(length)
    }

    /// `InitCriterionEstimations(Length, E1, E2, E3)` (`cxx:2256-2395`).
    fn init_criterion_estimations(&self, length: f64) -> (f64, f64, f64) {
        let eps1 = CONFUSION * 0.01;
        let dim = self.dimension as usize;

        let mut e1 = length * length;

        let mut vtang1 = MathVector::new(1, dim);
        let mut vtang2 = MathVector::new(1, dim);
        let mut vtang3 = MathVector::new(1, dim);
        let mut vscnd1 = MathVector::new(1, dim);
        let mut vscnd2 = MathVector::new(1, dim);
        let mut vscnd3 = MathVector::new(1, dim);

        // First point (`cxx:2268-2295`).
        let mut ipnt = self.first_point;
        self.est_tangent(ipnt, &mut vtang1);
        ipnt += 1;
        self.est_tangent(ipnt, &mut vtang2);
        ipnt += 1;
        self.est_tangent(ipnt, &mut vtang3);

        ipnt = self.first_point;
        self.est_secnd(ipnt, &vtang1, &vtang2, length, &mut vscnd1);
        ipnt += 1;
        self.est_secnd(ipnt, &vtang1, &vtang3, length, &mut vscnd2);

        let an_ind = ipnt;
        let mut delta = 0.5 * (self.param(an_ind) - self.param(ipnt - 1));
        if delta <= eps1 {
            delta = 1.0;
        }

        let mut e2 = vscnd1.norm2() * delta;
        let mut e3 = if delta > eps1 {
            vscnd2.subtracted(&vscnd1).norm2() / (4.0 * delta)
        } else {
            0.0
        };

        // Internal points (`cxx:2298-2357`).
        let mut curr_point = 2;
        for ipnt in (self.first_point + 1)..self.last_point {
            delta = 0.5 * (self.param(ipnt + 1) - self.param(ipnt - 1));

            if curr_point == 1 {
                if ipnt + 1 != self.last_point {
                    self.est_tangent(ipnt + 2, &mut vtang3);
                    self.est_secnd(ipnt + 1, &vtang1, &vtang3, length, &mut vscnd2);
                } else {
                    self.est_secnd(ipnt + 1, &vtang1, &vtang2, length, &mut vscnd2);
                }
                e2 += vscnd1.norm2() * delta;
                e3 += if delta > eps1 {
                    vscnd2.subtracted(&vscnd3).norm2() / (4.0 * delta)
                } else {
                    0.0
                };
            } else if curr_point == 2 {
                if ipnt + 1 != self.last_point {
                    self.est_tangent(ipnt + 2, &mut vtang1);
                    self.est_secnd(ipnt + 1, &vtang2, &vtang1, length, &mut vscnd3);
                } else {
                    self.est_secnd(ipnt + 1, &vtang2, &vtang3, length, &mut vscnd3);
                }
                e2 += vscnd2.norm2() * delta;
                e3 += if delta > eps1 {
                    vscnd3.subtracted(&vscnd1).norm2() / (4.0 * delta)
                } else {
                    0.0
                };
            } else {
                if ipnt + 1 != self.last_point {
                    self.est_tangent(ipnt + 2, &mut vtang2);
                    self.est_secnd(ipnt + 1, &vtang3, &vtang2, length, &mut vscnd1);
                } else {
                    self.est_secnd(ipnt + 1, &vtang3, &vtang1, length, &mut vscnd1);
                }
                e2 += vscnd3.norm2() * delta;
                e3 += if delta > eps1 {
                    vscnd1.subtracted(&vscnd2).norm2() / (4.0 * delta)
                } else {
                    0.0
                };
            }

            curr_point += 1;
            if curr_point == 4 {
                curr_point = 1;
            }
        }

        // Last point (`cxx:2358-2393`).
        delta = 0.5 * (self.param(self.last_point) - self.param(self.last_point - 1));
        if delta <= eps1 {
            delta = 1.0;
        }

        let aux0;
        if curr_point == 1 {
            e2 += vscnd1.norm2() * delta;
            aux0 = vscnd1.subtracted(&vscnd3).norm2();
            e3 += if delta > eps1 { aux0 / (4.0 * delta) } else { aux0 };
        } else if curr_point == 2 {
            e2 += vscnd2.norm2() * delta;
            aux0 = vscnd2.subtracted(&vscnd1).norm2();
            e3 += if delta > eps1 { aux0 / (4.0 * delta) } else { aux0 };
        } else {
            e2 += vscnd3.norm2() * delta;
            aux0 = vscnd3.subtracted(&vscnd2).norm2();
            e3 += if delta > eps1 { aux0 / (4.0 * delta) } else { aux0 };
        }

        let aux = length * length;
        e2 *= aux;
        e3 *= aux;
        e1 *= 1.0;
        (e1, e2, e3)
    }
}


/// `gp_Vec::IsNormal(other, tol)` (`gp_Vec.cxx`): `|Dot| <= |other| * 1 * tol`.
fn is_normal(vc_normalized: &GpVec, vt: &GpVec) -> bool {
    vc_normalized.dot(vt).abs() <= vt.magnitude() * ANGULAR
}
