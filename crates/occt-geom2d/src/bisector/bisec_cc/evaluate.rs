//! `Bisector_BisecCC` evaluation: `ValueAndDist` (`cxx:609-772`),
//! `ValueByInt` (`cxx:788-988`), `Values` (`cxx:1081-1207`) and `Extension`
//! (`cxx:1214-1323`), plus the `EvalD0..EvalDN` wrappers (`cxx:992-1059`).

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::math_fn::MathFunction;
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION, REAL_SMALL};
use occt_math::bissec_newton::MathBissecNewton;
use occt_math::function_set_root::MathFunctionRoot;

use crate::geom2d_int::Geom2dIntGInter;
use crate::line::Geom2dLine;

use super::super::bisec_pc::BisectorBisecPC;
use super::super::function_h::BisectorFunctionH;
use super::support::{curvature, discret_par, epsilon_of};
use super::BisectorBisecCC;

/// `gp_Dir2d(X, Y)` with the OCCT throw replaced by a zero direction.
fn dir2d_or_zero(x: f64, y: f64) -> GpDir2d {
    GpDir2d::new(x, y).unwrap_or(GpDir2d { x: 0.0, y: 0.0 })
}

impl BisectorBisecCC {
    /// `ValueAndDist(U, U1, U2, Dist)` (`cxx:609-772`).
    pub fn value_and_dist(&self, u: f64, u1: &mut f64, u2: &mut f64, dist: &mut f64) -> GpPnt2d {
        let mut t = GpVec2d::zero();

        if self.my_polygon.length() <= 1 {
            return self.extension(u, u1, u2, dist, &mut t);
        }
        if u < self.my_polygon.first().param_on_bis() {
            return self.extension(u, u1, u2, dist, &mut t);
        }
        if u > self.my_polygon.last().param_on_bis() {
            return self.extension(u, u1, u2, dist, &mut t);
        }

        let interval_index = self.my_polygon.interval(u);
        let u_min = self.my_polygon.value(interval_index).param_on_bis();
        let u_max = self.my_polygon.value(interval_index + 1).param_on_bis();
        let v_min_0 = self.my_polygon.value(interval_index).param_on_c2();
        let v_max_0 = self.my_polygon.value(interval_index + 1).param_on_c2();

        let v_init = if (u_max - u_min).abs() < REAL_SMALL {
            v_min_0
        } else {
            let alpha = (u - u_min) / (u_max - u_min);
            v_min_0 + alpha * (v_max_0 - v_min_0)
        };

        *u1 = self.link_bis_curve(u);
        let temp = v_min_0.min(v_max_0);
        let v_max = v_min_0.max(v_max_0);
        let v_min = temp;
        let mut valid = true;

        const EPS_H: f64 = 1.0e-9;
        const EPS_H100: f64 = 1.0e-7;
        let (p1, t1) = self.curve1.d1(*u1);
        let n1 = GpVec2d::new(t1.y(), -t1.x());

        if (v_max - v_min) < PCONFUSION {
            *u2 = v_init;
        } else {
            let mut h = BisectorFunctionH::new(
                self.curve2.clone(),
                p1,
                GpVec2d::new(self.sign1 * self.sign2 * t1.x(), self.sign1 * self.sign2 * t1.y()),
            );
            let mut f_init = 0.0f64;
            h.value(v_init, &mut f_init);
            if f_init.abs() < EPS_H {
                *u2 = v_init;
            } else {
                let mut a_new_solution = MathBissecNewton::new(EPS_H);
                a_new_solution.perform(&mut h, v_min - EPS_H100, v_max + EPS_H100, 10);
                if a_new_solution.is_done() {
                    *u2 = a_new_solution.root();
                } else {
                    let sol_root = MathFunctionRoot::new_with_bounds(
                        &mut h,
                        v_init,
                        EPS_H,
                        v_min - EPS_H100,
                        v_max + EPS_H100,
                        100,
                    );
                    if sol_root.is_done() {
                        *u2 = sol_root.root();
                    } else {
                        valid = false;
                    }
                }
            }
        }

        let mut p_bis = self.point_start;
        if valid {
            let p2 = self.curve2.value(*u2);
            let p2p1 = GpVec2d::new(p1.x() - p2.x(), p1.y() - p2.y());
            let square_p2p1 = p2p1.square_magnitude();
            let n1p2p1 = n1.dot(&p2p1);
            let an_eps = epsilon_of(1.0);

            if p1.is_equal(&p2, CONFUSION) {
                p_bis = p1;
                *dist = 0.0;
            } else if n1p2p1 * self.sign1 < an_eps {
                valid = false;
            } else {
                p_bis = p1.translated_vec(&n1.multiplied_scalar(-(0.5 * square_p2p1 / n1p2p1)));
                *dist = p1.square_distance(&p_bis);
            }
        }

        if !valid {
            let mut d_min = INFINITE;
            let bis_pc = BisectorBisecPC::from_interval(&*self.curve2, &p1, self.sign2, v_min, v_max);
            let nor_li = Geom2dLine::new(GpAx2d::new(p1, dir2d_or_zero(n1.x(), n1.y())));
            let mut intersect = Geom2dIntGInter::new();
            intersect.perform_curves(&bis_pc, &nor_li, CONFUSION, CONFUSION);

            if intersect.is_done() && intersect.nb_points() > 0 {
                for i in 1..=intersect.nb_points() {
                    let pt = intersect.point(i);
                    if pt.param_on_second() * self.sign1 < PCONFUSION {
                        let p = *pt.value();
                        if p.square_distance(&p1) < d_min {
                            d_min = p.square_distance(&p1);
                            p_bis = p;
                            *u2 = bis_pc.link_bis_curve(pt.param_on_first());
                            *dist = d_min;
                        }
                    }
                }
            }
        }
        p_bis
    }

    /// `ValueByInt(U, U1, U2, Dist)` (`cxx:788-988`).
    pub fn value_by_int(&self, u: f64, u1: &mut f64, u2: &mut f64, dist: &mut f64) -> GpPnt2d {
        *u1 = self.link_bis_curve(u);

        let (p1, tan1) = self.curve1.d1(*u1);
        let n1 = GpVec2d::new(tan1.y(), -tan1.x());
        let mut p2 = GpPnt2d::zero();
        let mut tan2 = GpVec2d::zero();

        if p1.distance(&self.curve2.value(self.curve2.first_parameter())) < CONFUSION {
            *u2 = self.curve2.first_parameter();
            let d = self.curve2.d1(*u2);
            p2 = d.0;
            tan2 = d.1;
            if self.is_convex1 && self.is_convex2 {
                *dist = 0.0;
                return p1;
            }
            if !tan1.is_parallel(&tan2, occt_core::precision::ANGULAR) {
                *dist = 0.0;
                return p1;
            }
        }
        if p1.distance(&self.curve2.value(self.curve2.last_parameter())) < CONFUSION {
            *u2 = self.curve2.last_parameter();
            let d = self.curve2.d1(*u2);
            p2 = d.0;
            tan2 = d.1;
            if self.is_convex1 && self.is_convex2 {
                *dist = 0.0;
                return p1;
            }
            if !tan1.is_parallel(&tan2, occt_core::precision::ANGULAR) {
                *dist = 0.0;
                return p1;
            }
        }
        let _ = (p2, tan2);

        let mut ya_sol = false;
        let mut d_min = INFINITE;
        let mut p_sol = p1;
        const EPS_MAX: f64 = 1.0e-6;
        const EPS_H: f64 = 1.0e-8;
        let nb_samples = 20;
        let mut u_first_on_c2 = self.curve2.first_parameter();
        let mut u_last_on_c2 = self.curve2.last_parameter();

        if !self.my_polygon.is_empty() {
            if self.sign1 == self.sign2 {
                u_last_on_c2 = self.my_polygon.last().param_on_c2();
            } else {
                u_first_on_c2 = self.my_polygon.last().param_on_c2();
            }
        }

        if (u_last_on_c2 - u_first_on_c2).abs() < PCONFUSION / 100.0 {
            *dist = INFINITE;
            return p1;
        }

        let (eps_x, _nb) = discret_par(
            (u_last_on_c2 - u_first_on_c2).abs(),
            EPS_H,
            EPS_MAX,
            2,
            20,
        );

        let mut h = BisectorFunctionH::new(
            self.curve2.clone(),
            p1,
            GpVec2d::new(self.sign1 * self.sign2 * tan1.x(), self.sign1 * self.sign2 * tan1.y()),
        );
        let _ = nb_samples;
        let sol_root = FunctionRoots::new(
            &mut h,
            u_first_on_c2,
            u_last_on_c2,
            20,
            eps_x,
            EPS_H,
            EPS_H,
            0.0,
        );
        if sol_root.is_done() {
            for j in 1..=sol_root.nb_solutions() {
                let u_sol = sol_root.value(j);
                let p2_curve2 = self.curve2.value(u_sol);
                let p2p1 = GpVec2d::new(p1.x() - p2_curve2.x(), p1.y() - p2_curve2.y());
                let square_p2p1 = p2p1.square_magnitude();
                let n1p2p1 = n1.dot(&p2p1);

                if n1p2p1 * self.sign1 > 0.0 {
                    let p = p1.translated_vec(&n1.multiplied_scalar(-(0.5 * square_p2p1 / n1p2p1)));
                    let dist_pp1 = p1.square_distance(&p);
                    if dist_pp1 < d_min {
                        d_min = dist_pp1;
                        p_sol = p;
                        *u2 = u_sol;
                        ya_sol = true;
                    }
                }
            }
        }

        if ya_sol {
            *dist = d_min;
            let p2 = self.curve2.value(*u2);
            let pp1 = GpVec2d::new(p1.x() - p_sol.x(), p1.y() - p_sol.y());
            let pp2 = GpVec2d::new(p2.x() - p_sol.x(), p2.y() - p_sol.y());

            if pp1.dot(&pp2) > (1.0 - occt_core::precision::ANGULAR) * *dist {
                ya_sol = false;
            } else {
                if !self.is_convex1 {
                    let k1 = curvature(&*self.curve1, *u1, CONFUSION);
                    if k1 != 0.0 && *dist > 1.0 / (k1 * k1) {
                        ya_sol = false;
                    }
                }
                if ya_sol && !self.is_convex2 {
                    let k2 = curvature(&*self.curve2, *u2, CONFUSION);
                    if k2 != 0.0 && *dist > 1.0 / (k2 * k2) {
                        ya_sol = false;
                    }
                }
            }
        }
        if !ya_sol {
            *dist = INFINITE;
            p_sol = p1;
        }
        p_sol
    }

    /// `Extension(U, U1, U2, Dist, T)` (`cxx:1214-1323`).
    pub fn extension(
        &self,
        u: f64,
        u1: &mut f64,
        u2: &mut f64,
        dist: &mut f64,
        t: &mut GpVec2d,
    ) -> GpPnt2d {
        let mut p_ref = super::super::point_on_bis::BisectorPointOnBis::new();
        let mut p = GpPnt2d::zero();
        let mut p1 = GpPnt2d::zero();
        let mut p2;
        let mut p_bis;
        let mut t1;
        let mut tang = GpVec2d::zero();
        let mut d_u = 0.0;
        let mut extension_tangent = false;

        if self.my_polygon.length() == 0 {
            d_u = u - self.start_first();
            p = self.point_start;
            p1 = self.point_end;
            *u1 = self.curve1.last_parameter();
            if self.sign1 == self.sign2 {
                *u2 = self.curve2.first_parameter();
            } else {
                *u2 = self.curve2.last_parameter();
            }
            tang.set_coord(p1.x() - p.x(), p1.y() - p.y());
        } else if u < self.my_polygon.first().param_on_bis() {
            p_ref = self.my_polygon.first();
            p = p_ref.point();
            d_u = u - p_ref.param_on_bis();
            if self.extension_start {
                p1 = self.point_start;
                *u1 = self.curve1.first_parameter();
                if self.sign1 == self.sign2 {
                    *u2 = self.curve2.last_parameter();
                } else {
                    *u2 = self.curve2.first_parameter();
                }
                tang.set_coord(p.x() - p1.x(), p.y() - p1.y());
            } else {
                extension_tangent = true;
            }
        } else if u > self.my_polygon.last().param_on_bis() {
            p_ref = self.my_polygon.last();
            p = p_ref.point();
            d_u = u - p_ref.param_on_bis();
            if self.extension_end {
                p1 = self.point_end;
                *u1 = self.curve1.last_parameter();
                if self.sign1 == self.sign2 {
                    *u2 = self.curve2.last_parameter();
                } else {
                    *u2 = self.curve2.first_parameter();
                }
                tang.set_coord(p1.x() - p.x(), p1.y() - p.y());
            } else {
                extension_tangent = true;
            }
        }

        if extension_tangent {
            *u1 = p_ref.param_on_c1();
            *u2 = p_ref.param_on_c2();
            p2 = self.curve2.value(*u2);
            let d = self.curve1.d1(*u1);
            p1 = d.0;
            t1 = d.1;
            tang.set_coord(2.0 * p.x() - p1.x() - p2.x(), 2.0 * p.y() - p1.y() - p2.y());
            if tang.magnitude() < CONFUSION {
                tang = t1;
            }
            if t1.dot(&tang) < 0.0 {
                tang = tang.reversed();
            }
        }

        *t = tang.normalized().unwrap_or(GpVec2d::zero());
        p_bis = GpPnt2d::new(p.x() + d_u * t.x(), p.y() + d_u * t.y());
        *dist = p1.distance(&p_bis);
        p_bis
    }

    /// `Values(U, N, P, V1, V2, V3)` (`cxx:1081-1207`).
    pub fn values(
        &self,
        u: f64,
        n: i32,
        p: &mut GpPnt2d,
        v1: &mut GpVec2d,
        v2: &mut GpVec2d,
        v3: &mut GpVec2d,
    ) {
        v1.set_coord(0.0, 0.0);
        v2.set_coord(0.0, 0.0);
        v3.set_coord(0.0, 0.0);

        let mut u0 = 0.0;
        let mut v0 = 0.0;
        let mut dist = 0.0;

        if self.my_polygon.length() <= 1 {
            *p = self.extension(u, &mut u0, &mut v0, &mut dist, v1);
        }
        if u < self.my_polygon.first().param_on_bis() {
            *p = self.extension(u, &mut u0, &mut v0, &mut dist, v1);
            return;
        }
        if u > self.my_polygon.last().param_on_bis() {
            *p = self.extension(u, &mut u0, &mut v0, &mut dist, v1);
            return;
        }
        *p = self.value_and_dist(u, &mut u0, &mut v0, &mut dist);

        if n == 0 {
            return;
        }

        let (p1, tu, tuu) = self.curve1.d2(u0);
        let nor = GpVec2d::new(-tu.y(), tu.x());
        let nu = GpVec2d::new(-tuu.y(), tuu.x());

        let (p2, tv, tvv) = self.curve2.d2(v0);

        let pu_pv = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());

        let tu_tu = tu.dot(&tu);
        let tv_tv = tv.dot(&tv);
        let tu_tv = tu.dot(&tv);
        let tu_pupv = tu.dot(&pu_pv);
        let tv_pupv = tv.dot(&pu_pv);
        let tuu_pupv = tuu.dot(&pu_pv);
        let tu_tuu = tu.dot(&tuu);
        let tvv_pupv = tvv.dot(&pu_pv);
        let tv_tvv = tv.dot(&tvv);

        let d_hdu = 2.0
            * (tu_pupv * (tuu_pupv - tu_tu) * tv_tv + tv_pupv * tu_tv * tu_tu
                - tu_tuu * tv_pupv * tv_pupv);
        let d_hdv = 2.0
            * (tu_pupv * tu_tv * tv_tv + tv_tvv * tu_pupv * tu_pupv
                - tv_pupv * (tvv_pupv + tv_tv) * tu_tu);

        let nor_pupv = nor.dot(&pu_pv);
        let nu_pupv = nu.dot(&pu_pv);
        let nor_tv = nor.dot(&tv);

        let a = 0.5 * pu_pv.square_magnitude();
        let b = -nor_pupv;
        let bb = b * b;
        let da_du = -tu_pupv;
        let db_du = -nu_pupv;
        let da_dv = tv_pupv;
        let db_dv = -nor_tv;

        if bb < REAL_SMALL {
            let sum = tu.normalized().unwrap_or(GpVec2d::zero()).added(
                &tv.normalized().unwrap_or(GpVec2d::zero()),
            );
            *v1 = sum.multiplied_scalar(0.5 * tu.square_magnitude());
        } else {
            let df_du = tu
                .subtracted(&nor.multiplied_scalar(da_du / b - db_du * a / bb))
                .subtracted(&nu.multiplied_scalar(a / b));
            let df_dv = nor.multiplied_scalar(-da_dv / b + db_dv * a / bb);

            if d_hdv.abs() > REAL_SMALL {
                *v1 = df_du.added(&df_dv.multiplied_scalar(-d_hdu / d_hdv));
            } else {
                *v1 = tu;
            }
        }
        if n == 1 {
            return;
        }
    }

    /// `EvalD0(U)` (`cxx:992-997`).
    pub fn eval_d0(&self, u: f64) -> GpPnt2d {
        let mut u1 = 0.0;
        let mut u2 = 0.0;
        let mut dist = 0.0;
        self.value_and_dist(u, &mut u1, &mut u2, &mut dist)
    }

    /// `EvalD1(U)` (`cxx:1002-1008`).
    pub fn eval_d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let mut p = GpPnt2d::zero();
        let mut d1 = GpVec2d::zero();
        let mut v2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, 1, &mut p, &mut d1, &mut v2, &mut v3);
        (p, d1)
    }

    /// `EvalD2(U)` (`cxx:1013-1021`).
    pub fn eval_d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let mut p = GpPnt2d::zero();
        let mut d1 = GpVec2d::zero();
        let mut d2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, 2, &mut p, &mut d1, &mut d2, &mut v3);
        (p, d1, d2)
    }

    /// `EvalD3(U)` (`cxx:1026-1033`).
    pub fn eval_d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        let mut p = GpPnt2d::zero();
        let mut d1 = GpVec2d::zero();
        let mut d2 = GpVec2d::zero();
        let mut d3 = GpVec2d::zero();
        self.values(u, 3, &mut p, &mut d1, &mut d2, &mut d3);
        (p, d1, d2, d3)
    }

    /// `EvalDN(U, N)` (`cxx:1038-1059`), for `N = 1, 2, 3`.
    pub fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        let mut p = GpPnt2d::zero();
        let mut v1 = GpVec2d::zero();
        let mut v2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, n, &mut p, &mut v1, &mut v2, &mut v3);
        match n {
            1 => v1,
            2 => v2,
            3 => v3,
            _ => GpVec2d::zero(),
        }
    }
}
