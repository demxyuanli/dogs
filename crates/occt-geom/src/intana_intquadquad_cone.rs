//! Cone-quadric `IntAna_IntQuadQuad::Perform`.

use occt_core::gp::GpCone;

use super::{IntQuadQuad, TrigFn, TrigPolyRoots, MAX_CURVES};
use crate::intana::IntAnaQuadric;

pub(super) fn perform(iqq: &mut IntQuadQuad, cone: &GpCone, quad: &IntAnaQuadric, _tol: f64) {
    let un_seul = false;
    let z_pos = true;
    let z_neg = false;
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
    iqq.done = true;
    iqq.identical = false;
    iqq.nb_curves = 0;
    iqq.nb_points = 0;
    iqq.reset_links();

    let q = {
        let mut ax = cone.position();
        ax.set_location(cone.apex());
        quad.new_coefficients(&ax)
    };
    let [qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1] = q;
    let tg = 1.0 / cone.semi_angle().tan();
    let eps = iqq.epsilon;

    let z2_cc = qxx;
    let z2_ss = qyy;
    let z2_cte = qzz * tg * tg;
    let z2_sc = qxy;
    let z2_c = tg * qxz;
    let z2_s = tg * qyz;
    let pol_z2 = TrigPolyRoots::new(
        z2_cc - z2_ss,
        z2_sc,
        z2_c + z2_c,
        z2_s + z2_s,
        z2_cte + z2_ss,
        0.0,
        two_pi,
    );
    if !pol_z2.is_done() {
        iqq.done = false;
        return;
    }
    let nbsol_z2 = pol_z2.nb_solutions();

    let z1_cte = 2.0 * tg * qz;
    let z1_s = qy;
    let z1_c = qx;
    let pol_z1 = TrigPolyRoots::new(0.0, 0.0, z1_c + z1_c, z1_s + z1_s, z1_cte, 0.0, two_pi);
    if !pol_z1.is_done() {
        iqq.done = false;
        return;
    }
    let mtf_z1 = TrigFn::new(0.0, 0.0, 0.0, z1_c, z1_s, z1_cte);
    let nbsol1 = pol_z1.nb_solutions();

    if pol_z2.infinite_roots() {
        if !pol_z1.infinite_roots() {
            if nbsol1 == 0 {
                iqq.curves[0].set_cone_quad_values(
                    cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, 0.0, two_pi, un_seul,
                    z_pos,
                );
                iqq.nb_curves = 1;
            }
        } else if q1.abs() <= iqq.epsilon {
            iqq.done = false;
        }
        return;
    }

    let c_1 = tg * tg * (qz * qz - qzz * q1);
    let c_ss = qy * qy - qyy * q1;
    let c_cc = qx * qx - qxx * q1;
    let c_s = tg * (qy * qz - qyz * q1);
    let c_c = tg * (qx * qz - qxz * q1);
    let c_sc = qx * qy - qxy * q1;
    let pol = TrigPolyRoots::new(c_cc - c_ss, c_sc, c_c + c_c, c_s + c_s, c_1 + c_ss, 0.0, two_pi);
    if !pol.is_done() {
        iqq.done = false;
        return;
    }
    let mut nbsol = pol.nb_solutions();
    let mtf = TrigFn::new(c_cc, c_ss, c_sc, c_c, c_s, c_1);

    if pol.infinite_roots() {
        iqq.curves[0].set_cone_quad_values(
            cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, 0.0, two_pi, un_seul, z_pos,
        );
        iqq.curves[1].set_cone_quad_values(
            cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, 0.0, two_pi, un_seul, z_neg,
        );
        iqq.nb_curves = 2;
        return;
    }

    if nbsol == 0 && mtf.value(std::f64::consts::PI) < 0.0 {
        return;
    }

    let disc_const_pos = nbsol == 0;
    if nbsol == 0 {
        nbsol = 1;
    }

    for i in 1..=nbsol {
        let (theta1, theta2) = if disc_const_pos {
            (0.0, two_pi - iqq.epsilon)
        } else {
            let t1 = pol.value(i);
            let t2 = if i < nbsol {
                pol.value(i + 1)
            } else {
                pol.value(1) + two_pi
            };
            (t1, t2)
        };
        if (theta2 - theta1).abs() <= iqq.epsilon {
            iqq.done = false;
            return;
        }
        let qwet = mtf.value(0.5 * (theta1 + theta2))
            + mtf.value(0.4 * theta1 + 0.6 * theta2)
            + mtf.value(0.6 * theta1 + 0.4 * theta2);
        if qwet < 0.0 {
            continue;
        }

        let mut z2_in = false;
        for i2 in 1..=nbsol_z2 {
            let r = pol_z2.value(i2);
            if r > theta1 && r < theta2 {
                z2_in = true;
            } else {
                let r2 = r + two_pi;
                if r2 > theta1 && r2 < theta2 {
                    z2_in = true;
                }
            }
        }

        if !z2_in {
            let n = iqq.nb_curves;
            if n + 1 < MAX_CURVES {
                iqq.curves[n].set_cone_quad_values(
                    cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, theta1, theta2, un_seul,
                    z_pos,
                );
                iqq.nb_curves += 1;
                let n = iqq.nb_curves;
                iqq.curves[n].set_cone_quad_values(
                    cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, theta1, theta2, un_seul,
                    z_neg,
                );
                iqq.nb_curves += 1;
            }
            continue;
        }

        let mut new_min = theta1;
        let mut no_changes = true;
        for i2 in 1..=(nbsol_z2 + nbsol_z2) {
            let to = if i2 > nbsol_z2 {
                pol_z2.value(i2 - nbsol_z2) + two_pi
            } else {
                pol_z2.value(i2)
            };
            if to < theta2 && to > new_min {
                no_changes = false;
                let n_neg = iqq.nb_curves;
                if n_neg + 1 >= MAX_CURVES {
                    break;
                }
                iqq.curves[n_neg].set_cone_quad_values(
                    cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, new_min, to, un_seul,
                    z_neg,
                );
                iqq.nb_curves += 1;
                let n_pos = iqq.nb_curves;
                iqq.curves[n_pos].set_cone_quad_values(
                    cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, new_min, to, un_seul,
                    z_pos,
                );
                if pol_z2.is_a_root(new_min) {
                    if mtf_z1.value(new_min) < 0.0 {
                        iqq.curves[n_pos].set_is_first_open(true);
                    } else {
                        iqq.curves[n_neg].set_is_first_open(true);
                    }
                }
                if mtf_z1.value(to) < 0.0 {
                    iqq.curves[n_pos].set_is_last_open(true);
                } else {
                    iqq.curves[n_neg].set_is_last_open(true);
                }
                iqq.nb_curves += 1;
                new_min = to;
            }
        }

        let n_neg = iqq.nb_curves;
        if n_neg + 1 >= MAX_CURVES {
            continue;
        }
        let a0 = if no_changes { theta1 } else { new_min };
        iqq.curves[n_neg].set_cone_quad_values(
            cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, a0, theta2, un_seul, z_neg,
        );
        iqq.nb_curves += 1;
        let n_pos = iqq.nb_curves;
        iqq.curves[n_pos].set_cone_quad_values(
            cone, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, a0, theta2, un_seul, z_pos,
        );
        let root_first = if no_changes { theta1 } else { new_min };
        if pol_z2.is_a_root(root_first) {
            if mtf_z1.value(root_first) < 0.0 {
                iqq.curves[n_pos].set_is_first_open(true);
            } else {
                iqq.curves[n_neg].set_is_first_open(true);
            }
        }
        if pol_z2.is_a_root(theta2) {
            if mtf_z1.value(theta2) < 0.0 {
                iqq.curves[n_pos].set_is_last_open(true);
            } else {
                iqq.curves[n_neg].set_is_last_open(true);
            }
        }
        iqq.nb_curves += 1;
    }
    iqq.internal_set_next_and_previous();
}
