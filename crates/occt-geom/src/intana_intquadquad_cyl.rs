//! Cylinder-quadric `IntAna_IntQuadQuad::Perform`.

use occt_core::gp::GpCylinder;

use super::{add_special_points_cyl, IntQuadQuad, TrigFn, TrigPolyRoots, MAX_CURVES};
use crate::intana::IntAnaQuadric;

pub(super) fn perform(iqq: &mut IntQuadQuad, cyl: &GpCylinder, quad: &IntAnaQuadric, _tol: f64) {
    iqq.done = true;
    iqq.identical = false;
    iqq.nb_curves = 0;
    iqq.nb_points = 0;
    iqq.reset_links();

    let un_seul = false;
    let deux = true;
    let z_pos = true;
    let z_indiff = true;
    let z_neg = false;
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
    let r_cyl = cyl.radius();
    let real_eps = f64::EPSILON;
    let eps = iqq.epsilon;

    let q = quad.new_coefficients(&cyl.position());
    let [qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1] = q;

    if qzz.abs() < iqq.epsilon_poly {
        iqq.done = false;
        return;
    }

    let r2 = r_cyl * r_cyl;
    let c_1 = qz * qz - qzz * q1;
    let c_ss = r2 * (qyz * qyz - qyy * qzz);
    let c_cc = r2 * (qxz * qxz - qxx * qzz);
    let c_s = r_cyl * (qyz * qz - qy * qzz);
    let c_c = r_cyl * (qxz * qz - qx * qzz);
    let c_sc = r2 * (qxz * qyz - qxy * qzz);
    let mtf = TrigFn::new(c_cc, c_ss, c_sc, c_c, c_s, c_1);
    let pol = TrigPolyRoots::new(c_cc - c_ss, c_sc, c_c + c_c, c_s + c_s, c_1 + c_ss, 0.0, two_pi);
    if !pol.is_done() {
        iqq.done = false;
        return;
    }

    let push2 = |iqq: &mut IntQuadQuad, t1: f64, t2: f64, two_z: bool, z_a: bool, z_b: bool| {
        let n = iqq.nb_curves;
        if n >= MAX_CURVES {
            return;
        }
        iqq.curves[n].set_cylinder_quad_values(
            cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, t1, t2, two_z, z_a,
        );
        iqq.nb_curves += 1;
        if !two_z {
            let n = iqq.nb_curves;
            if n >= MAX_CURVES {
                return;
            }
            iqq.curves[n].set_cylinder_quad_values(
                cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, t1, t2, two_z, z_b,
            );
            iqq.nb_curves += 1;
        }
    };

    if pol.infinite_roots() {
        push2(iqq, 0.0, two_pi, un_seul, z_pos, z_neg);
        return;
    }

    let nbsol = pol.nb_solutions();
    if nbsol == 0 {
        if mtf.value(std::f64::consts::PI) >= -real_eps {
            push2(iqq, 0.0, two_pi, un_seul, z_pos, z_neg);
        }
        return;
    }

    if nbsol == 1 {
        if mtf.value(pol.value(1) + std::f64::consts::PI) >= -real_eps {
            push2(iqq, 0.0, two_pi, un_seul, z_pos, z_neg);
        }
        return;
    }

    let mut un_pt_tg = false;
    if nbsol == 2 {
        for i in 1..=nbsol {
            let theta1 = pol.value(i);
            let theta2 = if i < nbsol {
                pol.value(i + 1)
            } else {
                pol.value(1) + two_pi
            };
            if (theta2 - theta1).abs() <= real_eps {
                un_pt_tg = true;
                let mut autre = theta1 - 0.1;
                if autre < 0.0 {
                    autre = theta1 + 0.1;
                }
                if mtf.value(autre) >= 0.0 {
                    let mut t1 = theta1;
                    let mut t2 = theta1 + two_pi;
                    add_special_points_cyl(quad, cyl, &mut t1, &mut t2);
                    push2(iqq, t1, t2, un_seul, z_pos, z_neg);
                }
            }
        }
    }

    if un_pt_tg {
        return;
    }
    for i in 1..=nbsol {
        let theta1 = pol.value(i);
        let mut theta2 = if i < nbsol {
            pol.value(i + 1)
        } else {
            pol.value(1) + two_pi
        };
        if (theta2 - theta1).abs() <= 1e-12 {
            continue;
        }
        let qwet = mtf.value(0.5 * (theta1 + theta2))
            + mtf.value(0.4 * theta1 + 0.6 * theta2)
            + mtf.value(0.6 * theta1 + 0.4 * theta2);
        if qwet < 0.0 {
            continue;
        }
        let theta3 = if i + 1 < nbsol {
            pol.value(i + 2)
        } else {
            pol.value(1) + two_pi
        };
        let close_double = (theta3 - theta2) < 5.0e-8;
        let mut t1 = theta1;
        add_special_points_cyl(quad, cyl, &mut t1, &mut theta2);
        if close_double {
            push2(iqq, t1, theta2, un_seul, z_pos, z_neg);
        } else {
            let n = iqq.nb_curves;
            if n < MAX_CURVES {
                iqq.curves[n].set_cylinder_quad_values(
                    cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1, eps, t1, theta2, deux,
                    z_indiff,
                );
                iqq.nb_curves += 1;
            }
        }
    }
}
