//! Adaptive U1 step from the linearized CylCyl system.
//! Source: `VBoundaryPrecise` / `DeltaU1Computing` / `StepComputing`.

use super::NUL_VALUE;

/// 3x5 matrix, 0-based, matching OCCT `math_Matrix(1,3,1,5)`.
pub(crate) type Mat35 = [[f64; 5]; 3];

fn col(m: &Mat35, j: usize) -> [f64; 3] {
    [m[0][j], m[1][j], m[2][j]]
}

fn det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn sub_scaled(a: [f64; 3], s: f64, b: [f64; 3]) -> [f64; 3] {
    [a[0] - s * b[0], a[1] - s * b[1], a[2] - s * b[2]]
}

/// `VBoundaryPrecise`: pick V1/V2 decrease when dV/dU1 < 0.
fn v_boundary_precise(m: &Mat35, v1_decr: f64, v2_decr: f64, v1_set: &mut f64, v2_set: &mut f64) {
    let c1 = col(m, 0);
    let c2 = col(m, 1);
    let du1 = col(m, 2);
    let du2 = col(m, 3);
    let det = det3(c1, c2, du2);
    let det1 = det3(du1, c2, du2);
    let det2 = det3(c1, du1, du2);
    if det * det1 > 0.0 {
        *v1_set = v1_decr;
    }
    if det * det2 > 0.0 {
        *v2_set = v2_decr;
    }
}

/// `DeltaU1Computing`. `syst` columns are (V-known partner, dU1, dU2).
fn delta_u1_computing(syst: [[f64; 3]; 3], free: [f64; 3]) -> Option<f64> {
    let det = det3(syst[0], syst[1], syst[2]);
    if det.abs() > NUL_VALUE {
        let det1 = det3(syst[0], free, syst[2]);
        return Some(det1.abs() / det.abs());
    }
    None
}

/// `StepComputing`. Returns `|dU1|` providing about `delta_v1` / `delta_v2`.
pub(crate) fn step_computing(
    m: &Mat35,
    v1_cur: f64,
    v2_cur: f64,
    delta_v1: f64,
    delta_v2: f64,
) -> Option<f64> {
    let mut found = f64::MAX;
    let mut v1_set = v1_cur + delta_v1;
    let mut v2_set = v2_cur + delta_v2;
    v_boundary_precise(m, v1_cur - delta_v1, v2_cur - delta_v2, &mut v1_set, &mut v2_set);

    let du1 = col(m, 2);
    let du2 = col(m, 3);
    let rhs = col(m, 4);
    let mut ok = false;
    for i in 0..2 {
        let (c0, free) = if i == 0 {
            (col(m, 1), sub_scaled(rhs, v1_set, col(m, 0)))
        } else {
            (col(m, 0), sub_scaled(rhs, v2_set, col(m, 1)))
        };
        let syst = [c0, du1, du2];
        if let Some(du) = delta_u1_computing(syst, free) {
            ok = true;
            if du < found {
                found = du;
            }
        }
    }
    if ok {
        return Some(found);
    }

    let free = sub_scaled(sub_scaled(rhs, v1_set, col(m, 0)), v2_set, col(m, 1));
    let det1 = m[0][2] * m[1][3] - m[1][2] * m[0][3];
    let det2 = m[0][2] * m[2][3] - m[2][2] * m[0][3];
    let det3 = m[1][2] * m[2][3] - m[2][2] * m[1][3];
    let a1 = det1.abs();
    let a2 = det2.abs();
    let a3 = det3.abs();
    if a1 >= a2 {
        if a1 >= a3 {
            if a1 <= NUL_VALUE {
                return None;
            }
            Some((free[0] * m[1][3] - free[1] * m[0][3]).abs() / a1)
        } else {
            if a3 <= NUL_VALUE {
                return None;
            }
            Some((free[1] * m[2][3] - free[2] * m[1][3]).abs() / a3)
        }
    } else if a2 >= a3 {
        if a2 <= NUL_VALUE {
            return None;
        }
        Some((free[0] * m[2][3] - free[2] * m[0][3]).abs() / a2)
    } else {
        if a3 <= NUL_VALUE {
            return None;
        }
        Some((free[1] * m[2][3] - free[2] * m[1][3]).abs() / a3)
    }
}

/// Fill the 3x5 Jacobian used by `StepComputing` (walker ~7403).
pub(crate) fn fill_step_matrix(c: &super::Coeffs, u1: f64, u2: f64) -> Mat35 {
    let s1 = u1.sin();
    let c1 = u1.cos();
    let s2 = u2.sin();
    let c2 = u2.cos();
    let mut m = [[0.0; 5]; 3];
    for r in 0..3 {
        let a1 = xyz(c.vec_a1, r);
        let a2 = xyz(c.vec_a2, r);
        let b1 = xyz(c.vec_b1, r);
        let b2 = xyz(c.vec_b2, r);
        let cc1 = xyz(c.vec_c1, r);
        let cc2 = xyz(c.vec_c2, r);
        let d = xyz(c.vec_d, r);
        m[r][0] = cc1;
        m[r][1] = cc2;
        m[r][2] = a1 * s1 - b1 * c1;
        m[r][3] = a2 * s2 - b2 * c2;
        m[r][4] = a1 * c1 + b1 * s1 + a2 * c2 + b2 * s2 + d;
    }
    m
}

fn xyz(v: occt_core::gp::GpXyz, i: usize) -> f64 {
    match i {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}
