//! `Convert_GridPolynomialToPoles` (the 1 x 1 patch ctor used by osculating).
//! Source: `Convert_GridPolynomialToPoles.cxx:80-291`.

use crate::bspl::banded_interp::{knot_sequence, schoenberg_points};
use crate::bspl::plib_eval::eval_poly2var;
use crate::bspl::surface_interpolate::interpolate_surface;
use crate::gp::GpPnt;

/// Result of `Convert_GridPolynomialToPoles`.
pub struct GridPoles {
    pub poles: Vec<Vec<GpPnt>>,
    pub u_knots: Vec<f64>,
    pub v_knots: Vec<f64>,
    pub u_mults: Vec<i32>,
    pub v_mults: Vec<i32>,
    pub u_degree: i32,
    pub v_degree: i32,
}

fn build_array(degree: i32, knots: &[f64], continuity: i32) -> (Vec<f64>, Vec<i32>, Vec<f64>) {
    let num_curves = knots.len() as i32 - 1;
    let multiplicities = degree - continuity;
    let mut mults = vec![multiplicities; knots.len()];
    if !mults.is_empty() {
        mults[0] = degree + 1;
        let last = knots.len() - 1;
        mults[last] = degree + 1;
    }
    let _ = num_curves;
    let flat = knot_sequence(knots, &mults, degree);
    let num_poles = flat.len() as i32 - degree - 1;
    let params = schoenberg_points(degree as usize, &flat, num_poles.max(0) as usize);
    (flat, mults, params)
}

/// `Convert_GridPolynomialToPoles` (`cxx:80-257`) for `NbUSurfaces=1`, `NbVSurfaces=1`.
pub fn grid_polynomial_to_poles(
    u_continuity: i32,
    v_continuity: i32,
    max_u_degree: i32,
    max_v_degree: i32,
    num_u: i32,
    num_v: i32,
    coefficients: &[f64],
    poly_u: &[f64],
    poly_v: &[f64],
    true_u: &[f64],
    true_v: &[f64],
) -> Option<GridPoles> {
    let real_u = max_u_degree.max(2 * u_continuity + 1);
    let real_v = max_v_degree.max(2 * v_continuity + 1);
    let u_degree = num_u - 1;
    let v_degree = num_v - 1;
    if u_degree > real_u || v_degree > real_v {
        return None;
    }
    if true_u.len() < 2 || true_v.len() < 2 || poly_u.len() < 2 || poly_v.len() < 2 {
        return None;
    }
    let (u_flat, u_mults, u_params) = build_array(u_degree, true_u, u_continuity);
    let (v_flat, v_mults, v_params) = build_array(v_degree, true_v, v_continuity);
    let siz_patch = 3 * (real_u + 1) * (real_v + 1);
    let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); v_params.len()]; u_params.len()];
    let mut patch = vec![0.0; ((u_degree + 1) * 3 * (v_degree + 1)) as usize];
    let mut u_index = 0usize;
    for (ii, &up) in u_params.iter().enumerate() {
        while u_index + 1 < true_u.len() - 1 && up > true_u[u_index + 1] {
            u_index += 1;
        }
        let denom_u = true_u[u_index + 1] - true_u[u_index];
        let n_u = if denom_u.abs() < 1e-30 {
            0.0
        } else {
            (up - true_u[u_index]) / denom_u
        };
        let u_value = (1.0 - n_u) * poly_u[0] + n_u * poly_u[1];
        let mut v_index = 0usize;
        for (jj, &vp) in v_params.iter().enumerate() {
            while v_index + 1 < true_v.len() - 1 && vp > true_v[v_index + 1] {
                v_index += 1;
            }
            let denom_v = true_v[v_index + 1] - true_v[v_index];
            let n_v = if denom_v.abs() < 1e-30 {
                0.0
            } else {
                (vp - true_v[v_index]) / denom_v
            };
            let v_value = (1.0 - n_v) * poly_v[0] + n_v * poly_v[1];
            let patch_indice = u_index + (true_u.len() - 1) * v_index;
            let mut ll = 0usize;
            for k1 in 1..=num_u {
                let mut pos = siz_patch * (patch_indice as i32) + 3 * (real_v + 1) * (k1 - 1);
                for _k2 in 1..=num_v {
                    patch[ll] = coefficients.get(pos as usize).copied().unwrap_or(0.0);
                    patch[ll + 1] = coefficients.get((pos + 1) as usize).copied().unwrap_or(0.0);
                    patch[ll + 2] = coefficients.get((pos + 2) as usize).copied().unwrap_or(0.0);
                    ll += 3;
                    pos += 3;
                }
            }
            let mut digit = [0.0; 3];
            eval_poly2var(
                u_value,
                v_value,
                num_u - 1,
                num_v - 1,
                3,
                &patch,
                &mut digit,
            );
            poles[ii][jj] = GpPnt::new(digit[0], digit[1], digit[2]);
        }
    }
    if !interpolate_surface(
        u_degree, v_degree, &u_flat, &v_flat, &u_params, &v_params, &mut poles,
    ) {
        return None;
    }
    Some(GridPoles {
        poles,
        u_knots: true_u.to_vec(),
        v_knots: true_v.to_vec(),
        u_mults,
        v_mults,
        u_degree,
        v_degree,
    })
}
