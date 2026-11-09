//! `BSplSLib::Interpolate` non-rational.
//! Source: `BSplSLib.cxx:3955-4038`.

use crate::bspl::banded_interp::interpolate;
use crate::gp::GpPnt;

/// `BSplSLib::Interpolate` (`cxx:3955-4038`). `poles[i][j]` is `Poles(i+1, j+1)`.
/// Returns `true` when `InversionProblem == 0`.
pub fn interpolate_surface(
    u_degree: i32,
    v_degree: i32,
    u_flat: &[f64],
    v_flat: &[f64],
    u_params: &[f64],
    v_params: &[f64],
    poles: &mut [Vec<GpPnt>],
) -> bool {
    let u_len = u_params.len();
    let v_len = v_params.len();
    if u_len == 0 || v_len == 0 {
        return false;
    }
    let dim_u = 3 * u_len;
    let mut points = vec![0.0; v_len * dim_u];
    for ii in 0..v_len {
        for jj in 0..u_len {
            let p = poles[jj][ii];
            let ll = jj * 3;
            points[ii * dim_u + ll] = p.x();
            points[ii * dim_u + ll + 1] = p.y();
            points[ii * dim_u + ll + 2] = p.z();
        }
    }
    if interpolate(v_degree as usize, v_flat, v_params, &mut points, dim_u).is_err() {
        return false;
    }
    let dim_v = 3 * v_len;
    let mut iso = vec![0.0; u_len * dim_v];
    for ii in 0..u_len {
        let kk = ii * 3;
        for jj in 0..v_len {
            let ll = jj * 3;
            iso[ii * dim_v + ll] = points[jj * dim_u + kk];
            iso[ii * dim_v + ll + 1] = points[jj * dim_u + kk + 1];
            iso[ii * dim_v + ll + 2] = points[jj * dim_u + kk + 2];
        }
    }
    if interpolate(u_degree as usize, u_flat, u_params, &mut iso, dim_v).is_err() {
        return false;
    }
    for ii in 0..u_len {
        for jj in 0..v_len {
            let ll = jj * 3;
            poles[ii][jj] = GpPnt::new(
                iso[ii * dim_v + ll],
                iso[ii * dim_v + ll + 1],
                iso[ii * dim_v + ll + 2],
            );
        }
    }
    true
}
