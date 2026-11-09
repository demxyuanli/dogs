//! `PrepareConvert` (`AdvApprox_ApproxAFunction.cxx:123-311`).

use super::prelude::*;

/// `AdvApprox_ApproxAFunction::PrepareConvert` (`cxx:123-311`).
///
/// Determines, for every interior node, the local continuity order actually
/// achieved by the piecewise polynomial result, using the same subspace
/// layout (`Num1DSS`/`Num2DSS`/`Num3DSS`) as `Perform`. The error bounds are
/// accumulated into `error_max` exactly as OCCT does. The returned array is
/// indexed like `Convert_CompPolynomialToPoles::Continuity`: entry `i` is the
/// continuity at the knot before span `i` (OCCT's 1-based `Continuity(i+1)`),
/// and entry 0 is never read by the converter.
#[allow(clippy::too_many_arguments)]
pub(crate) fn prepare_convert(
    num_curves: i32,
    max_degree: i32,
    continuity_order: i32,
    num_1d_ss: i32,
    num_2d_ss: i32,
    num_3d_ss: i32,
    num_coeff_per_curve: &[i32],
    coefficients: &[f64],
    poly_intervals: &[[f64; 2]],
    true_intervals: &[f64],
    local_tolerance: &[f64],
    error_max: &mut [f64],
) -> Vec<i32> {
    let dimension = (num_1d_ss + 2 * num_2d_ss + 3 * num_3d_ss) as usize;
    let nb_space = (num_1d_ss + num_2d_ss + num_3d_ss) as usize;
    let real_degree = (max_degree + 1).max(2 * continuity_order + 2) as usize;
    let mut tab_continuity = vec![0i32; num_curves.max(0) as usize];
    if continuity_order == 0 {
        return tab_continuity;
    }
    let block = (continuity_order + 1) as usize * dimension;
    let mut res = vec![0.0; 2 * block];
    let mut prec = vec![0.0; nb_space];
    let mut suivant = vec![0.0; nb_space];
    for icurve in 1..num_curves.max(0) as usize {
        let mut is_ci = true;
        let deg1 = num_coeff_per_curve[icurve - 1] - 1;
        let deg2 = num_coeff_per_curve[icurve] - 1;
        let coef1 = (icurve - 1) * dimension * real_degree;
        let coef2 = icurve * dimension * real_degree;
        eval_polynomial(
            poly_intervals[icurve - 1][1],
            continuity_order,
            deg1,
            dimension,
            &coefficients[coef1..],
            &mut res[..block],
        );
        eval_polynomial(
            poly_intervals[icurve][0],
            continuity_order,
            deg2,
            dimension,
            &coefficients[coef2..],
            &mut res[block..],
        );
        for iordre in 1..=continuity_order {
            if !is_ci {
                break;
            }
            const TOLER: f64 = 1.0e-5;
            let f1_dividend = poly_intervals[icurve - 1][1] - poly_intervals[icurve - 1][0];
            let f2_dividend = poly_intervals[icurve][1] - poly_intervals[icurve][0];
            let f1_divizor = true_intervals[icurve] - true_intervals[icurve - 1];
            let f2_divizor = true_intervals[icurve + 1] - true_intervals[icurve];
            let facteur1 = if f1_divizor.abs() < TOLER {
                0.0
            } else {
                (f1_dividend / f1_divizor).powi(iordre)
            };
            let facteur2 = if f2_divizor.abs() < TOLER {
                0.0
            } else {
                (f2_dividend / f2_divizor).powi(iordre)
            };
            let normal1 = f1_divizor.powi(iordre);
            let normal2 = f2_divizor.powi(iordre);
            let off = iordre as usize * dimension;
            let mut idim = 0usize;
            // 1D subspaces.
            for ii in 0..num_1d_ss.max(0) as usize {
                let v1 = res[off + ii] * facteur1;
                let v2 = res[block + off + ii] * facteur2;
                let eps = local_tolerance[idim] * 0.01;
                let diff = (v1 - v2).abs();
                let moy = (v1 + v2).abs();
                if diff > moy * 1.0e-9 {
                    prec[idim] = diff * normal1;
                    suivant[idim] = diff * normal2;
                    if prec[idim] > eps || suivant[idim] > eps {
                        is_ci = false;
                    }
                } else {
                    prec[idim] = 0.0;
                    suivant[idim] = 0.0;
                }
                idim += 1;
            }
            // 2D subspaces.
            for ii in 0..num_2d_ss.max(0) as usize {
                let idx = off + num_1d_ss.max(0) as usize + 2 * ii;
                let v1 = [res[idx] * facteur1, res[idx + 1] * facteur1];
                let v2 = [
                    res[block + idx] * facteur2,
                    res[block + idx + 1] * facteur2,
                ];
                let eps = local_tolerance[idim] * 0.01;
                let diff = (v1[0] - v2[0]).abs() + (v1[1] - v2[1]).abs();
                let moy = (v1[0] + v2[0]).abs() + (v1[1] + v2[1]).abs();
                if diff > moy * 1.0e-9 {
                    prec[idim] = diff * normal1;
                    suivant[idim] = diff * normal2;
                    if prec[idim] > eps || suivant[idim] > eps {
                        is_ci = false;
                    }
                } else {
                    prec[idim] = 0.0;
                    suivant[idim] = 0.0;
                }
                idim += 1;
            }
            // 3D subspaces.
            for ii in 0..num_3d_ss.max(0) as usize {
                let idx = off
                    + num_1d_ss.max(0) as usize
                    + 2 * num_2d_ss.max(0) as usize
                    + 3 * ii;
                let v1 = [
                    res[idx] * facteur1,
                    res[idx + 1] * facteur1,
                    res[idx + 2] * facteur1,
                ];
                let v2 = [
                    res[block + idx] * facteur2,
                    res[block + idx + 1] * facteur2,
                    res[block + idx + 2] * facteur2,
                ];
                let eps = local_tolerance[idim] * 0.01;
                let diff = (v1[0] - v2[0]).abs() + (v1[1] - v2[1]).abs() + (v1[2] - v2[2]).abs();
                let moy = (v1[0] + v2[0]).abs() + (v1[1] + v2[1]).abs() + (v1[2] + v2[2]).abs();
                if diff > moy * 1.0e-9 {
                    prec[idim] = diff * normal1;
                    suivant[idim] = diff * normal2;
                    if prec[idim] > eps || suivant[idim] > eps {
                        is_ci = false;
                    }
                } else {
                    prec[idim] = 0.0;
                    suivant[idim] = 0.0;
                }
                idim += 1;
            }
            if is_ci {
                tab_continuity[icurve] = iordre;
                let base = (icurve - 1) * nb_space;
                for idim in 0..nb_space {
                    error_max[base + idim] += prec[idim];
                    error_max[base + nb_space + idim] += suivant[idim];
                }
            }
        }
    }
    tab_continuity
}
