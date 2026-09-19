//! `AdvApprox_ApproxAFunction` for a single 3D subspace.
//! Source: `AdvApprox_ApproxAFunction.cxx` (constructor `cxx:603-630`,
//! `Perform` `cxx:663-956`, `Approximation` `cxx:364-599`).

use occt_core::bspl::comp_poly::CompPolynomialToPoles;
use occt_core::bspl::plib_eval::eval_polynomial;
use occt_core::bspl::plib_jacobi::{jacobi_parameters, JacobiPolynomial};
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;

use super::simple::{EvalFn, SimpleApprox};

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
fn prepare_convert(
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

/// Result of a 3D `ApproxAFunction` run.
pub struct ApproxAFunction3d {
    pub done: bool,
    pub has_result: bool,
    pub poles: Vec<GpPnt>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
}

/// Result of `AdvApprox_ApproxAFunction` with `Num1DSS=2` (`Approx_SameParameter`).
pub struct ApproxAFunction1dPair {
    pub done: bool,
    pub has_result: bool,
    pub poles_u: Vec<f64>,
    pub poles_v: Vec<f64>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
}

impl ApproxAFunction3d {
    /// `AdvApprox_ApproxAFunction` with default `DichoCutting`.
    /// Num1DSS=0, Num2DSS=0, Num3DSS=1, Continuity C1.
    pub fn approx_c1(
        first: f64,
        last: f64,
        max_deg: i32,
        max_seg: i32,
        tol: f64,
        eval: EvalFn<'_>,
    ) -> Result<Self, ()> {
        if last < first || max_deg < 1 || max_seg < 0 {
            return Err(());
        }
        let max_deg = max_deg.min(14);
        let continuity_order = 1i32;
        let num_max_coeffs = (max_deg + 1).max(2 * continuity_order + 2);
        let max_degree = num_max_coeffs - 1;
        let (nb_gauss, work_degree) = jacobi_parameters(Shape::C1, max_degree, 1)?;
        let jac = JacobiPolynomial::new(work_degree, Shape::C1)?;
        let mut approx = SimpleApprox::new(3, 1, Shape::C1, work_degree, nb_gauss, jac.clone())?;
        let local_dim = [3i32];
        let local_tol = [tol];
        let mut intervals = vec![0.0; (max_seg as usize) + 1];
        intervals[0] = first;
        intervals[1] = last;
        let mut nupil = 1i32;
        let mut num_curves = 0i32;
        let mut is_cut = false;
        let mut error_code = 0i32;
        let mut num_coeff = vec![0i32; max_seg as usize];
        let mut coeff = vec![0.0; (max_seg * num_max_coeffs * 3) as usize];
        let mut err_max = vec![0.0; max_seg as usize];
        let mut err_avg = vec![0.0; max_seg as usize];
        if max_seg < 1 || (last - first).abs() < 1.0e-9 {
            return Err(());
        }
        while nupil - num_curves != 0 {
            approx.perform(
                &local_dim,
                &local_tol,
                intervals[num_curves as usize],
                intervals[num_curves as usize + 1],
                max_degree,
                eval,
            );
            if !approx.is_done() {
                error_code = 1;
                break;
            }
            let ok_tol = approx.max_error(0) <= tol;
            if ok_tol {
                num_curves += 1;
            } else {
                let a = intervals[num_curves as usize];
                let b = intervals[num_curves as usize + 1];
                let large = (b - a).abs() >= 20.0 * PCONFUSION;
                let tmil = 0.5 * (a + b);
                if nupil < max_seg && large {
                    is_cut = true;
                    let from = num_curves as usize + 1;
                    for i in (from..=nupil as usize).rev() {
                        intervals[i + 1] = intervals[i];
                    }
                    intervals[from] = tmil;
                    nupil += 1;
                    continue;
                }
                num_curves += 1;
            }
            err_max[(num_curves - 1) as usize] = approx.max_error(0);
            err_avg[(num_curves - 1) as usize] = approx.average_error(0);
            let mut the_deg = approx.degree();
            if is_cut && the_deg < 2 * continuity_order + 1 {
                the_deg = 2 * continuity_order + 1;
            }
            num_coeff[(num_curves - 1) as usize] = the_deg + 1;
            let canon = jac.to_coefficients(3, the_deg, approx.coefficients());
            let f = ((the_deg + 1) * 3) as usize;
            let dest = ((num_curves - 1) * 3 * num_max_coeffs) as usize;
            for i in 0..f.min(canon.len()) {
                if dest + i < coeff.len() {
                    coeff[dest + i] = canon[i];
                }
            }
        }
        let _ = (&err_max, &err_avg);
        if error_code != 0 && error_code != -1 {
            return Ok(Self {
                done: false,
                has_result: false,
                poles: Vec::new(),
                knots: Vec::new(),
                mults: Vec::new(),
                degree: 0,
            });
        }
        for i in 0..num_curves as usize {
            num_coeff[i] = num_coeff[i].max(2);
        }
        let mut poly_iv = vec![[-1.0, 1.0]; num_curves as usize];
        for i in 0..num_curves as usize {
            poly_iv[i] = [-1.0, 1.0];
        }
        let true_iv = intervals[..=num_curves as usize].to_vec();
        let continuity = prepare_convert(
            num_curves,
            max_degree,
            continuity_order,
            0,
            0,
            1,
            &num_coeff[..num_curves as usize],
            &coeff,
            &poly_iv,
            &true_iv,
            &local_tol,
            &mut err_max,
        );
        let conv = CompPolynomialToPoles::new(
            num_curves,
            3,
            max_degree,
            &continuity,
            &num_coeff[..num_curves as usize],
            &coeff,
            &poly_iv,
            &true_iv,
        )?;
        if !conv.done {
            return Err(());
        }
        let poles = conv
            .poles
            .iter()
            .map(|p| GpPnt::new(p[0], p[1], p[2]))
            .collect();
        let done = error_code == 0;
        Ok(Self {
            done,
            has_result: true,
            poles,
            knots: conv.knots,
            mults: conv.mults,
            degree: conv.degree,
        })
    }
}

impl ApproxAFunction1dPair {
    /// `AdvApprox_ApproxAFunction(Num1DSS=2, Num2DSS=0, Num3DSS=0)` as used by
    /// `Approx_SameParameter` (`cxx:434-445`). Continuity C0 or C1.
    pub fn approx(
        first: f64,
        last: f64,
        continuity: Shape,
        max_deg: i32,
        max_seg: i32,
        tol_u: f64,
        tol_v: f64,
        eval: EvalFn<'_>,
    ) -> Result<Self, ()> {
        if last < first || max_deg < 1 || max_seg < 0 {
            return Err(());
        }
        let max_deg = max_deg.min(14);
        let continuity_order = match continuity {
            Shape::C0 => 0i32,
            Shape::C1 => 1i32,
            Shape::C2 => 2i32,
            _ => 1i32,
        };
        let num_max_coeffs = (max_deg + 1).max(2 * continuity_order + 2);
        let max_degree = num_max_coeffs - 1;
        let (nb_gauss, work_degree) = jacobi_parameters(continuity, max_degree, 1)?;
        let jac = JacobiPolynomial::new(work_degree, continuity)?;
        let mut approx = SimpleApprox::new(2, 2, continuity, work_degree, nb_gauss, jac.clone())?;
        let local_dim = [1i32, 1i32];
        let local_tol = [tol_u, tol_v];
        let mut intervals = vec![0.0; (max_seg as usize) + 1];
        intervals[0] = first;
        intervals[1] = last;
        let mut nupil = 1i32;
        let mut num_curves = 0i32;
        let mut is_cut = false;
        let mut error_code = 0i32;
        let mut num_coeff = vec![0i32; max_seg as usize];
        let mut coeff = vec![0.0; (max_seg * num_max_coeffs * 2) as usize];
        let mut err_max = vec![0.0; (max_seg * 2) as usize];
        let mut err_avg = vec![0.0; (max_seg * 2) as usize];
        if max_seg < 1 || (last - first).abs() < 1.0e-9 {
            return Err(());
        }
        while nupil - num_curves != 0 {
            approx.perform(
                &local_dim,
                &local_tol,
                intervals[num_curves as usize],
                intervals[num_curves as usize + 1],
                max_degree,
                eval,
            );
            if !approx.is_done() {
                error_code = 1;
                break;
            }
            let ok_tol = approx.max_error(0) <= tol_u && approx.max_error(1) <= tol_v;
            if ok_tol {
                num_curves += 1;
            } else {
                let a = intervals[num_curves as usize];
                let b = intervals[num_curves as usize + 1];
                let large = (b - a).abs() >= 20.0 * PCONFUSION;
                let tmil = 0.5 * (a + b);
                if nupil < max_seg && large {
                    is_cut = true;
                    let from = num_curves as usize + 1;
                    for i in (from..=nupil as usize).rev() {
                        intervals[i + 1] = intervals[i];
                    }
                    intervals[from] = tmil;
                    nupil += 1;
                    continue;
                }
                num_curves += 1;
            }
            err_max[((num_curves - 1) * 2) as usize] = approx.max_error(0);
            err_max[((num_curves - 1) * 2 + 1) as usize] = approx.max_error(1);
            err_avg[((num_curves - 1) * 2) as usize] = approx.average_error(0);
            err_avg[((num_curves - 1) * 2 + 1) as usize] = approx.average_error(1);
            let mut the_deg = approx.degree();
            if is_cut && the_deg < 2 * continuity_order + 1 {
                the_deg = 2 * continuity_order + 1;
            }
            num_coeff[(num_curves - 1) as usize] = the_deg + 1;
            let canon = jac.to_coefficients(2, the_deg, approx.coefficients());
            let f = ((the_deg + 1) * 2) as usize;
            let dest = ((num_curves - 1) * 2 * num_max_coeffs) as usize;
            for i in 0..f.min(canon.len()) {
                if dest + i < coeff.len() {
                    coeff[dest + i] = canon[i];
                }
            }
        }
        let _ = (&err_max, &err_avg);
        if error_code != 0 && error_code != -1 {
            return Ok(Self {
                done: false,
                has_result: false,
                poles_u: Vec::new(),
                poles_v: Vec::new(),
                knots: Vec::new(),
                mults: Vec::new(),
                degree: 0,
            });
        }
        for i in 0..num_curves as usize {
            num_coeff[i] = num_coeff[i].max(2);
        }
        let mut poly_iv = vec![[-1.0, 1.0]; num_curves as usize];
        for i in 0..num_curves as usize {
            poly_iv[i] = [-1.0, 1.0];
        }
        let true_iv = intervals[..=num_curves as usize].to_vec();
        let continuity_v = prepare_convert(
            num_curves,
            max_degree,
            continuity_order,
            2,
            0,
            0,
            &num_coeff[..num_curves as usize],
            &coeff,
            &poly_iv,
            &true_iv,
            &local_tol,
            &mut err_max,
        );
        let conv = CompPolynomialToPoles::new(
            num_curves,
            2,
            max_degree,
            &continuity_v,
            &num_coeff[..num_curves as usize],
            &coeff,
            &poly_iv,
            &true_iv,
        )?;
        if !conv.done {
            return Err(());
        }
        let mut poles_u = Vec::with_capacity(conv.poles.len());
        let mut poles_v = Vec::with_capacity(conv.poles.len());
        for p in &conv.poles {
            poles_u.push(*p.first().unwrap_or(&0.0));
            poles_v.push(*p.get(1).unwrap_or(&0.0));
        }
        let done = error_code == 0;
        Ok(Self {
            done,
            has_result: true,
            poles_u,
            poles_v,
            knots: conv.knots,
            mults: conv.mults,
            degree: conv.degree,
        })
    }
}
