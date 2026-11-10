//! `ApproxAFunction1d3d` and its `impl`.

use super::prelude::*;
use super::prepare_convert::*;

/// `AdvApprox_ApproxAFunction` with `Num1DSS=2`, `Num2DSS=0`, `Num3DSS=1`
/// (`Approx_CurvilinearParameter.cxx:340-353`): evaluator dimension 5.
pub struct ApproxAFunction1d3d {
    pub done: bool,
    pub has_result: bool,
    pub poles_1d_u: Vec<f64>,
    pub poles_1d_v: Vec<f64>,
    pub poles_3d: Vec<GpPnt>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
    pub max_error_1d: f64,
    pub max_error_3d: f64,
}

impl ApproxAFunction1d3d {
    /// PrefAndRec cutting, continuity as requested.
    #[allow(clippy::too_many_arguments)]
    pub fn approx(
        first: f64,
        last: f64,
        continuity: Shape,
        max_deg: i32,
        max_seg: i32,
        tol_1d: f64,
        tol_3d: f64,
        cut: &dyn Cutting,
        eval: EvalFn<'_>,
    ) -> Result<Self, ()> {
        if last < first || max_deg < 1 || max_seg < 0 {
            return Err(());
        }
        let continuity_order = match continuity {
            Shape::C0 => 0i32,
            Shape::C1 => 1i32,
            Shape::C2 => 2i32,
            _ => return Err(()),
        };
        let max_deg = max_deg.min(14);
        let num_max_coeffs = (max_deg + 1).max(2 * continuity_order + 2);
        let max_degree = num_max_coeffs - 1;
        let (nb_gauss, work_degree) = jacobi_parameters(continuity, max_degree, 1)?;
        let jac = JacobiPolynomial::new(work_degree, continuity)?;
        let mut approx = SimpleApprox::new(5, 3, continuity, work_degree, nb_gauss, jac.clone())?;
        let local_dim = [1i32, 1i32, 3i32];
        let local_tol = [tol_1d, tol_1d, tol_3d];
        let mut intervals = vec![0.0; (max_seg as usize) + 1];
        intervals[0] = first;
        intervals[1] = last;
        let mut nupil = 1i32;
        let mut num_curves = 0i32;
        let mut is_cut = false;
        let mut error_code = 0i32;
        let mut num_coeff = vec![0i32; max_seg as usize];
        let mut coeff = vec![0.0; (max_seg * num_max_coeffs * 5) as usize];
        let mut err_max = vec![0.0; (max_seg * 3) as usize];
        let mut err_avg = vec![0.0; (max_seg * 3) as usize];
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
            let ok_tol = approx.max_error(0) <= tol_1d
                && approx.max_error(1) <= tol_1d
                && approx.max_error(2) <= tol_3d;
            if ok_tol {
                num_curves += 1;
            } else {
                let a = intervals[num_curves as usize];
                let b = intervals[num_curves as usize + 1];
                let mut tmil = 0.0;
                let large = cut.value(a, b, &mut tmil);
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
            let base = ((num_curves - 1) * 3) as usize;
            err_max[base] = approx.max_error(0);
            err_max[base + 1] = approx.max_error(1);
            err_max[base + 2] = approx.max_error(2);
            err_avg[base] = approx.average_error(0);
            err_avg[base + 1] = approx.average_error(1);
            err_avg[base + 2] = approx.average_error(2);
            let mut deg = approx.degree();
            if is_cut && deg < 2 * continuity_order + 1 {
                deg = 2 * continuity_order + 1;
            }
            num_coeff[(num_curves - 1) as usize] = deg + 1;
            let canon = jac.to_coefficients(5, deg, approx.coefficients());
            let f = ((deg + 1) * 5) as usize;
            let dest = ((num_curves - 1) * 5 * num_max_coeffs) as usize;
            for i in 0..f.min(canon.len()) {
                if dest + i < coeff.len() {
                    coeff[dest + i] = canon[i];
                }
            }
        }
        if error_code != 0 && error_code != -1 {
            return Ok(Self {
                done: false,
                has_result: false,
                poles_1d_u: Vec::new(),
                poles_1d_v: Vec::new(),
                poles_3d: Vec::new(),
                knots: Vec::new(),
                mults: Vec::new(),
                degree: 0,
                max_error_1d: 0.0,
                max_error_3d: 0.0,
            });
        }
        for i in 0..num_curves as usize {
            num_coeff[i] = num_coeff[i].max(2);
        }
        let poly_iv = vec![[-1.0, 1.0]; num_curves as usize];
        let true_iv = intervals[..=num_curves as usize].to_vec();
        let continuity_v = prepare_convert(
            num_curves,
            max_degree,
            continuity_order,
            2,
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
            5,
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
        let mut poles_1d_u = Vec::with_capacity(conv.poles.len());
        let mut poles_1d_v = Vec::with_capacity(conv.poles.len());
        let mut poles_3d = Vec::with_capacity(conv.poles.len());
        for p in &conv.poles {
            poles_1d_u.push(*p.first().unwrap_or(&0.0));
            poles_1d_v.push(*p.get(1).unwrap_or(&0.0));
            poles_3d.push(GpPnt::new(
                *p.get(2).unwrap_or(&0.0),
                *p.get(3).unwrap_or(&0.0),
                *p.get(4).unwrap_or(&0.0),
            ));
        }
        let mut max_error_1d = 0.0_f64;
        let mut max_error_3d = 0.0_f64;
        for i in 0..num_curves as usize {
            max_error_1d = max_error_1d.max(err_max[i * 3]).max(err_max[i * 3 + 1]);
            max_error_3d = max_error_3d.max(err_max[i * 3 + 2]);
        }
        let _ = (&err_avg, &poly_iv);
        Ok(Self {
            done: error_code == 0,
            has_result: true,
            poles_1d_u,
            poles_1d_v,
            poles_3d,
            knots: conv.knots,
            mults: conv.mults,
            degree: conv.degree,
            max_error_1d,
            max_error_3d,
        })
    }
}
