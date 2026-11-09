//! `impl ApproxAFunction1dPair`.

use super::prelude::*;
use super::prepare_convert::*;
use super::ApproxAFunction1dPair;

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
            let mut deg = approx.degree();
            if is_cut && deg < 2 * continuity_order + 1 {
                deg = 2 * continuity_order + 1;
            }
            num_coeff[(num_curves - 1) as usize] = deg + 1;
            let canon = jac.to_coefficients(2, deg, approx.coefficients());
            let f = ((deg + 1) * 2) as usize;
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
