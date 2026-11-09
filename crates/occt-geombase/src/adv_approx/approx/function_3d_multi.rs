//! `ApproxAFunction3dMulti` and its `impl`.

use super::prelude::*;
use super::prepare_convert::*;

/// `GeomFill_ConstrainedFilling.cxx:714-751`, passes `N = nbd3` 1 or 2).
///
/// Same control flow as [`ApproxAFunction3d::approx`] with the single 3D
/// subspace replaced by `tols.len()` of them: the evaluator writes
/// `3 * N` values, `MaxError(3, i)` must satisfy `tols[i]` on every subspace
/// (`cxx:505-519`) and `Poles(3, i)` is the `i`-th pole array
/// (`cxx:823-860`).
pub struct ApproxAFunction3dMulti {
    pub done: bool,
    pub has_result: bool,
    /// `Poles(3, i)` for every subspace, 0-based (`i` is OCCT's 1-based index).
    pub poles: Vec<Vec<GpPnt>>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
    pub max_error: f64,
    pub average_error: f64,
}

impl ApproxAFunction3dMulti {
    /// `AdvApprox_ApproxAFunction` (`cxx:632-659`) / `Perform`
    /// (`cxx:663-956`) / `Approximation` (`cxx:364-599`) for
    /// `Num1DSS = Num2DSS = 0`, `Num3DSS = tols.len()`.
    #[allow(clippy::too_many_arguments)]
    pub fn approx(
        first: f64,
        last: f64,
        continuity: Shape,
        max_deg: i32,
        max_seg: i32,
        tols: &[f64],
        cut: &dyn Cutting,
        eval: EvalFn<'_>,
    ) -> Result<Self, ()> {
        if last < first || max_deg < 1 || max_seg < 0 || tols.is_empty() {
            return Err(());
        }
        let num_3d_ss = tols.len() as i32;
        let total_dim = 3 * tols.len();
        let continuity_order = match continuity {
            Shape::C0 => 0i32,
            Shape::C1 => 1i32,
            Shape::C2 => 2i32,
            // `cxx:699-700` throws `Standard_ConstructionError`.
            _ => return Err(()),
        };
        let max_deg = max_deg.min(14);
        let num_max_coeffs = (max_deg + 1).max(2 * continuity_order + 2);
        let max_degree = num_max_coeffs - 1;
        let (nb_gauss, work_degree) = jacobi_parameters(continuity, max_degree, 1)?;
        let jac = JacobiPolynomial::new(work_degree, continuity)?;
        let mut approx = SimpleApprox::new(
            total_dim,
            num_3d_ss,
            continuity,
            work_degree,
            nb_gauss,
            jac.clone(),
        )?;
        let local_dim = vec![3i32; tols.len()];
        let local_tol = tols.to_vec();
        let mut intervals = vec![0.0; (max_seg as usize) + 1];
        intervals[0] = first;
        intervals[1] = last;
        let mut nupil = 1i32;
        let mut num_curves = 0i32;
        let mut is_cut = false;
        let mut error_code = 0i32;
        let mut num_coeff = vec![0i32; max_seg as usize];
        let mut coeff = vec![0.0; (max_seg * num_max_coeffs * total_dim as i32) as usize];
        let mut err_max = vec![0.0; max_seg as usize];
        let mut err_avg = vec![0.0; max_seg as usize];
        // `cxx:415-419`: `MaxSegments < 1 || |Last - First| < 1e-9` -> ErrorCode 1.
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
            // `cxx:505-519`: every subspace must satisfy its own tolerance.
            let ok_tol = (0..tols.len()).all(|i| approx.max_error(i) <= tols[i]);
            if ok_tol {
                num_curves += 1;
            } else {
                let a = intervals[num_curves as usize];
                let b = intervals[num_curves as usize + 1];
                // `cxx:526`: `Large = CutTool.Value(a, b, TMIL)`.
                let mut tmil = 0.0;
                let large = cut.value(a, b, &mut tmil);
                if nupil < max_seg && large {
                    // `cxx:532-545`: insert TMIL as the new right end of the
                    // current interval and shift the remaining stack.
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
            // `cxx:561-566`.
            err_max[(num_curves - 1) as usize] = (0..tols.len())
                .map(|i| approx.max_error(i))
                .fold(0.0f64, f64::max);
            err_avg[(num_curves - 1) as usize] =
                (0..tols.len()).map(|i| approx.average_error(i)).sum::<f64>() / tols.len() as f64;
            let mut deg = approx.degree();
            if is_cut && deg < 2 * continuity_order + 1 {
                deg = 2 * continuity_order + 1;
            }
            // `cxx:577-584`.
            num_coeff[(num_curves - 1) as usize] = deg + 1;
            let canon = jac.to_coefficients(total_dim as i32, deg, approx.coefficients());
            let f = ((deg + 1) * total_dim as i32) as usize;
            let dest = ((num_curves - 1) * total_dim as i32 * num_max_coeffs) as usize;
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
                poles: Vec::new(),
                knots: Vec::new(),
                mults: Vec::new(),
                degree: 0,
                max_error: 0.0,
                average_error: 0.0,
            });
        }
        if num_curves <= 0 {
            return Err(());
        }
        // `cxx:782-788`: force a minimum degree of 1 (PRO5474).
        for i in 0..num_curves as usize {
            num_coeff[i] = num_coeff[i].max(2);
        }
        let poly_iv = vec![[-1.0, 1.0]; num_curves as usize];
        let true_iv = intervals[..=num_curves as usize].to_vec();
        // `cxx:790-802`.
        let continuity = prepare_convert(
            num_curves,
            max_degree,
            continuity_order,
            0,
            0,
            num_3d_ss,
            &num_coeff[..num_curves as usize],
            &coeff,
            &poly_iv,
            &true_iv,
            &local_tol,
            &mut err_max,
        );
        // `cxx:804-811`.
        let conv = CompPolynomialToPoles::new(
            num_curves,
            total_dim as i32,
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
        // `cxx:923-943`: the reported errors are the maxima over the spans.
        let max_error = (0..num_curves as usize)
            .map(|i| err_max[i])
            .fold(0.0f64, f64::max);
        let average_error =
            (0..num_curves as usize).map(|i| err_avg[i]).sum::<f64>() / num_curves as f64;
        // `Poles(3, i)`: each row holds the `3 * N` coordinates of every
        // subspace, subspace-major.
        let mut poles = vec![Vec::with_capacity(conv.poles.len()); tols.len()];
        for p in &conv.poles {
            for (i, slot) in poles.iter_mut().enumerate() {
                let base = 3 * i;
                slot.push(GpPnt::new(p[base], p[base + 1], p[base + 2]));
            }
        }
        let done = error_code == 0;
        Ok(Self {
            done,
            has_result: true,
            poles,
            knots: conv.knots,
            mults: conv.mults,
            degree: conv.degree,
            max_error,
            average_error,
        })
    }
}
