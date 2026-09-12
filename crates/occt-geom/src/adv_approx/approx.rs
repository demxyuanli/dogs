//! `AdvApprox_ApproxAFunction` for a single 3D subspace.
//! Source: `AdvApprox_ApproxAFunction.cxx` (constructor `cxx:603-630`,
//! `Perform` `cxx:663-956`, `Approximation` `cxx:364-599`).

use occt_core::bspl::comp_poly::CompPolynomialToPoles;
use occt_core::bspl::plib_jacobi::{jacobi_parameters, JacobiPolynomial};
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;

use super::simple::{EvalFn, SimpleApprox};

/// Result of a 3D `ApproxAFunction` run.
pub struct ApproxAFunction3d {
    pub done: bool,
    pub has_result: bool,
    pub poles: Vec<GpPnt>,
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
        let continuity = vec![continuity_order; num_curves as usize];
        let mut poly_iv = vec![[-1.0, 1.0]; num_curves as usize];
        for i in 0..num_curves as usize {
            poly_iv[i] = [-1.0, 1.0];
        }
        let true_iv = intervals[..=num_curves as usize].to_vec();
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
