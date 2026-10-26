//! `AdvApprox_ApproxAFunction` with `Num1DSS=2`, `Num2DSS=0`, `Num3DSS=0` or `1`
//! (`Approx_CurveOnSurface.cxx:403-518`): evaluator dimension 5 (or 2 when the
//! 3D subspace is absent), one tolerance per subspace.
//!
//! Source: `ModelingData/TKG3d/AdvApprox/AdvApprox_ApproxAFunction.cxx`
//! (`Approximation` `cxx:364-599`, `Perform` `cxx:663-956`).
//!
//! The control flow is the same as [`super::approx::ApproxAFunction1d3d`], which
//! serves `Approx_CurvilinearParameter` with one tolerance shared by the two
//! 1D subspaces. `Perform` keeps `LocalTolerances` per subspace
//! (`cxx:719-735`), and `Approx_CurveOnSurface` passes `TolU` and `TolV`
//! separately (`Approx_CurveOnSurface.cxx:428-460`), so the pair is kept apart
//! here: the per-curve error table is indexed by subspace
//! (`ErrorMaxArray(IS + (NumCurves-1)*TotalNumSS)`, `cxx:562-567`) and each
//! entry is compared against its own tolerance (`cxx:505-519`).
//!
//! `Approx_CurveOnSurface::Perform` reaches this layout three ways
//! (`cxx:403-518`): `Only3d` (`Num3DSS=1`, `Num1DSS=0`, served by
//! [`super::approx::ApproxAFunction3d`]), `Only2d` (`Num1DSS=2`, `Num3DSS=0`)
//! and the default (both). The last two share `Perform` and differ only in the
//! subspace tables, which is what `tol_3d: Option<f64>` selects.

use crate::bspl::comp_poly::CompPolynomialToPoles;
use crate::bspl::plib_jacobi::{jacobi_parameters, JacobiPolynomial};
use crate::gp::GpPnt;
use crate::kernel::geomabs::Shape;

use super::approx::prepare_convert;
use super::cutting::Cutting;
use super::simple::{EvalFn, SimpleApprox};

/// Result of `AdvApprox_ApproxAFunction` with `Num1DSS=2` and `Num3DSS` 0 or 1
/// (`Approx_CurveOnSurface`).
pub struct ApproxAFunction1dPair3d {
    pub done: bool,
    pub has_result: bool,
    /// `Poles1d(1, P)` (`cxx:995-1000`): the first 1D subspace (U).
    pub poles_1d_u: Vec<f64>,
    /// `Poles1d(2, P)`: the second 1D subspace (V).
    pub poles_1d_v: Vec<f64>,
    /// `Poles(1, P)` (`cxx:960-966`): the single 3D subspace; `None` when
    /// `tol_3d` was `None` (`Only2d`).
    pub poles_3d: Option<Vec<GpPnt>>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
    /// `MaxError(1, 1)` (`cxx:826-841`).
    pub max_error_u: f64,
    /// `MaxError(1, 2)`.
    pub max_error_v: f64,
    /// `MaxError(3, 1)` (`cxx:921-931`); `None` without the 3D subspace.
    pub max_error_3d: Option<f64>,
}

impl ApproxAFunction1dPair3d {
    /// `AdvApprox_ApproxAFunction` (`cxx:632-659`) / `Perform`
    /// (`cxx:663-956`) / `Approximation` (`cxx:364-599`) for
    /// `Num1DSS = 2`, `Num2DSS = 0`, `Num3DSS = tol_3d.is_some() as`.
    #[allow(clippy::too_many_arguments)]
    pub fn approx(
        first: f64,
        last: f64,
        continuity: Shape,
        max_deg: i32,
        max_seg: i32,
        tol_u: f64,
        tol_v: f64,
        tol_3d: Option<f64>,
        cut: &dyn Cutting,
        eval: EvalFn<'_>,
    ) -> Result<Self, ()> {
        let has_3d = tol_3d.is_some();
        let total_ss: i32 = if has_3d { 3 } else { 2 };
        let dimension: i32 = if has_3d { 5 } else { 2 };
        // `cxx:669-672`.
        if last < first || max_deg < 1 || max_seg < 0 {
            return Err(());
        }
        let continuity_order = match continuity {
            Shape::C0 => 0i32,
            Shape::C1 => 1i32,
            Shape::C2 => 2i32,
            // `cxx:699-700` throws `Standard_ConstructionError`.
            _ => return Err(()),
        };
        // `cxx:673-676`, `cxx:704-705`.
        let max_deg = max_deg.min(14);
        let num_max_coeffs = (max_deg + 1).max(2 * continuity_order + 2);
        let max_degree = num_max_coeffs - 1;
        let (nb_gauss, work_degree) = jacobi_parameters(continuity, max_degree, 1)?;
        let jac = JacobiPolynomial::new(work_degree, continuity)?;
        let mut approx = SimpleApprox::new(
            dimension as usize,
            total_ss,
            continuity,
            work_degree,
            nb_gauss,
            jac.clone(),
        )?;
        // `cxx:716-735`: LocalDimension = {1, 1[, 3]}, LocalTolerances per subspace.
        let local_dim: Vec<i32> = if has_3d { vec![1, 1, 3] } else { vec![1, 1] };
        let local_tol: Vec<f64> = match tol_3d {
            Some(t3) => vec![tol_u, tol_v, t3],
            None => vec![tol_u, tol_v],
        };
        let mut intervals = vec![0.0; (max_seg as usize) + 1];
        intervals[0] = first;
        intervals[1] = last;
        let mut nupil = 1i32;
        let mut num_curves = 0i32;
        let mut is_cut = false;
        let mut error_code = 0i32;
        let mut num_coeff = vec![0i32; max_seg as usize];
        let mut coeff = vec![0.0; (max_seg * num_max_coeffs * dimension) as usize];
        // `cxx:750-752`: one error slot per subspace per curve.
        let mut err_max = vec![0.0; (max_seg * total_ss) as usize];
        let mut err_avg = vec![0.0; (max_seg * total_ss) as usize];
        // `cxx:415-419`.
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
            let mut ok_tol = true;
            for (is, tol) in local_tol.iter().enumerate() {
                if approx.max_error(is) > *tol {
                    ok_tol = false;
                    break;
                }
            }
            if ok_tol {
                num_curves += 1;
            } else {
                let a = intervals[num_curves as usize];
                let b = intervals[num_curves as usize + 1];
                // `cxx:526`: `Large = CutTool.Value(a, b, TMIL)`.
                let mut tmil = 0.0;
                let large = cut.value(a, b, &mut tmil);
                if nupil < max_seg && large {
                    // `cxx:532-545`.
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
            // `cxx:561-567`.
            let base = ((num_curves - 1) * total_ss) as usize;
            for is in 0..total_ss as usize {
                err_max[base + is] = approx.max_error(is);
                err_avg[base + is] = approx.average_error(is);
            }
            let mut deg = approx.degree();
            // `cxx:570-575`.
            if is_cut && deg < 2 * continuity_order + 1 {
                deg = 2 * continuity_order + 1;
            }
            num_coeff[(num_curves - 1) as usize] = deg + 1;
            let canon = jac.to_coefficients(dimension, deg, approx.coefficients());
            let f = ((deg + 1) * dimension) as usize;
            let dest = ((num_curves - 1) * dimension * num_max_coeffs) as usize;
            for i in 0..f.min(canon.len()) {
                if dest + i < coeff.len() {
                    coeff[dest + i] = canon[i];
                }
            }
        }
        // `cxx:777-787`: ErrorCode 1 means no result; -1 keeps the best result.
        if error_code != 0 && error_code != -1 {
            return Ok(Self {
                done: false,
                has_result: false,
                poles_1d_u: Vec::new(),
                poles_1d_v: Vec::new(),
                poles_3d: None,
                knots: Vec::new(),
                mults: Vec::new(),
                degree: 0,
                max_error_u: 0.0,
                max_error_v: 0.0,
                max_error_3d: None,
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
        let continuity_v = prepare_convert(
            num_curves,
            max_degree,
            continuity_order,
            2,
            0,
            if has_3d { 1 } else { 0 },
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
            dimension,
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
        // `cxx:820-856` / `cxx:906-931`: per-subspace pole arrays and errors.
        let mut poles_1d_u = Vec::with_capacity(conv.poles.len());
        let mut poles_1d_v = Vec::with_capacity(conv.poles.len());
        let mut poles_3d = if has_3d {
            Some(Vec::with_capacity(conv.poles.len()))
        } else {
            None
        };
        for p in &conv.poles {
            poles_1d_u.push(*p.first().unwrap_or(&0.0));
            poles_1d_v.push(*p.get(1).unwrap_or(&0.0));
            if let Some(out) = poles_3d.as_mut() {
                out.push(GpPnt::new(
                    *p.get(2).unwrap_or(&0.0),
                    *p.get(3).unwrap_or(&0.0),
                    *p.get(4).unwrap_or(&0.0),
                ));
            }
        }
        let mut max_error_u = 0.0_f64;
        let mut max_error_v = 0.0_f64;
        let mut max_error_3d = has_3d.then_some(0.0_f64);
        for i in 0..num_curves as usize {
            let base = i * total_ss as usize;
            max_error_u = max_error_u.max(err_max[base]);
            max_error_v = max_error_v.max(err_max[base + 1]);
            if let Some(e) = max_error_3d.as_mut() {
                *e = e.max(err_max[base + 2]);
            }
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
            max_error_u,
            max_error_v,
            max_error_3d,
        })
    }
}
