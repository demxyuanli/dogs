//! `AdvApprox_ApproxAFunction` for a mixed set of 1D, 2D and 3D subspaces.
//! Source: `AdvApprox_ApproxAFunction.cxx` (constructor `cxx:603-659`,
//! `Perform` `cxx:663-956`, `Approximation` `cxx:364-599`, `Poles`
//! `cxx:960-969`, `Poles2d` `cxx:982-991`, `Poles1d` `cxx:993-1002`,
//! `MaxError` `cxx:1004-1037`, `AverageError` `cxx:1039-1070`).
//!
//! The other layouts in this module ([`super::ApproxAFunction3d`],
//! [`super::ApproxAFunction2d`], [`super::ApproxAFunction1dPair`],
//! [`super::ApproxAFunction1d3d`], [`super::ApproxAFunction3dMulti`]) are the
//! fixed subspace counts the curve approximations use. `Approx_SweepApproximation`
//! (`Approx_SweepApproximation.cxx:292-307`) needs the general layout:
//! `Num1DSS = Num3DSS = NbPoles` (weights plus poles, in both subspaces) and
//! `Num2DSS = Nb2dCurves`, so this module threads every subspace index.
//!
//! Translation notes:
//! - `NCollection_Array1` / `_Array2` become `Vec`; `Poles(Index, P)` is
//!   [`ApproxAFunction::poles`], `Poles2d` / `Poles1d` likewise, all 0-based.
//! - `MaxError(D, Index)` / `AverageError(D, Index)` keep OCCT's `(D, Index)`
//!   pair.
//! - The `Standard_ConstructionError` raises of `Perform` (`cxx:667-672`,
//!   `cxx:699-700`) become `Err(())`, like the other ports in this module. The
//!   `Approximation` error codes are not raises: a non-zero code yields the
//!   `done = false` result, exactly like `cxx:776-779`.

use occt_gp::bspl::comp_poly::CompPolynomialToPoles;
use occt_gp::bspl::plib_eval::eval_polynomial;
use occt_gp::bspl::plib_jacobi::{jacobi_parameters, JacobiPolynomial};
use occt_gp::gp::{GpPnt, GpPnt2d};
use occt_foundation::kernel::geomabs::Shape;

use super::cutting::Cutting;
use super::simple::{EvalFn, SimpleApprox};

/// Result of a mixed-dimension `AdvApprox_ApproxAFunction` run.
#[derive(Default)]
pub struct ApproxAFunction {
    pub done: bool,
    pub has_result: bool,
    /// `Degree()`.
    pub degree: i32,
    /// `Knots()`.
    pub knots: Vec<f64>,
    /// `Multiplicities()`.
    pub mults: Vec<i32>,
    /// Number of poles, i.e. `Poles(Index, P).Length()`.
    nb_poles: usize,
    /// `Poles1d(1, Index)` for every 1D subspace, 0-based (`cxx:993-1002`).
    poles_1d: Vec<Vec<f64>>,
    /// `Poles2d(Index, P)` for every 2D subspace, 0-based (`cxx:982-991`).
    poles_2d: Vec<Vec<GpPnt2d>>,
    /// `Poles(Index, P)` for every 3D subspace, 0-based (`cxx:960-969`).
    poles_3d: Vec<Vec<GpPnt>>,
    max_error_1d: Vec<f64>,
    average_error_1d: Vec<f64>,
    max_error_2d: Vec<f64>,
    average_error_2d: Vec<f64>,
    max_error_3d: Vec<f64>,
    average_error_3d: Vec<f64>,
}

/// Outputs of `AdvApprox_ApproxAFunction::Approximation` (`cxx:364-599`).
#[derive(Default)]
struct ApproximationResult {
    error_code: i32,
    num_curves: i32,
    num_coeff_per_curve: Vec<i32>,
    coefficients: Vec<f64>,
    intervals: Vec<f64>,
    error_max: Vec<f64>,
    average_error: Vec<f64>,
}

/// `AdvApprox_ApproxAFunction::PrepareConvert` (`cxx:123-311`).
///
/// Fills the achieved continuity order per interior knot and accumulates the
/// error bounds into `error_max`, indexed `(icurve - 1) * nb_space + idim` and
/// `+ nb_space` for the "next" term — OCCT's `ErrorMax(index)` /
/// `ErrorMax(index + NbSpace)`.
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
            // `cxx:190`: the divisor guard tolerance.
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
            // `cxx:226-247`: 1D subspaces.
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
            // `cxx:251-275`: 2D subspaces.
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
            // `cxx:277-299`: 3D subspaces.
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
            // `cxx:300-310`.
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

/// `AdvApprox_ApproxAFunction::Approximation` (`cxx:364-599`): the cutting loop.
#[allow(clippy::too_many_arguments)]
fn approximation(
    total_dim: usize,
    total_num_ss: i32,
    local_dim: &[i32],
    first: f64,
    last: f64,
    cut: &dyn Cutting,
    continuity_order: i32,
    num_max_coeffs: i32,
    max_seg: i32,
    local_tol: &[f64],
    eval: EvalFn<'_>,
) -> Result<ApproximationResult, ()> {
    let mut out = ApproximationResult::default();
    // `cxx:415-419`.
    if max_seg < 1 || (last - first).abs() < 1.0e-9 {
        out.error_code = 1;
        return Ok(out);
    }
    // `cxx:421-431`: the total dimension must be the sum of the subspaces.
    let idim: i32 = local_dim.iter().sum();
    if idim != total_dim as i32 {
        out.error_code = 1;
        return Ok(out);
    }
    // `cxx:432-445` (`GeomAbs_Shape` of the continuity order).
    let continuity = match continuity_order {
        0 => Shape::C0,
        1 => Shape::C1,
        2 => Shape::C2,
        _ => return Err(()),
    };
    // `cxx:452`: `PLib::JacobiParameters`.
    let (nb_gauss, work_degree) = jacobi_parameters(continuity, num_max_coeffs - 1, 1)?;
    // `cxx:463-477`.
    let jac = JacobiPolynomial::new(work_degree, continuity)?;
    let mut approx = SimpleApprox::new(total_dim, total_num_ss, continuity, work_degree, nb_gauss, jac.clone())?;

    let mut intervals = vec![0.0f64; max_seg as usize + 1];
    // `cxx:459-462`.
    intervals[0] = first;
    intervals[1] = last;
    let mut nupil = 1i32;
    let mut num_curves = 0i32;
    let mut is_cut = false;
    let max_degree = num_max_coeffs - 1;

    out.num_coeff_per_curve = vec![0i32; max_seg as usize];
    out.coefficients = vec![0.0; (max_seg * num_max_coeffs * total_dim as i32) as usize];
    out.error_max = vec![0.0; (max_seg * total_num_ss) as usize];
    out.average_error = vec![0.0; (max_seg * total_num_ss) as usize];

    // `cxx:479`: loop while the stack of intervals has not been consumed.
    while nupil - num_curves != 0 {
        approx.perform(
            local_dim,
            local_tol,
            intervals[num_curves as usize],
            intervals[num_curves as usize + 1],
            max_degree,
            eval,
        );
        if !approx.is_done() {
            out.error_code = 1;
            out.intervals = intervals;
            out.num_curves = num_curves;
            return Ok(out);
        }
        out.num_coeff_per_curve[num_curves as usize] = 0;

        // `cxx:504-519`: every subspace must satisfy its own tolerance.
        let mut max_max_err = true;
        for is in 0..total_num_ss as usize {
            if approx.max_error(is) > local_tol[is] {
                max_max_err = false;
                break;
            }
        }
        if max_max_err {
            num_curves += 1;
        } else {
            // `cxx:521-553`: split the current interval in two.
            let mut tmil = 0.0;
            let large = cut.value(
                intervals[num_curves as usize],
                intervals[num_curves as usize + 1],
                &mut tmil,
            );
            if nupil < max_seg && large {
                is_cut = true;
                // Shift the tail of the stack and insert TMIL.
                let from = num_curves as usize + 1;
                for i in (from..=nupil as usize).rev() {
                    intervals[i + 1] = intervals[i];
                }
                intervals[from] = tmil;
                nupil += 1;
                // `cxx:557-560`: start over with the shortened interval.
                continue;
            }
            num_curves += 1;
        }
        // `cxx:561-566`.
        for is in 0..total_num_ss as usize {
            let idx = is + (num_curves - 1) as usize * total_num_ss as usize;
            out.error_max[idx] = approx.max_error(is);
            out.average_error[idx] = approx.average_error(is);
        }
        let mut deg = approx.degree();
        // `cxx:570-575`.
        if is_cut && deg < 2 * continuity_order + 1 {
            deg = 2 * continuity_order + 1;
        }
        out.num_coeff_per_curve[num_curves as usize] = deg + 1;
        // `cxx:579-586`.
        let canon = jac.to_coefficients(total_dim as i32, deg, approx.coefficients());
        let f = ((deg + 1) * total_dim as i32) as usize;
        let dest = ((num_curves - 1) * total_dim as i32 * num_max_coeffs) as usize;
        for i in 0..f.min(canon.len()) {
            if dest + i < out.coefficients.len() {
                out.coefficients[dest + i] = canon[i];
            }
        }
    }
    out.num_curves = num_curves;
    out.intervals = intervals;
    Ok(out)
}

impl ApproxAFunction {
    /// `AdvApprox_ApproxAFunction(Num1DSS, Num2DSS, Num3DSS, OneDTol, TwoDTol,
    /// ThreeDTol, First, Last, Continuity, MaxDeg, MaxSeg, Func, CutTool)`
    /// (`cxx:631-659`) followed by `Perform` (`cxx:663-956`).
    #[allow(clippy::too_many_arguments)]
    pub fn approx(
        num_1d_ss: i32,
        num_2d_ss: i32,
        num_3d_ss: i32,
        one_d_tol: &[f64],
        two_d_tol: &[f64],
        three_d_tol: &[f64],
        first: f64,
        last: f64,
        continuity: Shape,
        max_deg: i32,
        max_seg: i32,
        cut: &dyn Cutting,
        eval: EvalFn<'_>,
    ) -> Result<Self, ()> {
        // `cxx:667-672`.
        if num_1d_ss < 0
            || num_2d_ss < 0
            || num_3d_ss < 0
            || num_1d_ss + num_2d_ss + num_3d_ss <= 0
            || last < first
            || max_deg < 1
            || max_seg < 0
        {
            return Err(());
        }
        if one_d_tol.len() < num_1d_ss as usize
            || two_d_tol.len() < num_2d_ss as usize
            || three_d_tol.len() < num_3d_ss as usize
        {
            return Err(());
        }
        if (num_1d_ss + 2 * num_2d_ss + 3 * num_3d_ss) <= 0 {
            return Err(());
        }
        // `cxx:673-676`.
        let max_deg = max_deg.min(14);
        let total_num_ss = num_1d_ss + num_2d_ss + num_3d_ss;
        let total_dim = (num_1d_ss + 2 * num_2d_ss + 3 * num_3d_ss) as usize;
        // `cxx:687-701`.
        let continuity_order = match continuity {
            Shape::C0 => 0i32,
            Shape::C1 => 1i32,
            Shape::C2 => 2i32,
            _ => return Err(()),
        };
        let num_max_coeffs = (max_deg + 1).max(2 * continuity_order + 2);
        let max_degree = num_max_coeffs - 1;

        // `cxx:720-735`: 1D subspaces, then 2D, then 3D.
        let mut local_dim = Vec::with_capacity(total_num_ss as usize);
        let mut local_tol = Vec::with_capacity(total_num_ss as usize);
        for jj in 0..num_1d_ss as usize {
            local_tol.push(one_d_tol[jj]);
            local_dim.push(1i32);
        }
        for jj in 0..num_2d_ss as usize {
            local_tol.push(two_d_tol[jj]);
            local_dim.push(2i32);
        }
        for jj in 0..num_3d_ss as usize {
            local_tol.push(three_d_tol[jj]);
            local_dim.push(3i32);
        }

        // `cxx:404-599`.
        let run = approximation(
            total_dim,
            total_num_ss,
            &local_dim,
            first,
            last,
            cut,
            continuity_order,
            num_max_coeffs,
            max_seg,
            &local_tol,
            eval,
        )?;

        let mut me = Self::default();
        // `cxx:776-779`: only codes 0 and -1 produce a result.
        if run.error_code != 0 && run.error_code != -1 {
            return Ok(me);
        }
        let num_curves = run.num_curves;
        if num_curves <= 0 {
            return Err(());
        }
        let mut num_coeff = run.num_coeff_per_curve;
        // `cxx:782-788`: force a minimum degree of 1 (PRO5474).
        for nc in num_coeff.iter_mut().take(num_curves as usize) {
            *nc = (*nc).max(2);
        }
        let poly_iv = vec![[-1.0, 1.0]; num_curves as usize];
        let true_iv = run.intervals[..=num_curves as usize].to_vec();
        let mut error_max = run.error_max;
        // `cxx:790-802`.
        let tab_continuity = prepare_convert(
            num_curves,
            max_degree,
            continuity_order,
            num_1d_ss,
            num_2d_ss,
            num_3d_ss,
            &num_coeff,
            &run.coefficients,
            &poly_iv,
            &true_iv,
            &local_tol,
            &mut error_max,
        );
        // `cxx:804-811`.
        let conv = CompPolynomialToPoles::new(
            num_curves,
            total_dim as i32,
            max_degree,
            &tab_continuity,
            &num_coeff,
            &run.coefficients,
            &poly_iv,
            &true_iv,
        )?;
        if !conv.done {
            return Err(());
        }
        me.knots = conv.knots.clone();
        me.mults = conv.mults.clone();
        me.degree = conv.degree;
        me.nb_poles = conv.poles.len();

        let nb_poles = conv.poles.len();
        let total_ss = total_num_ss as usize;
        // `cxx:820-853`: `Poles1d`.
        let mut index = 0usize;
        if num_1d_ss > 0 {
            me.poles_1d = vec![Vec::with_capacity(nb_poles); num_1d_ss as usize];
            for p in &conv.poles {
                for jj in 0..num_1d_ss as usize {
                    me.poles_1d[jj].push(p[index + jj]);
                }
            }
            for jj in 0..num_1d_ss as usize {
                let mut value = 0.0f64;
                for ii in 0..num_curves as usize {
                    value = value.max(error_max[ii * total_ss + jj]);
                }
                me.max_error_1d.push(value);
            }
            for jj in 0..num_1d_ss as usize {
                let mut value = 0.0f64;
                for ii in 0..num_curves as usize {
                    value += run.average_error[ii * total_ss + jj];
                }
                me.average_error_1d.push(value / num_curves as f64);
            }
            index += num_1d_ss as usize;
        }

        let mut dim_index = index;
        // `cxx:859-901`: `Poles2d`.
        if num_2d_ss > 0 {
            me.poles_2d = vec![Vec::with_capacity(nb_poles); num_2d_ss as usize];
            for p in &conv.poles {
                for jj in 0..num_2d_ss as usize {
                    let local = index + jj * 2;
                    me.poles_2d[jj].push(GpPnt2d::new(p[local], p[local + 1]));
                }
            }
            for jj in 0..num_2d_ss as usize {
                let mut value = 0.0f64;
                for ii in 0..num_curves as usize {
                    value = value.max(error_max[ii * total_ss + index + jj]);
                }
                me.max_error_2d.push(value);
            }
            for jj in 0..num_2d_ss as usize {
                let mut value = 0.0f64;
                for ii in 0..num_curves as usize {
                    value += run.average_error[ii * total_ss + index + jj];
                }
                me.average_error_2d.push(value / num_curves as f64);
            }
            index += num_2d_ss as usize;
            // `cxx:900-901`: the pole offset must be doubled for 3D subspaces.
            dim_index = index + num_2d_ss as usize;
        }

        // `cxx:904-943`: `Poles`.
        if num_3d_ss > 0 {
            me.poles_3d = vec![Vec::with_capacity(nb_poles); num_3d_ss as usize];
            for p in &conv.poles {
                for jj in 0..num_3d_ss as usize {
                    let local = dim_index + jj * 3;
                    me.poles_3d[jj].push(GpPnt::new(p[local], p[local + 1], p[local + 2]));
                }
            }
            for jj in 0..num_3d_ss as usize {
                let mut value = 0.0f64;
                for ii in 0..num_curves as usize {
                    value = value.max(error_max[ii * total_ss + index + jj]);
                }
                me.max_error_3d.push(value);
            }
            for jj in 0..num_3d_ss as usize {
                let mut value = 0.0f64;
                for ii in 0..num_curves as usize {
                    value += run.average_error[ii * total_ss + index + jj];
                }
                me.average_error_3d.push(value / num_curves as f64);
            }
        }

        // `cxx:944-953`.
        if run.error_code == 0 {
            me.done = true;
            me.has_result = true;
        } else {
            me.has_result = true;
        }
        Ok(me)
    }

    /// `Poles(Index, P)` (`cxx:960-969`).
    pub fn poles(&self, index: usize) -> &[GpPnt] {
        &self.poles_3d[index]
    }

    /// `Poles2d(Index, P)` (`cxx:982-991`).
    pub fn poles_2d(&self, index: usize) -> &[GpPnt2d] {
        &self.poles_2d[index]
    }

    /// `Poles1d(Index, P)` (`cxx:993-1002`).
    pub fn poles_1d(&self, index: usize) -> &[f64] {
        &self.poles_1d[index]
    }

    /// `NbPoles()` (`cxx:971-980`).
    pub fn nb_poles(&self) -> usize {
        self.nb_poles
    }

    /// `MaxError(D, Index)` (`cxx:1030-1037`).
    pub fn max_error(&self, d: i32, index: usize) -> f64 {
        match d {
            1 => self.max_error_1d[index],
            2 => self.max_error_2d[index],
            _ => self.max_error_3d[index],
        }
    }

    /// `AverageError(D, Index)` (`cxx:1065-1070`).
    pub fn average_error(&self, d: i32, index: usize) -> f64 {
        match d {
            1 => self.average_error_1d[index],
            2 => self.average_error_2d[index],
            _ => self.average_error_3d[index],
        }
    }
}
