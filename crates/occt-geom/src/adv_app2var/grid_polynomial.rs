//! `Convert_GridPolynomialToPoles` (`Convert_GridPolynomialToPoles.cxx:29-257`).

use occt_core::bspl::banded_interp::{knot_sequence, schoenberg_points};
use occt_core::bspl::plib_eval::eval_poly2var;
use occt_core::bspl::surface_interpolate::interpolate_surface;
use occt_core::gp::GpPnt;

pub(super) struct GridSurface {
    pub poles: Vec<Vec<GpPnt>>,
    pub u_knots: Vec<f64>,
    pub v_knots: Vec<f64>,
    pub u_mults: Vec<i32>,
    pub v_mults: Vec<i32>,
    pub u_degree: i32,
    pub v_degree: i32,
}

/// One polynomial (`Convert_GridPolynomialToPoles.cxx:29-78`).
#[allow(dead_code)]
pub(super) fn one_polynomial_to_poles(
    max_u_degree: i32,
    max_v_degree: i32,
    num_coeff_u: i32,
    num_coeff_v: i32,
    coefficients: &[f64],
    polynomial_u: &[f64],
    polynomial_v: &[f64],
) -> Result<GridSurface, String> {
    let expect = checked_len(&[3, max_u_degree + 1, max_v_degree + 1], 46)?;
    if coefficients.len() != expect {
        return Err("Convert : Wrong Coefficients".into());
    }
    if num_coeff_u - 1 > max_u_degree {
        return Err("Convert : Incoherence between NumCoeffPerSurface and MaxUDegree".into());
    }
    if num_coeff_v - 1 > max_v_degree {
        return Err("Convert : Incoherence between NumCoeffPerSurface and MaxVDegree".into());
    }
    perform(
        0,
        0,
        max_u_degree,
        max_v_degree,
        &[(num_coeff_u, num_coeff_v)],
        coefficients,
        polynomial_u,
        polynomial_v,
        polynomial_u,
        polynomial_v,
        num_coeff_u - 1,
        num_coeff_v - 1,
    )
}

/// A grid of polynomials (`Convert_GridPolynomialToPoles.cxx:80-147`).
pub(super) fn grid_polynomial_to_poles(
    nb_u: i32,
    nb_v: i32,
    u_continuity: i32,
    v_continuity: i32,
    max_u_degree: i32,
    max_v_degree: i32,
    num_coeff: &[(i32, i32)],
    coefficients: &[f64],
    polynomial_u: &[f64],
    polynomial_v: &[f64],
    true_u: &[f64],
    true_v: &[f64],
) -> Result<GridSurface, String> {
    let real_u = max_u_degree.max(2 * u_continuity + 1);
    let real_v = max_v_degree.max(2 * v_continuity + 1);
    let n_patch = checked_len(&[nb_u, nb_v], 102)?;
    if num_coeff.len() != n_patch {
        return Err("Convert : Wrong NumCoeffPerSurface".into());
    }
    let expect = checked_len(&[3, nb_u, nb_v, real_u + 1, real_v + 1], 110)?;
    if coefficients.len() != expect {
        return Err("Convert : Wrong Coefficients".into());
    }
    let mut u_degree = 0i32;
    let mut v_degree = 0i32;
    for &(nu, nv) in num_coeff {
        if nu > u_degree + 1 {
            u_degree = nu - 1;
        }
        if nv > v_degree + 1 {
            v_degree = nv - 1;
        }
    }
    if u_degree > real_u {
        return Err("Convert : Incoherence between NumCoeffPerSurface and MaxUDegree".into());
    }
    if v_degree > real_v {
        return Err("Convert : Incoherence between NumCoeffPerSurface and MaxVDegree".into());
    }
    perform(
        u_continuity,
        v_continuity,
        real_u,
        real_v,
        num_coeff,
        coefficients,
        polynomial_u,
        polynomial_v,
        true_u,
        true_v,
        u_degree,
        v_degree,
    )
}

fn perform(
    u_continuity: i32,
    v_continuity: i32,
    max_u_degree: i32,
    max_v_degree: i32,
    num_coeff: &[(i32, i32)],
    coefficients: &[f64],
    polynomial_u: &[f64],
    polynomial_v: &[f64],
    true_u: &[f64],
    true_v: &[f64],
    u_degree: i32,
    v_degree: i32,
) -> Result<GridSurface, String> {
    if polynomial_u.len() < 2 || polynomial_v.len() < 2 || true_u.len() < 2 || true_v.len() < 2 {
        return Err(convert_line(187));
    }
    let (u_flat, u_mults, u_params) = build_array(u_degree, true_u, u_continuity)?;
    let (v_flat, v_mults, v_params) = build_array(v_degree, true_v, v_continuity)?;
    let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); v_params.len()]; u_params.len()];
    let mut patch_index = 0i32;
    let mut patch_coeffs: Vec<f64> = Vec::new();
    if u_params.is_empty() {
        return Err(convert_line(184));
    }
    for ii in 1..=u_params.len() as i32 {
        let mut u_span = 1i32;
        while u_params[(ii - 1) as usize] > true_u[u_span as usize] && u_span < true_u.len() as i32 - 1
        {
            u_span += 1;
        }
        let u_value = span_to_poly(u_params[(ii - 1) as usize], true_u, u_span, polynomial_u)?;
        for jj in 1..=v_params.len() as i32 {
            let mut v_span = 1i32;
            while v_params[(jj - 1) as usize] > true_v[v_span as usize]
                && v_span < true_v.len() as i32 - 1
            {
                v_span += 1;
            }
            let v_value = span_to_poly(v_params[(jj - 1) as usize], true_v, v_span, polynomial_v)?;
            let index = u_span + (true_u.len() as i32 - 1) * (v_span - 1);
            if patch_index != index {
                patch_coeffs = extract_patch(
                    index,
                    num_coeff,
                    coefficients,
                    max_u_degree,
                    max_v_degree,
                )?;
                patch_index = index;
            }
            let (nu, nv) = *num_coeff.get((index - 1) as usize).ok_or_else(|| convert_line(217))?;
            if nu < 1 || nv < 1 {
                return Err(convert_line(235));
            }
            let mut digit = [0.0; 3];
            eval_poly2var(u_value, v_value, nu - 1, nv - 1, 3, &patch_coeffs, &mut digit);
            poles[(ii - 1) as usize][(jj - 1) as usize] = GpPnt::new(digit[0], digit[1], digit[2]);
        }
    }
    if !interpolate_surface(
        u_degree,
        v_degree,
        &u_flat,
        &v_flat,
        &u_params,
        &v_params,
        &mut poles,
    ) {
        return Err("GridPolynomialToPoles".into());
    }
    Ok(GridSurface {
        poles,
        u_knots: true_u.to_vec(),
        v_knots: true_v.to_vec(),
        u_mults,
        v_mults,
        u_degree,
        v_degree,
    })
}

fn extract_patch(
    index: i32,
    num_coeff: &[(i32, i32)],
    coefficients: &[f64],
    max_u_degree: i32,
    max_v_degree: i32,
) -> Result<Vec<f64>, String> {
    let (nu, nv) = *num_coeff.get((index - 1) as usize).ok_or_else(|| convert_line(217))?;
    let siz = checked_len(&[3, max_u_degree + 1, max_v_degree + 1], 176)?;
    let v_stride = checked_len(&[3, max_v_degree + 1], 219)?;
    let mut packed = Vec::new();
    if nu >= 1 {
        for k1 in 1..=nu {
            let row = (index - 1) as usize * siz + (k1 - 1) as usize * v_stride;
            if nv >= 1 {
                for k2 in 1..=nv {
                    let pos = row + (k2 - 1) as usize * 3;
                    let x = *coefficients.get(pos).ok_or_else(|| convert_line(222))?;
                    let y = *coefficients.get(pos + 1).ok_or_else(|| convert_line(223))?;
                    let z = *coefficients.get(pos + 2).ok_or_else(|| convert_line(224))?;
                    packed.push(x);
                    packed.push(y);
                    packed.push(z);
                }
            }
        }
    }
    Ok(packed)
}

fn span_to_poly(parameter: f64, knots: &[f64], span: i32, poly: &[f64]) -> Result<f64, String> {
    let left = *knots.get((span - 1) as usize).ok_or_else(|| convert_line(193))?;
    let right = *knots.get(span as usize).ok_or_else(|| convert_line(194))?;
    let width = right - left;
    if width == 0.0 {
        return Err(convert_line(194));
    }
    let n_value = (parameter - left) / width;
    let a = *poly.first().ok_or_else(|| convert_line(196))?;
    let b = *poly.get(1).ok_or_else(|| convert_line(196))?;
    Ok((1.0 - n_value) * a + n_value * b)
}

/// `BuildArray` (`Convert_GridPolynomialToPoles.cxx:259-291`).
fn build_array(
    degree: i32,
    knots: &[f64],
    continuity: i32,
) -> Result<(Vec<f64>, Vec<i32>, Vec<f64>), String> {
    if degree < 1 {
        return Err("BSplCLib::BuildSchoenbergPoints: BSplCLib.cxx:3331".into());
    }
    let num_curves = knots.len() as i32 - 1;
    if num_curves < 1 {
        return Err(convert_line(266));
    }
    let multiplicities = degree - continuity;
    if multiplicities < 0 {
        return Err(convert_line(269));
    }
    let mut mults = vec![multiplicities; knots.len()];
    if let Some(first) = mults.first_mut() {
        *first = degree + 1;
    }
    if let Some(last) = mults.last_mut() {
        *last = degree + 1;
    }
    let num_flat = (multiplicities as usize)
        .checked_mul((num_curves - 1) as usize)
        .and_then(|n| n.checked_add((2 * degree + 2) as usize))
        .ok_or_else(|| convert_line(280))?;
    let flat = knot_sequence(knots, &mults, degree);
    if flat.len() != num_flat {
        return Err(convert_line(283));
    }
    let num_poles = num_flat - degree as usize - 1;
    let params = schoenberg_points(degree as usize, &flat, num_poles);
    Ok((flat, mults, params))
}

fn checked_len(parts: &[i32], line: i32) -> Result<usize, String> {
    let mut n = 1i32;
    for part in parts {
        n = n.checked_mul(*part).ok_or_else(|| convert_line(line))?;
    }
    usize::try_from(n).map_err(|_| convert_line(line))
}

fn convert_line(line: i32) -> String {
    format!("Convert_GridPolynomialToPoles: Convert_GridPolynomialToPoles.cxx:{line}")
}
