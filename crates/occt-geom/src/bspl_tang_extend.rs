//! `BSplCLib::TangExtendToConstraint` — tangent extension of a non-rational
//! N-dimensional B-spline curve up to a constraint point.
//!
//! Source: `FoundationClasses/TKMath/BSplCLib/BSplCLib.cxx:3794-4238`
//! (OCCT 8.0.0). Companions ported here because they only exist as private
//! helpers of OCCT:
//!
//! - `PLib::CoefficientsPoles(dim, Coefs, NoWeights, Poles, NoWeights)`
//!   (`PLib.cxx:1522-1603`), the Bezier-to-power-basis conversion
//!   ([`coefficients_poles_dim`]),
//! - `BSplCLib::Eval` (`BSplCLib.cxx:3640-3790`) is
//!   [`occt_core::bspl::eval_extrap::eval_dim`],
//! - `BSplCLib::KnotSequence` (`cxx:488-516`) is
//!   [`occt_core::bspl::banded_interp::knot_sequence`],
//! - `BSplCLib::IncreaseDegree` / `RemoveKnot` are the
//!   `occt_core::bspl::increase_degree` / `occt_core::bspl::remove_knot`
//!   dimension-major ports.
//!
//! Divergences that are not heuristics:
//! - OCCT passes raw `double&` buffers that alias `NCollection_Array1`
//!   storage and writes `NbPolesResult` / `NbKnotsResult` through references;
//!   the port returns them in [`TangExtendResult`]. When the removal loop
//!   leaves `M != CDegree` OCCT copies from the `RemoveKnot` result arrays
//!   (`cxx:4215-4238`), which the port mirrors,
//! - the `#ifdef OCCT_DEBUG` degree check (`cxx:3806-3812`) only prints.

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::eval_extrap::eval_dim;
use occt_core::bspl::increase_degree::{increase_degree_count_knots, increase_degree_dim};
use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::bspl::plib::bin;
use occt_core::bspl::remove_knot::remove_knot_dim;
use occt_core::precision::CONFUSION;
use occt_math::hermite_coefficients;

/// The four outputs of `BSplCLib::TangExtendToConstraint`.
#[derive(Debug, Clone, PartialEq)]
pub struct TangExtendResult {
    /// `NbPolesResult`
    pub nb_poles: usize,
    /// `NbKnotsResult` — the flat knot count.
    pub nb_knots: usize,
    /// `PolesResult`, dimension-major (`nb_poles * dimension` values).
    pub poles: Vec<f64>,
    /// `KnotsResult`, the flat knot vector (`nb_knots` values).
    pub flat_knots: Vec<f64>,
}

/// `PLib::CoefficientsPoles(dim, Coefs, WCoefs = nullptr, Poles, Weights = nullptr)`
/// (`PLib.cxx:1522-1603`).
///
/// `coefs` and `poles` hold `reflen * dim` values in the same dimension-major
/// order; the weights branches of the source are not ported because
/// `TangExtendToConstraint` only calls the `NoWeights()` overload.
pub fn coefficients_poles_dim(dim: usize, coefs: &[f64], poles: &mut [f64]) {
    let dim = dim.max(1);
    let reflen = coefs.len() / dim;
    if reflen == 0 || poles.len() < reflen * dim {
        return;
    }
    // `cxx:1548-1553`: the two extremities.
    for i in 0..dim {
        poles[i] = coefs[i];
        poles[reflen * dim - 1 - i] = coefs[reflen * dim - 1 - i];
    }
    // `cxx:1560-1572`: the interior coefficients divided by the binomial.
    for i in 2..reflen {
        let cnp = bin((reflen - 1) as i32, (i - 1) as i32);
        for j in 0..dim {
            poles[dim * (i - 1) + j] = coefs[dim * (i - 1) + j] / cnp;
        }
    }
    // `cxx:1575-1589`: the accumulation towards the power basis.
    for i in 1..reflen {
        for j in (i..=reflen - 1).rev() {
            for k in 0..dim {
                poles[dim * j + k] += poles[dim * (j - 1) + k];
            }
        }
    }
}

/// `BSplCLib::TangExtendToConstraint(FlatKnots, C1Coefficient, NumPoles,
/// Poles, CDimension, CDegree, ConstraintPoint, Continuity, After,
/// NbPolesResult, NbKnotsResult, KnotsResult, PolesResult)`
/// (`BSplCLib.cxx:3794-4238`).
///
/// `poles` is the dimension-major pole array of the curve being extended
/// (`num_poles * dimension` values) and `flat_knots` its flat knot vector
/// (`num_poles + degree + 1` values), both exactly as OCCT receives them.
pub fn tang_extend_to_constraint(
    flat_knots: &[f64],
    c1_coefficient: f64,
    num_poles: usize,
    poles: &[f64],
    dimension: usize,
    degree: usize,
    constraint_point: &[f64],
    continuity: usize,
    after: bool,
) -> Option<TangExtendResult> {
    let dim = dimension.max(1);
    let csize = continuity + 2;

    // Hermite matrix (`cxx:3826-3843`).
    let mat_coefs = if after {
        hermite_coefficients(0.0, 1.0, continuity as i32, 0)?
    } else {
        hermite_coefficients(0.0, 1.0, 0, continuity as i32)?
    };

    // Position at the node of connection (`cxx:3846-3853`).
    let tbord = if after {
        *flat_knots.get(flat_knots.len().checked_sub(1 + degree)?)?
    } else {
        *flat_knots.get(degree)?
    };
    let derivative_request = continuity.max(1) as i32;
    let extrap_mode = [degree as i32; 2];

    // `cxx:3857-3868`.
    let eval_bs = eval_dim(
        tbord,
        false,
        derivative_request,
        extrap_mode,
        degree as i32,
        flat_knots,
        dim,
        poles,
    )?;
    let l1 = norm_of(&eval_bs, dim);

    // Matrix of constraints (`cxx:3880-3920`).
    let mut contraintes = vec![vec![0.0f64; dim]; csize];
    let c1_2 = c1_coefficient * c1_coefficient;
    let c1_3 = c1_2 * c1_coefficient;
    let eval_at = |order: usize, ipos: usize| -> f64 { eval_bs[order * dim + (ipos - 1)] };
    if after {
        for ipos in 1..=dim {
            contraintes[0][ipos - 1] = eval_at(0, ipos);
            contraintes[1][ipos - 1] = c1_coefficient * eval_at(1, ipos);
            if continuity >= 2 {
                contraintes[2][ipos - 1] = eval_at(2, ipos) * c1_2;
            }
            if continuity >= 3 {
                contraintes[3][ipos - 1] = eval_at(3, ipos) * c1_3;
            }
            contraintes[continuity + 1][ipos - 1] = constraint_point[ipos - 1];
        }
    } else {
        for ipos in 1..=dim {
            contraintes[0][ipos - 1] = constraint_point[ipos - 1];
            contraintes[1][ipos - 1] = eval_at(0, ipos);
            if continuity >= 1 {
                contraintes[2][ipos - 1] = c1_coefficient * eval_at(1, ipos);
            }
            if continuity >= 2 {
                contraintes[3][ipos - 1] = eval_at(2, ipos) * c1_2;
            }
            if continuity >= 3 {
                contraintes[4][ipos - 1] = eval_at(3, ipos) * c1_3;
            }
        }
    }

    // Coefficients of the extension (`cxx:3925-3938`).
    let mut extra_coeffs = vec![0.0; csize * dim];
    for ii in 1..=csize {
        for jj in 1..=csize {
            for kk in 1..=dim {
                extra_coeffs[kk - 1 + (jj - 1) * dim] +=
                    mat_coefs[ii - 1][jj - 1] * contraintes[ii - 1][kk - 1];
            }
        }
    }

    // Poles of the extension (`cxx:3941-3947`).
    let mut extrap_poles = vec![0.0; csize * dim];
    coefficients_poles_dim(dim, &extra_coeffs, &mut extrap_poles);

    // Nodes of the extension (`cxx:3950-3959`).
    let extrap_knots = [0.0, 1.0];
    let extrap_mults = [csize as i32, csize as i32];
    let fk2 = knot_sequence(&extrap_knots, &extrap_mults, (csize - 1) as i32);

    // Norm of the tangent at the connection point of the extension
    // (`cxx:3962-3977`). OCCT evaluates into the `EvalBS` buffer, so `Tgte`
    // reads the extension's derivative, not the surface's.
    let border = if after { 0.0 } else { 1.0 };
    let eval_ext = eval_dim(
        border,
        false,
        1,
        extrap_mode,
        (csize - 1) as i32,
        &fk2,
        dim,
        &extrap_poles,
    )?;
    let l2 = norm_of(&eval_ext, dim);

    // Harmonisation of degrees (`cxx:3979-4005`).
    let from_degree = (csize - 1) as i32;
    let to_degree = degree as i32;
    let (new_p2, new_k2, new_m2) = if from_degree < to_degree {
        let n_orig = csize;
        let from_k1 = first_u_knot_index(from_degree, &extrap_mults);
        let to_k2 = last_u_knot_index(from_degree, &extrap_mults);
        let step = to_degree - from_degree;
        let nb_poles = (n_orig as i32 + step * (to_k2 - from_k1)).max(0) as usize;
        let nb_knots = increase_degree_count_knots(
            from_degree,
            to_degree,
            false,
            &extrap_mults,
        )
        .max(0) as usize;
        let mut new_poles = vec![0.0; nb_poles * dim];
        let mut new_knots = vec![0.0; nb_knots];
        let mut new_mults = vec![0i32; nb_knots];
        increase_degree_dim(
            from_degree,
            to_degree,
            false,
            dim as i32,
            &extrap_poles,
            &extrap_knots,
            &extrap_mults,
            &mut new_poles,
            &mut new_knots,
            &mut new_mults,
        );
        (new_poles, new_knots, new_mults)
    } else {
        (extrap_poles.clone(), extrap_knots.to_vec(), extrap_mults.to_vec())
    };
    let new_fk2 = knot_sequence(&new_k2, &new_m2, degree as i32);

    // Ratio of reparametrisation (`cxx:4013-4034`).
    let mut ratio = 1.0;
    if l1 > CONFUSION && l2 > CONFUSION {
        ratio = l2 / l1;
    }
    if !(1.0e-5..=1.0e5).contains(&ratio) {
        ratio = 1.0;
    }
    let delta = if after {
        ratio * new_fk2.first().copied().unwrap_or(0.0) - flat_knots.last().copied().unwrap_or(0.0)
    } else {
        ratio * new_fk2.last().copied().unwrap_or(0.0) - flat_knots.first().copied().unwrap_or(0.0)
    };

    // Result of the concatenation (`cxx:4036-4085`).
    let nb_p1 = num_poles;
    let nb_p2 = degree + 1;
    let nb_k1 = flat_knots.len();
    let nb_k2 = 2 * (degree + 1);
    let mut new_poles = vec![0.0; (nb_p1 + nb_p2 - 1) * dim];
    let mut new_flats = vec![0.0; nb_k1 + nb_k2 - degree - 2];

    if after {
        for ii in 1..=nb_p1 + nb_p2 - 1 {
            for jj in 1..=dim {
                let ind_np = (ii - 1) * dim + jj;
                let ind_p = (ii - 1) * dim + jj - 1;
                let ind_ep = (ii - nb_p1) * dim + jj;
                new_poles[ind_np - 1] = if ii < nb_p1 {
                    // `Padr[indP]` reads `Poles(indP + 1)`; `ind_p` is that
                    // 0-based slot.
                    *poles.get(ind_p).unwrap_or(&0.0)
                } else {
                    *new_p2.get(ind_ep - 1).unwrap_or(&0.0)
                };
            }
        }
    } else {
        for ii in 1..=nb_p1 + nb_p2 - 1 {
            for jj in 1..=dim {
                let ind_np = (ii - 1) * dim + jj;
                let ind_ep = (ii - 1) * dim + jj;
                let ind_p = (ii - nb_p2) * dim + jj - 1;
                new_poles[ind_np - 1] = if ii < nb_p2 {
                    *new_p2.get(ind_ep - 1).unwrap_or(&0.0)
                } else {
                    *poles.get(ind_p).unwrap_or(&0.0)
                };
            }
        }
    }

    // Flat nodes (`cxx:4087-4117`).
    if after {
        for ii in 1..nb_k1 {
            new_flats[ii - 1] = flat_knots[ii - 1];
        }
        for ii in 1..=nb_k2 - degree - 1 {
            new_flats[nb_k1 + ii - 2] = ratio * new_fk2[ii + degree] - delta;
        }
    } else {
        for ii in 1..nb_k2 - degree {
            new_flats[ii - 1] = ratio * new_fk2[ii - 1] - delta;
        }
        for ii in 2..=nb_k1 {
            new_flats[nb_k2 + ii - degree - 3] = flat_knots[ii - 1];
        }
    }

    // Flat nodes -> nodes + multiplicities (`cxx:4125-4153`).
    let n_flat = nb_k1 + nb_k2 - degree - 2;
    let mut k_length = 1usize;
    for ii in 2..=n_flat {
        if new_flats[ii - 1] != new_flats[ii - 2] {
            k_length += 1;
        }
    }
    let mut new_knots = vec![0.0f64; k_length];
    let mut new_mults = vec![1i32; k_length];
    let mut jj = 1usize;
    new_knots[0] = new_flats[0];
    for ii in 2..=n_flat {
        if new_flats[ii - 1] == new_flats[ii - 2] {
            new_mults[jj - 1] += 1;
        } else {
            jj += 1;
            new_knots[jj - 1] = new_flats[ii - 1];
        }
    }

    // Reduction of multiplicity at the connection node (`cxx:4155-4187`).
    let index = if after { k_length - 1 } else { 2 };
    let mut m = degree as i32;
    let mut ok = true;
    let mut result_poles = vec![0.0; (nb_p1 + nb_p2 - 1) * dim];
    let mut result_knots = vec![0.0f64; k_length];
    let mut result_mults = vec![0i32; k_length];
    while m > degree as i32 - continuity as i32 && ok {
        ok = remove_knot_dim(
            index as i32,
            m - 1,
            degree as i32,
            false,
            dim as i32,
            &new_poles,
            &new_knots,
            &new_mults,
            &mut result_poles,
            &mut result_knots,
            &mut result_mults,
            1.0e-6,
        );
        if ok {
            m -= 1;
        }
    }

    if m == degree as i32 {
        // `cxx:4188-4212`: the poles and the flat knots of the concatenation.
        let nb_poles_result = nb_p1 + nb_p2 - 1;
        let mut flat = Vec::new();
        for (i, knot) in new_knots.iter().enumerate() {
            for _ in 0..new_mults[i] {
                flat.push(*knot);
            }
        }
        Some(TangExtendResult {
            nb_poles: nb_poles_result,
            nb_knots: flat.len(),
            poles: new_poles,
            flat_knots: flat,
        })
    } else {
        // `cxx:4214-4238`: the poles and the flat knots of the last removal.
        let nb_poles_result = (nb_p1 + nb_p2) as i32 - 1 - degree as i32 + m;
        let nb_poles_result = nb_poles_result.max(0) as usize;
        let mut flat = Vec::new();
        for (i, knot) in result_knots.iter().enumerate() {
            for _ in 0..result_mults[i] {
                flat.push(*knot);
            }
        }
        Some(TangExtendResult {
            nb_poles: nb_poles_result,
            nb_knots: flat.len(),
            poles: result_poles,
            flat_knots: flat,
        })
    }
}

/// Norm of the first-derivative block of an `Eval` result.
fn norm_of(results: &[f64], dim: usize) -> f64 {
    let start = dim;
    results
        .get(start..start + dim)
        .map(|block| block.iter().map(|v| v * v).sum::<f64>().sqrt())
        .unwrap_or(0.0)
}
