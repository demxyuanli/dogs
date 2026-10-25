//! `ApproxInt_KnotTools` — knot placement for the approximation of a discrete
//! point set (3D point plus one or two 2D parameter pairs).
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_KnotTools.cxx:1-848`
//! and `ApproxInt_KnotTools.hxx:38-132` (OCCT 8.0.0).
//!
//! Algorithm (`hxx:38-47`, after Razdan, "Knot Placement for B-Spline curve
//! Approximation"):
//! 1. discrete curvature on the point set (`BuildCurvature`, `cxx:88-167`),
//! 2. draft knot indices from its extrema and monotone intervals
//!    (`ComputeKnotInds`, `cxx:173-327`, `InsKnotBefI`, `cxx:452-552`),
//! 3. filtration of that draft (`FilterKnots`, `cxx:332-448`).
//!
//! Translation notes:
//! - `NCollection_Sequence<int>` (`theInds`, 1-based, grown by `InsertBefore`)
//!   becomes a `Vec<i32>` addressed through the 1-based helpers below so the
//!   OCCT index arithmetic stays visible.
//! - `NCollection_LocalArray<double> theCoords` is the flat, point-major
//!   coordinate table (`dim` doubles per point, `cxx:589-616`).
//! - `NCollection_Array1<double> theCurv` keeps its `Lower()..Upper()`
//!   indexing as a [`MathVector`] over the same range as `thePars`.
//! - `PLib::EvalLagrange(..., Dim, ...)` is
//!   [`eval_lagrange_nd`](occt_core::bspl::plib_eval::eval_lagrange_nd).
//!   OCCT ignores its `ReturnCode` and reads an uninitialised `Results` array
//!   on failure (`cxx:118`, `:139`, `:160`); the port uses zeros there, which
//!   is the only value OCCT's caller can be relying on.
//! - `DefineParType` (`cxx:669-848`) is UNPORTED: it builds a
//!   `GeomInt_TheMultiLineOfWLApprox` over an `IntPatch_WLine` (`cxx:686-699`),
//!   fills it through that line's `Value` (`cxx:707-744`) and parametrises it
//!   with `GeomInt_WLApprox::Parameters` (`cxx:794`); none of those are ported
//!   yet. It is not on the `ApproxInt_Approx` path either — that driver only
//!   calls `BuildKnots` (`ApproxInt_Approx.gxx:610`).

use occt_core::bspl::plib_eval::eval_lagrange_nd;
use occt_core::precision::{epsilon, Precision};
use occt_math::MathVector;

/// `aSinCoeff2 = 0.09549150281252627` (`cxx:28`): `(3 - sqrt(5)) / 8`, the
/// squared sine of the angle the draft-knot criteria test against.
const A_SIN_COEFF2: f64 = 0.09549150281252627;

/// `aMaxPntCoeff = 15` (`cxx:29`): the largest knot interval, as a multiple of
/// `theMinNbPnts`.
const A_MAX_PNT_COEFF: i32 = 15;

// ---- `NCollection_Sequence<int>` accessors, 1-based as in OCCT -------------

/// `theInds(i)` (`1 <= i <= Length()`).
fn seq_at(v: &[i32], i: i32) -> i32 {
    v[(i - 1) as usize]
}

/// `theInds.Length()`.
fn seq_len(v: &[i32]) -> i32 {
    v.len() as i32
}

/// `theInds.Last()`.
fn seq_last(v: &[i32]) -> i32 {
    *v.last().expect("empty NCollection_Sequence")
}

/// `theInds.Append(value)`.
fn seq_append(v: &mut Vec<i32>, value: i32) {
    v.push(value);
}

/// `theInds.InsertBefore(i, value)`.
fn seq_insert_before(v: &mut Vec<i32>, i: i32, value: i32) {
    v.insert((i - 1) as usize, value);
}

/// `theInds.InsertAfter(i, value)`.
fn seq_insert_after(v: &mut Vec<i32>, i: i32, value: i32) {
    v.insert(i as usize, value);
}

/// `EvalCurv` (`cxx:36-85`): `|V1^V2| / |V1|^3` for `dim`-dimensional vectors.
///
/// `V1` is the first and `V2` the second derivative of the curve; `V1^V2` is
/// their outer product. A vanishing `V1` answers `Precision::Infinite()` so the
/// caller inserts a knot at the singularity (`cxx:63-77`).
fn eval_curv(dim: i32, v1: &[f64], v2: &[f64]) -> f64 {
    let dim = dim.max(0) as usize;
    let mut mp = 0.0f64;
    for i in 1..dim {
        for j in 0..i {
            let p = v1[i] * v2[j] - v1[j] * v2[i];
            mp += p * p;
        }
    }
    let mut q = 0.0f64;
    for i in 0..dim {
        q += v1[i] * v1[i];
    }
    if q < 1.0 / Precision::INFINITE {
        return Precision::INFINITE;
    }
    q = q.min(Precision::INFINITE);
    q *= q * q;
    (mp / q).sqrt()
}

/// `BuildCurvature` (`cxx:88-167`): discrete curvature at every point.
///
/// The interior points evaluate the Lagrange polynomial through their two
/// neighbours and themselves at their own parameter; the two end points use the
/// three points they can see (`cxx:104-120`, `cxx:145-161`).
pub fn build_curvature(
    coords: &[f64],
    dim: i32,
    pars: &MathVector,
    curv: &mut MathVector,
    max_curv: &mut f64,
) {
    // `Val` / `Res` hold the three points and their three derivatives, one
    // block of `dim` per point (`cxx:94`).
    let mut val = [0.0f64; 21];
    let mut par = [0.0f64; 3];

    let lower = curv.lower() as i32;
    let upper = curv.upper() as i32;
    let dim_usize = dim.max(0) as usize;
    *max_curv = 0.0;

    if curv.len() < 3 {
        curv.init(0.0);
        return;
    }

    let mut i = lower;
    for j in 0..3 {
        let k = i + j;
        let ic = (k - lower) * dim;
        let l = dim * j;
        for m in 0..dim {
            val[(l + m) as usize] = coords[(ic + m) as usize];
        }
        par[j as usize] = pars.value(k as usize);
    }
    let res = eval_lagrange_nd(par[0], 2, 2, dim_usize, &val, &par)
        .unwrap_or_else(|_| vec![0.0f64; 3 * dim_usize]);
    curv.set_value(
        i as usize,
        eval_curv(
            dim,
            &res[dim_usize..2 * dim_usize],
            &res[2 * dim_usize..3 * dim_usize],
        ),
    );
    if curv.value(i as usize) > *max_curv {
        *max_curv = curv.value(i as usize);
    }

    for i in (lower + 1)..upper {
        for j in 0..3 {
            let k = i + j - 1;
            let ic = (k - lower) * dim;
            let l = dim * j;
            for m in 0..dim {
                val[(l + m) as usize] = coords[(ic + m) as usize];
            }
            par[j as usize] = pars.value(k as usize);
        }
        let res = eval_lagrange_nd(par[1], 2, 2, dim_usize, &val, &par)
            .unwrap_or_else(|_| vec![0.0f64; 3 * dim_usize]);
        curv.set_value(
            i as usize,
            eval_curv(
                dim,
                &res[dim_usize..2 * dim_usize],
                &res[2 * dim_usize..3 * dim_usize],
            ),
        );
        if curv.value(i as usize) > *max_curv {
            *max_curv = curv.value(i as usize);
        }
    }

    i = upper;
    for j in 0..3 {
        let k = i + j - 2;
        let ic = (k - lower) * dim;
        let l = dim * j;
        for m in 0..dim {
            val[(l + m) as usize] = coords[(ic + m) as usize];
        }
        par[j as usize] = pars.value(k as usize);
    }
    let res = eval_lagrange_nd(par[2], 2, 2, dim_usize, &val, &par)
        .unwrap_or_else(|_| vec![0.0f64; 3 * dim_usize]);
    curv.set_value(
        i as usize,
        eval_curv(
            dim,
            &res[dim_usize..2 * dim_usize],
            &res[2 * dim_usize..3 * dim_usize],
        ),
    );
    if curv.value(i as usize) > *max_curv {
        *max_curv = curv.value(i as usize);
    }
}

/// `ComputeKnotInds` (`cxx:173-327`): the draft knot indices, 1-based into
/// `theCoords` / `thePars`.
pub fn compute_knot_inds(coords: &[f64], dim: i32, pars: &MathVector, inds: &mut Vec<i32>) {
    let mut a_feature_inds: Vec<i32> = Vec::new();
    let mut a_curv = MathVector::new(pars.lower(), pars.upper());
    let mut a_max_curv = 0.0f64;
    build_curvature(coords, dim, pars, &mut a_curv, &mut a_max_curv);

    let lower = a_curv.lower() as i32;
    let upper = a_curv.upper() as i32;

    seq_append(inds, lower);
    if a_max_curv <= Precision::CONFUSION {
        // Linear case (`cxx:193-198`).
        seq_append(inds, upper);
        return;
    }

    // II: extremas of curvature (`cxx:200-226`).
    let eps = 1.0e-9f64;
    let eps1 = 1.0e3 * eps;
    for i in (lower + 1)..upper {
        let d1 = a_curv.value(i as usize) - a_curv.value((i - 1) as usize);
        let d2 = a_curv.value(i as usize) - a_curv.value((i + 1) as usize);
        let ad1 = d1.abs();
        let ad2 = d2.abs();
        if d1 * d2 > 0.0 && ad1 > eps && ad2 > eps {
            if i != seq_last(inds) {
                seq_append(inds, i);
                seq_append(&mut a_feature_inds, i);
            }
        } else if (ad1 < eps && ad2 > eps1) || (ad1 > eps1 && ad2 < eps) {
            if i != seq_last(inds) {
                seq_append(inds, i);
                seq_append(&mut a_feature_inds, i);
            }
        }
    }
    if upper != seq_last(inds) {
        seq_append(inds, upper);
    }

    // III: knots in the monotone intervals of curvature (`cxx:240-253`).
    let mut i = 1i32;
    loop {
        i += 1;
        let ok = ins_knot_bef_i(i, &a_curv, coords, dim, inds, true);
        if ok {
            i -= 1;
        }
        if i >= seq_len(inds) {
            break;
        }
    }

    // IV: checking feature points (`cxx:256-326`).
    let mut j = 2i32;
    for idx in 1..=seq_len(&a_feature_inds) {
        let an_ind = seq_at(&a_feature_inds, idx);
        while j <= seq_len(inds) - 1 {
            if seq_at(inds, j) == an_ind {
                let an_ind_prev = seq_at(inds, j - 1);
                let an_ind_next = seq_at(inds, j + 1);
                let ici = (an_ind_prev - lower) * dim;
                let ici1 = (an_ind_next - lower) * dim;
                let icm = (an_ind - lower) * dim;
                let mut v1 = vec![0.0f64; dim.max(0) as usize];
                let mut v2 = vec![0.0f64; dim.max(0) as usize];
                let mut mp = 0.0f64;
                let mut m1 = 0.0f64;
                let mut m2 = 0.0f64;
                for k in 0..dim {
                    v1[k as usize] = coords[(icm + k) as usize] - coords[(ici + k) as usize];
                    m1 += v1[k as usize] * v1[k as usize];
                    v2[k as usize] = coords[(ici1 + k) as usize] - coords[(icm + k) as usize];
                    m2 += v2[k as usize] * v2[k as usize];
                }
                for k in 1..dim {
                    for l in 0..k {
                        let p = v1[k as usize] * v2[l as usize] - v1[l as usize] * v2[k as usize];
                        mp += p * p;
                    }
                }
                if mp > A_SIN_COEFF2 * m1 * m2 {
                    // Insert new knots (`cxx:291-313`).
                    let d1 = (a_curv.value(an_ind as usize) - a_curv.value(an_ind_prev as usize)).abs();
                    let d2 = (a_curv.value(an_ind as usize) - a_curv.value(an_ind_next as usize)).abs();
                    if d1 > d2 {
                        let ok = ins_knot_bef_i(j, &a_curv, coords, dim, inds, false);
                        if ok {
                            j += 1;
                        } else {
                            break;
                        }
                    } else {
                        let ok = ins_knot_bef_i(j + 1, &a_curv, coords, dim, inds, false);
                        if !ok {
                            break;
                        }
                    }
                } else {
                    j += 1;
                    break;
                }
            } else {
                j += 1;
            }
        }
    }
}

/// `FilterKnots` (`cxx:332-448`): drops intervals that are too dense or too
/// sparse and fills the last knot.
pub fn filter_knots(inds: &mut Vec<i32>, min_nb_pnts: i32, l_knots: &mut Vec<i32>) {
    // Maximum number of points per knot interval (`cxx:335-337`).
    let a_max_nb_pnts = A_MAX_PNT_COEFF * min_nb_pnts;
    let mut i = 1i32;
    let a_min_nb_step = min_nb_pnts / 2;

    // I: filter too big number of points per knot interval (`cxx:340-354`).
    while i < seq_len(inds) {
        let nbint = seq_at(inds, i + 1) - seq_at(inds, i) + 1;
        if nbint <= a_max_nb_pnts {
            i += 1;
            continue;
        } else {
            let ind = seq_at(inds, i) + nbint / 2;
            seq_insert_after(inds, i, ind);
        }
    }

    // II: filter points with too small amount of points per knot interval
    // (`cxx:357-434`).
    i = 1;
    l_knots.push(seq_at(inds, i));
    let mut an_inds_prev = seq_at(inds, i);
    i = 2;
    while i <= seq_len(inds) {
        if seq_at(inds, i) - an_inds_prev <= min_nb_pnts {
            if i != seq_len(inds) {
                let mut an_idx = i + 1;
                while an_idx <= seq_len(inds) {
                    if seq_at(inds, an_idx) - an_inds_prev >= min_nb_pnts {
                        break;
                    }
                    an_idx += 1;
                }
                an_idx -= 1;

                let a_mid_idx = (seq_at(inds, an_idx) + an_inds_prev) / 2;
                if a_mid_idx - an_inds_prev < min_nb_pnts
                    && a_mid_idx - seq_at(inds, an_idx) < min_nb_pnts
                    && seq_at(inds, an_idx) - an_inds_prev >= a_min_nb_step
                {
                    if seq_at(inds, an_idx) - an_inds_prev > 2 * min_nb_pnts {
                        // Bad distribution points merge into one knot interval.
                        l_knots.push(an_inds_prev + min_nb_pnts);
                        an_inds_prev += min_nb_pnts;
                        i = an_idx - 1;
                    } else if seq_at(inds, an_idx - 1) - an_inds_prev >= min_nb_pnts / 2 {
                        // Bad distribution points merge into one knot interval.
                        l_knots.push(seq_at(inds, an_idx - 1));
                        an_inds_prev = seq_at(inds, an_idx - 1);
                        i = an_idx - 1;
                        if seq_at(inds, an_idx) - seq_at(inds, an_idx - 1) <= min_nb_pnts / 2 {
                            *l_knots.last_mut().expect("empty knot list") = seq_at(inds, an_idx);
                            an_inds_prev = seq_at(inds, an_idx);
                            i = an_idx;
                        }
                    } else {
                        // Bad distribution points merge into one knot interval.
                        l_knots.push(seq_at(inds, an_idx));
                        an_inds_prev = seq_at(inds, an_idx);
                        i = an_idx;
                    }
                } else if an_idx == seq_len(inds) && l_knots.len() >= 2 {
                    // Last point obtained; the previous knot may still be moved
                    // (`cxx:411-425`).
                    let a_last_good_idx = l_knots[l_knots.len() - 2];
                    if seq_last(inds) - 2 * min_nb_pnts >= a_last_good_idx {
                        *l_knots.last_mut().expect("empty knot list") =
                            seq_last(inds) - min_nb_pnts;
                        l_knots.push(seq_last(inds));
                        an_inds_prev = seq_at(inds, an_idx);
                        i = an_idx;
                    }
                }
            }
        } else {
            l_knots.push(seq_at(inds, i));
            an_inds_prev = seq_at(inds, i);
        }
        i += 1;
    }

    // III: fill last knot (`cxx:436-447`).
    if l_knots.len() < 2 {
        l_knots.push(seq_last(inds));
    } else if *l_knots.last().expect("empty knot list") < seq_last(inds) {
        *l_knots.last_mut().expect("empty knot list") = seq_last(inds);
    }
}

/// `InsKnotBefI` (`cxx:452-552`): inserts a knot before position `the_i`.
///
/// With `chk_curv` the insertion additionally requires the turn angle at the
/// midpoint to exceed `aSinCoeff` (`cxx:536-541`).
fn ins_knot_bef_i(
    the_i: i32,
    the_curv: &MathVector,
    the_coords: &[f64],
    the_dim: i32,
    the_inds: &mut Vec<i32>,
    chk_curv: bool,
) -> bool {
    let an_ind1 = seq_at(the_inds, the_i);
    let an_ind = seq_at(the_inds, the_i - 1);
    if an_ind1 - an_ind == 1 {
        return false;
    }

    let curv = 0.5 * (the_curv.value(an_ind as usize) + the_curv.value(an_ind1 as usize));
    let a_limit_curvature_change = 3.0f64;
    let curv_lower = the_curv.lower() as i32;

    for j in (an_ind + 1)..an_ind1 {
        let mut mid = 0i32;

        // I: curvature change criteria (`cxx:472-486`).
        if the_curv.value(j as usize) > Precision::CONFUSION
            && the_curv.value(an_ind as usize) > Precision::CONFUSION
        {
            let ratio = the_curv.value(j as usize) / the_curv.value(an_ind as usize);
            if ratio > a_limit_curvature_change || ratio < 1.0 / a_limit_curvature_change {
                mid = j;
                seq_insert_before(the_inds, the_i, mid);
                return true;
            }
        }

        // II: angular criteria (`cxx:488-505`).
        let ac = the_curv.value((j - 1) as usize);
        let ac1 = the_curv.value(j as usize);
        if (curv >= ac && curv <= ac1) || (curv >= ac1 && curv <= ac) {
            if (curv - ac).abs() < (curv - ac1).abs() {
                mid = j - 1;
            } else {
                mid = j;
            }
        }
        if mid == an_ind {
            mid += 1;
        }
        if mid == an_ind1 {
            mid -= 1;
        }
        if mid > 0 {
            if chk_curv {
                let ici = (an_ind - curv_lower) * the_dim;
                let ici1 = (an_ind1 - curv_lower) * the_dim;
                let icm = (mid - curv_lower) * the_dim;
                let mut v1 = vec![0.0f64; the_dim.max(0) as usize];
                let mut v2 = vec![0.0f64; the_dim.max(0) as usize];
                let mut mp = 0.0f64;
                let mut m1 = 0.0f64;
                let mut m2 = 0.0f64;
                for i in 0..the_dim {
                    v1[i as usize] =
                        the_coords[(icm + i) as usize] - the_coords[(ici + i) as usize];
                    m1 += v1[i as usize] * v1[i as usize];
                    v2[i as usize] =
                        the_coords[(ici1 + i) as usize] - the_coords[(icm + i) as usize];
                    m2 += v2[i as usize] * v2[i as usize];
                }
                for i in 1..the_dim {
                    for jj in 0..i {
                        let p = v1[i as usize] * v2[jj as usize] - v1[jj as usize] * v2[i as usize];
                        mp += p * p;
                    }
                }
                if mp > A_SIN_COEFF2 * m1 * m2 {
                    seq_insert_before(the_inds, the_i, mid);
                    return true;
                }
            } else {
                seq_insert_before(the_inds, the_i, mid);
                return true;
            }
        }
    }

    false
}

/// `BuildKnots` (`cxx:556-641`): the knot sequence of the point set.
///
/// At least one of the three point sets must be present; those the flags name
/// are concatenated into the same flat coordinate table, so the curvature is
/// computed on the combined `dim`-dimensional curve.
#[allow(clippy::too_many_arguments)]
pub fn build_knots(
    the_pnts_xyz: &[occt_core::gp::GpPnt],
    the_pnts_u1v1: &[occt_core::gp::GpPnt2d],
    the_pnts_u2v2: &[occt_core::gp::GpPnt2d],
    the_pars: &MathVector,
    the_approx_xyz: bool,
    the_approx_u1v1: bool,
    the_approx_u2v2: bool,
    the_min_nb_pnts: i32,
    the_knots: &mut Vec<i32>,
) {
    let mut a_knots: Vec<i32> = Vec::new();
    let mut a_dim = 0i32;

    // I: convert input data to the corresponding format (`cxx:569-617`).
    if the_approx_xyz {
        a_dim += 3;
    }
    if the_approx_u1v1 {
        a_dim += 2;
    }
    if the_approx_u2v2 {
        a_dim += 2;
    }

    if a_dim == 0 || the_pars.len() < 2 {
        the_knots.clear();
        return;
    }

    let lower = the_pars.lower() as i32;
    let upper = the_pars.upper() as i32;
    let mut a_coords = vec![0.0f64; the_pars.len() * a_dim as usize];
    for i in lower..=upper {
        let mut j = (i - lower) * a_dim;
        // The arrays are 1-based in OCCT and share `thePars`' range; the port
        // keeps them as slices indexed from `lower`.
        let idx = (i - lower) as usize;
        if the_approx_xyz {
            a_coords[j as usize] = the_pnts_xyz[idx].x();
            j += 1;
            a_coords[j as usize] = the_pnts_xyz[idx].y();
            j += 1;
            a_coords[j as usize] = the_pnts_xyz[idx].z();
            j += 1;
        }
        if the_approx_u1v1 {
            a_coords[j as usize] = the_pnts_u1v1[idx].x();
            j += 1;
            a_coords[j as usize] = the_pnts_u1v1[idx].y();
            j += 1;
        }
        if the_approx_u2v2 {
            a_coords[j as usize] = the_pnts_u2v2[idx].x();
            j += 1;
            a_coords[j as usize] = the_pnts_u2v2[idx].y();
        }
    }

    // II: build draft knot sequence (`cxx:620`).
    compute_knot_inds(&a_coords, a_dim, the_pars, &mut a_knots);

    // III: build output knot sequence (`cxx:631`).
    filter_knots(&mut a_knots, the_min_nb_pnts, the_knots);
}

/// `Epsilon(1.)` (`Standard_Real.hxx:242-246`), used by the UNPORTED
/// `DefineParType`; kept so the constant stays named with its OCCT source.
#[allow(dead_code)]
fn define_par_type_epsilon() -> f64 {
    epsilon(1.0)
}
