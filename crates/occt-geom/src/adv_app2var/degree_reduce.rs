//! Degree reduction mma2fx6_ (AdvApp2Var_ApproxF2var.cxx:7285-7572)
//! and the Fortran-to-C++ transpose (AdvApp2Var_Patch.cxx:942-957).

pub(super) fn fx6(line: i32) -> String {
    format!("AdvApp2Var_ApproxF2var::mma2fx6_: AdvApp2Var_ApproxF2var.cxx:{line}")
}

pub(super) fn hpaux_len(ndjacu: i32, ndjacv: i32, ndimen: i32) -> Result<usize, String> {
    let oob = || "AdvApp2Var_Patch::MakeApprox: AdvApp2Var_Patch.cxx:834".to_string();
    let n = (1 + ndjacu).checked_mul(1 + ndjacv).ok_or_else(oob)?;
    let n = n.checked_mul(ndimen).ok_or_else(oob)?;
    let n = n.checked_mul(2).ok_or_else(oob)?;
    usize::try_from(n).map_err(|_| oob())
}

/// `HPAUX` layout (`AdvApp2Var_Patch.cxx:942-957`): dimension, then V, then U.
pub(super) fn transpose_can(
    ndimen: i32,
    ncfmxu: i32,
    ncfmxv: i32,
    patcan: &[f64],
    equation: &mut [f64],
) -> Result<(), String> {
    if ndimen < 1 || ncfmxu < 1 || ncfmxv < 1 {
        return Ok(());
    }
    for dim in 1..=ndimen {
        for ii in 1..=ncfmxu {
            for jj in 1..=ncfmxv {
                let src = can_src(ii, jj, dim, ncfmxu, ncfmxv)?;
                let dst = cpp_dst(dim, jj, ii, ndimen, ncfmxv)?;
                let value = *patcan.get(src).ok_or_else(transpose_oob)?;
                let slot = equation.get_mut(dst).ok_or_else(transpose_oob)?;
                *slot = value;
            }
        }
    }
    Ok(())
}

fn transpose_oob() -> String {
    "AdvApp2Var_Patch::MakeApprox: AdvApp2Var_Patch.cxx:952".into()
}

fn can_src(ii: i32, jj: i32, dim: i32, ncfmxu: i32, ncfmxv: i32) -> Result<usize, String> {
    let plane = ncfmxu.checked_mul(ncfmxv).ok_or_else(transpose_oob)?;
    let at = (ii - 1)
        .checked_add((jj - 1).checked_mul(ncfmxu).ok_or_else(transpose_oob)?)
        .ok_or_else(transpose_oob)?
        .checked_add((dim - 1).checked_mul(plane).ok_or_else(transpose_oob)?)
        .ok_or_else(transpose_oob)?;
    usize::try_from(at).map_err(|_| transpose_oob())
}

fn cpp_dst(dim: i32, jj: i32, ii: i32, ndimen: i32, ncfmxv: i32) -> Result<usize, String> {
    let at = (dim - 1)
        .checked_add((jj - 1).checked_mul(ndimen).ok_or_else(transpose_oob)?)
        .ok_or_else(transpose_oob)?
        .checked_add(
            (ii - 1)
                .checked_mul(ndimen)
                .ok_or_else(transpose_oob)?
                .checked_mul(ncfmxv)
                .ok_or_else(transpose_oob)?,
        )
        .ok_or_else(transpose_oob)?;
    usize::try_from(at).map_err(|_| transpose_oob())
}

/// `mma2fx6_` (`AdvApp2Var_ApproxF2var.cxx:7285-7572`). Coefficient counts only.
pub(super) fn mma2fx6(
    ncfmxu: i32,
    ncfmxv: i32,
    ndimen: i32,
    nbsesp: i32,
    ndimse: &[i32],
    nbupat: i32,
    nbvpat: i32,
    iordru: i32,
    iordrv: i32,
    epsapr: &[f64],
    epsfro: &[f64],
    patcan: &[f64],
    errmax: &[f64],
    ncoefu: &mut [i32],
    ncoefv: &mut [i32],
) -> Result<(), String> {
    if nbvpat < 1 {
        return Ok(());
    }
    for jj in 1..=nbvpat {
        if nbupat < 1 {
            continue;
        }
        for ii in 1..=nbupat {
            let slot = coef_index(ii, jj, nbupat)?;
            let mut ncfu = *ncoefu.get(slot).ok_or_else(|| fx6(7408))?;
            let mut ncfv = *ncoefv.get(slot).ok_or_else(|| fx6(7409))?;
            loop {
                if !(ncfu <= ((iordru + 1) << 1) && ncfu > 2) {
                    break;
                }
                if !subspaces_ok(
                    true, ncfu, ncfv, ii, jj, ncfmxu, ncfmxv, ndimen, nbsesp, ndimse, nbupat,
                    nbvpat, epsapr, epsfro, patcan, errmax, 7460, 7470,
                )? {
                    break;
                }
                ncfu -= 1;
            }
            loop {
                if !(ncfv <= ((iordrv + 1) << 1) && ncfv > 2) {
                    break;
                }
                if !subspaces_ok(
                    false, ncfu, ncfv, ii, jj, ncfmxu, ncfmxv, ndimen, nbsesp, ndimse, nbupat,
                    nbvpat, epsapr, epsfro, patcan, errmax, 7531, 7541,
                )? {
                    break;
                }
                ncfv -= 1;
            }
            *ncoefu.get_mut(slot).ok_or_else(|| fx6(7556))? = ncfu.max(2);
            *ncoefv.get_mut(slot).ok_or_else(|| fx6(7557))? = ncfv.max(2);
        }
    }
    Ok(())
}

fn subspaces_ok(
    along_u: bool,
    ncfu: i32,
    ncfv: i32,
    ii: i32,
    jj: i32,
    ncfmxu: i32,
    ncfmxv: i32,
    ndimen: i32,
    nbsesp: i32,
    ndimse: &[i32],
    nbupat: i32,
    nbvpat: i32,
    epsapr: &[f64],
    epsfro: &[f64],
    patcan: &[f64],
    errmax: &[f64],
    bid_line: i32,
    err_line: i32,
) -> Result<bool, String> {
    let mut idim = 0i32;
    if nbsesp >= 1 {
        for ns in 1..=nbsesp {
            let mut tol = *epsapr.get((ns - 1) as usize).ok_or_else(|| fx6(7423))?;
            for k in [9, 10, 11, 12] {
                tol = tol.min(*epsfro.get(fro_index(ns, k, nbsesp)?).ok_or_else(|| fx6(7425))?);
            }
            if ii == 1 || ii == nbupat || jj == 1 || jj == nbvpat {
                for k in [5, 6, 7, 8] {
                    tol = tol.min(*epsfro.get(fro_index(ns, k, nbsesp)?).ok_or_else(|| fx6(7439))?);
                }
            }
            let mut bid = 0.0;
            let width = *ndimse.get((ns - 1) as usize).ok_or_else(|| fx6(7453))?;
            if width >= 1 {
                for nd in 1..=width {
                    let id = idim + nd;
                    if along_u {
                        if ncfv >= 1 {
                            for kv in 1..=ncfv {
                                let at = pat_index(
                                    ncfu, kv, id, ii, jj, ncfmxu, ncfmxv, ndimen, nbupat, bid_line,
                                )?;
                                bid += patcan.get(at).ok_or_else(|| fx6(bid_line))?.abs();
                            }
                        }
                    } else if ncfu >= 1 {
                        for ku in 1..=ncfu {
                            let at = pat_index(
                                ku, ncfv, id, ii, jj, ncfmxu, ncfmxv, ndimen, nbupat, bid_line,
                            )?;
                            bid += patcan.get(at).ok_or_else(|| fx6(bid_line))?.abs();
                        }
                    }
                }
            }
            let emax = *errmax
                .get(err_index(ns, ii, jj, nbsesp, nbupat, err_line)?)
                .ok_or_else(|| fx6(err_line))?;
            if bid > tol * 1.0e-6 || bid > emax {
                return Ok(false);
            }
            idim = idim.checked_add(width).ok_or_else(|| fx6(7474))?;
        }
    }
    Ok(true)
}

fn pat_index(
    u: i32,
    v: i32,
    id: i32,
    ii: i32,
    jj: i32,
    ncfmxu: i32,
    ncfmxv: i32,
    ndimen: i32,
    nbupat: i32,
    line: i32,
) -> Result<usize, String> {
    let plane = ncfmxu.checked_mul(ncfmxv).ok_or_else(|| fx6(line))?;
    let base = (u - 1)
        .checked_add((v - 1).checked_mul(ncfmxu).ok_or_else(|| fx6(line))?)
        .ok_or_else(|| fx6(line))?
        .checked_add((id - 1).checked_mul(plane).ok_or_else(|| fx6(line))?)
        .ok_or_else(|| fx6(line))?;
    let patch_i = (ii - 1)
        .checked_add((jj - 1).checked_mul(nbupat).ok_or_else(|| fx6(line))?)
        .ok_or_else(|| fx6(line))?;
    let patch = patch_i
        .checked_mul(ndimen)
        .ok_or_else(|| fx6(line))?
        .checked_mul(plane)
        .ok_or_else(|| fx6(line))?;
    usize::try_from(base.checked_add(patch).ok_or_else(|| fx6(line))?).map_err(|_| fx6(line))
}

fn fro_index(ns: i32, k: i32, nbsesp: i32) -> Result<usize, String> {
    let at = (ns - 1)
        .checked_add(nbsesp.checked_mul(k - 5).ok_or_else(|| fx6(7425))?)
        .ok_or_else(|| fx6(7425))?;
    usize::try_from(at).map_err(|_| fx6(7425))
}

fn err_index(ns: i32, ii: i32, jj: i32, nbsesp: i32, nbupat: i32, line: i32) -> Result<usize, String> {
    let at = (ns - 1)
        .checked_add((ii - 1).checked_mul(nbsesp).ok_or_else(|| fx6(line))?)
        .ok_or_else(|| fx6(line))?
        .checked_add(
            (jj - 1)
                .checked_mul(nbupat)
                .ok_or_else(|| fx6(line))?
                .checked_mul(nbsesp)
                .ok_or_else(|| fx6(line))?,
        )
        .ok_or_else(|| fx6(line))?;
    usize::try_from(at).map_err(|_| fx6(line))
}

fn coef_index(ii: i32, jj: i32, nbupat: i32) -> Result<usize, String> {
    let at = (ii - 1)
        .checked_add((jj - 1).checked_mul(nbupat).ok_or_else(|| fx6(7408))?)
        .ok_or_else(|| fx6(7408))?;
    usize::try_from(at).map_err(|_| fx6(7408))
}

