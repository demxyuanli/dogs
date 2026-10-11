//! Degree error of one square (mma2er1_, mma2er2_, mma2moy_).

fn er1_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2er1_: AdvApp2Var_ApproxF2var.cxx:6342".into()
}

/// `mma2er1_` (`AdvApp2Var_ApproxF2var.cxx:6209-6370`).
pub(super) fn mma2er1(
    ndjacu: i32,
    ndjacv: i32,
    ndimen: i32,
    mindgu: i32,
    maxdgu: i32,
    mindgv: i32,
    maxdgv: i32,
    iordru: i32,
    iordrv: i32,
    xmaxju: &[f64],
    xmaxjv: &[f64],
    patjac: &[f64],
    erreur: f64,
) -> Result<f64, String> {
    let minu = (iordru + 1) << 1;
    let minv = (iordrv + 1) << 1;
    let pat_dim1 = ndjacu + 1;
    let pat_dim2 = ndjacv + 1;
    let pat = |ii: i32, jj: i32, nd: i32| {
        ii + jj * pat_dim1 + (nd - 1) * pat_dim1 * pat_dim2
    };
    let mut work = vec![0.0; ndimen.max(0) as usize];
    for nd in 1..=ndimen {
        let mut bid1 = 0.0;
        if mindgv <= maxdgv {
            for jj in mindgv..=maxdgv {
                let mut bid0 = 0.0;
                if mindgu <= maxdgu {
                    for ii in mindgu..=maxdgu {
                        let coeff = patjac
                            .get(pat(ii, jj, nd) as usize)
                            .copied()
                            .ok_or_else(er1_oob)?;
                        let scale = xmaxju.get((ii - minu) as usize).copied().ok_or_else(er1_oob)?;
                        bid0 += coeff.abs() * scale;
                    }
                }
                let scale = xmaxjv.get((jj - minv) as usize).copied().ok_or_else(er1_oob)?;
                bid1 += bid0 * scale;
            }
        }
        let slot = work.get_mut((nd - 1) as usize).ok_or_else(er1_oob)?;
        *slot = bid1;
    }
    let bid1 = super::border_error::mzsnorm(&work);
    Ok(super::border_error::mzsnorm(&[erreur, bid1]))
}

fn pat_at(ndjacu: i32, ndjacv: i32, ii: i32, jj: i32, nd: i32) -> i32 {
    let dim1 = ndjacu + 1;
    let dim2 = ndjacv + 1;
    ii + jj * dim1 + (nd - 1) * dim1 * dim2
}

fn er2_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2er2_: AdvApp2Var_ApproxF2var.cxx:6525".into()
}

fn moy_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2moy_: AdvApp2Var_ApproxF2var.cxx:7880".into()
}

/// `mma2er2_` (`AdvApp2Var_ApproxF2var.cxx:6374-6620`).
fn mma2er2(
    ndjacu: i32,
    ndjacv: i32,
    ndimen: i32,
    mindgu: i32,
    maxdgu: i32,
    mindgv: i32,
    maxdgv: i32,
    iordru: i32,
    iordrv: i32,
    xmaxju: &[f64],
    xmaxjv: &[f64],
    patjac: &[f64],
    epmscut: f64,
    mut erreur: f64,
) -> Result<(i32, i32, f64), String> {
    let i2rdu = (iordru + 1) << 1;
    let i2rdv = (iordrv + 1) << 1;
    let mut nu = maxdgu;
    let mut nv = maxdgv;
    let mut work = vec![0.0; ndimen.max(0) as usize];
    let coeff = |ii: i32, jj: i32, nd: i32| -> Result<f64, String> {
        patjac
            .get(pat_at(ndjacu, ndjacv, ii, jj, nd) as usize)
            .copied()
            .ok_or_else(er2_oob)
    };
    loop {
        if nv > mindgv {
            let bid0 = *xmaxjv.get((nv - i2rdv) as usize).ok_or_else(er2_oob)?;
            for nd in 1..=ndimen {
                let mut bid1 = 0.0;
                if i2rdu <= nu {
                    for ii in i2rdu..=nu {
                        let scale = *xmaxju.get((ii - i2rdu) as usize).ok_or_else(er2_oob)?;
                        bid1 += coeff(ii, nv, nd)?.abs() * scale * bid0;
                    }
                }
                let slot = work.get_mut((nd - 1) as usize).ok_or_else(er2_oob)?;
                *slot = bid1;
            }
        } else if ndimen >= 1 {
            work[0] = epmscut * 2.0;
        }
        let errnv_inc = super::border_error::mzsnorm(&work);
        if nu > mindgu {
            let bid0 = *xmaxju.get((nu - i2rdu) as usize).ok_or_else(er2_oob)?;
            for nd in 1..=ndimen {
                let mut bid1 = 0.0;
                if i2rdv <= nv {
                    for jj in i2rdv..=nv {
                        let scale = *xmaxjv.get((jj - i2rdv) as usize).ok_or_else(er2_oob)?;
                        bid1 += coeff(nu, jj, nd)?.abs() * scale * bid0;
                    }
                }
                let slot = work.get_mut((nd - 1) as usize).ok_or_else(er2_oob)?;
                *slot = bid1;
            }
        } else if ndimen >= 1 {
            work[0] = epmscut * 2.0;
        }
        let errnu_inc = super::border_error::mzsnorm(&work);
        let errnu = super::border_error::mzsnorm(&[erreur, errnu_inc]);
        let errnv = super::border_error::mzsnorm(&[erreur, errnv_inc]);
        if errnu > errnv {
            if errnv < epmscut {
                erreur = errnv;
                nv -= 1;
            } else {
                break;
            }
        } else if errnu < epmscut {
            erreur = errnu;
            nu -= 1;
        } else {
            break;
        }
    }
    Ok((nu.max(1), nv.max(1), erreur))
}

/// `mma2moy_` (`AdvApp2Var_ApproxF2var.cxx:7767-7921`).
fn mma2moy(
    ndjacu: i32,
    ndjacv: i32,
    ndimen: i32,
    mindgu: i32,
    maxdgu: i32,
    mindgv: i32,
    maxdgv: i32,
    iordru: i32,
    iordrv: i32,
    patjac: &[f64],
) -> Result<f64, String> {
    let idebu = (iordru + 1) << 1;
    let idebv = (iordrv + 1) << 1;
    let minu = idebu.max(mindgu);
    let minv = idebv.max(mindgv);
    let mut bid0 = 0.0;
    let coeff = |ii: i32, jj: i32, nd: i32| -> Result<f64, String> {
        patjac
            .get(pat_at(ndjacu, ndjacv, ii, jj, nd) as usize)
            .copied()
            .ok_or_else(moy_oob)
    };
    if minv <= maxdgv && idebu <= maxdgu {
        for nd in 1..=ndimen {
            for jj in minv..=maxdgv {
                for ii in idebu..=maxdgu {
                    let bid1 = coeff(ii, jj, nd)?;
                    bid0 += bid1 * bid1;
                }
            }
        }
    }
    if idebv <= minv - 1 && minu <= maxdgu {
        for nd in 1..=ndimen {
            for jj in idebv..=(minv - 1) {
                for ii in minu..=maxdgu {
                    let bid1 = coeff(ii, jj, nd)?;
                    bid0 += bid1 * bid1;
                }
            }
        }
    }
    Ok((bid0 / 4.0).sqrt())
}

fn zero_rejected(
    patjac: &mut [f64],
    nu1: i32,
    nv1: i32,
    ndjacu: i32,
    ndjacv: i32,
    ndimen: i32,
) -> Result<(), String> {
    for nd in 1..=ndimen {
        if nv1 > ndjacv {
            continue;
        }
        for jj in nv1..=ndjacv {
            if nu1 > ndjacu {
                continue;
            }
            for kk in nu1..=ndjacu {
                let at = pat_at(ndjacu, ndjacv, kk, jj, nd) as usize;
                let slot = patjac.get_mut(at).ok_or_else(|| {
                    "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4704".to_string()
                })?;
                *slot = 0.0;
            }
        }
    }
    Ok(())
}

/// One finished square: canonical polynomial, reduced coefficient counts, and
/// the errors stored by `AddErrors` (`AdvApp2Var_Patch.cxx:886-957`).
pub(super) struct FinishedPatch {
    pub equation: Vec<f64>,
    pub nb_coeff_u: i32,
    pub nb_coeff_v: i32,
    pub max_error: f64,
    pub moy_error: f64,
    pub iso_errors: [f64; 4],
}

/// Degree reduction (`AdvApp2Var_ApproxF2var.cxx:4639-4714` and `:4828-4915`),
/// then `mma2can_` (`AdvApp2Var_Patch.cxx:902`), `AddConstraints` (`:920`),
/// `AddErrors` (`:921`), `mma2fx6_` (`ApproxF2var.cxx:7285-7572`), and the
/// Fortran-to-C++ transpose (`Patch.cxx:942-957`). The loop then runs
/// `Compute3DErrors`, and `ConvertBS` returns the B-spline.
pub(super) fn canon_stop(
    ndjacu: i32,
    ndjacv: i32,
    ndimen: i32,
    iordru: i32,
    iordrv: i32,
    ndmin_u: i32,
    ndmin_v: i32,
    ndguli: i32,
    ndgvli: i32,
    reduce: bool,
    erreur: f64,
    nu_keep: i32,
    nv_keep: i32,
    eps: f64,
    xmaxju: &[f64],
    xmaxjv: &[f64],
    patjac: &mut [f64],
    frame: &super::framework::Framework,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    nbsesp: i32,
    epsapr: &[f64],
    epsfro: &[f64],
) -> Result<FinishedPatch, String> {
    let i2rdu = (iordru + 1) << 1;
    let i2rdv = (iordrv + 1) << 1;
    let (nu, nv, kept) = if reduce && ndguli >= i2rdu && ndgvli >= i2rdv {
        let floor_u = 1.max((iordru << 1) + 1).max(ndmin_u);
        let floor_v = 1.max((iordrv << 1) + 1).max(ndmin_v);
        mma2er2(
            ndjacu, ndjacv, ndimen, floor_u, ndguli, floor_v, ndgvli, iordru, iordrv, xmaxju,
            xmaxjv, patjac, eps, erreur,
        )?
    } else {
        (nu_keep, nv_keep, erreur)
    };
    let errmoy = mma2moy(
        ndjacu,
        ndjacv,
        ndimen,
        nu + 1,
        ndjacu,
        nv + 1,
        ndjacv,
        iordru,
        iordrv,
        patjac,
    )?;
    zero_rejected(patjac, nu + 1, nv + 1, ndjacu, ndjacv, ndimen)?;
    mmfmca9(ndjacu + 1, ndjacv + 1, nu + 1, nv + 1, ndimen, patjac)?;
    let mut patcan = mma2can(
        ndguli + 1,
        ndgvli + 1,
        ndimen,
        iordru,
        iordrv,
        nu + 1,
        nv + 1,
        patjac,
    )?;
    super::constraints::add_constraints(
        frame,
        u0,
        u1,
        v0,
        v1,
        ndimen,
        iordru,
        iordrv,
        ndguli + 1,
        ndgvli + 1,
        &mut patcan,
    )?;
    let (errmax, errmoy, iso_err) = super::constraints::add_errors(
        frame, u0, u1, v0, v1, iordru, iordrv, kept, errmoy,
    )?;
    let ncfmxu = ndguli + 1;
    let ncfmxv = ndgvli + 1;
    let mut ncoef_u = [nu + 1];
    let mut ncoef_v = [nv + 1];
    // `MakeApprox` passes one subspace of dimension 3 (`itrois`) and one patch.
    super::degree_reduce::mma2fx6(
        ncfmxu,
        ncfmxv,
        ndimen,
        nbsesp,
        &[3],
        1,
        1,
        iordru,
        iordrv,
        epsapr,
        epsfro,
        &patcan,
        &[errmax],
        &mut ncoef_u,
        &mut ncoef_v,
    )?;
    let mut equation = vec![0.0; super::degree_reduce::hpaux_len(ndjacu, ndjacv, ndimen)?];
    super::degree_reduce::transpose_can(ndimen, ncfmxu, ncfmxv, &patcan, &mut equation)?;
    Ok(FinishedPatch {
        equation,
        nb_coeff_u: *ncoef_u.first().ok_or_else(|| super::degree_reduce::fx6(7556))?,
        nb_coeff_v: *ncoef_v.first().ok_or_else(|| super::degree_reduce::fx6(7557))?,
        max_error: errmax,
        moy_error: errmoy,
        iso_errors: iso_err,
    })
}

fn mmfmca9_oob() -> String {
    "AdvApp2Var_MathBase::mmfmca9_: AdvApp2Var_MathBase.cxx:3966".into()
}

fn packed(i: i32, j: i32, k: i32, dim1: i32, dim2: i32) -> i32 {
    (i - 1) + (j - 1) * dim1 + (k - 1) * dim1 * dim2
}

/// `mmfmca9_` (`AdvApp2Var_MathBase.cxx:3873-4005`). In-place on `PATJAC`.
fn mmfmca9(
    ndimax: i32,
    ncfumx: i32,
    ndimen: i32,
    ncoefu: i32,
    ncoefv: i32,
    tab: &mut [f64],
) -> Result<(), String> {
    if ndimen != ndimax {
        if ncoefv >= 1 && ncoefu >= 1 && ndimen >= 1 {
            for k in 1..=ncoefv {
                for j in 1..=ncoefu {
                    for i in 1..=ndimen {
                        let src = packed(i, j, k, ndimax, ncfumx) as usize;
                        let dst = packed(i, j, k, ndimen, ncoefu) as usize;
                        let value = *tab.get(src).ok_or_else(mmfmca9_oob)?;
                        let slot = tab.get_mut(dst).ok_or_else(mmfmca9_oob)?;
                        *slot = value;
                    }
                }
            }
        }
        return Ok(());
    }
    if ncoefu != ncfumx {
        let count = (ndimen * ncoefu).max(0) as usize;
        if ncoefv >= 1 {
            for k in 1..=ncoefv {
                let src = ((k - 1) * ndimax * ncfumx).max(0) as usize;
                let dst = ((k - 1) * ndimen * ncoefu).max(0) as usize;
                copy_fill(tab, src, dst, count, &mmfmca9_oob())?;
            }
        }
        return Ok(());
    }
    let count = (ndimen * ncoefu * ncoefv).max(0) as usize;
    copy_fill(tab, 0, 0, count, &mmfmca9_oob())
}

/// `mcrfill_` (`AdvApp2Var_SysBase.cxx:2268-2295`), counted in doubles.
fn copy_fill(
    tab: &mut [f64],
    src: usize,
    dst: usize,
    count: usize,
    oob: &str,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    let src_end = src.checked_add(count).ok_or_else(|| oob.to_string())?;
    let dst_end = dst.checked_add(count).ok_or_else(|| oob.to_string())?;
    if src_end > tab.len() || dst_end > tab.len() {
        return Err(oob.to_string());
    }
    if src == dst {
        return Ok(());
    }
    if src.abs_diff(dst) >= count {
        let copied = tab[src..src + count].to_vec();
        tab[dst..dst + count].copy_from_slice(&copied);
    } else if src > dst {
        for n in 0..count {
            tab[dst + n] = tab[src + n];
        }
    } else {
        for n in (0..count).rev() {
            tab[dst + n] = tab[src + n];
        }
    }
    Ok(())
}

fn mma2can_oob() -> String {
    "AdvApp2Var_Patch::MakeApprox : Error in FORTRAN".into()
}

fn can_tail_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2can_: AdvApp2Var_ApproxF2var.cxx:2332".into()
}

fn fmca8_oob() -> String {
    "AdvApp2Var_MathBase::mmfmca8_: AdvApp2Var_MathBase.cxx:3830".into()
}

/// `mma2can_` (`AdvApp2Var_ApproxF2var.cxx:2205-2361`), called from
/// `AdvApp2Var_Patch.cxx:902`. `NCFLMU` is `ndguli + 1`.
fn mma2can(
    ncfmxu: i32,
    ncfmxv: i32,
    ndimen: i32,
    iordru: i32,
    iordrv: i32,
    ncoefu: i32,
    ncoefv: i32,
    patjac: &[f64],
) -> Result<Vec<f64>, String> {
    if !(-1..=2).contains(&iordru)
        || !(-1..=2).contains(&iordrv)
        || ncoefu > ncfmxu
        || ncoefv > ncfmxv
    {
        return Err(mma2can_oob());
    }
    let can_len = usize::try_from(ncfmxu)
        .ok()
        .and_then(|u| usize::try_from(ncfmxv).ok().and_then(|v| u.checked_mul(v)))
        .and_then(|n| usize::try_from(ndimen).ok().and_then(|d| n.checked_mul(d)))
        .ok_or_else(can_tail_oob)?;
    let mut patcan = vec![0.0; can_len];
    super::jacobi_canon::mmjacpt(ndimen, ncoefu, ncoefv, iordru, iordrv, patjac, &mut patcan)?;
    // `mma2can_` passes `(ncoefu, ncoefv, ndimen, ncfmxu, ncfmxv, ndimen)`.
    mmfmca8(ncoefu, ncoefv, ndimen, ncfmxu, ncfmxv, &mut patcan)?;
    let ilon1 = ncfmxu - ncoefu;
    let ilon2 = ncfmxu
        .checked_mul(ncfmxv - ncoefv)
        .ok_or_else(can_tail_oob)?;
    if ndimen >= 1 {
        for nd in 1..=ndimen {
            if ilon1 > 0 && ncoefv >= 1 {
                for ii in 1..=ncoefv {
                    let start = tail_at(ncoefu, ii - 1, nd - 1, ncfmxu, ncfmxv)?;
                    raz(&mut patcan, start, ilon1)?;
                }
            }
            if ilon2 > 0 {
                let start = tail_at(0, ncoefv, nd - 1, ncfmxu, ncfmxv)?;
                raz(&mut patcan, start, ilon2)?;
            }
        }
    }
    Ok(patcan)
}

/// `mmfmca8_` (`AdvApp2Var_MathBase.cxx:3745-3869`). In-place expand.
fn mmfmca8(
    ndimen: i32,
    ncoefu: i32,
    ncoefv: i32,
    ndimax: i32,
    ncfumx: i32,
    tab: &mut [f64],
) -> Result<(), String> {
    if ndimax == ndimen {
        if ncoefu == ncfumx {
            let count = (ndimen * ncoefu * ncoefv).max(0) as usize;
            return copy_fill(tab, 0, 0, count, &fmca8_oob());
        }
        let count = (ndimen * ncoefu).max(0) as usize;
        if ncoefv >= 1 {
            for k in (1..=ncoefv).rev() {
                let src = ((k - 1) * ndimen * ncoefu).max(0) as usize;
                let dst = ((k - 1) * ndimax * ncfumx).max(0) as usize;
                copy_fill(tab, src, dst, count, &fmca8_oob())?;
            }
        }
        return Ok(());
    }
    if ncoefv >= 1 && ncoefu >= 1 && ndimen >= 1 {
        for k in (1..=ncoefv).rev() {
            for j in (1..=ncoefu).rev() {
                for i in (1..=ndimen).rev() {
                    let src = packed(i, j, k, ndimen, ncoefu) as usize;
                    let dst = packed(i, j, k, ndimax, ncfumx) as usize;
                    let value = *tab.get(src).ok_or_else(fmca8_oob)?;
                    let slot = tab.get_mut(dst).ok_or_else(fmca8_oob)?;
                    *slot = value;
                }
            }
        }
    }
    Ok(())
}

fn tail_at(u: i32, v: i32, nd: i32, ncfmxu: i32, ncfmxv: i32) -> Result<i32, String> {
    let stride = ncfmxu.checked_mul(ncfmxv).ok_or_else(can_tail_oob)?;
    let row = v.checked_mul(ncfmxu).ok_or_else(can_tail_oob)?;
    let plane = nd.checked_mul(stride).ok_or_else(can_tail_oob)?;
    u.checked_add(row)
        .and_then(|n| n.checked_add(plane))
        .ok_or_else(can_tail_oob)
}

fn raz(tab: &mut [f64], start: i32, count: i32) -> Result<(), String> {
    if count <= 0 {
        return Ok(());
    }
    let start = usize::try_from(start).map_err(|_| can_tail_oob())?;
    let count = count as usize;
    let end = start.checked_add(count).ok_or_else(can_tail_oob)?;
    let slots = tab.get_mut(start..end).ok_or_else(can_tail_oob)?;
    for slot in slots {
        *slot = 0.0;
    }
    Ok(())
}
