//! `mmjacan_` (`AdvApp2Var_MathBase.cxx:5660-5788`).
//!
//! `mnfndeb_` returns 0, so the debug messages are omitted.

include!("plgcan.rs");

/// Jacobi coefficients, already split into even then odd blocks, become
/// canonical coefficients on [-1, 1].
pub(super) fn mmjacan(
    ideriv: i32,
    ndeg: i32,
    poljac: &[f64],
    polcan: &mut [f64],
) -> Result<(), String> {
    let even_last = ndeg / 2;
    for i in 0..=even_last {
        let mut bid = 0.0;
        let iptt = i * 31 - (i + 1) * i / 2 + 1;
        for j in i..=even_last {
            let idx = iptt + j + ideriv * 992 + 991;
            let src = j as usize;
            if idx < 0 || idx as usize >= PLGCAN.len() || src >= poljac.len() {
                return Err(
                    "UNPORTED AdvApp2Var_MathBase::mmjacan_: AdvApp2Var_MathBase.cxx:5750".into(),
                );
            }
            bid += PLGCAN[idx as usize] * poljac[src];
        }
        let dst = (i * 2) as usize;
        if dst >= polcan.len() {
            return Err(
                "UNPORTED AdvApp2Var_MathBase::mmjacan_: AdvApp2Var_MathBase.cxx:5753".into(),
            );
        }
        polcan[dst] = bid;
    }
    if ndeg == 0 {
        return Ok(());
    }
    let dim1 = ndeg / 2 + 1;
    let odd_last = (ndeg - 1) / 2;
    for i in 0..=odd_last {
        let mut bid = 0.0;
        let iptt = i * 31 - (i + 1) * i / 2 + 1;
        for j in i..=odd_last {
            let idx = iptt + j + ((ideriv << 1) + 1) * 496 + 991;
            let src = (j + dim1) as usize;
            if idx < 0 || idx as usize >= PLGCAN.len() || src >= poljac.len() {
                return Err(
                    "UNPORTED AdvApp2Var_MathBase::mmjacan_: AdvApp2Var_MathBase.cxx:5772".into(),
                );
            }
            bid += PLGCAN[idx as usize] * poljac[src];
        }
        let dst = ((i << 1) + 1) as usize;
        if dst >= polcan.len() {
            return Err(
                "UNPORTED AdvApp2Var_MathBase::mmjacan_: AdvApp2Var_MathBase.cxx:5776".into(),
            );
        }
        polcan[dst] = bid;
    }
    Ok(())
}

fn jaccv_oob() -> String {
    "AdvApp2Var_MathBase::mmjaccv_: AdvApp2Var_MathBase.cxx:5872".into()
}

fn jacpt_oob() -> String {
    "AdvApp2Var_ApproxF2var::mmjacpt_: AdvApp2Var_ApproxF2var.cxx:8744".into()
}

/// `mmjaccv_` (`AdvApp2Var_MathBase.cxx:5792-5895`).
/// Even Jacobi coefficients, then odd ones, then `mmjacan_`.
pub(super) fn mmjaccv(
    ncoef: i32,
    ndim: i32,
    ider: i32,
    crvlgd: &[f64],
    crvcan: &mut [f64],
) -> Result<(), String> {
    if ncoef < 1 {
        return Err(jaccv_oob());
    }
    let ndeg = ncoef - 1;
    let half = ndeg / 2 + 1;
    let mut polaux = vec![0.0; ncoef as usize];
    if ndim < 1 {
        return Ok(());
    }
    for nd in 1..=ndim {
        let row = ((nd - 1) * ncoef) as usize;
        let mut ii = 0i32;
        for i in 0..=(ndeg / 2) {
            let src = row + ii as usize;
            polaux[i as usize] = *crvlgd.get(src).ok_or_else(jaccv_oob)?;
            ii += 2;
        }
        if ndeg >= 1 {
            ii = 1;
            for i in 0..=((ndeg - 1) / 2) {
                let src = row + ii as usize;
                polaux[(i + half) as usize] = *crvlgd.get(src).ok_or_else(jaccv_oob)?;
                ii += 2;
            }
        }
        let end = row + ncoef as usize;
        let row_can = crvcan.get_mut(row..end).ok_or_else(jaccv_oob)?;
        mmjacan(ider, ndeg, &polaux, row_can)?;
    }
    Ok(())
}

fn square_at(ii: i32, jj: i32, nd: i32, dim1: i32, dim2: i32) -> usize {
    ((ii - 1) + (jj - 1) * dim1 + (nd - 1) * dim1 * dim2) as usize
}

/// `mmjacpt_` (`AdvApp2Var_ApproxF2var.cxx:8645-8792`).
/// Canonical by U, swap, canonical by V, swap back.
pub(super) fn mmjacpt(
    ndimen: i32,
    ncoefu: i32,
    ncoefv: i32,
    iordru: i32,
    iordrv: i32,
    ptclgd: &[f64],
    ptccan: &mut [f64],
) -> Result<(), String> {
    if ncoefu < 1 || ncoefv < 1 || ndimen < 1 {
        return Err(jacpt_oob());
    }
    let block = (ncoefu as usize)
        .checked_mul(ncoefv as usize)
        .and_then(|n| n.checked_mul(ndimen as usize))
        .ok_or_else(jacpt_oob)?;
    let mut ptcaux = vec![0.0; block.checked_mul(2).ok_or_else(jacpt_oob)?];
    let kdim_u = ndimen.checked_mul(ncoefv).ok_or_else(jacpt_oob)?;
    mmjaccv(ncoefu, kdim_u, iordru, ptclgd, ptccan)?;
    for nd in 1..=ndimen {
        for jj in 1..=ncoefv {
            for ii in 1..=ncoefu {
                let src = square_at(ii, jj, nd, ncoefu, ncoefv);
                let dst = square_at(jj, ii, nd, ncoefv, ncoefu);
                let value = *ptccan.get(src).ok_or_else(jacpt_oob)?;
                let slot = ptcaux.get_mut(dst).ok_or_else(jacpt_oob)?;
                *slot = value;
            }
        }
    }
    let kdim_v = ndimen.checked_mul(ncoefu).ok_or_else(jacpt_oob)?;
    let (head, tail) = ptcaux.split_at_mut(block);
    mmjaccv(ncoefv, kdim_v, iordrv, head, tail)?;
    for nd in 1..=ndimen {
        for jj in 1..=ncoefv {
            for ii in 1..=ncoefu {
                let src = block + square_at(jj, ii, nd, ncoefv, ncoefu);
                let dst = square_at(ii, jj, nd, ncoefu, ncoefv);
                let value = *ptcaux.get(src).ok_or_else(jacpt_oob)?;
                let slot = ptccan.get_mut(dst).ok_or_else(jacpt_oob)?;
                *slot = value;
            }
        }
    }
    Ok(())
}
