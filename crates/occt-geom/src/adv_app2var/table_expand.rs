//! `mmfmca8_` (`AdvApp2Var_MathBase.cxx:3745-3869`).
//!
//! Expands a packed coefficient table into the caller's larger leading dimensions.

fn idx(i: i32, j: i32, k: i32, dim1: i32, dim2: i32) -> usize {
    ((i - 1) + (j - 1) * dim1 + (k - 1) * dim1 * dim2) as usize
}

pub(super) fn mmfmca8(
    ndimen: i32,
    ncoefu: i32,
    ncoefv: i32,
    ndimax: i32,
    ncfumx: i32,
    _ncfvmx: i32,
    tabini: &[f64],
    tabres: &mut [f64],
) -> Result<(), String> {
    if ndimax == ndimen {
        if ncoefu == ncfumx {
            let count = (ndimen * ncoefu * ncoefv) as usize;
            if count > tabini.len() || count > tabres.len() {
                return Err("AdvApp2Var_MathBase::mmfmca8_: AdvApp2Var_MathBase.cxx:3861".into());
            }
            tabres[..count].copy_from_slice(&tabini[..count]);
            return Ok(());
        }
        let count = (ndimen * ncoefu) as usize;
        for k in (1..=ncoefv).rev() {
            let src = idx(1, 1, k, ndimen, ncoefu);
            let dst = idx(1, 1, k, ndimax, ncfumx);
            if src + count > tabini.len() || dst + count > tabres.len() {
                return Err("AdvApp2Var_MathBase::mmfmca8_: AdvApp2Var_MathBase.cxx:3850".into());
            }
            tabres[dst..dst + count].copy_from_slice(&tabini[src..src + count]);
        }
        return Ok(());
    }
    for k in (1..=ncoefv).rev() {
        for j in (1..=ncoefu).rev() {
            for i in (1..=ndimen).rev() {
                let src = idx(i, j, k, ndimen, ncoefu);
                let dst = idx(i, j, k, ndimax, ncfumx);
                if src >= tabini.len() || dst >= tabres.len() {
                    return Err("AdvApp2Var_MathBase::mmfmca8_: AdvApp2Var_MathBase.cxx:3830".into());
                }
                tabres[dst] = tabini[src];
            }
        }
    }
    Ok(())
}
