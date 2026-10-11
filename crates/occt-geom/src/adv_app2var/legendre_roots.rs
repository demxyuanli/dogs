//! `AdvApp2Var_MathBase::mmrtptt_` (`AdvApp2Var_MathBase.cxx:7921-8027`) and
//! `AdvApp2Var_ApproxF2var::mma2roo_` (`AdvApp2Var_ApproxF2var.cxx:7925-8024`).
//!
//! `mnfndeb_` returns 0, so the debug messages are not reached. `mcrfill_`
//! copies `nsur2` doubles (`AdvApp2Var_SysBase.cxx:2268-2295`).

include!("legendre_rootab.rs");

/// Positive Legendre roots of degree `ndglgd` (`mmrtptt_`).
fn positive_roots(ndglgd: i32) -> Result<Vec<f64>, String> {
    if ndglgd < 2 {
        return Ok(Vec::new());
    }
    let nsur2 = ndglgd / 2;
    let nmod2 = ndglgd % 2;
    let ideb = nsur2 * (nsur2 - 1) / 2 + 1;
    let start = (ideb + nmod2 * 465 - 1) as usize;
    let count = nsur2 as usize;
    let end = start + count;
    if end > ROOTAB.len() {
        return Err("AdvApp2Var_MathBase::mmrtptt_: AdvApp2Var_MathBase.cxx:8015".into());
    }
    Ok(ROOTAB[start..end].to_vec())
}

/// Roots on `(-1, 1)` for `nbpnt` interior samples (`mma2roo_`).
pub(super) fn legendre_roots(nbpnt: i32) -> Result<Vec<f64>, String> {
    let n = nbpnt.max(0) as usize;
    let mut roots = vec![0.0; n];
    let positive = positive_roots(nbpnt)?;
    let dest = ((nbpnt + 1) / 2).max(0) as usize;
    let end = dest + positive.len();
    if end > n {
        return Err("AdvApp2Var_ApproxF2var::mma2roo_: AdvApp2Var_ApproxF2var.cxx:7989".into());
    }
    roots[dest..end].copy_from_slice(&positive);
    let half = n / 2;
    for ii in 1..=half {
        roots[ii - 1] = -roots[n - ii];
    }
    if nbpnt % 2 == 1 && n > 0 {
        roots[n / 2] = 0.0;
    }
    Ok(roots)
}
