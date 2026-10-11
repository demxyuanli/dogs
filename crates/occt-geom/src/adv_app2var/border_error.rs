//! `mma1fer_` (`AdvApp2Var_ApproxF2var.cxx:985-1200`).
//!
//! Continuity C0, C1, and C2 use Jacobi orders 2, 4, and 6.
//! Jacobi order -1 (`ia == 0`) uses `mmtrpj0_` and `mmaper0_`.

use super::jacobi_max::{XMAX2, XMAX6};
use super::xmaxj4::XMAXJ4;

fn jac_at(coeff: i32, dimen: i32, stride: i32) -> usize {
    ((coeff - 1) + (dimen - 1) * stride) as usize
}

/// `mzsnorm_` (`AdvApp2Var_MathBase.cxx:10503-10612`).
pub(super) fn mzsnorm(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut irmax = 0usize;
    for (index, value) in values.iter().enumerate().skip(1) {
        if values[irmax].abs() < value.abs() {
            irmax = index;
        }
    }
    if values[irmax].abs() < 1.0 {
        let sum: f64 = values.iter().map(|value| value * value).sum();
        sum.sqrt()
    } else {
        let pivot = values[irmax];
        let mut sum = 0.0;
        for (index, value) in values.iter().enumerate() {
            if index == irmax {
                sum += 1.0;
            } else {
                let ratio = value / pivot;
                sum += ratio * ratio;
            }
        }
        pivot.abs() * sum.sqrt()
    }
}

/// `mmeps1_` (`AdvApp2Var_MathBase.cxx:3385-3463`).
///
/// `mmprcsn_.eps1` is a zero-initialized static (`:221-225`).
fn eps1() -> f64 {
    0.0
}

/// `mmtrpj4_` (`AdvApp2Var_MathBase.cxx:8600-8768`).
fn mmtrpj4(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    epsi3d: f64,
    crvjac: &[f64],
    idim: i32,
) -> Result<(i32, f64), String> {
    let ia = 4;
    let mut ncfnew = ia;
    let mut ycvmax = vec![0.0; ndimen.max(0) as usize];
    let mut epstrc = 0.0;
    let ncut = ia + 1;
    for coeff in (ncut..=ncoeff).rev() {
        let slot = (coeff - ncut) as usize;
        if slot >= XMAXJ4.len() {
            return Err("AdvApp2Var_MathBase::mmtrpj4_: AdvApp2Var_MathBase.cxx:8716".into());
        }
        let bidon = XMAXJ4[slot];
        for nd in 1..=ndimen {
            let at = jac_at(coeff, idim + nd - 1, ncofmx);
            ycvmax[(nd - 1) as usize] += crvjac[at].abs() * bidon;
        }
        let error = mzsnorm(&ycvmax);
        if error > epsi3d {
            ncfnew = coeff;
            break;
        }
        epstrc = error;
    }
    if ncfnew == ia {
        let eps = eps1();
        let mut reduced = 1;
        for coeff in (2..=ia).rev() {
            let mut bid = 0.0;
            for nd in 1..=ndimen {
                bid += crvjac[jac_at(coeff, idim + nd - 1, ncofmx)].abs();
            }
            if bid > eps {
                reduced = coeff;
                break;
            }
        }
        ncfnew = reduced;
    }
    Ok((ncfnew, epstrc))
}

fn legendre_peak(index: i32) -> f64 {
    ((((index - 1) as f64) * 2.0 + 1.0) / 2.0).sqrt()
}

/// `mmtrpj0_` (`AdvApp2Var_MathBase.cxx:8314-8423`). Stops at degree 1.
fn mmtrpj0(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    epsi3d: f64,
    crvjac: &[f64],
    idim: i32,
) -> Result<(i32, f64), String> {
    let mut ncfnew = 1;
    let mut ycvmax = vec![0.0; ndimen.max(0) as usize];
    let mut epstrc = 0.0;
    let ncut = 2;
    for coeff in (ncut..=ncoeff).rev() {
        let bidon = legendre_peak(coeff);
        for nd in 1..=ndimen {
            let at = jac_at(coeff, idim + nd - 1, ncofmx);
            let value = *crvjac.get(at).ok_or_else(|| {
                "AdvApp2Var_MathBase::mmtrpj0_: AdvApp2Var_MathBase.cxx:8401".to_string()
            })?;
            ycvmax[(nd - 1) as usize] += value.abs() * bidon;
        }
        let error = mzsnorm(&ycvmax);
        if error > epsi3d {
            ncfnew = coeff;
            break;
        }
        epstrc = error;
    }
    Ok((ncfnew, epstrc))
}

/// `mmtrpj2_` (`AdvApp2Var_MathBase.cxx:8427-8596`).
fn mmtrpj2(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    epsi3d: f64,
    crvjac: &[f64],
    idim: i32,
) -> Result<(i32, f64), String> {
    trunc_jacobi(2, ncofmx, ndimen, ncoeff, epsi3d, crvjac, idim, &XMAX2, 8544, 8548, 8578)
}

/// `mmtrpj6_` (`AdvApp2Var_MathBase.cxx:8772-8939`).
fn mmtrpj6(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    epsi3d: f64,
    crvjac: &[f64],
    idim: i32,
) -> Result<(i32, f64), String> {
    trunc_jacobi(6, ncofmx, ndimen, ncoeff, epsi3d, crvjac, idim, &XMAX6, 8887, 8891, 8921)
}

fn trunc_jacobi(
    ia: i32,
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    epsi3d: f64,
    crvjac: &[f64],
    idim: i32,
    xmax: &[f64],
    slot_line: i32,
    coeff_line: i32,
    zero_line: i32,
) -> Result<(i32, f64), String> {
    let name = if ia == 2 { "mmtrpj2_" } else { "mmtrpj6_" };
    let oob = |line: i32| format!("AdvApp2Var_MathBase::{name}: AdvApp2Var_MathBase.cxx:{line}");
    let mut ncfnew = ia;
    let mut ycvmax = vec![0.0; ndimen.max(0) as usize];
    let mut epstrc = 0.0;
    let ncut = ia + 1;
    for coeff in (ncut..=ncoeff).rev() {
        let slot = (coeff - ncut) as usize;
        if slot >= xmax.len() {
            return Err(oob(slot_line));
        }
        let bidon = xmax[slot];
        for nd in 1..=ndimen {
            let at = jac_at(coeff, idim + nd - 1, ncofmx);
            let value = *crvjac.get(at).ok_or_else(|| oob(coeff_line))?;
            ycvmax[(nd - 1) as usize] += value.abs() * bidon;
        }
        let error = mzsnorm(&ycvmax);
        if error > epsi3d {
            ncfnew = coeff;
            break;
        }
        epstrc = error;
    }
    if ncfnew == ia {
        let eps = eps1();
        let mut reduced = 1;
        for coeff in (2..=ia).rev() {
            let mut bid = 0.0;
            for nd in 1..=ndimen {
                let at = jac_at(coeff, idim + nd - 1, ncofmx);
                bid += crvjac.get(at).ok_or_else(|| oob(zero_line))?.abs();
            }
            if bid > eps {
                reduced = coeff;
                break;
            }
        }
        ncfnew = reduced;
    }
    Ok((ncfnew, epstrc))
}

/// `mmaper0_` (`AdvApp2Var_MathBase.cxx:487-591`).
fn mmaper0(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    crvjac: &[f64],
    idim: i32,
    ncfnew: i32,
) -> Result<f64, String> {
    let mut ycvmax = vec![0.0; ndimen.max(0) as usize];
    let mut ncut = 1;
    if ncfnew >= ncut {
        ncut = ncfnew + 1;
    }
    for coeff in ncut..=ncoeff {
        let bidon = legendre_peak(coeff);
        for nd in 1..=ndimen {
            let at = jac_at(coeff, idim + nd - 1, ncofmx);
            let value = *crvjac.get(at).ok_or_else(|| {
                "AdvApp2Var_MathBase::mmaper0_: AdvApp2Var_MathBase.cxx:578".to_string()
            })?;
            ycvmax[(nd - 1) as usize] += value.abs() * bidon;
        }
    }
    Ok(mzsnorm(&ycvmax))
}

/// `mmaper2_` (`AdvApp2Var_MathBase.cxx:595-728`).
fn mmaper2(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    crvjac: &[f64],
    idim: i32,
    ncfnew: i32,
) -> Result<f64, String> {
    max_error_jacobi(3, ncofmx, ndimen, ncoeff, crvjac, idim, ncfnew, &XMAX2, "mmaper2_", 711, 715)
}

/// `mmaper6_` (`AdvApp2Var_MathBase.cxx:874-1005`).
fn mmaper6(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    crvjac: &[f64],
    idim: i32,
    ncfnew: i32,
) -> Result<f64, String> {
    max_error_jacobi(7, ncofmx, ndimen, ncoeff, crvjac, idim, ncfnew, &XMAX6, "mmaper6_", 988, 992)
}

fn max_error_jacobi(
    idec: i32,
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    crvjac: &[f64],
    idim: i32,
    ncfnew: i32,
    xmax: &[f64],
    name: &str,
    slot_line: i32,
    coeff_line: i32,
) -> Result<f64, String> {
    let mut ycvmax = vec![0.0; ndimen.max(0) as usize];
    let ncut = idec.max(ncfnew + 1);
    for coeff in ncut..=ncoeff {
        let slot = (coeff - idec) as usize;
        if slot >= xmax.len() {
            return Err(format!(
                "AdvApp2Var_MathBase::{name}: AdvApp2Var_MathBase.cxx:{slot_line}"
            ));
        }
        let bidon = xmax[slot];
        for nd in 1..=ndimen {
            let at = jac_at(coeff, idim + nd - 1, ncofmx);
            let value = *crvjac.get(at).ok_or_else(|| {
                format!("AdvApp2Var_MathBase::{name}: AdvApp2Var_MathBase.cxx:{coeff_line}")
            })?;
            ycvmax[(nd - 1) as usize] += value.abs() * bidon;
        }
    }
    Ok(mzsnorm(&ycvmax))
}

/// `mmaper4_` (`AdvApp2Var_MathBase.cxx:738-870`).
fn mmaper4(
    ncofmx: i32,
    ndimen: i32,
    ncoeff: i32,
    crvjac: &[f64],
    idim: i32,
    ncfnew: i32,
) -> Result<f64, String> {
    let mut ycvmax = vec![0.0; ndimen.max(0) as usize];
    let idec = 5;
    let ncut = idec.max(ncfnew + 1);
    for coeff in ncut..=ncoeff {
        let slot = (coeff - idec) as usize;
        if slot >= XMAXJ4.len() {
            return Err("AdvApp2Var_MathBase::mmaper4_: AdvApp2Var_MathBase.cxx:853".into());
        }
        let bidon = XMAXJ4[slot];
        for nd in 1..=ndimen {
            let at = jac_at(coeff, idim + nd - 1, ncofmx);
            ycvmax[(nd - 1) as usize] += crvjac[at].abs() * bidon;
        }
    }
    Ok(mzsnorm(&ycvmax))
}

/// `mmaperm_` (`AdvApp2Var_ApproxF2var.cxx:8198-8309`).
fn mmaperm(
    ncofmx: i32,
    ndim: i32,
    ncoeff: i32,
    iordre: i32,
    crvjac: &[f64],
    idim: i32,
    ncfnew: i32,
) -> f64 {
    let ia = (iordre + 1) << 1;
    let mut ncfcut = ia + 1;
    if ncfnew + 1 > ncfcut {
        ncfcut = ncfnew + 1;
    }
    let mut bid = 0.0;
    for nd in 1..=ndim {
        for coeff in ncfcut..=ncoeff {
            let value = crvjac[jac_at(coeff, idim + nd - 1, ncofmx)];
            bid += value * value;
        }
    }
    (bid / 2.0).sqrt()
}

/// Degree and errors of one border (`AdvApp2Var_ApproxF2var.cxx:1078-1181`).
///
/// `Ok(0)` met the tolerance. `Ok(-1)` kept `ncflim` coefficients.
/// `ncoeff` is the number of significant coefficients.
pub(super) fn mma1fer(
    nbsesp: i32,
    ndimse: &[i32],
    iordre: i32,
    ndgjac: i32,
    crvjac: &mut [f64],
    ncflim: i32,
    epsapr: &[f64],
) -> Result<(i32, i32, f64, f64), String> {
    let ia = (iordre + 1) << 1;
    let ncfja = ndgjac + 1;
    let mut iercod = 0;
    let mut idim = 1;
    let mut ncoeff = 0;
    let mut errmax = 0.0;
    let mut errmoy = 0.0;
    for ii in 1..=nbsesp {
        let ndses = *ndimse.get((ii - 1) as usize).ok_or_else(|| {
            "AdvApp2Var_ApproxF2var::mma1fer_: AdvApp2Var_ApproxF2var.cxx:1097".to_string()
        })?;
        let eps = *epsapr.get((ii - 1) as usize).ok_or_else(|| {
            "AdvApp2Var_ApproxF2var::mma1fer_: AdvApp2Var_ApproxF2var.cxx:1105".to_string()
        })?;
        let (ncfnw, emax) = if ia == 0 {
            // `mmtrpjj_` (`AdvApp2Var_MathBase.cxx:9004-9006`).
            mmtrpj0(ncfja, ndses, ncfja, eps, crvjac, idim)?
        } else if ia == 2 {
            mmtrpj2(ncfja, ndses, ncfja, eps, crvjac, idim)?
        } else if ia == 4 {
            mmtrpj4(ncfja, ndses, ncfja, eps, crvjac, idim)?
        } else {
            // `mmtrpjj_` else (`AdvApp2Var_MathBase.cxx:9018`).
            mmtrpj6(ncfja, ndses, ncfja, eps, crvjac, idim)?
        };
        errmax = emax;
        if ncfnw <= ncflim {
            errmoy = mmaperm(ncfja, ndses, ncfja, iordre, crvjac, idim, ncfnw);
            if ncfnw > ncoeff {
                ncoeff = ncfnw;
            }
            let nbr0 = ncflim - ncfnw;
            if nbr0 > 0 {
                for kk in 1..=ndses {
                    let start = jac_at(ncfnw + 1, idim + kk - 1, ncfja);
                    for offset in 0..nbr0 {
                        let at = start + offset as usize;
                        if at >= crvjac.len() {
                            return Err(
                                "AdvApp2Var_ApproxF2var::mma1fer_: AdvApp2Var_ApproxF2var.cxx:1133"
                                    .into(),
                            );
                        }
                        crvjac[at] = 0.0;
                    }
                }
            }
        } else {
            iercod = -1;
            errmax = if ia == 0 {
                // `mmaperx_` (`AdvApp2Var_MathBase.cxx:1074-1076`).
                mmaper0(ncfja, ndses, ncfja, crvjac, idim, ncflim)?
            } else if ia == 2 {
                mmaper2(ncfja, ndses, ncfja, crvjac, idim, ncflim)?
            } else if ia == 4 {
                mmaper4(ncfja, ndses, ncfja, crvjac, idim, ncflim)?
            } else if ia == 6 {
                mmaper6(ncfja, ndses, ncfja, crvjac, idim, ncflim)?
            } else {
                return Err("AdvApp2Var_MathBase::mmaperx_: AdvApp2Var_MathBase.cxx:1092".into());
            };
            ncoeff = ncflim;
            errmoy = mmaperm(ncfja, ndses, ncfja, iordre, crvjac, idim, ncflim);
        }
        idim += ndses;
    }
    Ok((iercod, ncoeff, errmax, errmoy))
}

/// `mmapcmp_` (`AdvApp2Var_MathBase.cxx:374-483`).
///
/// Even degrees land in the first half of each dimension, odd degrees in the second.
pub(super) fn mmapcmp(
    ndim: i32,
    ncofmx: i32,
    ncoeff: i32,
    crvold: &[f64],
) -> Result<Vec<f64>, String> {
    if ncoeff < 1 || ncofmx < 1 || ndim < 1 {
        return Err("AdvApp2Var_MathBase::mmapcmp_: AdvApp2Var_MathBase.cxx:445".into());
    }
    let half = (ncoeff - 1) / 2 + 1;
    let mut crvnew = vec![0.0; (ndim * 2 * half) as usize];
    let ndegre = ncoeff - 1;
    for nd in 1..=ndim {
        let mut ipair = 0;
        for idg in 0..=(ndegre / 2) {
            let src = (ipair + (nd - 1) * ncofmx) as usize;
            let dst = (idg + (nd - 1) * 2 * half) as usize;
            if src >= crvold.len() || dst >= crvnew.len() {
                return Err("AdvApp2Var_MathBase::mmapcmp_: AdvApp2Var_MathBase.cxx:453".into());
            }
            crvnew[dst] = crvold[src];
            ipair += 2;
        }
        if ndegre >= 1 {
            let mut impair = 1;
            for idg in 0..=((ndegre - 1) / 2) {
                let src = (impair + (nd - 1) * ncofmx) as usize;
                let dst = (idg + (nd - 1) * 2 * half + half) as usize;
                if src >= crvold.len() || dst >= crvnew.len() {
                    return Err("AdvApp2Var_MathBase::mmapcmp_: AdvApp2Var_MathBase.cxx:465".into());
                }
                crvnew[dst] = crvold[src];
                impair += 2;
            }
        }
    }
    Ok(crvnew)
}
