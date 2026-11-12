//! `mma1cdi_` (`cxx:233-465`), `mma1cnt_`, `mma1jak_`, `mmmapcoe_`, `mmmpocur_`.

use super::*;

/// `GeomConvert_ApproxSurface_Eval::Evaluate`
/// (`GeomConvert_ApproxSurface.cxx:55-310`).
///
/// Orders 0, 1, and 2 use `D0` / `D1` / `D2`. `D3` and `DN` are not on the
/// C1 iso and stay unported.
/// `mma1cdi_` (`AdvApp2Var_ApproxF2var.cxx:233-465`).
pub(crate) fn mma1cdi(
    ndimen: i32,
    nbroot: i32,
    rootlg: &[f64],
    iordre: i32,
    contr1: &[f64],
    contr2: &[f64],
    somtab: &mut [f64],
    diftab: &mut [f64],
) -> Result<(), String> {
    let hermit = mma1her(iordre)?;
    let ncfhe = (iordre + 1) << 1;
    let mut fpntab = vec![0.0; (ncfhe * nbroot) as usize];
    for ii in 1..=ncfhe {
        let poly = &hermit[((ii - 1) * ncfhe) as usize..(ii * ncfhe) as usize];
        for kk in 1..=nbroot {
            let value = mmmpocur(ncfhe, poly, rootlg[(kk - 1) as usize]);
            fpntab[((ii - 1) * nbroot + (kk - 1)) as usize] = value;
        }
    }

    let nroo2 = nbroot / 2;
    let dim1 = nroo2 + 1;
    for nd in 1..=ndimen {
        for ii in 1..=(iordre + 1) {
            let bid1 = contr1[contr_at(ii - 1, nd, ndimen)];
            let bid2 = contr2[contr_at(ii - 1, nd, ndimen)];
            let poly_m = (ii << 1) - 1;
            let poly_p = ii << 1;
            for kk in 1..=nroo2 {
                let kkm = nroo2 - kk + 1;
                let bid3 = bid1 * hermite_sample(&fpntab, nbroot, kkm, poly_m)
                    + bid2 * hermite_sample(&fpntab, nbroot, kkm, poly_p);
                somtab[pair_at(kk, nd, dim1)] -= bid3;
                diftab[pair_at(kk, nd, dim1)] += bid3;
            }
            for kk in 1..=nroo2 {
                let kkp = (nbroot + 1) / 2 + kk;
                let bid3 = bid1 * hermite_sample(&fpntab, nbroot, kkp, poly_m)
                    + bid2 * hermite_sample(&fpntab, nbroot, kkp, poly_p);
                somtab[pair_at(kk, nd, dim1)] -= bid3;
                diftab[pair_at(kk, nd, dim1)] -= bid3;
            }
        }
    }

    // The odd-degree arm applies only the last `ii` (`:433-445`).
    if nbroot % 2 == 1 {
        for nd in 1..=ndimen {
            let mut bid3 = 0.0;
            for ii in 1..=(iordre + 1) {
                let poly_m = (ii << 1) - 1;
                let poly_p = ii << 1;
                bid3 = hermite_sample(&fpntab, nbroot, nroo2 + 1, poly_m)
                    * contr1[contr_at(ii - 1, nd, ndimen)]
                    + hermite_sample(&fpntab, nbroot, nroo2 + 1, poly_p)
                        * contr2[contr_at(ii - 1, nd, ndimen)];
            }
            somtab[pair_at(0, nd, dim1)] -= bid3;
            diftab[pair_at(0, nd, dim1)] -= bid3;
        }
    }
    Ok(())
}

pub(crate) fn hermite_sample(fpntab: &[f64], nbroot: i32, root: i32, poly: i32) -> f64 {
    fpntab[((poly - 1) * nbroot + (root - 1)) as usize]
}

/// `mma1cnt_` (`AdvApp2Var_ApproxF2var.cxx:469-592`).
pub(crate) fn mma1cnt(
    ndimen: i32,
    iordre: i32,
    contr1: &[f64],
    contr2: &[f64],
    hermit: &[f64],
    ndgjac: i32,
    crvjac: &mut [f64],
) {
    let hdim = (iordre << 1) + 2;
    let jac_dim = ndgjac + 1;
    for nd in 1..=ndimen {
        for ii in 0..=((iordre << 1) + 1) {
            let mut bid = 0.0;
            for jj in 1..=(iordre + 1) {
                let at_m = (((jj << 1) - 2) * hdim + ii) as usize;
                let at_p = (((jj << 1) - 1) * hdim + ii) as usize;
                bid += contr1[contr_at(jj - 1, nd, ndimen)] * hermit[at_m]
                    + contr2[contr_at(jj - 1, nd, ndimen)] * hermit[at_p];
            }
            crvjac[(ii + (nd - 1) * jac_dim) as usize] = bid;
        }
    }
}

/// `mma1jak_` (`AdvApp2Var_ApproxF2var.cxx:1368-1486`).
pub(crate) fn mma1jak(
    ndimen: i32,
    nbroot: i32,
    iordre: i32,
    ndgjac: i32,
    somtab: &[f64],
    diftab: &[f64],
) -> Result<Vec<f64>, String> {
    let columns = ndgjac - 2 * iordre - 1;
    let cgauss = super::super::gauss::gauss_weights(ndgjac, nbroot, iordre, columns)?;
    Ok(mmmapcoe(
        ndimen, ndgjac, iordre, nbroot, somtab, diftab, &cgauss,
    ))
}

/// `mmmapcoe_` (`AdvApp2Var_ApproxF2var.cxx:8028-8194`).
pub(crate) fn mmmapcoe(
    ndim: i32,
    ndgjac: i32,
    iordre: i32,
    nbpnts: i32,
    somtab: &[f64],
    diftab: &[f64],
    gsstab: &[f64],
) -> Vec<f64> {
    let dim1 = nbpnts / 2 + 1;
    let nbroot = nbpnts / 2;
    let ikdeb = (iordre + 1) << 1;
    let jac_dim = ndgjac + 1;
    let mut crvjac = vec![0.0; (jac_dim * ndim).max(0) as usize];
    for nd in 1..=ndim {
        let mut ik = ikdeb;
        while ik <= ndgjac {
            let igss = ik - ikdeb;
            let mut bidon = 0.0;
            for ir in 1..=nbroot {
                bidon += somtab[pair_at(ir, nd, dim1)] * gsstab[(ir + igss * dim1) as usize];
            }
            crvjac[(ik + (nd - 1) * jac_dim) as usize] = bidon;
            ik += 2;
        }
        ik = ikdeb + 1;
        while ik <= ndgjac {
            let igss = ik - ikdeb;
            let mut bidon = 0.0;
            for ir in 1..=nbroot {
                bidon += diftab[pair_at(ir, nd, dim1)] * gsstab[(ir + igss * dim1) as usize];
            }
            crvjac[(ik + (nd - 1) * jac_dim) as usize] = bidon;
            ik += 2;
        }
    }
    if nbpnts % 2 != 0 {
        for nd in 1..=ndim {
            let mut ik = ikdeb;
            while ik <= ndgjac {
                let igss = ik - ikdeb;
                crvjac[(ik + (nd - 1) * jac_dim) as usize] +=
                    somtab[pair_at(0, nd, dim1)] * gsstab[(igss * dim1) as usize];
                ik += 2;
            }
        }
    }
    crvjac
}

/// `mmmpocur_` (`AdvApp2Var_MathBase.cxx:6596-6677`) for `NDIM = 1`.
pub(crate) fn mmmpocur(ndeg: i32, courbe: &[f64], tparam: f64) -> f64 {
    if ndeg < 1 {
        return 0.0;
    }
    let mut fu = courbe[(ndeg - 1) as usize];
    for i in (1..ndeg).rev() {
        fu = fu * tparam + courbe[(i - 1) as usize];
    }
    fu
}
