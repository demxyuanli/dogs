//! `CSLib::DNNUV` and `CSLib::DNNormal`.
//! Source: `CSLib.cxx:391-568`.

use crate::bspl::plib;
use crate::gp::GpVec;

fn at(grid: &[Vec<GpVec>], i: i32, j: i32) -> GpVec {
    grid[i as usize][j as usize]
}

fn at_f64(grid: &[Vec<f64>], i: i32, j: i32) -> f64 {
    grid[i as usize][j as usize]
}

fn set_vec(grid: &mut [Vec<GpVec>], i: i32, j: i32, v: GpVec) {
    grid[i as usize][j as usize] = v;
}

fn set_f64(grid: &mut [Vec<f64>], i: i32, j: i32, v: f64) {
    grid[i as usize][j as usize] = v;
}

/// `CSLib::DNNUV(Nu, Nv, DerSurf)` (`CSLib.cxx:391-408`).
/// `der_surf[i][j]` is `theDerSurf(i, j)`.
pub fn dnnuv(nu: i32, nv: i32, der_surf: &[Vec<GpVec>]) -> GpVec {
    let mut result = GpVec::new(0.0, 0.0, 0.0);
    for i in 0..=nu {
        for j in 0..=nv {
            let vg = at(der_surf, i + 1, j);
            let vd = at(der_surf, nu - i, nv + 1 - j);
            let cross = vg.crossed(&vd);
            let bin = plib::bin(nu, i) * plib::bin(nv, j);
            result = result.added(&cross.multiplied_scalar(bin));
        }
    }
    result
}

/// `CSLib::DNNUV(Nu, Nv, DerSurf1, DerSurf2)` (`CSLib.cxx:412-432`).
pub fn dnnuv2(nu: i32, nv: i32, der_surf1: &[Vec<GpVec>], der_surf2: &[Vec<GpVec>]) -> GpVec {
    let mut result = GpVec::new(0.0, 0.0, 0.0);
    for i in 0..=nu {
        for j in 0..=nv {
            let vg = at(der_surf1, i + 1, j);
            let vd = at(der_surf2, nu - i, nv + 1 - j);
            let cross = vg.crossed(&vd);
            let bin = plib::bin(nu, i) * plib::bin(nv, j);
            result = result.added(&cross.multiplied_scalar(bin));
        }
    }
    result
}

/// `CSLib::DNNormal` (`CSLib.cxx:436-568`).
/// `der_nuv[i][j]` is `theDerNUV(i, j)`. Defaults `iduref=0`, `idvref=0`.
pub fn dn_normal(
    nu: i32,
    nv: i32,
    der_nuv: &[Vec<GpVec>],
    iduref: i32,
    idvref: i32,
) -> GpVec {
    let kderiv = nu + nv;
    let n = (kderiv as usize) + 1;
    let mut der_vec_nor = vec![vec![GpVec::new(0.0, 0.0, 0.0); n]; n];
    let mut tab_scal = vec![vec![0.0; n]; n];
    let mut tab_norm = vec![vec![0.0; n]; n];

    let der_nor0 = at(der_nuv, iduref, idvref).normalized();
    set_vec(&mut der_vec_nor, 0, 0, der_nor0);
    let dnorm0 = at(der_nuv, iduref, idvref).dot(&at(&der_vec_nor, 0, 0));
    set_f64(&mut tab_norm, 0, 0, dnorm0);
    set_f64(&mut tab_scal, 0, 0, 0.0);

    for mderiv in 1..=kderiv {
        for pderiv in 0..=mderiv {
            let qderiv = mderiv - pderiv;
            if pderiv > nu || qderiv > nv {
                continue;
            }
            let mut scal = 0.0;
            if pderiv > qderiv {
                for jderiv in 1..=qderiv {
                    scal -= plib::bin(qderiv, jderiv)
                        * at(&der_vec_nor, 0, jderiv)
                            .dot(&at(&der_vec_nor, pderiv, qderiv - jderiv));
                }
                for jderiv in 0..qderiv {
                    scal -= plib::bin(qderiv, jderiv)
                        * at(&der_vec_nor, pderiv, jderiv)
                            .dot(&at(&der_vec_nor, 0, qderiv - jderiv));
                }
                for ideriv in 1..pderiv {
                    for jderiv in 0..=qderiv {
                        scal -= plib::bin(pderiv, ideriv)
                            * plib::bin(qderiv, jderiv)
                            * at(&der_vec_nor, ideriv, jderiv)
                                .dot(&at(&der_vec_nor, pderiv - ideriv, qderiv - jderiv));
                    }
                }
            } else {
                for ideriv in 1..=pderiv {
                    scal -= plib::bin(pderiv, ideriv)
                        * at(&der_vec_nor, ideriv, 0)
                            .dot(&at(&der_vec_nor, pderiv - ideriv, qderiv));
                }
                for ideriv in 0..pderiv {
                    scal -= plib::bin(pderiv, ideriv)
                        * at(&der_vec_nor, ideriv, qderiv)
                            .dot(&at(&der_vec_nor, pderiv - ideriv, 0));
                }
                for ideriv in 0..=pderiv {
                    for jderiv in 1..qderiv {
                        scal -= plib::bin(pderiv, ideriv)
                            * plib::bin(qderiv, jderiv)
                            * at(&der_vec_nor, ideriv, jderiv)
                                .dot(&at(&der_vec_nor, pderiv - ideriv, qderiv - jderiv));
                    }
                }
            }
            set_f64(&mut tab_scal, pderiv, qderiv, scal / 2.0);

            let mut dnorm = at(der_nuv, pderiv + iduref, qderiv + idvref)
                .dot(&at(&der_vec_nor, 0, 0));
            for jderiv in 0..qderiv {
                dnorm -= plib::bin(qderiv + idvref, jderiv + idvref)
                    * at_f64(&tab_norm, pderiv, jderiv)
                    * at_f64(&tab_scal, 0, qderiv - jderiv);
            }
            for ideriv in 0..pderiv {
                for jderiv in 0..=qderiv {
                    dnorm -= plib::bin(pderiv + iduref, ideriv + iduref)
                        * plib::bin(qderiv + idvref, jderiv + idvref)
                        * at_f64(&tab_norm, ideriv, jderiv)
                        * at_f64(&tab_scal, pderiv - ideriv, qderiv - jderiv);
                }
            }
            set_f64(&mut tab_norm, pderiv, qderiv, dnorm);

            let mut der_nor = at(der_nuv, pderiv + iduref, qderiv + idvref);
            for jderiv in 1..=qderiv {
                der_nor = der_nor.subtracted(
                    &at(&der_vec_nor, pderiv, qderiv - jderiv).multiplied_scalar(
                        plib::bin(pderiv + iduref, iduref)
                            * plib::bin(qderiv + idvref, jderiv + idvref)
                            * at_f64(&tab_norm, 0, jderiv),
                    ),
                );
            }
            for ideriv in 1..=pderiv {
                for jderiv in 0..=qderiv {
                    der_nor = der_nor.subtracted(
                        &at(&der_vec_nor, pderiv - ideriv, qderiv - jderiv).multiplied_scalar(
                            plib::bin(pderiv + iduref, ideriv + iduref)
                                * plib::bin(qderiv + idvref, jderiv + idvref)
                                * at_f64(&tab_norm, ideriv, jderiv),
                        ),
                    );
                }
            }
            der_nor = der_nor.divided(
                plib::bin(pderiv + iduref, iduref)
                    * plib::bin(qderiv + idvref, idvref)
                    * at_f64(&tab_norm, 0, 0),
            );
            set_vec(&mut der_vec_nor, pderiv, qderiv, der_nor);
        }
    }
    at(&der_vec_nor, nu, nv)
}
