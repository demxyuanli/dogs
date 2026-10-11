//! `Patch::MakeApprox` through the `mma2ce2_` checks
//! (`AdvApp2Var_Patch.cxx:764-880`, `AdvApp2Var_ApproxF2var.cxx:4243-4285`).
//! A finished square stores its equation. `Compute3DErrors` runs, then
//! `ConvertBS` returns the B-spline.

use super::context::App2VarContext;
use super::degree_error::FinishedPatch;
use super::framework::Framework;

pub(super) struct MakeApprox {
    pub sense: i32,
    pub finished: Option<FinishedPatch>,
}

/// `MakeApprox` until `mma2ce2_` (`AdvApp2Var_ApproxF2var.cxx:3958`).
pub(super) fn make_approx(
    frame: &Framework,
    ctx: &App2VarContext,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    num_dec: i32,
    sosotb: &[f64],
    disotb: &[f64],
    soditb: &[f64],
    diditb: &[f64],
) -> Result<MakeApprox, String> {
    let ndimen = ctx.total_dimension();
    let nbsesp = ctx.total_ssp();
    let mut nbpnt_u = ctx.u_roots.len() as i32;
    if ctx.ord_u > -1 {
        nbpnt_u -= 2;
    }
    let mut nbpnt_v = ctx.v_roots.len() as i32;
    if ctx.ord_v > -1 {
        nbpnt_v -= 2;
    }
    // `AdvApp2Var_Patch.cxx:799-809`: the floor is the neighbor iso degree.
    let ndmin_u = degree_floor(
        ctx.ord_u,
        frame.iso_v(u0, u1, v0)?.nb_coeff,
        frame.iso_v(u0, u1, v1)?.nb_coeff,
    );
    let ndmin_v = degree_floor(
        ctx.ord_v,
        frame.iso_u(u0, v0, v1)?.nb_coeff,
        frame.iso_u(u1, v0, v1)?.nb_coeff,
    );
    let (epsapr, epsfro) = tolerances(ctx, nbsesp)?;
    let cols_u = ctx.j_deg_u - 2 * ctx.ord_u - 1;
    let cols_v = ctx.j_deg_v - 2 * ctx.ord_v - 1;
    let u_gauss = super::gauss::gauss_weights(ctx.j_deg_u, nbpnt_u, ctx.ord_u, cols_u)?;
    let v_gauss = super::gauss::gauss_weights(ctx.j_deg_v, nbpnt_v, ctx.ord_v, cols_v)?;
    let j_max_u = super::jacobi_max::jacobi_max(ctx.j_deg_u, ctx.ord_u)?;
    let j_max_v = super::jacobi_max::jacobi_max(ctx.j_deg_v, ctx.ord_v)?;
    mma2ce2(
        num_dec,
        ndimen,
        nbsesp,
        ndmin_u,
        ndmin_v,
        ctx.lim_u - 1,
        ctx.lim_v - 1,
        ctx.j_deg_u,
        ctx.j_deg_v,
        ctx.ord_u,
        ctx.ord_v,
        nbpnt_u,
        nbpnt_v,
        sosotb,
        disotb,
        soditb,
        diditb,
        &u_gauss,
        &v_gauss,
        &j_max_u,
        &j_max_v,
        &epsapr,
        &epsfro,
        frame,
        u0,
        u1,
        v0,
        v1,
    )
}

/// `mma2ce2_` input checks (`AdvApp2Var_ApproxF2var.cxx:4243-4285`).
fn mma2ce2(
    num_dec: i32,
    ndimen: i32,
    nbsesp: i32,
    ndmin_u: i32,
    ndmin_v: i32,
    ndguli: i32,
    ndgvli: i32,
    ndjacu: i32,
    ndjacv: i32,
    iordru: i32,
    iordrv: i32,
    nbpnt_u: i32,
    nbpnt_v: i32,
    sosotb: &[f64],
    disotb: &[f64],
    soditb: &[f64],
    diditb: &[f64],
    gssutb: &[f64],
    gssvtb: &[f64],
    xmaxju: &[f64],
    xmaxjv: &[f64],
    epsapr: &[f64],
    epsfro: &[f64],
    frame: &Framework,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
) -> Result<MakeApprox, String> {
    if !(0..=5).contains(&num_dec)
        || (iordru << 1) + 1 > ndmin_u
        || ndmin_u > ndguli
        || ndguli >= ndjacu
        || (iordrv << 1) + 1 > ndmin_v
        || ndmin_v > ndgvli
        || ndgvli >= ndjacv
    {
        return Err("AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4926".into());
    }
    let itydec = 0;
    let ndegpu = ndmin_u;
    let ndegpv = ndmin_v;
    let errmax = vec![0.0; nbsesp.max(0) as usize];
    let vecerr = vec![0.0; (ndimen << 2).max(0) as usize];
    let mut patjac = vec![0.0; ((ndjacu + 1) * (ndjacv + 1) * ndimen).max(0) as usize];
    let i2rdu = (iordru + 1) << 1;
    let i2rdv = (iordrv + 1) << 1;
    let cols = ndjacu - i2rdu + 1;
    let pair_dim = nbpnt_v / 2 + 1;
    let impr_dim = nbpnt_v / 2;
    let mut chpair = vec![0.0; (pair_dim * cols * ndimen).max(0) as usize];
    let mut chimpr = vec![0.0; (impr_dim * cols * ndimen).max(0) as usize];
    if (1..=5).contains(&num_dec) {
        let minu = ndguli + 1;
        let maxu = ndjacu;
        let minv = ndgvli + 1;
        let maxv = ndjacv;
        if minu > maxu || minv > maxv {
            return Err("AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4926".into());
        }
        fill_u(
            minu, maxu, ndimen, nbpnt_u, nbpnt_v, i2rdu, cols, sosotb, disotb, soditb, diditb,
            gssutb, &mut chpair, &mut chimpr,
        )?;
        if nbsesp != 1 {
            return Err("UNPORTED AdvApp2Var_Iso::MakeApprox: AdvApp2Var_Iso.cxx:122".into());
        }
        let eps = *epsapr.first().ok_or_else(|| {
            "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4387".to_string()
        })?;
        let mut zone4 = 0.0;
        for jj in minv..=maxv {
            let igsv = jj - i2rdv;
            for nd in 1..=ndimen {
                mma2cfv(
                    jj, minu, maxu, nbpnt_v, igsv, nd, i2rdu, cols, ndjacu, ndjacv, gssvtb,
                    &chpair, &chimpr, &mut patjac,
                )?;
            }
            zone4 = super::degree_error::mma2er1(
                ndjacu, ndjacv, ndimen, minu, maxu, jj, jj, iordru, iordrv, xmaxju, xmaxjv,
                &patjac, zone4,
            )?;
            if zone4 > eps {
                let itydec = cut_sense(num_dec, 3)?;
                let _ = (
                    itydec, ndegpu, ndegpv, errmax, vecerr, patjac, chpair, chimpr, epsfro,
                );
                return Ok(MakeApprox { sense: itydec, finished: None });
            }
        }
        let minu2 = i2rdu;
        let maxu2 = ndguli;
        let mut zone2 = 0.0;
        if minu2 <= maxu2 {
            fill_u(
                minu2, maxu2, ndimen, nbpnt_u, nbpnt_v, i2rdu, cols, sosotb, disotb, soditb, diditb,
                gssutb, &mut chpair, &mut chimpr,
            )?;
            for jj in minv..=maxv {
                let igsv = jj - i2rdv;
                for nd in 1..=ndimen {
                    mma2cfv(
                        jj, minu2, maxu2, nbpnt_v, igsv, nd, i2rdu, cols, ndjacu, ndjacv, gssvtb,
                        &chpair, &chimpr, &mut patjac,
                    )?;
                }
            }
            zone2 = super::degree_error::mma2er1(
                ndjacu, ndjacv, ndimen, minu2, maxu2, minv, maxv, iordru, iordrv, xmaxju, xmaxjv,
                &patjac, 0.0,
            )?;
        }
        let minv3 = i2rdv;
        let maxv3 = ndgvli;
        let mut zone3 = 0.0;
        if minv3 <= maxv3 {
            for jj in minv3..=maxv3 {
                let igsv = jj - i2rdv;
                for nd in 1..=ndimen {
                    mma2cfv(
                        jj, minu, maxu, nbpnt_v, igsv, nd, i2rdu, cols, ndjacu, ndjacv, gssvtb,
                        &chpair, &chimpr, &mut patjac,
                    )?;
                }
            }
            zone3 = super::degree_error::mma2er1(
                ndjacu, ndjacv, ndimen, minu, maxu, minv3, maxv3, iordru, iordrv, xmaxju, xmaxjv,
                &patjac, 0.0,
            )?;
        }
        let vaux = [zone2, zone4, zone3];
        let err = super::border_error::mzsnorm(&vaux);
        if err > eps {
            let zv = super::border_error::mzsnorm(&vaux[..2]);
            let zu = super::border_error::mzsnorm(&vaux[1..]);
            let which = if zu > eps && zv > eps {
                3
            } else if zu > zv {
                1
            } else {
                2
            };
            let itydec = cut_sense(num_dec, which)?;
            let _ = (ndegpu, ndegpv, errmax, vecerr, patjac, chpair, chimpr, err);
            return Ok(MakeApprox { sense: itydec, finished: None });
        }
        let _ = (ndegpu, ndegpv, errmax, vecerr, chpair, chimpr);
        let finished = super::degree_error::canon_stop(
            ndjacu, ndjacv, ndimen, iordru, iordrv, ndmin_u, ndmin_v, ndguli, ndgvli, true, err,
            ndguli, ndgvli, eps, xmaxju, xmaxjv, &mut patjac, frame, u0, u1, v0, v1, nbsesp,
            epsapr, epsfro,
        )?;
        return Ok(MakeApprox { sense: 0, finished: Some(finished) });
    }
    let minu = i2rdu;
    let maxu = ndjacu;
    let minv = i2rdv;
    let maxv = ndjacv;
    fill_u(
        minu, maxu, ndimen, nbpnt_u, nbpnt_v, i2rdu, cols, sosotb, disotb, soditb, diditb, gssutb,
        &mut chpair, &mut chimpr,
    )?;
    if minv <= maxv {
        for nd in 1..=ndimen {
            for jj in minv..=maxv {
                mma2cfv(
                    jj, minu, maxu, nbpnt_v, jj - i2rdv, nd, i2rdu, cols, ndjacu, ndjacv, gssvtb,
                    &chpair, &chimpr, &mut patjac,
                )?;
            }
        }
    }
    if nbsesp != 1 {
        return Err("UNPORTED AdvApp2Var_Iso::MakeApprox: AdvApp2Var_Iso.cxx:122".into());
    }
    let eps = *epsapr.first().ok_or_else(|| {
        "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4828".to_string()
    })?;
    let mut errmax_v = super::degree_error::mma2er1(
        ndjacu, ndjacv, ndimen, i2rdu, ndjacu, ndgvli + 1, ndjacv, iordru, iordrv, xmaxju, xmaxjv,
        &patjac, 0.0,
    )?;
    if i2rdv <= ndgvli {
        errmax_v = super::degree_error::mma2er1(
            ndjacu, ndjacv, ndimen, ndguli + 1, ndjacu, i2rdv, ndgvli, iordru, iordrv, xmaxju,
            xmaxjv, &patjac, errmax_v,
        )?;
    }
    let _ = (itydec, ndegpu, ndegpv, errmax, vecerr, chpair, chimpr);
    let finished = super::degree_error::canon_stop(
        ndjacu,
        ndjacv,
        ndimen,
        iordru,
        iordrv,
        ndmin_u,
        ndmin_v,
        ndguli,
        ndgvli,
        errmax_v <= eps,
        errmax_v,
        ndguli,
        ndgvli,
        eps,
        xmaxju,
        xmaxjv,
        &mut patjac,
        frame,
        u0,
        u1,
        v0,
        v1,
        nbsesp,
        epsapr,
        epsfro,
    )?;
    Ok(MakeApprox { sense: 0, finished: Some(finished) })
}

fn cut_sense(num_dec: i32, which: i32) -> Result<i32, String> {
    if !(1..=5).contains(&num_dec) {
        return Err("AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4926".into());
    }
    let itydec = match which {
        1 => {
            if num_dec != 2 {
                1
            } else {
                2
            }
        }
        2 => {
            if num_dec != 1 {
                2
            } else {
                1
            }
        }
        3 => match num_dec {
            5 => 3,
            2 | 4 => 2,
            1 | 3 => 1,
            _ => {
                return Err(
                    "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4926".into(),
                );
            }
        },
        _ => {
            return Err("AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4926".into());
        }
    };
    Ok(itydec)
}

fn fill_u(
    minu: i32,
    maxu: i32,
    ndimen: i32,
    nbpnt_u: i32,
    nbpnt_v: i32,
    i2rdu: i32,
    cols: i32,
    sosotb: &[f64],
    disotb: &[f64],
    soditb: &[f64],
    diditb: &[f64],
    gssutb: &[f64],
    chpair: &mut [f64],
    chimpr: &mut [f64],
) -> Result<(), String> {
    if minu > maxu {
        return Ok(());
    }
    for nd in 1..=ndimen {
        for kk in minu..=maxu {
            mma2cfu(
                kk,
                nbpnt_u,
                nbpnt_v,
                nd,
                kk - i2rdu,
                cols,
                sosotb,
                disotb,
                soditb,
                diditb,
                gssutb,
                chpair,
                chimpr,
            )?;
        }
    }
    Ok(())
}

fn cfu_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2cfu_: AdvApp2Var_ApproxF2var.cxx:5126".into()
}

fn take_cfu(buf: &[f64], index: i32) -> Result<f64, String> {
    buf.get(index as usize).copied().ok_or_else(cfu_oob)
}

fn store_cfu(buf: &mut [f64], index: i32, value: f64) -> Result<(), String> {
    let slot = buf.get_mut(index as usize).ok_or_else(cfu_oob)?;
    *slot = value;
    Ok(())
}

/// `mma2cfu_` (`AdvApp2Var_ApproxF2var.cxx:4998-5224`).
fn mma2cfu(
    ndujac: i32,
    nbpnt_u: i32,
    nbpnt_v: i32,
    nd: i32,
    igsu: i32,
    cols: i32,
    sosotb: &[f64],
    disotb: &[f64],
    soditb: &[f64],
    diditb: &[f64],
    gssutb: &[f64],
    chpair: &mut [f64],
    chimpr: &mut [f64],
) -> Result<(), String> {
    let nptu2 = nbpnt_u / 2;
    let nptv2 = nbpnt_v / 2;
    let soso_dim = nptu2 + 1;
    let soso_dim2 = nptv2 + 1;
    let diso_dim = nptu2;
    let pair_dim = nptv2 + 1;
    let impr_dim = nptv2;
    let soso = |ii: i32, jj: i32| ii + jj * soso_dim + (nd - 1) * soso_dim * soso_dim2;
    let diso = |ii: i32, jj: i32| (ii - 1) + (jj - 1) * diso_dim + (nd - 1) * diso_dim * nptv2;
    let gss = |ii: i32| igsu * soso_dim + ii;
    let pair = |jj: i32| jj + igsu * pair_dim + (nd - 1) * pair_dim * cols;
    let impr = |jj: i32| igsu * impr_dim + (nd - 1) * impr_dim * cols + (jj - 1);
    if ndujac % 2 == 0 {
        for jj in 1..=nptv2 {
            let mut bid1 = 0.0;
            let mut bid2 = 0.0;
            for ii in 1..=nptu2 {
                let bid0 = take_cfu(gssutb, gss(ii))?;
                bid1 += take_cfu(sosotb, soso(ii, jj))? * bid0;
                bid2 += take_cfu(soditb, diso(ii, jj))? * bid0;
            }
            store_cfu(chpair, pair(jj), bid1)?;
            store_cfu(chimpr, impr(jj), bid2)?;
        }
    } else {
        for jj in 1..=nptv2 {
            let mut bid1 = 0.0;
            let mut bid2 = 0.0;
            for ii in 1..=nptu2 {
                let bid0 = take_cfu(gssutb, gss(ii))?;
                bid1 += take_cfu(disotb, diso(ii, jj))? * bid0;
                bid2 += take_cfu(diditb, soso(ii, jj))? * bid0;
            }
            store_cfu(chpair, pair(jj), bid1)?;
            store_cfu(chimpr, impr(jj), bid2)?;
        }
    }
    if nbpnt_u % 2 != 0 && ndujac % 2 == 0 {
        let bid0 = take_cfu(gssutb, gss(0))?;
        for jj in 1..=nptv2 {
            let even = take_cfu(chpair, pair(jj))? + take_cfu(sosotb, soso(0, jj))? * bid0;
            let odd = take_cfu(chimpr, impr(jj))? + take_cfu(diditb, soso(0, jj))? * bid0;
            store_cfu(chpair, pair(jj), even)?;
            store_cfu(chimpr, impr(jj), odd)?;
        }
    }
    if nbpnt_v % 2 != 0 {
        let mut bid1 = 0.0;
        if ndujac % 2 == 0 {
            for ii in 1..=nptu2 {
                bid1 += take_cfu(sosotb, soso(ii, 0))? * take_cfu(gssutb, gss(ii))?;
            }
        } else {
            for ii in 1..=nptu2 {
                bid1 += take_cfu(diditb, soso(ii, 0))? * take_cfu(gssutb, gss(ii))?;
            }
        }
        if nbpnt_u % 2 != 0 {
            bid1 += take_cfu(sosotb, soso(0, 0))? * take_cfu(gssutb, gss(0))?;
        }
        store_cfu(chpair, pair(0), bid1)?;
    }
    Ok(())
}

fn cfv_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2cfv_: AdvApp2Var_ApproxF2var.cxx:5330".into()
}

fn take_cfv(buf: &[f64], index: i32) -> Result<f64, String> {
    buf.get(index as usize).copied().ok_or_else(cfv_oob)
}

fn store_cfv(buf: &mut [f64], index: i32, value: f64) -> Result<(), String> {
    let slot = buf.get_mut(index as usize).ok_or_else(cfv_oob)?;
    *slot = value;
    Ok(())
}

/// `mma2cfv_` (`AdvApp2Var_ApproxF2var.cxx:5228-5379`).
fn mma2cfv(
    ndvjac: i32,
    mindgu: i32,
    maxdgu: i32,
    nbpnt_v: i32,
    igsv: i32,
    nd: i32,
    i2rdu: i32,
    cols: i32,
    ndjacu: i32,
    ndjacv: i32,
    gssvtb: &[f64],
    chpair: &[f64],
    chimpr: &[f64],
    patjac: &mut [f64],
) -> Result<(), String> {
    let nptv2 = nbpnt_v / 2;
    let pair_dim = nptv2 + 1;
    let impr_dim = nptv2;
    let pat_dim1 = ndjacu + 1;
    let pat_dim2 = ndjacv + 1;
    let gss = |jj: i32| igsv * pair_dim + jj;
    let pair = |ii: i32, jj: i32| jj + (ii - i2rdu) * pair_dim + (nd - 1) * pair_dim * cols;
    let impr = |ii: i32, jj: i32| {
        (ii - i2rdu) * impr_dim + (nd - 1) * impr_dim * cols + (jj - 1)
    };
    let pat = |ii: i32| ii + ndvjac * pat_dim1 + (nd - 1) * pat_dim1 * pat_dim2;
    if mindgu <= maxdgu && nptv2 >= 1 {
        if ndvjac % 2 == 0 {
            for ii in mindgu..=maxdgu {
                let mut bid1 = 0.0;
                for jj in 1..=nptv2 {
                    bid1 += take_cfv(chpair, pair(ii, jj))? * take_cfv(gssvtb, gss(jj))?;
                }
                store_cfv(patjac, pat(ii), bid1)?;
            }
        } else {
            for ii in mindgu..=maxdgu {
                let mut bid1 = 0.0;
                for jj in 1..=nptv2 {
                    bid1 += take_cfv(chimpr, impr(ii, jj))? * take_cfv(gssvtb, gss(jj))?;
                }
                store_cfv(patjac, pat(ii), bid1)?;
            }
        }
    }
    if nbpnt_v % 2 != 0 && ndvjac % 2 == 0 && mindgu <= maxdgu {
        let bid1 = take_cfv(gssvtb, gss(0))?;
        for ii in mindgu..=maxdgu {
            let sum = take_cfv(patjac, pat(ii))? + bid1 * take_cfv(chpair, pair(ii, 0))?;
            store_cfv(patjac, pat(ii), sum)?;
        }
    }
    Ok(())
}


fn degree_floor(order: i32, first_coeff: i32, second_coeff: i32) -> i32 {
    let mut ndmin = 1.max(2 * order + 1);
    ndmin = ndmin.max(first_coeff - 1);
    ndmin.max(second_coeff - 1)
}

fn tolerances(ctx: &App2VarContext, nbsesp: i32) -> Result<(Vec<f64>, Vec<f64>), String> {
    let ssp = nbsesp as usize;
    if ctx.internal_tol.len() != ssp
        || ctx.frontier_tol.len() != ssp
        || ctx.cutting_tol.len() != ssp
    {
        return Err("AdvApp2Var_Patch::MakeApprox: AdvApp2Var_Patch.cxx:817".into());
    }
    let mut epsapr = vec![0.0; ssp];
    let mut epsfro = vec![0.0; ssp * 8];
    for iesp in 1..=nbsesp {
        let at = (iesp - 1) as usize;
        epsapr[at] = ctx.internal_tol[at];
        for side in 0..4 {
            epsfro[at + side * ssp] = ctx.frontier_tol[at][side];
            epsfro[at + (side + 4) * ssp] = ctx.cutting_tol[at][side];
        }
    }
    Ok((epsapr, epsfro))
}
