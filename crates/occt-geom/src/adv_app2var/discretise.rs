//! `AdvApp2Var_Patch::Discretise` through `mma2ds1_`
//! (`AdvApp2Var_Patch.cxx:96-355`). `MakeApprox` cuts with `UpdateInU` and
//! `UpdateInV`. A finished square is stored. `Compute3DErrors` runs, then
//! `ConvertBS` returns the B-spline.

use crate::Surface;

use super::context::App2VarContext;
use super::framework::Framework;

fn coin_at(ndimen: i32, ord_u: i32, iu: i32, iv: i32, axis: i32) -> usize {
    (ndimen * iu + ndimen * (ord_u + 2) * iv + axis) as usize
}

fn scale_orders(table: &mut [f64], order: i32, size: usize, span: f64) {
    if order < 1 {
        return;
    }
    for iu in 1..=order {
        let rho = span.powi(iu);
        let start = (iu as usize) * size;
        for value in &mut table[start..start + size] {
            *value *= rho;
        }
    }
}

/// Corner constraints and the four frontier discretisations.
pub(super) fn discretise_corners(
    frame: &Framework,
    ctx: &App2VarContext,
    surf: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    max_patches: i32,
    nb_u: i32,
    nb_v: i32,
) -> Result<(super::make_approx::MakeApprox, f64, f64), String> {
    let ndimen = ctx.total_dimension();
    let ord_u = ctx.ord_u;
    let ord_v = ctx.ord_v;
    // `Discretise` hard-codes X/Y/Z (`AdvApp2Var_Patch.cxx:104-105`).
    // `IORDRU` / `IORDRV` of -1 are valid (`:115`, `:123`, `:185-188`).
    if ndimen != 3 {
        return Err("UNPORTED AdvApp2Var_Patch::Discretise: AdvApp2Var_Patch.cxx:104".into());
    }
    let mut nbpnt_u = ctx.u_roots.len() as i32;
    if ord_u > -1 {
        nbpnt_u -= 2;
    }
    let mut nbpnt_v = ctx.v_roots.len() as i32;
    if ord_v > -1 {
        nbpnt_v -= 2;
    }
    let coin_size = (ndimen * (ord_u + 2) * (ord_v + 2)) as usize;
    let mut coins = vec![0.0; coin_size * 4];
    let du = (u1 - u0) / 2.0;
    let dv = (v1 - v0) / 2.0;
    let corners = [
        frame.node_at(u0, v0)?,
        frame.node_at(u1, v0)?,
        frame.node_at(u0, v1)?,
        frame.node_at(u1, v1)?,
    ];
    for iu in 0..=ord_u {
        for iv in 0..=ord_v {
            let rho = du.powi(iu) * dv.powi(iv);
            for (corner, node) in corners.iter().enumerate() {
                let point = node.point(iu, iv)?;
                for axis in 0..3 {
                    let at = corner * coin_size + coin_at(ndimen, ord_u, iu, iv, axis);
                    coins[at] = rho * point[axis as usize];
                }
            }
        }
    }

    let u_size = ((1 + nbpnt_v / 2) * ndimen) as usize;
    let u_len = u_size * (ord_u + 1) as usize;
    let mut su0 = frame.iso_u(u0, v0, v1)?.somtab.clone();
    let mut du0 = frame.iso_u(u0, v0, v1)?.diftab.clone();
    let mut su1 = frame.iso_u(u1, v0, v1)?.somtab.clone();
    let mut du1 = frame.iso_u(u1, v0, v1)?.diftab.clone();
    if su0.len() != u_len || du0.len() != u_len || su1.len() != u_len || du1.len() != u_len {
        return Err("AdvApp2Var_Patch::Discretise: AdvApp2Var_Patch.cxx:189".into());
    }
    scale_orders(&mut su0, ord_u, u_size, du);
    scale_orders(&mut du0, ord_u, u_size, du);
    scale_orders(&mut su1, ord_u, u_size, du);
    scale_orders(&mut du1, ord_u, u_size, du);

    let v_size = ((1 + nbpnt_u / 2) * ndimen) as usize;
    let v_len = v_size * (ord_v + 1) as usize;
    let mut sv0 = frame.iso_v(u0, u1, v0)?.somtab.clone();
    let mut dv0 = frame.iso_v(u0, u1, v0)?.diftab.clone();
    let mut sv1 = frame.iso_v(u0, u1, v1)?.somtab.clone();
    let mut dv1 = frame.iso_v(u0, u1, v1)?.diftab.clone();
    if sv0.len() != v_len || dv0.len() != v_len || sv1.len() != v_len || dv1.len() != v_len {
        return Err("AdvApp2Var_Patch::Discretise: AdvApp2Var_Patch.cxx:234".into());
    }
    scale_orders(&mut sv0, ord_v, v_size, dv);
    scale_orders(&mut dv0, ord_v, v_size, dv);
    scale_orders(&mut sv1, ord_v, v_size, dv);
    scale_orders(&mut dv1, ord_v, v_size, dv);

    if (ctx.u_roots.len() as i32) < nbpnt_u || (ctx.v_roots.len() as i32) < nbpnt_v {
        return Err("AdvApp2Var_Patch::Discretise: AdvApp2Var_Patch.cxx:297".into());
    }
    let mut tables = mma2cdi(
        ndimen,
        nbpnt_u,
        &ctx.u_roots[..nbpnt_u as usize],
        nbpnt_v,
        &ctx.v_roots[..nbpnt_v as usize],
        ord_u,
        ord_v,
        &su0,
        &su1,
        &du0,
        &du1,
        &sv0,
        &sv1,
        &dv0,
        &dv1,
        &coins,
    )?;
    super::square_disc::mma2ds1(
        surf,
        ctx.fav,
        ndimen,
        [u0, u1],
        [v0, v1],
        nbpnt_u,
        &ctx.u_roots[..nbpnt_u as usize],
        nbpnt_v,
        &ctx.v_roots[..nbpnt_v as usize],
        &mut tables.0,
        &mut tables.2,
        &mut tables.1,
        &mut tables.3,
    )?;
    let (num_dec, u_dec, v_dec) =
        super::compute_patches::authorize_cuts(surf, u0, u1, v0, v1, nb_u, nb_v, max_patches)?;
    let step = super::make_approx::make_approx(
        frame,
        ctx,
        u0,
        u1,
        v0,
        v1,
        num_dec,
        &tables.0,
        &tables.2,
        &tables.1,
        &tables.3,
    )?;
    Ok((step, u_dec, v_dec))
}

/// `mma2cdi_` through `mma2cd1_` (`AdvApp2Var_ApproxF2var.cxx:3552-3686`).
fn mma2cdi(
    ndimen: i32,
    nbpnt_u: i32,
    uroot: &[f64],
    nbpnt_v: i32,
    vroot: &[f64],
    iordru: i32,
    iordrv: i32,
    sotbu1: &[f64],
    sotbu2: &[f64],
    ditbu1: &[f64],
    ditbu2: &[f64],
    sotbv1: &[f64],
    sotbv2: &[f64],
    ditbv1: &[f64],
    ditbv2: &[f64],
    coins: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), String> {
    if !(-1..=2).contains(&iordru) || !(-1..=2).contains(&iordrv) {
        return Err("AdvApp2Var_ApproxF2var::mma2cdi_: AdvApp2Var_ApproxF2var.cxx:3561".into());
    }
    let nu = nbpnt_u / 2 + 1;
    let nv = nbpnt_v / 2 + 1;
    let mut sosotb = vec![0.0; (nu * nv * ndimen) as usize];
    let mut diditb = vec![0.0; (nu * nv * ndimen) as usize];
    let half_u = nbpnt_u / 2;
    let half_v = nbpnt_v / 2;
    let mut soditb = vec![0.0; (half_u * half_v * ndimen) as usize];
    let mut disotb = vec![0.0; (half_u * half_v * ndimen) as usize];
    if iordru == -1 && iordrv == -1 {
        return Ok((sosotb, soditb, disotb, diditb));
    }
    let mut hermit_u = Vec::new();
    if (0..=2).contains(&iordru) {
        hermit_u = super::iso_fnc::mma1her(iordru)?;
        mma2cd3(
            ndimen,
            nbpnt_u,
            uroot,
            nbpnt_v,
            iordru,
            sotbu1,
            sotbu2,
            ditbu1,
            ditbu2,
            &hermit_u,
            &mut sosotb,
            &mut soditb,
            &mut disotb,
            &mut diditb,
        )?;
    }
    let mut hermit_v = Vec::new();
    if (0..=2).contains(&iordrv) {
        hermit_v = super::iso_fnc::mma1her(iordrv)?;
        mma2cd2(
            ndimen,
            nbpnt_u,
            nbpnt_v,
            vroot,
            iordrv,
            sotbv1,
            sotbv2,
            ditbv1,
            ditbv2,
            &hermit_v,
            &mut sosotb,
            &mut soditb,
            &mut disotb,
            &mut diditb,
        )?;
    }
    if iordru >= 0 && iordrv >= 0 {
        mma2cd1(
            ndimen,
            nbpnt_u,
            uroot,
            nbpnt_v,
            vroot,
            iordru,
            iordrv,
            coins,
            &hermit_u,
            &hermit_v,
            &mut sosotb,
            &mut soditb,
            &mut disotb,
            &mut diditb,
        )?;
    }
    Ok((sosotb, soditb, disotb, diditb))
}

fn cd3_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2cd3_: AdvApp2Var_ApproxF2var.cxx:3243".into()
}

fn take(buf: &[f64], index: i32) -> Result<f64, String> {
    buf.get(index as usize).copied().ok_or_else(cd3_oob)
}

fn store(buf: &mut [f64], index: i32, value: f64) -> Result<(), String> {
    let slot = buf.get_mut(index as usize).ok_or_else(cd3_oob)?;
    *slot = value;
    Ok(())
}

/// `mma2cd3_` (`AdvApp2Var_ApproxF2var.cxx:3043-3362`).
fn mma2cd3(
    ndimen: i32,
    nbpnt_u: i32,
    uroot: &[f64],
    nbpnt_v: i32,
    iordru: i32,
    sotbu1: &[f64],
    sotbu2: &[f64],
    ditbu1: &[f64],
    ditbu2: &[f64],
    uhermt: &[f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
) -> Result<(), String> {
    let ncfhu = (iordru + 1) << 1;
    let mut fpntab = vec![0.0; (ncfhu * nbpnt_u) as usize];
    for ii in 1..=ncfhu {
        let row = ((ii - 1) * ncfhu) as usize;
        if row + ncfhu as usize > uhermt.len() {
            return Err(cd3_oob());
        }
        for kk in 1..=nbpnt_u {
            let root = *uroot.get((kk - 1) as usize).ok_or_else(cd3_oob)?;
            let value = super::iso_fnc::mmmpocur(ncfhu, &uhermt[row..], root);
            fpntab[((kk - 1) + (ii - 1) * nbpnt_u) as usize] = value;
        }
    }
    let nvroo = nbpnt_v / 2;
    let nuroo = nbpnt_u / 2;
    let border_dim = nvroo + 1;
    let soso_dim1 = nuroo + 1;
    let soso_dim2 = nvroo + 1;
    let diso_dim1 = nuroo;
    let diso_dim2 = nvroo;
    let border = |jj: i32, nd: i32, ii: i32| jj + (nd - 1) * border_dim + (ii - 1) * border_dim * ndimen;
    let soso = |kk: i32, jj: i32, nd: i32| kk + jj * soso_dim1 + (nd - 1) * soso_dim1 * soso_dim2;
    let diso = |kk: i32, jj: i32, nd: i32| {
        (kk - 1) + (jj - 1) * diso_dim1 + (nd - 1) * diso_dim1 * diso_dim2
    };
    let fp = |kk: i32, ii: i32| (kk - 1) + (ii - 1) * nbpnt_u;
    for nd in 1..=ndimen {
        for ii in 1..=iordru + 1 {
            for jj in 1..=nvroo {
                let bid1 = take(sotbu1, border(jj, nd, ii))?;
                let bid2 = take(sotbu2, border(jj, nd, ii))?;
                let bid3 = take(ditbu1, border(jj, nd, ii))?;
                let bid4 = take(ditbu2, border(jj, nd, ii))?;
                for kk in 1..=nuroo {
                    let kkp = (nbpnt_u + 1) / 2 + kk;
                    let kkm = nuroo - kk + 1;
                    let odd_p = take(&fpntab, fp(kkp, (ii << 1) - 1))?;
                    let odd_m = take(&fpntab, fp(kkm, (ii << 1) - 1))?;
                    let even_p = take(&fpntab, fp(kkp, ii << 1))?;
                    let even_m = take(&fpntab, fp(kkm, ii << 1))?;
                    let soso_i = soso(kk, jj, nd);
                    let diso_i = diso(kk, jj, nd);
                    store(sosotb, soso_i, take(sosotb, soso_i)? - bid1 * (odd_p + odd_m) - bid2 * (even_p + even_m))?;
                    store(disotb, diso_i, take(disotb, diso_i)? - bid1 * (odd_p - odd_m) - bid2 * (even_p - even_m))?;
                    store(soditb, diso_i, take(soditb, diso_i)? - bid3 * (odd_p + odd_m) - bid4 * (even_p + even_m))?;
                    store(diditb, soso_i, take(diditb, soso_i)? - bid3 * (odd_p - odd_m) - bid4 * (even_p - even_m))?;
                }
            }
        }
        if nbpnt_u % 2 == 1 {
            for ii in 1..=iordru + 1 {
                for jj in 1..=nvroo {
                    let bid1 = take(sotbu1, border(jj, nd, ii))? * take(&fpntab, fp(nuroo + 1, (ii << 1) - 1))?
                        + take(sotbu2, border(jj, nd, ii))? * take(&fpntab, fp(nuroo + 1, ii << 1))?;
                    let at = soso(0, jj, nd);
                    store(sosotb, at, take(sosotb, at)? - bid1)?;
                    let bid2 = take(ditbu1, border(jj, nd, ii))? * take(&fpntab, fp(nuroo + 1, (ii << 1) - 1))?
                        + take(ditbu2, border(jj, nd, ii))? * take(&fpntab, fp(nuroo + 1, ii << 1))?;
                    store(diditb, at, take(diditb, at)? - bid2)?;
                }
            }
        }
        if nbpnt_v % 2 == 1 {
            for ii in 1..=iordru + 1 {
                for kk in 1..=nuroo {
                    let kkp = (nbpnt_u + 1) / 2 + kk;
                    let kkm = nuroo - kk + 1;
                    let odd_p = take(&fpntab, fp(kkp, (ii << 1) - 1))?;
                    let odd_m = take(&fpntab, fp(kkm, (ii << 1) - 1))?;
                    let even_p = take(&fpntab, fp(kkp, ii << 1))?;
                    let even_m = take(&fpntab, fp(kkm, ii << 1))?;
                    let bid1 = take(sotbu1, border(0, nd, ii))? * (odd_p + odd_m)
                        + take(sotbu2, border(0, nd, ii))? * (even_p + even_m);
                    let at = soso(kk, 0, nd);
                    store(sosotb, at, take(sosotb, at)? - bid1)?;
                    let bid2 = take(sotbu1, border(0, nd, ii))? * (odd_p - odd_m)
                        + take(sotbu2, border(0, nd, ii))? * (even_p - even_m);
                    let at_d = soso(kk, 0, nd);
                    store(diditb, at_d, take(diditb, at_d)? - bid2)?;
                }
            }
        }
        if nbpnt_u % 2 == 1 && nbpnt_v % 2 == 1 {
            for ii in 1..=iordru + 1 {
                let bid1 = take(sotbu1, border(0, nd, ii))? * take(&fpntab, fp(nuroo + 1, (ii << 1) - 1))?
                    + take(sotbu2, border(0, nd, ii))? * take(&fpntab, fp(nuroo + 1, ii << 1))?;
                let at = soso(0, 0, nd);
                store(sosotb, at, take(sosotb, at)? - bid1)?;
            }
        }
    }
    Ok(())
}

fn cd2_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2cd2_: AdvApp2Var_ApproxF2var.cxx:2911".into()
}

fn take2(buf: &[f64], index: i32) -> Result<f64, String> {
    buf.get(index as usize).copied().ok_or_else(cd2_oob)
}

fn store2(buf: &mut [f64], index: i32, value: f64) -> Result<(), String> {
    let slot = buf.get_mut(index as usize).ok_or_else(cd2_oob)?;
    *slot = value;
    Ok(())
}

/// `mma2cd2_` (`AdvApp2Var_ApproxF2var.cxx:2723-3039`).
fn mma2cd2(
    ndimen: i32,
    nbpnt_u: i32,
    nbpnt_v: i32,
    vroot: &[f64],
    iordrv: i32,
    sotbv1: &[f64],
    sotbv2: &[f64],
    ditbv1: &[f64],
    ditbv2: &[f64],
    vhermt: &[f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
) -> Result<(), String> {
    let ncfhv = (iordrv + 1) << 1;
    let mut fpntab = vec![0.0; (ncfhv * nbpnt_v) as usize];
    for ii in 1..=ncfhv {
        let row = ((ii - 1) * ncfhv) as usize;
        if row + ncfhv as usize > vhermt.len() {
            return Err(cd2_oob());
        }
        for jj in 1..=nbpnt_v {
            let root = *vroot.get((jj - 1) as usize).ok_or_else(cd2_oob)?;
            let value = super::iso_fnc::mmmpocur(ncfhv, &vhermt[row..], root);
            fpntab[((jj - 1) + (ii - 1) * nbpnt_v) as usize] = value;
        }
    }
    let nuroo = nbpnt_u / 2;
    let nvroo = nbpnt_v / 2;
    let border_dim = nuroo + 1;
    let soso_dim1 = nuroo + 1;
    let soso_dim2 = nvroo + 1;
    let diso_dim1 = nuroo;
    let diso_dim2 = nvroo;
    let border =
        |kk: i32, nd: i32, ii: i32| kk + (nd - 1) * border_dim + (ii - 1) * border_dim * ndimen;
    let soso = |kk: i32, jj: i32, nd: i32| kk + jj * soso_dim1 + (nd - 1) * soso_dim1 * soso_dim2;
    let diso = |kk: i32, jj: i32, nd: i32| {
        (kk - 1) + (jj - 1) * diso_dim1 + (nd - 1) * diso_dim1 * diso_dim2
    };
    let fp = |jj: i32, ii: i32| (jj - 1) + (ii - 1) * nbpnt_v;
    for nd in 1..=ndimen {
        for ii in 1..=iordrv + 1 {
            for kk in 1..=nuroo {
                let bid1 = take2(sotbv1, border(kk, nd, ii))?;
                let bid2 = take2(sotbv2, border(kk, nd, ii))?;
                let bid3 = take2(ditbv1, border(kk, nd, ii))?;
                let bid4 = take2(ditbv2, border(kk, nd, ii))?;
                for jj in 1..=nvroo {
                    let jjp = (nbpnt_v + 1) / 2 + jj;
                    let jjm = nvroo - jj + 1;
                    let odd_p = take2(&fpntab, fp(jjp, (ii << 1) - 1))?;
                    let odd_m = take2(&fpntab, fp(jjm, (ii << 1) - 1))?;
                    let even_p = take2(&fpntab, fp(jjp, ii << 1))?;
                    let even_m = take2(&fpntab, fp(jjm, ii << 1))?;
                    let soso_i = soso(kk, jj, nd);
                    let diso_i = diso(kk, jj, nd);
                    store2(
                        sosotb,
                        soso_i,
                        take2(sosotb, soso_i)? - bid1 * (odd_p + odd_m) - bid2 * (even_p + even_m),
                    )?;
                    store2(
                        disotb,
                        diso_i,
                        take2(disotb, diso_i)? - bid3 * (odd_p + odd_m) - bid4 * (even_p + even_m),
                    )?;
                    store2(
                        soditb,
                        diso_i,
                        take2(soditb, diso_i)? - bid1 * (odd_p - odd_m) - bid2 * (even_p - even_m),
                    )?;
                    store2(
                        diditb,
                        soso_i,
                        take2(diditb, soso_i)? - bid3 * (odd_p - odd_m) - bid4 * (even_p - even_m),
                    )?;
                }
            }
        }
        if nbpnt_v % 2 == 1 {
            for ii in 1..=iordrv + 1 {
                for kk in 1..=nuroo {
                    let bid1 = take2(sotbv1, border(kk, nd, ii))?
                        * take2(&fpntab, fp(nvroo + 1, (ii << 1) - 1))?
                        + take2(sotbv2, border(kk, nd, ii))?
                            * take2(&fpntab, fp(nvroo + 1, ii << 1))?;
                    let at = soso(kk, 0, nd);
                    store2(sosotb, at, take2(sosotb, at)? - bid1)?;
                    let bid2 = take2(ditbv1, border(kk, nd, ii))?
                        * take2(&fpntab, fp(nvroo + 1, (ii << 1) - 1))?
                        + take2(ditbv2, border(kk, nd, ii))?
                            * take2(&fpntab, fp(nvroo + 1, ii << 1))?;
                    store2(diditb, at, take2(diditb, at)? - bid2)?;
                }
            }
        }
        if nbpnt_u % 2 == 1 {
            for ii in 1..=iordrv + 1 {
                for jj in 1..=nvroo {
                    let jjp = (nbpnt_v + 1) / 2 + jj;
                    let jjm = nvroo - jj + 1;
                    let odd_p = take2(&fpntab, fp(jjp, (ii << 1) - 1))?;
                    let odd_m = take2(&fpntab, fp(jjm, (ii << 1) - 1))?;
                    let even_p = take2(&fpntab, fp(jjp, ii << 1))?;
                    let even_m = take2(&fpntab, fp(jjm, ii << 1))?;
                    let bid1 = take2(sotbv1, border(0, nd, ii))? * (odd_p + odd_m)
                        + take2(sotbv2, border(0, nd, ii))? * (even_p + even_m);
                    let at = soso(0, jj, nd);
                    store2(sosotb, at, take2(sosotb, at)? - bid1)?;
                    let bid2 = take2(sotbv1, border(0, nd, ii))? * (odd_p - odd_m)
                        + take2(sotbv2, border(0, nd, ii))? * (even_p - even_m);
                    let at_d = soso(jj, 0, nd);
                    store2(diditb, at_d, take2(diditb, at_d)? - bid2)?;
                }
            }
        }
        if nbpnt_u % 2 == 1 && nbpnt_v % 2 == 1 {
            for ii in 1..=iordrv + 1 {
                let bid1 = take2(sotbv1, border(0, nd, ii))?
                    * take2(&fpntab, fp(nvroo + 1, (ii << 1) - 1))?
                    + take2(sotbv2, border(0, nd, ii))? * take2(&fpntab, fp(nvroo + 1, ii << 1))?;
                let at = soso(0, 0, nd);
                store2(sosotb, at, take2(sosotb, at)? - bid1)?;
            }
        }
    }
    Ok(())
}

fn cd1_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2cd1_: AdvApp2Var_ApproxF2var.cxx:2589".into()
}

fn take1(buf: &[f64], index: i32) -> Result<f64, String> {
    buf.get(index as usize).copied().ok_or_else(cd1_oob)
}

fn store1(buf: &mut [f64], index: i32, value: f64) -> Result<(), String> {
    let slot = buf.get_mut(index as usize).ok_or_else(cd1_oob)?;
    *slot = value;
    Ok(())
}

/// `mma2cd1_` (`AdvApp2Var_ApproxF2var.cxx:2365-2719`).
fn mma2cd1(
    ndimen: i32,
    nbpnt_u: i32,
    uroot: &[f64],
    nbpnt_v: i32,
    vroot: &[f64],
    iordru: i32,
    iordrv: i32,
    coins: &[f64],
    uhermt: &[f64],
    vhermt: &[f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
) -> Result<(), String> {
    let ncfhu = (iordru + 1) << 1;
    let mut fpntbu = vec![0.0; (ncfhu * nbpnt_u) as usize];
    for ii in 1..=ncfhu {
        let row = ((ii - 1) * ncfhu) as usize;
        if row + ncfhu as usize > uhermt.len() {
            return Err(cd1_oob());
        }
        for ll in 1..=nbpnt_u {
            let root = *uroot.get((ll - 1) as usize).ok_or_else(cd1_oob)?;
            let value = super::iso_fnc::mmmpocur(ncfhu, &uhermt[row..], root);
            fpntbu[((ll - 1) + (ii - 1) * nbpnt_u) as usize] = value;
        }
    }
    let ncfhv = (iordrv + 1) << 1;
    let mut fpntbv = vec![0.0; (ncfhv * nbpnt_v) as usize];
    for jj in 1..=ncfhv {
        let row = ((jj - 1) * ncfhv) as usize;
        if row + ncfhv as usize > vhermt.len() {
            return Err(cd1_oob());
        }
        for kk in 1..=nbpnt_v {
            let root = *vroot.get((kk - 1) as usize).ok_or_else(cd1_oob)?;
            let value = super::iso_fnc::mmmpocur(ncfhv, &vhermt[row..], root);
            fpntbv[((kk - 1) + (jj - 1) * nbpnt_v) as usize] = value;
        }
    }
    let nuroo = nbpnt_u / 2;
    let nvroo = nbpnt_v / 2;
    let corner = ndimen * (iordru + 2) * (iordrv + 2);
    let stride_v = ndimen * (iordru + 2);
    let soso_dim1 = nuroo + 1;
    let soso_dim2 = nvroo + 1;
    let diso_dim1 = nuroo;
    let diso_dim2 = nvroo;
    let contr = |corner_i: i32, ii: i32, jj: i32, nd: i32| {
        corner_i * corner + (nd - 1) + (ii - 1) * ndimen + (jj - 1) * stride_v
    };
    let soso = |ll: i32, kk: i32, nd: i32| ll + kk * soso_dim1 + (nd - 1) * soso_dim1 * soso_dim2;
    let diso = |ll: i32, kk: i32, nd: i32| {
        (ll - 1) + (kk - 1) * diso_dim1 + (nd - 1) * diso_dim1 * diso_dim2
    };
    let fpu = |ll: i32, ii: i32| (ll - 1) + (ii - 1) * nbpnt_u;
    let fpv = |kk: i32, jj: i32| (kk - 1) + (jj - 1) * nbpnt_v;
    for nd in 1..=ndimen {
        for jj in 1..=iordrv + 1 {
            for ii in 1..=iordru + 1 {
                let bid1 = take1(coins, contr(0, ii, jj, nd))?;
                let bid2 = take1(coins, contr(1, ii, jj, nd))?;
                let bid3 = take1(coins, contr(2, ii, jj, nd))?;
                let bid4 = take1(coins, contr(3, ii, jj, nd))?;
                for kk in 1..=nvroo {
                    let kkp = (nbpnt_v + 1) / 2 + kk;
                    let kkm = nvroo - kk + 1;
                    let sov1 = take1(&fpntbv, fpv(kkp, (jj << 1) - 1))?
                        + take1(&fpntbv, fpv(kkm, (jj << 1) - 1))?;
                    let div1 = take1(&fpntbv, fpv(kkp, (jj << 1) - 1))?
                        - take1(&fpntbv, fpv(kkm, (jj << 1) - 1))?;
                    let sov2 = take1(&fpntbv, fpv(kkp, jj << 1))? + take1(&fpntbv, fpv(kkm, jj << 1))?;
                    let div2 = take1(&fpntbv, fpv(kkp, jj << 1))? - take1(&fpntbv, fpv(kkm, jj << 1))?;
                    for ll in 1..=nuroo {
                        let llp = (nbpnt_u + 1) / 2 + ll;
                        let llm = nuroo - ll + 1;
                        let sou1 = take1(&fpntbu, fpu(llp, (ii << 1) - 1))?
                            + take1(&fpntbu, fpu(llm, (ii << 1) - 1))?;
                        let diu1 = take1(&fpntbu, fpu(llp, (ii << 1) - 1))?
                            - take1(&fpntbu, fpu(llm, (ii << 1) - 1))?;
                        let sou2 =
                            take1(&fpntbu, fpu(llp, ii << 1))? + take1(&fpntbu, fpu(llm, ii << 1))?;
                        let diu2 =
                            take1(&fpntbu, fpu(llp, ii << 1))? - take1(&fpntbu, fpu(llm, ii << 1))?;
                        let soso_i = soso(ll, kk, nd);
                        let diso_i = diso(ll, kk, nd);
                        store1(
                            sosotb,
                            soso_i,
                            take1(sosotb, soso_i)?
                                - bid1 * sou1 * sov1
                                - bid2 * sou2 * sov1
                                - bid3 * sou1 * sov2
                                - bid4 * sou2 * sov2,
                        )?;
                        store1(
                            soditb,
                            diso_i,
                            take1(soditb, diso_i)?
                                - bid1 * sou1 * div1
                                - bid2 * sou2 * div1
                                - bid3 * sou1 * div2
                                - bid4 * sou2 * div2,
                        )?;
                        store1(
                            disotb,
                            diso_i,
                            take1(disotb, diso_i)?
                                - bid1 * diu1 * sov1
                                - bid2 * diu2 * sov1
                                - bid3 * diu1 * sov2
                                - bid4 * diu2 * sov2,
                        )?;
                        store1(
                            diditb,
                            soso_i,
                            take1(diditb, soso_i)?
                                - bid1 * diu1 * div1
                                - bid2 * diu2 * div1
                                - bid3 * diu1 * div2
                                - bid4 * diu2 * div2,
                        )?;
                    }
                }
                if nbpnt_u % 2 == 1 {
                    let sou1 = take1(&fpntbu, fpu(nuroo + 1, (ii << 1) - 1))?;
                    let sou2 = take1(&fpntbu, fpu(nuroo + 1, ii << 1))?;
                    for kk in 1..=nvroo {
                        let kkp = (nbpnt_v + 1) / 2 + kk;
                        let kkm = nvroo - kk + 1;
                        let sov1 = take1(&fpntbv, fpv(kkp, (jj << 1) - 1))?
                            + take1(&fpntbv, fpv(kkm, (jj << 1) - 1))?;
                        let div1 = take1(&fpntbv, fpv(kkp, (jj << 1) - 1))?
                            - take1(&fpntbv, fpv(kkm, (jj << 1) - 1))?;
                        let sov2 =
                            take1(&fpntbv, fpv(kkp, jj << 1))? + take1(&fpntbv, fpv(kkm, jj << 1))?;
                        let div2 =
                            take1(&fpntbv, fpv(kkp, jj << 1))? - take1(&fpntbv, fpv(kkm, jj << 1))?;
                        let at = soso(0, kk, nd);
                        store1(
                            sosotb,
                            at,
                            take1(sosotb, at)?
                                - bid1 * sou1 * sov1
                                - bid2 * sou2 * sov1
                                - bid3 * sou1 * sov2
                                - bid4 * sou2 * sov2,
                        )?;
                        store1(
                            diditb,
                            at,
                            take1(diditb, at)?
                                - bid1 * sou1 * div1
                                - bid2 * sou2 * div1
                                - bid3 * sou1 * div2
                                - bid4 * sou2 * div2,
                        )?;
                    }
                }
                if nbpnt_v % 2 == 1 {
                    let sov1 = take1(&fpntbv, fpv(nvroo + 1, (jj << 1) - 1))?;
                    let sov2 = take1(&fpntbv, fpv(nvroo + 1, jj << 1))?;
                    for ll in 1..=nuroo {
                        let llp = (nbpnt_u + 1) / 2 + ll;
                        let llm = nuroo - ll + 1;
                        let sou1 = take1(&fpntbu, fpu(llp, (ii << 1) - 1))?
                            + take1(&fpntbu, fpu(llm, (ii << 1) - 1))?;
                        let diu1 = take1(&fpntbu, fpu(llp, (ii << 1) - 1))?
                            - take1(&fpntbu, fpu(llm, (ii << 1) - 1))?;
                        let sou2 =
                            take1(&fpntbu, fpu(llp, ii << 1))? + take1(&fpntbu, fpu(llm, ii << 1))?;
                        let diu2 =
                            take1(&fpntbu, fpu(llp, ii << 1))? - take1(&fpntbu, fpu(llm, ii << 1))?;
                        let at = soso(ll, 0, nd);
                        store1(
                            sosotb,
                            at,
                            take1(sosotb, at)?
                                - bid1 * sou1 * sov1
                                - bid2 * sou2 * sov1
                                - bid3 * sou1 * sov2
                                - bid4 * sou2 * sov2,
                        )?;
                        store1(
                            diditb,
                            at,
                            take1(diditb, at)?
                                - bid1 * diu1 * sov1
                                - bid2 * diu2 * sov1
                                - bid3 * diu1 * sov2
                                - bid4 * diu2 * sov2,
                        )?;
                    }
                }
                if nbpnt_u % 2 == 1 && nbpnt_v % 2 == 1 {
                    let sou1 = take1(&fpntbu, fpu(nuroo + 1, (ii << 1) - 1))?;
                    let sou2 = take1(&fpntbu, fpu(nuroo + 1, ii << 1))?;
                    let sov1 = take1(&fpntbv, fpv(nvroo + 1, (jj << 1) - 1))?;
                    let sov2 = take1(&fpntbv, fpv(nvroo + 1, jj << 1))?;
                    let term = bid1 * sou1 * sov1
                        + bid2 * sou2 * sov1
                        + bid3 * sou1 * sov2
                        + bid4 * sou2 * sov2;
                    let at = soso(0, 0, nd);
                    store1(sosotb, at, take1(sosotb, at)? - term)?;
                    store1(diditb, at, take1(diditb, at)? - term)?;
                }
            }
        }
    }
    Ok(())
}
