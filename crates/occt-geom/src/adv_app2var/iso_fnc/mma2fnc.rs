//! `mma2fnc_` (`AdvApp2Var_ApproxF2var.cxx:7251`) over `BorderCut`.

use crate::Surface;

use super::super::context::App2VarContext;

use super::*;

/// One derivative order of `mma2fnc_`.
pub(crate) struct BorderCut {
    pub(crate) iercod: i32,
    pub(crate) contr1: Vec<f64>,
    pub(crate) contr2: Vec<f64>,
    pub(crate) errmax: f64,
    pub(crate) errmoy: f64,
    pub(crate) somtab: Vec<f64>,
    pub(crate) diftab: Vec<f64>,
    pub(crate) courbe: Vec<f64>,
    pub(crate) ncoeff: i32,
}

/// One derivative order of `mma2fnc_` through `mma1noc_`.
pub(crate) fn mma2fnc(
    surf: &dyn Surface,
    ctx: &App2VarContext,
    uvfonc: [f64; 4],
    tconst: f64,
    isofav: i32,
    iordre: i32,
    ideriv: i32,
    position: i32,
) -> Result<BorderCut, String> {
    let ndimen = ctx.total_dimension();
    if (uvfonc[1] - uvfonc[0]).abs() < eps3() || (uvfonc[3] - uvfonc[2]).abs() < eps3() {
        return Err(incoherent());
    }
    if isofav != 1 && isofav != 2 {
        return Err(incoherent());
    }

    let (nbroot, roots, ndgjac, ncflim) = if isofav == 2 {
        let mut nbroot = ctx.u_roots.len() as i32;
        if iordre > -1 {
            nbroot -= 2;
        }
        (nbroot, &ctx.u_roots, ctx.j_deg_u, ctx.lim_u)
    } else {
        let mut nbroot = ctx.v_roots.len() as i32;
        if iordre > -1 {
            nbroot -= 2;
        }
        (nbroot, &ctx.v_roots, ctx.j_deg_v, ctx.lim_v)
    };
    if nbroot < 0 || (roots.len() as i32) < nbroot || ndimen < 1 || ncflim < 1 {
        return Err(incoherent());
    }
    let rootlg = &roots[..nbroot as usize];
    let epsapr = eps_column(ctx, position)?;

    // Work-buffer sizes (`AdvApp2Var_ApproxF2var.cxx:6912-6928`).
    let isz1 = nbroot + 2;
    let ibid_points = ndimen * (nbroot + 2);
    let ibid_order = ((iordre + 1) << 1) * nbroot;
    let mut isz2 = ibid_points.max(ibid_order);
    let ibid_curve = (((ncflim - 1) / 2 + 1) << 1) * ndimen;
    isz2 = ibid_curve.max(isz2);
    let isz3 = ((iordre + 1) << 2) * (iordre + 1);
    let isz4 = (nbroot / 2 + 1) * (ndgjac + 1 - ((iordre + 1) << 1));
    let isz5 = (ndgjac + 1) * ndimen;
    if isz1 < 0 || isz2 < 0 || isz3 < 0 || isz4 < 0 || isz5 < 0 {
        return Err(alloc_failed());
    }

    let mut uvpav = [0.0; 4];
    let mut tabdec = [0.0; 2];
    if isofav == 1 {
        uvpav[0] = uvfonc[0];
        uvpav[1] = uvfonc[1];
        tabdec[0] = uvfonc[2];
        tabdec[1] = uvfonc[3];
    } else {
        tabdec[0] = uvfonc[0];
        tabdec[1] = uvfonc[1];
        uvpav[2] = uvfonc[2];
        uvpav[3] = uvfonc[3];
    }
    // `MakeApprox` passes `NBCRMX = 1` (`AdvApp2Var_Iso.cxx:116`), so the
    // cut arm at `mma2fnc_` `:7227-7236` never runs: `ier == -1` accepts
    // because `nupil == nbcrmx`, and `ier > 0` is an error.
    if isofav == 1 {
        uvpav[2] = tabdec[0];
        uvpav[3] = tabdec[1];
    } else {
        uvpav[0] = tabdec[0];
        uvpav[1] = tabdec[1];
    }

    let mut ttable = vec![0.0; isz1 as usize];
    mma1nop(nbroot, rootlg, uvpav, isofav, &mut ttable)?;

    let dim1 = nbroot / 2 + 1;
    let mut somtab = vec![0.0; (dim1 * ndimen) as usize];
    let mut diftab = vec![0.0; (dim1 * ndimen) as usize];
    let contr_len = ((iordre + 2).max(1) * ndimen) as usize;
    let mut contr1 = vec![0.0; contr_len];
    let mut contr2 = vec![0.0; contr_len];
    mma1fdi(
        surf,
        ndimen,
        uvpav,
        isofav,
        tconst,
        nbroot,
        &ttable,
        iordre,
        ideriv,
        &mut somtab,
        &mut diftab,
        &mut contr1,
        &mut contr2,
    )?;
    if iordre >= 0 {
        mma1cdi(
            ndimen,
            nbroot,
            rootlg,
            iordre,
            &contr1,
            &contr2,
            &mut somtab,
            &mut diftab,
        )?;
    }
    let mut crvjac = mma1jak(ndimen, nbroot, iordre, ndgjac, &somtab, &diftab)?;
    if iordre >= 0 {
        let hermit = super::super::constraints::mma1her(iordre)?;
        mma1cnt(ndimen, iordre, &contr1, &contr2, &hermit, ndgjac, &mut crvjac);
    }
    // `MakeApprox` passes one 3D subspace (`AdvApp2Var_Iso.cxx:121-122`).
    if ctx.total_ssp() != 1 || ndimen != 3 {
        return Err("UNPORTED AdvApp2Var_Iso::MakeApprox: AdvApp2Var_Iso.cxx:122".into());
    }
    let (ier, ncoeff, errmax, errmoy) = super::super::border_error::mma1fer(
        1,
        &[ndimen],
        iordre,
        ndgjac,
        &mut crvjac,
        ncflim,
        &epsapr,
    )?;
    let ncfja = ndgjac + 1;
    let packed = super::super::border_error::mmapcmp(ndimen, ncfja, ncoeff, &crvjac)?;
    let half = (ncoeff - 1) / 2 + 1;
    let mut canon = vec![0.0; (ndimen * ncoeff) as usize];
    for nd in 1..=ndimen {
        let src = (((nd - 1) << 1) * half) as usize;
        let dst = ((nd - 1) * ncoeff) as usize;
        super::super::jacobi_canon::mmjacan(
            iordre,
            ncoeff - 1,
            &packed[src..src + (half * 2) as usize],
            &mut canon[dst..dst + ncoeff as usize],
        )?;
    }
    let mut courbe = vec![0.0; (ncflim * ndimen) as usize];
    // `mma2fnc_` passes `ncoeff` as `NDIMEN` and the geometric dimension as
    // `NCOEFU` (`AdvApp2Var_ApproxF2var.cxx:7156-7164`).
    super::super::table_expand::mmfmca8(ncoeff, ndimen, 1, ncflim, ndimen, 1, &canon, &mut courbe)?;
    let uv11 = [-1.0, 1.0, -1.0, 1.0];
    if iordre >= 0 {
        for ii in 0..=iordre {
            let at = contr_at(ii, 1, ndimen);
            super::super::param_normalize::mma1noc(uv11, ndimen, ii, uvpav, isofav, ideriv, &mut contr1[at..])?;
            super::super::param_normalize::mma1noc(uv11, ndimen, ii, uvpav, isofav, ideriv, &mut contr2[at..])?;
        }
    }
    let span = (nbroot / 2 + 1) * ndimen;
    super::super::param_normalize::mma1noc(uv11, span, 0, uvpav, isofav, ideriv, &mut somtab)?;
    super::super::param_normalize::mma1noc(uv11, span, 0, uvpav, isofav, ideriv, &mut diftab)?;
    for nd in 1..=ndimen {
        let at = ((nd - 1) * ncflim) as usize;
        super::super::param_normalize::mma1noc(uv11, ncoeff, 0, uvpav, isofav, ideriv, &mut courbe[at..])?;
    }
    let _ = (isz2, isz3, isz4, isz5, ier);
    Ok(BorderCut {
        iercod: ier,
        contr1,
        contr2,
        errmax,
        errmoy,
        somtab,
        diftab,
        courbe,
        ncoeff,
    })
}
