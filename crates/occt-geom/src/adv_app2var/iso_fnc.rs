//! `AdvApp2Var_Iso::MakeApprox` (`AdvApp2Var_Iso.cxx:106-376`) through
//! `mma2fnc_` (`AdvApp2Var_ApproxF2var.cxx:6624-7281`).
//!
//! The first grid iso is the V frontier built by `InitGrid`
//! (`AdvApp2Var_ApproxAFunc2Var.cxx:271-279`). `FirstNotApprox` returns it
//! (`AdvApp2Var_Framework.cxx:77-106`). `mma2fnc_` runs `mma1nop_`, `mma1fdi_`,
//! `mma1cdi_`, `mma1jak_`, `mma1cnt_`, `mma1fer_`, `mmapcmp_`, and `mmjacan_`.
//! `ChangeIso` stores each frontier that meets tolerance. `Discretise` fills
//! the corner and iso tables. `mma2cdi_` runs `mma2cd3_`, `mma2cd2_`, and
//! `mma2cd1_`. `mma2ds1_` discretises the square. `NumDec` is chosen from
//! `PrefAndRec::Value`. `MakeApprox` cuts with `UpdateInU` and `UpdateInV`,
//! then stores a finished square. `Compute3DErrors` runs when every square
//! is approximated, and `ConvertBS` returns the B-spline.

use crate::bspline_surface::GeomBSplineSurface;
use crate::Surface;

use super::context::App2VarContext;
use super::eval_surface::evaluate;

pub(super) use super::constraints::mma1her;

/// `InitGrid(1)` then `ComputeConstraints`, then the corner tables of
/// `Patch::Discretise`. Cuts call `UpdateInU` and `UpdateInV`. The next
/// `ConvertBS` returns the B-spline.
pub fn approximate_until_cdi(
    surf: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ctx: &App2VarContext,
    max_patches: i32,
) -> Result<GeomBSplineSurface, String> {
    let mut frame = super::framework::Framework::init_grid(u0, u1, v0, v1, ctx.ord_u, ctx.ord_v);
    let mut net = super::network::Network::init(u0, u1, v0, v1, ctx.ord_u, ctx.ord_v);
    // `GeomConvert_ApproxSurface.cxx:381-398`: C2 breaks are recommended,
    // C3 breaks are preferred, and `Weight` stays 5.
    let u_choice = crate::adv_approx::PrefAndRec::new(surf.u_intervals(4), surf.u_intervals(5), 5.0)
        .map_err(|_| "PrefAndRec : Weight is too small".to_string())?;
    let v_choice = crate::adv_approx::PrefAndRec::new(surf.v_intervals(4), surf.v_intervals(5), 5.0)
        .map_err(|_| "PrefAndRec : Weight is too small".to_string())?;
    loop {
        frame.compute_constraints(
            &mut net,
            max_patches,
            &u_choice,
            &v_choice,
            |iso, begin, end| approx_iso(surf, ctx, iso, begin, end),
        )?;
        let Some(index) = net.first_not_approx() else {
            let _done = net.compute_3d_errors(ctx.nb_3d, &ctx.tol_3d, &ctx.tol_3d_front)?;
            return super::convert_bs::convert_bs(&mut net, ctx);
        };
        let (pu0, pu1, pv0, pv1) = net.domain(index)?;
        let (step, u_dec, v_dec) = super::discretise::discretise_corners(
            &frame,
            ctx,
            surf,
            pu0,
            pu1,
            pv0,
            pv1,
            max_patches,
            net.nb_patch_in_u(),
            net.nb_patch_in_v(),
        )?;
        if let Some(done) = step.finished {
            net.store_result(index, done)?;
            continue;
        }
        match step.sense {
            1 => {
                net.update_in_u(u_dec)?;
                frame.update_in_u(u_dec)?;
            }
            2 => {
                net.update_in_v(v_dec)?;
                frame.update_in_v(v_dec)?;
            }
            3 => {
                net.update_in_u(u_dec)?;
                frame.update_in_u(u_dec)?;
                net.update_in_v(v_dec)?;
                frame.update_in_v(v_dec)?;
            }
            _ => {
                return Err(
                    "AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error".into(),
                );
            }
        }
    }
}

fn approx_iso(
    surf: &dyn Surface,
    ctx: &App2VarContext,
    iso: &super::framework::GridIso,
    begin: &mut super::node::App2VarNode,
    end: &mut super::node::App2VarNode,
) -> Result<super::framework::IsoStatus, String> {
    let (extrem_order, deriv_order) = if iso.isofav == 2 {
        (ctx.ord_u, ctx.ord_v)
    } else {
        (ctx.ord_v, ctx.ord_u)
    };
    let mut errors = Vec::new();
    let mut somtab = Vec::new();
    let mut diftab = Vec::new();
    let mut polynom = Vec::new();
    let mut ncoeff = 0;
    let mut iercod = 0;
    for iapp in 0..=deriv_order {
        let cut = mma2fnc(
            surf,
            ctx,
            [iso.u0, iso.u1, iso.v0, iso.v1],
            iso.constant,
            iso.isofav,
            extrem_order,
            iapp,
            iso.position,
        )?;
        if cut.iercod > 0 {
            return Err("AdvApp2Var_ApproxAFunc2Var : Curve Approximation Error".into());
        }
        if cut.iercod == -1 {
            iercod = -1;
        }
        super::node::store_ends(
            begin,
            end,
            iapp,
            extrem_order,
            iso.isofav,
            &cut.contr1,
            &cut.contr2,
        )?;
        errors.push((cut.errmax, cut.errmoy));
        somtab.extend_from_slice(&cut.somtab);
        diftab.extend_from_slice(&cut.diftab);
        polynom.extend_from_slice(&cut.courbe);
        if ncoeff <= cut.ncoeff {
            ncoeff = cut.ncoeff;
        }
    }
    Ok(super::framework::IsoStatus {
        approximated: iercod == 0,
        has_result: iercod == 0 || iercod == -1,
        errors,
        somtab,
        diftab,
        polynom,
        nb_coeff: ncoeff,
    })
}

fn incoherent() -> String {
    "AdvApp2Var_ApproxF2var::mma2fnc_: AdvApp2Var_ApproxF2var.cxx:7251".into()
}

fn alloc_failed() -> String {
    "AdvApp2Var_ApproxF2var::mma2fnc_: AdvApp2Var_ApproxF2var.cxx:7247".into()
}

/// `mmveps3_` (`AdvApp2Var_MathBase.cxx:9171-9256`).
///
/// `mmprcsn_` is a zero-initialized static (`:221-225`). This path never
/// writes `eps3`.
fn eps3() -> f64 {
    0.0
}

fn eps_column(ctx: &App2VarContext, position: i32) -> Result<Vec<f64>, String> {
    let count = ctx.total_ssp();
    if count < 0 {
        return Err(incoherent());
    }
    let mut eps = vec![0.0; count as usize];
    for iesp in 0..count as usize {
        eps[iesp] = match position {
            0 => *ctx.internal_tol.get(iesp).ok_or_else(incoherent)?,
            1 => ctx.frontier_tol.get(iesp).ok_or_else(incoherent)?[0],
            2 => ctx.frontier_tol.get(iesp).ok_or_else(incoherent)?[1],
            3 => ctx.frontier_tol.get(iesp).ok_or_else(incoherent)?[2],
            4 => ctx.frontier_tol.get(iesp).ok_or_else(incoherent)?[3],
            // `MakeApprox` leaves the fresh array at 0 (`AdvApp2Var_Iso.cxx:182-216`).
            _ => 0.0,
        };
    }
    Ok(eps)
}

/// One derivative order of `mma2fnc_`.
struct BorderCut {
    iercod: i32,
    contr1: Vec<f64>,
    contr2: Vec<f64>,
    errmax: f64,
    errmoy: f64,
    somtab: Vec<f64>,
    diftab: Vec<f64>,
    courbe: Vec<f64>,
    ncoeff: i32,
}

/// One derivative order of `mma2fnc_` through `mma1noc_`.
fn mma2fnc(
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
        let hermit = super::constraints::mma1her(iordre)?;
        mma1cnt(ndimen, iordre, &contr1, &contr2, &hermit, ndgjac, &mut crvjac);
    }
    // `MakeApprox` passes one 3D subspace (`AdvApp2Var_Iso.cxx:121-122`).
    if ctx.total_ssp() != 1 || ndimen != 3 {
        return Err("UNPORTED AdvApp2Var_Iso::MakeApprox: AdvApp2Var_Iso.cxx:122".into());
    }
    let (ier, ncoeff, errmax, errmoy) = super::border_error::mma1fer(
        1,
        &[ndimen],
        iordre,
        ndgjac,
        &mut crvjac,
        ncflim,
        &epsapr,
    )?;
    let ncfja = ndgjac + 1;
    let packed = super::border_error::mmapcmp(ndimen, ncfja, ncoeff, &crvjac)?;
    let half = (ncoeff - 1) / 2 + 1;
    let mut canon = vec![0.0; (ndimen * ncoeff) as usize];
    for nd in 1..=ndimen {
        let src = (((nd - 1) << 1) * half) as usize;
        let dst = ((nd - 1) * ncoeff) as usize;
        super::jacobi_canon::mmjacan(
            iordre,
            ncoeff - 1,
            &packed[src..src + (half * 2) as usize],
            &mut canon[dst..dst + ncoeff as usize],
        )?;
    }
    let mut courbe = vec![0.0; (ncflim * ndimen) as usize];
    // `mma2fnc_` passes `ncoeff` as `NDIMEN` and the geometric dimension as
    // `NCOEFU` (`AdvApp2Var_ApproxF2var.cxx:7156-7164`).
    super::table_expand::mmfmca8(ncoeff, ndimen, 1, ncflim, ndimen, 1, &canon, &mut courbe)?;
    let uv11 = [-1.0, 1.0, -1.0, 1.0];
    if iordre >= 0 {
        for ii in 0..=iordre {
            let at = contr_at(ii, 1, ndimen);
            super::param_normalize::mma1noc(uv11, ndimen, ii, uvpav, isofav, ideriv, &mut contr1[at..])?;
            super::param_normalize::mma1noc(uv11, ndimen, ii, uvpav, isofav, ideriv, &mut contr2[at..])?;
        }
    }
    let span = (nbroot / 2 + 1) * ndimen;
    super::param_normalize::mma1noc(uv11, span, 0, uvpav, isofav, ideriv, &mut somtab)?;
    super::param_normalize::mma1noc(uv11, span, 0, uvpav, isofav, ideriv, &mut diftab)?;
    for nd in 1..=ndimen {
        let at = ((nd - 1) * ncflim) as usize;
        super::param_normalize::mma1noc(uv11, ncoeff, 0, uvpav, isofav, ideriv, &mut courbe[at..])?;
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

/// `mma1nop_` (`AdvApp2Var_ApproxF2var.cxx:1627-1734`).
fn mma1nop(
    nbroot: i32,
    rootlg: &[f64],
    uv: [f64; 4],
    isofav: i32,
    ttable: &mut [f64],
) -> Result<(), String> {
    let alinu = (uv[1] - uv[0]) / 2.0;
    let blinu = (uv[1] + uv[0]) / 2.0;
    let alinv = (uv[3] - uv[2]) / 2.0;
    let blinv = (uv[3] + uv[2]) / 2.0;
    if isofav == 1 {
        ttable[0] = uv[2];
        for ii in 1..=nbroot {
            ttable[ii as usize] = alinv * rootlg[(ii - 1) as usize] + blinv;
        }
        ttable[(nbroot + 1) as usize] = uv[3];
        Ok(())
    } else if isofav == 2 {
        ttable[0] = uv[0];
        for ii in 1..=nbroot {
            ttable[ii as usize] = alinu * rootlg[(ii - 1) as usize] + blinu;
        }
        ttable[(nbroot + 1) as usize] = uv[1];
        Ok(())
    } else {
        Err(incoherent())
    }
}

fn sample_at(point: i32, nd: i32, dim: i32) -> usize {
    ((nd - 1) + point * dim) as usize
}

fn pair_at(slot: i32, nd: i32, dim1: i32) -> usize {
    (slot + (nd - 1) * dim1) as usize
}

fn contr_at(order: i32, nd: i32, dim: i32) -> usize {
    (order * dim + (nd - 1)) as usize
}

/// `mma1fdi_` (`AdvApp2Var_ApproxF2var.cxx:596-981`).
fn mma1fdi(
    surf: &dyn Surface,
    ndimen: i32,
    uvfonc: [f64; 4],
    isofav: i32,
    tconst: f64,
    nbroot: i32,
    ttable: &[f64],
    iordre: i32,
    ideriv: i32,
    somtab: &mut [f64],
    diftab: &mut [f64],
    contr1: &mut [f64],
    contr2: &mut [f64],
) -> Result<(), String> {
    let (ideb, ifin) = if iordre >= 0 {
        (0, nbroot + 1)
    } else {
        (1, nbroot)
    };
    let nbp = ifin - ideb + 1;
    let nroo2 = nbroot / 2;
    let (ideru, iderv, renor) = if isofav == 1 {
        let half = (uvfonc[1] - uvfonc[0]) / 2.0;
        (ideriv, 0, super::param_normalize::pow_di(half, ideriv))
    } else {
        let half = (uvfonc[3] - uvfonc[2]) / 2.0;
        (0, ideriv, super::param_normalize::pow_di(half, ideriv))
    };

    let mut fpntab = vec![0.0; (nbp * ndimen) as usize];
    let code = evaluate(
        surf,
        ndimen,
        [uvfonc[0], uvfonc[1]],
        [uvfonc[2], uvfonc[3]],
        isofav,
        tconst,
        &ttable[ideb as usize..(ideb + nbp) as usize],
        ideru,
        iderv,
        &mut fpntab,
    )?;
    if code > 0 {
        return Err(format!(
            "AdvApp2Var_ApproxF2var::mma1fdi_: AdvApp2Var_ApproxF2var.cxx:973 ({})",
            code + 100
        ));
    }

    let dim1 = nroo2 + 1;
    // Samples start at `ideb`. With constraints `ideb` is 0, so point `p`
    // in this buffer is fpntab point `p` in the FORTRAN table.
    let point = |p: i32, nd: i32| sample_at(p - ideb, nd, ndimen);
    for nd in 1..=ndimen {
        for ii in 1..=nroo2 {
            let iip = (nbroot + 1) / 2 + ii;
            let iim = nroo2 - ii + 1;
            let bid1 = fpntab[point(iim, nd)];
            let bid2 = fpntab[point(iip, nd)];
            somtab[pair_at(ii, nd, dim1)] = renor * (bid2 + bid1);
            diftab[pair_at(ii, nd, dim1)] = renor * (bid2 - bid1);
        }
    }
    if nbroot % 2 == 1 {
        for nd in 1..=ndimen {
            let mid = fpntab[point(nroo2 + 1, nd)];
            somtab[pair_at(0, nd, dim1)] = renor * mid;
            diftab[pair_at(0, nd, dim1)] = renor * mid;
        }
    } else {
        for nd in 1..=ndimen {
            somtab[pair_at(0, nd, dim1)] = 0.0;
            diftab[pair_at(0, nd, dim1)] = 0.0;
        }
    }

    if iordre >= 0 {
        for nd in 1..=ndimen {
            contr1[contr_at(0, nd, ndimen)] = renor * fpntab[point(0, nd)];
            contr2[contr_at(0, nd, ndimen)] = renor * fpntab[point(nbroot + 1, nd)];
        }
        let along = if isofav == 1 {
            (uvfonc[3] - uvfonc[2]) / 2.0
        } else {
            (uvfonc[1] - uvfonc[0]) / 2.0
        };
        let end = (nbroot + 1) as usize;
        for order in 1..=iordre {
            // Iso-U derives along V. Iso-V derives along U (`:852-946`).
            let (du, dv) = if isofav == 1 {
                (ideru, order)
            } else {
                (order, iderv)
            };
            fill_constraint(
                surf,
                ndimen,
                uvfonc,
                isofav,
                tconst,
                &[ttable[0]],
                du,
                dv,
                contr1,
                order,
            )?;
            fill_constraint(
                surf,
                ndimen,
                uvfonc,
                isofav,
                tconst,
                &[ttable[end]],
                du,
                dv,
                contr2,
                order,
            )?;
        }
        let mut bid2 = renor;
        for ii in 1..=iordre {
            bid2 *= along;
            for nd in 1..=ndimen {
                contr1[contr_at(ii, nd, ndimen)] *= bid2;
                contr2[contr_at(ii, nd, ndimen)] *= bid2;
            }
        }
    }
    Ok(())
}

fn fill_constraint(
    surf: &dyn Surface,
    ndimen: i32,
    uvfonc: [f64; 4],
    isofav: i32,
    tconst: f64,
    parameters: &[f64],
    ideru: i32,
    iderv: i32,
    contr: &mut [f64],
    order: i32,
) -> Result<(), String> {
    let mut sample = vec![0.0; ndimen as usize];
    let code = evaluate(
        surf,
        ndimen,
        [uvfonc[0], uvfonc[1]],
        [uvfonc[2], uvfonc[3]],
        isofav,
        tconst,
        parameters,
        ideru,
        iderv,
        &mut sample,
    )?;
    if code > 0 {
        return Err(format!(
            "AdvApp2Var_ApproxF2var::mma1fdi_: AdvApp2Var_ApproxF2var.cxx:973 ({})",
            code + 100
        ));
    }
    for nd in 1..=ndimen {
        contr[contr_at(order, nd, ndimen)] = sample[(nd - 1) as usize];
    }
    Ok(())
}

/// `GeomConvert_ApproxSurface_Eval::Evaluate`
/// (`GeomConvert_ApproxSurface.cxx:55-310`).
///
/// Orders 0, 1, and 2 use `D0` / `D1` / `D2`. `D3` and `DN` are not on the
/// C1 iso and stay unported.
/// `mma1cdi_` (`AdvApp2Var_ApproxF2var.cxx:233-465`).
fn mma1cdi(
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

fn hermite_sample(fpntab: &[f64], nbroot: i32, root: i32, poly: i32) -> f64 {
    fpntab[((poly - 1) * nbroot + (root - 1)) as usize]
}

/// `mma1cnt_` (`AdvApp2Var_ApproxF2var.cxx:469-592`).
fn mma1cnt(
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
fn mma1jak(
    ndimen: i32,
    nbroot: i32,
    iordre: i32,
    ndgjac: i32,
    somtab: &[f64],
    diftab: &[f64],
) -> Result<Vec<f64>, String> {
    let columns = ndgjac - 2 * iordre - 1;
    let cgauss = super::gauss::gauss_weights(ndgjac, nbroot, iordre, columns)?;
    Ok(mmmapcoe(
        ndimen, ndgjac, iordre, nbroot, somtab, diftab, &cgauss,
    ))
}

/// `mmmapcoe_` (`AdvApp2Var_ApproxF2var.cxx:8028-8194`).
fn mmmapcoe(
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
pub(super) fn mmmpocur(ndeg: i32, courbe: &[f64], tparam: f64) -> f64 {
    if ndeg < 1 {
        return 0.0;
    }
    let mut fu = courbe[(ndeg - 1) as usize];
    for i in (1..ndeg).rev() {
        fu = fu * tparam + courbe[(i - 1) as usize];
    }
    fu
}
