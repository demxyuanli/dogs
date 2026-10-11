//! `AdvApp2Var_Patch::AddConstraints` (`AdvApp2Var_Patch.cxx:388-651`).

use super::framework::Framework;

fn fortran() -> String {
    "AdvApp2Var_Patch::AddConstraints : Error in FORTRAN".into()
}

fn copy_oob() -> String {
    "AdvApp2Var_Patch::AddConstraints: AdvApp2Var_Patch.cxx:409".into()
}

fn ac1_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2ac1_: AdvApp2Var_ApproxF2var.cxx:1882".into()
}

fn ac2_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2ac2_: AdvApp2Var_ApproxF2var.cxx:2023".into()
}

fn ac3_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2ac3_: AdvApp2Var_ApproxF2var.cxx:2168".into()
}

fn drc_oob() -> String {
    "AdvApp2Var_MathBase::mmdrc11_: AdvApp2Var_MathBase.cxx:2968".into()
}

fn coin_oob() -> String {
    "AdvApp2Var_Patch::AddConstraints: AdvApp2Var_Patch.cxx:533".into()
}

/// `mma1her_` (`AdvApp2Var_ApproxF2var.cxx:1204-1364`).
pub(super) fn mma1her(iordre: i32) -> Result<Vec<f64>, String> {
    let dim1 = (iordre + 1) << 1;
    if dim1 < 0 {
        return Err("AdvApp2Var_ApproxF2var::mma2fnc_: AdvApp2Var_ApproxF2var.cxx:7251".into());
    }
    let mut hermit = vec![0.0; (dim1 * dim1) as usize];
    let put = |table: &mut [f64], poly: i32, coeff: i32, value: f64| {
        table[((poly - 1) * dim1 + (coeff - 1)) as usize] = value;
    };
    if iordre == 0 {
        put(&mut hermit, 1, 1, 0.5);
        put(&mut hermit, 1, 2, -0.5);
        put(&mut hermit, 2, 1, 0.5);
        put(&mut hermit, 2, 2, 0.5);
    } else if iordre == 1 {
        let rows = [
            [0.5, -0.75, 0.0, 0.25],
            [0.5, 0.75, 0.0, -0.25],
            [0.25, -0.25, -0.25, 0.25],
            [-0.25, -0.25, 0.25, 0.25],
        ];
        for (poly, row) in rows.iter().enumerate() {
            for (coeff, value) in row.iter().enumerate() {
                put(&mut hermit, (poly as i32) + 1, (coeff as i32) + 1, *value);
            }
        }
    } else if iordre == 2 {
        let rows = [
            [0.5, -0.9375, 0.0, 0.625, 0.0, -0.1875],
            [0.5, 0.9375, 0.0, -0.625, 0.0, 0.1875],
            [0.3125, -0.4375, -0.375, 0.625, 0.0625, -0.1875],
            [-0.3125, -0.4375, 0.375, 0.625, -0.0625, -0.1875],
            [0.0625, -0.0625, -0.125, 0.125, 0.0625, -0.0625],
            [0.0625, 0.0625, -0.125, -0.125, 0.0625, 0.0625],
        ];
        for (poly, row) in rows.iter().enumerate() {
            for (coeff, value) in row.iter().enumerate() {
                put(&mut hermit, (poly as i32) + 1, (coeff as i32) + 1, *value);
            }
        }
    } else {
        return Err("AdvApp2Var_ApproxF2var::mma2fnc_: AdvApp2Var_ApproxF2var.cxx:7251".into());
    }
    Ok(hermit)
}

fn herm(table: &[f64], dim1: i32, coeff: i32, poly: i32, oob: fn() -> String) -> Result<f64, String> {
    let at = ((poly - 1) * dim1 + coeff) as usize;
    table.get(at).copied().ok_or_else(oob)
}

fn crb(
    table: &[f64],
    kk: i32,
    nd: i32,
    ii: i32,
    nclim: i32,
    ndimen: i32,
    oob: fn() -> String,
) -> Result<f64, String> {
    let at = (kk + (nd - 1) * nclim + (ii - 1) * ndimen * nclim) as usize;
    table.get(at).copied().ok_or_else(oob)
}

fn add_pat(
    pat: &mut [f64],
    ku: i32,
    kv: i32,
    nd: i32,
    dim1: i32,
    dim2: i32,
    delta: f64,
    oob: fn() -> String,
) -> Result<(), String> {
    let at = (ku + kv * dim1 + (nd - 1) * dim1 * dim2) as usize;
    let slot = pat.get_mut(at).ok_or_else(oob)?;
    *slot += delta;
    Ok(())
}

/// `mma2ac2_` (`AdvApp2Var_ApproxF2var.cxx:1908-2050`).
fn mma2ac2(
    ndimen: i32,
    mxujac: i32,
    mxvjac: i32,
    iordrv: i32,
    nclimu: i32,
    ncf0: i32,
    crb0: &[f64],
    ncf1: i32,
    crb1: &[f64],
    vhermt: &[f64],
    patjac: &mut [f64],
) -> Result<(), String> {
    let dim1 = mxujac + 1;
    let dim2 = mxvjac + 1;
    let hdim = (iordrv << 1) + 2;
    if iordrv < 0 || ndimen < 1 {
        return Ok(());
    }
    for ii in 1..=(iordrv + 1) {
        let ndg0 = ncf0 - 1;
        let ndg1 = ncf1 - 1;
        for nd in 1..=ndimen {
            for jj in 0..=((iordrv << 1) + 1) {
                let bid1 = herm(vhermt, hdim, jj, (ii << 1) - 1, ac2_oob)?;
                if ndg0 >= 0 {
                    for kk in 0..=ndg0 {
                        let coef = crb(crb0, kk, nd, ii, nclimu, ndimen, ac2_oob)?;
                        add_pat(patjac, kk, jj, nd, dim1, dim2, bid1 * coef, ac2_oob)?;
                    }
                }
                let bid2 = herm(vhermt, hdim, jj, ii << 1, ac2_oob)?;
                if ndg1 >= 0 {
                    for kk in 0..=ndg1 {
                        let coef = crb(crb1, kk, nd, ii, nclimu, ndimen, ac2_oob)?;
                        add_pat(patjac, kk, jj, nd, dim1, dim2, bid2 * coef, ac2_oob)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// `mma2ac3_` (`AdvApp2Var_ApproxF2var.cxx:2054-2201`).
fn mma2ac3(
    ndimen: i32,
    mxujac: i32,
    mxvjac: i32,
    iordru: i32,
    nclimv: i32,
    ncf0: i32,
    crb0: &[f64],
    ncf1: i32,
    crb1: &[f64],
    uhermt: &[f64],
    patjac: &mut [f64],
) -> Result<(), String> {
    let dim1 = mxujac + 1;
    let dim2 = mxvjac + 1;
    let hdim = (iordru << 1) + 2;
    let ndgu = (iordru << 1) + 1;
    if iordru < 0 || ndimen < 1 {
        return Ok(());
    }
    for ii in 1..=(iordru + 1) {
        let ndg0 = ncf0 - 1;
        let ndg1 = ncf1 - 1;
        for nd in 1..=ndimen {
            if ndg0 >= 0 {
                for jj in 0..=ndg0 {
                    let bid1 = crb(crb0, jj, nd, ii, nclimv, ndimen, ac3_oob)?;
                    for kk in 0..=ndgu {
                        let weight = herm(uhermt, hdim, kk, (ii << 1) - 1, ac3_oob)?;
                        add_pat(patjac, kk, jj, nd, dim1, dim2, bid1 * weight, ac3_oob)?;
                    }
                }
            }
            if ndg1 >= 0 {
                for jj in 0..=ndg1 {
                    let bid2 = crb(crb1, jj, nd, ii, nclimv, ndimen, ac3_oob)?;
                    for kk in 0..=ndgu {
                        let weight = herm(uhermt, hdim, kk, ii << 1, ac3_oob)?;
                        add_pat(patjac, kk, jj, nd, dim1, dim2, bid2 * weight, ac3_oob)?;
                    }
                }
            }
        }
    }
    Ok(())
}

fn contr(coins: &[f64], corner: i32, nd: i32, ii: i32, jj: i32, ndimen: i32, dim2: i32, size: i32) -> Result<f64, String> {
    let at = ((corner * size) + (nd - 1) + (ii - 1) * ndimen + (jj - 1) * ndimen * dim2) as usize;
    coins.get(at).copied().ok_or_else(ac1_oob)
}

/// `mma2ac1_` (`AdvApp2Var_ApproxF2var.cxx:1738-1904`).
fn mma2ac1(
    ndimen: i32,
    mxujac: i32,
    mxvjac: i32,
    iordru: i32,
    iordrv: i32,
    coins: &[f64],
    uhermt: &[f64],
    vhermt: &[f64],
    patjac: &mut [f64],
) -> Result<(), String> {
    let dim1 = mxujac + 1;
    let dim2 = mxvjac + 1;
    let hudim = (iordru << 1) + 2;
    let hvdim = (iordrv << 1) + 2;
    let cdim = iordru + 2;
    let size = ndimen * cdim * (iordrv + 2);
    let ndgu = (iordru << 1) + 1;
    let ndgv = (iordrv << 1) + 1;
    if iordru < 0 || iordrv < 0 || ndimen < 1 {
        return Ok(());
    }
    for jj in 1..=(iordrv + 1) {
        for ii in 1..=(iordru + 1) {
            for nd in 1..=ndimen {
                let cnt1 = contr(coins, 0, nd, ii, jj, ndimen, cdim, size)?;
                let cnt2 = contr(coins, 1, nd, ii, jj, ndimen, cdim, size)?;
                let cnt3 = contr(coins, 2, nd, ii, jj, ndimen, cdim, size)?;
                let cnt4 = contr(coins, 3, nd, ii, jj, ndimen, cdim, size)?;
                for kv in 0..=ndgv {
                    let bidv1 = herm(vhermt, hvdim, kv, (jj << 1) - 1, ac1_oob)?;
                    let bidv2 = herm(vhermt, hvdim, kv, jj << 1, ac1_oob)?;
                    for ku in 0..=ndgu {
                        let bidu1 = herm(uhermt, hudim, ku, (ii << 1) - 1, ac1_oob)?;
                        let bidu2 = herm(uhermt, hudim, ku, ii << 1, ac1_oob)?;
                        let delta = -bidu1 * bidv1 * cnt1
                            - bidu2 * bidv1 * cnt2
                            - bidu1 * bidv2 * cnt3
                            - bidu2 * bidv2 * cnt4;
                        add_pat(patjac, ku, kv, nd, dim1, dim2, delta, ac1_oob)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// `mfac_` (`AdvApp2Var_MathBase.cxx:347-370`). `f[1] = 1`, then `f[i] = i * f[i-1]`.
fn mfac(fact: &mut [f64], n: i32) -> Result<(), String> {
    if n < 1 || fact.len() <= n as usize {
        return Err(drc_oob());
    }
    fact[1] = 1.0;
    if n >= 2 {
        for i in 2..=n {
            let prev = fact[(i - 1) as usize];
            fact[i as usize] = (i as f64) * prev;
        }
    }
    Ok(())
}

/// `mmdrc11_` (`AdvApp2Var_MathBase.cxx:2831-3013`).
fn mmdrc11(
    iordre: i32,
    ndimen: i32,
    ncoeff: i32,
    courbe: &[f64],
    points: &mut [f64],
    fact: &mut [f64],
) -> Result<(), String> {
    if iordre < 0 || ncoeff < 1 {
        return Ok(());
    }
    if ndimen < 1 {
        return Ok(());
    }
    let span = (iordre + 1) * 2;
    for nd in 1..=ndimen {
        let base = ((nd - 1) * span) as usize;
        let last = ((nd - 1) * ncoeff + ncoeff - 1) as usize;
        let lead = *courbe.get(last).ok_or_else(drc_oob)?;
        let end = base + span as usize;
        let row = points.get_mut(base..end).ok_or_else(drc_oob)?;
        row[0] = lead;
        row[1] = lead;
        if iordre >= 1 {
            for j in 1..=iordre {
                row[(2 * j) as usize] = 0.0;
                row[(2 * j + 1) as usize] = 0.0;
            }
        }
    }
    for nd in 1..=ndimen {
        let base = ((nd - 1) * span) as usize;
        let row0 = ((nd - 1) * ncoeff) as usize;
        for ndeg in 1..ncoeff {
            let coeff = *courbe.get(row0 + (ncoeff - 1 - ndeg) as usize).ok_or_else(drc_oob)?;
            if iordre >= 1 {
                for i in (1..=iordre).rev() {
                    let hi = base + (2 * i) as usize;
                    let lo = base + (2 * (i - 1)) as usize;
                    let left = points[lo];
                    let right = points[lo + 1];
                    let cur = points[hi];
                    points[hi] = -cur + left;
                    points[hi + 1] += right;
                }
            }
            let cur = points[base];
            points[base] = -cur + coeff;
            points[base + 1] += coeff;
        }
    }
    if iordre > 1 {
        mfac(fact, iordre)?;
        for nd in 1..=ndimen {
            let base = ((nd - 1) * span) as usize;
            for i in 2..=iordre {
                let weight = fact[i as usize];
                points[base + (2 * i) as usize] *= weight;
                points[base + (2 * i + 1) as usize] *= weight;
            }
        }
    }
    Ok(())
}

fn take_iso(poly: &[f64], blocks: i32, size: i32) -> Result<Vec<f64>, String> {
    if blocks <= 0 || size <= 0 {
        return Ok(Vec::new());
    }
    let needed = (blocks as usize)
        .checked_mul(size as usize)
        .ok_or_else(copy_oob)?;
    if poly.len() < needed {
        return Err(copy_oob());
    }
    Ok(poly[..needed].to_vec())
}

fn scale_block(iso: &mut [f64], order: i32, size: i32, half: f64) -> Result<(), String> {
    if order < 1 || size <= 0 {
        return Ok(());
    }
    for iu in 1..=order {
        let rho = half.powi(iu);
        let base = (iu * size) as usize;
        let end = base + size as usize;
        let block = iso.get_mut(base..end).ok_or_else(copy_oob)?;
        for value in block {
            *value *= rho;
        }
    }
    Ok(())
}

fn coin_at(corner: i32, idim: i32, iu: i32, iv: i32, ndimen: i32, span_u: i32, size: i32) -> Result<usize, String> {
    let at = (corner * size) + (idim - 1) + ndimen * iu + ndimen * span_u * iv;
    usize::try_from(at).map_err(|_| coin_oob())
}

fn add_coin(coins: &mut [f64], at: usize, delta: f64) -> Result<(), String> {
    let slot = coins.get_mut(at).ok_or_else(coin_oob)?;
    *slot += delta;
    Ok(())
}

/// `AddConstraints` (`AdvApp2Var_Patch.cxx:388-651`).
pub(super) fn add_constraints(
    frame: &Framework,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ndimen: i32,
    iordru: i32,
    iordrv: i32,
    ncfmxu: i32,
    ncfmxv: i32,
    patcan: &mut [f64],
) -> Result<(), String> {
    let deg_u = ncfmxu - 1;
    let deg_v = ncfmxv - 1;
    let size_u = ncfmxv * ndimen;
    let size_v = ncfmxu * ndimen;
    let mut iso_u0 = take_iso(&frame.iso_u(u0, v0, v1)?.polynom, iordru + 1, size_u)?;
    let mut iso_u1 = take_iso(&frame.iso_u(u1, v0, v1)?.polynom, iordru + 1, size_u)?;
    let ncf_u0 = frame.iso_u(u0, v0, v1)?.nb_coeff;
    let ncf_u1 = frame.iso_u(u1, v0, v1)?.nb_coeff;
    let du = (u1 - u0) / 2.0;
    let dv = (v1 - v0) / 2.0;
    scale_block(&mut iso_u0, iordru, size_u, du)?;
    scale_block(&mut iso_u1, iordru, size_u, du)?;

    let mut iso_v0 = take_iso(&frame.iso_v(u0, u1, v0)?.polynom, iordrv + 1, size_v)?;
    let mut iso_v1 = take_iso(&frame.iso_v(u0, u1, v1)?.polynom, iordrv + 1, size_v)?;
    let ncf_v0 = frame.iso_v(u0, u1, v0)?.nb_coeff;
    let ncf_v1 = frame.iso_v(u0, u1, v1)?.nb_coeff;
    scale_block(&mut iso_v0, iordrv, size_v, dv)?;
    scale_block(&mut iso_v1, iordrv, size_v, dv)?;

    let herm_v = if iordrv >= 0 {
        let herm = mma1her(iordrv).map_err(|_| fortran())?;
        mma2ac2(
            ndimen, deg_u, deg_v, iordrv, ncfmxu, ncf_v0, &iso_v0, ncf_v1, &iso_v1, &herm, patcan,
        )?;
        Some(herm)
    } else {
        None
    };
    let herm_u = if iordru >= 0 {
        let herm = mma1her(iordru).map_err(|_| fortran())?;
        mma2ac3(
            ndimen, deg_u, deg_v, iordru, ncfmxv, ncf_u0, &iso_u0, ncf_u1, &iso_u1, &herm, patcan,
        )?;
        Some(herm)
    } else {
        None
    };

    let span_u = iordru + 2;
    let span_v = iordrv + 2;
    let size = ndimen
        .checked_mul(span_u)
        .and_then(|n| n.checked_mul(span_v))
        .ok_or_else(coin_oob)?;
    let mut coins = vec![0.0; (size * 4).max(0) as usize];
    if iordru >= 0 && iordrv >= 0 {
        for iu in 0..=iordru {
            for iv in 0..=iordrv {
                let rho = du.powi(iu) * dv.powi(iv);
                let corners = [
                    frame.node_at(u0, v0)?.point(iu, iv)?,
                    frame.node_at(u1, v0)?.point(iu, iv)?,
                    frame.node_at(u0, v1)?.point(iu, iv)?,
                    frame.node_at(u1, v1)?.point(iu, iv)?,
                ];
                for (corner, point) in corners.iter().enumerate() {
                    for idim in 1..=ndimen {
                        let at = coin_at(corner as i32, idim, iu, iv, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, -rho * point[(idim - 1) as usize])?;
                    }
                }
            }
        }
    }

    let iord_max = iordru.max(iordrv);
    let mut extr = vec![0.0; (2 * iord_max + 2).max(0) as usize];
    let mut fact = vec![0.0; (iord_max + 1).max(0) as usize];
    if iordru >= 0 {
        for iu in 1..=(iordru + 1) {
            for idim in 1..=ndimen {
                let row = ((idim - 1) * ncfmxv + ncfmxv * ndimen * (iu - 1)) as usize;
                let ncf = ncf_u0;
                if ncf >= 1 {
                    let end = row + ncf as usize;
                    let curve = iso_u0.get(row..end).ok_or_else(drc_oob)?;
                    mmdrc11(iordrv, 1, ncf, curve, &mut extr, &mut fact)?;
                }
                if iordrv >= 0 {
                    for iv in 1..=(iordrv + 1) {
                        let at = coin_at(0, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iv - 1)) as usize])?;
                        let at = coin_at(2, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iv - 1) + 1) as usize])?;
                    }
                }
                let ncf = ncf_u1;
                if ncf >= 1 {
                    let end = row + ncf as usize;
                    let curve = iso_u1.get(row..end).ok_or_else(drc_oob)?;
                    mmdrc11(iordrv, 1, ncf, curve, &mut extr, &mut fact)?;
                }
                if iordrv >= 0 {
                    for iv in 1..=(iordrv + 1) {
                        let at = coin_at(1, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iv - 1)) as usize])?;
                        let at = coin_at(3, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iv - 1) + 1) as usize])?;
                    }
                }
            }
        }
    }
    if iordrv >= 0 {
        for iv in 1..=(iordrv + 1) {
            for idim in 1..=ndimen {
                let row = ((idim - 1) * ncfmxu + ncfmxu * ndimen * (iv - 1)) as usize;
                if ncf_v0 >= 1 {
                    let end = row + ncf_v0 as usize;
                    let curve = iso_v0.get(row..end).ok_or_else(drc_oob)?;
                    mmdrc11(iordru, 1, ncf_v0, curve, &mut extr, &mut fact)?;
                }
                if iordru >= 0 {
                    for iu in 1..=(iordru + 1) {
                        let at = coin_at(0, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iu - 1)) as usize])?;
                        let at = coin_at(1, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iu - 1) + 1) as usize])?;
                    }
                }
                if ncf_v1 >= 1 {
                    let end = row + ncf_v1 as usize;
                    let curve = iso_v1.get(row..end).ok_or_else(drc_oob)?;
                    mmdrc11(iordru, 1, ncf_v1, curve, &mut extr, &mut fact)?;
                }
                if iordru >= 0 {
                    for iu in 1..=(iordru + 1) {
                        let at = coin_at(2, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iu - 1)) as usize])?;
                        let at = coin_at(3, idim, iu - 1, iv - 1, ndimen, span_u, size)?;
                        add_coin(&mut coins, at, extr[(2 * (iu - 1) + 1) as usize])?;
                    }
                }
            }
        }
    }
    if iordru >= 0 && iordrv >= 0 {
        let uherm = herm_u.ok_or_else(fortran)?;
        let vherm = herm_v.ok_or_else(fortran)?;
        mma2ac1(
            ndimen, deg_u, deg_v, iordru, iordrv, &coins, &uherm, &vherm, patcan,
        )?;
    }
    Ok(())
}

fn err_oob() -> String {
    "AdvApp2Var_Patch::AddErrors: AdvApp2Var_Patch.cxx:672".into()
}

fn hmax(order: i32) -> Result<f64, String> {
    const TABLE: [f64; 4] = [0.0, 1.0, 1.5, 1.75];
    let at = order + 1;
    TABLE
        .get(usize::try_from(at).map_err(|_| err_oob())?)
        .copied()
        .ok_or_else(err_oob)
}

fn peak(iso: &super::framework::GridIso, order: i32, moy: bool) -> Result<f64, String> {
    let mut err: f64 = 0.0;
    if order >= 0 {
        for iv in 1..=(order + 1) {
            let pair = iso.errors.get((iv - 1) as usize).ok_or_else(err_oob)?;
            let error = if moy { pair.1 } else { pair.0 };
            err = err.max(error);
        }
    }
    Ok(err)
}

fn column1(iso: &super::framework::GridIso) -> Result<f64, String> {
    iso.errors.first().map(|pair| pair.0).ok_or_else(err_oob)
}

/// `AddErrors` (`AdvApp2Var_Patch.cxx:655-760`). `NBSESP` is 1.
pub(super) fn add_errors(
    frame: &Framework,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    iordru: i32,
    iordrv: i32,
    errmax: f64,
    errmoy: f64,
) -> Result<(f64, f64, [f64; 4]), String> {
    let iso_v0 = frame.iso_v(u0, u1, v0)?;
    let iso_v1 = frame.iso_v(u0, u1, v1)?;
    let iso_u0 = frame.iso_u(u0, v0, v1)?;
    let iso_u1 = frame.iso_u(u1, v0, v1)?;
    let hu = hmax(iordru)?;
    let hv = hmax(iordrv)?;
    let err_u = peak(iso_v0, iordrv, false)?.max(peak(iso_v1, iordrv, false)?);
    let err_v = peak(iso_u0, iordru, false)?.max(peak(iso_u1, iordru, false)?);
    let max_out = errmax + err_u * hv + err_v * hu;
    let moy_u = peak(iso_v0, iordrv, true)?.max(peak(iso_v1, iordrv, true)?);
    let moy_v = peak(iso_u0, iordru, true)?.max(peak(iso_u1, iordru, true)?);
    let moy_out = (errmoy * errmoy + moy_u * hv * moy_u * hv + moy_v * hu * moy_v * hu).sqrt();
    let mut heriso = [
        column1(iso_v0)?,
        column1(iso_v1)?,
        column1(iso_u0)?,
        column1(iso_u1)?,
    ];
    let mut emax: [f64; 4] = [0.0; 4];
    if iordru >= 0 && iordrv >= 0 {
        for iu in 0..=iordru {
            for iv in 0..=iordrv {
                let corners = [
                    frame.node_at(u0, v0)?.error(iu, iv)?,
                    frame.node_at(u1, v0)?.error(iu, iv)?,
                    frame.node_at(u0, v1)?.error(iu, iv)?,
                    frame.node_at(u1, v1)?.error(iu, iv)?,
                ];
                for (slot, error) in emax.iter_mut().zip(corners) {
                    *slot = (*slot).max(error);
                }
            }
        }
    }
    let err1 = emax[0].max(emax[1]);
    let err2 = emax[2].max(emax[3]);
    let err3 = emax[0].max(emax[2]);
    let err4 = emax[1].max(emax[3]);
    if iso_v0.position == 0 {
        heriso[0] += err1 * hu;
    }
    if iso_v1.position == 0 {
        heriso[1] += err2 * hu;
    }
    if iso_u0.position == 0 {
        heriso[2] += err3 * hv;
    }
    if iso_u1.position == 0 {
        heriso[3] += err4 * hv;
    }
    Ok((max_out, moy_out, heriso))
}
