//! `mma2ds1_` and `mma2ds2_` (`AdvApp2Var_ApproxF2var.cxx:5383-6205`).
//!
//! `iuouv == 2` discretises directly. `iuouv == 1` transposes with
//! `mmfmtb1_` (`AdvApp2Var_ApproxF2var.cxx:5634`), calls `mma2ds2_` with
//! U and V swapped, then transposes back (`:5707`).

use crate::Surface;

use super::eval_surface::evaluate;

/// `mma2ds1_` (`AdvApp2Var_ApproxF2var.cxx:5383-5785`).
pub(super) fn mma2ds1(
    surf: &dyn Surface,
    isofav: i32,
    ndimen: i32,
    u_bounds: [f64; 2],
    v_bounds: [f64; 2],
    nbpnt_u: i32,
    uroot: &[f64],
    nbpnt_v: i32,
    vroot: &[f64],
    sosotb: &mut [f64],
    disotb: &mut [f64],
    soditb: &mut [f64],
    diditb: &mut [f64],
) -> Result<(), String> {
    let iuouv = if (1..=2).contains(&isofav) { isofav } else { 2 };
    if iuouv == 2 {
        return mma2ds2(
            surf, ndimen, u_bounds, v_bounds, nbpnt_u, uroot, nbpnt_v, vroot, iuouv, sosotb,
            disotb, soditb, diditb,
        );
    }
    // `iuouv == 1` (`AdvApp2Var_ApproxF2var.cxx:5626-5768`).
    transpose_pairs(ndimen, nbpnt_u, nbpnt_v, sosotb, diditb, soditb, disotb)?;
    mma2ds2(
        surf, ndimen, v_bounds, u_bounds, nbpnt_v, vroot, nbpnt_u, uroot, iuouv, sosotb, soditb,
        disotb, diditb,
    )?;
    transpose_pairs(ndimen, nbpnt_v, nbpnt_u, sosotb, diditb, soditb, disotb)
}

fn transpose_pairs(
    ndimen: i32,
    nbpnt_u: i32,
    nbpnt_v: i32,
    sosotb: &mut [f64],
    diditb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
) -> Result<(), String> {
    let s1 = nbpnt_u / 2 + 1;
    let s2 = nbpnt_v / 2 + 1;
    let d1 = nbpnt_u / 2;
    let d2 = nbpnt_v / 2;
    for nd in 1..=ndimen {
        let soso_base = ((nd - 1) * s1 * s2) as usize;
        mmfmtb1(sosotb, soso_base, s1, s1, s2, s2)?;
        mmfmtb1(diditb, soso_base, s1, s1, s2, s2)?;
        let diso_base = ((nd - 1) * d1 * d2) as usize;
        mmfmtb1(soditb, diso_base, d1, d1, d2, d2)?;
        mmfmtb1(disotb, diso_base, d1, d1, d2, d2)?;
    }
    Ok(())
}

/// `AdvApp2Var_MathBase::mmfmtb1_` (`AdvApp2Var_MathBase.cxx:4412-4550`).
/// In-place transpose of one block: `T(i,j)` becomes `T(j,i)`.
fn mmfmtb1(
    table: &mut [f64],
    base: usize,
    maxsz1: i32,
    isize1: i32,
    jsize1: i32,
    maxsz2: i32,
) -> Result<(), String> {
    if isize1 > maxsz1 || jsize1 > maxsz2 {
        return Err("AdvApp2Var_MathBase::mmfmtb1_: AdvApp2Var_MathBase.cxx:4490".into());
    }
    let work_len = (maxsz2 as i64)
        .checked_mul(isize1 as i64)
        .filter(|n| *n >= 0)
        .and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| "AdvApp2Var_MathBase::mmfmtb1_: AdvApp2Var_MathBase.cxx:4498".to_string())?;
    let mut work = vec![0.0; work_len];
    for ii in 1..=isize1 {
        for jj in 1..=jsize1 {
            let src = base + ((ii - 1) + (jj - 1) * maxsz1) as usize;
            let dst = ((ii - 1) * maxsz2 + (jj - 1)) as usize;
            let value = *table.get(src).ok_or_else(|| {
                "AdvApp2Var_MathBase::mmfmtb1_: AdvApp2Var_MathBase.cxx:4513".to_string()
            })?;
            let slot = work.get_mut(dst).ok_or_else(|| {
                "AdvApp2Var_MathBase::mmfmtb1_: AdvApp2Var_MathBase.cxx:4513".to_string()
            })?;
            *slot = value;
        }
    }
    let dest = table.get_mut(base..base + work_len).ok_or_else(|| {
        "AdvApp2Var_MathBase::mmfmtb1_: AdvApp2Var_MathBase.cxx:4519".to_string()
    })?;
    dest.copy_from_slice(&work);
    Ok(())
}

fn ds2_oob() -> String {
    "AdvApp2Var_ApproxF2var::mma2ds2_: AdvApp2Var_ApproxF2var.cxx:6073".into()
}

fn take_at(buf: &[f64], index: i32) -> Result<f64, String> {
    buf.get(index as usize).copied().ok_or_else(ds2_oob)
}

fn add_at(buf: &mut [f64], index: i32, delta: f64) -> Result<(), String> {
    let slot = buf.get_mut(index as usize).ok_or_else(ds2_oob)?;
    *slot += delta;
    Ok(())
}

/// `mma2ds2_` (`AdvApp2Var_ApproxF2var.cxx:5789-6205`).
fn mma2ds2(
    surf: &dyn Surface,
    ndimen: i32,
    uintfn: [f64; 2],
    vintfn: [f64; 2],
    nbpnt_u: i32,
    uroot: &[f64],
    nbpnt_v: i32,
    vroot: &[f64],
    iiuouv: i32,
    sosotb: &mut [f64],
    disotb: &mut [f64],
    soditb: &mut [f64],
    diditb: &mut [f64],
) -> Result<(), String> {
    let alinu = (uintfn[1] - uintfn[0]) / 2.0;
    let blinu = (uintfn[1] + uintfn[0]) / 2.0;
    let alinv = (vintfn[1] - vintfn[0]) / 2.0;
    let blinv = (vintfn[1] + vintfn[0]) / 2.0;
    let (dbfn1, dbfn2) = if iiuouv == 1 {
        (vintfn, uintfn)
    } else {
        (uintfn, vintfn)
    };
    let nuroo = nbpnt_u / 2;
    let nvroo = nbpnt_v / 2;
    let jdec = (nbpnt_u + 1) / 2;
    let soso_dim1 = nuroo + 1;
    let soso_dim2 = nvroo + 1;
    let diso_dim1 = nuroo;
    let diso_dim2 = nvroo;
    let soso =
        |iu: i32, iv: i32, id: i32| iu + iv * soso_dim1 + (id - 1) * soso_dim1 * soso_dim2;
    let diso = |iu: i32, iv: i32, id: i32| {
        (iu - 1) + (iv - 1) * diso_dim1 + (id - 1) * diso_dim1 * diso_dim2
    };
    let fp = |id: i32, col: i32| (id - 1) + (col - 1) * ndimen;
    let mut ttable = vec![0.0; nbpnt_u.max(0) as usize];
    for iu in 1..=nbpnt_u {
        let root = *uroot.get((iu - 1) as usize).ok_or_else(ds2_oob)?;
        ttable[(iu - 1) as usize] = blinu + alinu * root;
    }
    let mut fpntab = vec![0.0; (ndimen * nbpnt_u.max(0)) as usize];
    let eval = |surf: &dyn Surface,
                tcons: f64,
                fpntab: &mut [f64]|
     -> Result<(), String> {
        let code = evaluate(
            surf, ndimen, dbfn1, dbfn2, iiuouv, tcons, &ttable, 0, 0, fpntab,
        )?;
        if code > 0 {
            return Err(format!(
                "AdvApp2Var_ApproxF2var::mma2ds2_: AdvApp2Var_ApproxF2var.cxx:6197 iercod={}",
                code + 100
            ));
        }
        Ok(())
    };
    for iv in 1..=nvroo {
        let root = *vroot.get((iv - 1) as usize).ok_or_else(ds2_oob)?;
        let tcons = blinv + alinv * root;
        let slot = nvroo - iv + 1;
        eval(surf, tcons, &mut fpntab)?;
        for id in 1..=ndimen {
            for iu in 1..=nuroo {
                let up = take_at(&fpntab, fp(id, iu + jdec))?;
                let um = take_at(&fpntab, fp(id, nuroo - iu + 1))?;
                let soso_i = soso(iu, slot, id);
                let diso_i = diso(iu, slot, id);
                add_at(sosotb, soso_i, up + um)?;
                add_at(disotb, diso_i, up - um)?;
                add_at(soditb, diso_i, -(up + um))?;
                add_at(diditb, soso_i, -up + um)?;
            }
            if nbpnt_u % 2 != 0 {
                let up = take_at(&fpntab, fp(id, jdec))?;
                let at = soso(0, slot, id);
                add_at(sosotb, at, up)?;
                add_at(diditb, at, -up)?;
            }
        }
    }
    if nbpnt_v % 2 != 0 {
        eval(surf, blinv, &mut fpntab)?;
        for id in 1..=ndimen {
            for iu in 1..=nuroo {
                let up = take_at(&fpntab, fp(id, jdec + iu))?;
                let um = take_at(&fpntab, fp(id, nuroo - iu + 1))?;
                let at = soso(iu, 0, id);
                add_at(sosotb, at, up + um)?;
                add_at(diditb, at, up - um)?;
            }
            if nbpnt_u % 2 != 0 {
                let up = take_at(&fpntab, fp(id, jdec))?;
                add_at(sosotb, soso(0, 0, id), up)?;
            }
        }
    }
    for iv in 1..=nvroo {
        let at = ((nbpnt_v + 1) / 2 + iv - 1) as usize;
        let root = *vroot.get(at).ok_or_else(ds2_oob)?;
        let tcons = alinv * root + blinv;
        eval(surf, tcons, &mut fpntab)?;
        for id in 1..=ndimen {
            for iu in 1..=nuroo {
                let up = take_at(&fpntab, fp(id, iu + jdec))?;
                let um = take_at(&fpntab, fp(id, nuroo - iu + 1))?;
                let soso_i = soso(iu, iv, id);
                let diso_i = diso(iu, iv, id);
                add_at(sosotb, soso_i, up + um)?;
                add_at(disotb, diso_i, up - um)?;
                add_at(soditb, diso_i, up + um)?;
                add_at(diditb, soso_i, up - um)?;
            }
            if nbpnt_u % 2 != 0 {
                let up = take_at(&fpntab, fp(id, jdec))?;
                let at = soso(0, iv, id);
                add_at(sosotb, at, up)?;
                add_at(diditb, at, up)?;
            }
        }
    }
    Ok(())
}
