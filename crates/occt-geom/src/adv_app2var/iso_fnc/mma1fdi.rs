//! `mma1fdi_` (`cxx:596-981`) and `fill_constraint`.

use crate::Surface;


use super::super::eval_surface::evaluate;

use super::*;

/// `mma1fdi_` (`AdvApp2Var_ApproxF2var.cxx:596-981`).
pub(crate) fn mma1fdi(
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
        (ideriv, 0, super::super::param_normalize::pow_di(half, ideriv))
    } else {
        let half = (uvfonc[3] - uvfonc[2]) / 2.0;
        (0, ideriv, super::super::param_normalize::pow_di(half, ideriv))
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

pub(crate) fn fill_constraint(
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
