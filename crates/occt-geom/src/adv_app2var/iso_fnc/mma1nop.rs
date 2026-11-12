//! `mma1nop_` (`cxx:1627-1734`) and its index helpers.

use super::*;

/// `mma1nop_` (`AdvApp2Var_ApproxF2var.cxx:1627-1734`).
pub(crate) fn mma1nop(
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

pub(crate) fn sample_at(point: i32, nd: i32, dim: i32) -> usize {
    ((nd - 1) + point * dim) as usize
}

pub(crate) fn pair_at(slot: i32, nd: i32, dim1: i32) -> usize {
    (slot + (nd - 1) * dim1) as usize
}

pub(crate) fn contr_at(order: i32, nd: i32, dim: i32) -> usize {
    (order * dim + (nd - 1)) as usize
}
