//! `AdvApp2Var_ApproxF2var::mma2jmx_` (`AdvApp2Var_ApproxF2var.cxx:7576-7763`).
//!
//! `AdvApp2Var_SysBase::mnfndeb_` returns 0 (`AdvApp2Var_SysBase.cxx:2887-2891`),
//! so the debug messages around the table copy are not reached.

include!("jacobi_tables.rs");

/// `mma2jmx_`: `numax = ndgjac - 2 * (iordre + 1)`, then copy that prefix.
pub(super) fn jacobi_max(ndgjac: i32, iordre: i32) -> Result<Vec<f64>, String> {
    let numax = ndgjac - ((iordre + 1) << 1);
    if numax < 0 {
        return Ok(Vec::new());
    }
    if iordre == -1 {
        return Ok((0..=numax)
            .map(|ii| (((ii as f64) * 2.0 + 1.0) / 2.0).sqrt())
            .collect());
    }
    let table: &[f64] = match iordre {
        0 => &XMAX2,
        1 => &XMAX4,
        2 => &XMAX6,
        _ => {
            return Err(
                "AdvApp2Var_ApproxF2var::mma2jmx_: AdvApp2Var_ApproxF2var.cxx:7717".into(),
            );
        }
    };
    let last = numax as usize;
    if last >= table.len() {
        return Err("AdvApp2Var_ApproxF2var::mma2jmx_: AdvApp2Var_ApproxF2var.cxx:7727".into());
    }
    Ok(table[..=last].to_vec())
}
