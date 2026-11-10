//! Error constructors, `eps3` and `eps_column`.

use super::super::context::App2VarContext;


pub(crate) fn incoherent() -> String {
    "AdvApp2Var_ApproxF2var::mma2fnc_: AdvApp2Var_ApproxF2var.cxx:7251".into()
}

pub(crate) fn alloc_failed() -> String {
    "AdvApp2Var_ApproxF2var::mma2fnc_: AdvApp2Var_ApproxF2var.cxx:7247".into()
}

/// `mmveps3_` (`AdvApp2Var_MathBase.cxx:9171-9256`).
///
/// `mmprcsn_` is a zero-initialized static (`:221-225`). This path never
/// writes `eps3`.
pub(crate) fn eps3() -> f64 {
    0.0
}

pub(crate) fn eps_column(ctx: &App2VarContext, position: i32) -> Result<Vec<f64>, String> {
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
