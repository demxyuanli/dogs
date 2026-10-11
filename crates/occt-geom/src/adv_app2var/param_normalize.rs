//! `pow__di` (`AdvApp2Var_MathBase.cxx:9530-9553`) and `mma1noc_`
//! (`AdvApp2Var_ApproxF2var.cxx:1490-1623`).

/// `pow__di` (`AdvApp2Var_MathBase.cxx:9530-9553`).
pub(super) fn pow_di(x: f64, n: i32) -> f64 {
    let mut result = 1.0;
    let absolute = if n > 0 { n } else { -n };
    for _ in 0..absolute {
        result *= x;
    }
    if n < 0 {
        result = 1.0 / result;
    }
    result
}

/// `mma1noc_` (`AdvApp2Var_ApproxF2var.cxx:1490-1623`).
pub(super) fn mma1noc(
    dfuvin: [f64; 4],
    ndimen: i32,
    iordre: i32,
    duvout: [f64; 4],
    isofav: i32,
    ideriv: i32,
    cntr: &mut [f64],
) -> Result<(), String> {
    let (num_der, den_der, num_ord, den_ord) = if isofav == 1 {
        (
            dfuvin[1] - dfuvin[0],
            duvout[1] - duvout[0],
            dfuvin[3] - dfuvin[2],
            duvout[3] - duvout[2],
        )
    } else {
        (
            dfuvin[3] - dfuvin[2],
            duvout[3] - duvout[2],
            dfuvin[1] - dfuvin[0],
            duvout[1] - duvout[0],
        )
    };
    let bid = pow_di(num_der / den_der, ideriv) * pow_di(num_ord / den_ord, iordre);
    if ndimen < 0 || ndimen as usize > cntr.len() {
        return Err("AdvApp2Var_ApproxF2var::mma1noc_: AdvApp2Var_ApproxF2var.cxx:1610".into());
    }
    for nd in 0..ndimen as usize {
        cntr[nd] *= bid;
    }
    Ok(())
}
