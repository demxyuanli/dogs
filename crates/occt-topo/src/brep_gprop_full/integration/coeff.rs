
pub(in crate::brep_gprop_full) fn s_coeff(eps: f64) -> f64 {
    if eps < 0.1 {
        -0.15 * (1.0 + eps.log10()) + 0.75
    } else {
        0.75
    }
}

#[allow(dead_code)]
pub(in crate::brep_gprop_full) fn l_coeff(eps: f64) -> f64 {
    if eps < 0.1 {
        -0.50 * (1.0 + eps.log10()) + 0.25
    } else {
        0.25
    }
}
