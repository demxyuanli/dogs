//! The three parameter helpers of `ChFi3d_Builder_CnCrn.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:150-201`
//! (OCCT 8.0.0): `Calcul_Param`, `Calcul_P2dOnSurf` and `Calcul_C2dOnFace`.
//!
//! All three read surf data `indice` (1-based) of a stripe and pick the
//! interference of `jfposit`. Note that the two selectors are written
//! differently in OCCT and the port keeps both: `Calcul_Param` passes
//! `jfposit == 2` to `S2` and everything else to `S1` (`cxx:156-164`), while
//! `Calcul_P2dOnSurf` / `Calcul_C2dOnFace` pass `jfposit == 1` to `S1` and
//! everything else to `S2` (`cxx:171-183`, `:188-201`).

use occt_core::gp::GpPnt2d;

use crate::chfi_ds::Stripe;

/// The surf data `indice` of `stripe`, `indice` 1-based (`Value(indice)`).
fn data(stripe: &Stripe, indice: i32) -> Result<&crate::chfi_ds::SurfData, String> {
    if indice < 1 {
        return Err(format!("Calcul: surf data index {indice} is out of range"));
    }
    stripe
        .set_of_surf_data()
        .get(indice as usize - 1)
        .map(|d| d.as_ref())
        .ok_or_else(|| format!("Calcul: surf data index {indice} is out of range"))
}

/// `Calcul_Param(stripe, jfposit, indice, isfirst, param)` (`cxx:150-165`): the
/// parameter of the common point of the interference of `jfposit`.
pub fn calcul_param(
    stripe: &Stripe,
    jfposit: i32,
    indice: i32,
    isfirst: bool,
) -> Result<f64, String> {
    let d = data(stripe, indice)?;
    let on_s = if jfposit == 2 { 2 } else { 1 };
    Ok(d.interference(on_s).parameter(isfirst))
}

/// `Calcul_P2dOnSurf(stripe, jfposit, indice, param, p2)` (`cxx:168-184`): the
/// UV point of the surface p-curve of the interference of `jfposit`.
///
/// OCCT dereferences `PCurveOnSurf()` unconditionally; the port answers `Err`
/// when the interference carries no surface p-curve.
pub fn calcul_p2d_on_surf(
    stripe: &Stripe,
    jfposit: i32,
    indice: i32,
    param: f64,
) -> Result<GpPnt2d, String> {
    let d = data(stripe, indice)?;
    let on_s = if jfposit == 1 { 1 } else { 2 };
    let curve = d
        .interference(on_s)
        .p_curve_on_surf()
        .ok_or_else(|| "Calcul_P2dOnSurf: the interference has no surface p-curve".to_string())?;
    Ok(curve.d0(param))
}

// `Calcul_C2dOnFace(stripe, jfposit, indice, c2d)` (`cxx:187-201`) has no
// caller before `cxx:2948` / `:2969` (the `GeomPlate` intermediary curves); it
// is ported next to those callers, in `plate_intermediary::calcul_c2d_on_face`.

