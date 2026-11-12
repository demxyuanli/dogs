//! Error accumulation over the packed `NDIMSE` subspace layout.

#[allow(clippy::too_many_arguments)]
pub(super) fn subspace_errors(
    ndjacu: i32,
    ndjacv: i32,
    dimensions: &[i32],
    mindgu: i32,
    maxdgu: i32,
    mindgv: i32,
    maxdgv: i32,
    iordru: i32,
    iordrv: i32,
    xmaxju: &[f64],
    xmaxjv: &[f64],
    patjac: &[f64],
    previous: &[f64],
) -> Result<Vec<f64>, String> {
    if dimensions.len() != previous.len() {
        return Err("AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4369".into());
    }
    let plane = usize::try_from((ndjacu + 1) * (ndjacv + 1))
        .map_err(|_| "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4369")?;
    let mut offset = 0usize;
    let mut errors = Vec::with_capacity(dimensions.len());
    for (subspace, &width) in dimensions.iter().enumerate() {
        let len = usize::try_from(width)
            .ok()
            .and_then(|n| n.checked_mul(plane))
            .ok_or_else(|| {
                "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4369".to_string()
            })?;
        let end = offset.checked_add(len).ok_or_else(|| {
            "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4369".to_string()
        })?;
        let coefficients = patjac.get(offset..end).ok_or_else(|| {
            "AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4369".to_string()
        })?;
        errors.push(super::degree_error::mma2er1(
            ndjacu,
            ndjacv,
            width,
            mindgu,
            maxdgu,
            mindgv,
            maxdgv,
            iordru,
            iordrv,
            xmaxju,
            xmaxjv,
            coefficients,
            previous[subspace],
        )?);
        offset = end;
    }
    Ok(errors)
}

pub(super) fn exceeds(errors: &[f64], tolerances: &[f64]) -> Result<bool, String> {
    if errors.len() != tolerances.len() {
        return Err("AdvApp2Var_ApproxF2var::mma2ce2_: AdvApp2Var_ApproxF2var.cxx:4565".into());
    }
    Ok(errors
        .iter()
        .zip(tolerances)
        .any(|(error, tolerance)| error > tolerance))
}
