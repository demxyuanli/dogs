//! The two index helpers of `ChFi3d_Builder_CnCrn.cxx`: `Indices` and
//! `SurfIndex`, with its `ChFi3d_SurfType` selector.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:128-146`
//! and `:1143-1170` (OCCT 8.0.0).

use crate::chfi_ds::StripeHandle;

/// `ChFi3d_SurfType` (`cxx:1143-1148`), the selector of [`surf_index`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfType {
    /// `ChFiSURFACE` — `SurfData::Surf()`. Read by the intersection / plate
    /// stages (`SurfIndex(CD, ic, i.Value(...), ChFiSURFACE)` at `cxx:1982`,
    /// `:2084`, `:2097`, `:2670-2730`).
    ChFiSurface,
    /// `FACE1` — `SurfData::IndexOfS1()`.
    Face1,
    /// `FACE2` — `SurfData::IndexOfS2()`.
    Face2,
}

/// `Indices(n, ic, icplus, icmoins)` (`cxx:128-146`): the two neighbours of
/// `ic` on the ring of `n` entries, wrapping at both ends.
pub fn indices(n: usize, ic: usize) -> (usize, usize) {
    let icplus = if ic == n - 1 { 0 } else { ic + 1 };
    let icmoins = if ic == 0 { n - 1 } else { ic - 1 };
    (icplus, icmoins)
}

/// `SurfIndex(StripeArray1, StripeIndex, SurfDataIndex, SurfType)`
/// (`cxx:1152-1170`): one of the three shape indices of a `SurfData`.
///
/// `surf_data_index` is 1-based like `NCollection_Sequence::Value`; OCCT raises
/// on an out-of-range index or a null stripe, the port answers `Result`.
pub fn surf_index(
    stripes: &[Option<StripeHandle>],
    stripe_index: usize,
    surf_data_index: usize,
    surf_type: SurfType,
) -> Result<i32, String> {
    let stripe = stripes
        .get(stripe_index)
        .and_then(|s| s.as_ref())
        .ok_or_else(|| format!("SurfIndex: CD({stripe_index}) is null"))?;
    let guard = stripe
        .read()
        .map_err(|_| "SurfIndex: poisoned stripe".to_string())?;
    let data = surf_data_index
        .checked_sub(1)
        .and_then(|i| guard.set_of_surf_data().get(i))
        .ok_or_else(|| format!("SurfIndex: SurfData {surf_data_index} is out of range"))?;
    Ok(match surf_type {
        SurfType::ChFiSurface => data.surf(),
        SurfType::Face1 => data.index_of_s1(),
        SurfType::Face2 => data.index_of_s2(),
    })
}
