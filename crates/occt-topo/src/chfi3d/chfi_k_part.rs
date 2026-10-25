//! `ChFiKPart_ComputeData` — **UNPORTED**.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData.cxx`
//! (OCCT 8.0.0), the package of the analytic corner surfaces (`ChFiKPart` is
//! 4,938 `.cxx` lines over 5 classes and belongs to the `ChFiKPart` cell of the
//! S3 wave, not to `ChFi3d`).
//!
//! `PerformSetOfKPart` calls
//! `ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2, Spine, iedge)`
//! once per edge that [`crate::chfi3d::builder_0_helpers::ch_fi3d_k_particular`]
//! accepted (`ChFi3d_Builder_2.cxx:3071-3073`).
//!
//! Until that package is ported this function answers `Ok(false)`, which is
//! exactly the value OCCT's own `if (!Compute(...)) { }` arm handles: no
//! `SurfData` is filled, so `SplitKPart` is not called and no `ElSpine` is
//! appended. The placeholder is therefore inert rather than wrong — the KPart
//! branch of `PerformSetOfKPart` produces nothing, and the general path of
//! `PerformSetOfKGen` is the only one with data.

use crate::chfi_ds::{SpineHandle, SurfData};
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2, Spine, IEdge)`
/// (`ChFiKPart_ComputeData.cxx`) — **UNPORTED**, answers `Ok(false)`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_compute_data(
    _ds: &mut DataStructure,
    _sd: &mut SurfData,
    _hs1: &SurfaceAdaptor,
    _hs2: &SurfaceAdaptor,
    _or1: crate::abs::Orientation,
    _or2: crate::abs::Orientation,
    _spine: &SpineHandle,
    _iedge: i32,
) -> Result<bool, String> {
    Ok(false)
}
