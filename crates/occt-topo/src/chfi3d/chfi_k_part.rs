//! The `PerformSetOfKPart` call-out to `ChFiKPart_ComputeData::Compute`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:3070-3080`
//! (OCCT 8.0.0).
//!
//! [`crate::chfi_k_part::ch_fi_kpart_compute_data`] is the port of `Compute`
//! itself. This wrapper is the *wire* the builder uses, and deliberately keeps
//! answering `Ok(false)` until `SplitKPart` is ported:
//!
//! ```text
//! if (!ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2, Spine, iedge)) { }
//! else if (!SplitKPart(SD, LSD, Spine, iedge, HS1, It1, HS2, It2, intf, intl)) {
//!   LSD.Clear();
//! }
//! else { iedgelastkpart = iedge; }
//! ```
//!
//! `Compute`'s answer is only ever read by that `else if`: a `true` from here
//! sends the `SurfData` into [`crate::chfi3d::builder::split_k_part`], which is
//! UNPORTED (`ChFi3d_Builder_SpKP.cxx:749-1160`) and answers `Err`. Wiring the
//! ported `Compute` now would therefore turn every plane/plane, plane/cylinder
//! and plane/cone fillet edge from the empty arm into a failure, which is not
//! what OCCT does. `Ok(false)` is the empty arm, and the corner is simply not
//! built yet.

use crate::chfi_ds::{SpineHandle, SurfData};
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2, Spine, IEdge)`
/// (`ChFiKPart_ComputeData.cxx:52-637`) — the builder-side wire, see the module
/// comment for why it is still cut.
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
