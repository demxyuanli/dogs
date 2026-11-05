//! The `PerformSetOfKPart` call-out to `ChFiKPart_ComputeData::Compute`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:3070-3081`
//! (OCCT 8.0.0).
//!
//! [`crate::chfi_k_part::ch_fi_kpart_compute_data`] is the port of `Compute`
//! itself; this module is the thin wire the builder uses so that the branch
//! below stays readable:
//!
//! ```text
//! if (!ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2, Spine, iedge)) { }
//! else if (!SplitKPart(SD, LSD, Spine, iedge, HS1, It1, HS2, It2, intf, intl)) {
//!   LSD.Clear();
//! }
//! else { iedgelastkpart = iedge; }
//! ```
//!
//! `Compute`'s answer is only ever read by that `else if`: a `true` sends the
//! `SurfData` into [`crate::chfi3d::builder::split_k_part`]
//! (`ChFi3d_Builder_SpKP.cxx:749-1515`). Both halves are ported, so the wire
//! forwards the answer unchanged.

use crate::chfi_ds::{SpineHandle, SurfData};
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2, Spine, IEdge)`
/// (`ChFiKPart_ComputeData.cxx:52-637`) — the dispatcher lives in
/// [`crate::chfi_k_part`].
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_compute_data(
    ds: &mut DataStructure,
    sd: &mut SurfData,
    hs1: &SurfaceAdaptor,
    hs2: &SurfaceAdaptor,
    or1: crate::abs::Orientation,
    or2: crate::abs::Orientation,
    spine: &SpineHandle,
    iedge: i32,
) -> Result<bool, String> {
    crate::chfi_k_part::ch_fi_kpart_compute_data(
        ds, sd, hs1, hs2, or1, or2, spine, iedge,
    )
}
