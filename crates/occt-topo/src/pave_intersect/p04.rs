use super::prelude::*;
use super::*;

/// `BOPAlgo_PaveFiller::PutSEInOtherFaces` (`BOPAlgo_PaveFiller_6.cxx:4277`).
/// Intersects every section pave block with faces that did not create it.
pub(crate) fn put_se_in_other_faces(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ff_se::put_se_in_other_faces(f)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Edge/Face intersection: intersect every interfering edge with its face,
/// creating new vertices at the piercing points and recording coincident
/// sub-ranges as face paves.
///
/// Source: `BOPAlgo_PaveFiller::PerformEF` + `IntersectEF`
/// (`BOPAlgo_PaveFiller_5.cxx`).
pub fn perform_ef(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ef_perform::perform_ef(f)
}

/// Face/Face intersection: intersect every interfering pair of faces and build
/// a section edge for each intersection curve, recording the edge on both
/// faces and the F/F interference.
///
/// Source: `BOPAlgo_PaveFiller::PerformFF` (`BOPAlgo_PaveFiller_6.cxx`).
pub fn perform_ff(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ff_perform::perform_ff(f)
}
