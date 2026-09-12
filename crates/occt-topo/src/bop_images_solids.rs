//! `BOPAlgo_Builder::FillImagesSolids` — the three-step GF solid rebuild.
//!
//! Source: `BOPAlgo_Builder_3.cxx:70-93`:
//!
//! ```text
//! FillIn3DParts(aDraftSolids, theRange);
//! BuildSplitSolids(aDraftSolids, theRange);
//! FillInternalShapes(theRange);
//! ```
//!
//! This is the General Fuse path. Operation-state filtering (IN/OUT of
//! objects vs tools) is **not** applied here; OCCT does that later in
//! `BOPAlgo_BOP::BuildShape` (`BuildRC` / `BuildSolid` / `BuildBOP`).
//! `BOPAlgo_BOP::BuildResult` is the same GF assembly as
//! `BOPAlgo_Builder::BuildResult` (every argument image, no state skip).
//! The leftover `build_split_solids_full` mixed the two, which is why
//! overlapping boxes could emit two solids instead of the one GF image.
//!
//! Host: [`crate::bop_occt_util::BopSolidHost`]. `BopBuilder` implements it.

use crate::bop_fill_in3d::fill_in_3d_parts_builder;
use crate::bop_fill_internals_occt::fill_internal_shapes_occt;
use crate::bop_occt_util::BopSolidHost;
use crate::bop_split_solids_occt::build_split_solids_occt;
use crate::int_tools_full::IntToolsContext;

/// `BOPAlgo_Builder::FillImagesSolids`.
pub fn fill_images_solids_occt<B: BopSolidHost>(f: &mut B) -> Result<(), String> {
    let mut ctx = IntToolsContext::new();
    let fill = fill_in_3d_parts_builder(f, &mut ctx);
    build_split_solids_occt(f, &fill)?;
    fill_internal_shapes_occt(f)?;
    Ok(())
}

/// True when the DS contains at least one source solid. FillImagesSolids is
/// a no-op otherwise (`_3.cxx` loops `NbSourceShapes` and returns when none
/// are SOLID).
pub fn has_source_solids<B: BopSolidHost>(f: &B) -> bool {
    let n = f.ds().nb_source_shapes();
    (0..n).any(|i| {
        f.ds()
            .shape_info(i)
            .map(|s| s.shape_type() == crate::abs::ShapeType::Solid)
            .unwrap_or(false)
    })
}

/// FillImagesSolids that skips the three-step when there are no solids.
pub fn fill_images_solids_if_needed<B: BopSolidHost>(f: &mut B) -> Result<(), String> {
    if !has_source_solids(f) {
        return Ok(());
    }
    fill_images_solids_occt(f)
}
