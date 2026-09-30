//! `BOPTools_AlgoTools::IsSplitToReverseWithWarn`.
//!
//! Source: `BOPTools_AlgoTools.cxx` (`IsSplitToReverse` plus the WithWarn
//! wrapper that reports `BOPAlgo_AlertUnableToOrientTheShape` when the
//! classifier cannot decide). BuildSplitFaces (`_2.cxx:458`) and
//! BuildDraftFace (`_2.cxx:1169`) and BuildDraftSolid (`_3.cxx:321`) all call
//! the WithWarn form so a failed orientation test is a warning, not a skip
//! of the split.
//!
//! The geometric test itself is [`crate::algo_tools_face::is_split_to_reverse`].
//! This module only adds the warning path and a few orientation helpers the
//! Builder stages share.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use crate::abs::Orientation;
use crate::algo_tools_face::is_split_to_reverse;
use crate::int_tools_full::IntToolsContext;
use crate::shape::TopoShape;

/// Warning payload when `IsSplitToReverse` returns a non-zero error code.
#[derive(Debug, Clone)]
pub struct ReverseWarn {
    /// OCCT error code from `IsSplitToReverse` (`iErr != 0`).
    pub error: i32,
    /// Human-readable alert matching `BOPAlgo_AlertUnableToOrientTheShape`.
    pub message: String,
}

impl ReverseWarn {
    fn from_code(err: i32, split: &TopoShape, origin: &TopoShape) -> Self {
        Self {
            error: err,
            message: format!(
                "BOPAlgo_AlertUnableToOrientTheShape: IsSplitToReverse failed ({err}) for {:?} vs {:?}",
                split.shape_type(),
                origin.shape_type()
            ),
        }
    }
}

/// `BOPTools_AlgoTools::IsSplitToReverseWithWarn`.
///
/// Returns `(to_reverse, warning)`. `to_reverse` is false when the classifier
/// errors; the warning is `Some` in that case so the Builder can
/// `AddWarning` and still keep the split at its current orientation
/// (OCCT does not drop the split).
pub fn is_split_to_reverse_with_warn(
    split: &TopoShape,
    origin: &TopoShape,
    ctx: &mut IntToolsContext,
) -> (bool, Option<ReverseWarn>) {
    match is_split_to_reverse(split, origin, ctx) {
        Ok(flag) => (flag, None),
        Err(err) => (false, Some(ReverseWarn::from_code(err, split, origin))),
    }
}

/// Apply WithWarn: reverse `split` in place when the test says so, and
/// return the optional warning for the caller to record.
pub fn orient_split_with_warn(
    split: &mut TopoShape,
    origin: &TopoShape,
    ctx: &mut IntToolsContext,
) -> Option<ReverseWarn> {
    let (to_reverse, warn) = is_split_to_reverse_with_warn(split, origin, ctx);
    if to_reverse {
        split.reverse();
    }
    warn
}

/// Copy `split`, apply WithWarn, return the oriented copy plus warning.
pub fn oriented_split_with_warn(
    split: &TopoShape,
    origin: &TopoShape,
    ctx: &mut IntToolsContext,
) -> (TopoShape, Option<ReverseWarn>) {
    let mut sp = split.clone();
    let warn = orient_split_with_warn(&mut sp, origin, ctx);
    (sp, warn)
}

/// WithWarn that also forces `split` to start at orientation `base` before
/// the reverse test. BuildSplitFaces sets the split to the original edge
/// orientation first (`aSp.Orientation(anOriE)`), then possibly reverses.
pub fn orient_split_from_base_with_warn(
    split: &TopoShape,
    origin: &TopoShape,
    base: Orientation,
    ctx: &mut IntToolsContext,
) -> (TopoShape, Option<ReverseWarn>) {
    let mut sp = split.clone();
    sp.set_orientation(base);
    let warn = orient_split_with_warn(&mut sp, origin, ctx);
    (sp, warn)
}

/// Collect warnings from a batch of (split, origin) pairs. Used by the
/// draft-face / draft-solid loops that process an image list.
pub fn orient_image_list_with_warn(
    images: &[TopoShape],
    origin: &TopoShape,
    ctx: &mut IntToolsContext,
) -> (Vec<TopoShape>, Vec<ReverseWarn>) {
    let mut out = Vec::with_capacity(images.len());
    let mut warns = Vec::new();
    for im in images {
        let (sp, w) = oriented_split_with_warn(im, origin, ctx);
        if let Some(w) = w {
            warns.push(w);
        }
        out.push(sp);
    }
    (out, warns)
}

/// True when `split` and `origin` have opposite face/edge orientation after
/// the geometric test (convenience for callers that only need the flag).
pub fn split_needs_reverse(
    split: &TopoShape,
    origin: &TopoShape,
    ctx: &mut IntToolsContext,
) -> bool {
    is_split_to_reverse_with_warn(split, origin, ctx).0
}

/// Format a warning list into Builder alert strings.
pub fn warn_messages(warns: &[ReverseWarn]) -> Vec<String> {
    warns.iter().map(|w| w.message.clone()).collect()
}

/// Push each warning message onto a host `add_warning` callback.
pub fn emit_warnings<F: FnMut(String)>(warns: &[ReverseWarn], mut add_warning: F) {
    for w in warns {
        add_warning(w.message.clone());
    }
}

/// Decide the stored orientation of a split given the original orientation
/// `an_ori` and the WithWarn result. Mirrors the BuildSplitFaces closed-edge
/// path where both FORWARD and REVERSED copies are appended without reverse
/// testing, versus the open-edge path that tests and possibly reverses.
pub fn split_orientation_for_open_edge(
    split: &TopoShape,
    origin: &TopoShape,
    an_ori: Orientation,
    ctx: &mut IntToolsContext,
) -> (Orientation, Option<ReverseWarn>) {
    let mut sp = split.clone();
    sp.set_orientation(an_ori);
    let (flag, warn) = is_split_to_reverse_with_warn(&sp, origin, ctx);
    let or = if flag { an_ori.reversed() } else { an_ori };
    (or, warn)
}
