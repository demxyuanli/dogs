//! Builder-level `BOPAlgo_Builder::FillIn3DParts`.
//!
//! Source: `BOPAlgo_Builder_3.cxx:97-263`. This is *not* the same as
//! `BOPAlgo_FillIn3DParts::Perform` in Tools.cxx (the per-solid classifier).
//! The Builder method:
//! 1. Collects every source FACE, replaced by its image splits (`aLFaces`);
//! 2. For every source SOLID, builds a draft via
//!    [`crate::bop_draft_solid_occt::build_draft_solid_occt`] and keeps the
//!    INTERNAL faces as `aSolidsIF`;
//! 3. Calls [`crate::bop_classify_occt::classify_faces_occt`];
//! 4. Binds `theDraftSolids` for solids that have IN faces or split shells,
//!    and writes `myInParts` = IN faces + INTERNAL faces.
//!
//! `BopBuilder::fill_images_solids` passes the two maps into
//! [`crate::bop_split_solids_occt::build_split_solids_occt`].
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;

use crate::abs::ShapeType;
use crate::bop_classify_occt::{classify_faces_occt, combine_in_and_internal, solid_needs_draft};
use crate::bop_draft_solid_occt::build_draft_solid_occt;
use crate::bop_occt_util::{
    collect_candidate_faces, nb_source, shape_box_of, shape_key, source_solids, BopSolidHost,
};
use crate::int_tools_full::IntToolsContext;
use crate::shape::TopoShape;

/// Output of Builder `FillIn3DParts`.
#[derive(Debug, Clone, Default)]
pub struct FillIn3dPartsResult {
    /// `theDraftSolids`: source solid TShape key → draft solid.
    pub draft_solids: HashMap<usize, TopoShape>,
    /// `myInParts`: source solid TShape key → IN + INTERNAL faces.
    pub in_parts: HashMap<usize, Vec<TopoShape>>,
    /// Draft solid TShape key → source solid (reverse of `aDraftSolid`).
    pub draft_to_source: HashMap<usize, TopoShape>,
}

/// `BOPAlgo_Builder::FillIn3DParts` (`_3.cxx:97`).
pub fn fill_in_3d_parts_builder<B: BopSolidHost>(
    f: &mut B,
    ctx: &mut IntToolsContext,
) -> FillIn3dPartsResult {
    let mut result = FillIn3dPartsResult::default();
    let mut a_shape_box_map: HashMap<usize, BndBox> = HashMap::new();

    // 1. All faces (images preferred).
    let a_l_faces = collect_candidate_faces(f.history(), f.ds());
    for face in &a_l_faces {
        if f.ds().index(face).is_some() {
            a_shape_box_map.insert(shape_key(face), shape_box_of(f.ds(), face));
        }
    }

    // 2. Draft solids + INTERNAL faces.
    let mut a_l_solids: Vec<TopoShape> = Vec::new();
    let mut a_solids_if: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    // source key → draft shape (OCCT `aDraftSolid` IndexedDataMap).
    let mut a_draft_solid: Vec<(TopoShape, TopoShape)> = Vec::new();

    let n = nb_source(f.ds());
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let a_s = si.shape().clone();
        let mut a_box_s = shape_box_of(f.ds(), &a_s);
        if a_box_s.is_void() {
            a_box_s = crate::bbox_from_geometry::shape_bbox(&a_s);
        }
        let mut warnings: Vec<String> = Vec::new();
        let draft = build_draft_solid_occt(f, &a_s, ctx, &mut |msg| warnings.push(msg));
        for w in warnings {
            f.add_warning(w);
        }
        let a_sd = draft.draft;
        a_l_solids.push(a_sd.clone());
        a_solids_if.insert(shape_key(&a_sd), draft.internal_faces);
        a_shape_box_map.insert(shape_key(&a_sd), a_box_s);
        a_draft_solid.push((a_s, a_sd));
    }

    // 3. Classify faces relatively the draft solids.
    let an_in_parts = classify_faces_occt(
        &a_l_faces,
        &a_l_solids,
        ctx,
        &a_shape_box_map,
        &a_solids_if,
    );

    // 4. Analyze: bind drafts that need splitting, fill myInParts.
    for (a_solid, a_s_draft) in a_draft_solid {
        let a_l_in = an_in_parts
            .get(&shape_key(&a_s_draft))
            .cloned()
            .unwrap_or_default();
        let a_l_internal = a_solids_if
            .get(&shape_key(&a_s_draft))
            .cloned()
            .unwrap_or_default();
        let a_nb_in = a_l_in.len();
        let need = solid_needs_draft(f.history(), &a_solid, a_nb_in);
        if !need {
            continue;
        }
        result
            .draft_solids
            .insert(shape_key(&a_solid), a_s_draft.clone());
        result
            .draft_to_source
            .insert(shape_key(&a_s_draft), a_solid.clone());
        if !a_l_internal.is_empty() || a_nb_in != 0 {
            result.in_parts.insert(
                shape_key(&a_solid),
                combine_in_and_internal(&a_l_in, &a_l_internal),
            );
        }
    }
    let _ = source_solids(f.ds());
    result
}

/// Faces of `in_parts` for one source solid, or empty.
pub fn in_parts_of<'a>(
    result: &'a FillIn3dPartsResult,
    source: &TopoShape,
) -> &'a [TopoShape] {
    result
        .in_parts
        .get(&shape_key(source))
        .map(|v| v.as_slice())
        .unwrap_or(&[])
}

/// Draft of `source`, when FillIn3DParts bound one.
pub fn draft_of<'a>(result: &'a FillIn3dPartsResult, source: &TopoShape) -> Option<&'a TopoShape> {
    result.draft_solids.get(&shape_key(source))
}

/// Unique TShape keys of every IN face across all solids.
pub fn all_in_face_keys(result: &FillIn3dPartsResult) -> HashSet<usize> {
    let mut s = HashSet::new();
    for faces in result.in_parts.values() {
        for f in faces {
            s.insert(shape_key(f));
        }
    }
    s
}
