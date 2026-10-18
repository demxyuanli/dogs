//! `BOPAlgo_Builder::BuildDraftSolid`.
//!
//! Source: `BOPAlgo_Builder_3.cxx:267-368`. For every shell of the source
//! solid, each face is replaced by its image splits. INTERNAL faces (original
//! or split) are collected into `theLIF` and are *not* added to the draft
//! shell. Splits that have a same-domain representative (`myShapesSD.IsBound`)
//! are re-oriented with `IsSplitToReverseWithWarn`; other splits keep the
//! original face orientation.
//!
//! The previous port in [`crate::bop_build_common::build_draft_solid`] dropped
//! INTERNAL faces without returning `theLIF`, so ClassifyFaces never saw the
//! solid's own internal faces as `theSolidsIF`. This module returns both the
//! draft solid and `theLIF`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bop_occt_util::{
    add_shell_to_solid, as_face, builder_add, images_bound, images_of, is_internal,
    iter_children, make_empty_shell, make_empty_solid, shape_key, shells_of_solid, BopSolidHost,
};
use crate::bop_split_to_reverse::orient_split_with_warn;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Shell, Solid, TopoShape};
use crate::topo_tools::ShapeMapKey;

/// Whether `fx` is bound in the builder same-domain map (`myShapesSD.IsBound`).
fn shapes_sd_bound<B: BopSolidHost>(f: &B, fx: &TopoShape) -> bool {
    if f.seek_shapes_sd(fx).is_some() {
        return true;
    }
    f.ds()
        .index(fx)
        .and_then(|i| f.ds().has_shape_sd(i))
        .is_some()
}

/// Add `face` to `shell` once, skipping duplicates.
///
/// `BOPAlgo_Builder::BuildDraftSolid` adds every image face unconditionally
/// (`BOPAlgo_Builder_3.cxx:329` / `:342` / `:357`: `iFlag = 1; aBB.Add(aShD,
/// aFx)`). The previous body also dropped faces with fewer than three vertices;
/// that predicate does not exist in OCCT and is removed (audit A11). A duplicate
/// `TShape` cannot occur in OCCT's `myImages` list, so the identity guard is kept
/// only as a cheap invariant check.
fn add_draft_face_once(shell: &mut Shell, face: &TopoShape) -> bool {
    if as_face(face).is_none() {
        return false;
    }
    if iter_children(&shell.0)
        .iter()
        .any(|c| c.same_tshape(face))
    {
        return false;
    }
    builder_add(&mut shell.0, face);
    true
}

/// Result of `BuildDraftSolid`.
#[derive(Debug, Clone)]
pub struct DraftSolidResult {
    /// The rebuilt solid (`theDraftSolid`).
    pub draft: TopoShape,
    /// INTERNAL faces collected from the source (`theLIF`).
    pub internal_faces: Vec<TopoShape>,
}

/// `BOPAlgo_Builder::BuildDraftSolid` (`_3.cxx:267`).
pub fn build_draft_solid_occt<B: BopSolidHost>(
    f: &B,
    solid: &TopoShape,
    ctx: &mut IntToolsContext,
    add_warning: &mut dyn FnMut(String),
) -> DraftSolidResult {
    let or_sd = solid.orientation();
    let mut draft_val = make_empty_solid(or_sd);
    let mut lif: Vec<TopoShape> = Vec::new();

    let shells = if solid.shape_type() == ShapeType::Solid {
        shells_of_solid(solid)
    } else if solid.shape_type() == ShapeType::Shell {
        vec![Shell(solid.clone())]
    } else {
        Vec::new()
    };

    for sh in shells {
        let or_sh = sh.0.orientation();
        let mut sh_d = make_empty_shell(or_sh);
        let mut i_flag = false;

        for f_val in iter_children(&sh.0) {
            if f_val.shape_type() != ShapeType::Face {
                continue;
            }
            let or_f = f_val.orientation();
            if images_bound(f.history(), &f_val) {
                for fx_src in images_of(f.history(), &f_val) {
                    let mut fx = fx_src.clone();
                    if shapes_sd_bound(f, &fx) {
                        if is_internal(or_f) {
                            fx.set_orientation(or_f);
                            lif.push(fx);
                        } else if let Some(w) = orient_split_with_warn(&mut fx, &f_val, ctx) {
                            add_warning(w.message);
                            if add_draft_face_once(&mut sh_d, &fx) {
                                i_flag = true;
                            }
                        } else if add_draft_face_once(&mut sh_d, &fx) {
                            i_flag = true;
                        }
                    } else {
                        fx.set_orientation(or_f);
                        if is_internal(or_f) {
                            lif.push(fx);
                        } else if add_draft_face_once(&mut sh_d, &fx) {
                            i_flag = true;
                        }
                    }
                }
            } else if is_internal(or_f) {
                lif.push(f_val);
            } else if add_draft_face_once(&mut sh_d, &f_val) {
                i_flag = true;
            }
        }

        if i_flag {
            // `BOPAlgo_Builder_3.cxx:364`: `aShD.Closed(BRep_Tool::IsClosed(aShD))`.
            sh_d
                .0
                .set_closed(crate::brep_tool::BRepTool::is_closed(&sh_d.0));
            add_shell_to_solid(&mut draft_val, &sh_d);
        }
    }

    DraftSolidResult {
        draft: draft_val.0,
        internal_faces: lif,
    }
}

/// Host-agnostic draft solid used when only history is available (no SD map).
/// Same-domain splits are treated as unbound (`else` branch of `_3.cxx:332`).
pub fn build_draft_solid_from_history(
    history: &crate::bop_hist::BopHistory,
    solid: &TopoShape,
) -> DraftSolidResult {
    let or_sd = solid.orientation();
    let mut draft_val = make_empty_solid(or_sd);
    let mut lif: Vec<TopoShape> = Vec::new();
    let shells = if solid.shape_type() == ShapeType::Solid {
        shells_of_solid(solid)
    } else {
        Vec::new()
    };
    for sh in shells {
        let mut sh_d = make_empty_shell(sh.0.orientation());
        let mut i_flag = false;
        for f in iter_children(&sh.0) {
            if f.shape_type() != ShapeType::Face {
                continue;
            }
            let or_f = f.orientation();
            if images_bound(history, &f) {
                for fx_src in images_of(history, &f) {
                    let mut fx = fx_src.clone();
                    fx.set_orientation(or_f);
                    if is_internal(or_f) {
                        lif.push(fx);
                    } else if add_draft_face_once(&mut sh_d, &fx) {
                        i_flag = true;
                    }
                }
            } else if is_internal(or_f) {
                lif.push(f);
            } else if add_draft_face_once(&mut sh_d, &f) {
                i_flag = true;
            }
        }
        if i_flag {
            // `BOPAlgo_Builder_3.cxx:364`: `aShD.Closed(BRep_Tool::IsClosed(aShD))`.
            sh_d
                .0
                .set_closed(crate::brep_tool::BRepTool::is_closed(&sh_d.0));
            add_shell_to_solid(&mut draft_val, &sh_d);
        }
    }
    DraftSolidResult {
        draft: draft_val.0,
        internal_faces: lif,
    }
}

/// Diagnostics: number of faces on the draft solid.
pub fn draft_solid_face_count(draft: &TopoShape) -> usize {
    crate::topo_tools::faces_of(draft).len()
}

/// Diagnostics: whether the draft has any shell.
pub fn draft_solid_has_shell(draft: &TopoShape) -> bool {
    !shells_of_solid(draft).is_empty()
}

/// Shape keys of `theLIF`.
pub fn internal_face_keys(lif: &[TopoShape]) -> Vec<ShapeMapKey> {
    lif.iter().map(shape_key).collect()
}

/// Empty-solid constructor used by FillIn3DParts before BuildDraftSolid
/// (`_3.cxx:188` `aBB.MakeSolid(aSD)`).
pub fn make_solid_for_draft(or: Orientation) -> Solid {
    make_empty_solid(or)
}
