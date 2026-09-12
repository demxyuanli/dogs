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

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bop_occt_util::{
    add_shell_to_solid, as_face, builder_add, images_bound, images_of, is_internal,
    iter_children, make_empty_shell, make_empty_solid, shape_key, shells_of_solid, BopSolidHost,
};
use crate::bop_split_to_reverse::orient_split_with_warn;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Face, Shell, Solid, TopoShape};

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

/// True when `face` has no area (sliver). The existing draft path skipped
/// those; OCCT `BRep_Builder::Add` would still add them, but a zero-area
/// face cannot bound a shell, so the skip is kept as the degenerate-face
/// guard already used by `add_draft_face`.
fn face_is_degenerate(face: &Face) -> bool {
    crate::topo_tools_full::vertices_of(&face.0).len() < 3
}

/// Add `face` to `shell` once, skipping degenerates and duplicates.
fn add_draft_face_once(shell: &mut Shell, face: &TopoShape) -> bool {
    let Some(fc) = as_face(face) else {
        return false;
    };
    if face_is_degenerate(&fc) {
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
    the_solid: &TopoShape,
    ctx: &mut IntToolsContext,
    add_warning: &mut dyn FnMut(String),
) -> DraftSolidResult {
    let a_or_sd = the_solid.orientation();
    let mut the_draft = make_empty_solid(a_or_sd);
    let mut the_lif: Vec<TopoShape> = Vec::new();

    let shells = if the_solid.shape_type() == ShapeType::Solid {
        shells_of_solid(the_solid)
    } else if the_solid.shape_type() == ShapeType::Shell {
        vec![Shell(the_solid.clone())]
    } else {
        Vec::new()
    };

    for a_sh in shells {
        let a_or_sh = a_sh.0.orientation();
        let mut a_sh_d = make_empty_shell(a_or_sh);
        let mut i_flag = false;

        for a_f in iter_children(&a_sh.0) {
            if a_f.shape_type() != ShapeType::Face {
                continue;
            }
            let a_or_f = a_f.orientation();
            if images_bound(f.history(), &a_f) {
                for a_fx_src in images_of(f.history(), &a_f) {
                    let mut a_fx = a_fx_src.clone();
                    if shapes_sd_bound(f, &a_fx) {
                        if is_internal(a_or_f) {
                            a_fx.set_orientation(a_or_f);
                            the_lif.push(a_fx);
                        } else if let Some(w) = orient_split_with_warn(&mut a_fx, &a_f, ctx) {
                            add_warning(w.message);
                            if add_draft_face_once(&mut a_sh_d, &a_fx) {
                                i_flag = true;
                            }
                        } else if add_draft_face_once(&mut a_sh_d, &a_fx) {
                            i_flag = true;
                        }
                    } else {
                        a_fx.set_orientation(a_or_f);
                        if is_internal(a_or_f) {
                            the_lif.push(a_fx);
                        } else if add_draft_face_once(&mut a_sh_d, &a_fx) {
                            i_flag = true;
                        }
                    }
                }
            } else if is_internal(a_or_f) {
                the_lif.push(a_f);
            } else if add_draft_face_once(&mut a_sh_d, &a_f) {
                i_flag = true;
            }
        }

        if i_flag {
            a_sh_d
                .0
                .set_closed(!AlgoTools::is_open_shell(&a_sh_d.0));
            add_shell_to_solid(&mut the_draft, &a_sh_d);
        }
    }

    DraftSolidResult {
        draft: the_draft.0,
        internal_faces: the_lif,
    }
}

/// Host-agnostic draft solid used when only history is available (no SD map).
/// Same-domain splits are treated as unbound (`else` branch of `_3.cxx:332`).
pub fn build_draft_solid_from_history(
    history: &crate::bop_hist::BopHistory,
    the_solid: &TopoShape,
) -> DraftSolidResult {
    let a_or_sd = the_solid.orientation();
    let mut the_draft = make_empty_solid(a_or_sd);
    let mut the_lif: Vec<TopoShape> = Vec::new();
    let shells = if the_solid.shape_type() == ShapeType::Solid {
        shells_of_solid(the_solid)
    } else {
        Vec::new()
    };
    for a_sh in shells {
        let mut a_sh_d = make_empty_shell(a_sh.0.orientation());
        let mut i_flag = false;
        for a_f in iter_children(&a_sh.0) {
            if a_f.shape_type() != ShapeType::Face {
                continue;
            }
            let a_or_f = a_f.orientation();
            if images_bound(history, &a_f) {
                for a_fx_src in images_of(history, &a_f) {
                    let mut a_fx = a_fx_src.clone();
                    a_fx.set_orientation(a_or_f);
                    if is_internal(a_or_f) {
                        the_lif.push(a_fx);
                    } else if add_draft_face_once(&mut a_sh_d, &a_fx) {
                        i_flag = true;
                    }
                }
            } else if is_internal(a_or_f) {
                the_lif.push(a_f);
            } else if add_draft_face_once(&mut a_sh_d, &a_f) {
                i_flag = true;
            }
        }
        if i_flag {
            a_sh_d
                .0
                .set_closed(!AlgoTools::is_open_shell(&a_sh_d.0));
            add_shell_to_solid(&mut the_draft, &a_sh_d);
        }
    }
    DraftSolidResult {
        draft: the_draft.0,
        internal_faces: the_lif,
    }
}

/// Diagnostics: number of faces on the draft solid.
pub fn draft_solid_face_count(draft: &TopoShape) -> usize {
    crate::topo_tools_full::faces_of(draft).len()
}

/// Diagnostics: whether the draft has any shell.
pub fn draft_solid_has_shell(draft: &TopoShape) -> bool {
    !shells_of_solid(draft).is_empty()
}

/// TShape keys of `theLIF`.
pub fn internal_face_keys(lif: &[TopoShape]) -> Vec<usize> {
    lif.iter().map(shape_key).collect()
}

/// Empty-solid constructor used by FillIn3DParts before BuildDraftSolid
/// (`_3.cxx:188` `aBB.MakeSolid(aSD)`).
pub fn make_solid_for_draft(or: Orientation) -> Solid {
    make_empty_solid(or)
}
