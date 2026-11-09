//! File-local statics of `BRepOffset_MakeOffset::MakeShells`
//! (`BRepOffset_MakeOffset.cxx`): `UpdateHistory` (`:4649-4670`) and
//! `BuildShellsCompleteInter` (`:5064-5227`).
//!
//! `BuildShellsCompleteInter` builds the shells of the planar `Complete` arm
//! with `BOPAlgo_MakerVolume`: one pass over the split faces with
//! `SetIntersect(true)` (`:5071-5078`), and - only when that leaves a
//! non-manifold result - a second pass over the faces that belong to a single
//! solid (`:5159-5174`) and a third over the faces whose normal still agrees
//! with one of their origins (`:5214-5226`). Every early return keeps the
//! reference's `bDone = false`, which makes `MakeShells` fall back to
//! `BRepTools_Quilt`.
//!
//! `UpdateHistory` gives every split face of `theLF` its image in `theImage`,
//! through `Add` when the face already has one and `Bind` otherwise.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `Message_ProgressScope` (`:5069`) and the three `aPS.Next()` slices
//!   (`:5076`, `:5164`, `:5218`): the port has no progress machinery, so no
//!   `!More()` arm can fire. `BOPAlgo_MakerVolume::Perform(theRange)` is called
//!   without a range, as everywhere else in the port.
//! * the `NCollection_IncAllocator` handles (`:5095-5097`, `:5104`, `:5114`,
//!   `:5124`, `:5131`, `:5135`, `:5147`, `:5186`): Rust owns the containers, so
//!   `Reset` and `Bound(..., allocator)` have no counterpart.
//! * the `Standard_ASSERT_INVOKE` of `:5196-5197` (a face of the second result
//!   missing from the origins map): the port skips that face, which is the only
//!   behaviour the reference has past the assert in a release build.
//!
//! `anOrigins` maps a split face to the roots whose last image it is. OCCT
//! keeps the list in append order and reads it in that order (`:5201-5210`),
//! so the port keeps a `Vec` per key and never deduplicates: the same root can
//! be appended twice when `LastImage` lists the split twice.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::bop_maker_volume::BopMakerVolume;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::BRepOffsetTool;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors, shape_map_key, ShapeMapKey};

/// `theTolAng` default of `BRepOffset_Tool::CheckPlanesNormals`
/// (`BRepOffset_Tool.cxx:4509`).
const THE_TOL_ANG: f64 = 1.0e-6;

/// `void UpdateHistory(const List& theLF, BOPAlgo_Builder& theGF,
/// BRepAlgo_Image& theImage)` (`cxx:4649-4670`).
pub fn update_history(
    the_lf: &[TopoShape],
    the_gf: &BopMakerVolume,
    the_image: &mut BRepAlgoImage,
) {
    for a_f in the_lf {
        let a_lf_im: Vec<TopoShape> = the_gf.modified(a_f);
        if a_lf_im.is_empty() {
            continue;
        }
        // `cxx:4660-4667`.
        if the_image.has_image(a_f) {
            let _ = the_image.add_list(a_f, &a_lf_im);
        } else {
            let _ = the_image.bind_list(a_f, &a_lf_im);
        }
    }
}

/// `bool BuildShellsCompleteInter(const List& theLF, BRepAlgo_Image& theImage,
/// TopoDS_Shape& theShells, const Message_ProgressRange& theRange)`
/// (`cxx:5064-5227`).
///
/// Answers `theShells` on success and `None` where the reference returns
/// `false` (the caller then quilts instead).
pub fn build_shells_complete_inter(
    the_lf: &[TopoShape],
    the_image: &mut BRepAlgoImage,
) -> Option<TopoShape> {
    // `cxx:5071-5078`: the first volume pass, intersecting the arguments so
    // that the tangential faces are processed.
    let mut a_mv1 = BopMakerVolume::new();
    a_mv1.set_arguments(the_lf);
    a_mv1.set_intersect(true);
    a_mv1.set_avoid_internal_shapes(true);
    let _ = a_mv1.perform();
    if a_mv1.has_errors() {
        return None;
    }

    // `cxx:5084`.
    update_history(the_lf, &a_mv1, the_image);

    // `cxx:5086-5091`: a single solid needs no rebuilding.
    let a_result1 = a_mv1.shape().clone();
    if a_result1.shape_type() == ShapeType::Solid {
        return super::shape_prep::get_sub_shapes(&a_result1, ShapeType::Shell);
    }

    // `cxx:5099-5128`: `theImage` allows one ancestor per split, so a local map
    // collects every root whose `LastImage` is that split.
    let mut an_origins: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let a_roots = the_image.roots().to_vec();
    for a_f_r in &a_roots {
        let mut a_lf_im: Vec<TopoShape> = Vec::new();
        the_image.last_image(a_f_r, &mut a_lf_im);
        for a_f_im in a_lf_im {
            an_origins
                .entry(shape_map_key(&a_f_im))
                .or_default()
                .push(a_f_r.clone());
        }
    }

    // `cxx:5134-5144`: face -> the solids it belongs to.
    let a_dmfs = map_shapes_and_ancestors(&a_result1, ShapeType::Face, ShapeType::Solid);
    if a_dmfs.is_empty() {
        return None;
    }

    // `cxx:5147-5156`: the faces attached to exactly one solid.
    let mut a_lf: Vec<TopoShape> = Vec::new();
    for a_f in explorer_shapes(&a_result1, ShapeType::Face) {
        if a_dmfs.get(&shape_map_key(&a_f)).map_or(0, |l| l.len()) == 1 {
            a_lf.push(a_f);
        }
    }

    // `cxx:5159-5169`: the second pass, without intersecting.
    let mut a_mv2 = BopMakerVolume::new();
    a_mv2.set_arguments(&a_lf);
    a_mv2.set_intersect(false);
    a_mv2.set_avoid_internal_shapes(true);
    let _ = a_mv2.perform();
    if a_mv2.has_errors() {
        return None;
    }

    // `cxx:5171-5175`.
    let a_result2 = a_mv2.shape().clone();
    if a_result2.shape_type() == ShapeType::Solid {
        return super::shape_prep::get_sub_shapes(&a_result2, ShapeType::Shell);
    }

    // `cxx:5177-5182`.
    let a_faces2 = explorer_shapes(&a_result2, ShapeType::Face);
    if a_faces2.is_empty() {
        return None;
    }

    // `cxx:5188-5211`: the result is non-manifold, so only the faces whose
    // normal still agrees with one of their origins are kept.
    a_lf.clear();
    for a_f_shape in &a_faces2 {
        let Some(a_f) = Face::wrap(a_f_shape.clone()) else {
            continue;
        };
        let Some(a_l_origin) = an_origins.get(&shape_map_key(&a_f.0)) else {
            continue;
        };
        for a_f_or_shape in a_l_origin {
            let Some(a_f_or) = Face::wrap(a_f_or_shape.clone()) else {
                continue;
            };
            if BRepOffsetTool::check_planes_normals(&a_f, &a_f_or, THE_TOL_ANG) {
                a_lf.push(a_f.0.clone());
                break;
            }
        }
    }

    // `cxx:5214-5226`: the third pass over the outer faces.
    let mut a_mv3 = BopMakerVolume::new();
    a_mv3.set_arguments(&a_lf);
    a_mv3.set_intersect(false);
    a_mv3.set_avoid_internal_shapes(true);
    let _ = a_mv3.perform();
    if a_mv3.has_errors() {
        return None;
    }
    let a_result3 = a_mv3.shape().clone();
    super::shape_prep::get_sub_shapes(&a_result3, ShapeType::Shell)
}
