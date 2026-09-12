//! `BOPAlgo_Builder::FillInternalShapes`.
//!
//! Source: `BOPAlgo_Builder_3.cxx:622-887` plus `OwnInternalShapes` at 891.
//! After split solids are recorded as images, vertices and edges that belong
//! to the arguments (or that sit inside a source solid as INTERNAL children)
//! are classified with `ComputeStateByOnePoint` and, when IN, added to the
//! split solid as INTERNAL. An original (unsplit) solid is copied first so
//! the input argument is not mutated (`aMSOr` branch, `_3.cxx:844-868`).
//!
//! The previous port lives in [`crate::bop_build_common::fill_internal_shapes`].
//! This module is the same algorithm expressed against [`BopSolidHost`] so
//! FillImagesSolids can run the OCCT three-step without changing
//! `BopBuildOps`.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::bop_occt_util::{
    images_bound, images_of, iter_children, map_shapes_and_ancestors, merge_ancestor_maps,
    nb_source, origins_append, own_internal_shapes, shape_key, treat_compound, BopSolidHost,
};
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Solid, TopoShape};
use crate::topo_tools_full::edges_of;

/// `BOPTools_AlgoTools::ComputeStateByOnePoint` stand-in used by
/// FillInternalShapes (`_3.cxx:836`, tolerance `1.e-11`).
fn compute_state_by_one_point(shape: &TopoShape, solid: &TopoShape, tol: f64) -> FaceState {
    crate::bopalgo_tools_class::compute_state_by_one_point(shape, solid, tol)
}

/// `BOPAlgo_Builder::FillInternalShapes` (`_3.cxx:622`).
pub fn fill_internal_shapes_occt<B: BopSolidHost>(f: &mut B) -> Result<(), String> {
    // 1. Shapes to process from the pure arguments.
    let mut a_lsc: Vec<TopoShape> = Vec::new();
    let mut a_fence: HashSet<usize> = HashSet::new();
    for a in f.arguments() {
        treat_compound(a, &mut a_lsc, &mut a_fence);
    }
    let mut a_largs: Vec<TopoShape> = Vec::new();
    a_fence.clear();
    for s in &a_lsc {
        match s.shape_type() {
            ShapeType::Wire => {
                for e in edges_of(s) {
                    if a_fence.insert(shape_key(&e.0)) {
                        a_largs.push(e.0);
                    }
                }
            }
            ShapeType::Vertex | ShapeType::Edge => a_largs.push(s.clone()),
            _ => {}
        }
    }
    a_fence.clear();
    let mut a_msi: Vec<TopoShape> = Vec::new();
    for s in &a_largs {
        let t = s.shape_type();
        if t != ShapeType::Vertex && t != ShapeType::Edge && t != ShapeType::Wire {
            continue;
        }
        if !a_fence.insert(shape_key(s)) {
            continue;
        }
        if images_bound(f.history(), s) {
            for im in images_of(f.history(), s) {
                a_msi.push(im.clone());
            }
        } else {
            a_msi.push(s.clone());
        }
    }

    // 2. Internal vertices/edges from source solids + ancestor map of splits.
    a_fence.clear();
    let mut a_lsd: Vec<TopoShape> = Vec::new();
    let mut a_msx: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    let mut a_msor: HashSet<usize> = HashSet::new();
    let n = nb_source(f.ds());
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let a_s = si.shape().clone();
        for a_si in own_internal_shapes(&a_s) {
            if images_bound(f.history(), &a_si) {
                for sp in images_of(f.history(), &a_si) {
                    a_msi.push(sp.clone());
                }
            } else {
                a_msi.push(a_si);
            }
        }
        if images_bound(f.history(), &a_s) {
            for a_sp in images_of(f.history(), &a_s) {
                if a_fence.insert(shape_key(a_sp)) {
                    merge_ancestor_maps(
                        &mut a_msx,
                        map_shapes_and_ancestors(a_sp, ShapeType::Vertex, ShapeType::Edge),
                    );
                    merge_ancestor_maps(
                        &mut a_msx,
                        map_shapes_and_ancestors(a_sp, ShapeType::Vertex, ShapeType::Face),
                    );
                    merge_ancestor_maps(
                        &mut a_msx,
                        map_shapes_and_ancestors(a_sp, ShapeType::Edge, ShapeType::Face),
                    );
                    a_lsd.push(a_sp.clone());
                }
            }
        } else if a_fence.insert(shape_key(&a_s)) {
            merge_ancestor_maps(
                &mut a_msx,
                map_shapes_and_ancestors(&a_s, ShapeType::Vertex, ShapeType::Edge),
            );
            merge_ancestor_maps(
                &mut a_msx,
                map_shapes_and_ancestors(&a_s, ShapeType::Vertex, ShapeType::Face),
            );
            merge_ancestor_maps(
                &mut a_msx,
                map_shapes_and_ancestors(&a_s, ShapeType::Edge, ShapeType::Face),
            );
            a_lsd.push(a_s.clone());
            a_msor.insert(shape_key(&a_s));
        }
    }

    // 3. Drop candidates already tied to a split-solid ancestor.
    let mut a_lsi: Vec<TopoShape> = Vec::new();
    for a_si in &a_msi {
        match a_msx.get(&shape_key(a_si)) {
            Some(a_lsx) if !a_lsx.is_empty() => {}
            _ => a_lsi.push(a_si.clone()),
        }
    }
    if a_lsi.is_empty() {
        return Ok(());
    }

    // 5. Settle internal vertices and edges into solids.
    let bld = TopoBuilder::new();
    for a_sd in &a_lsd {
        let mut a_sd_cur = a_sd.clone();
        let mut i = 0;
        while i < a_lsi.len() {
            let mut a_si = a_lsi[i].clone();
            a_si.set_orientation(Orientation::Internal);
            let a_state = compute_state_by_one_point(&a_si, &a_sd_cur, 1e-11);
            if a_state != FaceState::In {
                i += 1;
                continue;
            }
            if a_msor.contains(&shape_key(&a_sd_cur)) {
                let mut a_sdx = Solid::new();
                for a_sh in iter_children(&a_sd_cur) {
                    bld.add(&mut a_sdx.0, &a_sh);
                }
                bld.add(&mut a_sdx.0, &a_si);
                f.history_mut().add_image(&a_sd_cur, a_sdx.0.clone());
                origins_append(f.origins_mut(), &a_sdx.0, a_sd_cur.clone());
                a_msor.remove(&shape_key(&a_sd_cur));
                a_sd_cur = a_sdx.0;
            } else {
                bld.add(&mut a_sd_cur, &a_si);
            }
            a_lsi.remove(i);
        }
    }
    Ok(())
}

/// `BOPAlgo_Tools::FillInternals` (`BOPAlgo_Tools.cxx:1751`).
///
/// Classifies `the_parts` (and their images) against `the_solids`. Vertices
/// and edges classified IN are added as INTERNAL immediately. IN faces are
/// grouped into connexity blocks and each block becomes an INTERNAL shell
/// (`_cxx:1867-1907`).
pub fn fill_internals_tools(
    the_solids: &mut [TopoShape],
    the_parts: &[TopoShape],
    the_images: &HashMap<usize, Vec<TopoShape>>,
    _the_ctx: &crate::int_tools_full::IntToolsContext,
) {
    if the_solids.is_empty() || the_parts.is_empty() {
        return;
    }
    let mut a_ms_solids: HashSet<usize> = HashSet::new();
    for a_solid in the_solids.iter() {
        if a_solid.shape_type() != ShapeType::Solid {
            continue;
        }
        for v in crate::topo_tools_full::vertices_of(a_solid) {
            a_ms_solids.insert(shape_key(&v.0));
        }
        for e in edges_of(a_solid) {
            a_ms_solids.insert(shape_key(&e.0));
        }
        for f in crate::topo_tools_full::faces_of(a_solid) {
            a_ms_solids.insert(shape_key(&f.0));
        }
    }

    let mut a_l_parts: Vec<TopoShape> = Vec::new();
    let mut a_l_input: Vec<TopoShape> = the_parts.to_vec();
    let mut ii = 0;
    while ii < a_l_input.len() {
        let a_part = a_l_input[ii].clone();
        ii += 1;
        match a_part.shape_type() {
            ShapeType::Vertex | ShapeType::Edge | ShapeType::Face => {
                if let Some(p_im) = the_images.get(&shape_key(&a_part)) {
                    for a_part_im in p_im {
                        if !a_ms_solids.contains(&shape_key(a_part_im)) {
                            a_l_parts.push(a_part_im.clone());
                        }
                    }
                } else if !a_ms_solids.contains(&shape_key(&a_part)) {
                    a_l_parts.push(a_part);
                }
            }
            _ => {
                for c in iter_children(&a_part) {
                    a_l_input.push(c);
                }
            }
        }
    }

    let mut an_in_faces: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    let bld = TopoBuilder::new();
    for a_sd in the_solids.iter_mut() {
        if a_sd.shape_type() != ShapeType::Solid {
            continue;
        }
        let mut keep: Vec<TopoShape> = Vec::new();
        for a_part in a_l_parts.drain(..) {
            let a_state = compute_state_by_one_point(&a_part, a_sd, occt_core::precision::CONFUSION);
            if a_state == FaceState::In {
                if a_part.shape_type() == ShapeType::Face {
                    an_in_faces
                        .entry(shape_key(a_sd))
                        .or_default()
                        .push(a_part);
                } else {
                    let mut p = a_part;
                    p.set_orientation(Orientation::Internal);
                    bld.add(a_sd, &p);
                }
            } else {
                keep.push(a_part);
            }
        }
        a_l_parts = keep;
    }

    for (k, faces) in an_in_faces {
        let Some(a_sd) = the_solids.iter_mut().find(|s| shape_key(s) == k) else {
            continue;
        };
        crate::bop_connexity_faces::add_internal_face_shells(a_sd, &faces);
    }
}

/// Image list of `s` copied into the FillInternals map form.
pub fn images_of_shape(
    h: &crate::bop_hist::BopHistory,
    s: &TopoShape,
) -> Option<Vec<TopoShape>> {
    h.image(s).map(|v| v.to_vec())
}

/// Build the FillInternals `theImages` map from a list of source shapes.
pub fn images_map_from_shapes(
    h: &crate::bop_hist::BopHistory,
    sources: &[TopoShape],
) -> HashMap<usize, Vec<TopoShape>> {
    let mut m = HashMap::new();
    for s in sources {
        if let Some(im) = h.image(s) {
            m.insert(shape_key(s), im.to_vec());
        }
    }
    m
}
