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
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::bop_occt_util::{
    images_bound, images_of, iter_children, nb_source, origins_append, own_internal_shapes,
    shape_key, treat_compound, BopSolidHost,
};
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Solid, TopoShape};
use crate::topo_tools_full::{edges_of, map_shapes_and_ancestors_into, shape_map_key, ShapeMapKey};

/// `BOPTools_AlgoTools::ComputeStateByOnePoint` stand-in used by
/// FillInternalShapes (`_3.cxx:836`, tolerance `1.e-11`).
fn compute_state_by_one_point(shape: &TopoShape, solid: &TopoShape, tol: f64) -> FaceState {
    crate::bopalgo_tools_class::compute_state_by_one_point(shape, solid, tol)
}

/// `BOPAlgo_Builder::FillInternalShapes` (`_3.cxx:622`).
pub fn fill_internal_shapes_occt<B: BopSolidHost>(f: &mut B) -> Result<(), String> {
    // 1. Shapes to process from the pure arguments.
    let mut lsc: Vec<TopoShape> = Vec::new();
    let mut fence: HashSet<ShapeMapKey> = HashSet::new();
    for a in f.arguments() {
        treat_compound(a, &mut lsc, &mut fence);
    }
    let mut largs: Vec<TopoShape> = Vec::new();
    fence.clear();
    for s in &lsc {
        match s.shape_type() {
            ShapeType::Wire => {
                for e in edges_of(s) {
                    if fence.insert(shape_key(&e.0)) {
                        largs.push(e.0);
                    }
                }
            }
            ShapeType::Vertex | ShapeType::Edge => largs.push(s.clone()),
            _ => {}
        }
    }
    fence.clear();
    let mut msi: Vec<TopoShape> = Vec::new();
    for s in &largs {
        let t = s.shape_type();
        if t != ShapeType::Vertex && t != ShapeType::Edge && t != ShapeType::Wire {
            continue;
        }
        if !fence.insert(shape_key(s)) {
            continue;
        }
        if images_bound(f.history(), s) {
            for im in images_of(f.history(), s) {
                msi.push(im.clone());
            }
        } else {
            msi.push(s.clone());
        }
    }

    // 2. Internal vertices/edges from source solids + ancestor map of splits.
    fence.clear();
    let mut lsd: Vec<TopoShape> = Vec::new();
    let mut msx: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut msor: HashSet<ShapeMapKey> = HashSet::new();
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
                    msi.push(sp.clone());
                }
            } else {
                msi.push(a_si);
            }
        }
        if images_bound(f.history(), &a_s) {
            for a_sp in images_of(f.history(), &a_s) {
                if fence.insert(shape_key(a_sp)) {
                    map_shapes_and_ancestors_into(a_sp, ShapeType::Vertex, ShapeType::Edge, &mut msx);
                    map_shapes_and_ancestors_into(a_sp, ShapeType::Vertex, ShapeType::Face, &mut msx);
                    map_shapes_and_ancestors_into(a_sp, ShapeType::Edge, ShapeType::Face, &mut msx);
                    lsd.push(a_sp.clone());
                }
            }
        } else if fence.insert(shape_key(&a_s)) {
            map_shapes_and_ancestors_into(&a_s, ShapeType::Vertex, ShapeType::Edge, &mut msx);
            map_shapes_and_ancestors_into(&a_s, ShapeType::Vertex, ShapeType::Face, &mut msx);
            map_shapes_and_ancestors_into(&a_s, ShapeType::Edge, ShapeType::Face, &mut msx);
            lsd.push(a_s.clone());
            msor.insert(shape_key(&a_s));
        }
    }

    // 3. Drop candidates already tied to a split-solid ancestor.
    let mut lsi: Vec<TopoShape> = Vec::new();
    for a_si in &msi {
        match msx.get(&shape_map_key(a_si)) {
            Some(lsx) if !lsx.is_empty() => {}
            _ => lsi.push(a_si.clone()),
        }
    }
    if lsi.is_empty() {
        return Ok(());
    }

    // 5. Settle internal vertices and edges into solids.
    let bld = TopoBuilder::new();
    for sd in &lsd {
        let mut sd_cur = sd.clone();
        let mut i = 0;
        while i < lsi.len() {
            let mut a_si = lsi[i].clone();
            a_si.set_orientation(Orientation::Internal);
            let state = compute_state_by_one_point(&a_si, &sd_cur, 1e-11);
            if state != FaceState::In {
                i += 1;
                continue;
            }
            if msor.contains(&shape_key(&sd_cur)) {
                let mut sdx = Solid::new();
                for sh in iter_children(&sd_cur) {
                    bld.add(&mut sdx.0, &sh);
                }
                bld.add(&mut sdx.0, &a_si);
                f.history_mut().add_image(&sd_cur, sdx.0.clone());
                origins_append(f.origins_mut(), &sdx.0, sd_cur.clone());
                msor.remove(&shape_key(&sd_cur));
                sd_cur = sdx.0;
            } else {
                bld.add(&mut sd_cur, &a_si);
            }
            lsi.remove(i);
        }
    }
    Ok(())
}

/// `BOPAlgo_Tools::FillInternals` (`BOPAlgo_Tools.cxx:1751`).
///
/// Classifies `parts` (and their images) against `solids`. Vertices
/// and edges classified IN are added as INTERNAL immediately. IN faces are
/// grouped into connexity blocks and each block becomes an INTERNAL shell
/// (`_cxx:1867-1907`).
pub fn fill_internals_tools(
    solids: &mut [TopoShape],
    parts: &[TopoShape],
    images: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    _the_ctx: &crate::int_tools_full::IntToolsContext,
) {
    if solids.is_empty() || parts.is_empty() {
        return;
    }
    let mut ms_solids: HashSet<ShapeMapKey> = HashSet::new();
    for solid in solids.iter() {
        if solid.shape_type() != ShapeType::Solid {
            continue;
        }
        for v in crate::topo_tools_full::vertices_of(solid) {
            ms_solids.insert(shape_key(&v.0));
        }
        for e in edges_of(solid) {
            ms_solids.insert(shape_key(&e.0));
        }
        for f in crate::topo_tools_full::faces_of(solid) {
            ms_solids.insert(shape_key(&f.0));
        }
    }

    let mut l_parts: Vec<TopoShape> = Vec::new();
    let mut l_input: Vec<TopoShape> = parts.to_vec();
    let mut ii = 0;
    while ii < l_input.len() {
        let part = l_input[ii].clone();
        ii += 1;
        match part.shape_type() {
            ShapeType::Vertex | ShapeType::Edge | ShapeType::Face => {
                if let Some(p_im) = images.get(&shape_key(&part)) {
                    for part_im in p_im {
                        if !ms_solids.contains(&shape_key(part_im)) {
                            l_parts.push(part_im.clone());
                        }
                    }
                } else if !ms_solids.contains(&shape_key(&part)) {
                    l_parts.push(part);
                }
            }
            _ => {
                for c in iter_children(&part) {
                    l_input.push(c);
                }
            }
        }
    }

    let mut an_in_faces: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let bld = TopoBuilder::new();
    for sd in solids.iter_mut() {
        if sd.shape_type() != ShapeType::Solid {
            continue;
        }
        let mut keep: Vec<TopoShape> = Vec::new();
        for part in l_parts.drain(..) {
            let state = compute_state_by_one_point(&part, sd, occt_core::precision::CONFUSION);
            if state == FaceState::In {
                if part.shape_type() == ShapeType::Face {
                    an_in_faces
                        .entry(shape_key(sd))
                        .or_default()
                        .push(part);
                } else {
                    let mut p = part;
                    p.set_orientation(Orientation::Internal);
                    bld.add(sd, &p);
                }
            } else {
                keep.push(part);
            }
        }
        l_parts = keep;
    }

    for (k, faces) in an_in_faces {
        let Some(sd) = solids.iter_mut().find(|s| shape_key(s) == k) else {
            continue;
        };
        crate::bop_connexity_faces::add_internal_face_shells(sd, &faces);
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
) -> HashMap<ShapeMapKey, Vec<TopoShape>> {
    let mut m = HashMap::new();
    for s in sources {
        if let Some(im) = h.image(s) {
            m.insert(shape_key(s), im.to_vec());
        }
    }
    m
}
