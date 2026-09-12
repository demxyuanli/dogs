//! `BOPAlgo_Builder::FillSameDomainFaces`.
//!
//! Source: `BOPAlgo_Builder_2.cxx:580-925`, plus the local `AddEdgeSet` at
//! `:562-576`. This replaces the previous geometric-signature shortcut in
//! [`crate::bop_build_faces::fill_same_domain_faces`]: OCCT never compares
//! vertex-grid edge signatures. It:
//!
//! 1. Returns immediately when there are no F/F interference records
//!    (`aNbFFs == 0`).
//! 2. Builds `aFaceToParent`: each source SOLID's faces, then propagates the
//!    parent onto each face's images. Two faces of one solid cannot be
//!    same-domain.
//! 3. Collects DS indices of faces that appear in an F/F record **and** have
//!    FaceInfo, sorts them, and groups each (image or original) by
//!    [`crate::bop_tools_set::BopToolsSet`] of EDGE sub-shapes.
//! 4. Planar faces whose DS box is closed on every side skip
//!    `AreFacesSameDomain` and are linked directly when they share an edge
//!    set. Other pairs go through [`crate::bop_pair_sd::PairOfShapeBoolean`].
//! 5. `BOPAlgo_Tools::FillMap` + `MakeBlocks` produce the SD groups. The
//!    representative is the original (DS-indexed) face with the smallest
//!    index, else the first face of the group. Every member is bound in
//!    `myShapesSD`. Original members are appended to their own image list.
//! 6. Image lists of source faces are rewritten to the SD representative and
//!    `myOrigins` is filled.

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::bop_build_faces::BopBuilderLike;
use crate::bop_occt_util::{has_face_info, nb_source, origins_append, shape_key};
use crate::bop_pair_sd::{perform_pair_vector, PairOfShapeBoolean};
use crate::bop_tools_set::BopToolsSet;
use crate::bopalgo_tools_blocks::{fill_map_shapes, make_blocks_shapes, IndexedShapeAdj};
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Face, TopoShape};
use crate::topo_tools_full::faces_of;

/// `AddEdgeSet` (`_2.cxx:562-576`): intern `theS` under its EDGE
/// [`BopToolsSet`] and keep the connection from that set to the shape.
fn add_edge_set(the_s: &TopoShape, the_map: &mut HashMap<BopToolsSet, Vec<TopoShape>>) {
    let mut a_se = BopToolsSet::new();
    a_se.add(the_s, ShapeType::Edge);
    the_map.entry(a_se).or_default().push(the_s.clone());
}

/// Parent solid of a face, keyed by TShape (`aFaceToParent`).
///
/// `_2.cxx:595-649`: every source SOLID binds each of its faces that is not
/// already bound; then each image of a bound parent is bound to the same
/// solid when the image is not yet a key.
fn face_to_parent_map<B: BopBuilderLike>(f: &B) -> HashMap<usize, TopoShape> {
    let ds = f.ds();
    let n_src = nb_source(ds);
    let mut a_face_to_parent: HashMap<usize, TopoShape> = HashMap::new();
    for i_src in 0..n_src {
        let Some(a_si) = ds.shape_info(i_src) else { continue };
        if a_si.shape_type() != ShapeType::Solid {
            continue;
        }
        let a_solid = a_si.shape().clone();
        for a_f in faces_of(&a_solid) {
            a_face_to_parent
                .entry(shape_key(&a_f.0))
                .or_insert_with(|| a_solid.clone());
        }
    }
    // `_2.cxx:619-648`: propagate the parent onto images that are not yet keys.
    let mut a_propagation: HashMap<usize, TopoShape> = HashMap::new();
    for i_src in 0..n_src {
        let Some(a_si) = ds.shape_info(i_src) else { continue };
        if a_si.shape_type() != ShapeType::Face {
            continue;
        }
        let a_f = a_si.shape();
        let Some(parent) = a_face_to_parent.get(&shape_key(a_f)).cloned() else {
            continue;
        };
        let Some(imgs) = f.history().image(a_f) else { continue };
        for piece in imgs {
            if !a_face_to_parent.contains_key(&shape_key(piece)) {
                a_propagation.insert(shape_key(piece), parent.clone());
            }
        }
    }
    for (k, v) in a_propagation {
        a_face_to_parent.insert(k, v);
    }
    a_face_to_parent
}

/// Whether the DS box of face `n_f` is closed on every side
/// (`_2.cxx:713-716`). A missing box is treated as open: the planar
/// shortcut is then not taken and `AreFacesSameDomain` decides.
fn box_is_bounded(ds: &crate::bopds::BopdsDS, n_f: usize) -> bool {
    let Some(a_box) = ds.box_of(n_f) else {
        return false;
    };
    !(a_box.is_open_xmin()
        || a_box.is_open_xmax()
        || a_box.is_open_ymin()
        || a_box.is_open_ymax()
        || a_box.is_open_zmin()
        || a_box.is_open_zmax())
}

/// `myContext->SurfaceAdaptor(TopoDS::Face(aF)).GetType() == GeomAbs_Plane`.
fn face_is_plane(a_f: &TopoShape) -> bool {
    let Some(surf) = BRepTool::face_surface_world(&Face(a_f.clone())) else {
        return false;
    };
    classify_surface(surf.as_ref()) == SurfaceKind::Plane
}

/// Collect sorted unique DS indices of faces that participate in an F/F
/// record and have FaceInfo (`_2.cxx:658-687`).
fn face_indices_from_ff<B: BopBuilderLike>(f: &B) -> Vec<usize> {
    let ds = f.ds();
    let mut a_m_fence: HashSet<usize> = HashSet::new();
    let mut a_fi_vec: Vec<usize> = Vec::new();
    for a_ff in ds.interf_ff() {
        let (n0, n1) = a_ff.indices();
        for n_f in [n0, n1] {
            if !has_face_info(ds, n_f) {
                continue;
            }
            if !a_m_fence.insert(n_f) {
                continue;
            }
            a_fi_vec.push(n_f);
        }
    }
    a_fi_vec.sort_unstable();
    a_fi_vec
}

/// Group each participating face (or its images) by EDGE set
/// (`_2.cxx:696-741`). Planar bounded faces are remembered in `a_mf_planar`.
fn collect_edge_set_faces<B: BopBuilderLike>(
    f: &B,
    a_fi_vec: &[usize],
) -> (HashMap<BopToolsSet, Vec<TopoShape>>, HashSet<usize>) {
    let ds = f.ds();
    let mut an_e_set_faces: HashMap<BopToolsSet, Vec<TopoShape>> = HashMap::new();
    let mut a_mf_planar: HashSet<usize> = HashSet::new();
    for &n_f in a_fi_vec {
        let Some(a_si) = ds.shape_info(n_f) else { continue };
        let a_f = a_si.shape();
        let mut b_check_planar = false;
        if face_is_plane(a_f) {
            b_check_planar = box_is_bounded(ds, n_f);
        }
        if let Some(p_lf_sp) = f.history().image(a_f) {
            for a_it_lf in p_lf_sp {
                add_edge_set(a_it_lf, &mut an_e_set_faces);
                if b_check_planar {
                    a_mf_planar.insert(shape_key(a_it_lf));
                }
            }
        } else {
            add_edge_set(a_f, &mut an_e_set_faces);
            if b_check_planar {
                a_mf_planar.insert(shape_key(a_f));
            }
        }
    }
    (an_e_set_faces, a_mf_planar)
}

/// Pair every two faces that share an EDGE set (`_2.cxx:750-793`).
///
/// Same-parent solids are skipped. Two planar-bounded faces are linked
/// immediately. Remaining pairs are queued for `AreFacesSameDomain`.
fn collect_sd_pairs(
    an_e_set_faces: &HashMap<BopToolsSet, Vec<TopoShape>>,
    a_mf_planar: &HashSet<usize>,
    a_face_to_parent: &HashMap<usize, TopoShape>,
    fuzzy: f64,
    a_dmsls: &mut IndexedShapeAdj,
) -> Vec<PairOfShapeBoolean> {
    let mut a_vpsb: Vec<PairOfShapeBoolean> = Vec::new();
    for a_lf in an_e_set_faces.values() {
        if a_lf.len() < 2 {
            continue;
        }
        for i1 in 0..a_lf.len() {
            let a_f1 = &a_lf[i1];
            let b_check_planar = a_mf_planar.contains(&shape_key(a_f1));
            let p_parent1 = a_face_to_parent.get(&shape_key(a_f1));
            for a_f2 in a_lf.iter().skip(i1 + 1) {
                let p_parent2 = a_face_to_parent.get(&shape_key(a_f2));
                if let (Some(p1), Some(p2)) = (p_parent1, p_parent2) {
                    if p1.same_tshape(p2) {
                        continue;
                    }
                }
                if b_check_planar && a_mf_planar.contains(&shape_key(a_f2)) {
                    fill_map_shapes(a_f1, a_f2, a_dmsls);
                    continue;
                }
                a_vpsb.push(PairOfShapeBoolean::from_faces(
                    a_f1.clone(),
                    a_f2.clone(),
                    fuzzy,
                ));
            }
        }
    }
    a_vpsb
}

/// Representative of an SD group (`_2.cxx:837-873`): the original DS face
/// with the smallest index, else the first member.
fn group_representative<B: BopBuilderLike>(f: &B, a_lsd: &[TopoShape]) -> TopoShape {
    let mut p_fsd: Option<TopoShape> = None;
    let mut n_f_min = usize::MAX;
    for a_f in a_lsd {
        if let Some(n_f) = f.ds().index(a_f) {
            if n_f < n_f_min {
                n_f_min = n_f;
                p_fsd = Some(a_f.clone());
            }
        }
    }
    p_fsd.unwrap_or_else(|| a_lsd[0].clone())
}

/// Bind SD connections and mark original faces as split into themselves
/// (`_2.cxx:829-882`).
fn bind_sd_groups<B: BopBuilderLike>(f: &mut B, a_m_blocks: &[Vec<TopoShape>]) {
    for a_lsd in a_m_blocks {
        if a_lsd.is_empty() {
            continue;
        }
        let p_fsd = group_representative(f, a_lsd);
        for a_f in a_lsd {
            if f.ds().index(a_f).is_some() {
                // `_2.cxx:858`: Bound + Append the original onto its image list.
                f.history_mut().add_image(a_f, a_f.clone());
            }
        }
        for a_f in a_lsd {
            f.bind_shapes_sd(a_f.clone(), p_fsd.clone());
        }
    }
}

/// Rewrite source-face images to the SD representative and fill origins
/// (`_2.cxx:884-921`).
fn update_images_and_origins<B: BopBuilderLike>(f: &mut B) {
    let n = nb_source(f.ds());
    let mut jobs: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
    for i in 0..n {
        let Some(a_si) = f.ds().shape_info(i) else { continue };
        if a_si.shape_type() != ShapeType::Face {
            continue;
        }
        let a_f = a_si.shape().clone();
        let Some(p_lf_im) = f.history().image(&a_f) else { continue };
        let mut updated: Vec<TopoShape> = Vec::new();
        for a_f_im in p_lf_im {
            let resolved = f.seek_shapes_sd(a_f_im).unwrap_or_else(|| a_f_im.clone());
            updated.push(resolved);
        }
        jobs.push((a_f, updated));
    }
    for (a_f, updated) in jobs {
        for a_f_im in &updated {
            origins_append(f.origins_mut(), a_f_im, a_f.clone());
        }
        if let Some(list) = f.history_mut().image_list_mut(&a_f) {
            *list = updated;
        }
    }
}

/// `BOPAlgo_Builder::FillSameDomainFaces` (`_2.cxx:580-925`).
pub fn fill_same_domain_faces_occt<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    if f.ds().interf_ff().is_empty() {
        return Ok(());
    }
    let a_face_to_parent = face_to_parent_map(f);
    let a_fi_vec = face_indices_from_ff(f);
    let (an_e_set_faces, a_mf_planar) = collect_edge_set_faces(f, &a_fi_vec);
    let mut a_dmsls = IndexedShapeAdj::new();
    let mut a_vpsb = collect_sd_pairs(
        &an_e_set_faces,
        &a_mf_planar,
        &a_face_to_parent,
        f.fuzzy_value(),
        &mut a_dmsls,
    );
    let mut ctx = IntToolsContext::new();
    perform_pair_vector(&mut a_vpsb, &mut ctx);
    for a_psb in &a_vpsb {
        if a_psb.flag() {
            fill_map_shapes(a_psb.shape1(), a_psb.shape2(), &mut a_dmsls);
        }
    }
    let a_m_blocks = make_blocks_shapes(&a_dmsls);
    bind_sd_groups(f, &a_m_blocks);
    update_images_and_origins(f);
    Ok(())
}
