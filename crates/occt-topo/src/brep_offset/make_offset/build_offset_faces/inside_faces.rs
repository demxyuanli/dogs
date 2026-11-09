//! `RemoveInsideFaces` (`cxx:4291-4636`) and its two file-local helpers
//! `buildPairs` (`cxx:4073-4106`) and `buildIntersectionPairs`
//! (`cxx:4110-4287`).
//!
//! `RemoveInsideFaces` hands every group of face splits (plus the faces made
//! only of invalid edges) to `BOPAlgo_MakerVolume`; the volume builder's
//! solids say which splits bound the material and which are inside the offset
//! shell. Everything not on a solid boundary, or shared by two solids, is
//! dropped (`cxx:4461-4472`), solids made only of invalid splits are dropped
//! too (`cxx:4546-4605`), and the survivors are recorded in `mySolids`
//! (`cxx:4601-4603`). The removed edges drive `theMEInside` (`cxx:4616-4626`)
//! and the rebuilt intersection pairs (`cxx:4629-4635`).

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bop_maker_volume::BopMakerVolume;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::builder::TopoBuilder;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::statics::{find_shape, take_modified};
use super::BRepOffsetBuildOffsetFaces;

/// `static void buildPairs(const IndexedMap& theSMap, DataMap<Shape, Map>&
/// theIntPairs)` (`cxx:4073-4106`).
///
/// Binds an empty pair set for every shape, then adds every unordered pair
/// both ways (`cxx:4093-4105`).
pub(super) fn build_pairs(
    the_s_map: &ShapeIndexedMap,
    the_int_pairs: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
) {
    let a_nb_s = the_s_map.extent();
    // `cxx:4083-4086`.
    if a_nb_s < 2 {
        return;
    }
    // `cxx:4087-4093`.
    for it1 in 1..=a_nb_s {
        let Some(a_s) = the_s_map.value(it1).cloned() else {
            continue;
        };
        the_int_pairs.entry(shape_map_key(&a_s)).or_default();
    }
    // `cxx:4095-4105`.
    for it1 in 1..=a_nb_s {
        let Some(a_s1) = the_s_map.value(it1).cloned() else {
            continue;
        };
        let a_key1 = shape_map_key(&a_s1);
        for it2 in (it1 + 1)..=a_nb_s {
            let Some(a_s2) = the_s_map.value(it2).cloned() else {
                continue;
            };
            let a_key2 = shape_map_key(&a_s2);
            if let Some(a_map1) = the_int_pairs.get_mut(&a_key1) {
                a_map1.add(&a_s2);
            }
            if let Some(a_map2) = the_int_pairs.get_mut(&a_key2) {
                a_map2.add(&a_s1);
            }
        }
    }
}

/// `static void buildIntersectionPairs(const IndexedDataMap& myOFImages,
/// const IndexedDataMap& myInvalidFaces, const BOPAlgo_Builder& theBuilder,
/// const Map& theMFRemoved, const DataMap<Shape, Shape>& theFOrigins,
/// DataMap<Shape, DataMap<Shape, Map>>& theIntPairs)` (`cxx:4110-4287`).
///
/// For every invalid face the reference collects the splits still kept and the
/// splits already removed, splits them into connexity blocks and records the
/// pairs (in terms of the *original* faces) that still have to be intersected.
/// The kept splits are paired through the vertices they share with the volume
/// builder result (`cxx:4189-4227`); the removed splits only among the faces
/// that bound the same edge of the result (`cxx:4229-4283`).
#[allow(clippy::too_many_arguments)]
pub(super) fn build_intersection_pairs(
    my_of_images: &Dmvv,
    my_invalid_faces: &Dmvv,
    the_builder: &BopMakerVolume,
    the_mf_removed: &ShapeIndexedMap,
    the_f_origins: &HashMap<ShapeMapKey, TopoShape>,
    the_int_pairs: &mut HashMap<ShapeMapKey, HashMap<ShapeMapKey, ShapeIndexedMap>>,
) {
    // `cxx:4124-4130`: vertices of the builder's result -> faces.
    let a_dm_vf = map_shapes_and_ancestors(the_builder.shape(), ShapeType::Vertex, ShapeType::Face);

    // `cxx:4132-4137`.
    let an_images = the_builder.images();
    let an_origins = the_builder.origins();

    // `cxx:4140`.
    let a_nb_f = my_invalid_faces.extent();
    for i_f in 1..=a_nb_f {
        let Some(a_f_inv) = my_invalid_faces.keys().get(i_f - 1).cloned() else {
            continue;
        };

        // `cxx:4142-4148`: kept and removed splits.
        let mut a_cf_splits: Vec<TopoShape> = Vec::new();
        let mut a_cf_rem_splits: Vec<TopoShape> = Vec::new();

        for i_c in 0..2 {
            // `cxx:4150-4151`.
            let a_lf: Vec<TopoShape> = if i_c == 0 {
                my_invalid_faces.seek(&a_f_inv).map(|l| l.to_vec()).unwrap_or_default()
            } else {
                my_of_images.seek(&a_f_inv).map(|l| l.to_vec()).unwrap_or_default()
            };

            for it in &a_lf {
                // `cxx:4156-4157`.
                let mut a_lf_im: Vec<TopoShape> = Vec::new();
                take_modified(it, an_images, &mut a_lf_im, None);

                for it_im in &a_lf_im {
                    // `cxx:4160-4170`.
                    if the_mf_removed.contains(it_im) {
                        a_cf_rem_splits.push(it_im.clone());
                    } else {
                        a_cf_splits.push(it_im.clone());
                    }
                }
            }
        }

        // `cxx:4173-4180`.
        let a_cf = make_compound_from(a_cf_splits.iter());
        let a_lcb = AlgoTools::make_connexity_blocks_in(&a_cf, ShapeType::Edge, ShapeType::Face);
        if a_lcb.is_empty() {
            continue;
        }

        // `cxx:4182-4187`: `theIntPairs.Bound(aFInv, ...)`.
        let p_f_inter_map = the_int_pairs.entry(shape_map_key(&a_f_inv)).or_default();

        // `cxx:4190-4227`: pairs among the kept splits.
        for a_cb in &a_lcb {
            let mut a_mf_inter = ShapeIndexedMap::new();
            for a_cs in explorer_shapes(a_cb, ShapeType::Vertex) {
                let Some(p_lf_v) = a_dm_vf.get(&shape_map_key(&a_cs)) else {
                    continue;
                };
                for a_f_connected in p_lf_v.clone() {
                    // `cxx:4210-4218`.
                    let mut a_lf_or: Vec<TopoShape> = Vec::new();
                    take_modified(&a_f_connected, an_origins, &mut a_lf_or, None);
                    for it_or in &a_lf_or {
                        if let Some(p_f_or) = the_f_origins.get(&shape_map_key(it_or)) {
                            a_mf_inter.add(p_f_or);
                        }
                    }
                }
            }
            // `cxx:4221`.
            build_pairs(&a_mf_inter, p_f_inter_map);
        }

        // `cxx:4229-4231`: the same for the removed splits.
        let a_cf_rem = make_compound_from(a_cf_rem_splits.iter());
        let a_lcb_rem = AlgoTools::make_connexity_blocks_in(&a_cf_rem, ShapeType::Edge, ShapeType::Face);
        if a_lcb_rem.is_empty() {
            continue;
        }

        for a_cb in &a_lcb_rem {
            // `cxx:4237-4256`: edge -> faces of the result around the block.
            let mut a_dm_ef: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
            for a_cs in explorer_shapes(a_cb, ShapeType::Vertex) {
                let Some(p_lf_v) = a_dm_vf.get(&shape_map_key(&a_cs)) else {
                    continue;
                };
                for a_f_connected in p_lf_v.clone() {
                    crate::topo_tools::map_shapes_and_ancestors_into(
                        &a_f_connected,
                        ShapeType::Edge,
                        ShapeType::Face,
                        &mut a_dm_ef,
                    );
                }
            }

            // `cxx:4258-4284`: per edge of the result, pair the original faces
            // of the splits around it.
            let keys: Vec<ShapeMapKey> = a_dm_ef.keys().cloned().collect();
            for a_key in keys {
                let a_lf_connected = match a_dm_ef.get(&a_key) {
                    Some(l) => l.clone(),
                    None => continue,
                };
                if a_lf_connected.len() < 2 {
                    continue;
                }
                let mut a_mf_inter = ShapeIndexedMap::new();
                for a_f_connected in &a_lf_connected {
                    let mut a_lf_or: Vec<TopoShape> = Vec::new();
                    take_modified(a_f_connected, an_origins, &mut a_lf_or, None);
                    for it_or in &a_lf_or {
                        if let Some(p_f_or) = the_f_origins.get(&shape_map_key(it_or)) {
                            a_mf_inter.add(p_f_or);
                        }
                    }
                }
                build_pairs(&a_mf_inter, p_f_inter_map);
            }
        }
    }
}

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::RemoveInsideFaces(const List&
    /// theInvertedFaces, const IndexedMap& theMFToCheckInt, const IndexedMap&
    /// theMFInvInHole, const TopoDS_Shape& theFHoles, IndexedMap&
    /// theMERemoved, IndexedMap& theMEInside, const Message_ProgressRange&)`
    /// (`cxx:4291-4636`).
    ///
    /// The reference's `Message_ProgressScope` guards (`cxx:4313`, `:4409`,
    /// `:4550`) are the only use of `theRange`; the port has no progress
    /// reporting, so those arms are absent. Every other branch is kept.
    #[allow(clippy::too_many_arguments)]
    pub fn remove_inside_faces(
        &mut self,
        the_inverted_faces: &[TopoShape],
        the_mf_to_check_int: &ShapeIndexedMap,
        the_mf_inv_in_hole: &ShapeIndexedMap,
        the_f_holes: &TopoShape,
        the_me_removed: &mut ShapeIndexedMap,
        the_me_inside: &mut ShapeIndexedMap,
    ) {
        // `cxx:4306-4309`.
        let mut a_ls: Vec<TopoShape> = Vec::new();
        let mut a_mfence = ShapeIndexedMap::new();
        let mut a_mf_inv = ShapeIndexedMap::new();
        let mut a_dmf_imf: HashMap<ShapeMapKey, TopoShape> = HashMap::new();

        // `cxx:4311-4350`: one compound per offset face, holding its invalid
        // and valid splits.
        let a_nb = self.my_of_images.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let mut a_cf_imi_splits: Vec<TopoShape> = Vec::new();

            for j in 0..2 {
                // `cxx:4327-4330`.
                let a_splits: Vec<TopoShape> = if j == 0 {
                    self.my_invalid_faces.seek(&a_f).map(|l| l.to_vec()).unwrap_or_default()
                } else {
                    self.my_of_images.seek(&a_f).map(|l| l.to_vec()).unwrap_or_default()
                };

                for a_f_im in &a_splits {
                    // `cxx:4336-4346`.
                    if a_mfence.add(a_f_im) == a_mfence.extent() {
                        a_cf_imi_splits.push(a_f_im.clone());
                        a_dmf_imf.insert(shape_map_key(a_f_im), a_f.clone());
                        if j == 0 {
                            a_mf_inv.add(a_f_im);
                        }
                    }
                }
            }

            // `cxx:4349`.
            a_ls.push(make_compound_from(a_cf_imi_splits.iter()));
        }

        // `cxx:4352-4363`: the faces made only of invalid edges join the
        // arguments.
        let a_nb = the_mf_to_check_int.extent();
        for i in 1..=a_nb {
            let Some(a_f_sp) = the_mf_to_check_int.value(i).cloned() else {
                continue;
            };
            if a_mfence.add(&a_f_sp) == a_mfence.extent() {
                a_ls.push(a_f_sp);
            }
        }

        // `cxx:4365-4371`.
        let mut a_mv = BopMakerVolume::new();
        a_mv.set_arguments(&a_ls);
        a_mv.set_intersect(true);
        let _ = a_mv.perform();
        if a_mv.has_errors() {
            return;
        }

        // `cxx:4374-4378`: record the connections for the rebuilding.
        self.shapes_connections(&a_dmf_imf, &a_mv);

        // `cxx:4379-4390`.
        let a_sols = a_mv.shape().clone();
        let a_dm_fs = map_shapes_and_ancestors(&a_sols, ShapeType::Face, ShapeType::Solid);
        if a_dm_fs.is_empty() {
            return;
        }

        // `cxx:4392-4460`: when the volume builder dropped shapes, the solids
        // may be broken - every dropped split has to sit on a hole edge or be
        // connected to a solid.
        let mut a_mf_to_rem = ShapeIndexedMap::new();
        if a_mv.has_deleted() {
            // `cxx:4399-4407`.
            let mut a_me_holes = ShapeIndexedMap::new();
            for a_e in explorer_shapes(the_f_holes, ShapeType::Edge) {
                a_me_holes.add(&a_e);
            }
            let mut a_me_sols = ShapeIndexedMap::new();
            for a_e in explorer_shapes(&a_sols, ShapeType::Edge) {
                a_me_sols.add(&a_e);
            }

            // `cxx:4408-4459`.
            let a_nb = self.my_of_images.extent();
            for i in 1..=a_nb {
                let a_lf_im: Vec<TopoShape> = {
                    let Some(a_f) = self.my_of_images.keys().get(i - 1) else {
                        continue;
                    };
                    self.my_of_images.seek(a_f).map(|l| l.to_vec()).unwrap_or_default()
                };
                if a_lf_im.is_empty() {
                    continue;
                }

                let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                    continue;
                };
                let b_invalid = self.my_invalid_faces.seek(&a_f).is_some();
                let mut b_connected = false;
                let mut b_face_kept = false;

                // `cxx:4428-4452`.
                for a_f_im in &a_lf_im {
                    if !a_mv.is_deleted(a_f_im) {
                        b_face_kept = true;
                        continue;
                    }
                    for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                        if a_me_holes.contains(&a_e) {
                            b_face_kept = true;
                            a_mf_to_rem.add(a_f_im);
                            break;
                        }
                        if !b_face_kept && b_invalid && !b_connected {
                            b_connected = a_me_sols.contains(&a_e);
                        }
                    }
                }

                // `cxx:4455-4458`.
                if !b_face_kept && !b_connected {
                    return;
                }
            }
        }

        // `cxx:4461-4472`: faces on two solids, and the boundary edges of the
        // faces that sit on exactly one.
        let mut a_me_boundary = ShapeIndexedMap::new();
        let a_fs_order: Vec<ShapeMapKey> = a_dm_fs.keys().cloned().collect();
        for a_key in &a_fs_order {
            let a_l_sol = match a_dm_fs.get(a_key) {
                Some(l) => l,
                None => continue,
            };
            let a_f_im = match a_l_sol.first() {
                Some(f) => f,
                None => continue,
            };
            if a_l_sol.len() > 1 {
                a_mf_to_rem.add(a_f_im);
            } else if a_f_im.orientation() != Orientation::Internal {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    a_me_boundary.add(&a_e);
                }
            }
        }

        // `cxx:4475-4478`: the splits of the result solids.
        let a_mv_images = a_mv.images();

        // `cxx:4480-4487`: the invalid splits are replaced by their images.
        let a_nb = a_mf_inv.extent();
        for i in 1..=a_nb {
            let Some(a_f_inv) = a_mf_inv.value(i).cloned() else {
                continue;
            };
            take_modified(&a_f_inv, a_mv_images, &mut a_mf_inv, None);
        }

        // `cxx:4489-4494`: the faces invalidated by inverted edges.
        for it_lf in the_inverted_faces {
            take_modified(it_lf, a_mv_images, &mut a_mf_inv, None);
        }

        // `cxx:4496-4541`: an invalid split that fills a hole is only invalid
        // when its normal flipped.
        let a_nb_fh = the_mf_inv_in_hole.extent();
        for i in 1..=a_nb_fh {
            let Some(a_f_inv) = the_mf_inv_in_hole.value(i).cloned() else {
                continue;
            };
            let mut a_lf_inv_im = a_mv.modified(&a_f_inv);
            if a_lf_inv_im.is_empty() {
                a_lf_inv_im.push(a_f_inv.clone());
            }
            let Some(p_f_offset) = a_dmf_imf.get(&shape_map_key(&a_f_inv)).cloned() else {
                continue;
            };

            // `cxx:4516-4540`.
            for a_f_inv_im in &a_lf_inv_im {
                let Some(p_l_sols) = a_dm_fs.get(&shape_map_key(a_f_inv_im)) else {
                    continue;
                };
                if p_l_sols.len() != 1 {
                    continue;
                }
                let a_f_sol = p_l_sols[0].clone();
                let mut a_fx = TopoShape::new(ShapeType::Shape);
                if !find_shape(a_f_inv_im, &a_f_sol, None, &mut a_fx) {
                    continue;
                }
                // `cxx:4537-4540`.
                if crate::brep_offset::BRepOffsetTool::check_planes_normals(
                    &Face(a_fx),
                    &Face(p_f_offset.clone()),
                    1.0e-8,
                ) {
                    a_mf_to_rem.add(a_f_inv_im);
                }
            }
        }

        // `cxx:4543-4605`: a solid made only of invalid splits is dropped; the
        // splits already marked for removal are flipped back into `aMFKeep` so
        // the solid's own boundary survives.
        let b = TopoBuilder::new();
        let mut a_solids = b.make_compound_of(&[]);
        let mut a_mf_keep = ShapeIndexedMap::new();

        for a_sol in explorer_shapes(&a_sols, ShapeType::Solid) {
            let mut b_all_inv = true;
            let mut b_all_removed = true;

            for a_fs in explorer_shapes(&a_sol, ShapeType::Face) {
                // `cxx:4562-4566`.
                if a_fs.orientation() == Orientation::Internal {
                    a_mf_to_rem.add(&a_fs);
                    continue;
                }
                if a_mf_to_rem.contains(&a_fs) {
                    continue;
                }
                b_all_removed = false;
                b_all_inv &= a_mf_inv.contains(&a_fs);
            }

            if b_all_inv && !b_all_removed {
                // `cxx:4578-4595`: keep the splits already marked, remove the
                // rest.
                for a_fs in explorer_shapes(&a_sol, ShapeType::Face) {
                    if a_mf_to_rem.contains(&a_fs) {
                        if a_mf_keep.add(&a_fs) != a_mf_keep.extent() {
                            a_mf_keep.remove_key(&a_fs);
                        }
                    } else {
                        a_mf_to_rem.add(&a_fs);
                    }
                }
            } else {
                // `cxx:4597-4603`.
                let a_bb = TopoBuilder::new();
                a_bb.add(&mut a_solids.0, &a_sol);
                self.my_solids = Some(a_solids.0.clone());
            }
        }

        // `cxx:4607-4609`.
        let a_keep: Vec<TopoShape> = a_mf_keep.iter().cloned().collect();
        for a_fs in &a_keep {
            a_mf_to_rem.remove_key(a_fs);
        }

        // `cxx:4611-4614`: the isolated hanging parts.
        self.remove_hanging_parts(&a_mv, &a_dmf_imf, &a_mf_inv, &mut a_mf_to_rem);

        // `cxx:4616-4618`.
        self.remove_valid_splits(&a_mf_to_rem, &a_mv, the_me_removed);
        self.remove_invalid_splits(&a_mf_to_rem, &a_mv, the_me_removed);

        // `cxx:4619-4627`: the removed edges not on a solid boundary are
        // inside faces.
        the_me_inside.clear();
        let a_nb = the_me_removed.extent();
        for i in 1..=a_nb {
            let Some(a_e) = the_me_removed.value(i).cloned() else {
                continue;
            };
            if !a_me_boundary.contains(&a_e) {
                the_me_inside.add(&a_e);
            }
        }

        // `cxx:4629-4635`.
        if a_mf_to_rem.extent() != 0 {
            build_intersection_pairs(
                &self.my_of_images,
                &self.my_invalid_faces,
                &a_mv,
                &a_mf_to_rem,
                &a_dmf_imf,
                &mut self.my_intersection_pairs,
            );
        }
    }
}
