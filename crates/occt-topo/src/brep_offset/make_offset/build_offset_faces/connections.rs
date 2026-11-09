//! `ShapesConnections` (`cxx:4638-4964`): the face/face intersections the
//! volume builder found that the offset algorithm must remember.
//!
//! The method walks the interference records of `BOPAlgo_MakerVolume` and
//! records, per original offset face, the faces it has to be intersected with
//! during the rebuilding:
//!
//! * `mySSInterfs` (`cxx:4817-4837`) - for a face/face intersection involving
//!   at least one invalid face, when a split of one face carries both an
//!   invalid edge and an intersection edge (`cxx:4779-4803`), the opposite
//!   face is appended for every invalid edge of that split and for the offset
//!   face itself;
//! * `mySSInterfsArt` (`cxx:4946-4959`) - for the artificial-invalid case, the
//!   face/face and edge/face connections plus the intersection vertices, so the
//!   rebuilding can trim the intersection edge.
//!
//! The file-local `AddToContainer` overloads are the port's
//! [`super::statics::AddToContainer`] (`Vec` for lists, `ShapeIndexedMap` for
//! maps).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bop_maker_volume::BopMakerVolume;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{children_keep_orientation, explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::statics::{add_to_container_of_map, append_to_list};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::ShapesConnections(const
    /// DataMap<Shape, Shape>& theDMFOr, BOPAlgo_Builder& theBuilder)`
    /// (`cxx:4638-4964`).
    ///
    /// `the_dmf_or` is `aDMFImF`: split -> the offset face it came from. The
    /// builder is read only (`Modified`, `PDS`), so it is taken by reference.
    pub fn shapes_connections(
        &mut self,
        the_dmf_or: &HashMap<ShapeMapKey, TopoShape>,
        the_builder: &BopMakerVolume,
    ) {
        // `cxx:4648-4655`: the connexity blocks of the invalid edges, so a whole
        // block is used instead of the single edge.
        let mut a_ce_inv_edges: Vec<TopoShape> = Vec::new();
        let a_nb = self.my_invalid_edges.extent();
        for i in 1..=a_nb {
            if let Some(a_e) = self.my_invalid_edges.value(i) {
                a_ce_inv_edges.push(a_e.clone());
            }
        }
        let a_ce_inv = make_compound_from(a_ce_inv_edges.iter());

        let a_lcb = AlgoTools::make_connexity_blocks_in(&a_ce_inv, ShapeType::Vertex, ShapeType::Edge);

        // `cxx:4659-4668`: edge -> the block it belongs to.
        let mut a_ecb_map: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
        for a_cb in &a_lcb {
            // `cxx:4663-4667`: `TopoDS_Iterator itE(itCB.Value())`.
            for it_e in children_keep_orientation(a_cb) {
                a_ecb_map.insert(shape_map_key(&it_e), a_cb.clone());
            }
        }

        // `cxx:4670-4688`: invalid edge -> its images, or itself when the
        // builder did not touch it.
        let mut a_dme_or: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let a_nb = self.my_invalid_edges.extent();
        for i in 1..=a_nb {
            let Some(a_e_inv) = self.my_invalid_edges.value(i).cloned() else {
                continue;
            };
            let a_le_im = the_builder.modified(&a_e_inv);
            if a_le_im.is_empty() {
                // `cxx:4678-4680`.
                add_to_container_of_map(&mut a_dme_or, &a_e_inv, &a_e_inv);
                continue;
            }
            for a_e_im in &a_le_im {
                add_to_container_of_map(&mut a_dme_or, a_e_im, &a_e_inv);
            }
        }

        // `cxx:4690-4695`.
        let p_ds = the_builder.ds();
        let a_ffs = p_ds.interf_ff();

        // `cxx:4697-4962`: every face/face interference.
        for a_ff in a_ffs {
            let a_vnc = a_ff.curves();
            // `cxx:4700-4703`.
            if a_vnc.is_empty() {
                continue;
            }

            // `cxx:4706-4717`.
            let Some(a_f_im1) = p_ds.shape(a_ff.index1()).cloned() else {
                continue;
            };
            let Some(a_f_im2) = p_ds.shape(a_ff.index2()).cloned() else {
                continue;
            };
            let Some(p_f1) = the_dmf_or.get(&shape_map_key(&a_f_im1)).cloned() else {
                continue;
            };
            let Some(p_f2) = the_dmf_or.get(&shape_map_key(&a_f_im2)).cloned() else {
                continue;
            };
            if is_same(&p_f1, &p_f2) {
                continue;
            }

            // `cxx:4722-4728`.
            let b_inv1 = self.my_invalid_faces.seek(&p_f1).is_some();
            let b_inv2 = self.my_invalid_faces.seek(&p_f2).is_some();
            if !b_inv1 && !b_inv2 {
                continue;
            }

            // `cxx:4730-4750`: the intersection edges of the interference.
            let mut a_me_int = ShapeIndexedMap::new();
            for a_nc in a_vnc {
                for a_pb in a_nc.pave_blocks() {
                    if a_pb.has_edge() {
                        if let Some(a_e_int) = p_ds.shape(a_pb.edge()) {
                            a_me_int.add(a_e_int);
                        }
                    }
                }
            }
            if a_me_int.extent() == 0 {
                continue;
            }

            // `cxx:4753-4838`: check whether the invalid edges of the face sit
            // in the same split as the intersection edges.
            for i in 0..2usize {
                if (i == 0 && !b_inv1) || (i == 1 && !b_inv2) {
                    continue;
                }
                let a_f = if i == 0 { p_f1.clone() } else { p_f2.clone() };
                let a_f_op = if i == 0 { p_f2.clone() } else { p_f1.clone() };
                let a_f_im = if i == 0 { a_f_im1.clone() } else { a_f_im2.clone() };

                let mut b_found = false;

                // `cxx:4767-4771`.
                let mut a_lf_im = the_builder.modified(&a_f_im);
                if a_lf_im.is_empty() {
                    a_lf_im.push(a_f_im.clone());
                }

                for a_f_im_im in &a_lf_im {
                    // `cxx:4778-4803`.
                    let mut b_inv = false;
                    let mut b_int = false;
                    for a_e in explorer_shapes(a_f_im_im, ShapeType::Edge) {
                        if !b_inv {
                            b_inv = a_dme_or.contains_key(&shape_map_key(&a_e));
                        }
                        if !b_int {
                            b_int = a_me_int.contains(&a_e);
                        }
                        if b_inv && b_int {
                            break;
                        }
                    }
                    if !b_int || !b_inv {
                        continue;
                    }
                    b_found = true;

                    // `cxx:4804-4827`: the opposite face goes with every invalid
                    // edge of the split.
                    for a_e in explorer_shapes(a_f_im_im, ShapeType::Edge) {
                        let Some(p_le_or) = a_dme_or.get(&shape_map_key(&a_e)).cloned() else {
                            continue;
                        };
                        for a_e_or in &p_le_or {
                            append_to_list(
                                self.my_ss_interfs
                                    .entry(shape_map_key(a_e_or))
                                    .or_default(),
                                &a_f_op,
                            );
                        }
                    }
                }

                // `cxx:4829-4838`.
                if b_found {
                    append_to_list(
                        self.my_ss_interfs.entry(shape_map_key(&a_f)).or_default(),
                        &a_f_op,
                    );
                }
            }

            // `cxx:4840-4961`: the artificial case.
            for i_f in 0..2usize {
                let a_f_art = if i_f == 0 { p_f1.clone() } else { p_f2.clone() };
                let a_f_opposite = if i_f == 0 { p_f2.clone() } else { p_f1.clone() };

                // `cxx:4847-4853`.
                if !self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f_art)) {
                    continue;
                }
                if self.my_invalid_faces.seek(&a_f_opposite).is_some()
                    && !self
                        .my_art_invalid_faces
                        .contains_key(&shape_map_key(&a_f_opposite))
                {
                    continue;
                }

                // `cxx:4857-4882`: the invalid edges to avoid, taken from the
                // whole block of every own invalid edge of the opposite face.
                let mut a_me_avoid = ShapeIndexedMap::new();
                if let Some(p_fe_inv) = self
                    .my_art_invalid_faces
                    .get(&shape_map_key(&a_f_opposite))
                {
                    let a_fe_inv = p_fe_inv.clone();
                    let a_nb_e = a_fe_inv.extent();
                    for i_e in 1..=a_nb_e {
                        let Some(a_e) = a_fe_inv.value(i_e) else {
                            continue;
                        };
                        if let Some(p_cb) = a_ecb_map.get(&shape_map_key(a_e)) {
                            for a_ee in explorer_shapes(p_cb, ShapeType::Edge) {
                                a_me_avoid.add(&a_ee);
                            }
                        }
                    }
                } else if let Some(p_lf_im) = self
                    .my_of_images
                    .seek(&a_f_opposite)
                    .map(|l| l.to_vec())
                {
                    for it_lf_im in &p_lf_im {
                        for exp_e in explorer_shapes(it_lf_im, ShapeType::Edge) {
                            if let Some(p_cb) = a_ecb_map.get(&shape_map_key(&exp_e)) {
                                for a_ee in explorer_shapes(p_cb, ShapeType::Edge) {
                                    a_me_avoid.add(&a_ee);
                                }
                            }
                        }
                    }
                }

                let a_f_art_im = if i_f == 0 { a_f_im1.clone() } else { a_f_im2.clone() };
                let a_f_opposite_im = if i_f == 0 { a_f_im2.clone() } else { a_f_im1.clone() };

                // `cxx:4889`.
                let Some(n_f_op) = p_ds.index(&a_f_opposite_im) else {
                    continue;
                };

                // `cxx:4891-4959`.
                for a_e in explorer_shapes(&a_f_art_im, ShapeType::Edge) {
                    if !self.my_invalid_edges.contains(&a_e)
                        || self.my_inverted_edges.contains(&a_e)
                        || a_me_avoid.contains(&a_e)
                    {
                        continue;
                    }
                    let Some(n_e) = p_ds.index(&a_e) else {
                        continue;
                    };
                    // `cxx:4903-4906`.
                    if !p_ds.has_interf_pair(n_e, n_f_op) {
                        continue;
                    }

                    // `cxx:4908-4922`: the vertex parts of the edge/face
                    // interferences.
                    let mut a_lv: Vec<TopoShape> = Vec::new();
                    for a_ef in p_ds.interf_ef() {
                        if a_ef.contains(n_e) && a_ef.contains(n_f_op) {
                            let is_vertex_part = a_ef
                                .common_part()
                                .map(|cp| cp.is_vertex_part())
                                .unwrap_or(false);
                            if is_vertex_part {
                                if let Some(n_new) = a_ef.new_index() {
                                    if let Some(a_v) = p_ds.shape(n_new) {
                                        a_lv.push(a_v.clone());
                                    }
                                }
                            }
                        }
                    }
                    if a_lv.is_empty() {
                        continue;
                    }

                    // `cxx:4930-4944`: an opposite intersection must exist.
                    let Some(n_f_art) = p_ds.index(&a_f_art_im) else {
                        continue;
                    };
                    let mut b_opposite_found = false;
                    for a_e_op in explorer_shapes(&a_f_opposite_im, ShapeType::Edge) {
                        let Some(n_e_op) = p_ds.index(&a_e_op) else {
                            continue;
                        };
                        if p_ds.has_interf_pair(n_e_op, n_f_art) {
                            b_opposite_found = true;
                            break;
                        }
                    }
                    if !b_opposite_found {
                        continue;
                    }

                    // `cxx:4946-4959`.
                    add_to_container_of_map(&mut self.my_ss_interfs_art, &a_f_art, &a_f_opposite);
                    add_to_container_of_map(&mut self.my_ss_interfs_art, &a_e, &a_f_opposite);
                    for a_v in &a_lv {
                        add_to_container_of_map(&mut self.my_ss_interfs_art, &a_e, a_v);
                    }
                }
            }
        }
    }
}
