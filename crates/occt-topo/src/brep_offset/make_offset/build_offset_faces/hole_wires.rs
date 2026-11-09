//! Hole-wire classification of `BRepOffset_BuildOffsetFaces` -
//! `FindFacesInsideHoleWires` (`cxx:3089-3310`).
//!
//! While offsetting an original face with holes, the offset of the face has no
//! holes: the hole wires are simply absent from its splits. The splits that
//! *fill* a hole therefore have to be recognized and removed, otherwise the
//! offset solid would be built with the holes filled.
//!
//! The method has three phases, all of them kept here in the reference order:
//!
//! 1. build the "new holes" of the original face - the images of its hole
//!    wires assembled into planar faces (`cxx:3125-3234`);
//! 2. classify every split of the offset face by a point inside it against the
//!    new holes and record the ones that fall inside (`cxx:3244-3272`);
//! 3. drop the holes whose boundary is connected to splits of another offset
//!    face, because removing them would leave an uncovered edge
//!    (`cxx:3276-3308`).

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools3d::point_in_face;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_tools::outer_wire;
use crate::brep_tool::BRepTool;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::bop_builder_face_perform::perform_face_builder;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, is_same, map_shapes_and_ancestors,
    map_shapes_and_ancestors_into, shape_map_key, wires_of_face, ShapeMapKey,
};

use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::FindFacesInsideHoleWires(const
    /// TopoDS_Face& theFOrigin, const TopoDS_Face& theFOffset, const List&
    /// theLFImages, const DataMap<Shape, List>& theDMEOrLEIm, const
    /// IndexedDataMap<Shape, List>& theEFMap, Map& theMFHoles)`
    /// (`cxx:3089-3310`).
    ///
    /// * `the_dme_or_le_im` is the reference's original edge -> its images map
    ///   (`theDMEOrLEIm`), carried by `HashMap<ShapeMapKey, Vec<TopoShape>>`.
    /// * `the_ef_map` is the reference's `IndexedDataMap` edge -> faces over
    ///   every split of every offset face; only `FindFromKey` is used on it
    ///   (`cxx:3304`), so the key order is not observable.
    /// * `the_mf_holes` is the `theMFHoles` map, filled and then pruned in
    ///   place.
    pub fn find_faces_inside_hole_wires(
        &mut self,
        the_f_origin: &Face,
        the_f_offset: &Face,
        the_lf_images: &[TopoShape],
        the_dme_or_le_im: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_ef_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_mf_holes: &mut ShapeIndexedMap,
    ) {
        // `cxx:3100-3103`.
        if the_lf_images.is_empty() {
            return;
        }

        // `cxx:3105-3117`: the hole wires of the original face.
        let mut a_l_hole_wires: Vec<TopoShape> = Vec::new();
        let an_outer_wire = outer_wire(the_f_origin);
        for a_hole_wire in wires_of_face(the_f_origin) {
            let is_outer = match &an_outer_wire {
                Some(w) => is_same(&a_hole_wire.0, &w.0),
                None => false,
            };
            if !is_outer && a_hole_wire.0.orientation() != Orientation::Internal {
                a_l_hole_wires.push(a_hole_wire.0);
            }
        }
        if a_l_hole_wires.is_empty() {
            return;
        }

        // `cxx:3124-3130`: `ChangeSeek` then `Bound` - the entry is created
        // empty when absent, and stays bound even if no hole face is found.
        let a_f_origin_key = shape_map_key(&the_f_origin.0);
        let mut a_lf_new_holes: Vec<TopoShape> = self
            .my_f_new_holes
            .get(&a_f_origin_key)
            .cloned()
            .unwrap_or_default();

        if a_lf_new_holes.is_empty() {
            // `cxx:3136-3143`: the edges of every split of the offset face.
            let mut a_me_splits = ShapeIndexedMap::new();
            for a_f_im in the_lf_images {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    a_me_splits.add(&a_e);
                }
            }

            // `cxx:3146-3233`: one planar face per original hole wire whose
            // images are all present among the splits.
            for a_hole_wire in &a_l_hole_wires {
                // `cxx:3151-3170`: the images of the wire's edges.
                let mut a_me_im_wire = ShapeIndexedMap::new();
                for a_e_or in children_keep_orientation(a_hole_wire) {
                    let Some(a_le_im) = the_dme_or_le_im.get(&shape_map_key(&a_e_or)) else {
                        continue;
                    };
                    for a_e_im in a_le_im {
                        if a_me_splits.contains(a_e_im) {
                            a_me_im_wire.add(a_e_im);
                        }
                    }
                }
                if a_me_im_wire.extent() == 0 {
                    continue;
                }

                // `cxx:3177-3187`: feed both orientations of every edge so
                // `BOPAlgo_BuilderFace` can assemble the wire in either sense.
                let mut a_le: Vec<TopoShape> = Vec::new();
                for i in 1..=a_me_im_wire.extent() {
                    let Some(a_e) = a_me_im_wire.value(i) else {
                        continue;
                    };
                    let mut a_ef = a_e.clone();
                    a_ef.set_orientation(Orientation::Forward);
                    a_le.push(a_ef);
                    let mut a_er = a_e.clone();
                    a_er.set_orientation(Orientation::Reversed);
                    a_le.push(a_er);
                }

                let mut a_bf = FaceBuilder::new();
                let mut a_f_offset_forward = the_f_offset.clone();
                a_f_offset_forward.0.set_orientation(Orientation::Forward);
                a_bf.set_face(&a_f_offset_forward);
                a_bf.set_shapes(&a_le);
                // `cxx:3188`: `aBF.Perform()`.
                let _ = perform_face_builder(&mut a_bf);

                let a_lf_new: Vec<TopoShape> = a_bf.areas().to_vec();
                if a_lf_new.is_empty() {
                    continue;
                }

                // `cxx:3194-3223`: an outer edge of the new faces that is
                // inverted means the hole was filled during the offset, so
                // there is nothing to remove.
                let mut a_dmef_new: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
                for a_f_new in &a_lf_new {
                    for (k, mut v) in map_shapes_and_ancestors(a_f_new, ShapeType::Edge, ShapeType::Face) {
                        a_dmef_new.entry(k).or_default().append(&mut v);
                    }
                }
                let b_has_inverted = a_dmef_new.values().any(|a_lf| {
                    a_lf.len() == 1 && self.my_inverted_edges.contains(&a_lf[0])
                });
                if b_has_inverted {
                    continue;
                }

                // `cxx:3225-3233`.
                for a_f_new in &a_lf_new {
                    a_lf_new_holes.push(a_f_new.clone());
                }
            }

            // `cxx:3126-3130`: the pointer stays bound with what was found
            // (possibly nothing).
            self.my_f_new_holes.insert(a_f_origin_key, a_lf_new_holes.clone());
        }

        // `cxx:3236-3242`: edge -> face maps of the splits and of the holes.
        let mut an_ef_splits_map: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut an_ef_holes_map: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:3245-3272`: classify every split against the new holes.
        for a_f_im_shape in the_lf_images {
            let Some(a_f_im) = Face::wrap(a_f_im_shape.clone()) else {
                continue;
            };
            map_shapes_and_ancestors_into(
                &a_f_im.0,
                ShapeType::Edge,
                ShapeType::Face,
                &mut an_ef_splits_map,
            );

            // `cxx:3250-3257`: a point inside the split.
            let Ok((a_p3d, _a_p2d)) = point_in_face(&a_f_im, &mut self.my_context) else {
                continue;
            };
            let a_tol = BRepTool::face_tolerance(&a_f_im);

            // `cxx:3260-3271`.
            for a_f_new_shape in &a_lf_new_holes {
                let Some(a_f_new) = Face::wrap(a_f_new_shape.clone()) else {
                    continue;
                };
                if self
                    .my_context
                    .is_point_in_on_face(&a_f_new, &a_p3d, None, a_tol)
                    .unwrap_or(false)
                {
                    the_mf_holes.add(&a_f_im.0);
                    map_shapes_and_ancestors_into(
                        &a_f_im.0,
                        ShapeType::Edge,
                        ShapeType::Face,
                        &mut an_ef_holes_map,
                    );
                    break;
                }
            }
        }

        // `cxx:3276-3308`: keep the holes whose removal would uncover an edge
        // connected to another offset face.
        //
        // `IndexedDataMap` key order is invisible here: every entry is tested
        // independently and only `theMFHoles` is mutated.
        let a_hole_entries: Vec<(ShapeMapKey, Vec<TopoShape>)> = an_ef_holes_map
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (an_edge, a_lf_holes) in a_hole_entries {
            // `cxx:3283-3287`: an outer edge of the hole.
            if a_lf_holes.len() != 1 {
                continue;
            }
            let a_f_hole = &a_lf_holes[0];
            // `cxx:3289-3293`.
            if !the_mf_holes.contains(a_f_hole) {
                continue;
            }
            // `cxx:3296-3299`: the edge must belong to exactly one split.
            let Some(a_l_splits) = an_ef_splits_map.get(&an_edge) else {
                continue;
            };
            if a_l_splits.len() == 1 {
                continue;
            }

            // `cxx:3302-3308`.
            let Some(a_lf_all) = the_ef_map.get(&an_edge) else {
                continue;
            };
            if a_lf_all.len() == 2 {
                the_mf_holes.remove_key(a_f_hole);
            }
        }
    }
}
