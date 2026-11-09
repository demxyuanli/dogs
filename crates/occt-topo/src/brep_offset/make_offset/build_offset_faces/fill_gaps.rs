//! `FillGaps` (`cxx:9297-9407`): fills the gaps (holes) the splitting left in
//! the images of the offset faces, so that the volume builder can still close
//! a solid.
//!
//! Two passes over `myOFImages`, exactly as the reference:
//!
//! 1. `cxx:9311-9325`: `anEFMap` collects the `(EDGE, FACE)` ancestors of
//!    every split of every offset face. The map is accumulated across the
//!    face loop (`TopExp::MapShapesAndAncestors` never resets it), so an edge
//!    shared by two splits carries two ancestors.
//! 2. `cxx:9329-9406`: per offset face, the edges of its splits are explored;
//!    a FORWARD/REVERSED edge whose ancestor list holds a single face is
//!    *free* (`cxx:9362-9366`) and every edge of the splits goes into the
//!    `anEdges` compound (`cxx:9368`). With at least one free edge the face is
//!    split again over that compound (`cxx:9387`), and every new split that
//!    carries a free edge is appended to the face's image list
//!    (`cxx:9393-9405`) - those pieces are the ones filling the holes.
//!
//! `Message_ProgressScope` (`cxx:9305`, `:9315`, `:9331`) is the only use of
//! `theRange`; the port has no progress reporting, so the two `return` arms it
//! guards are absent. Every other branch is kept.
//!
//! `anEFMap.FindFromKey(aE)` (`cxx:9362`) throws on an unbound key in OCCT;
//! every edge of the splits is bound by the first pass, so the port reads the
//! entry and treats a missing one as the empty list (extent 0, never free).
//! `aDummy` (`cxx:9385`) only receives the face-origin bindings of
//! `BuildSplitsOfFace` and is dropped, like the reference's local map.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors_into, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::splits::build_splits_of_face;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::FillGaps(const
    /// Message_ProgressRange&)` (`cxx:9297-9407`).
    pub fn fill_gaps(&mut self) {
        // `cxx:9299-9303`.
        let a_nb_f = self.my_of_images.extent();
        if a_nb_f == 0 {
            return;
        }

        // `cxx:9311-9325`: the (EDGE, FACE) ancestors of all splits.
        let mut an_ef_map: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        for i in 1..=a_nb_f {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_images: Vec<TopoShape> =
                self.my_of_images.seek(&a_f).map(|l| l.to_vec()).unwrap_or_default();
            for a_f_im in &a_lf_images {
                map_shapes_and_ancestors_into(
                    a_f_im,
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut an_ef_map,
                );
            }
        }

        // `cxx:9329-9406`: fill the holes of each offset face.
        for i in 1..=a_nb_f {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_images: Vec<TopoShape> =
                self.my_of_images.seek(&a_f).map(|l| l.to_vec()).unwrap_or_default();
            if a_lf_images.is_empty() {
                continue;
            }

            // `cxx:9343-9370`: the free edges and the `anEdges` compound.
            let mut a_free_edges_map = ShapeIndexedMap::new();
            let mut a_compound_edges: Vec<TopoShape> = Vec::new();
            for a_f_im in &a_lf_images {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    if a_e.orientation() != Orientation::Forward
                        && a_e.orientation() != Orientation::Reversed
                    {
                        // `cxx:9356-9360`: skip internals.
                        continue;
                    }

                    // `cxx:9362-9366`: a single ancestor makes the edge free;
                    // the map ignores the orientation, so the reversed copy
                    // answers `Contains` for the original edge too.
                    let a_lf_extent = an_ef_map.get(&shape_map_key(&a_e)).map_or(0, |l| l.len());
                    if a_lf_extent == 1 {
                        let a_reversed = a_e.oriented(a_e.orientation().reversed());
                        a_free_edges_map.add(&a_reversed);
                    }

                    a_compound_edges.push(a_e);
                }
            }

            // `cxx:9372-9376`.
            if a_free_edges_map.extent() == 0 {
                continue;
            }

            // `cxx:9380-9387`: split the face over all kept edges.
            let Some(a_face) = Face::wrap(a_f.clone()) else {
                continue;
            };
            let an_edges = make_compound_from(a_compound_edges.iter());
            let mut a_lf_new: Vec<TopoShape> = Vec::new();
            let mut a_dummy: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
            build_splits_of_face(&a_face, &an_edges, &mut a_dummy, &mut a_lf_new);

            // `cxx:9389-9405`: the splits carrying a free edge fill the holes.
            let mut a_to_add: Vec<TopoShape> = Vec::new();
            for a_f_new in &a_lf_new {
                if explorer_shapes(a_f_new, ShapeType::Edge)
                    .iter()
                    .any(|a_e| a_free_edges_map.contains(a_e))
                {
                    a_to_add.push(a_f_new.clone());
                }
            }

            if !a_to_add.is_empty() {
                if let Some(a_list) = self.my_of_images.change_seek(&a_f) {
                    a_list.extend(a_to_add);
                }
            }
        }
    }
}
