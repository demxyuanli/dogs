//! Bounds of the offset-face splits and the split filtering that consumes
//! them: `GetBounds` (`cxx:8551-8587`), `GetBoundsToUpdate` (`cxx:8591-8667`)
//! and `FilterSplits` (`cxx:8994-9042`).
//!
//! `GetBounds` collects the edges of the splits of the given offset faces into
//! one compound, skipping the edges already recorded as bounding (`theMEB`) and
//! keeping the first occurrence of every edge (`aMFence`).
//!
//! `GetBoundsToUpdate` does the same for a set of faces whose *descendant*
//! edges are taken through `myOEImages` (an untrimmed edge maps to its splits,
//! an edge without images is used as is), and additionally answers the offset
//! edges to re-image (`theLAValid`, the origins of the split edge when it has
//! any, the edge itself otherwise).
//!
//! `FilterSplits` drops from `theEImages` every split whose membership in
//! `theMEFilter` equals `theIsInv`, unbinding the key once its list runs empty,
//! and appends the survivors to the answered compound.
//!
//! `GetInvalidEdgesByBounds` (`cxx:8672-8990`), the only other consumer of the
//! bounds and the reason `theMVOld`/`theMENew`/`theDMEOr` exist, is UNPORTED: it
//! needs `BOPAlgo_Section` (`cxx:8695-8700`), which the port does not have yet.

use crate::abs::ShapeType;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void GetBounds(const List& theLFaces, const Map& theMEB, Shape&
    /// theBounds)` (`cxx:8551-8587`): the compound of the split edges of
    /// `theLFaces`.
    pub fn get_bounds(
        &self,
        the_lfaces: &[TopoShape],
        the_meb: &ShapeIndexedMap,
    ) -> TopoShape {
        // `cxx:8561`: the fence of the collected edges.
        let mut a_mfence = ShapeIndexedMap::new();
        // `cxx:8559-8560`.
        let mut a_bounds_edges: Vec<TopoShape> = Vec::new();

        for a_f in the_lfaces {
            // `cxx:8566-8570`: a face without splits contributes nothing.
            let Some(a_lf_im) = self.my_of_images.seek(a_f) else {
                continue;
            };
            for a_f_im in a_lf_im {
                // `cxx:8576-8585`.
                for a_e_im in explorer_shapes(a_f_im, ShapeType::Edge) {
                    if !the_meb.contains(&a_e_im) && a_mfence.add(&a_e_im) == a_mfence.extent() {
                        a_bounds_edges.push(a_e_im);
                    }
                }
            }
        }

        // `cxx:8587`.
        make_compound_from(a_bounds_edges.iter())
    }

    /// `void GetBoundsToUpdate(const List& theLF, const Map& theMEB, List&
    /// theLABounds, List& theLAValid, Shape& theBounds)` (`cxx:8591-8667`).
    ///
    /// Answers `(theLABounds, theLAValid, theBounds)`.
    pub fn get_bounds_to_update(
        &self,
        the_lf: &[TopoShape],
        the_meb: &ShapeIndexedMap,
    ) -> (Vec<TopoShape>, Vec<TopoShape>, TopoShape) {
        // `cxx:8605-8606`.
        let mut a_ma_valid = ShapeIndexedMap::new();
        let mut a_mfence = ShapeIndexedMap::new();
        let mut a_bounds_edges: Vec<TopoShape> = Vec::new();
        let mut a_labounds: Vec<TopoShape> = Vec::new();
        let mut a_lavalid: Vec<TopoShape> = Vec::new();

        for a_f in the_lf {
            // `cxx:8614-8634`: the edges of the face's splits, keyed by the
            // face's descendants.
            let mut a_mde = ShapeIndexedMap::new();
            let a_lf_des: Vec<TopoShape> = self.my_as_des.descendant(a_f).to_vec();
            for a_ed in &a_lf_des {
                match self.my_oe_images.seek(a_ed) {
                    // `cxx:8620-8624`: an edge without images is used as is.
                    None => {
                        a_mde.add(a_ed);
                    }
                    Some(a_led_im) => {
                        for a_ed_im in a_led_im {
                            a_mde.add(a_ed_im);
                        }
                    }
                }
            }

            // `cxx:8636-8665`.
            let a_nb_e = a_mde.extent();
            for j in 1..=a_nb_e {
                let Some(a_e_im) = a_mde.value(j).cloned() else {
                    continue;
                };
                if !the_meb.contains(&a_e_im) && a_mfence.add(&a_e_im) == a_mfence.extent() {
                    a_bounds_edges.push(a_e_im.clone());
                    a_labounds.push(a_e_im.clone());
                }
                // `cxx:8646-8665`: the offset edges to re-image.
                match self.my_oe_origins.get(&shape_map_key(&a_e_im)) {
                    Some(a_lo) => {
                        for a_eo in a_lo {
                            if a_ma_valid.add(a_eo) == a_ma_valid.extent() {
                                a_lavalid.push(a_eo.clone());
                            }
                        }
                    }
                    None => {
                        if a_ma_valid.add(&a_e_im) == a_ma_valid.extent() {
                            a_lavalid.push(a_e_im);
                        }
                    }
                }
            }
        }

        (
            a_labounds,
            a_lavalid,
            make_compound_from(a_bounds_edges.iter()),
        )
    }

    /// `void FilterSplits(const List& theLE, const Map& theMEFilter, const bool
    /// theIsInv, DataMap& theEImages, Shape& theSplits)` (`cxx:8994-9042`).
    ///
    /// Removes from `theEImages` the splits whose `Contains` in
    /// `theMEFilter` equals `theIsInv`, unbinds the keys whose list empties,
    /// and answers the compound of the surviving splits.
    pub fn filter_splits(
        &self,
        the_le: &[TopoShape],
        the_me_filter: &ShapeIndexedMap,
        the_is_inv: bool,
        the_e_images: &mut std::collections::HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) -> TopoShape {
        // `cxx:9006`.
        let mut a_mfence = ShapeIndexedMap::new();
        let mut a_splits_edges: Vec<TopoShape> = Vec::new();

        for a_e in the_le {
            // `cxx:9012-9016`.
            let Some(a_le_im) = the_e_images.get(&shape_map_key(a_e)) else {
                continue;
            };
            let a_le_im = a_le_im.clone();

            // `cxx:9018-9037`: the survivors of this key; removal and `Next()`
            // interleave, so the walk keeps the order of the kept splits.
            let mut a_kept: Vec<TopoShape> = Vec::with_capacity(a_le_im.len());
            for a_e_im in a_le_im {
                if the_me_filter.contains(&a_e_im) == the_is_inv {
                    continue;
                }
                if a_mfence.add(&a_e_im) == a_mfence.extent() {
                    a_splits_edges.push(a_e_im.clone());
                }
                a_kept.push(a_e_im);
            }

            // `cxx:9035-9038`.
            if a_kept.is_empty() {
                the_e_images.remove(&shape_map_key(a_e));
            } else {
                the_e_images.insert(shape_map_key(a_e), a_kept);
            }
        }

        // `cxx:9041`.
        make_compound_from(a_splits_edges.iter())
    }
}
