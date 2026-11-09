//! `BRepFeat_RibSlot::UpdateDescendants(const BRepAlgoAPI_BooleanOperation&,
//! const TopoDS_Shape&, bool)` (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:2672-2736`).
//!
//! Rewrites every `myMap` entry through the boolean history of `aBOP`: a face
//! the result kept stays, a face the operation changed is replaced by
//! `Modified(face)`, and the survivors are filtered again on membership in the
//! result. `NoSlidingProfile` calls it after the `BndFace - profile` cut
//! (`cxx:2644`).

use crate::brep_feat::algo_api::BooleanOperation;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

use crate::abs::ShapeType;

use super::RibSlot;

/// `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` add (`IsSame` key).
fn add_unique(list: &mut Vec<TopoShape>, s: &TopoShape) {
    if !list.iter().any(|e| is_same(e, s)) {
        list.push(s.clone());
    }
}

impl RibSlot {
    /// `UpdateDescendants(aBOP, S, SkipFace)` (`cxx:2672-2736`).
    pub fn update_descendants_bop(
        &mut self,
        a_bop: &BooleanOperation,
        s: &TopoShape,
        skip_face: bool,
    ) {
        // `cxx:2682`: `TopExp_Explorer exp` over the result, queried per face.
        let result_faces = explorer_shapes(s, ShapeType::Face);
        let entries: Vec<(TopoShape, Vec<TopoShape>)> =
            self.my_map.iter().cloned().collect();
        for (orig, value) in entries {
            // `cxx:2685-2688`.
            if skip_face && orig.shape_type() == ShapeType::Face {
                continue;
            }
            // `cxx:2689-2720`.
            let mut newdsc: Vec<TopoShape> = Vec::new();
            for sh in &value {
                // `cxx:2697-2699`: only faces are rewritten.
                if sh.shape_type() != ShapeType::Face {
                    continue;
                }
                let fdsc = Face(sh.clone());
                // `cxx:2701-2709`: a face the result preserved stays.
                if result_faces.iter().any(|f| is_same(f, &fdsc.0)) {
                    add_unique(&mut newdsc, &fdsc.0);
                    continue;
                }
                // `cxx:2710-2719`: otherwise the `Modified` faces replace it.
                for a_s in a_bop.modified(&fdsc.0) {
                    add_unique(&mut newdsc, &a_s);
                }
            }
            // `cxx:2721-2733`.
            if let Some(slot) = self.my_map.change(&orig) {
                slot.clear();
            }
            for d in &newdsc {
                if result_faces.iter().any(|f| is_same(f, d)) {
                    if let Some(slot) = self.my_map.change(&orig) {
                        slot.push(d.clone());
                    }
                }
            }
        }
    }
}
