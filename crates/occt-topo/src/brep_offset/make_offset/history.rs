//! `BRepOffset_MakeOffset::Generated` (`BRepOffset_MakeOffset.cxx:5231-5330`),
//! `Modified` (`:5334-5365`) and `IsDeleted` (`:5369-5382`): the history
//! queries of the offset driver.
//!
//! `Generated` switches on the shape type and falls through (`:5237-5307`):
//!
//! * `VERTEX` (`:5237-5266`): when the vertex has ancestors in `myAnalyse`,
//!   walk them, take `myInitOffsetEdge.LastImage` of the first edge that has
//!   one, and append the first vertex the `TopoDS_Iterator` sees *twice* (the
//!   `aMFence.Add` guard at `:5258`).
//! * `EDGE` (`:5269-5274`): `myInitOffsetEdge.LastImage(theS, myGenerated)`.
//! * `FACE` (`:5276-5299`): the plan-face substitution from
//!   `myFacePlanfaceMap` (`:5277-5282`), then `myInitOffsetFace.LastImage` when
//!   the face is not a closing face (`:5284-5298`); with closing faces every
//!   image is reversed in place (`:5292-5296`).
//! * `SOLID` (`:5301-5306`): `theS.IsSame(myShape)` appends `myOffsetShape`.
//!
//! The fallthrough means a `VERTEX` also runs the `EDGE` and `FACE` arms with
//! the same `theS`, and an `EDGE` also runs the `FACE` arm; both tests are
//! normally false for those types, and the port keeps them anyway instead of
//! short-circuiting, so the reference control flow stays visible.
//!
//! `myResMap` caches `TopExp::MapShapes(myOffsetShape)` (`:5312-5315`,
//! `:5371-5374`) and every generated shape outside it is dropped
//! (`:5317-5327`). `Modified` is the FACE-only counterpart whose
//! `myFaces.Contains(aS)` test is *not* negated (`:5347`). `IsDeleted` is
//! `!myResMap.Contains(S) && Generated(S).IsEmpty() && Modified(S).IsEmpty()`
//! (`:5376-5381`).
//!
//! Port notes:
//!
//! * `it.ChangeValue().Reverse()` (`:5295`, `:5358`) mutates the list entry,
//!   so the port reverses the stored list ([`TopoShape::reverse`]).
//! * `TopoDS_Iterator(itLEIm.Value())` is [`children_keep_orientation`]
//!   (`TopoDS_Iterator(S, false)`).
//! * `myGenerated` is kept as the scratch list it is in OCCT; the accessors
//!   answer the list they just built, which is the same content.

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{
    all_subshapes, children_keep_orientation, is_same, shape_map_key, ShapeMapKey,
};

use super::BRepOffsetMakeOffset;

/// `TopExp::MapShapes(theShape, theMap)` (`TopExp.cxx`, includes the shape
/// itself) into the key set `myResMap` holds (`BRepOffset_MakeOffset.hxx:283`).
fn map_shapes_into(shape: &TopoShape, out: &mut std::collections::HashSet<ShapeMapKey>) {
    out.insert(shape_map_key(shape));
    for sub in all_subshapes(shape) {
        out.insert(shape_map_key(&sub));
    }
}

/// `it.ChangeValue().Reverse()` over a list (`:5295`, `:5358`).
fn reverse_all(list: &mut [TopoShape]) {
    for s in list.iter_mut() {
        s.reverse();
    }
}

impl BRepOffsetMakeOffset {
    /// `Generated(theS)` (`cxx:5231-5330`).
    pub fn generated(&mut self, the_s: &TopoShape) -> Vec<TopoShape> {
        self.my_generated.clear();
        let mut out: Vec<TopoShape> = Vec::new();
        let ty = the_s.shape_type();

        // `cxx:5237-5266`: case TopAbs_VERTEX.
        if ty == ShapeType::Vertex && self.my_analyse.has_ancestor(the_s) {
            let mut a_m_fence: Vec<ShapeMapKey> = Vec::new();
            let ancestors = self.my_analyse.ancestors(the_s).to_vec();
            for a_e in &ancestors {
                if !out.is_empty() {
                    break;
                }
                if !self.my_init_offset_edge.has_image(a_e) {
                    continue;
                }
                let mut a_le_im: Vec<TopoShape> = Vec::new();
                self.my_init_offset_edge.last_image(a_e, &mut a_le_im);
                for the_edge_image in &a_le_im {
                    if !out.is_empty() {
                        break;
                    }
                    for a_v in children_keep_orientation(the_edge_image) {
                        let key = shape_map_key(&a_v);
                        if a_m_fence.contains(&key) {
                            // `!aMFence.Add(itV.Value())`: already fenced.
                            out.push(a_v);
                            break;
                        }
                        a_m_fence.push(key);
                    }
                }
            }
        }

        // `cxx:5269-5274`: case TopAbs_EDGE, reached by the VERTEX fallthrough
        // as well.
        if matches!(ty, ShapeType::Vertex | ShapeType::Edge)
            && self.my_init_offset_edge.has_image(the_s)
        {
            self.my_init_offset_edge.last_image(the_s, &mut out);
        }

        // `cxx:5276-5299`: case TopAbs_FACE, reached by the VERTEX and EDGE
        // fallthroughs as well.
        if matches!(
            ty,
            ShapeType::Vertex | ShapeType::Edge | ShapeType::Face
        ) {
            let mut a_s = the_s.clone();
            if let Some(planface) = self.my_face_planface_map.get(&shape_map_key(the_s)) {
                a_s = planface.clone();
            }

            if !self.my_faces.contains(&a_s) && self.my_init_offset_face.has_image(&a_s) {
                let mut list = Vec::new();
                self.my_init_offset_face.last_image(&a_s, &mut list);

                if self.my_faces.extent() != 0 {
                    // Reverse generated shapes in case of small solids.
                    reverse_all(&mut list);
                }
                out.extend(list);
            }
        }

        // `cxx:5301-5306`: case TopAbs_SOLID.
        if ty == ShapeType::Solid {
            if let Some(shape) = &self.my_shape {
                if is_same(shape, the_s) {
                    if let Some(offset_shape) = &self.my_offset_shape {
                        out.push(offset_shape.clone());
                    }
                }
            }
        }

        // `cxx:5312-5315`.
        if self.my_res_map.is_empty() {
            if let Some(offset_shape) = self.my_offset_shape.clone() {
                map_shapes_into(&offset_shape, &mut self.my_res_map);
            }
        }

        // `cxx:5317-5327`.
        out.retain(|s| self.my_res_map.contains(&shape_map_key(s)));

        self.my_generated = out.clone();
        out
    }

    /// `Modified(theShape)` (`cxx:5334-5365`): the FACE arm only, and the
    /// `myFaces.Contains(aS)` test is *not* negated.
    pub fn modified(&mut self, the_shape: &TopoShape) -> Vec<TopoShape> {
        self.my_generated.clear();
        let mut out: Vec<TopoShape> = Vec::new();

        if the_shape.shape_type() == ShapeType::Face {
            // `cxx:5340-5345`.
            let mut a_s = the_shape.clone();
            if let Some(planface) = self.my_face_planface_map.get(&shape_map_key(the_shape)) {
                a_s = planface.clone();
            }

            // `cxx:5347-5361`.
            if self.my_faces.contains(&a_s) && self.my_init_offset_face.has_image(&a_s) {
                let mut list = Vec::new();
                self.my_init_offset_face.last_image(&a_s, &mut list);
                if self.my_faces.extent() != 0 {
                    reverse_all(&mut list);
                }
                out.extend(list);
            }
        }

        self.my_generated = out.clone();
        out
    }

    /// `IsDeleted(theS)` (`cxx:5369-5382`).
    pub fn is_deleted(&mut self, the_s: &TopoShape) -> bool {
        // `cxx:5371-5374`.
        if self.my_res_map.is_empty() {
            if let Some(offset_shape) = self.my_offset_shape.clone() {
                map_shapes_into(&offset_shape, &mut self.my_res_map);
            }
        }

        // `cxx:5376-5379`.
        if self.my_res_map.contains(&shape_map_key(the_s)) {
            return false;
        }

        // `cxx:5381`.
        self.generated(the_s).is_empty() && self.modified(the_s).is_empty()
    }
}
