//! `CompletInt` (`BRepOffset_Inter3d.cxx:92-148`) and its single pair worker
//! `FaceInter` (`:155-254`).

use crate::abs::ShapeType;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_bnd_lib::add as bnd_add;
use crate::brep_offset::BRepOffsetTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, is_same};
use occt_core::bnd::BndBox;

use super::BRepOffsetInter3d;

impl BRepOffsetInter3d {
    /// `CompletInt(SetOfFaces, InitOffsetFace, theRange)` (`cxx:92-148`): every
    /// pair of faces with overlapping bounding boxes is intersected by
    /// [`Self::face_inter`].
    ///
    /// The reference builds a `BOPTools_BoxTree` over the boxes and runs
    /// `BOPTools_BoxPairSelector` with `SetSame(true)`, i.e. it keeps the pairs
    /// `(i, j)` with `i < j` (insertion indices) whose boxes are not
    /// [`Bnd_Box::IsOut`] of each other, then `Sort()`s them by
    /// `(i, j)` (`BOPTools_PairSelector.hxx:40-104`). Both steps are only a
    /// spatial filter over the same predicate, so the nested `i < j` loop below
    /// visits exactly the same pairs in exactly the same order.
    pub fn complet_int(&mut self, set_of_faces: &[TopoShape], init_offset_face: &BRepAlgoImage) {
        // `cxx:105-124`: the ordered faces and their boxes.
        let mut faces: Vec<Face> = Vec::new();
        let mut boxes: Vec<BndBox> = Vec::new();
        for s in set_of_faces {
            if s.shape_type() != ShapeType::Face {
                continue;
            }
            let f = Face(s.clone());
            let mut b = BndBox::new();
            // `BRepBndLib::Add(aF, aBoxF)`.
            bnd_add(&f.0, &mut b, false);
            faces.push(f);
            boxes.push(b);
        }

        // `cxx:126-141`: `SetSame(true)` + `Sort()`.
        for i in 0..faces.len() {
            for j in (i + 1)..faces.len() {
                // `BVH_Box::IsOut` / `Bnd_Box::IsOut`.
                if boxes[i].is_out_box(&boxes[j]) {
                    continue;
                }
                self.face_inter(&faces[i], &faces[j], init_offset_face);
            }
        }
    }

    /// `FaceInter(F1, F2, InitOffsetFace)` (`cxx:155-254`).
    pub fn face_inter(&mut self, f1: &Face, f2: &Face, init_offset_face: &BRepAlgoImage) {
        let mut l_int1: Vec<TopoShape> = Vec::new();
        let mut l_int2: Vec<TopoShape> = Vec::new();

        // `cxx:162-170`.
        if is_same(&f1.0, &f2.0) {
            return;
        }
        if self.is_done(f1, f2) {
            return;
        }

        // `cxx:172-177`: `ImageFrom` raises for an unbound face in OCCT; the
        // port returns.
        let (Ok(init_f1), Ok(init_f2)) = (
            init_offset_face.image_from(&f1.0),
            init_offset_face.image_from(&f2.0),
        ) else {
            return;
        };
        if is_same(&init_f1, &init_f2) {
            return;
        }

        // `cxx:179-180`.
        let inter_pipes = init_f2.shape_type() == ShapeType::Edge
            && init_f1.shape_type() == ShapeType::Edge;
        let inter_faces = init_f1.shape_type() == ShapeType::Face
            && init_f2.shape_type() == ShapeType::Face;

        // `cxx:182-183`.
        let mut le: Vec<TopoShape> = Vec::new();
        let mut lv: Vec<TopoShape> = Vec::new();

        // `cxx:185`: `FindCommonShapes(F1, F2, LE, LV)` short-circuits the
        // `HasCommonDescendant`, which clears `LE` before filling it.
        if BRepOffsetTool::find_common_shapes(f1, f2, &mut le, &mut lv)
            || self.as_des.has_common_descendant(&f1.0, &f2.0, &mut le)
        {
            // `cxx:189-191`.
            if le.is_empty() && !lv.is_empty() {
                if inter_pipes {
                    // `cxx:192-216`: the two tubes share a vertex. When that
                    // vertex has no image (no sphere) the reference intersects
                    // the tubes through `PipeInter`.
                    let ee1 = Edge(init_f1.clone());
                    let ee2 = Edge(init_f2.clone());
                    let v1s = edge_vertices(&ee1);
                    let v2s = edge_vertices(&ee2);
                    // `cxx:202-210`: no `break`, the last match wins.
                    let mut shared: Option<TopoShape> = None;
                    for a in [v1s.0.as_ref(), v1s.1.as_ref()].into_iter().flatten() {
                        for b in [v2s.0.as_ref(), v2s.1.as_ref()].into_iter().flatten() {
                            if is_same(&a.0, &b.0) {
                                shared = Some(a.0.clone());
                            }
                        }
                    }
                    if let Some(v) = shared {
                        // `if (!InitOffsetFace.HasImage(V)) { // no sphere`.
                        if !init_offset_face.has_image(&v) {
                            // `cxx:214`: `PipeInter(F1, F2, LInt1, LInt2,
                            // mySide)`.
                            BRepOffsetTool::pipe_inter(
                                f1,
                                f2,
                                &mut l_int1,
                                &mut l_int2,
                                self.side,
                            );
                        }
                    }
                } else if inter_faces {
                    // `cxx:218-239`.
                    let init_f1_face = Face(init_f1.clone());
                    let init_f2_face = Face(init_f2.clone());
                    if BRepOffsetTool::find_common_shapes(
                        &init_f1_face,
                        &init_f2_face,
                        &mut le,
                        &mut lv,
                    ) {
                        if !le.is_empty() {
                            BRepOffsetTool::inter3d(
                                f1, f2, &mut l_int1, &mut l_int2, self.side, None, None, None,
                            );
                        }
                    } else {
                        BRepOffsetTool::inter3d(
                            f1, f2, &mut l_int1, &mut l_int2, self.side, None, None, None,
                        );
                    }
                }
            }
        } else {
            // `cxx:243-253`.
            if inter_pipes {
                // `cxx:246`: `PipeInter(F1, F2, LInt1, LInt2, mySide)`.
                BRepOffsetTool::pipe_inter(f1, f2, &mut l_int1, &mut l_int2, self.side);
            } else {
                BRepOffsetTool::inter3d(
                    f1, f2, &mut l_int1, &mut l_int2, self.side, None, None, None,
                );
            }
        }

        // `cxx:254`.
        self.store(f1, f2, &l_int1, &l_int2);
    }
}
