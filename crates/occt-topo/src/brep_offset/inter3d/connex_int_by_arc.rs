//! `ConnexIntByArc` (`BRepOffset_Inter3d.cxx:258-443`): the intersections of
//! the offset faces of the initial faces that meet at a sharp edge (stage 1),
//! then of the tubes of edges sharing a vertex that has no sphere (stage 2).

use crate::abs::{ShapeType, State};
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::{Analyse, BRepOffsetTool};
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::BRepOffsetInter3d;

impl BRepOffsetInter3d {
    /// `ConnexIntByArc(SetOfFaces, ShapeInit, Analyse, InitOffsetFace,
    /// theRange)` (`cxx:258-443`). `SetOfFaces` is unused in the reference
    /// (`/*SetOfFaces*/`).
    pub fn connex_int_by_arc(
        &mut self,
        shape_init: &TopoShape,
        analyse: &Analyse,
        init_offset_face: &BRepAlgoImage,
    ) {
        // `cxx:264-268`: the concavity of the processed side.
        let ot = if self.side == State::Out {
            TypeOfConcavity::Convex
        } else {
            TypeOfConcavity::Concave
        };

        let mut l_int1: Vec<TopoShape> = Vec::new();
        let mut l_int2: Vec<TopoShape> = Vec::new();

        //----------------------------------------------
        // etape 1 : intersection of the offset faces of the initial faces
        //           separated by an edge of type `ot` (`cxx:280-308`).
        //----------------------------------------------
        for e_shape in explorer_shapes(shape_init, ShapeType::Edge) {
            let e = Edge(e_shape);
            let l = analyse.type_of(&e);
            if l.is_empty() || l[0].kind() != ot {
                continue;
            }
            // `cxx:293-298`: `Analyse.Ancestors(E)`; `Anc.Extent() == 2`.
            let anc = analyse.ancestors(&e.0);
            if anc.len() != 2 {
                continue;
            }
            let init_f1 = Face(anc[0].clone());
            let init_f2 = Face(anc[anc.len() - 1].clone());

            // `F1 = TopoDS::Face(InitOffsetFace.Image(InitF1).First())`; an
            // empty image list raises in OCCT, the port skips the edge.
            let Some(f1) = first_image(init_offset_face, &init_f1.0) else {
                continue;
            };
            let Some(f2) = first_image(init_offset_face, &init_f2.0) else {
                continue;
            };
            if !self.is_done(&f1, &f2) {
                // `cxx:304-305`: the reference edge is `E`.
                BRepOffsetTool::inter3d(
                    &f1,
                    &f2,
                    &mut l_int1,
                    &mut l_int2,
                    self.side,
                    Some(&e),
                    Some(&init_f1),
                    Some(&init_f2),
                );
                self.store(&f1, &f2, &l_int1, &l_int2);
            }
        }

        //----------------------------------------------
        // etape 2 : intersections of the tubes sharing a vertex without sphere
        //           (`cxx:314-441`).
        //----------------------------------------------
        for e1_shape in explorer_shapes(shape_init, ShapeType::Edge) {
            let e1 = Edge(e1_shape);
            // `if (InitOffsetFace.HasImage(E1))`.
            if !init_offset_face.has_image(&e1.0) {
                continue;
            }
            // `F1 = TopoDS::Face(InitOffsetFace.Image(E1).First())`.
            let Some(f1) = first_image(init_offset_face, &e1.0) else {
                continue;
            };
            // `TopExp::Vertices(E1, V[0], V[1])` uses the default
            // `CumOri = false`.
            let v = edge_vertices(&e1);
            let v = [v.0, v.1];
            // `const NCollection_List<TopoDS_Shape>& AncE1 = Analyse.Ancestors(E1);`
            let anc_e1: Vec<TopoShape> = analyse.ancestors(&e1.0).to_vec();

            for i in 0..2 {
                let Some(vi) = v[i].as_ref() else {
                    continue;
                };
                // `if (!InitOffsetFace.HasImage(V[i]))`: the vertex has no
                // sphere.
                if init_offset_face.has_image(&vi.0) {
                    continue;
                }
                let anc: Vec<TopoShape> = analyse.ancestors(&vi.0).to_vec();
                let mut tang_on_v: Vec<TopoShape> = Vec::new();
                analyse.tangent_edges(&e1, vi, &mut tang_on_v);
                // `NCollection_Map<TopoDS_Shape> MTEV` (`IsSame` set).
                let mtev: std::collections::HashSet<ShapeMapKey> =
                    tang_on_v.iter().map(shape_map_key).collect();

                for e2_shape in &anc {
                    if e2_shape.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    let e2 = Edge(e2_shape.clone());
                    // `cxx:353-366`: the OCC4455 skip rule.
                    let mut is_to_skip = false;
                    if !is_same(&e1.0, &e2.0) {
                        let al = analyse.type_of(&e2);
                        is_to_skip = mtev.contains(&shape_map_key(e2_shape))
                            && (al.is_empty() || al[0].kind() != ot);
                    }
                    if is_same(&e1.0, &e2.0) || is_to_skip {
                        continue;
                    }

                    if init_offset_face.has_image(&e2.0) {
                        // `cxx:369-383`: `E2` generated a tube.
                        let Some(f2) = first_image(init_offset_face, &e2.0) else {
                            continue;
                        };
                        if !self.is_done(&f1, &f2) {
                            // `cxx:379`: `PipeInter(F1, F2, LInt1, LInt2,
                            // mySide)`.
                            BRepOffsetTool::pipe_inter(
                                &f1,
                                &f2,
                                &mut l_int1,
                                &mut l_int2,
                                self.side,
                            );
                            self.store(&f1, &f2, &l_int1, &l_int2);
                        }
                    } else {
                        // `cxx:388-431`: the tube of `E1` against the offset
                        // faces of the faces carrying `E2`.
                        let l = analyse.type_of(&e2);
                        if !l.is_empty() && l[0].kind() == TypeOfConcavity::Tangential {
                            continue;
                        }
                        let anc_e2: Vec<TopoShape> = analyse.ancestors(&e2.0).to_vec();
                        if anc_e2.len() != 2 {
                            continue;
                        }
                        let first = Face(anc_e2[0].clone());
                        let last = Face(anc_e2[anc_e2.len() - 1].clone());
                        for init_f2 in [first, last] {
                            // `TangentFaces = InitF2.IsSame(AncE1.First()) ||
                            // InitF2.IsSame(AncE1.Last())`.
                            let tangent_faces = anc_e1
                                .first()
                                .map(|s| is_same(&init_f2.0, s))
                                .unwrap_or(false)
                                || anc_e1
                                    .last()
                                    .map(|s| is_same(&init_f2.0, s))
                                    .unwrap_or(false);
                            if tangent_faces {
                                continue;
                            }
                            let Some(f2) = first_image(init_offset_face, &init_f2.0) else {
                                continue;
                            };
                            if !self.is_done(&f1, &f2) {
                                BRepOffsetTool::inter3d(
                                    &f1,
                                    &f2,
                                    &mut l_int1,
                                    &mut l_int2,
                                    self.side,
                                    None,
                                    None,
                                    None,
                                );
                                self.store(&f1, &f2, &l_int1, &l_int2);
                            }
                        }
                    }
                }
            }
        }
    }
}

/// `TopoDS::Face(InitOffsetFace.Image(S).First())`: the head of the image
/// list, `None` where OCCT's `First()` on an empty list raises.
fn first_image(image: &BRepAlgoImage, s: &TopoShape) -> Option<Face> {
    image.image(s).first().map(|i| Face(i.clone()))
}
