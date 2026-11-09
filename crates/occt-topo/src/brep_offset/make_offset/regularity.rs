//! `BRepOffset_MakeOffset::EncodeRegularity` (`BRepOffset_MakeOffset.cxx:
//! 3962-4194`): code the `GeomAbs_G1` / `GeomAbs_CN` regularity of the result's
//! edges.
//!
//! Every edge of `myOffsetShape` is walked once (`MS`, `:3984-3992`) and gets
//! a 3D curve first (`BRepLib::BuildCurve3d(OE, myTol)`, `:3987`). The edge is
//! then mapped back to its root image (`:3994-3997`) whose `AsDes` ascendants
//! must be exactly two faces (`:3999-4008`). The two roots in
//! `myInitOffsetFace` (`:4017-4018`) decide the branch:
//!
//! * same face (`:4024-4066`): a seam edge. In `GeomAbs_Intersection` mode a
//!   cylindrical / conical / spherical / toroidal face codes `GeomAbs_CN`
//!   (`:4030-4041`); otherwise the initial edge image's own continuity decides
//!   `GeomAbs_CN` vs `GeomAbs_G1` (`:4043-4063`).
//! * `TYPE1 == EDGE` and `TYPE2 == VERTEX` (or swapped): the vertex is an
//!   endpoint of that edge -> `GeomAbs_G1` (`:4077-4092`).
//! * face vs edge: the edge is contained in the face -> `GeomAbs_G1`
//!   (`:4093-4117`).
//! * face vs face: the two roots share exactly one edge which
//!   `BRepOffset_Analyse` types `ChFiDS_Tangential` -> `GeomAbs_G1`
//!   (`:4119-4137`).
//! * edge vs edge: the two roots share exactly one vertex on which
//!   `BRepOffset_Analyse::TangentEdges` reaches the other root ->
//!   `GeomAbs_G1` (`:4138-4156`).
//!
//! Port notes:
//!
//! * `BRep_Builder::Continuity(E, F1, F2, C)` (`:4040`, `:4056`, `:4060`,
//!   `:4084`, `:4091`, `:4103`, `:4112`, `:4134`, `:4151`) is
//!   `brep_lib_continuity::set_stored_continuity`, with `GeomAbs_CN = 6` and
//!   `GeomAbs_G1 = 1` in the crate's encoding.
//! * `BRep_Tool::Continuity(EI, FRoot, FRoot)` (`:4055`) is
//!   `brep_lib_continuity::stored_continuity`.
//! * `BRep_Tool::IsClosed(OE, F1)` (`:4026`) is
//!   `BRepTool::is_closed_edge_face`.
//! * `BRepAdaptor_Surface(F1, false).GetType()` (`:4033-4035`) is
//!   `pcurve_full::classify_surface_kind` on the registered surface; the
//!   `GeomAbs_Plane` / `GeomAbs_OtherCurve` arms are outside the test.
//! * `TopoDS::Face` / `TopoDS::Edge` / `TopoDS::Vertex` (`:4009`, `:4010`,
//!   `:4051`, `:4052`, `:4080`, `:4088`, `:4097`, `:4110`, `:4127`, `:4144`)
//!   throw on a type mismatch; the port guards with the `wrap` constructors.
//! * `TopoDS_Shape::IsSame` (`:4024`, `:4081`, `:4090`, `:4100`, `:4112`,
//!   `:4149`) is the `ShapeMapKey` comparison.
//! * `BRepAlgo_Image::Root` / `ImageFrom` (`:3996`, `:4053`) raise
//!   `Standard_NoSuchObject` when unbound; the port guards with `Result`.
//! * the `OCCT_DEBUG` / `OCCT_DEBUG_VERB` blocks (`:3964-3972`, `:4004-4006`,
//!   `:4158-4163`).

use std::collections::HashSet;

use occt_core::kernel::geomabs::{JoinType, Shape};
use crate::abs::ShapeType;
use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_lib_continuity::{set_stored_continuity, stored_continuity};
use crate::brep_offset::BRepOffsetTool as _Tool;
use crate::brep_surface::SurfaceKind;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::TypeOfConcavity;
use crate::pcurve_full::classify_surface_kind;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};

use super::BRepOffsetMakeOffset;

/// `GeomAbs_CN` in the crate's `GeomAbs_Shape` encoding.
const GEOM_ABS_CN: u8 = 6;
/// `GeomAbs_G1` in the crate's `GeomAbs_Shape` encoding.
const GEOM_ABS_G1: u8 = 1;

/// `BRepAdaptor_Surface(aFace, false).GetType()` (`cxx:4033`): `None` when the
/// face carries no registered surface (OCCT would throw).
fn surface_kind(face: &Face) -> Option<SurfaceKind> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    Some(classify_surface_kind(surf.as_ref()))
}

impl BRepOffsetMakeOffset {
    /// `EncodeRegularity()` (`cxx:3962-4194`).
    pub fn encode_regularity(&mut self) {
        let Some(shape) = self.my_offset_shape.clone() else {
            return;
        };

        // `cxx:3984`: `NCollection_Map MS`.
        let mut ms: HashSet<ShapeMapKey> = HashSet::new();

        for oe_shape in explorer_shapes(&shape, ShapeType::Edge) {
            let Some(oe) = Edge::wrap(oe_shape) else {
                continue;
            };
            // `cxx:3987`: `BRepLib::BuildCurve3d(OE, myTol)`, defaults
            // `GeomAbs_C1`, `MaxDegree = 14`, `MaxSegment = 0`.
            let _ = build_curve3d(&oe, self.my_tol, Shape::C1, 14, 0);

            // `cxx:3988-3992`.
            if !ms.insert(shape_map_key(&oe.0)) {
                continue;
            }

            // `cxx:3994-3997`.
            let mut roe = oe.clone();
            if self.my_image_offset.is_image(&oe.0) {
                if let Ok(root) = self.my_image_offset.root(&oe.0) {
                    if let Some(e) = Edge::wrap(root) {
                        roe = e;
                    }
                }
            }

            // `cxx:3999-4008`.
            let lof = self.my_as_des.ascendant(&roe.0);
            if lof.len() != 2 {
                continue;
            }
            let (Some(f1), Some(f2)) = (Face::wrap(lof[0].clone()), Face::wrap(lof[1].clone()))
            else {
                continue;
            };

            // `cxx:4017-4018`.
            let (Ok(root1), Ok(root2)) = (
                self.my_init_offset_face.root(&f1.0),
                self.my_init_offset_face.root(&f2.0),
            ) else {
                continue;
            };
            let type1 = root1.shape_type();
            let type2 = root2.shape_type();

            // `cxx:4024-4066`: the edge is a seam of one face.
            if shape_map_key(&f1.0) == shape_map_key(&f2.0) {
                if BRepTool::is_closed_edge_face(&oe, &f1) {
                    if self.my_join == JoinType::Intersection {
                        // `cxx:4030-4041`.
                        if let Some(kind) = surface_kind(&f1) {
                            if matches!(
                                kind,
                                SurfaceKind::Cylinder
                                    | SurfaceKind::Cone
                                    | SurfaceKind::Sphere
                                    | SurfaceKind::Torus
                            ) {
                                set_stored_continuity(&oe, &f1, &f1, GEOM_ABS_CN);
                            }
                        }
                    } else if self.my_init_offset_edge.is_image(&roe.0) {
                        // `cxx:4043-4063`.
                        if type1 == ShapeType::Face && type2 == ShapeType::Face {
                            let ei = self.my_init_offset_edge.image_from(&roe.0).ok();
                            if let (Some(fr), Some(ei)) = (Face::wrap(root1.clone()), ei.and_then(|e| Edge::wrap(e)))
                            {
                                let conti = stored_continuity(&ei, &fr, &fr);
                                if conti == GEOM_ABS_CN {
                                    set_stored_continuity(&oe, &f1, &f1, GEOM_ABS_CN);
                                } else if conti > 0 {
                                    set_stored_continuity(&oe, &f1, &f1, GEOM_ABS_G1);
                                }
                            }
                        }
                    }
                }
                continue;
            }

            // `cxx:4077-4092`: edge root vs vertex root.
            if type1 == ShapeType::Edge && type2 == ShapeType::Vertex {
                if let (Some(edge), Some(v)) = (Edge::wrap(root1.clone()), Vertex::wrap(root2.clone()))
                {
                    let (v1, v2) = edge_vertices(&edge);
                    let k = shape_map_key(&v.0);
                    let (k1, k2) = (
                        v1.map(|x| shape_map_key(&x.0)),
                        v2.map(|x| shape_map_key(&x.0)),
                    );
                    if Some(k.clone()) == k1 || Some(k) == k2 {
                        set_stored_continuity(&oe, &f1, &f2, GEOM_ABS_G1);
                    }
                }
            } else if type1 == ShapeType::Vertex && type2 == ShapeType::Edge {
                if let (Some(edge), Some(v)) = (Edge::wrap(root2.clone()), Vertex::wrap(root1.clone()))
                {
                    let (v1, v2) = edge_vertices(&edge);
                    let k = shape_map_key(&v.0);
                    let (k1, k2) = (
                        v1.map(|x| shape_map_key(&x.0)),
                        v2.map(|x| shape_map_key(&x.0)),
                    );
                    if Some(k.clone()) == k1 || Some(k) == k2 {
                        set_stored_continuity(&oe, &f1, &f2, GEOM_ABS_G1);
                    }
                }
            } else if type1 == ShapeType::Face && type2 == ShapeType::Edge {
                // `cxx:4093-4105`: the edge root is an edge of the face root.
                let k = shape_map_key(&root2);
                if explorer_shapes(&root1, ShapeType::Edge)
                    .iter()
                    .any(|e| shape_map_key(e) == k)
                {
                    set_stored_continuity(&oe, &f1, &f2, GEOM_ABS_G1);
                }
            } else if type1 == ShapeType::Edge && type2 == ShapeType::Face {
                // `cxx:4106-4117`.
                let k = shape_map_key(&root1);
                if explorer_shapes(&root2, ShapeType::Edge)
                    .iter()
                    .any(|e| shape_map_key(e) == k)
                {
                    set_stored_continuity(&oe, &f1, &f2, GEOM_ABS_G1);
                }
            } else if type1 == ShapeType::Face && type2 == ShapeType::Face {
                // `cxx:4119-4137`.
                let mut le: Vec<TopoShape> = Vec::new();
                let _ = _Tool::find_common_shapes_of_type(
                    &root1,
                    &root2,
                    ShapeType::Edge,
                    &mut le,
                );
                if le.len() == 1 {
                    if let Some(ed) = Edge::wrap(le[0].clone()) {
                        if self.my_analyse.has_ancestor(&ed.0) {
                            let li = self.my_analyse.type_of(&ed);
                            if li.len() == 1 && li[0].kind() == TypeOfConcavity::Tangential {
                                set_stored_continuity(&oe, &f1, &f2, GEOM_ABS_G1);
                            }
                        }
                    }
                }
            } else if type1 == ShapeType::Edge && type2 == ShapeType::Edge {
                // `cxx:4138-4156`.
                let mut lv: Vec<TopoShape> = Vec::new();
                let _ = _Tool::find_common_shapes_of_type(
                    &root1,
                    &root2,
                    ShapeType::Vertex,
                    &mut lv,
                );
                if lv.len() == 1 {
                    if let (Some(e_root), Some(v)) =
                        (Edge::wrap(root1.clone()), Vertex::wrap(lv[0].clone()))
                    {
                        let mut led_tg: Vec<TopoShape> = Vec::new();
                        self.my_analyse.tangent_edges(&e_root, &v, &mut led_tg);
                        let k = shape_map_key(&root2);
                        if led_tg.iter().any(|x| shape_map_key(x) == k) {
                            set_stored_continuity(&oe, &f1, &f2, GEOM_ABS_G1);
                        }
                    }
                }
            }
        }
    }
}
