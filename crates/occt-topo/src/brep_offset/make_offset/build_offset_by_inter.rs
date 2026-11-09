//! `BRepOffset_MakeOffset::BuildOffsetByInter` (`BRepOffset_MakeOffset.cxx:
//! 1277-1823`), the intersection join construction of the Complete mode.
//!
//! The method runs the `.cxx` blocks in order:
//!
//! 1. `:1333-1339` `MakeOffsetFaces(MapSF)`.
//! 2. `:1360-1388` `BRepOffset_Inter3d(AsDes, TopAbs_IN, myTol)` followed by
//!    `ConnexIntByInt(myFaceComp, ...)` (the parallel faces) and
//!    `ContextIntByInt(myFaces, ExtentContext, ...)`, where `ExtentContext` is
//!    `myOffset > 0` (`:1353-1358`).
//! 3. `:1390-1416` `IntersectEdges(aLFaces, ...)` over the faces of
//!    `myFaceComp` plus `myAnalyse.NewFaces()`.
//! 4. `:1417-1440` `TrimEdges(...)`; a `false` answer is
//!    `BRepOffset_CannotTrimEdges`.
//! 5. `:1443-1498` `GetEnlargedFaces` then `BRepOffset_Inter2d::Compute` on
//!    every enlarged face and every cap of `myFaces`, closed by
//!    `FuseVertices(aDMVV, AsDes, myImageVV)`.
//! 6. `:1505-1541` the loop rebuild: `BuildSplitsOfExtendedFaces` for
//!    Complete + Intersection + planar, `myMakeLoops.Build` otherwise.
//! 7. `:1550-1817` the history fill of `myInitOffsetFace`, `myImageOffset`,
//!    `myAsDes`, `myInitOffsetEdge`.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `analyzeProgress` (`:5410-5441`) and every `Message_ProgressScope` /
//!   `aSteps` (`:1286-1326`, `:1338`, `:1370`, `:1384`, `:1411`, `:1455-1458`,
//!   `:1473-1481`, `:1512`, `:1535`, `:1551-1554`, `:1765-1773`): the port has
//!   no progress machinery anywhere in `BRepOffset`, so no `!More()`
//!   `BRepOffset_UserBreak` arm can fire.
//!
//! Port guards, each replacing an OCCT exception or null-handle dereference:
//!
//! * `MapSF(FI)` (`:1560`, `:1739`) raises `Standard_NoSuchObject` on an
//!   unbound key; the port skips that face / edge.
//! * `TopoDS::Face(...)` (`:1562`, `:1467`, `:1487`) throws on a type mismatch;
//!   the port wraps through [`Face::wrap`] and skips the shape.
//! * `BRepAlgo_Image::Bind` (`:1566`, `:1603`, `:1650`, `:1705`, `:1750`) and
//!   `Add` (`:1616`, `:1744`) throw on an already-bound root; the port ignores
//!   the failure (`Err`), like every other `BRepOffset` module.

use std::collections::{HashMap, HashSet};

use occt_core::kernel::geomabs::JoinType;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::inter2d::{BRepOffsetInter2d, Dmvv};
use crate::brep_offset::{fuse_vertices, BRepOffsetInter3d, Error, MapSf, ShapeMap};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};

use super::intersect_trim::{get_enlarged_faces, intersect_edges, trim_edges};
use super::make_offset_faces::make_offset_faces;
use super::BRepOffsetMakeOffset;

/// `BuildOffsetByInter(theRange)` (`cxx:1277-1823`).
pub fn build_offset_by_inter(this: &mut BRepOffsetMakeOffset) {
    // `cxx:1333-1339`: the offset faces.
    let mut map_sf: MapSf = HashMap::new();
    make_offset_faces(this, &mut map_sf);
    if this.my_error != Error::NoError {
        return;
    }

    // `cxx:1344-1348`: `MES` (offset shape -> extended shape), `Build`
    // (initial shape -> the shape built by inter), `Failed`, and the `AsDes`
    // handle the 3D phase fills.
    let mut mes: ShapeMap = HashMap::new();
    let mut build: ShapeMap = HashMap::new();
    let mut failed: Vec<TopoShape> = Vec::new();

    // `cxx:1347`: `TopAbs_State Side = TopAbs_IN`.
    let side = State::In;

    // `cxx:1353-1358`.
    let extent_context = this.my_offset > 0.0;

    // `cxx:1360`: `BRepOffset_Inter3d Inter3(AsDes, Side, myTol)`.
    let mut inter3 = BRepOffsetInter3d::new(AsDes::new(), side, this.my_tol);

    // `cxx:1362-1369`: `Inter3.ConnexIntByInt(myFaceComp, MapSF, myAnalyse,
    // MES, Build, Failed, ..., myIsPlanar)`.
    inter3.connex_int_by_int(
        &this.my_face_comp.0,
        &map_sf,
        &this.my_analyse,
        &mut mes,
        &mut build,
        &mut failed,
        this.my_is_planar,
    );

    // `cxx:1375-1388`: `Inter3.ContextIntByInt(myFaces, ExtentContext, MapSF,
    // myAnalyse, MES, Build, Failed, ..., myIsPlanar)`.
    let context_faces: Vec<TopoShape> = this.my_faces.iter().cloned().collect();
    inter3.context_int_by_int(
        &context_faces,
        extent_context,
        &map_sf,
        &this.my_analyse,
        &mut mes,
        &mut build,
        &mut failed,
        this.my_is_planar,
    );

    // `cxx:1390-1399`: `aLFaces` = the faces of `myFaceComp` then the faces
    // `BRepOffset_Analyse` created.
    let mut a_l_faces: Vec<TopoShape> = explorer_shapes(&this.my_face_comp.0, ShapeType::Face);
    a_l_faces.extend(this.my_analyse.new_faces().iter().cloned());

    // The `Inter3` object owns the `AsDes` handle; the rest of the method is
    // the only consumer left.
    let mut as_des = inter3.into_as_des();

    // `cxx:1402-1416`: `IntersectEdges(aLFaces, MapSF, MES, Build, AsDes,
    // AsDes2d, ...)`, then `myError != BRepOffset_NoError` returns.
    let mut as_des2d = AsDes::new();
    intersect_edges(
        this,
        &a_l_faces,
        &map_sf,
        &mut mes,
        &build,
        &as_des,
        &mut as_des2d,
    );
    if this.my_error != Error::NoError {
        return;
    }

    // `cxx:1417-1420`: `anEdgesOrigins` (offset edge -> initial edges),
    // `NewEdges` (`IsSame`-unique) and `aETrimEInf` (trimmed -> not trimmed).
    let mut edges_origins: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut new_edges: HashSet<ShapeMapKey> = HashSet::new();
    let mut e_trim_e_inf: ShapeMap = HashMap::new();

    // `cxx:1429-1440`.
    if !trim_edges(
        this,
        &map_sf,
        &mes,
        &build,
        &mut as_des,
        &as_des2d,
        &mut new_edges,
        &mut e_trim_e_inf,
        &mut edges_origins,
    ) {
        this.my_error = Error::CannotTrimEdges;
        return;
    }

    // `cxx:1443-1452`: the 2D stage.
    let mut dmvv = Dmvv::new();
    let mut faces_origins: ShapeMap = HashMap::new();
    let mut imoe = BRepAlgoImage::new();
    let mut lfe: Vec<TopoShape> = Vec::new();
    get_enlarged_faces(
        &a_l_faces,
        &map_sf,
        &mes,
        &mut faces_origins,
        &mut imoe,
        &mut lfe,
    );

    // `cxx:1457-1477`: `BRepOffset_Inter2d::Compute` on every enlarged face.
    for nef_shape in &lfe {
        let nef = Face(nef_shape.clone());
        let a_curr_face_tol = BRepTool::face_tolerance(&nef);
        BRepOffsetInter2d::compute(
            &mut as_des,
            &nef,
            &new_edges,
            a_curr_face_tol,
            &this.my_edge_int_edges,
            &mut dmvv,
        );
    }

    // `cxx:1480-1497`: the same on every cap of `myFaces`.
    for i in 1..=this.my_faces.extent() {
        let Some(cork_shape) = this.my_faces.value(i).cloned() else {
            continue;
        };
        let cork = Face(cork_shape);
        let a_curr_face_tol = BRepTool::face_tolerance(&cork);
        BRepOffsetInter2d::compute(
            &mut as_des,
            &cork,
            &new_edges,
            a_curr_face_tol,
            &this.my_edge_int_edges,
            &mut dmvv,
        );
    }

    // `cxx:1499`: `FuseVertices(aDMVV, AsDes, myImageVV)`.
    let _ = fuse_vertices(&dmvv, &mut as_des, &mut this.my_image_vv);

    // `cxx:1504`: `aMFDone`.
    let mut mf_done: HashSet<ShapeMapKey> = HashSet::new();

    // `cxx:1506-1540`.
    if this.my_join == JoinType::Intersection && this.my_inter && this.my_is_planar {
        // `cxx:1508-1519`.
        super::build_offset_faces::build_splits_of_extended_faces(
            &lfe,
            &this.my_analyse,
            &as_des,
            &mut edges_origins,
            &mut faces_origins,
            &e_trim_e_inf,
            &mut imoe,
        );
        if this.my_error != Error::NoError {
            return;
        }

        // `cxx:1521-1526`.
        for a_s in &lfe {
            mf_done.insert(shape_map_key(a_s));
        }
    } else {
        // `cxx:1529-1534`: the `else` arm.
        this.my_make_loops
            .build(&lfe, &mut as_des, &mut imoe, &this.my_image_vv);
    }

    // `cxx:1550-1676`: history for the offset faces.
    for fi in &a_l_faces {
        let Some(bof) = map_sf.get(&shape_map_key(fi)) else {
            continue;
        };
        this.my_init_offset_face.set_root(fi);
        let Some(of) = bof.face().cloned() else {
            continue;
        };
        let mut of = of;

        if let Some(m) = mes.get(&shape_map_key(&of.0)).cloned() {
            // `cxx:1562-1563`: `OF = TopoDS::Face(MES(OF))`.
            let Some(mf) = Face::wrap(m) else {
                continue;
            };
            of = mf;

            if imoe.has_image(&of.0) {
                // `cxx:1565-1595`.
                let lofe = imoe.image(&of.0);
                let _ = this.my_init_offset_face.bind_list(fi, &lofe);
                for ofe in &lofe {
                    this.my_image_offset.set_root(ofe);
                    let mut view: HashSet<ShapeMapKey> = HashSet::new();
                    for coe in explorer_shapes(&ofe.oriented(Orientation::Forward), ShapeType::Edge)
                    {
                        this.my_as_des.add(ofe, &coe);
                        if view.insert(shape_map_key(&coe)) {
                            if !this.my_as_des.has_descendant(&coe) {
                                let (cv1, cv2) = edge_vertices(&Edge(coe.clone()));
                                if let Some(cv1) = cv1 {
                                    this.my_as_des
                                        .add(&coe, &cv1.0.oriented(Orientation::Forward));
                                }
                                if let Some(cv2) = cv2 {
                                    this.my_as_des
                                        .add(&coe, &cv2.0.oriented(Orientation::Reversed));
                                }
                            }
                        }
                    }
                }
            } else {
                // `cxx:1597-1648`.
                if mf_done.contains(&shape_map_key(&of.0)) {
                    continue;
                }
                let _ = this.my_init_offset_face.bind(fi, &of.0);
                this.my_image_offset.set_root(&of.0);
                let le = as_des.descendant(&of.0).to_vec();
                for oe in le {
                    if imoe.has_image(&oe) {
                        let loe = imoe.image(&oe);
                        for loe_i in loe {
                            // `cxx:1616-1619`: the image carries the edge's
                            // orientation.
                            let a_local_shape = loe_i.oriented(oe.orientation());
                            this.my_as_des.add(&of.0, &a_local_shape);
                            if !this.my_as_des.has_descendant(&a_local_shape) {
                                let (cv1, cv2) = edge_vertices(&Edge(a_local_shape.clone()));
                                if let Some(cv1) = cv1 {
                                    this.my_as_des
                                        .add(&a_local_shape, &cv1.0.oriented(Orientation::Forward));
                                }
                                if let Some(cv2) = cv2 {
                                    this.my_as_des.add(
                                        &a_local_shape,
                                        &cv2.0.oriented(Orientation::Reversed),
                                    );
                                }
                            }
                        }
                    } else {
                        // `cxx:1638-1643`.
                        this.my_as_des.add(&of.0, &oe);
                        let lv = as_des.descendant(&oe).to_vec();
                        this.my_as_des.add_list(&oe, &lv);
                    }
                }
            }
        } else {
            // `cxx:1647-1675`: no extended image; every edge of the offset
            // face is recorded with its endpoints.
            let _ = this.my_init_offset_face.bind(fi, &of.0);
            this.my_image_offset.set_root(&of.0);
            let mut view: HashSet<ShapeMapKey> = HashSet::new();
            for coe in explorer_shapes(&of.0.oriented(Orientation::Forward), ShapeType::Edge) {
                this.my_as_des.add(&of.0, &coe);
                if view.insert(shape_map_key(&coe)) {
                    if !this.my_as_des.has_descendant(&coe) {
                        let (cv1, cv2) = edge_vertices(&Edge(coe.clone()));
                        if let Some(cv1) = cv1 {
                            this.my_as_des
                                .add(&coe, &cv1.0.oriented(Orientation::Forward));
                        }
                        if let Some(cv2) = cv2 {
                            this.my_as_des
                                .add(&coe, &cv2.0.oriented(Orientation::Reversed));
                        }
                    }
                }
            }
        }
    }

    // `cxx:1679-1758`: history for the offset edges.
    let mut a_map_edges: HashSet<ShapeMapKey> = HashSet::new();
    for a_face_ref in &a_l_faces {
        for an_edge_ref in explorer_shapes(&a_face_ref.oriented(Orientation::Forward), ShapeType::Edge)
        {
            if !a_map_edges.insert(shape_map_key(&an_edge_ref)) {
                continue;
            }
            this.my_init_offset_edge.set_root(&an_edge_ref);
            if let Some(a_new_shape) = build.get(&shape_map_key(&an_edge_ref)).cloned() {
                if a_new_shape.shape_type() == ShapeType::Edge {
                    // `cxx:1695-1708`.
                    if imoe.has_image(&a_new_shape) {
                        let a_list_new_e = imoe.image(&a_new_shape);
                        let _ = this.my_init_offset_edge.bind_list(&an_edge_ref, &a_list_new_e);
                    } else {
                        let _ = this.my_init_offset_edge.bind(&an_edge_ref, &a_new_shape);
                    }
                } else {
                    // `cxx:1710-1734`: the intersection landed on several edges.
                    let mut a_list_new_edge: Vec<TopoShape> = Vec::new();
                    for a_res_edge in explorer_shapes(&a_new_shape, ShapeType::Edge) {
                        if imoe.has_image(&a_res_edge) {
                            a_list_new_edge.extend(imoe.image(&a_res_edge));
                        } else {
                            a_list_new_edge.push(a_res_edge.clone());
                        }
                    }
                    let _ = this
                        .my_init_offset_edge
                        .bind_list(&an_edge_ref, &a_list_new_edge);
                }
            } else {
                // `cxx:1737-1756`: a free boundary.
                let Some(bof) = map_sf.get(&shape_map_key(a_face_ref)) else {
                    continue;
                };
                let Some(mut a_new_edge) = bof.generated(&an_edge_ref) else {
                    continue;
                };
                if let Some(m) = mes.get(&shape_map_key(&a_new_edge)).cloned() {
                    a_new_edge = m;
                }
                if imoe.has_image(&a_new_edge) {
                    let a_list_new_e = imoe.image(&a_new_edge);
                    let _ = this.my_init_offset_edge.bind_list(&an_edge_ref, &a_list_new_e);
                } else {
                    let _ = this.my_init_offset_edge.bind(&an_edge_ref, &a_new_edge);
                }
            }
        }
    }

    // `cxx:1761-1817`: history for the caps.
    for i in 1..=this.my_faces.extent() {
        let Some(cork) = this.my_faces.value(i).cloned() else {
            continue;
        };
        let le = as_des.descendant(&cork).to_vec();
        for oe in le {
            if imoe.has_image(&oe) {
                let loe = imoe.image(&oe);
                for coe in loe {
                    this.my_as_des.add(&cork, &coe.oriented(oe.orientation()));
                    if !this.my_as_des.has_descendant(&coe) {
                        let (cv1, cv2) = edge_vertices(&Edge(coe.clone()));
                        if let Some(cv1) = cv1 {
                            this.my_as_des
                                .add(&coe, &cv1.0.oriented(Orientation::Forward));
                        }
                        if let Some(cv2) = cv2 {
                            this.my_as_des
                                .add(&coe, &cv2.0.oriented(Orientation::Reversed));
                        }
                    }
                }
            } else {
                this.my_as_des.add(&cork, &oe);
                if as_des.has_descendant(&oe) {
                    let lv = as_des.descendant(&oe).to_vec();
                    this.my_as_des.add_list(&oe, &lv);
                }
            }
        }
    }
}
