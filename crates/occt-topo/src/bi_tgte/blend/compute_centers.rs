//! `BiTgte_Blend::ComputeCenters` (`BiTgte_Blend.cxx:1286-1757`): build the
//! offset of every analysed face / free-border edge, intersect them, and cut
//! the resulting centre edges into loops.
//!
//! Control flow, all from the `.cxx`:
//!
//! * `:1288-1291` `myAnalyse.Perform(myShape, 2 * asin(myTol / |myRadius / 2|))`.
//! * `:1295-1297` the faces touched by the caps (`Touched`, a no-op in the
//!   reference).
//! * `:1301-1308` `BRepOffset_Inter3d(myAsDes, Side, myTol)` with
//!   `Side = myRadius < 0 ? TopAbs_OUT : TopAbs_IN`.
//! * `:1323-1400` the outer `JenRajoute` loop: locate the faces connected to
//!   `myEdges` (`:1335-1369`), then build one `BRepOffset_Offset` per new face
//!   (`:1382-1403`), enlarging the face when it is touched by a cap.
//! * `:1405-1448` fill `EdgeTgt`, the `Created` map handed to the next
//!   `Init(Face, Radius, Created)`, with the images of the tangent edges and of
//!   their vertices.
//! * `:1450-1475` a reversed / degenerated offset is dropped; otherwise
//!   `myInitOffsetFace` and `MapSBox` are updated and [`Blend::intersect`] runs.
//! * `:1480-1572` tubes on the convex edges of `Co` (the compound of all
//!   analysed faces): `BRepOffset_Offset(E, EOn1, EOn2, myRadius, E1f, E1l)`,
//!   then [`Blend::intersect`] again; `JenRajoute` is set when an intersection
//!   reached a restriction.
//! * `:1578-1580` `myEdges = Inter.NewEdges()`.
//! * `:1587-1660` `BRepOffset_Inter2d::Compute` on every generated face, after
//!   `myAsDes->Add` of the restrictions that do not disappear.
//! * `:1663-1707` the same for the tube faces of `myMapSF` (faces whose key is
//!   not a face), with every wire edge added to `myAsDes` first.
//! * `:1708-1717` `FuseVertices` then `BRepOffset_MakeLoops::Build`.
//! * `:1726-1750` the unwound faces of `myImageOffset` are added to
//!   `myResult`, skipping the faces that are their own image.
//!
//! Deviations:
//! - `myAncestors` is [`AncestorMap`], and `TopExp::MapShapesAndAncestors` is
//!   called twice for the `Co` map of `:1494-1497` (edges -> faces and
//!   vertices -> edges) because the port's map holds one inner type; the two
//!   key sets are disjoint (edges vs vertices), so every lookup is unchanged.
//! - `EdgeTgt` cannot store a null image (`OF1.Generated(Cur)` answering a null
//!   shape): the port skips that binding instead of binding a null shape.
//! - `BRepTools::Update` has no counterpart (no cached UV points).

use std::collections::{HashMap, HashSet};

use occt_core::bnd::box3d::BndBox;
use occt_core::kernel::geomabs::{JoinType, Shape};

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_bnd_lib;
use crate::brep_offset::ShapeIndexedMap;
use crate::brep_offset::{EnLargeFaceOptions, ShapeMap};
use crate::brep_offset::{
    fuse_vertices, BRepOffsetInter2d, BRepOffsetInter3d, BRepOffsetMakeLoops, BRepOffsetOffset,
    BRepOffsetStatus, BRepOffsetTool, Dmvv,
};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, shape_map_key, AncestorMap, ShapeMapKey,
};

use super::statics::touched;
use super::Blend;

impl Blend {
    /// `void BiTgte_Blend::ComputeCenters()` (`cxx:1286-1757`).
    pub(super) fn compute_centers(&mut self) {
        let builder = TopoBuilder::new();
        let Some(shape) = self.my_shape.clone() else {
            return;
        };

        // `cxx:1290-1291`.
        let tol_angle = 2.0 * (self.my_tol / (self.my_radius * 0.5).abs()).asin();
        self.my_analyse
            .perform(&shape, tol_angle)
            .expect("BRepOffset_Analyse::Perform");

        // `cxx:1295-1297`.
        let mut touched_by_cork: HashSet<ShapeMapKey> = HashSet::new();
        touched(
            &self.my_analyse,
            &self.my_stop_faces,
            &shape,
            &mut touched_by_cork,
        );

        // `cxx:1301-1305`.
        let side = if self.my_radius < 0.0 {
            State::Out
        } else {
            State::In
        };
        let mut inter =
            BRepOffsetInter3d::new(std::mem::take(&mut self.my_as_des), side, self.my_tol);

        // `cxx:1309-1317`.
        let mut map_sbox: HashMap<ShapeMapKey, BndBox> = HashMap::new();
        let mut co = builder.make_shape(ShapeType::Compound);

        // `cxx:1320-1576`.
        let mut jen_rajoute = true;
        while jen_rajoute {
            jen_rajoute = false;
            let mut fini = false;
            let mut edge_tgt: ShapeMap = ShapeMap::new();

            while !fini {
                fini = true;

                // `cxx:1335-1369`: locate in `myFaces` the faces connected to
                // `myEdges`.
                let mut i = 0;
                while i < self.my_edges.extent() {
                    let e = Edge(self.my_edges.value(i).expect("myEdges(i)").clone());
                    i += 1;
                    if BRepTool::is_degenerated(&e) {
                        continue;
                    }
                    let l = self
                        .my_ancestors
                        .find_from_key(&e.0)
                        .unwrap_or(&[])
                        .to_vec();
                    if l.len() == 1 {
                        // `cxx:1349-1354`: a free border on which the ball
                        // should roll.
                        self.my_faces.add(&e.0);
                        self.my_stop_faces.insert(shape_map_key(&e.0));
                    } else {
                        for sh in l {
                            if !self.my_stop_faces.contains(&shape_map_key(&sh)) {
                                self.my_faces.add(&sh);
                            }
                        }
                    }
                }
                self.my_edges = ShapeIndexedMap::new();

                // `cxx:1375-1476`.
                let mut i = 0;
                while i < self.my_faces.extent() {
                    let as_ = self.my_faces.value(i).expect("myFaces(i)").clone();
                    i += 1;

                    if self.my_map_sf.is_bound(&as_) {
                        continue;
                    }

                    let mut of1 = BRepOffsetOffset::new();

                    if as_.shape_type() == ShapeType::Face {
                        let f = Face(as_.clone());
                        if touched_by_cork.contains(&shape_map_key(&as_)) {
                            // `cxx:1390-1393`: `EnLargeFace(F, BigF, true)`.
                            let opts = EnLargeFaceOptions {
                                can_extent_surface: true,
                                ..EnLargeFaceOptions::default()
                            };
                            let (big_f, _) =
                                BRepOffsetTool::enlarge_face(&f, &opts).expect("EnLargeFace");
                            of1.init_face(&big_f, self.my_radius, &edge_tgt, true, JoinType::Arc)
                                .expect("BRepOffset_Offset::Init(EnlargedFace)");
                        } else {
                            of1.init_face(&f, self.my_radius, &edge_tgt, true, JoinType::Arc)
                                .expect("BRepOffset_Offset::Init(Face)");
                        }
                    } else {
                        // `cxx:1400-1403`: a free border edge on which the ball
                        // rolls.
                        of1.init_edge(&Edge(as_.clone()), self.my_radius)
                            .expect("BRepOffset_Offset::Init(Edge)");
                    }

                    // `cxx:1405-1448`: increment the map of created tangents.
                    let let_ = if as_.shape_type() == ShapeType::Face {
                        self.my_analyse
                            .edges_of_face(&Face(as_.clone()), TypeOfConcavity::Tangential)
                    } else {
                        Vec::new()
                    };
                    for cur in let_ {
                        if edge_tgt.contains_key(&shape_map_key(&cur.0)) {
                            continue;
                        }
                        let Some(gen) = of1.generated(&cur.0) else {
                            continue;
                        };
                        let ote = Edge(gen.clone());
                        edge_tgt.insert(shape_map_key(&cur.0), gen);
                        let (v1, v2) = edge_vertices(&cur);
                        let (ov1, ov2) = edge_vertices(&ote);
                        if let Some(v1) = v1 {
                            if !edge_tgt.contains_key(&shape_map_key(&v1.0)) {
                                let le = self
                                    .my_analyse
                                    .edges_of_vertex(&v1, TypeOfConcavity::Tangential);
                                let la = self.my_analyse.ancestors(&v1.0);
                                if le.len() == la.len() {
                                    if let Some(ov1) = ov1.clone() {
                                        edge_tgt.insert(shape_map_key(&v1.0), ov1.0);
                                    }
                                }
                            }
                        }
                        if let Some(v2) = v2 {
                            if !edge_tgt.contains_key(&shape_map_key(&v2.0)) {
                                let le = self
                                    .my_analyse
                                    .edges_of_vertex(&v2, TypeOfConcavity::Tangential);
                                let la = self.my_analyse.ancestors(&v2.0);
                                if le.len() == la.len() {
                                    if let Some(ov2) = ov2.clone() {
                                        edge_tgt.insert(shape_map_key(&v2.0), ov2.0);
                                    }
                                }
                            }
                        }
                    }

                    // `cxx:1450-1453`.
                    if matches!(
                        of1.status(),
                        BRepOffsetStatus::Reversed | BRepOffsetStatus::Degenerated
                    ) {
                        continue;
                    }

                    // `cxx:1455-1463`.
                    let Some(f1) = of1.face().cloned() else {
                        continue;
                    };
                    self.my_init_offset_face.set_root(&as_);
                    let _ = self.my_init_offset_face.bind(&as_, &f1.0);
                    let mut box1 = BndBox::new();
                    brep_bnd_lib::add(&f1.0, &mut box1, false);
                    map_sbox.insert(shape_map_key(&f1.0), box1);

                    // `cxx:1468-1473`.
                    fini = !self.intersect(&as_, &f1, &map_sbox, &of1, &mut inter);
                    if as_.shape_type() == ShapeType::Face {
                        builder.add(&mut co, &as_);
                    }
                    self.my_map_sf.bind(&as_, of1);
                }
            }

            // `cxx:1487-1490`.
            let ot = if self.my_radius < 0.0 {
                TypeOfConcavity::Concave
            } else {
                TypeOfConcavity::Convex
            };

            // `cxx:1494-1497`: `TopExp::MapShapesAndAncestors` into one map,
            // edges first then vertices. Only the edge entries are read here
            // (`:1507`); the vertex entries are never queried, so the port
            // builds the edge -> face map alone.
            let map_ef = AncestorMap::map_shapes_and_ancestors(&co, ShapeType::Edge, ShapeType::Face);

            // `cxx:1499-1574`.
            for occ in explorer_shapes(&co, ShapeType::Edge) {
                let e = Edge(occ.clone());
                if self.my_map_sf.is_bound(&e.0) {
                    continue;
                }
                let Some(anc) = map_ef.find_from_key(&e.0) else {
                    continue;
                };
                if anc.len() != 2 {
                    continue;
                }
                let l = self.my_analyse.type_of(&e);
                if l.is_empty() || l[0].kind() != ot {
                    continue;
                }

                let e_on1 = self
                    .my_map_sf
                    .get(&anc[0])
                    .and_then(|of| of.generated(&e.0));
                let e_on2 = self
                    .my_map_sf
                    .get(&anc[1])
                    .and_then(|of| of.generated(&e.0));
                let (Some(e_on1), Some(e_on2)) = (e_on1, e_on2) else {
                    continue;
                };
                let mut e_on1 = Edge(e_on1);
                let mut e_on2 = Edge(e_on2);

                // `cxx:1520-1552`: the tangent edges of both end vertices.
                let (v1f, v1l) = edge_vertices(&e);
                let mut e1f: Option<Edge> = None;
                let mut e1l: Option<Edge> = None;
                if let Some(v1f) = v1f.as_ref() {
                    let mut tang_e: Vec<TopoShape> = Vec::new();
                    self.my_analyse.tangent_edges(&e, v1f, &mut tang_e);
                    for t in tang_e {
                        if self.my_map_sf.is_bound(&t) {
                            if let Some(gen) = self.my_map_sf.get(&t).and_then(|of| of.generated(&v1f.0)) {
                                e1f = Some(Edge(gen));
                                break;
                            }
                        }
                    }
                }
                if let Some(v1l) = v1l.as_ref() {
                    let mut tang_e: Vec<TopoShape> = Vec::new();
                    self.my_analyse.tangent_edges(&e, v1l, &mut tang_e);
                    for t in tang_e {
                        if self.my_map_sf.is_bound(&t) {
                            if let Some(gen) = self.my_map_sf.get(&t).and_then(|of| of.generated(&v1l.0)) {
                                e1l = Some(Edge(gen));
                                break;
                            }
                        }
                    }
                }

                // `cxx:1553-1554`: the `FirstEdge` / `LastEdge` overload with
                // the header defaults `Polynomial = false`, `Tol = 1.0e-4`,
                // `Conti = GeomAbs_C1`.
                let mut of1 = BRepOffsetOffset::new();
                of1.init_pipe_first_last(
                    &e,
                    &mut e_on1,
                    &mut e_on2,
                    self.my_radius,
                    e1f.as_ref(),
                    e1l.as_ref(),
                    false,
                    1.0e-4,
                    Shape::C1,
                )
                .expect("BRepOffset_Offset::Init(Path, Edge1, Edge2, First, Last)");

                // `cxx:1555-1569`.
                let Some(f1) = of1.face().cloned() else {
                    continue;
                };
                self.my_init_offset_face.set_root(&e.0);
                let _ = self.my_init_offset_face.bind(&e.0, &f1.0);
                let mut box1 = BndBox::new();
                brep_bnd_lib::add(&f1.0, &mut box1, false);
                map_sbox.insert(shape_map_key(&f1.0), box1);
                let is_on_rest = self.intersect(&e.0, &f1, &map_sbox, &of1, &mut inter);
                jen_rajoute = jen_rajoute || is_on_rest;
                self.my_map_sf.bind(&e.0, of1);
            }
        }

        // `cxx:1578-1580`.
        self.my_edges = ShapeIndexedMap::new();
        self.my_edges = inter.new_edges().clone();

        // `cxx:1587-1590`: the concavity whose edges do *not* disappear.
        let ot = if self.my_radius < 0.0 {
            TypeOfConcavity::Convex
        } else {
            TypeOfConcavity::Concave
        };

        // `cxx:1593-1660`.
        let mut dm_vv = Dmvv::new();
        let mut lof: Vec<TopoShape> = Vec::new();
        let empty_map: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        let mut i = 0;
        while i < self.my_faces.extent() {
            let cur_s = self.my_faces.value(i).expect("myFaces(i)").clone();
            i += 1;

            // `cxx:1601-1604`: an undesirable tube on a free border.
            if self.my_stop_faces.contains(&shape_map_key(&cur_s)) {
                continue;
            }
            if !self.my_map_sf.is_bound(&cur_s) {
                continue; // inverted or degenerated
            }
            let Some(cur_of) = self.my_map_sf.get(&cur_s).and_then(|of| of.face().cloned()) else {
                continue;
            };
            lof.push(cur_of.0.clone());

            if cur_s.shape_type() == ShapeType::Face {
                let cur_f = Face(cur_s.clone());
                let cur_f_fwd = cur_f.0.oriented(Orientation::Forward);
                for expe in explorer_shapes(&cur_f_fwd, ShapeType::Edge) {
                    let cur_e = Edge(expe.clone());
                    let l = self.my_analyse.type_of(&cur_e);
                    // `cxx:1628-1647`: the edges that will disappear are only
                    // set when one of their ancestors leaves `myFaces`.
                    let set = if !l.is_empty() && l[0].kind() != ot {
                        true
                    } else {
                        let lanc = self.my_analyse.ancestors(&cur_e.0);
                        let first = lanc.first();
                        let last = lanc.last();
                        let contains = |f: Option<&TopoShape>| {
                            f.map(|f| self.my_faces.contains(f)).unwrap_or(false)
                        };
                        let stop = |f: Option<&TopoShape>| {
                            f.map(|f| self.my_stop_faces.contains(&shape_map_key(f)))
                                .unwrap_or(false)
                        };
                        !contains(first) || !contains(last) || stop(first) || stop(last)
                    };
                    if !set {
                        continue;
                    }
                    if let Some(cur_oe) = self
                        .my_map_sf
                        .get(&cur_f.0)
                        .and_then(|of| of.generated(&cur_e.0))
                    {
                        inter
                            .as_des_mut()
                            .add(&cur_of.0, &cur_oe.oriented(cur_e.0.orientation()));
                    }
                }

                let new_edges: HashSet<ShapeMapKey> =
                    self.my_edges.iter().map(shape_map_key).collect();
                BRepOffsetInter2d::compute(
                    inter.as_des_mut(),
                    &cur_f,
                    &new_edges,
                    self.my_tol,
                    &empty_map,
                    &mut dm_vv,
                );
            }
        }

        // `cxx:1663-1707`: also intersect the tubes built on free borders.
        let entries: Vec<(TopoShape, BRepOffsetOffset)> = self
            .my_map_sf
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (cur_s, of) in entries {
            if cur_s.shape_type() == ShapeType::Face {
                continue;
            }
            let Some(cur_of) = of.face().cloned() else {
                continue;
            };
            // `cxx:1679-1682`: no unwinding by tubes on a free border.
            if self.my_stop_faces.contains(&shape_map_key(&cur_s)) {
                continue;
            }
            lof.push(cur_of.0.clone());

            // `cxx:1688-1697`: `myAsDes->Add(CurOF, CurOE)` for every edge.
            let cur_of_fwd = cur_of.0.oriented(Orientation::Forward);
            for expe in explorer_shapes(&cur_of_fwd, ShapeType::Edge) {
                inter.as_des_mut().add(&cur_of.0, &expe);
            }

            let new_edges: HashSet<ShapeMapKey> =
                self.my_edges.iter().map(shape_map_key).collect();
            BRepOffsetInter2d::compute(
                inter.as_des_mut(),
                &cur_of,
                &new_edges,
                self.my_tol,
                &empty_map,
                &mut dm_vv,
            );
        }

        // `cxx:1708-1717`.
        let mut an_empty_image = BRepAlgoImage::new();
        fuse_vertices(&dm_vv, inter.as_des_mut(), &mut an_empty_image);
        let mut make_loops = BRepOffsetMakeLoops::new();
        make_loops.build(
            &lof,
            inter.as_des_mut(),
            &mut self.my_image_offset,
            &an_empty_image,
        );

        // `cxx:1726-1751`.
        self.my_result = Some(builder.make_shape(ShapeType::Compound));
        for cur_lof in lof.iter() {
            if !self.my_image_offset.has_image(cur_lof) {
                continue;
            }
            let mut lim: Vec<TopoShape> = Vec::new();
            self.my_image_offset.last_image(cur_lof, &mut lim);
            for cur_lim in lim {
                // `:1741-1747`: a face that is its own image is not added.
                if is_same(&cur_lim, cur_lof) {
                    break;
                }
                if let Some(result) = self.my_result.as_mut() {
                    builder.add(result, &cur_lim);
                }
            }
        }

        // `cxx:1304-1308`: hand the shared `AsDes` back to the blend.
        self.my_as_des = inter.into_as_des();
    }
}



