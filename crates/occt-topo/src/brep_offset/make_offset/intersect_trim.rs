//! `BRepOffset_MakeOffset`'s 2D stage: `IntersectEdges` (`cxx:4674-4747`),
//! the file-local `TrimEdges` (`:4752-4923`) / `TrimEdge` (`:4921-5033`) and
//! `GetEnlargedFaces` (`:5036-5056`).
//!
//! `IntersectEdges` is the driver of the 2D arm of `BuildOffsetByInter`
//! (`:1404-1411`): it runs `BRepOffset_Inter2d::ConnexIntByInt` on every face
//! of `myFaceComp` + `myAnalyse.NewFaces()`, then
//! `ConnexIntByIntInVert` on the faces that collected vertices, and finally
//! `FuseVertices`.
//!
//! `TrimEdges` restricts every new edge to its intersection parameters:
//!
//! * `:4783-4788` the faces of `myAnalyse.NewFaces()` join the face list and
//!   get an edge -> face ancestry map (`aDMEF`);
//! * `:4790-4802` per face, `NF` is `MapSF(FI).Face()`, replaced by
//!   `MES(NF)` when bound;
//! * `:4802-4806` the edge list is `MapShapes(FI, EDGE)` followed by
//!   `MapShapes(FI, VERTEX)` in one `IsSame`-unique map, so edges come first;
//! * `:4815-4867` a shape the 3D phase rebuilt (`Build(aS)`) is trimmed once
//!   per `NewEdges` entry; an edge arm trims the whole shape, a compound arm
//!   trims each piece only when `AsDes2d` has descendants for it and
//!   otherwise drops it from `AsDes`;
//! * `:4868-4915` a free edge is skipped when it is a vertex, or an edge of a
//!   generated face with a single face ancestor; otherwise its offset image
//!   (orientation carried over) is trimmed, or - when `MES` has no image -
//!   only extended by `BRepOffset_Inter2d::ExtentEdge` if it is a line; then
//!   it becomes a descendant of `NF`.
//!
//! `TrimEdge` (`:4920-5033`) is the per-edge trim: it splits the edge to get
//! a "source" copy, takes the parameter of every `AsDes2d` descendant vertex
//! (`FindParameter`, then projection when that fails), and when the two
//! extreme vertices differ it frees the edge, replaces its endpoints, writes
//! `[UMin, UMax]`, re-applies `SameParameter` and records the source edge in
//! `theETrimEInf`. The untrimmed arm drops a line edge from `AsDes`.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `Message_ProgressScope` (`:4687`, `:4689`, `:4720`, `:4752-5033`,
//!   `:4763`): the port has no progress machinery anywhere in `BRepOffset`.
//!
//! Port guards, each replacing an OCCT exception or null-handle dereference:
//!
//! * `theMapSF(FI)` (`:4791`, `:4880`, `:1739` in `BuildOffsetByInter`)
//!   raises `Standard_NoSuchObject` on an unbound key; the port fails the
//!   call (`TrimEdges` answers `false`, i.e. `BRepOffset_CannotTrimEdges`).
//! * `BOPTools_AlgoTools::MakeSplitEdge` (`:4939`) throws when the edge has no
//!   3D curve; the port's [`AlgoTools::make_split_edge`] answers `Err` and the
//!   trim fails the same way.
//! * `GeomAPI_ProjectPointOnCurve` with `NbPoints() == 0` is the `false`
//!   return of `FindParameter` (`:4967-4972`), kept.
//! * `GeomAdaptor_Curve(NE).GetType()` on an edge without a 3D curve (`:4905`,
//!   `:5014`) throws in OCCT; the port treats it as "not a line", which only
//!   skips the `ExtentEdge` / `AsDes::Remove` arms of a geometry-less edge.

use std::collections::{HashMap, HashSet};

use occt_core::precision::CONFUSION;
use occt_geom::geom_api::project_point_on_curve;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_offset::inter2d::{
    connex_int_by_int, connex_int_by_int_in_vert, extent_edge, fuse_vertices, ConnexByIntData,
    Dmvv,
};
use crate::brep_offset::{Error, MapSf, ShapeMap};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, map_shapes_and_unique_ancestors, shape_map_key,
    ShapeMapKey,
};

use super::BRepOffsetMakeOffset;

/// `IntersectEdges(theFaces, theMapSF, theMES, theBuild, theAsDes,
/// theAsDes2d, theRange)` (`cxx:4674-4749`).
#[allow(clippy::too_many_arguments)]
pub fn intersect_edges(
    this: &mut BRepOffsetMakeOffset,
    faces: &[TopoShape],
    map_sf: &MapSf,
    mes: &mut ShapeMap,
    build: &ShapeMap,
    as_des: &AsDes,
    as_des2d: &mut AsDes,
) {
    // `aDMVV` (`cxx:4684`) survives both passes and reaches `FuseVertices`.
    let mut dmvv = Dmvv::new();
    // `aMFV` (`cxx:4686`): the `IsSame`-unique `IndexedMap` of faces that
    // collected vertices.
    let mut amfv: Vec<TopoShape> = Vec::new();

    // `cxx:4689-4716`.
    for it in faces {
        let Some(a_f) = Face::wrap(it.clone()) else {
            continue;
        };
        let a_tol_f = BRepTool::face_tolerance(&a_f);
        let Some(ofi) = map_sf.get(&shape_map_key(&a_f.0)) else {
            // `theMapSF(aF)` throws `Standard_NoSuchObject`.
            this.my_error = Error::UnknownError;
            return;
        };
        let mut data = ConnexByIntData {
            faces_with_verts: &mut amfv,
            image_vv: &mut this.my_image_vv,
            edge_int_edges: &mut this.my_edge_int_edges,
            dmvv: &mut dmvv,
        };
        let done = connex_int_by_int(
            &a_f,
            ofi,
            mes,
            build,
            as_des,
            as_des2d,
            this.my_offset,
            a_tol_f,
            &this.my_analyse,
            &mut data,
        );
        if !done {
            this.my_error = Error::CannotExtentEdge;
            return;
        }
    }

    // `cxx:4718-4741`: the faces that collected vertices.
    let a_nb_f = amfv.len();
    for i in 0..a_nb_f {
        let Some(a_f) = Face::wrap(amfv[i].clone()) else {
            continue;
        };
        let a_tol_f = BRepTool::face_tolerance(&a_f);
        let Some(ofi) = map_sf.get(&shape_map_key(&a_f.0)) else {
            continue;
        };
        connex_int_by_int_in_vert(
            &a_f,
            ofi,
            mes,
            build,
            as_des,
            as_des2d,
            a_tol_f,
            &this.my_analyse,
            &mut dmvv,
        );
    }

    // `cxx:4743-4748`.
    if !fuse_vertices(&dmvv, as_des2d, &mut this.my_image_vv) {
        this.my_error = Error::CannotFuseVertices;
    }
}

/// `TrimEdges(theShape, theOffset, Analyse, theMapSF, theMES, theBuild,
/// theAsDes, theAsDes2d, theNewEdges, theETrimEInf, theEdgesOrigins)`
/// (`cxx:4752-4918`).
///
/// `false` is the reference's `return false` (`:4835`, `:4854`, `:4898`),
/// which `BuildOffsetByInter` reports as `BRepOffset_CannotTrimEdges`.
#[allow(clippy::too_many_arguments)]
pub fn trim_edges(
    this: &mut BRepOffsetMakeOffset,
    map_sf: &MapSf,
    mes: &ShapeMap,
    build: &ShapeMap,
    as_des: &mut AsDes,
    as_des2d: &AsDes,
    new_edges: &mut HashSet<ShapeMapKey>,
    e_trim_e_inf: &mut ShapeMap,
    edges_origins: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> bool {
    // `cxx:4774-4778`: the faces of `theShape`.
    let mut a_l_faces: Vec<TopoShape> =
        explorer_shapes(&this.my_face_comp.0, ShapeType::Face);

    // `cxx:4779-4788`: the generated faces and their edge -> face ancestry.
    let mut a_mf_generated: HashSet<ShapeMapKey> = HashSet::new();
    let mut admef: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    for a_fg in this.my_analyse.new_faces() {
        a_l_faces.push(a_fg.clone());
        a_mf_generated.insert(shape_map_key(a_fg));
        // `TopExp::MapShapesAndUniqueAncestors(aFG, EDGE, FACE, aDMEF)`: the
        // map accumulates across the generated faces, so an edge shared by two
        // of them keeps both ancestors (`TopExp.cxx:124-177`).
        let m = map_shapes_and_unique_ancestors(a_fg, ShapeType::Edge, ShapeType::Face, false);
        for (k, v) in m {
            let list = admef.entry(k).or_default();
            for anc in v {
                if !list.iter().any(|x| is_same(x, &anc)) {
                    list.push(anc);
                }
            }
        }
    }

    // `cxx:4789-4917`.
    for face_shape in &a_l_faces {
        let Some(fi) = Face::wrap(face_shape.clone()) else {
            continue;
        };
        // `NF = theMapSF(FI).Face()` (`:4791`).
        let Some(of) = map_sf.get(&shape_map_key(&fi.0)) else {
            return false;
        };
        let Some(of_face) = of.face() else {
            return false;
        };
        let mut nf = Face(of_face.0.clone());
        if let Some(m) = mes.get(&shape_map_key(&nf.0)) {
            nf = Face(m.clone());
        }

        // `cxx:4797-4806`: `View` (dedup) over `VEmap` = the face's edges then
        // its vertices, both `IsSame`-unique.
        let mut view: HashSet<ShapeMapKey> = HashSet::new();
        let mut ve_map: Vec<TopoShape> = Vec::new();
        let mut ve_seen: HashSet<ShapeMapKey> = HashSet::new();
        for s in explorer_shapes(&fi.0.oriented(Orientation::Forward), ShapeType::Edge) {
            if ve_seen.insert(shape_map_key(&s)) {
                ve_map.push(s);
            }
        }
        for s in explorer_shapes(&fi.0.oriented(Orientation::Forward), ShapeType::Vertex) {
            if ve_seen.insert(shape_map_key(&s)) {
                ve_map.push(s);
            }
        }

        for a_s in ve_map {
            if !view.insert(shape_map_key(&a_s)) {
                continue;
            }
            // `cxx:4815-4867`: a shape the 3D phase rebuilt.
            if let Some(ne) = build.get(&shape_map_key(&a_s)).cloned() {
                // `cxx:4817-4827`: keep the connection to the original edges.
                for nec in explorer_shapes(&ne, ShapeType::Edge) {
                    edges_origins
                        .entry(shape_map_key(&nec))
                        .or_default()
                        .push(a_s.clone());
                }
                if ne.shape_type() == ShapeType::Edge {
                    // `cxx:4829-4838`.
                    if new_edges.insert(shape_map_key(&ne)) {
                        if !trim_edge(&mut Edge(ne.clone()), as_des2d, as_des, e_trim_e_inf) {
                            return false;
                        }
                    }
                } else {
                    // `cxx:4839-4867`: the intersection lands on several edges.
                    for nec in explorer_shapes(&ne, ShapeType::Edge) {
                        if !new_edges.insert(shape_map_key(&nec)) {
                            continue;
                        }
                        if !as_des2d.descendant(&nec).is_empty() {
                            if !trim_edge(&mut Edge(nec.clone()), as_des2d, as_des, e_trim_e_inf) {
                                return false;
                            }
                        } else if as_des.has_ascendant(&nec) {
                            let _ = as_des.remove(&nec);
                        }
                    }
                }
                continue;
            }

            // `cxx:4868-4878`: the free arm.
            if a_s.shape_type() != ShapeType::Edge {
                continue;
            }
            if a_mf_generated.contains(&shape_map_key(&fi.0))
                && admef
                    .get(&shape_map_key(&a_s))
                    .map(|v| v.len())
                    .unwrap_or(0)
                    == 1
            {
                continue;
            }

            // `cxx:4880-4888`: `NE = theMapSF(FI).Generated(aS)`, orientation
            // carried over (`OCC4455`).
            let Some(of) = map_sf.get(&shape_map_key(&fi.0)) else {
                return false;
            };
            let Some(mut ne) = of.generated(&a_s) else {
                continue;
            };
            ne.set_orientation(a_s.orientation());
            edges_origins
                .entry(shape_map_key(&ne))
                .or_default()
                .push(a_s.clone());

            // `cxx:4890-4903`: trim the extended image when `MES` has one.
            if let Some(m) = mes.get(&shape_map_key(&ne)).cloned() {
                ne = m;
                ne.set_orientation(a_s.orientation());
                if new_edges.insert(shape_map_key(&ne)) {
                    if !trim_edge(&mut Edge(ne.clone()), as_des2d, as_des, e_trim_e_inf) {
                        return false;
                    }
                }
            } else {
                // `cxx:4904-4913`: only a line edge is extended here.
                let an_edge = Edge(ne.clone());
                let is_line = BRepTool::edge_curve(&an_edge)
                    .map(|c| curve_kind(c.as_ref()) == GeomAbsCurveType::Line)
                    .unwrap_or(false);
                if is_line {
                    if let Some(new_edge) = extent_edge(&an_edge, this.my_offset) {
                        e_trim_e_inf
                            .insert(shape_map_key(&an_edge.0), new_edge.0.clone());
                    }
                }
            }
            // `cxx:4915`.
            as_des.add(&nf.0, &ne);
        }
    }
    true
}

/// `TrimEdge(NE, AsDes2d, AsDes, theETrimEInf)` (`cxx:4920-5033`).
///
/// `NE` is rewritten in place (`Free`, endpoint `Remove`/`Add`, `Range`,
/// orientation) exactly like the reference's `TopoDS_Edge&`.
///
/// `false` is one of the reference's `return false` statements (`:4970`,
/// `:4988`), plus the port guard for a failing `MakeSplitEdge`.
pub fn trim_edge(
    ne: &mut Edge,
    as_des2d: &AsDes,
    as_des: &mut AsDes,
    e_trim_e_inf: &mut ShapeMap,
) -> bool {
    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();

    // `cxx:4936-4939`: the vertices, the range and the source edge.
    let (v1, v2) = edge_vertices(ne);
    let (a_t1, a_t2) = BRepTool::edge_parameters(ne);
    let v1_shape = v1.as_ref().map(|v| v.0.clone());
    let v2_shape = v2.as_ref().map(|v| v.0.clone());
    let Ok(a_source_edge) = AlgoTools::make_split_edge(
        ne,
        v1_shape.as_ref(),
        a_t1,
        v2_shape.as_ref(),
        a_t2,
    ) else {
        return false;
    };

    // `cxx:4943-4947`.
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut first = v1.clone();
    let mut last = v2.clone();

    let le: Vec<TopoShape> = as_des2d.descendant(&ne.0).to_vec();
    let mut b_trim = false;
    if le.len() > 1 {
        // `cxx:4956-4979`.
        for v_shape in &le {
            let mut v = Vertex(v_shape.clone());
            if ne.orientation() == Orientation::Reversed {
                v.0.reverse();
            }
            let u = match BRepTool::parameter_of_vertex_checked(&ne.0, &v.0) {
                Some(u) => u,
                None => {
                    // `cxx:4963-4972`: project the vertex on the 3D curve.
                    let Some(the_curve) = BRepTool::edge_curve(ne) else {
                        return false;
                    };
                    let the_point = BRepTool::vertex_point_world(&v);
                    match project_point_on_curve(the_curve.as_ref(), &the_point, CONFUSION) {
                        Some(p) => p.parameter,
                        None => return false,
                    }
                }
            };
            if u < u_min {
                u_min = u;
                first = Some(v.clone());
            }
            if u > u_max {
                u_max = u;
                last = Some(v.clone());
            }
        }

        // `cxx:4985-4991`.
        let (Some(first_v), Some(last_v)) = (first.as_ref(), last.as_ref()) else {
            return false;
        };
        if !is_same(&first_v.0, &last_v.0) {
            // `cxx:4992-5010`.
            ne.0.set_free(true);
            let or = ne.orientation();
            ne.0.set_orientation(Orientation::Forward);
            let (vf, vl) = edge_vertices(ne);
            if let Some(vf) = vf.as_ref() {
                builder.remove(&mut ne.0, &vf.0);
            }
            if let Some(vl) = vl.as_ref() {
                builder.remove(&mut ne.0, &vl.0);
            }
            builder.add(&mut ne.0, &first_v.0.oriented(Orientation::Forward));
            builder.add(&mut ne.0, &last_v.0.oriented(Orientation::Reversed));
            reg.set_edge_range(&ne.0, u_min, u_max);
            ne.0.set_orientation(or);
            as_des.add(&ne.0, &first_v.0.oriented(Orientation::Forward));
            as_des.add(&ne.0, &last_v.0.oriented(Orientation::Reversed));
            same_parameter_edge_inplace(ne, CONFUSION);
            b_trim = true;
        }
    }

    // `cxx:5012-5030`.
    if !b_trim {
        let is_line = BRepTool::edge_curve(ne)
            .map(|c| curve_kind(c.as_ref()) == GeomAbsCurveType::Line)
            .unwrap_or(false);
        if is_line {
            if as_des.has_ascendant(&ne.0) {
                let _ = as_des.remove(&ne.0);
            }
        }
    } else if !e_trim_e_inf.contains_key(&shape_map_key(&ne.0)) {
        e_trim_e_inf.insert(shape_map_key(&ne.0), a_source_edge.0.clone());
    }
    true
}

/// `GetEnlargedFaces(theFaces, theMapSF, theMES, theFacesOrigins, theImage,
/// theLSF)` (`cxx:5036-5056`): the enlarged face of every face whose offset
/// image `MES` rebuilt.
pub fn get_enlarged_faces(
    faces: &[TopoShape],
    map_sf: &MapSf,
    mes: &ShapeMap,
    faces_origins: &mut ShapeMap,
    image: &mut BRepAlgoImage,
    lsf: &mut Vec<TopoShape>,
) {
    for fi in faces {
        let Some(ofi) = map_sf.get(&shape_map_key(fi)) else {
            continue;
        };
        let Some(of_face) = ofi.face() else {
            continue;
        };
        if let Some(m) = mes.get(&shape_map_key(&of_face.0)).cloned() {
            // `cxx:5050-5054`.
            lsf.push(m.clone());
            image.set_root(&m);
            faces_origins.insert(shape_map_key(&m), fi.clone());
        }
    }
}

/// The `myEdgeIntEdges` map type (`NCollection_DataMap<TopoDS_Shape,
/// NCollection_List<TopoDS_Shape>>`).
pub type EdgeIntEdges = HashMap<ShapeMapKey, Vec<TopoShape>>;

/// The `theAnEdgesOrigins` value type: the initial edges an offset edge came
/// from (`NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`).
pub type EdgesOrigins = HashMap<ShapeMapKey, Vec<TopoShape>>;
