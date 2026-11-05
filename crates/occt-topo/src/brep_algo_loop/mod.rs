//! `BRepAlgo_Loop` — builds the loops (wires) of a set of edges on one face.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_Loop.cxx:1-1151` and
//! `.hxx:32-112` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepAlgo_Loop()` | `cxx:57-60` | [`BRepAlgoLoop::new`] |
//! | `Init` | `cxx:64-74` | [`BRepAlgoLoop::init`] |
//! | `AddEdge` | `cxx:126-130` | [`BRepAlgoLoop::add_edge`] |
//! | `AddConstEdge` | `cxx:134-137` | [`BRepAlgoLoop::add_const_edge`] |
//! | `AddConstEdges` | `cxx:140-147` | [`BRepAlgoLoop::add_const_edges`] |
//! | `SetImageVV` | `cxx:151-154` | [`BRepAlgoLoop::set_image_vv`] |
//! | `Perform` | `cxx:581-787` | [`BRepAlgoLoop::perform`] |
//! | `CutEdge` | `cxx:790-947` | [`cut_edge`] (`cut_edge.rs`) |
//! | `NewWires` | `cxx:951-954` | [`BRepAlgoLoop::new_wires`] |
//! | `WiresToFaces` | `cxx:965-992` | [`BRepAlgoLoop::wires_to_faces`] |
//! | `NewFaces` | `cxx:958-961` | [`BRepAlgoLoop::new_faces`] |
//! | `NewEdges` | `cxx:996-999` | [`BRepAlgoLoop::new_edges`] |
//! | `GetVerticesForSubstitute` / `VerticesForSubstitute` | `cxx:1003-1014` | [`BRepAlgoLoop::vertices_for_substitute_map`] / [`BRepAlgoLoop::set_vertices_for_substitute`] |
//! | `SetTolConf` / `GetTolConf` | `hxx:94-98` | [`BRepAlgoLoop::set_tol_conf`] / [`BRepAlgoLoop::tol_conf`] |
//! | `UpdateVEmap` | `cxx:1018-1151` | [`update_ve_map`] (`update_ve_map.rs`) |
//! | `Bubble`, `UpdateClosedEdge` | `cxx:80-227` | `cut_edge.rs` |
//! | `RemovePendingEdges`, `SamePnt2d`, `SelectEdge`, `PurgeNewEdges`, `StoreInMVE` | `cxx:230-578` | `mve.rs` |
//!
//! Translation notes:
//!
//! * `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>` over
//!   `TopTools_ShapeMapHasher` is `HashMap<ShapeMapKey, Vec<TopoShape>>`
//!   ([`shape_map_key`], `IsSame` identity) with `Vec` in insertion order.
//! * `NCollection_IndexedDataMap` keeps insertion order and is indexed 1-based
//!   in OCCT; [`VertexEdgeMap`] reproduces both, because `Perform` starts each
//!   wire from `MVE(1).First()` (`cxx:682`).
//! * `TopoDS_Shape::Oriented(o)` is a copy carrying `o`; the mutable
//!   `TopoDS_Shape::Orientation(o)` calls of the reference (`cxx:509`, `:1234`,
//!   `:1237`) modify a shape *view*, which the port keeps by mutating the local
//!   clone it passes on.
//! * `BRep_Tool::Pnt(V)` is [`BRepTool::vertex_point`], the convention the rest
//!   of the port uses.
//! * `BRep_Builder::UpdateVertex(V, Tol)` is [`Vertex::set_tolerance`] and
//!   `UpdateVertex(V, P, Tol)` is `set_point` followed by `set_tolerance`, which
//!   is what the registry-backed `BRep_TVertex` stores.
//! * `BRep_Builder::Range(E, U1, U2)` is
//!   [`GeometryRegistry::set_ranges_all`] (3D curve range plus every p-curve
//!   range, like `BRep_TEdge::Range`).
//! * The `OCCT_DEBUG_ALGO` block (`cxx:47-54`, `:347-357`, `:586-602`,
//!   `:764-783`) is not ported.

mod cut_edge;
mod mve;
mod update_ve_map;

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo_face_restrictor::BRepAlgoFaceRestrictor;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, ShapeMapKey};

use cut_edge::cut_edge;
use mve::{
    purge_new_edges, remove_pending_edges, same_pnt2d, select_edge, store_in_mve, VertexEdgeMap,
};
use update_ve_map::update_ve_map;

/// `class BRepAlgo_Loop` (`hxx:33-112`).
pub struct BRepAlgoLoop {
    /// `myFace` (`hxx:101`).
    my_face: Option<Face>,
    /// `myConstEdges` (`hxx:102`).
    my_const_edges: Vec<TopoShape>,
    /// `myEdges` (`hxx:103`).
    my_edges: Vec<TopoShape>,
    /// `myVerOnEdges` (`hxx:104-105`).
    my_ver_on_edges: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myNewWires` (`hxx:106`).
    my_new_wires: Vec<TopoShape>,
    /// `myNewFaces` (`hxx:107`).
    my_new_faces: Vec<TopoShape>,
    /// `myCutEdges` (`hxx:108-109`).
    my_cut_edges: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myVerticesForSubstitute` (`hxx:110`).
    my_vertices_for_substitute: HashMap<ShapeMapKey, TopoShape>,
    /// `myImageVV` (`hxx:111`).
    my_image_vv: BRepAlgoImage,
    /// `myTolConf` (`hxx:112`), initialised to `0.001` (`cxx:57-60`).
    my_tol_conf: f64,
}

impl Default for BRepAlgoLoop {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAlgoLoop {
    /// `BRepAlgo_Loop()` (`cxx:57-60`).
    pub fn new() -> Self {
        Self {
            my_face: None,
            my_const_edges: Vec::new(),
            my_edges: Vec::new(),
            my_ver_on_edges: HashMap::new(),
            my_new_wires: Vec::new(),
            my_new_faces: Vec::new(),
            my_cut_edges: HashMap::new(),
            my_vertices_for_substitute: HashMap::new(),
            my_image_vv: BRepAlgoImage::new(),
            my_tol_conf: 0.001,
        }
    }

    /// `Init(F)` (`cxx:64-74`).
    pub fn init(&mut self, f: &Face) {
        self.my_const_edges.clear();
        self.my_edges.clear();
        self.my_ver_on_edges.clear();
        self.my_new_wires.clear();
        self.my_new_faces.clear();
        self.my_cut_edges.clear();
        self.my_face = Some(f.clone());
    }

    /// `AddEdge(E, LV)` (`cxx:126-130`).
    pub fn add_edge(&mut self, e: &Edge, lv: &[TopoShape]) {
        self.my_edges.push(e.0.clone());
        self.my_ver_on_edges.insert(shape_map_key(&e.0), lv.to_vec());
    }

    /// `AddConstEdge(E)` (`cxx:134-137`).
    pub fn add_const_edge(&mut self, e: &Edge) {
        self.my_const_edges.push(e.0.clone());
    }

    /// `AddConstEdges(LE)` (`cxx:140-147`).
    pub fn add_const_edges(&mut self, le: &[TopoShape]) {
        for s in le {
            self.my_const_edges.push(s.clone());
        }
    }

    /// `SetImageVV` (`cxx:151-154`).
    pub fn set_image_vv(&mut self, the_image_vv: BRepAlgoImage) {
        self.my_image_vv = the_image_vv;
    }

    /// `SetTolConf` (`hxx:94`).
    pub fn set_tol_conf(&mut self, tol: f64) {
        self.my_tol_conf = tol;
    }

    /// `GetTolConf` (`hxx:98`).
    pub fn tol_conf(&self) -> f64 {
        self.my_tol_conf
    }

    /// `Perform()` (`cxx:581-787`).
    pub fn perform(&mut self) {
        let face = match self.my_face.clone() {
            Some(f) => f,
            None => return,
        };
        let mut ya_couture = false;

        // `cxx:604-616`: cut every edge that carries vertices.
        for e in self.my_edges.clone() {
            if let Some(verts) = self.my_ver_on_edges.get(&shape_map_key(&e)).cloned() {
                let mut lce = Vec::new();
                cut_edge(&Edge(e.clone()), &verts, &mut lce, &face, self.my_tol_conf);
                self.my_cut_edges.insert(shape_map_key(&e), lce);
            }
        }

        // `cxx:619-639`: the cut edges, one `StoreInMVE` per distinct edge.
        let mut mve = VertexEdgeMap::new();
        let mut emap: HashSet<ShapeMapKey> = HashSet::new();
        for e in self.my_edges.clone() {
            if let Some(lce) = self.my_cut_edges.get(&shape_map_key(&e)).cloned() {
                for mut ne in lce {
                    if !emap.insert(shape_map_key(&ne)) {
                        continue;
                    }
                    store_in_mve(
                        &face,
                        &mut ne,
                        &mut mve,
                        &mut ya_couture,
                        &mut self.my_vertices_for_substitute,
                        self.my_tol_conf,
                    );
                }
            }
        }

        // `cxx:642-652`: the const edges, once each.
        let mut deja_vu: HashSet<ShapeMapKey> = HashSet::new();
        for mut e in self.my_const_edges.clone() {
            if !deja_vu.insert(shape_map_key(&e)) {
                continue;
            }
            store_in_mve(
                &face,
                &mut e,
                &mut mve,
                &mut ya_couture,
                &mut self.my_vertices_for_substitute,
                self.my_tol_conf,
            );
        }

        let builder = TopoBuilder::new();
        update_ve_map(&mut mve, &mut self.my_image_vv, &self.my_const_edges);

        let mut used_edges: HashSet<ShapeMapKey> = HashSet::new();

        while mve.extent() > 0 {
            let mut nw = TopoShape::new(ShapeType::Wire);

            // `cxx:670-677`.
            remove_pending_edges(&mut mve);
            if mve.extent() == 0 {
                break;
            }

            // `cxx:680-695`: the starting edge and the vertex the wire starts at.
            let ef = match mve.find_from_index(1).and_then(|l| l.first().cloned()) {
                Some(s) => s,
                None => break,
            };
            let mut ce = ef.clone();
            let (v1, v2) = edge_vertices(&Edge(ce.clone()));
            let start = if ce.orientation() == Orientation::Forward {
                v1
            } else {
                v2
            };
            let vf = start.clone();
            let Some(cv_start) = start else {
                break;
            };
            let mut cv = Some(cv_start.clone());

            if !mve.contains(&cv_start) {
                continue;
            }
            if let Some(list) = mve.change_from_key(&cv_start) {
                if let Some(i) = list.iter().position(|x| {
                    crate::topo_tools::is_equal(x, &ce)
                }) {
                    list.remove(i);
                }
            }
            let mut end = false;

            while !end {
                // `cxx:715-723`: the other end of the edge.
                let (a, b) = edge_vertices(&Edge(ce.clone()));
                let same_as_v1 = match (&cv, &a) {
                    (Some(c), Some(a)) => is_same(&c.0, &a.0),
                    (None, None) => true,
                    _ => false,
                };
                cv = if !same_as_v1 { a } else { b };

                builder.add(&mut nw, &ce);
                used_edges.insert(shape_map_key(&ce));

                let Some(cur) = cv.clone() else {
                    end = true;
                    continue;
                };
                let empty = mve.find_from_key(&cur).map(|l| l.is_empty()).unwrap_or(true);
                if empty {
                    end = true;
                } else {
                    let mut ne: Option<TopoShape> = None;
                    let found = match mve.change_from_key(&cur) {
                        Some(list) => select_edge(&face, &ce, &cur, &mut ne, list),
                        None => false,
                    };
                    end = !found;
                    if !end {
                        if let Some(ne) = ne {
                            ce = ne;
                        }
                        if mve.find_from_key(&cur).map(|l| l.is_empty()).unwrap_or(false) {
                            mve.remove_key(&cur);
                        }
                    }
                }
            }

            // `cxx:746-784`.
            let mut closed = false;
            if let (Some(a), Some(b)) = (&vf, &cv) {
                if is_same(&a.0, &b.0) {
                    if same_pnt2d(a, &ef, &ce, &face) {
                        closed = true;
                    } else if BRepTool::vertex_tolerance(a) < self.my_tol_conf {
                        a.set_tolerance(self.my_tol_conf);
                        closed = same_pnt2d(a, &ef, &ce, &face);
                    }
                }
            }
            if closed {
                nw.set_closed(true);
                self.my_new_wires.push(nw);
            }
        }

        purge_new_edges(&mut self.my_cut_edges, &used_edges);
    }

    /// `NewWires()` (`cxx:951-954`).
    pub fn new_wires(&self) -> &[TopoShape] {
        &self.my_new_wires
    }

    /// `WiresToFaces()` (`cxx:965-992`).
    pub fn wires_to_faces(&mut self) {
        if self.my_new_wires.is_empty() {
            return;
        }
        let Some(face) = self.my_face.clone() else {
            return;
        };
        let mut fr = BRepAlgoFaceRestrictor::new();
        fr.init(Some(Face(face.clone().oriented(Orientation::Forward))), false, false);
        for w in &self.my_new_wires {
            fr.add(w);
        }
        fr.perform();
        if fr.is_done() {
            let ori_f = face.orientation();
            for f in fr.faces().to_vec() {
                self.my_new_faces.push(f.oriented(ori_f));
            }
        }
    }

    /// `NewFaces()` (`cxx:958-961`).
    pub fn new_faces(&self) -> &[TopoShape] {
        &self.my_new_faces
    }

    /// `NewEdges(E)` (`cxx:996-999`): `myCutEdges(E)`, which raises when `E`
    /// was not cut (`NCollection_DataMap::operator()`).
    pub fn new_edges(&self, e: &Edge) -> Option<&Vec<TopoShape>> {
        self.my_cut_edges.get(&shape_map_key(&e.0))
    }

    /// `GetVerticesForSubstitute` (`cxx:1003-1007`).
    pub fn vertices_for_substitute(&self) -> &HashMap<ShapeMapKey, TopoShape> {
        &self.my_vertices_for_substitute
    }

    /// `VerticesForSubstitute` (`cxx:1012-1015`).
    pub fn set_vertices_for_substitute(&mut self, map: HashMap<ShapeMapKey, TopoShape>) {
        self.my_vertices_for_substitute = map;
    }
}

/// `BRep_Builder::UpdateVertex(V, theTol)` on a vertex that may be null.
#[allow(dead_code)]
pub(crate) fn update_vertex_tolerance(v: Option<&Vertex>, tol: f64) {
    if let Some(v) = v {
        v.set_tolerance(tol);
    }
}

/// `BRep_Tool::Pnt(V)`.
#[allow(dead_code)]
pub(crate) fn vertex_pnt(v: &Vertex) -> occt_core::gp::GpPnt {
    BRepTool::vertex_point(v)
}

/// `TopoDS::Wire(shape)`; `None` when the shape is not a wire.
#[allow(dead_code)]
pub(crate) fn as_wire(s: &TopoShape) -> Option<Wire> {
    Wire::wrap(s.clone())
}
