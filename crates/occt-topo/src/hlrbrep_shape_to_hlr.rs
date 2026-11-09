//! `HLRBRep_ShapeToHLR` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_ShapeToHLR.hxx:32-86`,
//! `.cxx:40-313`): turns one [`HlrTopoBRepOutLiner`] into a filled
//! [`HlrBRepData`].
//!
//! `Load` (`cxx:40-215`) runs `S->Fill(P, MST, nbIso)` (`cxx:46`), merges the
//! outlined faces/edges into two index maps (`cxx:58-91`), counts the vertices
//! through `VerticesToEdges` (`cxx:93-101`), creates the data structure
//! (`cxx:116`) and loads every edge (`cxx:127-199`) before `ExploreShape`
//! (`cxx:201`). `ExploreFace` (`cxx:221-283`) fills one `HLRBRep_FaceData`;
//! `ExploreShape` (`cxx:289-313`) walks the shells of the *original* shape and
//! decides the `closed` flag of each one.
//!
//! Deviations:
//! * `NCollection_IndexedMap<TopoDS_Shape>` is a `Vec<TopoShape>` plus a
//!   `HashMap<ShapeMapKey, usize>` index, reproducing `Add` / `FindIndex` /
//!   `Extent` (`FM`, `EM`, `VerticesToEdges`, `EdgesToFaces`).
//! * `TopExp::MapShapesAndAncestors(EM(i), VERTEX, EDGE, VerticesToEdges)` is
//!   the loop below: it keeps the first-appearance order OCCT's indexed map
//!   uses, which is what `FindIndex(V)` returns.
//! * `BRep_Tool::Continuity(E, F1, F2)` (`BRep_Tool.cxx:1180-1188`) is
//!   [`stored_continuity`] (`brep_lib_continuity.rs:67-73`); `GeomAbs_G1` /
//!   `GeomAbs_G2` are the `1` / `3` values of that store.
//! * `HLRBRep_EdgeData::Set` / `HLRBRep_FaceData::Set` / `SetWire` / `SetWEdge`
//!   already exist in the port and are called with the same arguments.
//! * `EDataArray().ChangeValue(1)` + pointer increment (`cxx:119-124`) becomes
//!   `e_data_mut(i)` inside the loop.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::precision;

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_continuity::stored_continuity;
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::hlr_algo::{
    add_min_max, copy_min_max, decode_min_max, encode_min_max, MinMaxIndices,
};
use crate::hlrbrep_data::HlrBRepData;
use crate::hlrtopobrep::{HlrTopoBRepOutLiner, ShapeToTopAdaptorToolMap};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey,
};

/// `RealFirst()` (`Standard_Real.hxx`).
const REAL_FIRST: f64 = -f64::MAX;
/// `RealLast()` (`Standard_Real.hxx`).
const REAL_LAST: f64 = f64::MAX;

/// An `NCollection_IndexedMap<TopoDS_Shape>`: insertion order plus `FindIndex`.
#[derive(Debug, Default, Clone)]
struct IndexedMap {
    /// Slots `1..=len`, the OCCT indexing.
    items: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
}

impl IndexedMap {
    fn new() -> Self {
        Self::default()
    }

    /// `Add(shape)`: appends when absent and answers `true`, as OCCT's
    /// `NCollection_Map::Add` does for the `ShapeMap` of `cxx:294`.
    fn add(&mut self, shape: &TopoShape) -> bool {
        let key = shape_map_key(shape);
        if self.index.contains_key(&key) {
            return false;
        }
        self.items.push(shape.clone());
        self.index.insert(key, self.items.len());
        true
    }

    /// `FindIndex(shape)`: the 1-based index, `0` when absent.
    fn find_index(&self, shape: &TopoShape) -> i32 {
        self.index
            .get(&shape_map_key(shape))
            .map(|i| *i as i32)
            .unwrap_or(0)
    }

    /// `Extent()`.
    fn extent(&self) -> i32 {
        self.items.len() as i32
    }

    /// `operator()(i)`, 1-based.
    fn get(&self, i: i32) -> Option<&TopoShape> {
        if i < 1 {
            return None;
        }
        self.items.get(i as usize - 1)
    }
}

/// The `TopTools_ShapeMapHasher` key of the ancestors of one sub-shape: an
/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// in first-appearance order.
#[derive(Debug, Default, Clone)]
struct IndexedAncestorMap {
    /// Slot 0 unused; slot `k` is the ancestor list of `keys[k-1]`.
    keys: Vec<TopoShape>,
    ancestors: Vec<Vec<TopoShape>>,
    index: HashMap<ShapeMapKey, usize>,
}

impl IndexedAncestorMap {
    fn new() -> Self {
        Self::default()
    }

    /// `TopExp::MapShapesAndAncestors(ancestor, TS, TA, theMap)`: appends
    /// `ancestor` to every `TS` sub-shape of `ancestor`.
    fn map_shapes_and_ancestors(&mut self, ancestor: &TopoShape, ts: ShapeType) {
        for sub in explorer_shapes(ancestor, ts) {
            let key = shape_map_key(&sub);
            let slot = match self.index.get(&key) {
                Some(slot) => *slot,
                None => {
                    self.keys.push(sub.clone());
                    self.ancestors.push(Vec::new());
                    let slot = self.ancestors.len();
                    self.index.insert(key, slot);
                    slot
                }
            };
            self.ancestors[slot - 1].push(ancestor.clone());
        }
    }

    /// `FindIndex(shape)`: the 1-based index, `0` when absent.
    fn find_index(&self, shape: &TopoShape) -> i32 {
        self.index
            .get(&shape_map_key(shape))
            .map(|i| *i as i32)
            .unwrap_or(0)
    }

    /// `Extent()`.
    fn extent(&self) -> i32 {
        self.keys.len() as i32
    }

    /// `operator()(i)`, 1-based: the ancestor list.
    fn get(&self, i: i32) -> Option<&[TopoShape]> {
        if i < 1 {
            return None;
        }
        self.ancestors.get(i as usize - 1).map(|v| v.as_slice())
    }
}

/// `HLRBRep_ShapeToHLR::Load(S, P, MST, nbIso)` (`cxx:40-215`).
///
/// `S` is the outlined shape; OCCT reaches it through a `handle` and writes
/// into it, so the caller's `Arc` must be unique here (the port's callers pass
/// a freshly created outliner).
pub fn load(
    s: &mut HlrTopoBRepOutLiner,
    p: &crate::hlr_algo::HlrAlgoProjector,
    mst: &mut ShapeToTopAdaptorToolMap,
    nb_iso: i32,
) -> Result<HlrBRepData, String> {
    // `S->Fill(P, MST, nbIso);` (`cxx:46`).
    s.fill(p, mst, nb_iso)?;

    // `HLRTopoBRep_Data& TopDS = S->DataStructure();` (`cxx:48`).
    let top_ds = s.data_structure().clone();

    // `cxx:58-72`: faces in a shell.
    let mut fm = IndexedMap::new();
    let out_lined = s
        .out_lined_shape()
        .cloned()
        .ok_or_else(|| "HLRBRep_ShapeToHLR::Load: OutLinedShape is null".to_string())?;
    for shell_shape in explorer_shapes(&out_lined, ShapeType::Shell) {
        for face_shape in explorer_shapes(&shell_shape, ShapeType::Face) {
            fm.add(&face_shape);
        }
    }
    // `cxx:74-81`: faces not in a shell.
    for face_shape in explorer_shapes_avoiding(&out_lined, ShapeType::Face, ShapeType::Shell) {
        fm.add(&face_shape);
    }

    // `TopExp::MapShapes(S->OutLinedShape(), TopAbs_EDGE, EM);` (`cxx:84`).
    let mut em = IndexedMap::new();
    for edge_shape in explorer_shapes(&out_lined, ShapeType::Edge) {
        em.add(&edge_shape);
    }
    let nb_edge = em.extent();

    // `cxx:88-91`: `VerticesToEdges`.
    let mut vertices_to_edges = IndexedAncestorMap::new();
    for i in 1..=nb_edge {
        if let Some(e) = em.get(i) {
            vertices_to_edges.map_shapes_and_ancestors(e, ShapeType::Vertex);
        }
    }
    let nb_vert = vertices_to_edges.extent();
    let nb_face = fm.extent();

    // `occ::handle<HLRBRep_Data> DS = new HLRBRep_Data(nbVert, nbEdge, nbFace);`
    // (`cxx:116`).
    let mut ds = HlrBRepData::new(nb_vert, nb_edge, nb_face);

    // `cxx:125-127`: `EdgesToFaces`.
    let mut edges_to_faces = IndexedAncestorMap::new();
    for i in 1..=nb_face {
        if let Some(f) = fm.get(i) {
            edges_to_faces.map_shapes_and_ancestors(f, ShapeType::Edge);
        }
    }

    // `cxx:130-199`: load the edges.
    for i in 1..=nb_edge {
        let Some(edg_shape) = em.get(i) else {
            continue;
        };
        let Some(edg) = Edge::wrap(edg_shape.clone()) else {
            continue;
        };
        // `TopExp::Vertices(Edg, VF, VL);` (`cxx:132`).
        let (vf, vl) = crate::topo_tools::edge_vertices(&edg);
        // `BRep_Tool::Range(Edg, pf, pl);` (`cxx:133`).
        let (mut pf, mut pl) = BRepTool::edge_parameters(&edg);

        // `cxx:134-150`: regularity from the two faces sharing the edge.
        let mut reg1 = false;
        let mut regn = false;
        let inde = edges_to_faces.find_index(&edg);
        if inde > 0 {
            if let Some(faces) = edges_to_faces.get(inde) {
                if faces.len() == 2 {
                    let f1 = Face::wrap(faces[0].clone());
                    let f2 = Face::wrap(faces[1].clone());
                    if let (Some(f1), Some(f2)) = (f1, f2) {
                        // `GeomAbs_Shape rg = BRep_Tool::Continuity(Edg, F1, F2);`
                        let rg = stored_continuity(&edg, &f1, &f2);
                        reg1 = rg >= 1; // `GeomAbs_G1`
                        regn = rg >= 3; // `GeomAbs_G2`
                    }
                }
            }
        }

        // `cxx:152-186`: the two vertices.
        // `cxx:152-170`: the first vertex.
        let (i1, o1, c1, tf) = match &vf {
            None => {
                pf = REAL_FIRST;
                (0, false, false, precision::epsilon(pf) as f32)
            }
            Some(v) => (
                vertices_to_edges.find_index(v.shape()),
                top_ds.is_out_v(v),
                top_ds.is_int_v(v),
                BRepTool::vertex_tolerance(v) as f32,
            ),
        };
        // `cxx:172-186`: the last vertex.
        let (i2, o2, c2, tl) = match &vl {
            None => {
                pl = REAL_LAST;
                (0, false, false, precision::epsilon(pl) as f32)
            }
            Some(v) => (
                vertices_to_edges.find_index(v.shape()),
                top_ds.is_out_v(v),
                top_ds.is_int_v(v),
                BRepTool::vertex_tolerance(v) as f32,
            ),
        };

        // `ed->Set(reg1, regn, Edg, i1, i2, o1, o2, c1, c2, pf, tf, pl, tl);`
        // (`cxx:188`).
        ds.e_data_mut(i).set(
            reg1, regn, &edg, i1, i2, o1, o2, c1, c2, pf, tf, pl, tl,
        );
        // `DS->EdgeMap().Add(Edg);` (`cxx:189`).
        ds.edge_map_push(edg.shape());
    }

    // `ExploreShape(S, DS, FM, EM); return DS;` (`cxx:201-202`).
    explore_shape(s, &mut ds, &fm, &em)?;
    Ok(ds)
}

/// `HLRBRep_ShapeToHLR::ExploreFace(S, DS, FM, EM, i, F, closed)`
/// (`cxx:221-283`).
fn explore_face(
    s: &HlrTopoBRepOutLiner,
    ds: &mut HlrBRepData,
    fm: &IndexedMap,
    em: &IndexedMap,
    i: &mut i32,
    f: &Face,
    closed: bool,
) {
    *i += 1;
    let top_ds = s.data_structure();

    // `TopAbs_Orientation orient = FM(i).Orientation();`
    // `TopoDS_Face theFace = TopoDS::Face(FM(i)); theFace.Orientation(FORWARD);`
    // (`cxx:231-234`).
    let Some(fm_face) = fm.get(*i) else {
        return;
    };
    let orient = fm_face.orientation();
    let mut the_face_shape = fm_face.clone();
    the_face_shape.set_orientation(Orientation::Forward);
    let Some(the_face) = Face::wrap(the_face_shape) else {
        return;
    };

    // `int nw = 0; for (Ex1.Init(theFace, WIRE); ...) nw++;` (`cxx:238-245`).
    let wires = explorer_shapes(the_face.shape(), ShapeType::Wire);
    // `fd.Set(theFace, orient, closed, nw);` (`cxx:247`).
    ds.f_data_mut(*i).set(&the_face, orient, closed, wires.len());

    // `cxx:250-279`: one wire at a time.
    for (wi, wire_shape) in wires.iter().enumerate() {
        let nw = wi + 1;
        let wire_edges = explorer_shapes(wire_shape, ShapeType::Edge);
        let mut ne = 0usize;
        for e_shape in &wire_edges {
            if let Some(e) = Edge::wrap(e_shape.clone()) {
                if !BRepTool::is_degenerated(&e) {
                    ne += 1;
                }
            }
        }
        // `fd.SetWire(nw, ne);` (`cxx:262`).
        ds.f_data_mut(*i).set_wire(nw, ne);

        ne = 0;
        for e_shape in &wire_edges {
            let Some(e) = Edge::wrap(e_shape.clone()) else {
                continue;
            };
            if BRepTool::is_degenerated(&e) {
                continue; // `cxx:269-272`
            }
            ne += 1;
            let ie = em.find_index(&e);
            let an_orient_e = e.shape().orientation();
            let inte = top_ds.is_int_l_face_edge(f, &e);
            let iso = top_ds.is_iso_l_face_edge(f, &e);
            let out = top_ds.is_out_l_face_edge(f, &e);
            let dbl = is_really_closed(&e, &the_face);
            ds.f_data_mut(*i)
                .set_wedge(nw, ne, ie, an_orient_e, out, inte, dbl, iso);
        }
    }
    // `DS->FaceMap().Add(theFace);` (`cxx:282`).
    ds.face_map_push(the_face.shape());
}

/// `HLRBRep_ShapeToHLR::ExploreShape(S, DS, FM, EM)` (`cxx:289-313`).
fn explore_shape(
    s: &HlrTopoBRepOutLiner,
    ds: &mut HlrBRepData,
    fm: &IndexedMap,
    em: &IndexedMap,
) -> Result<(), String> {
    let original = s
        .original_shape()
        .cloned()
        .ok_or_else(|| "HLRBRep_ShapeToHLR::ExploreShape: OriginalShape is null".to_string())?;

    // `NCollection_Map<TopoDS_Shape> ShapeMap; int i = 0;` (`cxx:292-293`).
    let mut shape_map = IndexedMap::new();
    let mut i = 0;

    // `cxx:295-308`: faces in a shell.
    for shell_shape in explorer_shapes(&original, ShapeType::Shell) {
        // `bool closed = exshell.Current().Closed();` (`cxx:298`).
        let mut closed = shell_shape.is_closed();

        // `cxx:300-326`: an open shell is closed when every edge is used once
        // forward and once reversed.
        if !closed {
            let nb_edge = em.extent();
            let mut flag = vec![0i32; nb_edge as usize + 1];
            for e_shape in explorer_shapes(&shell_shape, ShapeType::Edge) {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let ie = em.find_index(&e);
                if ie < 0 || ie > nb_edge {
                    continue;
                }
                let orient = e.shape().orientation();
                if !BRepTool::is_degenerated(&e) {
                    if orient == Orientation::Forward {
                        flag[ie as usize] += 1;
                    } else if orient == Orientation::Reversed {
                        flag[ie as usize] -= 1;
                    }
                }
            }
            closed = true;
            for ie in 1..=nb_edge {
                if flag[ie as usize] != 0 {
                    closed = false;
                    break;
                }
            }
        }

        // `cxx:310-318`.
        for face_shape in explorer_shapes(&shell_shape, ShapeType::Face) {
            if !shape_map.add(&face_shape) {
                continue; // `ShapeMap.Add(...)`
            }
            let Some(face) = Face::wrap(face_shape) else {
                continue;
            };
            explore_face(s, ds, fm, em, &mut i, &face, closed);
        }
    }

    // `cxx:320-326`: faces not in a shell, always `closed = false`.
    let free_faces = explorer_shapes_avoiding(&original, ShapeType::Face, ShapeType::Shell);
    for face_shape in free_faces {
        if !shape_map.add(&face_shape) {
            continue;
        }
        let Some(face) = Face::wrap(face_shape) else {
            continue;
        };
        explore_face(s, ds, fm, em, &mut i, &face, false);
    }
    Ok(())
}

/// The `MinMaxIndices` accumulator type re-exported for
/// `HLRBRep_InternalAlgo::Update` (`HLRBRep_InternalAlgo.cxx:150-184`).
pub type ShapeMinMax = MinMaxIndices;

/// `HLRBRep_InternalAlgo::Update`'s per-shape bounds fold
/// (`HLRBRep_InternalAlgo.cxx:150-184`): `DecodeMinMax` the edge bounds of
/// `e1..=e2`, the wire bounds of `f1..=f2`, then `EncodeMinMax` into the block.
pub fn update_min_max_of_shape(
    ds: &HlrBRepData,
    e1: i32,
    e2: i32,
    f1: i32,
    f2: i32,
    out: &mut MinMaxIndices,
) {
    let mut shap_min = MinMaxIndices::default();
    let mut shap_max = MinMaxIndices::default();
    let mut first_time = true;
    let mut the_min = MinMaxIndices::default();
    let mut the_max = MinMaxIndices::default();

    for e in e1..=e2 {
        decode_min_max(ds.e_data(e).min_max(), &mut the_min, &mut the_max);
        if first_time {
            first_time = false;
            let (mn, mx) = copy_min_max(&the_min, &the_max);
            shap_min = mn;
            shap_max = mx;
        } else {
            add_min_max(&the_min, &the_max, &mut shap_min, &mut shap_max);
        }
    }

    for f in f1..=f2 {
        if let Some(wires) = ds.f_data(f).wires() {
            decode_min_max(wires.min_max(), &mut the_min, &mut the_max);
            add_min_max(&the_min, &the_max, &mut shap_min, &mut shap_max);
        }
    }
    encode_min_max(&shap_min, &shap_max, out);
}

/// Convenience wrapper for `const handle(HLRTopoBRep_OutLiner)&` callers: clones
/// the outliner out of an `Arc`, fills it and answers the data structure. The
/// caller must rebind the filled outliner, as OCCT's handle mutation would.
pub fn load_shared(
    s: &Arc<HlrTopoBRepOutLiner>,
    p: &crate::hlr_algo::HlrAlgoProjector,
    mst: &mut ShapeToTopAdaptorToolMap,
    nb_iso: i32,
) -> Result<(HlrBRepData, HlrTopoBRepOutLiner), String> {
    let mut owned = (**s).clone();
    let ds = load(&mut owned, p, mst, nb_iso)?;
    Ok((ds, owned))
}
