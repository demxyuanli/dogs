//! `HLRTopoBRep_Data` (`ModelingAlgorithms/TKHLR/HLRTopoBRep/HLRTopoBRep_Data.hxx:39-144`,
//! `HLRTopoBRep_Data.cxx:28-303`, `HLRTopoBRep_Data.lxx:26-106`): the result
//! store of the outline (`HLRTopoBRep_OutLiner`) and iso-line
//! (`HLRTopoBRep_FaceIsoLiner`) passes.
//!
//! The OCCT maps are `NCollection_DataMap` / `NCollection_Map` keyed by
//! `TopoDS_Shape` with `TopTools_ShapeMapHasher`, i.e. keyed by
//! `(TShape, Location)` - exactly [`ShapeMapKey`]. Iteration order of
//! `NCollection_DataMap` follows hash buckets and is not part of the contract;
//! the port keeps insertion order for `myEdgesVertices` (the only map that is
//! iterated through `InitEdge` / `NextEdge`) so the traversal is deterministic,
//! and notes it here instead of reproducing bucket order.

use std::collections::{HashMap, HashSet};

use occt_topo::shape::{Edge, Face, TopoShape, Vertex};
use occt_topo::topo_tools::{is_same_occurrence, shape_map_key, ShapeMapKey};

use super::face_data::HlrTopoBRepFaceData;
use super::v_data::HlrTopoBRepVData;

/// One `myEdgesVertices` entry: the edge itself (OCCT recovers it from the map
/// key, `lxx:113-116`) plus its vertex list.
#[derive(Debug, Clone, Default)]
struct EdgeVertices {
    edge: Option<TopoShape>,
    vertices: Vec<HlrTopoBRepVData>,
}

/// `HLRTopoBRep_Data`.
#[derive(Debug, Clone, Default)]
pub struct HlrTopoBRepData {
    /// `myOldS` (`hxx:133`): new shape -> original shape.
    my_old_s: HashMap<ShapeMapKey, TopoShape>,
    /// `mySplE` (`hxx:134`): split edge -> the pieces.
    my_spl_e: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myData` (`hxx:135`): face -> its three outline lists.
    my_data: HashMap<ShapeMapKey, HlrTopoBRepFaceData>,
    /// `myOutV` (`hxx:136`): outline vertices on restriction.
    my_out_v: HashSet<ShapeMapKey>,
    /// `myIntV` (`hxx:137`): internal outline vertices.
    my_int_v: HashSet<ShapeMapKey>,
    /// `myEdgesVertices` (`hxx:138-139`).
    my_edges_vertices: HashMap<ShapeMapKey, EdgeVertices>,
    /// Insertion order of `myEdgesVertices` (see the module note).
    edge_order: Vec<ShapeMapKey>,
    /// `myEIterator` (`hxx:140-141`) as a key snapshot + cursor.
    edge_iter_keys: Vec<ShapeMapKey>,
    edge_iter_pos: usize,
    /// `myVIterator` (`hxx:142`) as an index into the current edge's list.
    current_edge: Option<ShapeMapKey>,
    vertex_iter_pos: usize,
}

impl HlrTopoBRepData {
    /// `HLRTopoBRep_Data()` (`cxx:28-31`): `Clear()`.
    pub fn new() -> Self {
        let mut self_ = Self::default();
        self_.clear();
        self_
    }

    /// `Clear()` (`cxx:35-45`).
    pub fn clear(&mut self) {
        self.my_old_s.clear();
        self.my_spl_e.clear();
        self.my_data.clear();
        self.my_out_v.clear();
        self.my_int_v.clear();
        self.my_edges_vertices.clear();
        self.edge_order.clear();
        self.edge_iter_keys.clear();
        self.edge_iter_pos = 0;
        self.current_edge = None;
        self.vertex_iter_pos = 0;
    }

    /// `Clean()` (`cxx:47`): empty body in OCCT.
    pub fn clean(&mut self) {}

    /// `EdgeHasSplE(E)` (`cxx:51-60`).
    pub fn edge_has_spl_e(&self, e: &Edge) -> bool {
        match self.my_spl_e.get(&shape_map_key(e.shape())) {
            Some(list) => !list.is_empty(),
            None => false,
        }
    }

    /// `FaceHasIntL(F)` (`cxx:62-71`).
    pub fn face_has_int_l(&self, f: &Face) -> bool {
        match self.my_data.get(&shape_map_key(f.shape())) {
            Some(data) => !data.face_int_l().is_empty(),
            None => false,
        }
    }

    /// `FaceHasOutL(F)` (`cxx:73-82`).
    pub fn face_has_out_l(&self, f: &Face) -> bool {
        match self.my_data.get(&shape_map_key(f.shape())) {
            Some(data) => !data.face_out_l().is_empty(),
            None => false,
        }
    }

    /// `FaceHasIsoL(F)` (`cxx:84-93`).
    pub fn face_has_iso_l(&self, f: &Face) -> bool {
        match self.my_data.get(&shape_map_key(f.shape())) {
            Some(data) => !data.face_iso_l().is_empty(),
            None => false,
        }
    }

    /// `IsSplEEdgeEdge(E1, E2)` (`cxx:95-114`).
    pub fn is_spl_e_edge_edge(&self, e1: &Edge, e2: &Edge) -> bool {
        let mut found = false;
        if self.edge_has_spl_e(e1) {
            if let Some(list) = self.my_spl_e.get(&shape_map_key(e1.shape())) {
                for s in list {
                    found = is_same_occurrence(s, e2.shape());
                    if found {
                        break;
                    }
                }
            }
        } else {
            found = is_same_occurrence(e1.shape(), e2.shape());
        }
        found
    }

    /// `IsIntLFaceEdge(F, E)` (`cxx:116-131`).
    pub fn is_int_l_face_edge(&self, f: &Face, e: &Edge) -> bool {
        let mut found = false;
        if self.face_has_int_l(f) {
            if let Some(data) = self.my_data.get(&shape_map_key(f.shape())) {
                for s in data.face_int_l() {
                    found = self.is_shape_edge_split(s, e);
                    if found {
                        break;
                    }
                }
            }
        }
        found
    }

    /// `IsOutLFaceEdge(F, E)` (`cxx:133-148`).
    pub fn is_out_l_face_edge(&self, f: &Face, e: &Edge) -> bool {
        let mut found = false;
        if self.face_has_out_l(f) {
            if let Some(data) = self.my_data.get(&shape_map_key(f.shape())) {
                for s in data.face_out_l() {
                    found = self.is_shape_edge_split(s, e);
                    if found {
                        break;
                    }
                }
            }
        }
        found
    }

    /// `IsIsoLFaceEdge(F, E)` (`cxx:150-165`).
    pub fn is_iso_l_face_edge(&self, f: &Face, e: &Edge) -> bool {
        let mut found = false;
        if self.face_has_iso_l(f) {
            if let Some(data) = self.my_data.get(&shape_map_key(f.shape())) {
                for s in data.face_iso_l() {
                    found = self.is_shape_edge_split(s, e);
                    if found {
                        break;
                    }
                }
            }
        }
        found
    }

    /// `IsSplEEdgeEdge(TopoDS::Edge(list item), E)`: the list items are stored
    /// as `TopoDS_Shape` and down-cast with `TopoDS::Edge`.
    fn is_shape_edge_split(&self, shape: &TopoShape, e: &Edge) -> bool {
        match Edge::wrap(shape.clone()) {
            Some(edge) => self.is_spl_e_edge_edge(&edge, e),
            None => false,
        }
    }

    /// `EdgeSplE(E)` (`lxx:24-27`): the pieces of a split edge.
    ///
    /// OCCT writes `return mySplE(E);`, i.e. an unguarded `NCollection_DataMap`
    /// lookup that throws for an unbound edge; every caller tests
    /// [`Self::edge_has_spl_e`] first (`HLRTopoBRep_OutLiner.cxx:136`, `:228`,
    /// `:269`). The port answers an empty slice for the unbound case.
    pub fn edge_spl_e(&self, e: &Edge) -> &[TopoShape] {
        match self.my_spl_e.get(&shape_map_key(e.shape())) {
            Some(list) => list,
            None => &[],
        }
    }

    /// `FaceIntL(F)` (`lxx:31-34`): the internal outlines of a face. Same
    /// unguarded-lookup note as [`Self::edge_spl_e`]; callers test
    /// [`Self::face_has_int_l`] first (`HLRTopoBRep_OutLiner.cxx:157`).
    pub fn face_int_l(&self, f: &Face) -> &[TopoShape] {
        match self.my_data.get(&shape_map_key(f.shape())) {
            Some(data) => data.face_int_l(),
            None => &[],
        }
    }

    /// `FaceOutL(F)` (`lxx:38-41`): the outlines on the face's restriction.
    pub fn face_out_l(&self, f: &Face) -> &[TopoShape] {
        match self.my_data.get(&shape_map_key(f.shape())) {
            Some(data) => data.face_out_l(),
            None => &[],
        }
    }

    /// `FaceIsoL(F)` (`lxx:45-48`): the iso lines of a face. Callers test
    /// [`Self::face_has_iso_l`] first (`HLRTopoBRep_OutLiner.cxx:260`).
    pub fn face_iso_l(&self, f: &Face) -> &[TopoShape] {
        match self.my_data.get(&shape_map_key(f.shape())) {
            Some(data) => data.face_iso_l(),
            None => &[],
        }
    }

    /// `NewSOldS(NewS)` (`cxx:167-179`).
    pub fn new_s_old_s(&self, new_s: &TopoShape) -> TopoShape {
        match self.my_old_s.get(&shape_map_key(new_s)) {
            Some(old) => old.clone(),
            None => new_s.clone(),
        }
    }

    /// `AddOldS(NewS, OldS)` (`cxx:181-189`).
    pub fn add_old_s(&mut self, new_s: &TopoShape, old_s: &TopoShape) {
        self.my_old_s
            .entry(shape_map_key(new_s))
            .or_insert_with(|| old_s.clone());
    }

    /// `AddSplE(E)` (`cxx:191-201`).
    pub fn add_spl_e(&mut self, e: &Edge) -> &mut Vec<TopoShape> {
        self.my_spl_e
            .entry(shape_map_key(e.shape()))
            .or_default()
    }

    /// `AddIntL(F)` (`cxx:203-213`).
    pub fn add_int_l(&mut self, f: &Face) -> &mut Vec<TopoShape> {
        self.my_data
            .entry(shape_map_key(f.shape()))
            .or_default()
            .add_int_l()
    }

    /// `AddOutL(F)` (`cxx:215-225`).
    pub fn add_out_l(&mut self, f: &Face) -> &mut Vec<TopoShape> {
        self.my_data
            .entry(shape_map_key(f.shape()))
            .or_default()
            .add_out_l()
    }

    /// `AddIsoL(F)` (`cxx:227-237`).
    pub fn add_iso_l(&mut self, f: &Face) -> &mut Vec<TopoShape> {
        self.my_data
            .entry(shape_map_key(f.shape()))
            .or_default()
            .add_iso_l()
    }

    /// `IsOutV(V)` (`lxx:54-57`).
    pub fn is_out_v(&self, v: &Vertex) -> bool {
        self.my_out_v.contains(&shape_map_key(v.shape()))
    }

    /// `IsIntV(V)` (`lxx:61-64`).
    pub fn is_int_v(&self, v: &Vertex) -> bool {
        self.my_int_v.contains(&shape_map_key(v.shape()))
    }

    /// `AddOutV(V)` (`lxx:68-71`).
    pub fn add_out_v(&mut self, v: &Vertex) {
        self.my_out_v.insert(shape_map_key(v.shape()));
    }

    /// `AddIntV(V)` (`lxx:75-78`).
    pub fn add_int_v(&mut self, v: &Vertex) {
        self.my_int_v.insert(shape_map_key(v.shape()));
    }

    /// `InitEdge()` (`cxx:239-249`).
    pub fn init_edge(&mut self) {
        self.edge_iter_keys = self.edge_order.clone();
        self.edge_iter_pos = 0;
        self.skip_empty_edges();
    }

    /// `NextEdge()` (`cxx:251-261`).
    pub fn next_edge(&mut self) {
        self.edge_iter_pos += 1;
        self.skip_empty_edges();
    }

    fn current_edge_key(&self) -> Option<ShapeMapKey> {
        self.edge_iter_keys.get(self.edge_iter_pos).cloned()
    }

    fn skip_empty_edges(&mut self) {
        while let Some(key) = self.current_edge_key() {
            let empty = self
                .my_edges_vertices
                .get(&key)
                .map(|e| e.vertices.is_empty())
                .unwrap_or(true);
            if !empty {
                break;
            }
            self.edge_iter_pos += 1;
        }
    }

    /// `MoreEdge()` (`lxx:82-85`).
    pub fn more_edge(&self) -> bool {
        self.edge_iter_pos < self.edge_iter_keys.len()
    }

    /// `Edge()` (`lxx:89-92`): `TopoDS::Edge(myEIterator.Key())`.
    pub fn edge(&self) -> Option<&TopoShape> {
        self.current_edge_key()
            .and_then(|key| self.my_edges_vertices.get(&key))
            .and_then(|e| e.edge.as_ref())
    }

    /// `InitVertex(E)` (`cxx:263-275`).
    pub fn init_vertex(&mut self, e: &Edge) {
        let key = shape_map_key(e.shape());
        if !self.my_edges_vertices.contains_key(&key) {
            self.edge_order.push(key.clone());
            self.my_edges_vertices.insert(
                key.clone(),
                EdgeVertices {
                    edge: Some(e.shape().clone()),
                    vertices: Vec::new(),
                },
            );
        }
        self.current_edge = Some(key);
        self.vertex_iter_pos = 0;
    }

    fn current_vertices(&self) -> Option<&Vec<HlrTopoBRepVData>> {
        self.current_edge
            .as_ref()
            .and_then(|key| self.my_edges_vertices.get(key))
            .map(|e| &e.vertices)
    }

    /// `MoreVertex()` (`lxx:96-99`).
    pub fn more_vertex(&self) -> bool {
        match self.current_vertices() {
            Some(list) => self.vertex_iter_pos < list.len(),
            None => false,
        }
    }

    /// `NextVertex()` (`lxx:103-106`).
    pub fn next_vertex(&mut self) {
        self.vertex_iter_pos += 1;
    }

    /// `Vertex()` (`cxx:277-282`): `TopoDS::Vertex(myVIterator.Value().Vertex())`.
    pub fn vertex(&self) -> Option<&TopoShape> {
        let list = self.current_vertices()?;
        list.get(self.vertex_iter_pos)?.vertex()
    }

    /// `Parameter()` (`cxx:284-289`).
    pub fn parameter(&self) -> f64 {
        match self.current_vertices().and_then(|l| l.get(self.vertex_iter_pos)) {
            Some(vd) => vd.parameter(),
            None => 0.0,
        }
    }

    /// `InsertBefore(V, P)` (`cxx:291-297`): the OCCT list iterator keeps
    /// pointing at the same node, which is now one slot further on.
    pub fn insert_before(&mut self, v: &TopoShape, p: f64) {
        let Some(key) = self.current_edge.clone() else {
            return;
        };
        if let Some(entry) = self.my_edges_vertices.get_mut(&key) {
            let pos = self.vertex_iter_pos.min(entry.vertices.len());
            entry
                .vertices
                .insert(pos, HlrTopoBRepVData::with_vertex(p, v));
            self.vertex_iter_pos = pos + 1;
        }
    }

    /// `Append(V, P)` (`cxx:299-303`).
    pub fn append(&mut self, v: &TopoShape, p: f64) {
        let Some(key) = self.current_edge.clone() else {
            return;
        };
        if let Some(entry) = self.my_edges_vertices.get_mut(&key) {
            entry.vertices.push(HlrTopoBRepVData::with_vertex(p, v));
        }
    }
}
