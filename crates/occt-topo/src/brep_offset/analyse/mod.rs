//! `BRepOffset_Analyse` (`BRepOffset_Analyse.hxx`, `BRepOffset_Analyse.cxx`).
//!
//! Analyses a shape to find the parts of edges connecting convex, concave or
//! tangent faces. The connectivity of every edge is stored as a list of
//! [`Interval`]s; the ancestors map (vertex -> edges, edge -> faces) is the
//! input of the three followers of this class (`BRepOffset_Tool`,
//! `BRepOffset_Inter2d` / `Inter3d`, `BRepOffset_MakeLoops`).
//!
//! Split over three files:
//! - this one: the object, its maps and the result queries
//!   (`cxx:810-1131`);
//! - [`perform`]: `Perform` / `EdgeAnalyse` / `CheckMixedContinuity`
//!   (`cxx:81-367`);
//! - [`tangent`]: `TreatTangentFaces` (`cxx:371-806`).

use std::cell::RefCell;
use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::brep_offset::{BRepOffsetTool, Interval};
use crate::brep_tool::BRepTool;
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::{Compound, Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{
    explorer_shapes, is_equal, is_same, shape_map_key, ShapeMapKey,
};
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use occt_core::gp::GpVec;
use occt_core::precision::REAL_SMALL;
use crate::builder::TopoBuilder;

mod perform;
mod tangent;

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>`, the type of `myAncestors`.
///
/// The indexed map keeps the insertion order of the keys (`FindKey(i)` /
/// `operator(i)`, used by `Descendants` at `cxx:843-865`) and answers `Add`
/// with *insert if absent* (`NCollection_IndexedDataMap.hxx:1055-1076`: a bound
/// key keeps its item and only the index is returned).
#[derive(Debug, Default, Clone)]
pub(crate) struct ShapeListMap {
    keys: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
    values: Vec<Vec<TopoShape>>,
}

impl ShapeListMap {
    fn clear(&mut self) {
        self.keys.clear();
        self.index.clear();
        self.values.clear();
    }

    /// `NCollection_IndexedDataMap::Add(Key, Item)` (`hxx:1055-1076`): the
    /// item is stored only for a new key; answers the 1-based index.
    fn add(&mut self, key: &TopoShape, items: &[TopoShape]) -> usize {
        let k = shape_map_key(key);
        if let Some(&i) = self.index.get(&k) {
            return i + 1;
        }
        let idx = self.values.len();
        self.keys.push(key.clone());
        self.values.push(items.to_vec());
        self.index.insert(k, idx);
        idx + 1
    }

    /// `operator(i)` / `ChangeFromKey` with the key created when absent
    /// (`NCollection_IndexedDataMap::operator()`, `hxx:384-402`).
    fn change_or_insert(&mut self, key: &TopoShape) -> &mut Vec<TopoShape> {
        let k = shape_map_key(key);
        let idx = match self.index.get(&k) {
            Some(&i) => i,
            None => {
                let i = self.values.len();
                self.keys.push(key.clone());
                self.values.push(Vec::new());
                self.index.insert(k, i);
                i
            }
        };
        &mut self.values[idx]
    }

    /// `operator(i)` of a 1-based index.
    fn value_at(&mut self, i: usize) -> Option<&mut Vec<TopoShape>> {
        if i == 0 {
            return None;
        }
        self.values.get_mut(i - 1)
    }

    /// `FindKey(i)` of a 1-based index.
    fn key_at(&self, i: usize) -> Option<&TopoShape> {
        if i == 0 {
            return None;
        }
        self.keys.get(i - 1)
    }

    fn seek(&self, key: &TopoShape) -> Option<&Vec<TopoShape>> {
        let k = shape_map_key(key);
        self.index.get(&k).and_then(|&i| self.values.get(i))
    }

    fn extent(&self) -> usize {
        self.values.len()
    }
}

/// `BRepOffset_Analyse` (`hxx:52-181`).
#[derive(Debug, Default, Clone)]
pub struct Analyse {
    /// `myShape`
    shape: Option<TopoShape>,
    /// `myAngle`
    angle: f64,
    /// `myOffset`
    offset: f64,
    /// `myFaceOffsetMap`
    face_offset_map: HashMap<ShapeMapKey, f64>,
    /// `myDone`
    done: bool,
    /// `myMapEdgeType`
    map_edge_type: HashMap<ShapeMapKey, Vec<Interval>>,
    /// `myAncestors`
    ancestors: ShapeListMap,
    /// `myReplacement`: face -> (edge -> replacement edge).
    replacement: HashMap<ShapeMapKey, HashMap<ShapeMapKey, Edge>>,
    /// `myDescendants`, filled on the first query (`hxx:154-156`).
    descendants: RefCell<HashMap<ShapeMapKey, Vec<TopoShape>>>,
    /// `myNewFaces`
    new_faces: Vec<TopoShape>,
    /// `myGenerated`: edge / vertex -> generated face / edge.
    generated: HashMap<ShapeMapKey, TopoShape>,
}

impl Analyse {
    /// Empty constructor (`cxx:64-68`): `myOffset = 0`, `myDone = false`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `BRepOffset_Analyse(S, Angle)` (`cxx:74-77`).
    pub fn with_shape(shape: &TopoShape, angle: f64) -> Result<Self, String> {
        let mut a = Self::new();
        a.perform(shape, angle)?;
        Ok(a)
    }

    /// `IsDone()` (`hxx:83`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `SetOffsetValue(theOffset)` (`hxx:137`).
    pub fn set_offset_value(&mut self, offset: f64) {
        self.offset = offset;
    }

    /// `SetFaceOffsetMap(theMap)` (`hxx:140-146`).
    pub fn set_face_offset_map(&mut self, map: HashMap<ShapeMapKey, f64>) {
        self.face_offset_map = map;
    }

    /// `NewFaces()` (`hxx:152`).
    pub fn new_faces(&self) -> &[TopoShape] {
        &self.new_faces
    }

    /// `HasAncestor(theS)` (`hxx:123`).
    pub fn has_ancestor(&self, s: &TopoShape) -> bool {
        self.ancestors.seek(s).is_some()
    }

    /// `Ancestors(theS)` (`hxx:126-129`).
    ///
    /// `FindFromKey` asserts on an absent key in OCCT; the port answers the
    /// empty list instead.
    pub fn ancestors(&self, s: &TopoShape) -> &[TopoShape] {
        self.ancestors.seek(s).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// `Type(theE)` (`cxx:891-894`).
    ///
    /// `myMapEdgeType(E)` asserts on an absent key in OCCT. The port answers
    /// the empty slice; every caller either tests `IsEmpty` first
    /// (`cxx:1073`, `:1112`) or only reaches bound edges.
    pub fn type_of(&self, e: &Edge) -> &[Interval] {
        self.map_edge_type
            .get(&shape_map_key(&e.0))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// `Edges(V, T, LE)` (`cxx:898-931`): the edges of type `T` on vertex `V`.
    ///
    /// `V1` (the first vertex of the edge) reads `Last()` of the interval list
    /// and `V2` reads `First()`, exactly as `cxx:916` / `:923` do; a closed
    /// edge whose two endpoints are `V` is appended twice.
    pub fn edges_of_vertex(&self, v: &Vertex, t: TypeOfConcavity) -> Vec<Edge> {
        let mut out = Vec::new();
        for anc in self.ancestors(&v.0) {
            if anc.shape_type() != ShapeType::Edge {
                continue;
            }
            let e = Edge(anc.clone());
            let intervals = self.type_of(&e);
            if intervals.is_empty() {
                continue;
            }
            let (v1, v2) = BRepOffsetTool::edge_vertices(&e);
            if v1.as_ref().is_some_and(|x| is_same(&x.0, &v.0))
                && intervals[intervals.len() - 1].kind() == t
            {
                out.push(e.clone());
            }
            if v2.as_ref().is_some_and(|x| is_same(&x.0, &v.0)) && intervals[0].kind() == t {
                out.push(e);
            }
        }
        out
    }

    /// `Edges(F, T, LE)` (`cxx:934-955`): the edges of type `T` on face `F`,
    /// once per matching interval (`cxx:945-952` has no `break`).
    pub fn edges_of_face(&self, f: &Face, t: TypeOfConcavity) -> Vec<Edge> {
        let mut out = Vec::new();
        for cur in explorer_shapes(&f.0, ShapeType::Edge) {
            let e = Edge(cur);
            for interval in self.type_of(&e) {
                if interval.kind() == t {
                    out.push(e.clone());
                }
            }
        }
        out
    }

    /// `TangentEdges(Edge, Vertex, Edges)` (`cxx:959-1001`): the edges of the
    /// shape that are tangent to `edge` at `vertex`.
    ///
    /// The early return of `cxx:972-975` happens *before* `Edges.Clear()`
    /// (`cxx:977`), so the caller's list is left untouched then - the port
    /// keeps that by returning without touching `edges`.
    pub fn tangent_edges(&self, edge: &Edge, vertex: &Vertex, edges: &mut Vec<TopoShape>) {
        let u_ref = BRepTool::parameter_of_vertex(&edge.0, &vertex.0);
        let c3d_ref = CurveAdaptor::of_edge(edge);
        let Some((_, mut v_ref)) = c3d_ref.d1(u_ref) else {
            return;
        };
        correct_orientation_of_tangent(&mut v_ref, vertex, edge);
        if v_ref.square_magnitude() < REAL_SMALL {
            return;
        }

        edges.clear();

        for anc in self.ancestors(&vertex.0) {
            if anc.shape_type() != ShapeType::Edge {
                continue;
            }
            let cur_e = Edge(anc.clone());
            if is_same(&cur_e.0, &edge.0) {
                continue;
            }
            let u = BRepTool::parameter_of_vertex(&cur_e.0, &vertex.0);
            let c3d = CurveAdaptor::of_edge(&cur_e);
            let Some((_, mut v)) = c3d.d1(u) else {
                continue;
            };
            correct_orientation_of_tangent(&mut v, vertex, &cur_e);
            if v.square_magnitude() < REAL_SMALL {
                continue;
            }
            if v.is_opposite_ang(&v_ref, self.angle) {
                edges.push(cur_e.0);
            }
        }
    }

    /// `Explode(List, T)` (`cxx:1005-1027`): one compound per cloud of faces
    /// whose connex edges are all of type `T`.
    pub fn explode(&self, t: TypeOfConcavity) -> Vec<TopoShape> {
        let mut list = Vec::new();
        let builder = TopoBuilder::new();
        let mut map: Vec<ShapeMapKey> = Vec::new();
        let Some(shape) = self.shape.as_ref() else {
            return list;
        };
        for cur in explorer_shapes(shape, ShapeType::Face) {
            if map.contains(&shape_map_key(&cur)) {
                continue;
            }
            map.push(shape_map_key(&cur));
            let face = Face(cur);
            let mut co = builder.make_compound_of(&[]);
            builder.add_compound(&mut co, &face.0);
            self.add_faces(&face, &mut co, &mut map, t);
            list.push(co.0);
        }
        list
    }

    /// `Explode(List, T1, T2)` (`cxx:1031-1054`).
    pub fn explode_two(&self, t1: TypeOfConcavity, t2: TypeOfConcavity) -> Vec<TopoShape> {
        let mut list = Vec::new();
        let builder = TopoBuilder::new();
        let mut map: Vec<ShapeMapKey> = Vec::new();
        let Some(shape) = self.shape.as_ref() else {
            return list;
        };
        for cur in explorer_shapes(shape, ShapeType::Face) {
            if map.contains(&shape_map_key(&cur)) {
                continue;
            }
            map.push(shape_map_key(&cur));
            let face = Face(cur);
            let mut co = builder.make_compound_of(&[]);
            builder.add_compound(&mut co, &face.0);
            self.add_faces_two(&face, &mut co, &mut map, t1, t2);
            list.push(co.0);
        }
        list
    }

    /// `AddFaces(Face, Co, Map, T)` (`cxx:1058-1092`): adds to `co` the faces
    /// reachable from `face` through edges whose first interval is of type `T`.
    ///
    /// `BRepOffset_MakeOffset::UpdateFaceOffset` (`BRepOffset_MakeOffset.cxx:2370`)
    /// is the public caller: it passes a fresh `Map` so the recursion only
    /// dedups within one call.
    pub fn add_faces(
        &self,
        face: &Face,
        co: &mut Compound,
        map: &mut Vec<ShapeMapKey>,
        t: TypeOfConcavity,
    ) {
        let builder = TopoBuilder::new();
        let Some(descendants) = self.descendants(&face.0, false) else {
            return;
        };
        for anc in descendants {
            if anc.shape_type() != ShapeType::Edge {
                continue;
            }
            let e = Edge(anc);
            let intervals = self.type_of(&e);
            if intervals.is_empty() || intervals[0].kind() != t {
                continue;
            }
            let l = self.ancestors(&e.0);
            if l.len() != 2 {
                continue;
            }
            let mut f1 = Face(l[0].clone());
            if is_same(&f1.0, &face.0) {
                f1 = Face(l[1].clone());
            }
            if map.contains(&shape_map_key(&f1.0)) {
                continue;
            }
            map.push(shape_map_key(&f1.0));
            builder.add_compound(co, &f1.0);
            self.add_faces(&f1, co, map, t);
        }
    }

    /// `AddFaces(Face, Co, Map, T1, T2)` (`cxx:1096-1131`).
    pub fn add_faces_two(
        &self,
        face: &Face,
        co: &mut Compound,
        map: &mut Vec<ShapeMapKey>,
        t1: TypeOfConcavity,
        t2: TypeOfConcavity,
    ) {
        let builder = TopoBuilder::new();
        let Some(descendants) = self.descendants(&face.0, false) else {
            return;
        };
        for anc in descendants {
            if anc.shape_type() != ShapeType::Edge {
                continue;
            }
            let e = Edge(anc);
            let intervals = self.type_of(&e);
            if intervals.is_empty() || (intervals[0].kind() != t1 && intervals[0].kind() != t2) {
                continue;
            }
            let l = self.ancestors(&e.0);
            if l.len() != 2 {
                continue;
            }
            let mut f1 = Face(l[0].clone());
            if is_same(&f1.0, &face.0) {
                f1 = Face(l[1].clone());
            }
            if map.contains(&shape_map_key(&f1.0)) {
                continue;
            }
            map.push(shape_map_key(&f1.0));
            builder.add_compound(co, &f1.0);
            self.add_faces_two(&f1, co, map, t1, t2);
        }
    }

    /// `Generated(theS)` (`cxx:831-835`): the face generated from an edge, or
    /// the edge generated from a vertex. `None` is the static null shape.
    pub fn generated(&self, s: &TopoShape) -> Option<TopoShape> {
        self.generated.get(&shape_map_key(s)).cloned()
    }

    /// `HasGenerated(theS)` (`hxx:163`).
    pub fn has_generated(&self, s: &TopoShape) -> bool {
        self.generated.contains_key(&shape_map_key(s))
    }

    /// `EdgeReplacement(theFace, theEdge)` (`cxx:810-827`): the replacement of
    /// the edge in the face, or the edge itself when there is none.
    pub fn edge_replacement(&self, face: &Face, edge: &Edge) -> Edge {
        match self.replacement.get(&shape_map_key(&face.0)) {
            Some(m) => match m.get(&shape_map_key(&edge.0)) {
                Some(e) => e.clone(),
                None => edge.clone(),
            },
            None => edge.clone(),
        }
    }

    /// `Descendants(theS, theUpdate)` (`cxx:840-870`): the shapes that have
    /// `theS` in their ancestor list, built from `myAncestors` on the first
    /// query (or every query when `theUpdate`).
    ///
    /// `None` is the null pointer `myDescendants.Seek` answers for a shape no
    /// ancestor list points at.
    pub fn descendants(&self, s: &TopoShape, update: bool) -> Option<Vec<TopoShape>> {
        {
            let mut cache = self.descendants.borrow_mut();
            if cache.is_empty() || update {
                cache.clear();
                let n = self.ancestors.extent();
                for i in 1..=n {
                    let Some(a_ss) = self.ancestors.key_at(i) else {
                        continue;
                    };
                    let a_l = self.ancestors.seek(a_ss).cloned().unwrap_or_default();
                    for a_sa in a_l {
                        let k = shape_map_key(&a_sa);
                        let list = cache.entry(k).or_default();
                        if !list.iter().any(|x| is_same(x, a_ss)) {
                            list.push(a_ss.clone());
                        }
                    }
                }
            }
            cache.get(&shape_map_key(s)).cloned()
        }
    }

    /// `Clear()` (`cxx:874-885`).
    pub fn clear(&mut self) {
        self.done = false;
        self.shape = None;
        self.map_edge_type.clear();
        self.ancestors.clear();
        self.face_offset_map.clear();
        self.replacement.clear();
        self.descendants.borrow_mut().clear();
        self.new_faces.clear();
        self.generated.clear();
    }

    // -----------------------------------------------------------------------
    // helpers shared with `perform` / `tangent`
    // -----------------------------------------------------------------------

    /// `myMapEdgeType(aE)` (`NCollection_DataMap::operator()`, insertion on
    /// demand).
    pub(super) fn edge_types_mut(&mut self, e: &Edge) -> &mut Vec<Interval> {
        self.map_edge_type
            .entry(shape_map_key(&e.0))
            .or_default()
    }

    /// `myMapEdgeType.Seek(aE)`.
    pub(super) fn edge_types(&self, e: &Edge) -> Option<&Vec<Interval>> {
        self.map_edge_type.get(&shape_map_key(&e.0))
    }

    /// `myMapEdgeType.Bind(aE, List())`: replaces the stored list.
    pub(super) fn bind_edge_types(&mut self, e: &Edge, list: Vec<Interval>) {
        self.map_edge_type.insert(shape_map_key(&e.0), list);
    }

    /// `myAncestors.ChangeFromKey(...)`.
    pub(super) fn ancestors_mut(&mut self, s: &TopoShape) -> &mut Vec<TopoShape> {
        self.ancestors.change_or_insert(s)
    }

    /// `myAncestors.Add(key, items)` followed by `myAncestors(index)`.
    pub(super) fn ancestors_add(&mut self, key: &TopoShape, items: &[TopoShape]) -> usize {
        self.ancestors.add(key, items)
    }

    /// `myAncestors(index)`.
    pub(super) fn ancestors_at(&mut self, i: usize) -> Option<&mut Vec<TopoShape>> {
        self.ancestors.value_at(i)
    }

    /// `myGenerated.Bind(key, value)`.
    pub(super) fn bind_generated(&mut self, key: &TopoShape, value: TopoShape) {
        self.generated.insert(shape_map_key(key), value);
    }

    /// `myReplacement.ChangeSeek(face)` / `myReplacement.Bound(face, ...)`.
    pub(super) fn replacement_map_mut(&mut self, face: &Face) -> &mut HashMap<ShapeMapKey, Edge> {
        self.replacement
            .entry(shape_map_key(&face.0))
            .or_default()
    }

    /// `myNewFaces.Append(face)`.
    pub(super) fn push_new_face(&mut self, face: &Face) {
        self.new_faces.push(face.0.clone());
    }
}

/// `CorrectOrientationOfTangent(TangVec, aVertex, anEdge)` (`cxx:46-55`): the
/// tangent is reversed when `aVertex` is the *last* vertex of the edge.
fn correct_orientation_of_tangent(tang: &mut GpVec, vertex: &Vertex, edge: &Edge) {
    let (_, last) = crate::topo_tools::edge_vertices(edge);
    if last.as_ref().is_some_and(|v| is_same(&v.0, &vertex.0)) {
        tang.reverse();
    }
}
