//! `BRepOffset_Inter3d` (`BRepOffset_Inter3d.hxx` / `BRepOffset_Inter3d.cxx`):
//! the 3D intersections between offset faces (and between offset faces and
//! context faces), recorded in `BRepAlgo_AsDes`.
//!
//! Ported here so far (control flow taken from
//! `src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_Inter3d.cxx`):
//!
//! - file-local `ExtentEdge` (`:62-88`) - `extent_edge.rs`
//! - `CompletInt` (`:92-148`) - `complet_int.rs`
//! - `FaceInter` (`:155-254`) - `complet_int.rs`
//! - `ConnexIntByArc` (`:258-443`) - `connex_int_by_arc.rs`
//! - `ConnexIntByInt` (`:446-1041`) - `connex_int_by_int.rs`
//! - `ContextIntByInt` (`:1045-1272`) - `context_int_by_int.rs`
//! - `ContextIntByArc` (`:1276-1494`) - `context_int_by_arc.rs`
//! - `SetDone` (`:1500-1516`), `IsDone` (`:1519-1533`) and
//!   `Store` (`:1536-1554`) - `store.rs`
//!
//! `MapSF` is the `NCollection_DataMap<TopoDS_Shape, BRepOffset_Offset>` the
//! two integer-mode arms read; [`MapSf`] is its port, keyed on
//! `TopTools_ShapeMapHasher` identity and valued by
//! [`BRepOffsetOffset`](crate::brep_offset::BRepOffsetOffset). `MES` and
//! `Build` are plain `DataMap<TopoDS_Shape, TopoDS_Shape>` and shared
//! [`ShapeMap`]s.
//!
//! The 2D half of the integer mode - `BRepOffset_Inter2d::ConnexIntByInt`
//! (`BRepOffset_Inter2d.cxx:1830-2094`) and `ConnexIntByIntInVert`
//! (`:2100-2308`) - lives with that class
//! ([`connex_int_by_int`](crate::brep_offset::connex_int_by_int)); the
//! `ConnexIntByInt` here calls it after the 3D arm.
//!
//! Two reference arms call `BRepOffset_Tool::PipeInter`
//! (`BRepOffset_Tool.cxx:747-804`, ported as `BRepOffsetTool::pipe_inter` in
//! `brep_offset/tool/pipe_inter.rs`): the `InterPipes` arm of `FaceInter`
//! (`:191-217`, `:246-249`) and the tube/tube arm of `ConnexIntByArc`
//! (`:369-383`). Both are wired to that port.
//!
//! `myAsDes` is the shared `BRepAlgo_AsDes` handle of the class. The port owns
//! it by value ([`BRepOffsetInter3d::as_des`] / [`into_as_des`]) because the
//! Rust caller cannot alias a `&mut AsDes` for the whole lifetime of the
//! object; the semantics of every `Add` are unchanged.
//!
//! [`into_as_des`]: BRepOffsetInter3d::into_as_des

use std::collections::HashMap;

use crate::abs::State;
use crate::brep_algo::AsDes;
use crate::brep_offset::BRepOffsetOffset;
use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

mod complet_int;
mod connex_int_by_arc;
mod connex_int_by_int;
mod context_int_by_arc;
mod context_int_by_int;
mod extent_edge;
mod store;

/// `NCollection_DataMap<TopoDS_Shape, BRepOffset_Offset, TopTools_ShapeMapHasher>`
/// `MapSF` (`BRepOffset_Inter3d.cxx:447`): the offset built for every initial
/// face (`BRepOffset_MakeOffset::myMapSF`).
pub type MapSf = HashMap<ShapeMapKey, BRepOffsetOffset>;

/// `BRepOffset_Inter3d` (`BRepOffset_Inter3d.hxx:41-135`).
pub struct BRepOffsetInter3d {
    /// `myAsDes`.
    as_des: AsDes,
    /// `myTouched`: `NCollection_IndexedMap<TopoDS_Shape>` of the faces the
    /// intersections touched, insertion-ordered and `IsSame`-unique.
    touched: ShapeIndexedMap,
    /// `myDone`: `NCollection_DataMap<TopoDS_Shape, NCollection_List<
    /// TopoDS_Shape>>`, bound both ways by `SetDone`.
    done: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myNewEdges`: `NCollection_IndexedMap<TopoDS_Shape>` of the
    /// intersection edges `Store` recorded.
    new_edges: ShapeIndexedMap,
    /// `mySide`.
    side: State,
    /// `myTol`.
    tol: f64,
}

impl BRepOffsetInter3d {
    /// `BRepOffset_Inter3d(AsDes, Side, Tol)` (`cxx:53-59`).
    pub fn new(as_des: AsDes, side: State, tol: f64) -> Self {
        Self {
            as_des,
            touched: ShapeIndexedMap::new(),
            done: HashMap::new(),
            new_edges: ShapeIndexedMap::new(),
            side,
            tol,
        }
    }

    /// `AsDes()` (`hxx:113-115`).
    pub fn as_des(&self) -> &AsDes {
        &self.as_des
    }

    /// Mutable view of the owned `AsDes` handle.
    pub fn as_des_mut(&mut self) -> &mut AsDes {
        &mut self.as_des
    }

    /// Takes the `AsDes` back out (the Rust stand-in for dropping the last
    /// reference to the shared handle).
    pub fn into_as_des(self) -> AsDes {
        self.as_des
    }

    /// `TouchedFaces()` (`hxx:109-112`).
    pub fn touched_faces(&self) -> &ShapeIndexedMap {
        &self.touched
    }

    /// `NewEdges()` (`hxx:118`).
    pub fn new_edges(&self) -> &ShapeIndexedMap {
        &self.new_edges
    }
}

/// `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>`: insertion
/// ordered, `IsSame`-unique, 1-based lookup as in the reference.
#[derive(Debug, Clone, Default)]
pub struct ShapeIndexedMap {
    items: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
}

impl ShapeIndexedMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// `NCollection_IndexedMap::Add(S)`: appends and answers the 1-based index,
    /// or the existing index when `S` is already present.
    pub fn add(&mut self, s: &TopoShape) -> usize {
        let key = shape_map_key(s);
        if let Some(i) = self.index.get(&key) {
            return *i + 1;
        }
        self.items.push(s.clone());
        self.index.insert(key, self.items.len() - 1);
        self.items.len()
    }

    pub fn contains(&self, s: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(s))
    }

    /// `NCollection_Map::Remove(S)` / `NCollection_IndexedMap::Remove(S)`:
    /// drops `S` and answers whether it was present. The later indices shift
    /// down like [`ShapeIndexedMap::remove_key`] so `value(i)` stays dense.
    pub fn remove(&mut self, s: &TopoShape) -> bool {
        let key = shape_map_key(s);
        let Some(i) = self.index.remove(&key) else {
            return false;
        };
        self.items.remove(i);
        for idx in self.index.values_mut() {
            if *idx > i {
                *idx -= 1;
            }
        }
        true
    }

    /// `NCollection_IndexedMap::RemoveKey(S)` (`NCollection_IndexedMap.hxx`):
    /// drops `S` and shifts every later index down by one, so `value(i)` stays
    /// 1-based and dense. `RemoveCorks` (`BRepOffset_MakeOffset.cxx:732`) and
    /// `IsPlanar` (`cxx:5534-5537`) both rely on that shift: they remove a key
    /// and immediately `Add` it back so it moves to the end.
    pub fn remove_key(&mut self, s: &TopoShape) {
        let key = shape_map_key(s);
        let Some(i) = self.index.remove(&key) else {
            return;
        };
        self.items.remove(i);
        for idx in self.index.values_mut() {
            if *idx > i {
                *idx -= 1;
            }
        }
    }

    pub fn extent(&self) -> usize {
        self.items.len()
    }

    /// `NCollection_IndexedMap::Clear()` (`NCollection_IndexedMap.hxx`): drops
    /// every entry, after which `Add` numbers them from one again.
    pub fn clear(&mut self) {
        self.items.clear();
        self.index.clear();
    }

    /// `FindKey(i)` / `operator()(i)`, 1-based.
    pub fn value(&self, i: usize) -> Option<&TopoShape> {
        self.items.get(i.checked_sub(1)?)
    }

    pub fn iter(&self) -> impl Iterator<Item = &TopoShape> {
        self.items.iter()
    }
}
