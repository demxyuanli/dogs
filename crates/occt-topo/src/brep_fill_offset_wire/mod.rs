//! `BRepFill_OffsetWire` — the offset of a planar wire (or of a planar face's
//! wires).
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_OffsetWire.cxx:17-2940`
//! and `.hxx:36-137` (OCCT 8.0.0).
//!
//! The class computes the locus of the points at distance `Offset` from the
//! spine, on the `MAT_Left` side (or on the outer side for a closed wire), by
//! walking the bisecting locus of the spine and offsetting each edge by the
//! 2D bisector between it and its neighbours:
//!
//! * `Init` (`cxx:327-366`): copy the forward spine, cut it at curvature
//!   extrema ([`spine::prepare_spine`], `cxx:1195-1290`) and compute the
//!   bisecting locus plus the topological link.
//! * `Perform` (`cxx:454-639`): the `KPartCircle` shortcut for complete
//!   circles, then the retry of the locus on a subdivided spine when
//!   `CheckBadEdges` finds edges the locus cannot see.
//! * `PerformWithBiLo` (`cxx:680-1184`): the body — circles and offset curves
//!   on the spine (`MakeCircle` / `MakeOffset`), vertices on the parallels from
//!   the bisectors (`TrimEdge`, `VertexFromNode`), then `MakeWires` and
//!   `FixHoles`.
//!
//! Ported layout:
//! * this file — the type, the `myMap` container, the constructors and the
//!   accessors,
//! * [`spine`] — `PrepareSpine` (`cxx:1195-1290`) with `CutEdge`
//!   (`cxx:1930-2070`) and `CutCurve` (`cxx:2076-2125`).
//!
//! Deviations (handle correspondence only, no missing branch):
//! * `TopoDS_Shape` may be null in OCCT and is `Option<TopoShape>` here
//!   (`shape.rs:69-73`: every `TopoShape` is non-null), so `mySpine`,
//!   `myWorkSpine` and `myShape` are `Option`.
//! * `myMap` keeps OCCT's `NCollection_IndexedDataMap` semantics in
//!   [`IndexedShapeMap`]: insertion keeps the first value bound to a key and
//!   the slot order is the insertion order, which `MakeWires` / `MakeOffset`
//!   rely on (`FindKey(1)`, `RemoveFromIndex`, ...).
//! * `BRepFill_OffsetWire::Add` (`hxx:103`) has no definition in
//!   `BRepFill_OffsetWire.cxx` — it is declared and never defined or called, so
//!   it has no port.
//! * `GeneratedShapes` returns an owned `Vec` instead of OCCT's reference to a
//!   static empty `NCollection_List`.

use std::collections::HashMap;

use occt_core::kernel::geomabs::JoinType;

use crate::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use crate::brep_mat2d_link_topo_bilo::BRepMat2dLinkTopoBilo;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{shape_map_key, ShapeMapKey};

pub mod curves;
pub mod detromp;
pub mod fix_holes;
pub mod init;
pub mod make_wires;
pub mod maps;
pub mod perform;
pub mod perform_with_bilo;
pub mod spine;
pub mod trim_edge;
pub mod vertices;

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// (`BRepFill_OffsetWire.hxx:128`).
///
/// `NCollection_IndexedDataMap::Add` keeps the value already bound to a key and
/// only appends a slot for a new key (`NCollection_IndexedDataMap.hxx`,
/// `addImpl`), which is what [`IndexedShapeMap::add`] does; every other method
/// mirrors the one OCCT calls.
#[derive(Clone, Default)]
pub struct IndexedShapeMap {
    /// `(key, key shape, list)` in slot order; slot `i` is `items[i - 1]`
    /// (1-based like OCCT).
    items: Vec<(ShapeMapKey, TopoShape, Vec<TopoShape>)>,
}

impl IndexedShapeMap {
    /// An empty map.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.items.clear();
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.items.len()
    }

    /// `Add(K, T)`: append `key` unless it is already bound, in which case the
    /// existing value is kept.
    pub fn add(&mut self, key: &TopoShape, value: Vec<TopoShape>) {
        let k = shape_map_key(key);
        if !self.items.iter().any(|(existing, _, _)| *existing == k) {
            self.items.push((k, key.clone(), value));
        }
    }

    /// `Contains(K)`.
    pub fn contains(&self, key: &TopoShape) -> bool {
        let k = shape_map_key(key);
        self.items.iter().any(|(existing, _, _)| *existing == k)
    }

    /// `FindFromKey(K)`; `None` where OCCT raises `Standard_NoSuchObject`.
    pub fn find_from_key(&self, key: &TopoShape) -> Option<&Vec<TopoShape>> {
        let k = shape_map_key(key);
        self.items
            .iter()
            .find(|(existing, _, _)| *existing == k)
            .map(|(_, _, v)| v)
    }

    /// `ChangeFromKey(K)`; `None` for an unbound key.
    pub fn change_from_key(&mut self, key: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        let k = shape_map_key(key);
        self.items
            .iter_mut()
            .find(|(existing, _, _)| *existing == k)
            .map(|(_, _, v)| v)
    }

    /// `RemoveKey(K)`.
    pub fn remove_key(&mut self, key: &TopoShape) {
        let k = shape_map_key(key);
        self.items.retain(|(existing, _, _)| *existing != k);
    }

    /// `FindKey(Index)` (`index` is 1-based).
    pub fn find_key(&self, index: usize) -> Option<&TopoShape> {
        self.items.get(index.checked_sub(1)?).map(|(_, s, _)| s)
    }

    /// `operator()(Index)` (`index` is 1-based).
    pub fn value(&self, index: usize) -> Option<&Vec<TopoShape>> {
        self.items.get(index.checked_sub(1)?).map(|(_, _, v)| v)
    }

    /// `ChangeFromIndex(Index)` (`index` is 1-based).
    pub fn value_mut(&mut self, index: usize) -> Option<&mut Vec<TopoShape>> {
        self.items.get_mut(index.checked_sub(1)?).map(|(_, _, v)| v)
    }

    /// `RemoveFromIndex(Index)` (`index` is 1-based).
    pub fn remove_from_index(&mut self, index: usize) {
        if index >= 1 && index <= self.items.len() {
            self.items.remove(index - 1);
        }
    }

    /// `RemoveLast()`.
    pub fn remove_last(&mut self) {
        self.items.pop();
    }

    /// The slot keys in slot order, for the drivers that walk the map by index.
    pub fn keys(&self) -> impl Iterator<Item = &TopoShape> {
        self.items.iter().map(|(_, s, _)| s)
    }
}

/// `class BRepFill_OffsetWire` (`BRepFill_OffsetWire.hxx:48-137`).
pub struct BRepFillOffsetWire {
    /// `mySpine` (`.hxx:121`).
    my_spine: Option<Face>,
    /// `myWorkSpine` (`.hxx:122`).
    my_work_spine: Option<Face>,
    /// `myOffset` (`.hxx:123`).
    my_offset: f64,
    /// `myIsOpenResult` (`.hxx:124`).
    my_is_open_result: bool,
    /// `myShape` (`.hxx:125`).
    my_shape: Option<TopoShape>,
    /// `myIsDone` (`.hxx:126`).
    my_is_done: bool,
    /// `myJoinType` (`.hxx:127`).
    my_join_type: JoinType,
    /// `myMap` (`.hxx:128`).
    my_map: IndexedShapeMap,
    /// `myBilo` (`.hxx:129`).
    my_bilo: BRepMat2dBisectingLocus,
    /// `myLink` (`.hxx:130`).
    my_link: BRepMat2dLinkTopoBilo,
    /// `myMapSpine` (`.hxx:131`): `NCollection_DataMap<TopoDS_Shape,
    /// TopoDS_Shape, TopTools_ShapeMapHasher>`, kept as `(key shape, value)`.
    my_map_spine: HashMap<ShapeMapKey, (TopoShape, TopoShape)>,
    /// `myCallGen` (`.hxx:132`).
    my_call_gen: bool,
}

impl BRepFillOffsetWire {
    /// `BRepFill_OffsetWire()` (`cxx:310-314`): `myIsOpenResult` and
    /// `myIsDone` are `false`; `myJoinType` is left indeterminate by OCCT, the
    /// port uses the `GeomAbs_Arc` default that every caller passes
    /// (`hxx:54-55`).
    pub fn new() -> Self {
        Self {
            my_spine: None,
            my_work_spine: None,
            my_offset: 0.0,
            my_is_open_result: false,
            my_shape: None,
            my_is_done: false,
            my_join_type: JoinType::Arc,
            my_map: IndexedShapeMap::new(),
            my_bilo: BRepMat2dBisectingLocus::new(),
            my_link: BRepMat2dLinkTopoBilo::new(),
            my_map_spine: HashMap::new(),
            my_call_gen: false,
        }
    }

    /// `IsDone()` (`cxx:369-372`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `Spine()` (`cxx:376-379`).
    pub fn spine(&self) -> Option<&Face> {
        self.my_spine.as_ref()
    }

    /// `Shape()` (`cxx:383-386`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `JoinType()` (`cxx:447-450`).
    pub fn join_type(&self) -> JoinType {
        self.my_join_type
    }

    /// `GeneratedShapes(SpineShape)` (`cxx:390-444`).
    ///
    /// The first call folds `myMapSpine` into `myMap` so that the shapes
    /// reported for the cut spine are reported for the original spine edges
    /// (`cxx:395-431`).
    pub fn generated_shapes(&mut self, spine_shape: &TopoShape) -> Vec<TopoShape> {
        if !self.my_call_gen {
            if !self.my_map_spine.is_empty() {
                // `cxx:398-430`, in `myMapSpine` iteration order. OCCT walks a
                // `TopTools_ShapeMapHasher` DataMap, the port sorts nothing and
                // folds the pairs one at a time, which reaches the same folds
                // because each pair only moves one slot.
                let pairs: Vec<(TopoShape, TopoShape)> =
                    self.my_map_spine.values().cloned().collect();
                for (key, value) in pairs {
                    if self.my_map.contains(&key) {
                        if !self.my_map.contains(&value) {
                            self.my_map.add(&value, Vec::new());
                        }
                        if !crate::topo_tools::is_same(&value, &key) {
                            let moved = self.my_map.find_from_key(&key).cloned().unwrap_or_default();
                            if let Some(slot) = self.my_map.change_from_key(&value) {
                                slot.extend(moved);
                            }
                            self.my_map.remove_key(&key);
                        }
                    }
                    let key_rev = key.oriented(key.orientation().reversed());
                    let value_rev = value.oriented(value.orientation().reversed());
                    if self.my_map.contains(&key_rev) {
                        if !self.my_map.contains(&value_rev) {
                            self.my_map.add(&value_rev, Vec::new());
                        }
                        if !crate::topo_tools::is_same(&value, &key) {
                            let moved = self
                                .my_map
                                .find_from_key(&key_rev)
                                .cloned()
                                .unwrap_or_default();
                            if let Some(slot) = self.my_map.change_from_key(&value_rev) {
                                slot.extend(moved);
                            }
                            self.my_map.remove_key(&key_rev);
                        }
                    }
                }
            }
            self.my_call_gen = true;
        }

        self.my_map
            .find_from_key(spine_shape)
            .cloned()
            .unwrap_or_default()
    }

    /// `Generated()` (`cxx:1187-1191`).
    pub fn generated(&mut self) -> &mut IndexedShapeMap {
        &mut self.my_map
    }
}

