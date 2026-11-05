//! The `TopTools_*` containers and the two spatial indexes used by
//! `BRepBuilderAPI_Sewing`.
//!
//! `TopTools_IndexedDataMapOfShapeShape` / `TopTools_DataMapOfShapeListOfShape`
//! / `TopTools_MapOfShape` all key by `TopoDS_Shape::IsSame`
//! (`TopTools_ShapeMapHasher.cxx:22-26`: TShape address + location), which is
//! [`shape_map_key`]. OCCT's indexed maps answer lookups by position with the
//! insertion order; the data maps iterate in bucket order, for which this port
//! keeps the insertion order as a deterministic stand-in.
//!
//! `NCollection_UBTree<Standard_Integer, Bnd_Box>`
//! (`BRepBuilderAPI_BndBoxTreeSelector.hxx:36-38`) is only an acceleration
//! structure for "all boxes not `IsOut` of the query box", so a scan answers
//! the same set; `BRepBuilderAPI_CellFilter` + `BRepBuilderAPI_VertexInspector`
//! (`BRepBuilderAPI_VertexInspector.cxx:5305-5325`) likewise reject a point
//! when any axis exceeds the tolerance, which the scan reproduces exactly.
//! Only the enumeration *order* differs, and it is fixed to insertion order.

use std::collections::HashMap;

use occt_core::bnd::BndBox;
use occt_core::gp::GpXyz;

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `TopTools_DataMapOfShape*` / `TopTools_IndexedDataMapOfShape*`.
#[derive(Clone, Default)]
pub(crate) struct ShapeKeyedMap<V> {
    items: Vec<(TopoShape, V)>,
    index: HashMap<ShapeMapKey, usize>,
}

impl<V> ShapeKeyedMap<V> {
    pub(crate) fn new() -> Self {
        Self {
            items: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// `NCollection_Map::Add`: insert when absent, answer whether it was added.
    pub(crate) fn add(&mut self, key: &TopoShape, value: V) -> bool {
        let k = shape_map_key(key);
        if self.index.contains_key(&k) {
            return false;
        }
        self.index.insert(k, self.items.len());
        self.items.push((key.clone(), value));
        true
    }

    /// `NCollection_DataMap::Bind`: insert or overwrite in place.
    pub(crate) fn bind(&mut self, key: &TopoShape, value: V) {
        let k = shape_map_key(key);
        match self.index.get(&k) {
            Some(&i) => self.items[i].1 = value,
            None => {
                self.index.insert(k, self.items.len());
                self.items.push((key.clone(), value));
            }
        }
    }

    pub(crate) fn contains(&self, key: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(key))
    }

    pub(crate) fn find(&self, key: &TopoShape) -> Option<&V> {
        let i = *self.index.get(&shape_map_key(key))?;
        self.items.get(i).map(|(_, v)| v)
    }

    /// `ChangeFind` / `ChangeFromKey`: answer the stored value, inserting the
    /// default first when the key is absent.
    pub(crate) fn change(&mut self, key: &TopoShape) -> &mut V
    where
        V: Default,
    {
        let k = shape_map_key(key);
        if !self.index.contains_key(&k) {
            let pos = self.items.len();
            self.index.insert(shape_map_key(key), pos);
            self.items.push((key.clone(), V::default()));
        }
        let i = *self.index.get(&shape_map_key(key)).unwrap();
        &mut self.items[i].1
    }

    /// `IsBound` + `(key)` read.
    pub(crate) fn get(&self, key: &TopoShape) -> Option<V>
    where
        V: Clone,
    {
        self.find(key).cloned()
    }

    /// `UnBind` / `RemoveKey`. Removal keeps the relative order of the rest.
    pub(crate) fn remove(&mut self, key: &TopoShape) {
        let k = shape_map_key(key);
        if let Some(i) = self.index.remove(&k) {
            self.items.remove(i);
            self.reindex_from(i);
        }
    }

    fn reindex_from(&mut self, from: usize) {
        self.index.clear();
        for (i, (key, _)) in self.items.iter().enumerate() {
            if i >= from {
                self.index.insert(shape_map_key(key), i);
            }
        }
        // Entries before `from` keep their position.
        for i in 0..from.min(self.items.len()) {
            self.index.insert(shape_map_key(&self.items[i].0), i);
        }
    }

    pub(crate) fn extent(&self) -> usize {
        self.items.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub(crate) fn clear(&mut self) {
        self.items.clear();
        self.index.clear();
    }

    /// `FindKey(i)` (1-based).
    pub(crate) fn key_at(&self, i: usize) -> Option<&TopoShape> {
        self.items.get(i.checked_sub(1)?).map(|(k, _)| k)
    }

    /// `FindFromIndex(i)` (1-based).
    pub(crate) fn value_at(&self, i: usize) -> Option<&V> {
        self.items.get(i.checked_sub(1)?).map(|(_, v)| v)
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &(TopoShape, V)> {
        self.items.iter()
    }

    /// Snapshot of the keys in insertion order (safe while mutating `self`).
    pub(crate) fn keys(&self) -> Vec<TopoShape> {
        self.items.iter().map(|(k, _)| k.clone()).collect()
    }
}

/// `TopTools_MapOfShape` / `TopTools_MapOfOrientedShape`: a set of shapes
/// (`TopTools_MapOfOrientedShape` additionally compares the orientation, which
/// the caller encodes in the key).
#[derive(Clone, Default)]
pub(crate) struct ShapeSet {
    items: Vec<(ShapeMapKey, TopoShape)>,
}

impl ShapeSet {
    pub(crate) fn new() -> Self {
        Self { items: Vec::new() }
    }

    pub(crate) fn add(&mut self, shape: &TopoShape) -> bool {
        let k = shape_map_key(shape);
        if self.items.iter().any(|(kk, _)| *kk == k) {
            return false;
        }
        self.items.push((k, shape.clone()));
        true
    }

    pub(crate) fn contains(&self, shape: &TopoShape) -> bool {
        let k = shape_map_key(shape);
        self.items.iter().any(|(kk, _)| *kk == k)
    }

    pub(crate) fn remove(&mut self, shape: &TopoShape) {
        let k = shape_map_key(shape);
        self.items.retain(|(kk, _)| *kk != k);
    }

    pub(crate) fn extent(&self) -> usize {
        self.items.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub(crate) fn clear(&mut self) {
        self.items.clear();
    }

    pub(crate) fn key_at(&self, i: usize) -> Option<&TopoShape> {
        self.items.get(i.checked_sub(1)?).map(|(_, s)| s)
    }

    pub(crate) fn keys(&self) -> Vec<TopoShape> {
        self.items.iter().map(|(_, s)| s.clone()).collect()
    }
}

/// `TColStd_SequenceOfShape` plus the parallel `TColStd_SequenceOfBoolean` /
/// `TColStd_SequenceOfReal` carried next to it in the merging stage.
#[derive(Clone, Default)]
pub(crate) struct ShapeSequence<T> {
    pub(crate) items: Vec<(TopoShape, T)>,
}

#[allow(dead_code)]
impl<T> ShapeSequence<T> {
    pub(crate) fn new() -> Self {
        Self { items: Vec::new() }
    }
    pub(crate) fn len(&self) -> usize {
        self.items.len()
    }
    pub(crate) fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    pub(crate) fn append(&mut self, shape: &TopoShape, extra: T) {
        self.items.push((shape.clone(), extra));
    }
    pub(crate) fn prepend(&mut self, shape: &TopoShape, extra: T) {
        self.items.insert(0, (shape.clone(), extra));
    }
    pub(crate) fn insert_before(&mut self, index: usize, shape: &TopoShape, extra: T) {
        let i = index.max(1) - 1;
        self.items.insert(i.min(self.items.len()), (shape.clone(), extra));
    }
    pub(crate) fn insert_after(&mut self, index: usize, shape: &TopoShape, extra: T) {
        let i = index.min(self.items.len());
        self.items.insert(i, (shape.clone(), extra));
    }
    /// `Remove(i)` (1-based).
    pub(crate) fn remove(&mut self, index: usize) {
        if index >= 1 && index <= self.items.len() {
            self.items.remove(index - 1);
        }
    }
    /// `Value(i)` (1-based).
    pub(crate) fn value(&self, index: usize) -> Option<&TopoShape> {
        self.items.get(index.checked_sub(1)?).map(|(s, _)| s)
    }
    pub(crate) fn extra(&self, index: usize) -> Option<&T> {
        self.items.get(index.checked_sub(1)?).map(|(_, t)| t)
    }
    pub(crate) fn set_extra(&mut self, index: usize, value: T) {
        if index >= 1 && index <= self.items.len() {
            self.items[index - 1].1 = value;
        }
    }
    pub(crate) fn first(&self) -> Option<&TopoShape> {
        self.items.first().map(|(s, _)| s)
    }
    pub(crate) fn clone_shapes(&self) -> Vec<TopoShape> {
        self.items.iter().map(|(s, _)| s.clone()).collect()
    }
}

/// `TColStd_SequenceOfReal`.
#[derive(Clone, Default)]
pub(crate) struct RealSequence {
    pub(crate) items: Vec<f64>,
}

#[allow(dead_code)]
impl RealSequence {
    pub(crate) fn new() -> Self {
        Self { items: Vec::new() }
    }
    pub(crate) fn len(&self) -> usize {
        self.items.len()
    }
    pub(crate) fn append(&mut self, v: f64) {
        self.items.push(v);
    }
    pub(crate) fn prepend(&mut self, v: f64) {
        self.items.insert(0, v);
    }
    pub(crate) fn insert_before(&mut self, index: usize, v: f64) {
        let i = index.max(1) - 1;
        self.items.insert(i.min(self.items.len()), v);
    }
    pub(crate) fn remove(&mut self, index: usize) {
        if index >= 1 && index <= self.items.len() {
            self.items.remove(index - 1);
        }
    }
    pub(crate) fn value(&self, index: usize) -> Option<f64> {
        self.items.get(index.checked_sub(1)?).copied()
    }
    pub(crate) fn set_value(&mut self, index: usize, v: f64) {
        if index >= 1 && index <= self.items.len() {
            self.items[index - 1] = v;
        }
    }
}

/// `TColStd_SequenceOfInteger`.
#[derive(Clone, Default)]
pub(crate) struct IntSequence {
    pub(crate) items: Vec<i32>,
}

#[allow(dead_code)]
impl IntSequence {
    pub(crate) fn new() -> Self {
        Self { items: Vec::new() }
    }
    pub(crate) fn len(&self) -> usize {
        self.items.len()
    }
    pub(crate) fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    pub(crate) fn append(&mut self, v: i32) {
        self.items.push(v);
    }
    pub(crate) fn prepend(&mut self, v: i32) {
        self.items.insert(0, v);
    }
    pub(crate) fn prepend_seq(&mut self, other: &IntSequence) {
        for v in other.items.iter().rev() {
            self.items.insert(0, *v);
        }
    }
    pub(crate) fn insert_before(&mut self, index: usize, v: i32) {
        let i = index.max(1) - 1;
        self.items.insert(i.min(self.items.len()), v);
    }
    pub(crate) fn remove(&mut self, index: usize) {
        if index >= 1 && index <= self.items.len() {
            self.items.remove(index - 1);
        }
    }
    pub(crate) fn value(&self, index: usize) -> Option<i32> {
        self.items.get(index.checked_sub(1)?).copied()
    }
    pub(crate) fn set(&mut self, index: usize, v: i32) {
        if index >= 1 && index <= self.items.len() {
            self.items[index - 1] = v;
        }
    }
    pub(crate) fn first(&self) -> Option<i32> {
        self.items.first().copied()
    }
}

/// `TColStd_IndexedMapOfInteger` (1-based, insertion order).
#[derive(Clone, Default)]
pub(crate) struct IntIndexedMap {
    pub(crate) items: Vec<i32>,
}

#[allow(dead_code)]
impl IntIndexedMap {
    pub(crate) fn new() -> Self {
        Self { items: Vec::new() }
    }
    pub(crate) fn add(&mut self, v: i32) -> bool {
        if self.items.contains(&v) {
            return false;
        }
        self.items.push(v);
        true
    }
    pub(crate) fn contains(&self, v: i32) -> bool {
        self.items.contains(&v)
    }
    pub(crate) fn extent(&self) -> usize {
        self.items.len()
    }
    pub(crate) fn last(&self) -> Option<i32> {
        self.items.last().copied()
    }
    pub(crate) fn value(&self, index: usize) -> Option<i32> {
        self.items.get(index.checked_sub(1)?).copied()
    }
    pub(crate) fn remove_last(&mut self) {
        self.items.pop();
    }
}

/// `NCollection_UBTree<Standard_Integer, Bnd_Box>`: boxes with their 1-based
/// insertion index (the index `BRepBuilderAPI_Sewing::Cutting` stores in
/// `myVertexNode`).
#[derive(Clone, Default)]
pub(crate) struct BoxTree {
    boxes: Vec<(usize, BndBox)>,
}

impl BoxTree {
    pub(crate) fn new() -> Self {
        Self { boxes: Vec::new() }
    }
    /// `NCollection_UBTreeFiller::Add(theIndex, theBox)`.
    pub(crate) fn add(&mut self, index: usize, b: BndBox) {
        self.boxes.push((index, b));
    }
    /// `UBTree::Select(selector)` with `Reject(box) = myBox.IsOut(box)`
    /// (`BRepBuilderAPI_BndBoxTreeSelector.hxx:44-50`).
    pub(crate) fn select(&self, current: &BndBox) -> Vec<usize> {
        self.boxes
            .iter()
            .filter(|(_, b)| !current.is_out_box(b))
            .map(|(i, _)| *i)
            .collect()
    }
}

/// `BRepBuilderAPI_VertexInspector` + `BRepBuilderAPI_CellFilter`
/// (`BRepBuilderAPI_VertexInspector.hxx:38-65`, `.cxx:5305-5325`).
///
/// The inspector keeps the searched points, the current point and the squared
/// tolerance; `Inspect` keeps a target when `|d| <= tol` on every axis.
pub(crate) struct VertexInspector {
    tol2: f64,
    points: Vec<GpXyz>,
    current: GpXyz,
    res: Vec<usize>,
}

impl VertexInspector {
    pub(crate) fn new(tol: f64) -> Self {
        Self {
            tol2: tol * tol,
            points: Vec::new(),
            current: GpXyz::new(0., 0., 0.),
            res: Vec::new(),
        }
    }

    /// `Add(theXYZ)`.
    pub(crate) fn add(&mut self, p: &GpXyz) {
        self.points.push(*p);
    }

    pub(crate) fn clear_res_list(&mut self) {
        self.res.clear();
    }

    pub(crate) fn set_current(&mut self, p: &GpXyz) {
        self.current = *p;
    }

    pub(crate) fn res_ind(&self) -> &[usize] {
        &self.res
    }

    /// `CellFilter::Inspect(pMin, pMax, inspector)`: the points whose cell
    /// range covers the query range, in stored order.
    pub(crate) fn inspect_range(&mut self, p_min: &GpXyz, p_max: &GpXyz) {
        for (i, q) in self.points.iter().enumerate() {
            if q.x >= p_min.x && q.x <= p_max.x && q.y >= p_min.y && q.y <= p_max.y
                && q.z >= p_min.z && q.z <= p_max.z
            {
                // `Inspect(theTarget)` keeps or rejects on the inspector test.
                let dx = self.current.x - q.x;
                let dy = self.current.y - q.y;
                let dz = self.current.z - q.z;
                if dx * dx <= self.tol2 && dy * dy <= self.tol2 && dz * dz <= self.tol2 {
                    self.res.push(i + 1);
                }
            }
        }
    }
}


