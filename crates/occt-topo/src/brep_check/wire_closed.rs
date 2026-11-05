//! `BRepCheck_Wire::Closed` (`BRepCheck_Wire.cxx:283-425`).
//!
//! The oriented-edge / oriented-vertex closure control:
//!
//! * every oriented edge of the wire is collected together with its oriented
//!   vertices (`cxx:320-351`);
//! * the edges reached from the first oriented edge by shared vertices must
//!   cover all collected edges, otherwise `NotConnected` (`cxx:353-373`);
//! * an edge used three times, or twice with the same orientation, is
//!   `RedundantEdge` (`cxx:375-404`);
//! * a vertex incident to an odd number of oriented edge occurrences is
//!   `NotClosed` (`cxx:406-417`).
//!
//! UNPORTED, and why: `myCdone` / `myCstat` / `myMap` keep the answer cached per
//! measured shape and `Closed` returns an already recorded non-`NoError` status
//! before re-measuring (`cxx:296-311`). This port has no `BRepCheck_Result`
//! cache — the analysis is redone on every call and only the status of this
//! measurement is returned, which is what a fresh `BRepCheck_Wire(S)` answers.
//! `myNbori` (`cxx:323`, `327`, commented out in OCCT) is not modelled.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_check::Status;
use crate::shape::{TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

/// `IsOriented(S)` (`BRepCheck_Wire.cxx:95-98`): `FORWARD` or `REVERSED`.
fn is_oriented(s: &TopoShape) -> bool {
    matches!(s.orientation(), Orientation::Forward | Orientation::Reversed)
}

/// An insertion-ordered set of shape keys — `NCollection_Map<TopoDS_Shape>`.
struct ShapeSet {
    order: Vec<ShapeMapKey>,
    seen: HashMap<ShapeMapKey, ()>,
}

impl ShapeSet {
    fn new() -> Self {
        Self { order: Vec::new(), seen: HashMap::new() }
    }

    fn add(&mut self, s: &TopoShape) {
        let k = shape_map_key(s);
        if self.seen.insert(k.clone(), ()).is_none() {
            self.order.push(k);
        }
    }

    fn contains(&self, s: &TopoShape) -> bool {
        self.seen.contains_key(&shape_map_key(s))
    }

    fn extent(&self) -> usize {
        self.order.len()
    }

    fn clear(&mut self) {
        self.order.clear();
        self.seen.clear();
    }
}

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// keyed by vertex, valued by the oriented edge occurrences touching it
/// (`myMapVE`, `BRepCheck_Wire.hxx:120`): index-based like `NCollection_*`
/// (`FindIndex` is 0 when absent, `Add` appends).
struct VertexEdges {
    keys: Vec<ShapeMapKey>,
    index: HashMap<ShapeMapKey, usize>,
    lists: Vec<Vec<TopoShape>>,
}

impl VertexEdges {
    fn new() -> Self {
        Self { keys: Vec::new(), index: HashMap::new(), lists: Vec::new() }
    }

    fn index_of(&self, v: &TopoShape) -> Option<usize> {
        self.index.get(&shape_map_key(v)).copied()
    }

    fn add(&mut self, v: &TopoShape, e: &TopoShape) -> usize {
        let k = shape_map_key(v);
        if let Some(i) = self.index.get(&k) {
            self.lists[*i].push(e.clone());
            return *i;
        }
        let i = self.lists.len();
        self.keys.push(k.clone());
        self.index.insert(k, i);
        self.lists.push(vec![e.clone()]);
        i
    }

    fn clear(&mut self) {
        self.keys.clear();
        self.index.clear();
        self.lists.clear();
    }
}

/// `Propagate(mapVE, edg, mapE)` (`BRepCheck_Wire.cxx:1777-1823`): flood the
/// edges reachable from `edg` through vertices listed in `map_ve`, adding every
/// reached edge to `map_e`.
fn propagate(ve: &VertexEdges, edg: &TopoShape, map_e: &mut ShapeSet) {
    let mut current: Vec<TopoShape> = vec![edg.clone()];
    while !current.is_empty() {
        let mut next: Vec<TopoShape> = Vec::new();
        for edge in &current {
            map_e.add(edge);
            for v in explorer_shapes(edge, ShapeType::Vertex) {
                let Some(i) = ve.index_of(&v) else { continue };
                for e in &ve.lists[i] {
                    if !crate::topo_tools::is_same_occurrence(e, edge) && !map_e.contains(e) {
                        map_e.add(e);
                        next.push(e.clone());
                    }
                }
            }
        }
        current = next;
    }
}

/// `BRepCheck_Wire::Closed(bool Update = false)` (`BRepCheck_Wire.cxx:283-425`).
pub fn wire_closed(w: &Wire) -> Status {
    let mut map_s = ShapeSet::new();
    // `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>> Cradoc`
    // (`cxx:317`): every occurrence of an oriented edge, keyed by shape.
    let mut cradoc: HashMap<ShapeMapKey, Vec<Orientation>> = HashMap::new();
    let mut map_ve = VertexEdges::new();

    let edges = explorer_shapes(&w.0, ShapeType::Edge);
    for e in &edges {
        if !is_oriented(e) {
            continue;
        }
        // `Cradoc(Current()).Append(Current())` (`cxx:325-334`).
        cradoc.entry(shape_map_key(e)).or_default().push(e.orientation());
        map_s.add(e);
        // `for (expv.Init(Current(), TopAbs_VERTEX); ...)` (`cxx:337-348`).
        for v in explorer_shapes(e, ShapeType::Vertex) {
            if is_oriented(&v) {
                map_ve.add(&v, e);
            }
        }
    }

    let the_nbori = map_s.extent();
    if the_nbori >= 2 {
        map_s.clear();
        // `for (exp.ReInit(); exp.More(); exp.Next()) if (IsOriented) break;`
        // (`cxx:354-361`).
        if let Some(first) = edges.iter().find(|e| is_oriented(e)) {
            propagate(&map_ve, first, &mut map_s);
        }
    }
    if the_nbori != map_s.extent() {
        return Status::NotConnected;
    }

    // `cxx:375-404`: an edge may appear at most twice, once FORWARD and once
    // REVERSED.
    for orients in cradoc.values() {
        if orients.len() >= 3 {
            return Status::RedundantEdge;
        }
        if orients.len() == 2 && orients[0] == orients[1] {
            return Status::RedundantEdge;
        }
    }

    // `cxx:406-417`: every vertex must be met an even number of times.
    for list in &map_ve.lists {
        if list.len() % 2 != 0 {
            return Status::NotClosed;
        }
    }

    Status::NoError
}
