//! `ShapeAnalysis_Shell` (`ShapeAnalysis_Shell.cxx:26-321`,
//! `ShapeAnalysis_Shell.hxx:26-131`) and `ShapeAnalysis_FreeBounds`
//! (`ShapeAnalysis_FreeBounds.cxx:56-648`,
//! `ShapeAnalysis_FreeBounds.hxx:58-233`), plus the two `ShapeExtend` classes
//! they are built on: `ShapeExtend_WireData`
//! (`ShapeExtend_WireData.cxx:40-700`) and the `ShapeAnalysis_Wire` slice the
//! free-bound chaining uses (`ShapeAnalysis_Wire.cxx:201-204`, `:693-761`).
//!
//! `ShapeAnalysis_Shell` splits the edges of a shell into three buckets
//! (`cxx:66-69`): an edge seen twice in the same orientation is *bad*, an edge
//! seen once is *free*, and an edge seen in both orientations is a *connected*
//! (internal) edge. `ShapeAnalysis_FreeBounds` then chains the free edges into
//! wires and separates the closed from the open ones.
//!
//! Port notes:
//!
//! * `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>` keeps its
//!   keys in insertion order and compares them with `IsSame`
//!   (`TopTools_ShapeMapHasher::IsEqual`, orientation included in neither the
//!   hash nor the equality). [`IndexedShapeMap`] is that container: a `Vec`
//!   searched with [`shape_map_key`].
//! * `TopoDS_Iterator(S)` (the default `cumOri = cumLoc = true`) is
//!   [`crate::iterator::cumulated_children`]; the `CumOri = false` form used by
//!   `TopExp::Vertices(E, V1, V2)` (`TopExp.cxx:214-251`) reads the stored
//!   children instead.
//! * `TopExp::Vertices(W, V1, V2)` (`TopExp.cxx:255-319`) is
//!   [`topexp_vertices_of_wire`]; `TopExp::FirstVertex` / `LastVertex`
//!   (`TopExp.cxx:182-210`) are [`super::first_vertex`] /
//!   [`super::last_vertex`], which compose the edge orientation the way
//!   `ShapeAnalysis_Edge::FirstVertex` / `LastVertex` do.
//! * `ShapeExtend_WireData` is [`WireData`]. The `mySeamF` / `mySeamR` /
//!   `mySeams` fields only feed `ShapeFix_Wire` and are not modelled; the
//!   `mySeamF = -1` invalidation the mutators perform has no consumer here.
//! * `ShapeAnalysis_Wire::CheckConnected` (`ShapeAnalysis_Wire.cxx:693-761`) is
//!   [`WireAnalyzer::check_connected`]. `myStatus` is `i32` with the
//!   `ShapeExtend::EncodeStatus` bits.
//! * `NCollection_UBTree` / `NCollection_UBTreeFiller` /
//!   `ShapeAnalysis_BoxBndTreeSelector` are
//!   [`super::box_bnd_tree`]; `mySeq` is passed into
//!   [`UbTree::select`] instead of being a selector field, because the OCCT
//!   array is mutated (a wire is reversed at `cxx:308`) while the selector
//!   holds it.
//! * `ShapeBuild_Vertex::CombineVertex` is [`super::combine_vertex`];
//!   `ShapeBuild_Edge::CopyReplaceVertices` is
//!   [`super::copy_replace_vertices_with`].
//! * `NCollection_HSequence<TopoDS_Shape>` becomes `Vec`; the 1-based OCCT
//!   indexing is kept in the comments (`Value(i)` -> `[i - 1]`).
//! * `ShapeExtend_Explorer::SeqFromCompound` / `CompoundFromSeq`
//!   (`ShapeExtend_Explorer.cxx:32-86`) are [`seq_from_compound`] and
//!   [`compound_from_seq`].
//! * UNPORTED: nothing in the two classes' reachable control flow. The
//!   `BRepBuilderAPI_Sewing` constructor (`cxx:62-95`) is ported through
//!   [`Sewing`]; `myTolerance` / `myShared` / `mySplitClosed` / `mySplitOpen`
//!   keep their OCCT names as `tolerance` / `shared` / `split_closed` /
//!   `split_open`.

use std::collections::HashMap;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_api_sewing::Sewing;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::iterator::cumulated_children;
use crate::shape::{Compound, Edge, Shell, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, ShapeMapKey};

use super::box_bnd_tree::{
    ShapeAnalysisBoxBndTreeSelector, UbTree, UbTreeFiller, SE_DONE1, SE_FAIL, SE_OK,
};
use super::combine_vertex;
use super::copy_replace_vertices_with;

/// `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>`
/// (`NCollection_IndexedMap.hxx`): insertion ordered, `IsSame` keyed.
#[derive(Default)]
struct IndexedShapeMap {
    items: Vec<TopoShape>,
}

impl IndexedShapeMap {
    fn new() -> Self {
        Self { items: Vec::new() }
    }

    /// `FindIndex(S)`: the 1-based position, `0` when absent.
    fn find_index(&self, shape: &TopoShape) -> usize {
        let key = shape_map_key(shape);
        match self.items.iter().position(|s| shape_map_key(s) == key) {
            Some(i) => i + 1,
            None => 0,
        }
    }

    /// `Add(S)`: `true` when the key was new.
    fn add(&mut self, shape: &TopoShape) -> bool {
        if self.find_index(shape) != 0 {
            return false;
        }
        self.items.push(shape.clone());
        true
    }

    /// `Contains(S)`.
    fn contains(&self, shape: &TopoShape) -> bool {
        self.find_index(shape) != 0
    }

    /// `FindKey(i)`: 1-based.
    fn find_key(&self, i: usize) -> Option<&TopoShape> {
        self.items.get(i.checked_sub(1)?)
    }

    /// `Extent()`.
    fn extent(&self) -> usize {
        self.items.len()
    }
}

/// `TopoDS_Iterator(S, false, false).Value()`: the stored children.
fn stored_children(s: &TopoShape) -> Vec<TopoShape> {
    s.tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone()
}

/// `TopExp::Vertices(W, Vfirst, Vlast)` (`TopExp.cxx:255-319`).
///
/// The wire is walked with `TopoDS_Iterator it(W)` (cumulated orientation), a
/// REVERSED edge swaps its two endpoints (`cxx:266-273`), and the endpoints go
/// through an `IsSame` toggle map so that a closed wire leaves the map empty
/// (`cxx:288-296`) and an open one leaves exactly two entries
/// (`cxx:297-318`).
pub fn topexp_vertices_of_wire(wire: &Wire) -> (Option<Vertex>, Option<Vertex>) {
    // `NCollection_Map<TopoDS_Shape>` keyed by `IsSame`; the stored value
    // keeps the orientation it was first inserted with.
    let mut vmap: Vec<Vertex> = Vec::new();
    let mut last_v2: Option<Vertex> = None;

    for child in cumulated_children(&wire.0) {
        if child.shape_type() != ShapeType::Edge {
            continue;
        }
        let edge = Edge(child);
        let mut v1;
        let mut v2;
        if edge.0.orientation() == Orientation::Reversed {
            // `cxx:268`: `TopExp::Vertices(E, V2, V1)`.
            let (a, b) = edge_vertices(&edge);
            v1 = b;
            v2 = a;
        } else {
            let (a, b) = edge_vertices(&edge);
            v1 = a;
            v2 = b;
        }
        // `cxx:275-284`: force the orientation role and toggle the map.
        let Some(mut key1) = v1.take() else {
            continue;
        };
        key1.0.set_orientation(Orientation::Forward);
        toggle_vertex(&mut vmap, &key1);
        if let Some(mut key2) = v2.take() {
            key2.0.set_orientation(Orientation::Reversed);
            toggle_vertex(&mut vmap, &key2);
            last_v2 = Some(key2);
        }
    }

    if vmap.is_empty() {
        // `cxx:288-296`: closed, both answers are the last `V2`.
        let Some(v2) = last_v2 else {
            return (None, None);
        };
        return (
            Some(Vertex(v2.0.oriented(Orientation::Forward))),
            Some(Vertex(v2.0.oriented(Orientation::Reversed))),
        );
    }
    if vmap.len() == 2 {
        // `cxx:297-318`: open, the FORWARD entry is the first vertex and the
        // REVERSED entry the last one.
        let first = vmap
            .iter()
            .find(|v| v.0.orientation() == Orientation::Forward)
            .cloned();
        let last = vmap
            .iter()
            .find(|v| v.0.orientation() == Orientation::Reversed)
            .cloned();
        return (first, last);
    }
    (None, None)
}

/// `NCollection_Map::Add` / `Remove` on an `IsSame` keyed sequence
/// (`TopExp.cxx:277-284`): a key already present is removed, a missing key is
/// appended.
fn toggle_vertex(map: &mut Vec<Vertex>, v: &Vertex) {
    let key = shape_map_key(&v.0);
    match map.iter().position(|s| shape_map_key(&s.0) == key) {
        Some(pos) => {
            map.remove(pos);
        }
        None => map.push(v.clone()),
    }
}

/// `ShapeAnalysis::FindBounds(shape, V1, V2)` (`ShapeAnalysis.cxx:73-103`).
pub fn find_bounds(shape: &TopoShape) -> (Option<Vertex>, Option<Vertex>) {
    match shape.shape_type() {
        ShapeType::Wire => topexp_vertices_of_wire(&Wire(shape.clone())),
        ShapeType::Edge => {
            let edge = Edge(shape.clone());
            (super::first_vertex(&edge), super::last_vertex(&edge))
        }
        ShapeType::Vertex => {
            let v = Vertex(shape.clone());
            (Some(v.clone()), Some(v))
        }
        _ => (None, None),
    }
}

/// `ShapeExtend_Explorer::SeqFromCompound(comp, expcomp)`
/// (`ShapeExtend_Explorer.cxx:70-86`).
fn seq_from_compound(comp: &TopoShape, expcomp: bool) -> Vec<TopoShape> {
    let mut list = Vec::new();
    if comp.shape_type() != ShapeType::Compound {
        list.push(comp.clone());
        return list;
    }
    fill_list(&mut list, comp, expcomp);
    list
}

/// `FillList` (`ShapeExtend_Explorer.cxx:48-68`).
fn fill_list(list: &mut Vec<TopoShape>, comp: &TopoShape, expcomp: bool) {
    for sub in stored_children(comp) {
        if sub.shape_type() != ShapeType::Compound || !expcomp {
            list.push(sub);
        } else {
            fill_list(list, &sub, expcomp);
        }
    }
}

/// `ShapeExtend_Explorer::CompoundFromSeq(seqval)`
/// (`ShapeExtend_Explorer.cxx:32-44`).
fn compound_from_seq(seq: &[TopoShape]) -> Compound {
    let builder = crate::builder::TopoBuilder::new();
    let mut comp = Compound::new();
    for s in seq {
        builder.add_compound(&mut comp, s);
    }
    comp
}

mod connect_wires;
mod free_bounds_class;
mod shape_analysis_shell;
mod split_wires;
mod wire_data;

pub use free_bounds_class::*;
pub use shape_analysis_shell::*;
pub(crate) use connect_wires::{connect_edges_to_wires, connect_wires_to_wires};
pub(crate) use split_wires::split_wires_impl;
pub(crate) use wire_data::{WireAnalyzer, WireData};

