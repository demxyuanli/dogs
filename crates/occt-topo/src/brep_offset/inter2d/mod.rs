//! `BRepOffset_Inter2d` (`BRepOffset_Inter2d.hxx` / `BRepOffset_Inter2d.cxx`):
//! the 2D (pcurve) intersections of the edges lying on one face, and the
//! intersection vertices they produce.
//!
//! Ported here so far (control flow taken from
//! `src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_Inter2d.cxx`):
//!
//! - file-local `CommonVertex` (`:83-97`) - `orient_edges.rs`
//! - file-local `DefineClosedness` (`:101-124`) - `orient_edges.rs`
//! - file-local `GetEdgesOrientedInFace` (`:128-266`) - `orient_edges.rs`
//! - file-local `Store` (`:276-390`) and its two-edge overload (`:396-412`)
//!   - `store_vertices.rs`
//! - file-local `EdgeInter` (`:416-699`) - `edge_inter.rs`
//! - file-local `RefEdgeInter` (`:703-1070`) - `ref_edge_inter.rs`
//! - file-local `UpdateVertex` (`:1680-1721`) - `vertex_update.rs`
//! - `Compute` (`:1725-1827`) - `compute.rs`
//! - `ExtentEdge` (`:1193-1676`) with its helpers `ExtendPCurve`
//!   (`:1100-1187`) and `evaluateMaxSegment` (`:1078-1097`) - `extent_edge.rs`
//! - `ConnexIntByInt` (`:1830-2094`) and `ConnexIntByIntInVert` (`:2100-2308`)
//!   - `connex_int_by_int.rs`
//! - file-local `MakeChain` (`:2311-2332`) and `FuseVertices` (`:2335-2393`)
//!   - `fuse_vertices.rs`
//!
//! The whole `.cxx` of this class is ported; nothing is UNPORTED.
//!
//! `theDMVV` ([`Dmvv`]) is the `NCollection_IndexedDataMap<TopoDS_Shape,
//! NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>` of
//! coincident-vertex chains `Store` fills and `FuseVertices` consumes: keyed
//! by [`ShapeMapKey`] (`IsSame`), with order-preserving lists and an
//! order-preserving key walk (`FindKey(i)`).

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use occt_geom2d::Curve2d;

use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{shape_map_key, ShapeMapKey};

mod compute;
mod connex_int_by_int;
mod edge_inter;
mod extent_edge;
mod fuse_vertices;
mod orient_edges;
mod ref_edge_inter;
mod store_vertices;
mod vertex_update;

pub use connex_int_by_int::{connex_int_by_int, connex_int_by_int_in_vert, ConnexByIntData};
pub use extent_edge::extent_edge;
pub use fuse_vertices::fuse_vertices;

/// `theDMVV` (`NCollection_IndexedDataMap<TopoDS_Shape,
/// NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>`) of every `Store`
/// / `Compute` / `ConnexIntByInt` / `FuseVertices` signature.
///
/// The map keys are vertices, `IsSame` identity (TShape plus location, see
/// [`ShapeMapKey`]), and `Add(key, value)` binds `key` - appending it to the
/// key order on first use - then appends `value` to that key's list. The key
/// order is what `FindKey(i)` walks in `FuseVertices` (`cxx:2345-2391`).
#[derive(Default)]
pub struct Dmvv {
    /// `myData1` of the `IndexedDataMap`: the keys in insertion order.
    order: Vec<TopoShape>,
    /// `myIndex`: key membership.
    keys: HashSet<ShapeMapKey>,
    /// `myData2`: key -> the list appended under it.
    lists: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl Dmvv {
    pub fn new() -> Self {
        Self::default()
    }

    /// `theDMVV.Add(theKey, theValue)`: bind the key if absent, then append.
    pub fn add(&mut self, key: &TopoShape, value: &TopoShape) {
        let k = shape_map_key(key);
        if self.keys.insert(k.clone()) {
            self.order.push(key.clone());
        }
        self.lists.entry(k).or_default().push(value.clone());
    }

    /// `theDMVV.Add(theKey, theList)`.
    pub fn add_list(&mut self, key: &TopoShape, values: &[TopoShape]) {
        for value in values {
            self.add(key, value);
        }
    }

    /// `NCollection_IndexedDataMap::Add(theKey, theItem)` where `theItem` is a
    /// list: the key keeps its position on the first use, and its list is
    /// **replaced** by `values` (`NCollection_IndexedDataMap.hxx`: `Add`
    /// rewrites `myData2` for an existing key). Used by
    /// `BRepOffset_BuildOffsetFaces` (`_1.cxx:1075`, `:1080`) where the
    /// reference does
    /// `myOFImages(myOFImages.Add(aF, NCollection_List<...>())).Append(...)`.
    pub fn bind_list(&mut self, key: &TopoShape, values: &[TopoShape]) {
        let k = shape_map_key(key);
        if self.keys.insert(k.clone()) {
            self.order.push(key.clone());
        }
        self.lists.insert(k, values.to_vec());
    }

    /// `theDMVV.ChangeSeek(theKey)`: the mutable list under `key`.
    pub fn change_seek(&mut self, key: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        self.lists.get_mut(&shape_map_key(key))
    }

    /// `theDMVV.Seek(theKey)`: the list under `key`, `None` for an unbound
    /// key (the null pointer the reference tests).
    pub fn seek(&self, key: &TopoShape) -> Option<&[TopoShape]> {
        self.lists.get(&shape_map_key(key)).map(|v| v.as_slice())
    }

    /// `theDMVV.FindKey(i)` (`1 <= i <= Extent()`): the keys in insertion
    /// order.
    pub fn keys(&self) -> &[TopoShape] {
        &self.order
    }

    /// `theDMVV.Extent()`.
    pub fn extent(&self) -> usize {
        self.order.len()
    }

    /// `NCollection_DataMap::IsBound` / `IndexedDataMap::Contains`.
    pub fn contains(&self, key: &TopoShape) -> bool {
        self.keys.contains(&shape_map_key(key))
    }

    /// `NCollection_DataMap::UnBind(theKey)` (`FilterSplits`,
    /// `_1.cxx:9038`; `TrimNewIntersectionEdges`, `:8289`/`:8371`): drops the
    /// key and its list. The reference's `IsBound` then answers false, so the
    /// key must leave the order too.
    pub fn unbind(&mut self, key: &TopoShape) -> bool {
        let k = shape_map_key(key);
        if !self.keys.remove(&k) {
            return false;
        }
        self.lists.remove(&k);
        self.order.retain(|s| shape_map_key(s) != k);
        true
    }
}

/// `BRepOffset_Inter2d` (the class has no data members; every method is
/// static).
pub struct BRepOffsetInter2d;

/// `BRepTool::Curve(E, F)` evaluated at `u`, i.e. what
/// `BRepAdaptor_Curve(E, F)` returns: `BRep_Tool::CurveOnSurface(E, F, f, l)`
/// maps the pcurve onto the face surface (`Geom2d_Curve::To3d(S, E.Location(),
/// First, Last)`) and `BRepAdaptor_Curve` applies the edge location.
pub(crate) fn pcurve_value_on_face(
    basurf: &crate::topopebrep::tool::adaptor::SurfaceAdaptor,
    pc: &dyn Curve2d,
    u: f64,
    e: &Edge,
) -> Option<occt_core::gp::GpPnt> {
    let p2d = pc.d0(u);
    let p = basurf.value(p2d.x(), p2d.y())?;
    let loc = e.0.location();
    if loc.is_identity() {
        Some(p)
    } else {
        Some(p.transformed(&loc.transformation()))
    }
}


