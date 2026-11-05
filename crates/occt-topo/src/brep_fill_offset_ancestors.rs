//! `BRepFill_OffsetAncestors` — which spine sub-shape each offset shape came
//! from.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_OffsetAncestors.cxx:26-94`
//! and `.hxx:25-66` (OCCT 8.0.0).
//!
//! `Perform` (`cxx:42-67`) walks every edge and then every vertex of the
//! offset wire's spine and binds each shape `BRepFill_OffsetWire` generated
//! from it to that sub-shape, so `Ancestor` answers the generating spine shape
//! of an offset edge. `IsDone` (`cxx:73-76`) reports whether `Perform` ran and
//! `Ancestor` (`cxx:86-94`) raises `StdFail_NotDone` before that.
//!
//! Deviation: `myMap` is an `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>`
//! with the default `TopTools_ShapeMapHasher` (keyed by `IsSame`), i.e. a
//! `HashMap<ShapeMapKey, TopoShape>` here; `Bind` replaces an existing binding
//! as `NCollection_DataMap::Bind` does.

use std::collections::HashMap;

use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::shape::TopoShape;
use crate::topo_tools::{edges_of, shape_map_key, vertices_of, ShapeMapKey};

/// `class BRepFill_OffsetAncestors` (`BRepFill_OffsetAncestors.hxx:25-66`).
#[derive(Clone, Default)]
pub struct BRepFillOffsetAncestors {
    /// `myMap` (`.hxx:62`): offset shape -> generating spine sub-shape.
    my_map: HashMap<ShapeMapKey, TopoShape>,
    /// `myIsPerform` (`.hxx:63`).
    my_is_perform: bool,
}

impl BRepFillOffsetAncestors {
    /// `BRepFill_OffsetAncestors()` (`cxx:26-30`).
    pub fn new() -> Self {
        Self {
            my_map: HashMap::new(),
            my_is_perform: false,
        }
    }

    /// `BRepFill_OffsetAncestors(Paral)` (`cxx:33-37`).
    pub fn with_offset_wire(paral: &mut BRepFillOffsetWire) -> Self {
        let mut ancestors = Self::new();
        ancestors.perform(paral);
        ancestors
    }

    /// `Perform(Paral)` (`cxx:42-67`).
    pub fn perform(&mut self, paral: &mut BRepFillOffsetWire) {
        let spine = paral
            .spine()
            .expect("BRepFill_OffsetAncestors::Perform: null spine")
            .0
            .clone();

        // `cxx:47-56`: iterate the edges of the spine.
        for edge in edges_of(&spine) {
            for generated in paral.generated_shapes(&edge.0) {
                self.my_map.insert(shape_map_key(&generated), edge.0.clone());
            }
        }

        // `cxx:58-65`: then the vertices.
        for vertex in vertices_of(&spine) {
            for generated in paral.generated_shapes(&vertex.0) {
                self.my_map.insert(shape_map_key(&generated), vertex.0.clone());
            }
        }

        self.my_is_perform = true;
    }

    /// `IsDone()` (`cxx:73-76`).
    pub fn is_done(&self) -> bool {
        self.my_is_perform
    }

    /// `HasAncestor(S1)` (`cxx:78-81`): `myMap.IsBound(S1)`.
    pub fn has_ancestor(&self, s1: &TopoShape) -> bool {
        self.my_map.contains_key(&shape_map_key(s1))
    }

    /// `Ancestor(S1)` (`cxx:87-94`); OCCT raises `StdFail_NotDone` when
    /// `Perform` has not run and `Standard_NoSuchObject` for an unbound shape,
    /// both of which panic here.
    pub fn ancestor(&self, s1: &TopoShape) -> &TopoShape {
        assert!(
            self.my_is_perform,
            "BRepFill_OffsetAncestors::Ancestor() - Perform() should be called before accessing results"
        );
        self.my_map
            .get(&shape_map_key(s1))
            .unwrap_or_else(|| panic!("BRepFill_OffsetAncestors::Ancestor: shape is not bound"))
    }
}
