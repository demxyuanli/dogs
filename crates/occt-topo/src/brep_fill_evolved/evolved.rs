//! The `BRepFill_Evolved` object itself (`BRepFill_Evolved.hxx:38-200`): the
//! state the sweep accumulates, and the accessors over it
//! (`BRepFill_Evolved.cxx:186-218`, `:1783-1800`, `:1959-1966`, `:2057-2078`).
//!
//! The performing half of the class lives in the sibling modules:
//! [`super::make_pipe`] (`MakePipe`, `cxx:2266-2327`), [`super::make_revol`]
//! (`MakeRevol`, `cxx:2329-2400`), and the `Perform` chain that drives them.
//!
//! OCCT's null `TopoDS_*` members become `Option<...>`, the convention the
//! rest of the port uses for a null `TShape` ([`TopoShape::is_null`] is always
//! false here).

use occt_core::kernel::geomabs::JoinType;

use crate::shape::{Face, TopoShape, Wire};

use super::shape_map::ShapeMap;

/// `BRepFill_Evolved` (`BRepFill_Evolved.hxx:38-200`).
#[derive(Clone)]
pub struct BRepFillEvolved {
    /// `mySpine` (`hxx:186`).
    pub(crate) my_spine: Option<Face>,
    /// `myProfile` (`hxx:187`).
    pub(crate) my_profile: Option<Wire>,
    /// `myShape` (`hxx:188`).
    pub(crate) my_shape: Option<TopoShape>,
    /// `myIsDone` (`hxx:189`).
    pub(crate) my_is_done: bool,
    /// `mySpineType` (`hxx:190`): `true` for the default constructor, then
    /// `false` for a wire spine and `true` for a face spine (`cxx:317`,
    /// `cxx:330`).
    pub(crate) my_spine_type: bool,
    /// `myJoinType` (`hxx:191`).
    pub(crate) my_join_type: JoinType,
    /// `myMap` (`hxx:192-197`):
    /// `DataMap<Shape, DataMap<Shape, List<Shape>>>`.
    pub(crate) my_map: ShapeMap<ShapeMap<Vec<TopoShape>>>,
    /// `myTop` (`hxx:198`).
    pub(crate) my_top: Option<TopoShape>,
    /// `myBottom` (`hxx:199`).
    pub(crate) my_bottom: Option<TopoShape>,
}

impl Default for BRepFillEvolved {
    /// `BRepFill_Evolved()` (`cxx:186-190`): `myIsDone = false`,
    /// `mySpineType = true`; every other member keeps its `TopoDS_*` /
    /// `GeomAbs_*` default.
    fn default() -> Self {
        Self {
            my_spine: None,
            my_profile: None,
            my_shape: None,
            my_is_done: false,
            my_spine_type: true,
            my_join_type: JoinType::Arc,
            my_map: ShapeMap::new(),
            my_top: None,
            my_bottom: None,
        }
    }
}

impl BRepFillEvolved {
    /// `BRepFill_Evolved()` (`cxx:186-190`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IsDone()` (`cxx:2057-2062`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `Shape()` (`cxx:2064-2076`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `JoinType()` (`BRepFill_Evolved.hxx:96`): `myJoinType` as set by the
    /// last `PrivatePerform`.
    pub fn join_type(&self) -> JoinType {
        self.my_join_type
    }

    /// `Top()` (`cxx:1783-1788`).
    pub fn top(&self) -> Option<&TopoShape> {
        self.my_top.as_ref()
    }

    /// `Bottom()` (`cxx:1790-1795`).
    pub fn bottom(&self) -> Option<&TopoShape> {
        self.my_bottom.as_ref()
    }

    /// `ChangeShape()` (`cxx:1959-1964`).
    pub fn change_shape(&mut self) -> Option<&mut TopoShape> {
        self.my_shape.as_mut()
    }

    /// `Generated()` (`BRepFill_Evolved.hxx:130-134`): the mutable `myMap`.
    pub fn generated(&mut self) -> &mut ShapeMap<ShapeMap<Vec<TopoShape>>> {
        &mut self.my_map
    }

    /// `GeneratedShapes(SpineShape, ProfShape)` (`cxx:1797-1810`): the list of
    /// shapes generated from a spine sub-shape by a profile sub-shape. The
    /// reference's static empty list is an empty `Vec` here.
    pub fn generated_shapes(&self, spine_shape: &TopoShape, prof_shape: &TopoShape) -> Vec<TopoShape> {
        self.my_map
            .list(spine_shape, prof_shape)
            .cloned()
            .unwrap_or_default()
    }
}
