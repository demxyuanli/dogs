//! `TopOpeBRepDS_DataStructure` — the shape-with-state maps and the rejected
//! shape maps.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_DataStructure.cxx:1287-1355`.
//!
//! OCCT has two maps of the same type and picks one by probing it, returning
//! a `static` dummy when the shape is in neither. A `&mut` cannot be handed
//! out for a static dummy, so [`DataStructure::change_map_of_shape_with_state`]
//! reports which map owns the shape instead (and `None` for the dummy case,
//! i.e. OCCT's `aFlag == false`), and the caller then picks the map through
//! [`DataStructure::change_map_of_shape_with_state_obj`] /
//! [`DataStructure::change_map_of_shape_with_state_tool`].

use crate::shape::TopoShape;

use super::super::shape_with_state::ShapeWithState;
use super::{DataStructure, IndexedShapeData, IndexedShapeSet};

/// Which of OCCT's two `myMapOfShapeWithState*` maps owns a shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeWithStateMap {
    /// `myMapOfShapeWithStateObj`
    Obj,
    /// `myMapOfShapeWithStateTool`
    Tool,
}

impl DataStructure {
    /// `ChangeMapOfShapeWithStateObj()` (`cxx:1287-1290`).
    pub fn change_map_of_shape_with_state_obj(&mut self) -> &mut IndexedShapeData<ShapeWithState> {
        &mut self.map_of_shape_with_state_obj
    }

    /// `ChangeMapOfShapeWithStateTool()` (`cxx:1292-1296`).
    pub fn change_map_of_shape_with_state_tool(
        &mut self,
    ) -> &mut IndexedShapeData<ShapeWithState> {
        &mut self.map_of_shape_with_state_tool
    }

    /// `myMapOfShapeWithStateObj` (read-only). OCCT reaches it inline through
    /// `ChangeMapOfShapeWithStateObj()`; the port also needs the immutable view
    /// where the caller already holds `&DataStructure` (e.g.
    /// `TopOpeBRepBuild_Builder1::PerformPieceOn2D` `cxx:1550-1553`).
    pub fn map_of_shape_with_state_obj(&self) -> &IndexedShapeData<ShapeWithState> {
        &self.map_of_shape_with_state_obj
    }

    /// `myMapOfShapeWithStateTool` (read-only), see
    /// [`DataStructure::map_of_shape_with_state_obj`].
    pub fn map_of_shape_with_state_tool(&self) -> &IndexedShapeData<ShapeWithState> {
        &self.map_of_shape_with_state_tool
    }

    /// `ChangeMapOfShapeWithState(aShape, aFlag)` (`cxx:1303-1321`).
    ///
    /// `Some((map, which))` is OCCT's `aFlag == true` plus the returned map;
    /// `None` is `aFlag == false` plus the `dummy` map.
    pub fn change_map_of_shape_with_state(
        &mut self,
        a_shape: &TopoShape,
    ) -> Option<(&mut IndexedShapeData<ShapeWithState>, ShapeWithStateMap)> {
        if self.map_of_shape_with_state_obj.contains(a_shape) {
            return Some((
                &mut self.map_of_shape_with_state_obj,
                ShapeWithStateMap::Obj,
            ));
        }
        if self.map_of_shape_with_state_tool.contains(a_shape) {
            return Some((
                &mut self.map_of_shape_with_state_tool,
                ShapeWithStateMap::Tool,
            ));
        }
        None
    }

    /// `GetShapeWithState(aShape)` (`cxx:1327-1341`). OCCT returns a `static`
    /// default-constructed dummy when the shape is in neither map: `None`.
    pub fn get_shape_with_state(&self, a_shape: &TopoShape) -> Option<&ShapeWithState> {
        if self.map_of_shape_with_state_obj.contains(a_shape) {
            return self.map_of_shape_with_state_obj.find_from_key(a_shape);
        }
        if self.map_of_shape_with_state_tool.contains(a_shape) {
            return self.map_of_shape_with_state_tool.find_from_key(a_shape);
        }
        None
    }

    /// `ChangeMapOfRejectedShapesObj()` (`cxx:1345-1348`).
    pub fn change_map_of_rejected_shapes_obj(&mut self) -> &mut IndexedShapeSet {
        &mut self.map_of_rejected_shapes_obj
    }

    /// `ChangeMapOfRejectedShapesTool()` (`cxx:1352-1355`).
    pub fn change_map_of_rejected_shapes_tool(&mut self) -> &mut IndexedShapeSet {
        &mut self.map_of_rejected_shapes_tool
    }
}
