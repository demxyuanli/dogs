//! `LocOpe_SplitDrafts` (`ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitDrafts`).
//!
//! The class splits one face of a shape along a wire lying on it and puts a
//! draft angle on each side of the wire. `Perform` builds two planes (`NewSg`
//! for the "left" side, `NewSd` for the "right") that envelope the drafted
//! walls, splits the wire with the line where those planes meet, splits the
//! face along the wire with `LocOpe_Spliter`, then rebuilds every piece of the
//! face on the corresponding draft plane and substitutes the new edges back
//! into the rest of the shape with `BRepTools_Substitution`.
//!
//! Module layout (each mirrors one block of the reference):
//!
//! * [`state`] - the locals `Perform` carries between stages.
//! * [`perform`] - `Perform` (`cxx:105-1436`), the stage driver.
//! * [`wire_split`] - the `theW` split by the pipe surface (`cxx:180-351`).
//! * [`split_face`] - `LocOpe_Spliter` pass and wire bookkeeping (`cxx:353-403`).
//! * [`generated`] - the `MapEV` "generated face per edge" pass (`cxx:405-470`).
//! * [`vertices`] - the per-vertex `Choice` 1 / 2 / 3 pass (`cxx:471-866`).
//! * [`wire_edges`] - `MapW` / `MapSg` / `MapSd` per wire edge (`cxx:868-1100`).
//! * [`draft_faces`] - the drafted faces and `MapSonS` (`cxx:1102-1358`).
//! * [`descendants`] - `myMap` fix-ups and `LocOpe_BuildShape` (`cxx:1360-1435`).
//! * [`statics`] - `NewPlane` / `MakeFace` / `Contains` / `NewEdge`
//!   (`cxx:1461-1976`).

mod descendants;
mod draft_faces;
mod generated;
mod perform;
mod split_face;
mod state;
mod statics;
mod vertices;
mod wire_edges;
mod wire_split;

use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::TopoShape;

pub use perform::PerformArgs;

/// `LocOpe_SplitDrafts` (`hxx:37-111`).
#[derive(Debug, Clone, Default)]
pub struct SplitDrafts {
    /// `myShape` (`hxx:108`). `None` is the reference's null `TopoDS_Shape`.
    my_shape: Option<TopoShape>,
    /// `myResult` (`hxx:109`). `None` is the null result of a failed run.
    my_result: Option<TopoShape>,
    /// `myMap` (`hxx:110`):
    /// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`.
    my_map: ShapeValueMap<Vec<TopoShape>>,
}

impl SplitDrafts {
    /// Empty constructor (`hxx:43`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_SplitDrafts(const TopoDS_Shape& S)` (`hxx:48-51`).
    pub fn with_shape(s: &TopoShape) -> Self {
        let mut this = Self::new();
        this.init(s);
        this
    }

    /// `Init(const TopoDS_Shape& S)` (`cxx:85-91`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_shape = Some(s.clone());
        self.my_result = None;
        self.my_map.clear();
    }

    /// `IsDone()` (`hxx:95`).
    pub fn is_done(&self) -> bool {
        self.my_result.is_some()
    }

    /// `OriginalShape()` (`hxx:97`).
    pub fn original_shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `Shape()` (`cxx:1438-1445`); `None` is the `StdFail_NotDone` raise.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_result.as_ref()
    }

    /// `ShapesFromShape(S)` (`cxx:1449-1457`); `None` is the `StdFail_NotDone`
    /// raise that precedes the `myMap(S)` lookup.
    pub fn shapes_from_shape(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        if self.my_result.is_none() {
            return None;
        }
        self.my_map.find(s).map(|v| v.as_slice())
    }
}
