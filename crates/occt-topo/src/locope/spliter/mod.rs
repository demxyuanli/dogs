//! `LocOpe_Spliter` (`LocOpe_Spliter.hxx:30-73`, `.cxx:63-577`, `.lxx:22-70`).
//!
//! Splits a shape by the wires/edges projected onto it by
//! [`WiresOnShape`](crate::locope::wires_on_shape::WiresOnShape), keeping the
//! faces on the "left" of those sections.
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | ctor / `Init` / `IsDone` / `Shape` / `ResultingShape` (`lxx:22-70`) | this file |
//! | `Perform` stages 1-6: vertex substitution, section bookkeeping, wire rebuild, descendant table, result shape (`cxx:63-266`) | [`perform`] |
//! | `Perform` stages 7-9: edge-on-edge substitution, degenerated-edge vertex merge, `DirectLeft` / `Left` (`cxx:268-577`) | [`finish`] |
//! | `Put` / `RebuildWires` / `Select` (file-local statics, `cxx:620-719`) | [`statics`] |
//!
//! Prerequisites already in the tree: `LocOpe_WiresOnShape` (T-493),
//! `LocOpe_SplitShape` (T-494), `LocOpe_BuildWires` (T-495) and
//! `BRepTools_Substitution` (`brep_tools_substitution.rs`).
//!
//! Deviations, all forced by the language or by an OCCT raise:
//! - `myRes` / `myShape` have no null state (`TopoShape` has no null shape), so
//!   they are `Option`s; `Perform` on a null shape returns the
//!   `Standard_NullObject` raise as `Err`.
//! - `InitEdgeIterator` is replayed by snapshotting the `PW` pairs up front
//!   (`myMapEF` cannot change while `Perform` runs), which is the same loop the
//!   reference writes.
//! - `NCollection_Map` iteration is insertion-ordered here instead of hash
//!   ordered, so `DirectLeft` / `Left` answers can differ in *order* only.

mod finish;
mod perform;
mod statics;

use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::TopoShape;

/// `LocOpe_Spliter` (`hxx:30-73`).
#[derive(Debug, Clone, Default)]
pub struct Spliter {
    /// `TopoDS_Shape myShape`
    my_shape: Option<TopoShape>,
    /// `bool myDone`
    my_done: bool,
    /// `TopoDS_Shape myRes`
    my_res: Option<TopoShape>,
    /// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>> myMap`
    my_map: ShapeValueMap<Vec<TopoShape>>,
    /// `NCollection_List<TopoDS_Shape> myDLeft`
    my_d_left: Vec<TopoShape>,
    /// `NCollection_List<TopoDS_Shape> myLeft`
    my_left: Vec<TopoShape>,
}

impl Spliter {
    /// `LocOpe_Spliter()` (`lxx:22-25`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_Spliter(S)` (`lxx:29-33`).
    pub fn with_shape(s: &TopoShape) -> Self {
        let mut this = Self::default();
        this.init(s);
        this
    }

    /// `Init(S)` (`lxx:37-41`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_shape = Some(s.clone());
        self.my_done = false;
    }

    /// `IsDone()` (`lxx:45-48`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` (`lxx:58-61`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `ResultingShape()` (`lxx:65-70`); `None` is the `StdFail_NotDone` raise.
    pub fn resulting_shape(&self) -> Option<&TopoShape> {
        if self.my_done {
            self.my_res.as_ref()
        } else {
            None
        }
    }

    /// `DirectLeft()` (`cxx:599-605`); `None` is the `StdFail_NotDone` raise.
    pub fn direct_left(&self) -> Option<&[TopoShape]> {
        if self.my_done {
            Some(&self.my_d_left)
        } else {
            None
        }
    }

    /// `Left()` (`cxx:610-616`); `None` is the `StdFail_NotDone` raise.
    pub fn left(&self) -> Option<&[TopoShape]> {
        if self.my_done {
            Some(&self.my_left)
        } else {
            None
        }
    }

    /// `DescendantShapes(S)` (`cxx:580-594`): the reference's `static empty`
    /// list for an unbound key, and the `StdFail_NotDone` raise for a
    /// not-yet-performed run, are both the empty answer here.
    pub fn descendant_shapes(&self, s: &TopoShape) -> Vec<TopoShape> {
        if !self.my_done {
            return Vec::new();
        }
        self.my_map.find(s).cloned().unwrap_or_default()
    }
}
