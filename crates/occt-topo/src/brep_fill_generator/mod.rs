//! `BRepFill_Generator` — the shell ruled through a sequence of section wires.
//!
//! Source: `BRepFill_Generator.cxx:598-1238`, `.hxx:31-72`,
//! `BRepFill_ThruSectionErrorStatus.hxx:18-30` (`TKBool/BRepFill`).
//!
//! Each pair of consecutive wires is walked in parallel by two
//! `BRepTools_WireExplorer`s (`cxx:629`, `cxx:667-671`); every edge pair gets a
//! ruled patch from `GeomFill_Generator` (`cxx:851-856`) or, for the six
//! geometric particular cases of `DetectKPart` (`cxx:62-355`), an analytic
//! quadric from `CreateKPart` (`cxx:362-569`). The four boundary edges of the
//! patch get their pcurves written as iso lines of the patch surface
//! (`cxx:1008-1077`), and the two connector edges are shared between adjacent
//! patches through the per-wire-pair vertex map (`cxx:884-972`).
//!
//! UNPORTED: the `myMutableInput == false` path (`cxx:974-1006`, `:1136-1157`)
//! needs `BRepTools_ReShape` (`Replace` / `Value` / `Apply`) and
//! `BRepLib::SameParameter(theShape, theReshaper)`, neither of which is ported.
//! `set_mutable_input(false)` is therefore accepted but `perform` keeps the
//! caller's edges in place (it attaches the pcurves to them) instead of
//! pre-copying them through `CreateNewEdge` (`cxx:573-595`).

use crate::shape::{Face, Shell, TopoShape, Wire};
use crate::topo_tools::ShapeMapKey;

pub mod detect_kpart;
pub mod perform;

pub use detect_kpart::{create_kpart, detect_kpart};

/// `BRepFill_ThruSectionErrorStatus` (`BRepFill_ThruSectionErrorStatus.hxx:18-30`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ThruSectionErrorStatus {
    /// `BRepFill_ThruSectionErrorStatus_Done`.
    Done,
    /// `BRepFill_ThruSectionErrorStatus_NotDone` — the constructor state
    /// (`cxx:600`).
    #[default]
    NotDone,
    /// `BRepFill_ThruSectionErrorStatus_NotSameTopology`.
    NotSameTopology,
    /// `BRepFill_ThruSectionErrorStatus_ProfilesInconsistent`.
    ProfilesInconsistent,
    /// `BRepFill_ThruSectionErrorStatus_WrongUsage`.
    WrongUsage,
    /// `BRepFill_ThruSectionErrorStatus_Null3DCurve`.
    Null3DCurve,
    /// `BRepFill_ThruSectionErrorStatus_Failed`.
    Failed,
}

/// `BRepFill_Generator` (`BRepFill_Generator.hxx:31-72`).
#[derive(Debug, Clone, Default)]
pub struct Generator {
    /// `myWires` — the profile wires in section order.
    pub(crate) wires: Vec<Wire>,
    /// `myShell` — null until `Perform` runs (`cxx:624`).
    pub(crate) shell: Option<Shell>,
    /// `myMap` — input edge to the faces ruled from it (`cxx:1113-1120`).
    pub(crate) map: Vec<(ShapeMapKey, Vec<Face>)>,
    /// `myOldNewShapes` — input edge to its copy (`cxx:696`, `:712`).
    pub(crate) old_new: Vec<(ShapeMapKey, TopoShape)>,
    /// `myMutableInput` — default `true` (`cxx:599`).
    pub(crate) mutable_input: bool,
    /// `myStatus` — default `NotDone` (`cxx:600`).
    pub(crate) status: ThruSectionErrorStatus,
}

impl Generator {
    /// `BRepFill_Generator()` (`cxx:598-602`): `myMutableInput = true`,
    /// `myStatus = NotDone`.
    pub fn new() -> Self {
        Self {
            wires: Vec::new(),
            shell: None,
            map: Vec::new(),
            old_new: Vec::new(),
            mutable_input: true,
            status: ThruSectionErrorStatus::NotDone,
        }
    }

    /// `AddWire(Wire)` (`cxx:606-609`).
    pub fn add_wire(&mut self, wire: &Wire) {
        self.wires.push(wire.clone());
    }

    /// `Perform()` (`cxx:613-1186`). The status carries the outcome, exactly as
    /// `GetStatus()` does in OCCT.
    pub fn perform(&mut self) {
        perform::perform(self);
    }

    /// `Shell()` (`.lxx:17-22`): null before `Perform`.
    pub fn shell(&self) -> Option<&Shell> {
        self.shell.as_ref()
    }

    /// `Generated()` (`cxx:1206-1210`).
    pub fn generated(&self) -> &[(ShapeMapKey, Vec<Face>)] {
        &self.map
    }

    /// `GeneratedShapes(SSection)` (`cxx:1190-1201`): an unbound section
    /// answers the static empty list.
    pub fn generated_shapes(&self, section: &TopoShape) -> Vec<TopoShape> {
        let key = crate::topo_tools::shape_map_key(section);
        self.map
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, faces)| faces.iter().map(|f| f.0.clone()).collect())
            .unwrap_or_default()
    }

    /// `ResultShape(theShape)` (`cxx:1214-1225`): the copy bound in
    /// `myOldNewShapes`, following the `myReshaper` chain (UNPORTED, see the
    /// module docs — the chain collapses to the `myOldNewShapes` lookup here).
    pub fn result_shape(&self, the_shape: &TopoShape) -> TopoShape {
        let key = crate::topo_tools::shape_map_key(the_shape);
        self.old_new
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, s)| s.clone())
            .unwrap_or_else(|| the_shape.clone())
    }

    /// `SetMutableInput(theIsMutableInput)` (`cxx:1229-1232`).
    pub fn set_mutable_input(&mut self, is_mutable: bool) {
        self.mutable_input = is_mutable;
    }

    /// `IsMutableInput()` (`cxx:1236-1239`).
    pub fn is_mutable_input(&self) -> bool {
        self.mutable_input
    }

    /// `GetStatus()` (`.hxx:63`).
    pub fn status(&self) -> ThruSectionErrorStatus {
        self.status
    }

    /// `myOldNewShapes.Bind(anOld, aNew)` (`cxx:696`, `:712`).
    pub(crate) fn bind_old_new(&mut self, old: &TopoShape, new: TopoShape) {
        let key = crate::topo_tools::shape_map_key(old);
        if let Some(slot) = self.old_new.iter_mut().find(|(k, _)| *k == key) {
            slot.1 = new;
        } else {
            self.old_new.push((key, new));
        }
    }

    /// `myMap(aREd).Append(Face)` (`cxx:1120`).
    pub(crate) fn append_generated(&mut self, edge: &TopoShape, face: &Face) {
        let key = crate::topo_tools::shape_map_key(edge);
        if let Some(slot) = self.map.iter_mut().find(|(k, _)| *k == key) {
            slot.1.push(face.clone());
        } else {
            self.map.push((key, vec![face.clone()]));
        }
    }
}
