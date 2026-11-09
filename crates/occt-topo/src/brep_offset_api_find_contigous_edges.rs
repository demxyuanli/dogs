//! `BRepOffsetAPI_FindContigousEdges` - the seam/contiguity report of a
//! set of shapes, read out of a `BRepBuilderAPI_Sewing` run.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_FindContigousEdges.cxx:22-134`
//! and `.hxx:33-157` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepOffsetAPI_FindContigousEdges(Tol, Option)` | `cxx:25-30` | [`BRepOffsetApiFindContigousEdges::new`] |
//! | `Init(Tol, Option)` | `cxx:34-37` | [`BRepOffsetApiFindContigousEdges::init`] |
//! | `Add(Shape)` | `cxx:41-44` | [`BRepOffsetApiFindContigousEdges::add`] |
//! | `Perform()` | `cxx:48-51` | [`BRepOffsetApiFindContigousEdges::perform`] |
//! | `NbContigousEdges()` | `cxx:55-58` | [`BRepOffsetApiFindContigousEdges::nb_contigous_edges`] |
//! | `ContigousEdge(Index)` | `cxx:62-67` | [`BRepOffsetApiFindContigousEdges::contigous_edge`] |
//! | `ContigousEdgeCouple(Index)` | `cxx:72-77` | [`BRepOffsetApiFindContigousEdges::contigous_edge_couple`] |
//! | `SectionToBoundary(Section)` | `cxx:82-87` | [`BRepOffsetApiFindContigousEdges::section_to_boundary`] |
//! | `NbDegeneratedShapes()` | `cxx:92-95` | [`BRepOffsetApiFindContigousEdges::nb_degenerated_shapes`] |
//! | `DegeneratedShape(Index)` | `cxx:99-104` | [`BRepOffsetApiFindContigousEdges::degenerated_shape`] |
//! | `IsDegenerated(Shape)` | `cxx:108-111` | [`BRepOffsetApiFindContigousEdges::is_degenerated`] |
//! | `IsModified(Shape)` | `cxx:115-118` | [`BRepOffsetApiFindContigousEdges::is_modified`] |
//! | `Modified(Shape)` | `cxx:122-127` | [`BRepOffsetApiFindContigousEdges::modified`] |
//! | `Dump()` | `cxx:131-134` | [`BRepOffsetApiFindContigousEdges::dump`] |
//!
//! # Translation notes
//!
//! * The class is a pure forwarder (`hxx:154`: `mySewing` is the only member)
//!   and the port keeps it that way: [`Sewing`] is the whole state.
//! * `mySewing = new BRepBuilderAPI_Sewing` (`cxx:28`) is the OCCT default
//!   constructor, i.e. `Init(1.0e-06, true, true, true, false)`
//!   (`BRepBuilderAPI_Sewing.cxx:2090-2103`), which is [`Sewing::default`].
//! * `Init(tolerance, option)` forwards `Init(tolerance, option, false, true,
//!   false)` (`cxx:36`): the fourth argument keeps the `option4` default
//!   `false` of `BRepBuilderAPI_Sewing::Init` (`BRepBuilderAPI_Sewing.cxx:2104-2133`).
//! * The `Standard_OutOfRange` guards of `ContigousEdge` /
//!   `ContigousEdgeCouple` / `DegeneratedShape` are written `index < 0 || index
//!   > Nb*()` (`cxx:64`, `:74`, `:101`) - the lower bound is 0, not 1, and
//!   index 0 is then rejected by the `NCollection_IndexedMap` lookup itself.
//!   The port answers `None` for both cases instead of raising.
//! * `SectionToBoundary` raises `Standard_NoSuchObject` when
//!   `IsSectionBound(section)` is false (`cxx:84-86`); `Modified` does the
//!   same under `!IsModified()` (`cxx:124-126`). Both answer `None` here.
//! * `Perform` and `Dump` are `void` in the reference. `Sewing::perform`
//!   answers `Result` (it never fails today) and `Sewing::dump` builds the
//!   text the reference prints, so the port returns that text.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `NbEdges()` is declared (`hxx:89`) but has no definition anywhere in
//!   `BRepOffsetAPI_FindContigousEdges.cxx` (the file ends at `:134` with
//!   `Dump`), so the reference symbol itself is unimplemented; the port
//!   registers it as UNPORTED rather than inventing a count.

use crate::brep_builder_api_sewing::Sewing;
use crate::shape::{Edge, TopoShape};

/// `class BRepOffsetAPI_FindContigousEdges` (`hxx:44-157`).
pub struct BRepOffsetApiFindContigousEdges {
    /// `mySewing` (`hxx:154`): the only member.
    my_sewing: Sewing,
}

impl Default for BRepOffsetApiFindContigousEdges {
    fn default() -> Self {
        Self::new(1.0e-06, true)
    }
}

impl BRepOffsetApiFindContigousEdges {
    /// `BRepOffsetAPI_FindContigousEdges(const double tolerance = 1.0e-06,
    /// const bool option = true)` (`cxx:25-30`): `new BRepBuilderAPI_Sewing`
    /// then `Init(tolerance, option)`.
    pub fn new(tolerance: f64, option: bool) -> Self {
        // `cxx:28`: `mySewing = new BRepBuilderAPI_Sewing`.
        let mut this = Self {
            my_sewing: Sewing::default(),
        };
        // `cxx:29`.
        this.init(tolerance, option);
        this
    }

    /// `Init(const double tolerance, const bool option)` (`cxx:34-37`):
    /// `mySewing->Init(tolerance, option, false, true)`.
    pub fn init(&mut self, tolerance: f64, option: bool) {
        // `cxx:36`: the fifth argument is `BRepBuilderAPI_Sewing::Init`'s
        // `option4 = false` default.
        self.my_sewing.init(tolerance, option, false, true, false);
    }

    /// `Add(const TopoDS_Shape& aShape)` (`cxx:41-44`).
    pub fn add(&mut self, shape: &TopoShape) {
        self.my_sewing.add(shape);
    }

    /// `Perform()` (`cxx:48-51`).
    pub fn perform(&mut self) {
        let _ = self.my_sewing.perform();
    }

    /// `NbContigousEdges()` (`cxx:55-58`).
    pub fn nb_contigous_edges(&self) -> usize {
        self.my_sewing.nb_contigous_edges()
    }

    /// `ContigousEdge(const int index)` (`cxx:62-67`). The reference's guard is
    /// `index < 0 || index > NbContigousEdges()`; a miss answers `None`.
    pub fn contigous_edge(&self, index: usize) -> Option<Edge> {
        self.my_sewing.contigous_edge(index)
    }

    /// `ContigousEdgeCouple(const int index)` (`cxx:72-77`).
    pub fn contigous_edge_couple(&self, index: usize) -> Option<&Vec<TopoShape>> {
        self.my_sewing.contigous_edge_couple(index)
    }

    /// `SectionToBoundary(const TopoDS_Edge& section)` (`cxx:82-87`):
    /// `Standard_NoSuchObject` when the section is not bound.
    pub fn section_to_boundary(&self, section: &TopoShape) -> Option<Edge> {
        self.my_sewing.section_to_boundary(section)
    }

    /// `NbDegeneratedShapes()` (`cxx:92-95`).
    pub fn nb_degenerated_shapes(&self) -> usize {
        self.my_sewing.nb_degenerated_shapes()
    }

    /// `DegeneratedShape(const int index)` (`cxx:99-104`).
    pub fn degenerated_shape(&self, index: usize) -> Option<TopoShape> {
        self.my_sewing.degenerated_shape(index)
    }

    /// `IsDegenerated(const TopoDS_Shape& aShape)` (`cxx:108-111`).
    /// `BRepBuilderAPI_Sewing::IsDegenerated` is non-const in the port
    /// (`cxx:2465-2494` fills `myContigSecBound` while it walks), so the
    /// forwarder needs `&mut self`, like every other `IsDegenerated` caller.
    pub fn is_degenerated(&mut self, shape: &TopoShape) -> bool {
        self.my_sewing.is_degenerated(shape)
    }

    /// `IsModified(const TopoDS_Shape& aShape)` (`cxx:115-118`).
    pub fn is_modified(&self, shape: &TopoShape) -> bool {
        self.my_sewing.is_modified(shape)
    }

    /// `Modified(const TopoDS_Shape& aShape)` (`cxx:122-127`):
    /// `Standard_NoSuchObject` when `IsModified(aShape)` is false.
    pub fn modified<'a>(&'a self, shape: &'a TopoShape) -> Option<&'a TopoShape> {
        if !self.is_modified(shape) {
            return None;
        }
        Some(self.my_sewing.modified(shape))
    }

    /// `Dump()` (`cxx:131-134`): the reference prints the sewing report; the
    /// port answers it.
    pub fn dump(&mut self) -> String {
        self.my_sewing.dump()
    }

    /// The `BRepBuilderAPI_Sewing` object the reference hides in `mySewing`
    /// (`hxx:154`).
    pub fn sewing(&self) -> &Sewing {
        &self.my_sewing
    }
}
