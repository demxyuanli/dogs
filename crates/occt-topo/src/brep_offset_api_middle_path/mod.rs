//! `BRepOffsetAPI_MiddlePath` - the middle path of a pipe-like shape.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MiddlePath.cxx:296-1148`
//! and `.hxx:31-57` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepOffsetAPI_MiddlePath(Shape, StartShape, EndShape)` | `cxx:296-331` | [`BRepOffsetApiMiddlePath::new`] |
//! | `Build(Range)` | `cxx:336-1148` | [`BRepOffsetApiMiddlePath::build`] |
//! | - wire walks / `myPaths` init | `cxx:339-421` | `paths::collect` |
//! | - `myPaths` filling | `cxx:424-497` | `paths::fill_paths` |
//! | - sections | `cxx:522-843` | `sections::build_sections` |
//! | - middle path | `cxx:845-1147` | `middle_edges::build_middle_path` |
//!
//! # Translation notes
//!
//! * The constructor unifies the input shape before replacing every start/end
//!   wire edge by its `History()->Modified` images (`cxx:266-303`).
//! * `myClosedSection = myStartWire.Closed()` (`cxx:329`) and
//!   `myClosedRing = myStartWire.IsSame(myEndWire)` (`cxx:330`) are stored
//!   verbatim; `myClosedRing` is never read by OCCT.
//! * `StartShape`/`EndShape` that are neither a face nor a wire are a
//!   `TopoDS::Wire` cast (`cxx:313`, `:322`), i.e. `Standard_TypeMismatch` in
//!   OCCT; the port answers `Err`.
//! * `Build` is `void` in the reference but has two `return` arms
//!   (`cxx:416-421`) that leave the object `!IsDone()`, and the `TopoDS::Edge`
//!   casts inside it raise on a non-edge. The port returns `Result` and keeps
//!   both silent-return arms as `Ok` without `Done()`.
//!
//! UNPORTED:
//!
//! * `Message_ProgressRange` (`cxx:336`) - there is no progress port.
//! * The `std::cout << "Interpolation failed"` (`cxx:1097`) and the two
//!   commented-out "Temporary" debug blocks (`cxx:500-519`, `cxx:1117-1145`).

mod middle_edges;
mod paths;
mod sections;
mod statics;

use std::collections::{HashMap, HashSet};

use crate::brep_offset_api_middle_path::statics::get_unified_wire;
use crate::shape::{Face, TopoShape, Wire};
use crate::shhealing::UnifySameDomain;
use crate::topo_tools::{is_same, AncestorMap, ShapeMapKey};
/// The locals of `Build` that more than one phase reads. Everything else stays
/// inside the phase that uses it.
#[derive(Default)]
pub(super) struct Ctx {
    /// `EndVertices` (`cxx:340`).
    pub end_vertices: HashSet<ShapeMapKey>,
    /// `EndEdges` (`cxx:341`).
    pub end_edges: HashSet<ShapeMapKey>,
    /// `SectionsEdges` (`cxx:342`), 1-based in the reference.
    pub sections_edges: Vec<Vec<Option<TopoShape>>>,
    /// `EdgeSeq` (`cxx:345`), 1-based in the reference.
    pub edge_seq: Vec<Option<TopoShape>>,
    /// `VEmap` (`cxx:379-381`).
    pub vemap: AncestorMap,
    /// `EFmap` (`cxx:379-382`); `EFmap.Add(NewEdge, ListOneFace)` writes it.
    pub efmap: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `CurVertices` (`cxx:384`).
    pub cur_vertices: HashSet<ShapeMapKey>,
}

/// `class BRepOffsetAPI_MiddlePath` (`hxx:31-57`).
pub struct BRepOffsetApiMiddlePath {
    /// `myInitialShape` (`hxx:46`).
    my_initial_shape: TopoShape,
    /// `myStartWire` (`hxx:47`).
    my_start_wire: Wire,
    /// `myEndWire` (`hxx:48`).
    my_end_wire: Wire,
    /// `myClosedSection` (`hxx:49`).
    my_closed_section: bool,
    /// `myClosedRing` (`hxx:50`): written, never read (as in OCCT).
    #[allow(dead_code)]
    my_closed_ring: bool,
    /// `myStartWireEdges` (`hxx:51`).
    my_start_wire_edges: HashSet<ShapeMapKey>,
    /// `myEndWireEdges` (`hxx:52`).
    my_end_wire_edges: HashSet<ShapeMapKey>,
    /// `myPaths` (`hxx:53`), 1-based in the reference.
    my_paths: Vec<Vec<TopoShape>>,
    /// `BRepBuilderAPI_MakeShape::myShape`.
    my_shape: Option<TopoShape>,
    /// `BRepBuilderAPI_MakeShape::myDone`.
    my_done_flag: bool,
}

impl BRepOffsetApiMiddlePath {
    /// `BRepOffsetAPI_MiddlePath(const TopoDS_Shape& aShape, const
    /// TopoDS_Shape& StartShape, const TopoDS_Shape& EndShape)` (`cxx:296-331`).
    pub fn new(
        a_shape: &TopoShape,
        start_shape: &TopoShape,
        end_shape: &TopoShape,
    ) -> Result<Self, String> {
        // `cxx:301-303`.
        let mut unifier = UnifySameDomain::new(a_shape.clone());
        unifier.build()?;
        let my_initial_shape = unifier.shape().clone();

        // `cxx:305-325`.
        let a_start_wire = wire_of(start_shape)?;
        let an_end_wire = wire_of(end_shape)?;

        // `cxx:327-328`.
        let my_start_wire = get_unified_wire(&a_start_wire, unifier.history())?;
        let my_end_wire = get_unified_wire(&an_end_wire, unifier.history())?;

        // `cxx:329-330`.
        let my_closed_section = my_start_wire.0.is_closed();
        let my_closed_ring = is_same(&my_start_wire.0, &my_end_wire.0);

        Ok(Self {
            my_initial_shape,
            my_start_wire,
            my_end_wire,
            my_closed_section,
            my_closed_ring,
            my_start_wire_edges: HashSet::new(),
            my_end_wire_edges: HashSet::new(),
            my_paths: Vec::new(),
            my_shape: None,
            my_done_flag: false,
        })
    }

    /// `void Build(const Message_ProgressRange&)` (`cxx:336-1148`).
    pub fn build(&mut self) -> Result<(), String> {
        let mut ctx = Ctx::default();
        // `cxx:339-497`.
        if !paths::collect(self, &mut ctx)? {
            // `cxx:416-421`: a start vertex without a free starting edge.
            return Ok(());
        }
        // `cxx:522-843`.
        sections::build_sections(self, &mut ctx)?;
        // `cxx:845-1147`.
        middle_edges::build_middle_path(self, &mut ctx)?;
        // `cxx:1147`.
        self.my_done_flag = true;
        Ok(())
    }

    /// `IsDone()` (`BRepBuilderAPI_Command`, inherited).
    pub fn is_done(&self) -> bool {
        self.my_done_flag
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`, inherited): the final wire.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `myPaths` (`hxx:53`) after `Build`.
    pub fn paths(&self) -> &[Vec<TopoShape>] {
        &self.my_paths
    }

    /// `myInitialShape` (`hxx:46`).
    pub fn initial_shape(&self) -> &TopoShape {
        &self.my_initial_shape
    }
}

/// `TopoDS::Wire(StartShape)` / `BRepTools::OuterWire(TopoDS::Face(StartShape))`
/// (`cxx:305-325`).
fn wire_of(shape: &TopoShape) -> Result<Wire, String> {
    if shape.shape_type() == crate::abs::ShapeType::Face {
        let face = Face::wrap(shape.clone())
            .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
        // `cxx:308-309` / `:318-319`.
        crate::brep_tools::outer_wire(&face)
            .ok_or_else(|| "Standard_NullObject: BRepTools::OuterWire".to_string())
    } else {
        // `cxx:313` / `:322`.
        Wire::wrap(shape.clone()).ok_or_else(|| "Standard_TypeMismatch: TopoDS::Wire".to_string())
    }
}
