//! `BRepFill_SectionLaw` — the shared state of every section law of a sweep.
//!
//! Source: `BRepFill_SectionLaw.cxx:33-111`, `.hxx:39-88`.
//!
//! OCCT declares this as an abstract class and lets each concrete law
//! (`BRepFill_ShapeLaw`, `BRepFill_NSections`, `BRepFill_EvolvedSection`, ...)
//! inherit its data. The port keeps that data in [`SectionLawBase`] (held by
//! each concrete law as a field) and the pure virtuals in the [`SectionLaw`]
//! trait, whose `base()` accessor reaches the shared state.
//!
//! The one piece of real work here is [`SectionLawBase::current_edge`]: it walks
//! the section wire through `BRepTools_WireExplorer` and *skips the degenerated
//! edges* (`cxx:86-111`), so the laws only ever see edges that carry a curve.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::toploc::LocationKey;
use occt_geom::geom_fill::SectionLaw as GeomFillSectionLaw;

use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, TopoShape, Vertex, Wire};

/// `TopTools_ShapeMapHasher` key: `TShape` identity plus location
/// (`myIndices`, `hxx:83`).
type ShapeKey = (usize, LocationKey);

/// `BRepFill_SectionLaw` (`BRepFill_SectionLaw.hxx:39-88`).
#[derive(Default)]
pub struct SectionLawBase {
    /// `myLaws` (`hxx:79`) — 1-based in OCCT, so index `i` is `laws[i - 1]`.
    /// `None` is a null handle, which OCCT leaves behind when a law is built
    /// with `Build = false` (`BRepFill_ShapeLaw.cxx:133-184`).
    laws: Vec<Option<Arc<dyn GeomFillSectionLaw>>>,
    /// `uclosed` (`hxx:80`).
    u_closed: bool,
    /// `vclosed` (`hxx:81`).
    v_closed: bool,
    /// `myDone` (`hxx:82`).
    done: bool,
    /// `myIndices` (`hxx:83`): edge -> its 1-based law index. Signed, because
    /// `BRepFill_NSections::Init` binds `-ii` for a `REVERSED` section edge
    /// (`BRepFill_NSections.cxx:489-495`) and `BRepFill_PipeShell::BuildHistory`
    /// reads that sign back (`BRepFill_PipeShell.cxx:1362-1366`).
    indices: HashMap<ShapeKey, i32>,
    /// `myIterator` (`hxx:86`). `RefCell` keeps `Init` / `CurrentEdge`
    /// callable through the `Arc<dyn SectionLaw>` handle that `BRepFill_Sweep`
    /// holds, which OCCT does with a plain mutable member.
    iterator: std::cell::RefCell<WireExplorer>,
}

impl SectionLawBase {
    pub fn new() -> Self {
        Self::default()
    }

    /// `NbLaw` (`cxx:33-37`): `myLaws->Length()`.
    pub fn nb_law(&self) -> usize {
        self.laws.len()
    }

    /// `Law(Index)` (`cxx:41-44`), 1-based as in OCCT.
    pub fn law(&self, index: usize) -> Option<Arc<dyn GeomFillSectionLaw>> {
        index
            .checked_sub(1)
            .and_then(|i| self.laws.get(i))
            .and_then(|slot| slot.clone())
    }

    /// `IndexOfEdge` (`cxx:48-51`): `myIndices(anEdge)`, i.e.
    /// `NCollection_DataMap::operator()`, which throws `Standard_NoSuchObject`
    /// for an edge the law does not know. The answer is signed: a section edge
    /// stored `REVERSED` binds `-ii` (`BRepFill_NSections.cxx:489-495`).
    pub fn index_of_edge(&self, an_edge: &TopoShape) -> Result<i32, String> {
        self.indices
            .get(&(an_edge.shape_addr(), an_edge.location().key()))
            .copied()
            .ok_or_else(|| "Standard_NoSuchObject: BRepFill_SectionLaw::IndexOfEdge".to_string())
    }

    /// `IsUClosed` (`cxx:55-58`).
    pub fn is_u_closed(&self) -> bool {
        self.u_closed
    }

    /// `IsVClosed` (`cxx:61-64`).
    pub fn is_v_closed(&self) -> bool {
        self.v_closed
    }

    /// `IsDone` (`cxx:68-71`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Init(W)` (`cxx:76-79`).
    pub fn init(&self, w: &Wire) {
        self.iterator.borrow_mut().init(w, None);
    }

    /// `CurrentEdge()` (`cxx:86-111`): the current non-degenerated edge, then
    /// advance past it.
    pub fn current_edge(&self) -> Option<Edge> {
        let mut iterator = self.iterator.borrow_mut();
        let mut e = None;
        // `cxx:92-96`.
        let mut suivant = false;
        if iterator.more() {
            e = iterator.current_edge();
            suivant = e.as_ref().map(|edge| BRepTool::is_degenerated(edge)).unwrap_or(false);
        }
        // `cxx:98-105`.
        while suivant {
            iterator.next();
            e = iterator.current_edge();
            suivant = e
                .as_ref()
                .map(|edge| BRepTool::is_degenerated(edge) && iterator.more())
                .unwrap_or(false);
        }
        // `cxx:107-110`.
        if iterator.more() {
            iterator.next();
        }
        e
    }

    // ---- state written by the concrete laws ----

    /// `myLaws` (1-based `Value(Index)` assignment).
    pub fn set_laws(&mut self, laws: Vec<Arc<dyn GeomFillSectionLaw>>) {
        self.laws = laws.into_iter().map(Some).collect();
    }

    /// `myLaws = new NCollection_HArray1<...>(1, Count)` (`BRepFill_ShapeLaw.cxx:129`):
    /// every slot starts as a null handle.
    pub fn reset_laws(&mut self, count: usize) {
        self.laws = vec![None; count];
    }

    /// `myLaws->ChangeValue(Index) = Law` (1-based), including the null
    /// assignment a `Build = false` law leaves behind.
    pub fn set_law_at(&mut self, index: usize, law: Option<Arc<dyn GeomFillSectionLaw>>) {
        if let Some(slot) = index.checked_sub(1).and_then(|i| self.laws.get_mut(i)) {
            *slot = law;
        }
    }

    /// `myIndices.Bind(anEdge, Index)`.
    pub fn bind_index(&mut self, an_edge: &TopoShape, index: i32) {
        self.indices.insert((an_edge.shape_addr(), an_edge.location().key()), index);
    }

    pub fn set_u_closed(&mut self, v: bool) {
        self.u_closed = v;
    }

    pub fn set_v_closed(&mut self, v: bool) {
        self.v_closed = v;
    }

    pub fn set_done(&mut self, v: bool) {
        self.done = v;
    }
}

/// The pure virtual part of `BRepFill_SectionLaw` (`hxx:50-70`).
pub trait SectionLaw {
    /// The shared state this law was built on.
    fn base(&self) -> &SectionLawBase;

    /// `IsConstant` (`hxx:50`).
    fn is_constant(&self) -> bool;

    /// `IsVertex` (`hxx:59`): the input shape is a vertex.
    fn is_vertex(&self) -> bool;

    /// `ConcatenedLaw` (`hxx:61`).
    fn concatened_law(&self) -> Option<Arc<dyn GeomFillSectionLaw>>;

    /// `Continuity` (`hxx:63`): `GeomAbs_Shape` as the port encodes it.
    fn continuity(&self, index: usize, tol_angular: f64) -> u8;

    /// `VertexTol` (`hxx:66`): the gap between two consecutive section edges.
    fn vertex_tol(&self, index: usize, param: f64) -> f64;

    /// `Vertex` (`hxx:68`).
    fn vertex(&self, index: usize, param: f64) -> Option<Vertex>;

    /// `D0(U, S)` (`hxx:70`): the section shape at parameter `u`.
    fn d0(&self, u: f64) -> Option<TopoShape>;
}
