//! `Draft_VertexInfo` (`Draft_VertexInfo.cxx:27-117`): the new point of a
//! modified vertex plus the parameter it takes on every edge that meets there.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | ctor | `cxx:27` | [`DraftVertexInfo::new`] |
//! | `Add(E)` | `cxx:29-43` | [`DraftVertexInfo::add`] |
//! | `Geometry()` | `cxx:47-50` | [`DraftVertexInfo::geometry`] |
//! | `ChangeGeometry()` | `cxx:54-57` | [`DraftVertexInfo::change_geometry`] |
//! | `Parameter(E)` | `cxx:62-73` | [`DraftVertexInfo::parameter`] |
//! | `ChangeParameter(E)` | `cxx:77-89` | [`DraftVertexInfo::change_parameter`] |
//! | `InitEdgeIterator()` | `cxx:93-96` | [`DraftVertexInfo::init_edge_iterator`] |
//! | `Edge()` | `cxx:100-103` | [`DraftVertexInfo::edge`] |
//! | `MoreEdge()` | `cxx:107-110` | [`DraftVertexInfo::more_edge`] |
//! | `NextEdge()` | `cxx:114-117` | [`DraftVertexInfo::next_edge`] |
//!
//! # Translation notes
//!
//! * The two `NCollection_List`s (`hxx:57-58`) stay parallel vectors and the
//!   member iterator (`hxx:59`) is an index into `my_edges`, so
//!   `InitEdgeIterator` / `MoreEdge` / `NextEdge` keep OCCT's cursor semantics.
//! * `myParams` starts at `RealLast()` for every edge the vertex does not
//!   already carry (`cxx:41`), i.e. `f64::MAX`.
//! * `Parameter` / `ChangeParameter` throw `Standard_DomainError` when the edge
//!   is absent (`cxx:72`, `cxx:88`); the port panics with that name.

use occt_core::gp::GpPnt;

use crate::shape::{Edge, TopoShape};
use crate::topo_tools::is_same;

/// `class Draft_VertexInfo` (`Draft_VertexInfo.hxx:31-61`).
#[derive(Debug, Clone)]
pub struct DraftVertexInfo {
    /// `myGeom` (`hxx:56`): the new point, set by `Perform`.
    my_geom: GpPnt,
    /// `myEdges` (`hxx:57`).
    my_edges: Vec<TopoShape>,
    /// `myParams` (`hxx:58`), parallel to `my_edges`.
    my_params: Vec<f64>,
    /// `myItEd` (`hxx:59`).
    my_it_ed: usize,
}

impl Default for DraftVertexInfo {
    fn default() -> Self {
        Self::new()
    }
}

impl DraftVertexInfo {
    /// `Draft_VertexInfo()` (`cxx:27`).
    pub fn new() -> Self {
        Self {
            my_geom: GpPnt::new(0.0, 0.0, 0.0),
            my_edges: Vec::new(),
            my_params: Vec::new(),
            my_it_ed: 0,
        }
    }

    /// `Add(E)` (`cxx:29-43`): appends `E` with `RealLast()` unless the vertex
    /// already carries the same edge.
    pub fn add(&mut self, e: &Edge) {
        // `cxx:31-37`
        let found = self.my_edges.iter().any(|w| is_same(w, &e.0));
        // `cxx:38-42`
        if !found {
            self.my_edges.push(e.0.clone());
            self.my_params.push(f64::MAX);
        }
    }

    /// `Geometry()` (`cxx:47-50`).
    pub fn geometry(&self) -> &GpPnt {
        &self.my_geom
    }

    /// `ChangeGeometry()` (`cxx:54-57`).
    pub fn change_geometry(&mut self) -> &mut GpPnt {
        &mut self.my_geom
    }

    /// `Parameter(E)` (`cxx:62-73`): `Standard_DomainError` when `E` is absent.
    pub fn parameter(&self, e: &Edge) -> f64 {
        match self.index_of(e) {
            Some(i) => self.my_params[i],
            // `cxx:72`
            None => panic!("Standard_DomainError: Draft_VertexInfo::Parameter"),
        }
    }

    /// `ChangeParameter(E)` (`cxx:77-89`).
    pub fn change_parameter(&mut self, e: &Edge) -> &mut f64 {
        match self.index_of(e) {
            Some(i) => &mut self.my_params[i],
            // `cxx:88`
            None => panic!("Standard_DomainError: Draft_VertexInfo::ChangeParameter"),
        }
    }

    /// `InitEdgeIterator()` (`cxx:93-96`).
    pub fn init_edge_iterator(&mut self) {
        self.my_it_ed = 0;
    }

    /// `Edge()` (`cxx:100-103`): `TopoDS::Edge` on the cursor's value.
    pub fn edge(&self) -> Edge {
        Edge::wrap(self.my_edges[self.my_it_ed].clone())
            .expect("Standard_TypeMismatch: Draft_VertexInfo::Edge")
    }

    /// `MoreEdge()` (`cxx:107-110`).
    pub fn more_edge(&self) -> bool {
        self.my_it_ed < self.my_edges.len()
    }

    /// `NextEdge()` (`cxx:114-117`).
    pub fn next_edge(&mut self) {
        self.my_it_ed += 1;
    }

    /// Position of `E` in the parallel lists (`cxx:64-71`).
    fn index_of(&self, e: &Edge) -> Option<usize> {
        self.my_edges.iter().position(|w| is_same(w, &e.0))
    }
}
