//! `BRepFilletAPI_LocalOperation` — the abstract contour / simulation
//! interface shared by `BRepFilletAPI_MakeFillet` and
//! `BRepFilletAPI_MakeChamfer`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepFilletAPI/BRepFilletAPI_LocalOperation.hxx:33-99`
//! (OCCT 8.0.0). Every member is `virtual ... = 0` there, so the class has no
//! body of its own: it only fixes the signatures the two subclasses answer
//! (`BRepFilletAPI_MakeFillet.cxx:276-432`,
//! `BRepFilletAPI_MakeChamfer.cxx:180-371`).
//!
//! Translation notes:
//! - OCCT's `void` members that answer by throwing (`Add`, `ResetContour`,
//!   `Remove`, `Simulate` reach `ChFi3d_Builder` / `ChFiDS_Spine` failures)
//!   answer `Result<(), String>` here, the port's standing spelling of
//!   "void + `Standard_Failure`",
//! - `const TopoDS_Edge& Edge(I, J)` becomes `Option<TopoShape>`: `I` and `J`
//!   index `Value(I)->Edges(J)` and OCCT raises `Standard_OutOfRange` on a
//!   missing contour or edge; `None` is that arm without a panic. The edge is
//!   returned by value because the port's spine lives behind a lock guard,
//!   while a `TopoShape` clone shares the same `TShape`,
//! - `occ::handle<NCollection_HArray1<ChFiDS_CircSection>> Sect(IC, IS)` is
//!   `Option<Arc<CircSectionArray>>`, the port's null-handle spelling,
//! - `TopoDS_Vertex FirstVertex/LastVertex(IC)` answer `Option<TopoShape>` for
//!   OCCT's null vertex.

use std::sync::Arc;

use crate::chfi_ds::CircSectionArray;
use crate::shape::{Edge, TopoShape, Vertex};

/// `BRepFilletAPI_LocalOperation` (`hxx:33-99`).
pub trait LocalOperation {
    /// `Add(E)` (`hxx:40`): builds a contour of tangent edges.
    fn add(&mut self, e: &Edge) -> Result<(), String>;

    /// `ResetContour(IC)` (`hxx:44`): empties the contour of index `IC`.
    fn reset_contour(&mut self, ic: usize) -> Result<(), String>;

    /// `NbContours()` (`hxx:47`).
    fn nb_contours(&self) -> usize;

    /// `Contour(E)` (`hxx:51`): the index of the contour holding `E`, `0` when
    /// `E` belongs to none.
    fn contour(&self, e: &Edge) -> usize;

    /// `NbEdges(I)` (`hxx:55`).
    fn nb_edges(&self, i: usize) -> usize;

    /// `Edge(I, J)` (`hxx:58`).
    fn edge(&self, i: usize, j: usize) -> Option<TopoShape>;

    /// `Remove(E)` (`hxx:61`): drops the contour holding `E`.
    fn remove(&mut self, e: &Edge) -> Result<(), String>;

    /// `Length(IC)` (`hxx:64`).
    fn length(&self, ic: usize) -> f64;

    /// `FirstVertex(IC)` (`hxx:67`).
    fn first_vertex(&self, ic: usize) -> Option<TopoShape>;

    /// `LastVertex(IC)` (`hxx:70`).
    fn last_vertex(&self, ic: usize) -> Option<TopoShape>;

    /// `Abscissa(IC, V)` (`hxx:74`).
    fn abscissa(&self, ic: usize, v: &Vertex) -> f64;

    /// `RelativeAbscissa(IC, V)` (`hxx:77`).
    fn relative_abscissa(&self, ic: usize, v: &Vertex) -> f64;

    /// `ClosedAndTangent(IC)` (`hxx:81`).
    fn closed_and_tangent(&self, ic: usize) -> bool;

    /// `Closed(IC)` (`hxx:84`).
    fn closed(&self, ic: usize) -> bool;

    /// `Reset()` (`hxx:90`): undoes `Build` and leaves the contours and radii
    /// editable again.
    fn reset(&mut self);

    /// `Simulate(IC)` (`hxx:92`).
    fn simulate(&mut self, ic: usize) -> Result<(), String>;

    /// `NbSurf(IC)` (`hxx:94`).
    fn nb_surf(&self, ic: usize) -> usize;

    /// `Sect(IC, IS)` (`hxx:96-98`).
    fn sect(&self, ic: usize, is: usize) -> Option<Arc<CircSectionArray>>;
}
