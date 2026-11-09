//! `LocOpe_GeneratedShape` (`LocOpe_GeneratedShape.hxx:29-52`): the abstract
//! query interface a generator (`LocOpe_Generator`) drives.
//!
//! A concrete shape tells the generator which edges of the base shape bound the
//! generated feature (`GeneratingEdges`), which edge a base vertex produces
//! (`Generated(Vertex)`), which face a base edge produces
//! (`Generated(Edge)`), and which of those faces come back correctly oriented
//! (`OrientedFaces`). Both `Generated` overloads answer a *null* shape when the
//! sub-shape produces nothing (`hxx:35-42`); the port spells that `None`.
//!
//! The reference's protected `myGEdges` / `myList` members (`hxx:50-51`) are
//! caches its subclasses fill (`LocOpe_Prism`, `LocOpe_Revol`, the two Forms);
//! they are not part of the interface and do not appear here.

use crate::shape::TopoShape;

/// `LocOpe_GeneratedShape` (`hxx:29-46`).
pub trait GeneratedShape {
    /// `GeneratingEdges()` (`hxx:32`): the edges of the base shape that bound
    /// the generated feature.
    fn generating_edges(&self) -> &[TopoShape];

    /// `Generated(const TopoDS_Vertex&)` (`hxx:36`): the edge the vertex
    /// produces, or `None` for the reference's null shape.
    fn generated_edge(&self, v: &TopoShape) -> Option<TopoShape>;

    /// `Generated(const TopoDS_Edge&)` (`hxx:41`): the face the edge
    /// produces, or `None` for the reference's null shape.
    fn generated_face(&self, e: &TopoShape) -> Option<TopoShape>;

    /// `OrientedFaces()` (`hxx:45`): the generated faces, correctly oriented.
    fn oriented_faces(&self) -> &[TopoShape];
}
