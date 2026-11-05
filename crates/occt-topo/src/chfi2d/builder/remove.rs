//! `RemoveFillet` (`ChFi2d_Builder.cxx:290-497`) and `RemoveChamfer`
//! (`ChFi2d_Builder_0.cxx:310-519`).
//!
//! The two bodies are the same apart from one line: the rebuilt edges take
//! their orientation and location from `E1` / `E2` in `RemoveFillet`
//! (`Builder.cxx:401-402`, `:414-415`) and from `basisEdge1` / `basisEdge2` in
//! `RemoveChamfer` (`_0.cxx:423-424`, `:436-437`). That difference is
//! [`OrientSource`].
//!
//! Translation notes:
//! - `BRepLib_MakeEdge::Init(curve, V1, V2)` is
//!   [`crate::brep_lib_make_edge::init_3d_projecting`], which is the same call
//!   [`super::wire_build`] uses for `BuildNewEdge`,
//! - `newFace.Nullify(); newFace = mFace;` (`Builder.cxx:490-492`) is
//!   [`ChFi2dBuilder::rebuild_face`],
//! - OCCT leaves `newEdge1` / `newEdge2` null when neither
//!   `firstVertex.IsSame(connectionE1Fillet)` nor
//!   `lastVertex.IsSame(connectionE1Fillet)` holds, and then adds that null
//!   edge to the wire (`Builder.cxx:474-486`). A null shape cannot be stored
//!   here, so that arm answers `None` — the same null vertex the other early
//!   returns answer.

use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_equal, is_same};
use crate::{brep_lib_make_edge, builder::TopoBuilder};

use crate::brep_tool::BRepTool;
use crate::abs::ShapeType;
use crate::chfi2d::{common_vertex, find_connected_edges, ConstructionError};
use crate::topo_tools::explorer_shapes;

use super::ChFi2dBuilder;

/// Whether the rebuilt edges inherit their orientation and location from the
/// trimmed edge (fillets) or from its basis edge (chamfers).
#[derive(Clone, Copy, PartialEq, Eq)]
enum OrientSource {
    /// `newEdge.Orientation(E.Orientation()); newEdge.Location(E.Location());`
    Trimmed,
    /// `newEdge.Orientation(basisEdge.Orientation());
    ///  newEdge.Location(basisEdge.Location());`
    Basis,
}

impl ChFi2dBuilder {
    /// `RemoveFillet(Fillet)` (`Builder.cxx:290-497`): drop the fillet from the
    /// face and restore the two adjacent edges over their basis edges.
    ///
    /// Returns the vertex the fillet was attached to, `None` where OCCT returns
    /// a null `TopoDS_Vertex`.
    pub fn remove_fillet(&mut self, fillet: &Edge) -> Option<TopoShape> {
        self.remove(fillet, true, OrientSource::Trimmed)
    }

    /// `RemoveChamfer(Chamfer)` (`_0.cxx:310-519`).
    pub fn remove_chamfer(&mut self, chamfer: &Edge) -> Option<TopoShape> {
        self.remove(chamfer, false, OrientSource::Basis)
    }

    /// The shared body of the two removals.
    fn remove(&mut self, target: &Edge, is_fillet: bool, source: OrientSource) -> Option<TopoShape> {
        // `Builder.cxx:293-308` / `_0.cxx:314-331`: drop the edge from the
        // list, if it is one of them.
        let list = if is_fillet { &mut self.fillets } else { &mut self.chamfers };
        let position = list.iter().position(|s| is_same(s, &target.0))?;
        list.remove(position);

        let face = self.new_face.clone()?;
        // `TopExp::Vertices(Fillet, firstVertex, lastVertex)`
        // (`Builder.cxx:311-312`).
        let (first_vertex, last_vertex) = edge_vertices(target);
        let (first_vertex, last_vertex) = (first_vertex?, last_vertex?);

        // `Builder.cxx:314-333`: `E1` is the neighbour of the first vertex that
        // is not the fillet itself.
        let (adj1, adj2) = match find_connected_edges(&face, &first_vertex.0) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        let e1 = if is_same(&adj1.0, &target.0) { adj2 } else { adj1 };
        let basis_edge1 = self.basis_edge(&e1);

        // `Builder.cxx:334-346`: the same at the last vertex.
        let (adj1, adj2) = match find_connected_edges(&face, &last_vertex.0) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        let e2 = if is_same(&adj1.0, &target.0) { adj2 } else { adj1 };
        let basis_edge2 = self.basis_edge(&e2);

        // `Builder.cxx:347-366`: the three connection checks.
        let common = match common_vertex(&basis_edge1, &basis_edge2) {
            Some(v) => v,
            None => {
                self.status = ConstructionError::ConnexionError;
                return None;
            }
        };
        let connection1 = match common_vertex(&e1, target) {
            Some(v) => v,
            None => {
                self.status = ConstructionError::ConnexionError;
                return None;
            }
        };
        let connection2 = match common_vertex(&e2, target) {
            Some(v) => v,
            None => {
                self.status = ConstructionError::ConnexionError;
                return None;
            }
        };

        // `Builder.cxx:368-466`: rebuild both edges.
        let new_edge1 = self.rebuild_without(&e1, &basis_edge1, &common, &connection1, source)?;
        let new_edge2 = self.rebuild_without(&e2, &basis_edge2, &common, &connection2, source)?;

        // `Builder.cxx:467-488`: the wire with `E1` / `E2` replaced and the
        // fillet dropped. `theEdge == E1` is `TopoDS_Shape::operator==`, i.e.
        // `IsSame` plus orientation and location, which is
        // [`crate::topo_tools::is_equal`].
        let mut edges: Vec<Edge> = Vec::new();
        for child in explorer_shapes(&face.0, ShapeType::Edge) {
            let Some(the_edge) = Edge::wrap(child) else { continue };
            if !is_same(&the_edge.0, &e1.0)
                && !is_same(&the_edge.0, &e2.0)
                && !is_same(&the_edge.0, &target.0)
            {
                edges.push(the_edge);
            } else if is_equal(&the_edge.0, &e1.0) {
                edges.push(new_edge1.clone());
            } else if is_equal(&the_edge.0, &e2.0) {
                edges.push(new_edge2.clone());
            }
        }
        // `Builder.cxx:489-493`.
        let wire = TopoBuilder::new().make_wire(&edges);
        self.rebuild_face(&wire);

        // `Builder.cxx:495`: `UpDateHistory(basisEdge1, basisEdge2, newEdge1,
        //  newEdge2)`.
        self.replace_history(&basis_edge1, &basis_edge2, &new_edge1, &new_edge2);

        // `Builder.cxx:496`: `return commonVertex;`.
        Some(common)
    }

    /// One of the two rebuilds of `RemoveFillet` (`Builder.cxx:369-416`,
    /// `:418-463`): the basis edge when the trimmed edge carries only this one
    /// fillet, otherwise an edge over `E`'s curve trimmed to the far side of
    /// the fillet.
    fn rebuild_without(
        &self,
        e: &Edge,
        basis: &Edge,
        common: &TopoShape,
        connection: &TopoShape,
        source: OrientSource,
    ) -> Option<Edge> {
        let (first_vertex, last_vertex) = edge_vertices(e);
        let (first_vertex, last_vertex) = (first_vertex?, last_vertex?);
        // `TopExp::Vertices(basisEdge1, v1, v2); v = the one that is not the
        //  common vertex` (`Builder.cxx:372-386`).
        let (v1, v2) = edge_vertices(basis);
        let (v1, v2) = (v1?, v2?);
        let v = if is_same(&v1.0, common) { v2 } else { v1 };

        // `if (firstVertex.IsSame(v) || lastVertex.IsSame(v)) newEdge1 =
        //  basisEdge1;` (`Builder.cxx:388-393`).
        if is_same(&first_vertex.0, &v.0) || is_same(&last_vertex.0, &v.0) {
            return Some(basis.clone());
        }
        // `Builder.cxx:396-416`: the edge supports one fillet on each end.
        let (from, to): (&TopoShape, &TopoShape) = if is_same(&first_vertex.0, connection) {
            (common, &last_vertex.0)
        } else if is_same(&last_vertex.0, connection) {
            (&first_vertex.0, common)
        } else {
            // See the module note: OCCT keeps a null `newEdge1` here.
            return None;
        };
        // `BRep_Tool::Curve(E1, loc, first, last); makeEdge.Init(curve, from,
        //  to);` (`Builder.cxx:399-402`).
        let curve = BRepTool::edge_curve(e)?;
        let build = brep_lib_make_edge::init_3d_projecting(
            &curve,
            Some(&Vertex(from.clone())),
            Some(&Vertex(to.clone())),
        )?;
        let builder = TopoBuilder::new();
        let mut edge = builder.make_edge(curve, build.first, build.last);
        brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
        // `newEdge1.Orientation(...); newEdge1.Location(...);`.
        let model = match source {
            OrientSource::Trimmed => e,
            OrientSource::Basis => basis,
        };
        edge.0.set_orientation(model.0.orientation());
        edge.0.set_location(model.0.location());
        Some(edge)
    }
}

