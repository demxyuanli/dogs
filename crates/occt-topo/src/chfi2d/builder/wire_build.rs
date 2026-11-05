//! Wire and edge rebuilding: `BuildNewWire` (`ChFi2d_Builder.cxx:533-600`) and
//! both `BuildNewEdge` overloads (`cxx:604-681`).

use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_projecting};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi2d::ConstructionError;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_equal, is_same};

use super::ChFi2dBuilder;

impl ChFi2dBuilder {
    /// `BuildNewWire(OldE1, OldE2, E1, Fillet, E2)` (`cxx:535-600`): rebuild
    /// `newFace`'s wire with the two old edges replaced by their trimmed
    /// versions and the new fillet/chamfer edge in between.
    ///
    /// The fillet is appended after the first of the two replaced edges
    /// (`filletIsAdded` latch, `cxx:571-584`), and a degenerated side is simply
    /// dropped (`cxx:569`, `:581`).
    pub(super) fn build_new_wire(
        &mut self,
        old_e1: &Edge,
        old_e2: &Edge,
        e1: &Edge,
        fillet: &Edge,
        e2: &Edge,
    ) {
        // `cxx:542-549`: the closed flag comes from the reference face's wire.
        let closed_status = self.ref_wire_closed();

        let Some(face) = self.new_face.clone() else { return };
        // `cxx:553-557`: `Ex.Init(newFace, TopAbs_EDGE); B.MakeWire(newWire)`.
        let mut fillet_is_added = false;
        let mut edges: Vec<Edge> = Vec::new();
        for child in explorer_shapes(&face.0, crate::abs::ShapeType::Edge) {
            let Some(the_edge) = Edge::wrap(child) else { continue };
            if !is_same(&the_edge.0, &old_e1.0) && !is_same(&the_edge.0, &old_e2.0) {
                // `cxx:561-563`.
                edges.push(the_edge);
                continue;
            }
            if is_equal(&the_edge.0, &old_e1.0) {
                // `cxx:568-578`.
                if self.status != ConstructionError::FirstEdgeDegenerated
                    && self.status != ConstructionError::BothEdgesDegenerated
                {
                    edges.push(e1.clone());
                }
                if !fillet_is_added {
                    edges.push(fillet.clone());
                    fillet_is_added = true;
                }
            } else {
                // `cxx:579-591`.
                if self.status != ConstructionError::LastEdgeDegenerated
                    && self.status != ConstructionError::BothEdgesDegenerated
                {
                    edges.push(e2.clone());
                }
                if !fillet_is_added {
                    edges.push(fillet.clone());
                    fillet_is_added = true;
                }
            }
        }

        // `cxx:594-599`: `newWire.Closed(aClosedStatus)` then
        // `BRepLib_MakeFace(Adaptor3dSurface.Plane(), newWire)`.
        let wire = TopoBuilder::new().make_wire(&edges);
        wire.0.set_closed(closed_status);
        self.rebuild_face(&wire);
    }

    /// `BuildNewEdge(E1, OldExtr, NewExtr)` (`cxx:604-629`): the edge over
    /// `E1`'s curve with `OldExtr` moved to `NewExtr`.
    ///
    /// `None` is `BRepLib_MakeEdge::Init` failing to project a vertex
    /// (`BRepLib_PointProjectionFailed`), where OCCT returns the null edge of
    /// `makeEdge`.
    ///
    /// `Standard_EXPORT` in `ChFi2d_Builder.hxx:240-245`, so the port keeps it
    /// `pub`; OCCT itself has no live caller — both use sites are commented out
    /// (`_0.cxx:210-211`, `:246-247`) in favour of the `IsDegenerated` overload.
    pub fn build_new_edge(
        &self,
        e1: &Edge,
        old_extr: &TopoShape,
        new_extr: &TopoShape,
    ) -> Option<Edge> {
        // `cxx:611-617`.
        let curve = BRepTool::edge_curve(e1)?;
        let (first_vertex, last_vertex) = edge_vertices(e1);
        let first_vertex = crate::shape::Vertex(first_vertex?.0);
        let last_vertex = crate::shape::Vertex(last_vertex?.0);
        let new_vertex = crate::shape::Vertex(new_extr.clone());
        // `cxx:618-625`.
        let build = if is_same(&first_vertex.0, old_extr) {
            init_3d_projecting(&curve, Some(&new_vertex), Some(&last_vertex))
        } else {
            init_3d_projecting(&curve, Some(&first_vertex), Some(&new_vertex))
        }?;
        let b = TopoBuilder::new();
        let mut edge = b.make_edge(curve, build.first, build.last);
        add_endpoint_vertices(&b, &mut edge, &build);
        // `cxx:626`: `anEdge.Orientation(E1.Orientation())`.
        edge.0.set_orientation(e1.0.orientation());
        Some(edge)
    }

    /// `BuildNewEdge(E1, OldExtr, NewExtr, IsDegenerated)` (`cxx:636-680`): the
    /// same rebuild, but reporting an edge that collapses onto its other
    /// endpoint — `IsDegenerated = true` and `anEdge = E1` (`cxx:669-673`).
    pub(super) fn build_new_edge_or_degenerated(
        &self,
        e1: &Edge,
        old_extr: &TopoShape,
        new_extr: &TopoShape,
    ) -> Option<(Edge, bool)> {
        // `cxx:643-647`.
        let curve = BRepTool::edge_curve(e1)?;
        let (first_vertex, last_vertex) = edge_vertices(e1);
        let first_vertex = crate::shape::Vertex(first_vertex?.0);
        let last_vertex = crate::shape::Vertex(last_vertex?.0);
        let new_vertex = crate::shape::Vertex(new_extr.clone());
        // `cxx:648-649`: `Tol = Precision::Confusion()`.
        let tol = occt_core::precision::CONFUSION;
        let p_new = BRepTool::vertex_point(&new_vertex);

        // `cxx:650-665`: `PonctualEdge` compares `Pnew` with the point of the
        // *other* vertex of `E1`.
        let (build, other_point) = if is_same(&first_vertex.0, old_extr) {
            (
                init_3d_projecting(&curve, Some(&new_vertex), Some(&last_vertex)),
                BRepTool::vertex_point(&last_vertex),
            )
        } else {
            (
                init_3d_projecting(&curve, Some(&first_vertex), Some(&new_vertex)),
                BRepTool::vertex_point(&first_vertex),
            )
        };
        let ponctual_edge = p_new.distance(&other_point) < tol;

        // `cxx:667-679`: `error == BRepLib_LineThroughIdenticPoints ||
        // PonctualEdge` means degenerated (`anEdge = E1`); any other `Init`
        // failure (`BRepLib_PointProjectionFailed`) leaves OCCT with the null
        // shape of `makeEdge`, which the port answers with `None`.
        let Some(build) = build else {
            return if ponctual_edge { Some((e1.clone(), true)) } else { None };
        };
        if ponctual_edge {
            return Some((e1.clone(), true));
        }
        let b = TopoBuilder::new();
        let mut edge = b.make_edge(curve, build.first, build.last);
        add_endpoint_vertices(&b, &mut edge, &build);
        // `cxx:676`: `anEdge.Orientation(E1.Orientation())`.
        edge.0.set_orientation(e1.0.orientation());
        Some((edge, false))
    }
}
