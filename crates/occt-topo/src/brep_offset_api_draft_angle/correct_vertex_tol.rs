//! `BRepOffsetAPI_DraftAngle::CorrectVertexTol()`
//! (`BRepOffsetAPI_DraftAngle.cxx:913-993`, OCCT 8.0.0).
//!
//! A draft moves the vertices of the tuned faces, so some of the new edges are
//! now met by a vertex whose tolerance does not cover the gap. `CorrectVertexTol`
//! walks the edges of `myShape` that are not part of `myInitialShape` (and that
//! have not been seen yet), and for every vertex of such an edge that existed
//! in the initial shape:
//!
//! * if the vertex was already scheduled for replacement, raises the
//!   replacement's tolerance (`cxx:954-956`);
//! * else, when the vertex tolerance is below the edge tolerance, creates a
//!   new vertex at the same point with tolerance `ETol + Epsilon(ETol)` and
//!   records it in `myVtxToReplace` (`cxx:958-969`);
//! * a vertex that was *not* in the initial shape only has its tolerance
//!   raised (`cxx:972-974`).
//!
//! When nothing was scheduled the method returns with `myShape` untouched
//! (`cxx:978-982`); otherwise the replacements are pushed into `mySubs`, the
//! whole shape is re-applied through it and `myShape` becomes
//! `mySubs.Value(myShape)` (`cxx:984-992`).
//!
//! # Translation notes
//!
//! * `NCollection_Map<TopoDS_Shape>` (`cxx:915`) is a `HashSet` over
//!   [`shape_map_key`], the `TopTools_ShapeMapHasher` key.
//! * `TopoDS_Iterator(anEdge)` (`cxx:920`, `:949`) walks the edge's *stored*
//!   vertex children, which is [`children_keep_orientation`].
//! * `BRep_Builder::UpdateVertex(V, Tol)` (`cxx:956`, `:973`) is
//!   `TV->UpdateTolerance` (`BRep_TBuilder` raise-only), i.e.
//!   [`GeometryRegistry::raise_vertex_tolerance`].
//! * `BRep_Builder::MakeVertex(V, P, Tol)` is [`TopoBuilder::make_vertex`];
//!   `V.Orientation(anOri)` is `set_orientation`.
//! * `myVtxToReplace.IsBound(V)` / `Bind(V, NewV)` are
//!   [`BRepOffsetApiDraftAngle::vtx_replacement`] /
//!   [`BRepOffsetApiDraftAngle::bind_vtx_replacement`] (an insertion-ordered
//!   vector keyed by `IsSame`).
//! * `mySubs.Clear()` (`cxx:985`) has no `MapReShape::clear` in the port, so a
//!   fresh [`MapReShape`] replaces it; `mySubs.Apply(myShape)` answers the
//!   rebuilt shape, then `mySubs.Value(myShape)` gives the final one.

use std::collections::HashSet;

use occt_core::precision::epsilon;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex};
use crate::shape_fix_compose_shell::MapReShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{children_keep_orientation, explorer_shapes, shape_map_key};

use super::BRepOffsetApiDraftAngle;

/// `CorrectVertexTol()` (`cxx:913-993`).
impl BRepOffsetApiDraftAngle {
    pub(super) fn correct_vertex_tol(&mut self) {
        let bb = TopoBuilder::new();

        // `cxx:915-926`: the vertices and edges of the initial shape.
        let mut an_init_vertices: HashSet<_> = HashSet::new();
        let mut an_init_edges: HashSet<_> = HashSet::new();
        if let Some(initial) = self.my_initial_shape.clone() {
            for an_e in explorer_shapes(&initial, ShapeType::Edge) {
                an_init_edges.insert(shape_map_key(&an_e));
                // `cxx:920-924`: the edge's stored vertex children.
                for v in children_keep_orientation(&an_e) {
                    an_init_vertices.insert(shape_map_key(&v));
                }
            }
        }

        // `cxx:928-929`
        self.my_vtx_to_replace.clear();

        // `cxx:930-977`
        let mut a_new_edges: HashSet<_> = HashSet::new();
        let Some(shape) = self.my_shape.clone() else {
            return;
        };
        for an_e in explorer_shapes(&shape, ShapeType::Edge) {
            // `cxx:934-937`: skip the edges the draft did not touch.
            if an_init_edges.contains(&shape_map_key(&an_e)) {
                continue;
            }
            // `cxx:940-944`: skip edges already processed.
            if !a_new_edges.insert(shape_map_key(&an_e)) {
                continue;
            }
            // `cxx:946`
            let Some(edge) = Edge::wrap(an_e.clone()) else {
                continue;
            };
            let an_e_tol = BRepTool::edge_tolerance(&edge);
            let tol = an_e_tol + epsilon(an_e_tol);
            for v in children_keep_orientation(&an_e) {
                let Some(a_vtx) = Vertex::wrap(v) else {
                    continue;
                };
                // `cxx:953-971`
                if an_init_vertices.contains(&shape_map_key(&a_vtx.0)) {
                    if let Some(new_vtx) = self.vtx_replacement(&a_vtx.0) {
                        // `cxx:954-956`
                        GeometryRegistry::global().raise_vertex_tolerance(&new_vtx, tol);
                    } else {
                        // `cxx:958-969`
                        let a_v_tol = BRepTool::vertex_tolerance(&a_vtx);
                        if a_v_tol < an_e_tol {
                            let a_v_pnt = BRepTool::vertex_point(&a_vtx);
                            let mut new_vtx = bb.make_vertex(a_v_pnt, tol);
                            new_vtx.0.set_orientation(a_vtx.0.orientation());
                            self.bind_vtx_replacement(&a_vtx.0, &new_vtx.0);
                        }
                    }
                } else {
                    // `cxx:971-975`
                    GeometryRegistry::global().raise_vertex_tolerance(&a_vtx.0, tol);
                }
            }
        }

        // `cxx:979-982`
        if self.my_vtx_to_replace.is_empty() {
            return;
        }

        // `cxx:984-990`: `mySubs.Clear()` then replace every recorded vertex.
        self.my_subs = MapReShape::new();
        for (old, new) in self.my_vtx_to_replace.iter() {
            self.my_subs.replace(old, new);
        }

        // `cxx:991`: `mySubs.Apply(myShape); myShape = mySubs.Value(myShape);`
        let applied = self.my_subs.apply(&shape);
        let replaced = self.my_subs.value(&applied).unwrap_or(applied);
        self.my_shape = Some(replaced);
    }
}
