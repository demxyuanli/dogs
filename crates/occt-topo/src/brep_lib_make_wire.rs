//! `BRepLib_MakeWire` — the wire builder that **decides the orientation of every
//! added edge**; it is what `BRepBuilderAPI_MakeWire` (`BRepBuilderAPI_MakeWire.cxx:117-125`)
//! delegates to.
//!
//! Source: `BRepLib_MakeWire.cxx:35-38` (ctor, `myError = EmptyWire`),
//! `:42-77` (ctors from edges), `:96-106` (`Add(Wire)`),
//! `:110-113` (`Add(Edge)` → `Add(Edge, true)`),
//! `:123-453` (`Add(Edge, IsCheckGeometryProximity)` — the whole control flow),
//! `:457-460` (`Wire`), `:464-467` (`Edge`), `:471-474` (`Vertex`),
//! `:485-488` (`Error`); the error enum is `BRepLib_WireError`.
//!
//! Supporting OCCT control flow reproduced here:
//! - `TopoDS_Iterator(shape, CumOri = true)` — children in stored order with the
//!   parent orientation composed in; the algorithm iterates
//!   `E.Oriented(TopAbs_FORWARD)`, i.e. the **stored** child orientations.
//! - `TopExp::Vertices(E, V1, V2, CumOri = true)` (`TopExp.cxx`) — `V1`/`V2` are
//!   the children seen FORWARD/REVERSED, so the edge's own orientation swaps
//!   them ([`edge_first_last_vertices`]).
//! - `TopTools_ShapeMapHasher` (`TopTools_ShapeMapHasher.hxx:35-38`) — equality is
//!   `TopoDS_Shape::IsSame` (`TopoDS_Shape.hxx`: same TShape handle **and** same
//!   location; two null handles compare equal). [`is_same`] reproduces that,
//!   with the location part reduced to "identity" as the rest of this port does
//!   (`bop_occt_util.rs:31`: "IsSame without location"), which holds on every
//!   wire-building path here.
//!
//! **Why this module exists (task T-32)**: `TopoBuilder::make_wire` is the
//! `BRep_Builder` level append (`BRep_Builder::MakeWire` + `Add`, no orientation
//! control). The shared edges of a triangle mesh are cached with the direction
//! of whichever triangle created them first, so appending them naively yields a
//! wire whose edges are **not chained** (board round 8: "edges all Forward,
//! source/sink vertices"). Building such wires through this class reproduces
//! OCCT: an edge **occurrence** is reversed when needed, leaving the shared
//! `TShape` untouched so neighbouring faces are unaffected.
//!
//! **Boundary of the port**: OCCT's `TopoDS_Shape::EmptyCopied` +
//! `BRep_Builder::Transfert` in the copy-edge branch are reproduced by rebuilding
//! an edge from the same curve and parameter range (`TopoBuilder::make_edge`) and
//! adding the resolved vertices with `TopoBuilder::add_edge_vertices` — exactly
//! the `BRepLib_MakeEdge` convention (first vertex FORWARD, last REVERSED) that
//! `B.Add(myEdge, myVertex)` reproduces for an edge whose children are in
//! storage order.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex, Wire};
use crate::tgeometry::{GeometryRegistry, VertexGeom};
use crate::topo_tools_full::edge_vertices;

/// `BRepLib_WireError` (`BRepLib_MakeWire.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WireError {
    /// No edge added yet (`BRepLib_MakeWire.cxx:36`).
    EmptyWire,
    /// The added edge shares no vertex with the wire (`:290`).
    DisconnectedWire,
    /// Ambiguous / non-manifold connection (`:167`, `:189`, `:218`, `:274`, `:410`, `:436`).
    NonManifoldWire,
    /// The last `Add` succeeded (`:451`).
    WireDone,
}

/// `TopoDS_Shape::IsSame` (`TopoDS_Shape.hxx`) on optional (possibly null)
/// handles: same `TShape`, and two null handles compare equal.
fn is_same(a: Option<&Vertex>, b: Option<&Vertex>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => x.0.same_tshape(&y.0),
        _ => false,
    }
}

/// Children of `E.Oriented(TopAbs_FORWARD)` as `TopoDS_Iterator(shape, true)`
/// sees them: stored order, stored orientations.
fn edge_forward_vertices(e: &Edge) -> Vec<Vertex> {
    let stored = e
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    stored
        .into_iter()
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .map(Vertex)
        .collect()
}

/// `TopExp::Vertices(E, V1, V2, CumOri = true)` (`TopExp.cxx`): `V1`/`V2` are the
/// stored FORWARD/REVERSED children, swapped when the edge occurrence is
/// REVERSED.
fn edge_first_last_vertices(e: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (first, last) = edge_vertices(e);
    if e.0.orientation() == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}

/// Port of `BRepLib_MakeWire` (`BRepLib_MakeWire.cxx`).
pub struct MakeWire {
    /// `myShape` — the wire under construction.
    shape: Wire,
    /// `myEdge` — the edge added last.
    edge: Option<Edge>,
    /// `myVertices` — `TopTools_IndexedMapOfShape` of the wire's vertices, in
    /// insertion order (`FindKey(i)` walks it; `Add` does not duplicate).
    vertices: Vec<Vertex>,
    /// `myVertex` — the vertex the last edge was connected by.
    vertex: Option<Vertex>,
    /// `VF` / `VL` — the wire's first / last vertex (`TopExp::Vertices`).
    first_vertex: Option<Vertex>,
    last_vertex: Option<Vertex>,
    /// `myError` — persistent across `Add` calls.
    error: WireError,
    /// `IsDone()` (`BRepLib_MakeShape`).
    done: bool,
}

impl Default for MakeWire {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeWire {
    /// `BRepLib_MakeWire::BRepLib_MakeWire()` (`cxx:35-38`).
    pub fn new() -> Self {
        Self {
            shape: Wire::new(),
            edge: None,
            vertices: Vec::new(),
            vertex: None,
            first_vertex: None,
            last_vertex: None,
            error: WireError::EmptyWire,
            done: false,
        }
    }

    /// `BRepLib_MakeWire::Add(const TopoDS_Edge&)` (`cxx:110-113`).
    pub fn add(&mut self, e: &Edge) -> Result<(), WireError> {
        self.add_checked(e, true)
    }

    /// `BRepLib_MakeWire::Add(const TopoDS_Edge&, bool IsCheckGeometryProximity)`
    /// (`cxx:123-453`).
    ///
    /// `Ok(())` mirrors reaching `Done()` at `:452`; `Err` carries the `myError`
    /// value OCCT leaves on the early `return` at `:291-292` (in
    /// `BRepBuilderAPI_MakeWire` terms: `IsDone() == false`).
    pub fn add_checked(
        &mut self,
        e: &Edge,
        is_check_geometry_proximity: bool,
    ) -> Result<(), WireError> {
        let builder = TopoBuilder::new();

        // to tell if it has been decided to add forward / reversed (`cxx:126-130`)
        let mut forward = false;
        let mut reverse = false;
        let init;

        if self.edge.is_none() {
            // First edge, create the wire (`cxx:135-149`). `myEdge = E` keeps the
            // **input** orientation (`cxx:142`).
            init = true;
            self.edge = Some(e.clone());
            for v in edge_forward_vertices(e) {
                self.push_vertex(v);
            }
        } else {
            // `init = myShape.Closed();` (`cxx:153`)
            init = self.shape.closed();
            // `TopoDS_Shape aLocalShape = E.Oriented(TopAbs_FORWARD);` (`cxx:154-155`)
            let ee = Edge(e.oriented(Orientation::Forward));

            let mut connected = false;
            let mut copyedge = false;

            // `if (VF.IsNull() || VL.IsNull()) myError = NonManifoldWire;` (`cxx:165-169`)
            if self.error != WireError::NonManifoldWire
                && (self.first_vertex.is_none() || self.last_vertex.is_none())
            {
                self.error = WireError::NonManifoldWire;
            }

            for ve in edge_forward_vertices(&ee) {
                // `if (myVertices.Contains(VE))` (`cxx:177`)
                if self.contains_vertex(&ve) {
                    connected = true;
                    self.vertex = Some(ve.clone());
                    if self.error != WireError::NonManifoldWire {
                        if is_same(self.first_vertex.as_ref(), self.last_vertex.as_ref()) {
                            // Orientation indetermined (in 3d): preserve the initial (`cxx:184-191`)
                            if !is_same(self.first_vertex.as_ref(), Some(&ve)) {
                                self.error = WireError::NonManifoldWire;
                            }
                        } else if is_same(self.first_vertex.as_ref(), Some(&ve)) {
                            if ve.0.orientation() == Orientation::Forward {
                                reverse = true;
                            } else {
                                forward = true;
                            }
                        } else if is_same(self.last_vertex.as_ref(), Some(&ve)) {
                            if ve.0.orientation() == Orientation::Reversed {
                                reverse = true;
                            } else {
                                forward = true;
                            }
                        } else {
                            self.error = WireError::NonManifoldWire;
                        }
                    }
                } else if is_check_geometry_proximity {
                    // Search a similar vertex in the edge (`cxx:223-285`).
                    let pe = BRepTool::vertex_point(&ve);
                    for i in 0..self.vertices.len() {
                        let vw = self.vertices[i].clone();
                        let pw = BRepTool::vertex_point(&vw);
                        let l = pe.distance(&pw);
                        if l < BRepTool::vertex_tolerance(&ve)
                            || l < BRepTool::vertex_tolerance(&vw)
                        {
                            copyedge = true;
                            if self.error != WireError::NonManifoldWire {
                                if is_same(self.first_vertex.as_ref(), self.last_vertex.as_ref()) {
                                    if !is_same(self.first_vertex.as_ref(), Some(&vw)) {
                                        self.error = WireError::NonManifoldWire;
                                    }
                                } else if is_same(self.first_vertex.as_ref(), Some(&vw)) {
                                    if ve.0.orientation() == Orientation::Forward {
                                        reverse = true;
                                    } else {
                                        forward = true;
                                    }
                                } else if is_same(self.last_vertex.as_ref(), Some(&vw)) {
                                    if ve.0.orientation() == Orientation::Reversed {
                                        reverse = true;
                                    } else {
                                        forward = true;
                                    }
                                } else {
                                    self.error = WireError::NonManifoldWire;
                                }
                            }
                            break;
                        }
                    }
                    if copyedge {
                        connected = true;
                    }
                }
            }

            if !connected {
                // `myError = BRepLib_DisconnectedWire; NotDone(); return;` (`cxx:288-293`)
                self.error = WireError::DisconnectedWire;
                self.done = false;
                return Err(WireError::DisconnectedWire);
            }

            if !copyedge {
                // `myEdge = EE; myVertices.Add(it.Value());` (`cxx:296-303`)
                self.edge = Some(ee.clone());
                for v in edge_forward_vertices(&ee) {
                    self.push_vertex(v);
                }
            } else {
                // Copy the edge (`cxx:304-368`); `myVertex = VW` inside (`cxx:353-356`).
                let (copied, merged) = self.copy_edge_with_merged_vertices(&ee);
                self.edge = Some(copied);
                if merged.is_some() {
                    self.vertex = merged;
                }
            }

            // Make a decision about the orientation of the edge (`cxx:370-379`):
            // an ambiguous case (nothing decided) preserves the input orientation.
            let e_is_reversed = e.0.orientation() == Orientation::Reversed;
            if ((forward == reverse) && e_is_reversed) || (reverse && !forward) {
                if let Some(me) = self.edge.as_mut() {
                    me.0.reverse();
                }
            }
        }

        // Add myEdge to myShape (`cxx:382-384`).
        let me = self.edge.clone().expect("myEdge is set in both branches");
        builder.add_edge(&mut self.shape, &me);
        self.shape.set_closed(false);

        // Initialize / update VF, VL (`cxx:386-444`).
        if init {
            // `TopExp::Vertices(TopoDS::Wire(myShape), VF, VL);` (`cxx:387-390`)
            let (v1, v2) = edge_first_last_vertices(&me);
            self.first_vertex = v1;
            self.last_vertex = v2;
        } else if self.error == WireError::WireDone {
            // `TopoDS_Vertex V1, V2, VRef; TopExp::Vertices(myEdge, V1, V2);` (`cxx:393-411`)
            let (vv1, vv2) = edge_first_last_vertices(&me);
            let v_ref = if is_same(vv1.as_ref(), self.vertex.as_ref()) {
                vv2
            } else if is_same(vv2.as_ref(), self.vertex.as_ref()) {
                vv1
            } else {
                self.error = WireError::NonManifoldWire;
                None
            };

            if is_same(self.first_vertex.as_ref(), self.last_vertex.as_ref()) {
                // Particular case: it is required to control the orientation
                // (`cxx:413-420`; OCCT only reports `VF == myVertex` in debug).
            } else if is_same(self.first_vertex.as_ref(), self.vertex.as_ref()) {
                self.first_vertex = v_ref;
            } else if is_same(self.last_vertex.as_ref(), self.vertex.as_ref()) {
                self.last_vertex = v_ref;
            } else {
                self.error = WireError::NonManifoldWire;
            }
        }
        if self.error == WireError::NonManifoldWire {
            // `VF = VL = TopoDS_Vertex();` (`cxx:440-443`)
            self.first_vertex = None;
            self.last_vertex = None;
        }

        // `if (!VF.IsNull() && !VL.IsNull() && VF.IsSame(VL)) myShape.Closed(true);` (`cxx:445-449`)
        if self.first_vertex.is_some()
            && self.last_vertex.is_some()
            && is_same(self.first_vertex.as_ref(), self.last_vertex.as_ref())
        {
            self.shape.set_closed(true);
        }

        self.error = WireError::WireDone;
        self.done = true;
        Ok(())
    }

    /// `myVertices.Add(...)` on a `TopTools_IndexedMapOfShape` (`cxx:147`, `:301`):
    /// insertion, no duplicates.
    fn push_vertex(&mut self, v: Vertex) {
        if !self.contains_vertex(&v) {
            self.vertices.push(v);
        }
    }

    /// `myVertices.Contains(VE)` (`TopTools_ShapeMapHasher` = `IsSame`).
    fn contains_vertex(&self, v: &Vertex) -> bool {
        self.vertices.iter().any(|w| w.0.same_tshape(&v.0))
    }

    /// The copy-edge branch of `Add` (`BRepLib_MakeWire.cxx:304-368`):
    /// `EE.EmptyCopied()` and then, for every child vertex `VE`, either the
    /// matching wire vertex `VW` — updated to the weighted point with
    /// `maxtol = 0.5 * (tolW + tolE + l)` (`B.UpdateVertex`) — or `VE` itself is
    /// added to the copy with `VE`'s orientation (`B.Add` + `B.Transfert`).
    ///
    /// Returns the copy and the connection vertex `myVertex` (`cxx:353-356`).
    fn copy_edge_with_merged_vertices(&self, ee: &Edge) -> (Edge, Option<Vertex>) {
        let builder = TopoBuilder::new();
        let curve = BRepTool::edge_curve(ee).expect("copy_edge: edge has no curve");
        let (first, last) = BRepTool::edge_parameters(ee);
        let mut copy = builder.make_edge(curve, first, last);

        let mut resolved: Vec<Vertex> = Vec::new();
        let mut merged_vertex: Option<Vertex> = None;
        for ve in edge_forward_vertices(ee) {
            let pe = BRepTool::vertex_point(&ve);
            let tol_e = BRepTool::vertex_tolerance(&ve);
            let mut merged = None;
            for vw in &self.vertices {
                let pw = BRepTool::vertex_point(vw);
                let l = pe.distance(&pw);
                let tol_w = BRepTool::vertex_tolerance(vw);
                if l < tol_e || l < tol_w {
                    // `maxtol = .5 * (tolW + tolE + l)`, then the weighted point
                    // (`cxx:328-349`).
                    let maxtol = 0.5 * (tol_w + tol_e + l);
                    let (c_w, c_e, tol) = if maxtol > tol_w && maxtol > tol_e {
                        let c_w = (maxtol - tol_e) / l;
                        (c_w, 1.0 - c_w, maxtol)
                    } else if maxtol > tol_w {
                        (0.0, 1.0, tol_e)
                    } else {
                        (1.0, 0.0, tol_w)
                    };
                    let pc = occt_core::gp::GpPnt::new(
                        c_w * pw.x() + c_e * pe.x(),
                        c_w * pw.y() + c_e * pe.y(),
                        c_w * pw.z() + c_e * pe.z(),
                    );
                    // `B.UpdateVertex(VW, PC, maxtol)` — the wire vertex is shared,
                    // so this updates it for the whole wire (as OCCT does).
                    GeometryRegistry::global()
                        .set_vertex(&vw.0, VertexGeom { point: pc, tolerance: tol });
                    let mut v = vw.clone();
                    v.0.set_orientation(ve.0.orientation());
                    merged = Some(v);
                    break;
                }
            }
            match merged {
                Some(v) => {
                    merged_vertex = Some(v.clone());
                    resolved.push(v);
                }
                // `myVertices.Add(VE); B.Add(myEdge, VE); B.Transfert(EE, myEdge, VE, VE);`
                // (`cxx:361-366`) — the port's geometry side table already holds
                // `VE`'s data, so the vertex is used as is.
                None => resolved.push(ve),
            }
        }

        // `B.Add(myEdge, myVertex)` in child order; `add_edge_vertices` writes
        // exactly the `BRepLib_MakeEdge` convention (first FORWARD, last
        // REVERSED) that the loop above restored.
        if resolved.len() == 2 {
            builder.add_edge_vertices(&mut copy, &resolved[0], &resolved[1]);
        } else {
            for v in &resolved {
                builder.add(&mut copy.0, &v.0);
            }
        }
        (copy, merged_vertex)
    }

    /// `IsDone()` (`BRepLib_MakeShape`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Error()` (`BRepLib_MakeWire.cxx:485-488`).
    pub fn error(&self) -> WireError {
        self.error
    }

    /// `Wire()` (`BRepLib_MakeWire.cxx:457-460`).
    pub fn wire(&self) -> Wire {
        self.shape.clone()
    }

    /// `Edge()` (`BRepLib_MakeWire.cxx:464-467`).
    pub fn last_edge(&self) -> Option<&Edge> {
        self.edge.as_ref()
    }

    /// `Vertex()` (`BRepLib_MakeWire.cxx:471-474`).
    pub fn last_vertex(&self) -> Option<&Vertex> {
        self.vertex.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
    use occt_geom::GeomLine;
    use std::sync::Arc;

    fn segment(b: &TopoBuilder, p1: &GpPnt, p2: &GpPnt, v1: &Vertex, v2: &Vertex) -> Edge {
        let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).unwrap();
        let mut e = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(*p1, dir))),
            0.0,
            p1.distance(p2),
        );
        b.add_edge_vertices(&mut e, v1, v2);
        e
    }

    /// Shared edges created in "foreign" directions (as a mesh's edge cache
    /// does) must be reoriented so that the wire chains head-to-tail. Like
    /// `mesh_to_brep::ensure_edge`, the three edges share the **same** vertex
    /// shapes.
    #[test]
    fn reorients_shared_edges_into_a_chain() {
        let b = TopoBuilder::new();
        let pa = GpPnt::new(0.0, 0.0, 0.0);
        let pb = GpPnt::new(1.0, 0.0, 0.0);
        let pc = GpPnt::new(0.0, 1.0, 0.0);
        let (va, vb, vc) = (
            b.make_vertex(pa, 0.0),
            b.make_vertex(pb, 0.0),
            b.make_vertex(pc, 0.0),
        );
        let e_ab = segment(&b, &pa, &pb, &va, &vb); // a->b (as needed)
        let e_cb = segment(&b, &pc, &pb, &vc, &vb); // stored c->b, needed b->c
        let e_ca = segment(&b, &pc, &pa, &vc, &va); // stored c->a (as needed)

        let mut mw = MakeWire::new();
        mw.add(&e_ab).unwrap();
        mw.add(&e_cb).unwrap();
        mw.add(&e_ca).unwrap();
        assert!(mw.is_done());
        assert_eq!(mw.error(), WireError::WireDone);

        let wire = mw.wire();
        let edges: Vec<Edge> = wire
            .0
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .filter(|s| s.shape_type() == ShapeType::Edge)
            .map(|s| Edge(s.clone()))
            .collect();
        assert_eq!(edges.len(), 3);
        let mut end: Option<Vertex> = None;
        for e in &edges {
            let (v1, v2) = edge_first_last_vertices(e);
            let (v1, v2) = (v1.unwrap(), v2.unwrap());
            if let Some(prev) = &end {
                assert!(
                    prev.0.same_tshape(&v1.0),
                    "wire is not chained: previous end != next start"
                );
            }
            end = Some(v2);
        }
        assert!(wire.closed(), "a triangle wire closes back onto its first vertex");
    }

    #[test]
    fn disconnected_edge_is_reported() {
        let b = TopoBuilder::new();
        let (p0, p1) = (GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0));
        let (p4, p5) = (GpPnt::new(5.0, 0.0, 0.0), GpPnt::new(6.0, 0.0, 0.0));
        let (v0, v1, v4, v5) = (
            b.make_vertex(p0, 0.0),
            b.make_vertex(p1, 0.0),
            b.make_vertex(p4, 0.0),
            b.make_vertex(p5, 0.0),
        );
        let e1 = segment(&b, &p0, &p1, &v0, &v1);
        let e2 = segment(&b, &p4, &p5, &v4, &v5);
        let mut mw = MakeWire::new();
        mw.add(&e1).unwrap();
        assert_eq!(mw.add(&e2), Err(WireError::DisconnectedWire));
        assert!(!mw.is_done());
    }
}
