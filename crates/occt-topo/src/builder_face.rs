//! Face reconstruction — a Rust port of `BOPAlgo_BuilderFace`.
//!
//! Rebuilds a closed planar `Face` from an unordered set of boundary `Edge`s:
//!
//! * [`make_face_from_wire`] orders the edges into a single closed loop and
//!   trims a (given or plane-derived) surface to it.
//! * [`build_face_with_holes`] builds a face with one outer loop and any
//!   number of inner hole loops, checking that the hole rings wind opposite
//!   the outer ring.
//! * [`FaceBuilder`] implements the [`AreaBuilder`] pipeline over a generatrix
//!   face (`BOPAlgo_BuilderFace::SetFace`).
//!
//! The loop-closing is the vertex-adjacency cycle tracer from
//! [`crate::builder_area`]; no `BOPAlgo_WireSplitter` dependency is pulled in.
//!
//! Orientation note: this port's flat `TShape` child tree does not retain the
//! per-wire orientation of an edge (see `TopoBuilder::add`), so a wire is a
//! multiset of edges. The *geometric* winding of a loop is therefore computed
//! separately by [`loop_signed_area`], which walks the true 3-D boundary and
//! is what the hole-orientation check uses.

use std::sync::Arc;

use occt_core::gp::{GpPln, GpVec};
use occt_geom::{GeomPlane, Surface};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::builder_area::{
    AreaBuilder, AreaBuilderBase, build_loops, edges_form_closed_loop, make_wire_from_edges,
    plane_from_loop,
};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, vertex_position};

/// A face-oriented area builder (`BOPAlgo_BuilderFace`).
///
/// Sets a generatrix face whose surface is reused for the area faces. The
/// four-phase pipeline is the `AreaBuilder` one; `perform_areas` builds each
/// loop into a face on the generatrix surface (or a plane derived from the
/// loop when no face is set).
#[derive(Debug)]
pub struct FaceBuilder {
    /// Shared area-builder state.
    pub base: AreaBuilderBase,
    /// The generatrix face (`BOPAlgo_BuilderFace::myFace`).
    pub face: Option<Face>,
    /// Orientation of the generatrix face preserved for the result.
    pub orientation: Orientation,
}

impl Default for FaceBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl FaceBuilder {
    pub fn new() -> Self {
        Self {
            base: AreaBuilderBase::default(),
            face: None,
            orientation: Orientation::External,
        }
    }

    /// Set the generatrix face (`BOPAlgo_BuilderFace::SetFace`). Its
    /// orientation is recorded and `myFace` is stored `FORWARD`.
    pub fn set_face(&mut self, f: &Face) {
        self.orientation = f.0.orientation();
        let mut ff = f.clone();
        ff.0.set_orientation(Orientation::Forward);
        self.face = Some(ff);
    }

    /// The generatrix face, if set.
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }
}

impl AreaBuilder for FaceBuilder {
    fn set_shapes(&mut self, shapes: &[TopoShape]) {
        self.base.set_shapes(shapes);
    }

    fn perform_shapes_to_avoid(&mut self) -> Result<(), String> {
        crate::builder_face_occt::perform_shapes_to_avoid(self)
    }

    fn perform_loops(&mut self) -> Result<(), String> {
        crate::builder_face_occt::perform_loops(self)
    }

    fn perform_areas(&mut self) -> Result<(), String> {
        if self.face.is_some() {
            return crate::builder_face_occt::perform_areas(self);
        }
        self.base.areas.clear();
        let loops = self.base.loops.clone();
        let b = TopoBuilder::new();
        let base_surface = self.face.as_ref().and_then(|f| BRepTool::face_surface(f));
        for l in loops {
            let wire = Wire(l);
            let edges = edges_of_wire(&wire);
            let surf: Arc<dyn Surface> = match &base_surface {
                Some(s) => s.clone(),
                None => Arc::new(GeomPlane::new(plane_from_loop(&edges)?)),
            };
            let mut face = b.make_face(surf, &[wire]);
            if let Some(gf) = &self.face {
                face.0.set_orientation(gf.0.orientation());
            }
            self.base.areas.push(face.0);
        }
        Ok(())
    }

    fn perform_internal_shapes(&mut self) -> Result<(), String> {
        crate::builder_face_occt::perform_internal_shapes(self)
    }

    fn areas(&self) -> &[TopoShape] {
        &self.base.areas
    }
}

/// Rebuild a closed face from an unordered set of boundary edges.
///
/// The edges are ordered into a single closed loop by vertex adjacency
/// ([`build_loops`]); when `surface` is `None` a plane through the loop
/// vertices is derived. Mirrors `BRepBuilderAPI_MakeFace(Wire)` after the
/// `BOPAlgo_BuilderFace` loop construction.
pub fn make_face_from_wire(edges: &[Edge], surface: Option<Arc<dyn Surface>>) -> Result<Face, String> {
    if edges.len() < 3 {
        return Err("make_face_from_wire: need at least 3 edges".into());
    }
    let loops = build_loops(edges)?;
    if loops.len() != 1 {
        return Err(format!(
            "make_face_from_wire: expected a single closed loop, found {}",
            loops.len()
        ));
    }
    let loop_edges = &loops[0];
    if !edges_form_closed_loop(loop_edges) {
        return Err("make_face_from_wire: edges do not form a closed loop".into());
    }
    let wire = make_wire_from_edges(loop_edges);
    let surf: Arc<dyn Surface> = match surface {
        Some(s) => s,
        None => Arc::new(GeomPlane::new(plane_from_loop(loop_edges)?)),
    };
    let face = TopoBuilder::new().make_face(surf, &[wire]);
    Ok(face)
}

/// Build a face from one outer loop and any number of inner hole loops
/// (`BRepBuilderAPI_MakeFace` with multiple wires).
///
/// Each loop is closed by vertex adjacency. The outer loop provides the
/// supporting plane; every hole ring must be *closed* and wind **opposite** to
/// the outer ring (a hole has negative signed area relative to the outer on
/// the same plane). A hole that winds the same way is reported as an error —
/// this is the "loop orientation check". The returned face carries one wire
/// per loop.
pub fn build_face_with_holes(outer: &[Edge], holes: &[Vec<Edge>]) -> Result<Face, String> {
    let outer_loop = single_loop(outer, "outer")?;
    if outer_loop.len() < 3 {
        return Err("build_face_with_holes: outer loop has fewer than 3 edges".into());
    }
    let pln = plane_from_loop(&outer_loop)?;
    let outer_signed = loop_signed_area(&outer_loop, &pln);
    if outer_signed.abs() < 1e-9 {
        return Err("build_face_with_holes: degenerate outer loop".into());
    }

    let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
    let b = TopoBuilder::new();
    let mut wires = vec![make_wire_from_edges(&outer_loop)];

    for (hi, h) in holes.iter().enumerate() {
        let hole_loop = single_loop(h, &format!("hole {hi}"))?;
        let hole_signed = loop_signed_area(&hole_loop, &pln);
        if hole_signed.abs() < 1e-9 {
            return Err(format!("build_face_with_holes: degenerate hole {hi} loop"));
        }
        // Orientation check: the hole ring must wind opposite the outer ring.
        if outer_signed * hole_signed > 0.0 {
            return Err(format!(
                "build_face_with_holes: hole {hi} winds the same way as the outer loop"
            ));
        }
        wires.push(make_wire_from_edges(&hole_loop));
    }

    let face = b.make_face(surf, &wires);
    Ok(face)
}

/// Order an edge set into a single closed loop, erroring otherwise.
fn single_loop(edges: &[Edge], name: &str) -> Result<Vec<Edge>, String> {
    if edges.len() < 3 {
        return Err(format!("build_face_with_holes: {name} loop has fewer than 3 edges"));
    }
    let loops = build_loops(edges)?;
    if loops.len() != 1 {
        return Err(format!(
            "build_face_with_holes: {name} edges form {} loops, expected 1",
            loops.len()
        ));
    }
    Ok(loops.into_iter().next().expect("one loop"))
}

/// Signed area of a planar loop traced by following the edges, projected onto
/// `pln`.
///
/// The walk starts at the first edge's start vertex and, for every following
/// edge, continues from whichever of its two endpoints matches the current
/// point — mentally reversing edges as needed. This yields the *true*
/// geometric winding of the loop, independent of the flat wire's inability to
/// store per-edge orientation. Positive for CCW, negative for CW (with respect
/// to the plane's `(u, v)` frame). Returns `0.0` when the edges do not chain.
pub fn loop_signed_area(edges: &[Edge], pln: &GpPln) -> f64 {
    let Some(first) = edges.first() else { return 0.0 };
    let (a0, b0) = edge_vertices(first);
    let (Some(va0), Some(vb0)) = (a0, b0) else { return 0.0 };
    let start = vertex_position(&va0);
    let mut cur = vertex_position(&vb0);
    let mut pts = vec![start, cur];
    for e in &edges[1..] {
        let (a, b) = edge_vertices(e);
        let (Some(va), Some(vb)) = (a, b) else { return 0.0 };
        let (pa, pb) = (vertex_position(&va), vertex_position(&vb));
        if pa.distance(&cur) < 1e-9 {
            cur = pb;
        } else if pb.distance(&cur) < 1e-9 {
            cur = pa;
        } else {
            return 0.0; // not a chain
        }
        pts.push(cur);
    }

    let pos = pln.position();
    let xd = *pos.x_direction().xyz();
    let yd = *pos.y_direction().xyz();
    let loc = pos.location();
    let n = pts.len();
    let mut acc = 0.0;
    for i in 0..n {
        let d1 = GpVec::from_pnts(&loc, &pts[i]);
        let d2 = GpVec::from_pnts(&loc, &pts[(i + 1) % n]);
        let (u1, v1) = (d1.dot(&GpVec::from_xyz(&xd)), d1.dot(&GpVec::from_xyz(&yd)));
        let (u2, v2) = (d2.dot(&GpVec::from_xyz(&xd)), d2.dot(&GpVec::from_xyz(&yd)));
        acc += u1 * v2 - u2 * v1;
    }
    0.5 * acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_faces::face_area_exact;
    use crate::builder::TopoBuilder;
    use crate::builder_area::build_loops;
    use crate::shape::TopoShape;
    use crate::topo_tools_full::edges_of;
    use occt_core::gp::GpPnt;

    /// Four segment edges of a polygon; each edge carries its own vertex
    /// children (positions still chain by [`vertex_key`]).
    fn square_edges(b: &TopoBuilder, pts: &[GpPnt; 4]) -> Vec<Edge> {
        (0..4).map(|i| b.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect()
    }

    #[test]
    fn box_faces_rebuild_closed_with_area_one() {
        let ub = unit_box();
        for (i, face) in ub.faces.iter().enumerate() {
            let edges = edges_of(&face.0);
            assert_eq!(edges.len(), 4, "face {i} has 4 distinct boundary edges");
            let rebuilt = make_face_from_wire(&edges, None).expect("rebuild box face");
            // closed: a single closed loop through the 4 edges.
            let loops = build_loops(&edges).unwrap();
            assert_eq!(loops.len(), 1, "face {i} edges form one loop");
            assert!(edges_form_closed_loop(&loops[0]), "face {i} loop is closed");
            // one wire with 4 edges.
            let wires = crate::topo_tools_full::wires_of_face(&rebuilt);
            assert_eq!(wires.len(), 1);
            assert_eq!(crate::topo_tools_full::edges_of_wire(&wires[0]).len(), 4);
            // area 1.0 for a unit quad.
            let area = face_area_exact(&rebuilt, 16, 16);
            assert!((area - 1.0).abs() < 1e-6, "face {i} area = {area}");
        }
    }

    #[test]
    fn face_builder_pipeline_rebuilds_box_face() {
        let ub = unit_box();
        let face = &ub.faces[0];
        let shapes: Vec<TopoShape> = edges_of(&face.0).into_iter().map(|e| e.0.clone()).collect();
        let mut fb = FaceBuilder::new();
        fb.set_face(face);
        fb.set_shapes(&shapes);
        fb.perform().expect("face pipeline");
        assert_eq!(fb.areas().len(), 1);
        let area_face = Face::wrap(fb.areas()[0].clone()).expect("area is a face");
        let area = face_area_exact(&area_face, 16, 16);
        assert!((area - 1.0).abs() < 1e-6, "rebuilt area = {area}");
    }

    #[test]
    fn face_with_hole_area_is_outer_minus_hole() {
        let b = TopoBuilder::new();
        // Outer unit square, CCW when viewed from +z.
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let outer_edges = square_edges(&b, &outer);
        // Hole quarter square, CW (right → down → left → up).
        let hole = [
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let hole_edges = square_edges(&b, &hole);

        let face = build_face_with_holes(&outer_edges, &[hole_edges.clone()])
            .expect("face with a square hole");
        assert_eq!(crate::topo_tools_full::wires_of_face(&face).len(), 2);

        let pln = plane_from_loop(&build_loops(&outer_edges).unwrap()[0]).unwrap();
        let outer_loop = &build_loops(&outer_edges).unwrap()[0];
        let hole_loop = &build_loops(&hole_edges).unwrap()[0];
        let a_outer = loop_signed_area(outer_loop, &pln).abs();
        let a_hole = loop_signed_area(hole_loop, &pln).abs();
        assert!((a_outer - 1.0).abs() < 1e-9, "outer area {a_outer}");
        assert!((a_hole - 0.25).abs() < 1e-9, "hole area {a_hole}");
        assert!((a_outer - a_hole - 0.75).abs() < 1e-9, "face area = outer − hole");
    }

    #[test]
    fn hole_winding_same_as_outer_is_rejected() {
        let b = TopoBuilder::new();
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let outer_edges = square_edges(&b, &outer);
        // Hole wound CCW — same as the outer ring → must be rejected.
        let hole = [
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let hole_edges = square_edges(&b, &hole);
        let r = build_face_with_holes(&outer_edges, &[hole_edges]);
        assert!(r.is_err(), "same-winding hole must be rejected");
    }

    #[test]
    fn degenerate_open_edges_rejected() {
        let b = TopoBuilder::new();
        // Two segments forming an open chain cannot make a wire face.
        let e1 = b.make_edge_segment(&GpPnt::new(0., 0., 0.), &GpPnt::new(1., 0., 0.));
        let e2 = b.make_edge_segment(&GpPnt::new(1., 0., 0.), &GpPnt::new(1., 1., 0.));
        assert!(make_face_from_wire(&[e1, e2], None).is_err());
    }
}
