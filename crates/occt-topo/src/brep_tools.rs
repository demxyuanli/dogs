//! `BRepTools` remainder: UV bounds, outer wire, triangulation dump.
//!
//! Source: `ModelingData/TKBRep/BRepTools`. [`add_uv_bounds`] delegates to
//! [`crate::brep_uv_bounds`] (including the periodic B-spline extra probe).
//! [`write_triangulation`] dumps the incremental-mesh triangles of a shape
//! (`BRepTools::Write` triangulation block).

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_geom::check_curve_on_surface::CheckCurveOnSurface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::brep_uv_bounds::{add_uv_bounds_edge, add_uv_bounds_face, add_uv_bounds_wire};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same, wires_of_face};

/// `BRepTools::OriEdgeInFace(E, F)` (`BRepTools.cxx:1346-1362`): the
/// orientation of the occurrence of `E` under `F` seen FORWARD. OCCT raises
/// `Standard_ConstructionError` when the edge is not in the face; that becomes
/// `Err`.
pub fn ori_edge_in_face(edge: &Edge, face: &Face) -> Result<Orientation, String> {
    let forward = face.0.oriented(Orientation::Forward);
    for e in explorer_shapes(&forward, ShapeType::Edge) {
        if is_same(&e, &edge.0) {
            return Ok(e.orientation());
        }
    }
    Err("Standard_ConstructionError: BRepTools::OriEdgeInFace".to_string())
}

/// `BRepTools::IsReallyClosed(E, F)` (`BRepTools.cxx:1204-1220`): `E` is
/// closed on `F` and occurs exactly twice in the face's edge exploration.
///
/// The occurrence count is `TopExp_Explorer(F, TopAbs_EDGE)` with
/// `IsSame`, so the two seam occurrences of one `TShape` count as two.
pub fn is_really_closed(edge: &Edge, face: &Face) -> bool {
    if !BRepTool::is_closed_edge_face(edge, face) {
        return false;
    }
    let nb_occ = explorer_shapes(&face.0, ShapeType::Edge)
        .iter()
        .filter(|e| is_same(e, &edge.0))
        .count();
    nb_occ == 2
}

/// `BRepTools::AddUVBounds(Face, Box2d)`.
pub fn add_uv_bounds(face: &Face, box2d: &mut BndBox2d) {
    add_uv_bounds_face(face, box2d);
}

/// `BRepTools::AddUVBounds(Face, Wire, Box2d)`.
pub fn add_uv_bounds_on_wire(face: &Face, wire: &Wire, box2d: &mut BndBox2d) {
    add_uv_bounds_wire(face, wire, box2d);
}

/// `BRepTools::AddUVBounds(Face, Edge, Box2d)`.
pub fn add_uv_bounds_on_edge(face: &Face, edge: &Edge, box2d: &mut BndBox2d) {
    add_uv_bounds_edge(face, edge, box2d);
}

/// `BRepTools::UVBounds(Face)` → `(u_min, u_max, v_min, v_max)`.
///
/// When the accumulated box is void OCCT sets all four bounds to zero
/// (`BRepTools.cxx:76-79`), not to the unit range.
pub fn uv_bounds(face: &Face) -> (f64, f64, f64, f64) {
    let mut b = BndBox2d::new();
    add_uv_bounds_face(face, &mut b);
    match b.get() {
        Some((u0, v0, u1, v1)) => (u0, u1, v0, v1),
        None => (0.0, 0.0, 0.0, 0.0),
    }
}

/// `BRepTools::DetectClosedness(F, uclosed, vclosed)`
/// (`BRepTools.cxx:1224-1251`): which of the two surface directions the face
/// closes on, from its seam edges.
///
/// For every edge that is closed on `F` *and* really closed (`IsReallyClosed`),
/// the two pcurves are evaluated at the shared `FirstParameter` of the
/// `CurveOnSurface` representation: the seam runs along U when the two points
/// differ more in X than in Y (`cxx:1240`), along V otherwise.
///
/// Returns `(u_closed, v_closed)`.
///
/// `BRep_Tool::CurveOnSurface(E, F, f, l)` answers the forward pcurve and the
/// reversed-edge overload the reversed one; the port reads both slots of the
/// registry and keeps OCCT's single shared `fpar` (both COS ranges are the
/// edge range, `BRep_Builder::UpdateEdge(E, C1, C2, F)`).
pub fn detect_closedness(face: &Face) -> (bool, bool) {
    let mut u_closed = false;
    let mut v_closed = false;
    let face_key = GeometryRegistry::shape_key(&face.0);
    for shape in explorer_shapes(&face.0, ShapeType::Edge) {
        let edge = Edge(shape);
        if !is_really_closed(&edge, face) {
            continue;
        }
        let pcurves = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
        // Port guard: OCCT dereferences both `PCurve1` and `PCurve2` at
        // `cxx:1238-1239`; a closed edge without the second representation
        // faults there.
        if pcurves.len() < 2 {
            continue;
        }
        // `CurveOnSurface` writes the COS range into `fpar`/`lpar`; the second
        // call overwrites them with the same values, so `Point1` and `Point2`
        // are both read at that single `fpar`.
        let fpar = GeometryRegistry::global()
            .pcurve_range(&edge.0, face_key)
            .map(|(f, _)| f)
            .unwrap_or_else(|| BRepTool::edge_parameters(&edge).0);
        let point1 = pcurves[0].d0(fpar);
        let point2 = pcurves[1].d0(fpar);
        // `cxx:1240`.
        let is_u_iso = (point1.x() - point2.x()).abs() > (point1.y() - point2.y()).abs();
        if is_u_iso {
            u_closed = true;
        } else {
            v_closed = true;
        }
    }
    (u_closed, v_closed)
}

/// `BRepTools::OuterWire`. The wire whose UV box contains every other wire
/// of the face; if none, the first wire.
pub fn outer_wire(face: &Face) -> Option<Wire> {
    let wires = wires_of_face(face);
    if wires.is_empty() {
        return None;
    }
    let mut best = 0usize;
    let mut best_span = -1.0;
    for (i, w) in wires.iter().enumerate() {
        let mut b = BndBox2d::new();
        add_uv_bounds_wire(face, w, &mut b);
        let span = match b.get() {
            Some((u0, v0, u1, v1)) => (u1 - u0).abs() * (v1 - v0).abs(),
            None => 0.0,
        };
        if span > best_span {
            best_span = span;
            best = i;
        }
    }
    Some(wires[best].clone())
}

/// `BRepTools::Map3DEdges` — unique 3-D edges of `shape`.
pub fn map_3d_edges(shape: &TopoShape) -> Vec<Edge> {
    crate::topo_tools::edges_of(shape)
}

/// Dump the incremental-mesh triangulation of `shape` as ASCII
/// (`BRepTools` triangulation write). Each triangle is `v0 v1 v2` after a
/// vertex list. `deflection` is the mesher linear deflection.
pub fn write_triangulation(shape: &TopoShape, deflection: f64) -> String {
    // One mesher only (`BRepMesh_IncrementalMesh`, via `meshing::incremental_mesh`).
    // A shape whose triangulation fails contributes no triangles; OCCT has no
    // second mesher to fall back to (audit A17/A18/A19).
    let mesh = crate::meshing::incremental_mesh::incremental_mesh_to_shape_mesh(shape, deflection)
        .unwrap_or_else(|_| crate::mesh::ShapeMesh {
            vertices: Vec::new(),
            triangles: Vec::new(),
            source_shape: crate::abs::ShapeType::Shape,
        });
    let mut out = String::new();
    out.push_str(&format!("Triangulation {}\n", mesh.vertices.len()));
    for p in &mesh.vertices {
        out.push_str(&format!("{} {} {}\n", p.x(), p.y(), p.z()));
    }
    out.push_str(&format!("Triangles {}\n", mesh.triangles.len()));
    for t in &mesh.triangles {
        out.push_str(&format!("{} {} {}\n", t.n0, t.n1, t.n2));
    }
    out
}

/// `BRepTools::Compare` vertices: same point within the larger tolerance.
pub fn compare_vertices(a: &crate::shape::Vertex, b: &crate::shape::Vertex) -> bool {
    let pa = BRepTool::vertex_point(a);
    let pb = BRepTool::vertex_point(b);
    let tol = BRepTool::vertex_tolerance(a).max(BRepTool::vertex_tolerance(b));
    pa.distance(&pb) <= tol
}

/// `BRepTools::EvalAndUpdateTol(E, C3d, C2d, S, f, l)`
/// (`BRepTools.cxx:1255-1342`): the max distance between the 3D curve and the
/// pcurve lifted to `S` over `[f, l]`, written back onto the edge when it
/// exceeds the stored tolerance.
///
/// `first` / `last` are clamped to the curves' own ranges unless the curve is
/// periodic (`cxx:1266-1278`). `GeomLib_CheckCurveOnSurface`
/// (`cxx:1281-1289`) is [`CheckCurveOnSurface`] over the
/// `Adaptor3d_CurveOnSurface` of the pcurve; the `ErrorStatus() == 3` /
/// `== 2 with a periodic input` fallback (`cxx:1292-1330`) samples 22
/// intervals instead and answers `1.2 * sqrt(dmax2)`.
pub fn eval_and_update_tol(
    e: &Edge,
    c3d: &Arc<dyn occt_geom::Curve>,
    c2d: &Arc<dyn Curve2d>,
    s: &Arc<dyn occt_geom::Surface>,
    f: f64,
    l: f64,
) -> f64 {
    let mut newtol = 0.0_f64;
    let mut first = f;
    let mut last = l;
    // `cxx:1266-1278`.
    if !c3d.is_periodic() {
        first = first.max(c3d.first_parameter());
        last = last.min(c3d.last_parameter());
    }
    if !c2d.is_periodic() {
        first = first.max(c2d.first_parameter());
        last = last.min(c2d.last_parameter());
    }

    // `cxx:1281-1289`: `GeomLib_CheckCurveOnSurface` over the curve-on-surface.
    let cos = occt_geom::curve_on_surface::CurveOnSurface::new(c2d.clone(), s.clone());
    let mut check = CheckCurveOnSurface::with_curve(c3d.as_ref(), occt_core::precision::PCONFUSION);
    check.perform(&cos, c2d.as_ref());
    if check.is_done() {
        newtol = check.max_distance();
    } else if check.error_status() == 3
        || (check.error_status() == 2 && (c3d.is_periodic() || c2d.is_periodic()))
    {
        // `cxx:1294-1330`: estimate by sample points.
        let nbint = 22;
        let mut dt = (last - first) / nbint as f64;
        dt = dt.max(occt_core::precision::CONFUSION);
        let mut dmax = 0.0_f64;
        let mut cnt = 0;
        let mut t = first;
        while t <= last {
            cnt += 1;
            let p2d = c2d.d0(t);
            let pc = c3d.d0(t);
            let ps = s.d0(p2d.x(), p2d.y());
            let d = ps.square_distance(&pc);
            if d > dmax {
                dmax = d;
            }
            t += dt;
        }
        if cnt < nbint + 1 {
            let t = last;
            let p2d = c2d.d0(t);
            let pc = c3d.d0(t);
            let ps = s.d0(p2d.x(), p2d.y());
            let d = ps.square_distance(&pc);
            if d > dmax {
                dmax = d;
            }
        }
        newtol = 1.2 * dmax.sqrt();
    }

    // `cxx:1334-1339`.
    let mut tol = BRepTool::edge_tolerance(e);
    if newtol > tol {
        tol = newtol;
        GeometryRegistry::global().set_edge_tolerance(&e.0, tol);
    }
    tol
}
