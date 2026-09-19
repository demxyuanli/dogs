//! Port of `BRepBndLib::Add(const TopoDS_Shape& S, Bnd_Box& B, bool useTriangulation)`
//! (`BRepBndLib.cxx:81-215`).
//!
//! The `useTriangulation == false` control flow is the one the mesh pipeline
//! needs (`Prs3d::GetDeflection`, `Prs3d.hxx:66-72`): every face whose surface
//! is not a plane contributes `BndLib_AddSurface::Add` (restricted to the
//! face's `BRepTools::UVBounds`), a plane face contributes the
//! `BndLib_Add3dCurve::Add` bounds of its boundary edges plus the face
//! tolerance, and every edge not owned by a face contributes
//! `BndLib_Add3dCurve::Add`. The box is built from the control net / analytic
//! extrema OCCT uses, not from a sampled envelope.
//!
//! PARK (no `Poly_Triangulation` / `Poly_Polygon3D` in this port's BRep
//! registry; reported, not fudged):
//! - the triangulation arm of the face loop (`BRepBndLib.cxx:95-101`) and of
//!   the free-edge loop (`:157-179`), both reachable only when the shape
//!   already carries a triangulation. A face with a surface never reaches
//!   them here, so the geometry arms below are the whole `false` path.
//! - edges that carry only a pcurve (`BRepAdaptor_Curve::Initialize` fallback,
//!   `BRepAdaptor_Curve.cxx:93-106`): the registry stores 3d curves only, so
//!   such an edge is skipped.
//! - `TopoDS_Compound` / `CompSolid` need no special case in `Add`: the shape
//!   is traversed by `TopExp_Explorer`, and this port walks the same children
//!   tree. Vertices are never added by `Add`.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::bnd::BndBox;

use crate::brep_tool::BRepTool;
use crate::brep_uv_bounds::uv_box_of_face;
use crate::geom_bnd_lib_curve3d::box_curve;
use crate::geom_bnd_lib_surface3d::{add_surface, surface_kind, SurfaceKind};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools_full::{edges_of, faces_of, vertices_of};

/// `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)` (`BRepTools.cxx:110-125`):
/// scalar form of `BRepTools::AddUVBounds`, with the void/whole fallbacks.
fn uv_bounds(f: &Face) -> (f64, f64, f64, f64) {
    let b = uv_box_of_face(f);
    if b.is_void() {
        (1.0, -1.0, 1.0, -1.0)
    } else if b.is_whole() {
        (-1.0, 1.0, -1.0, 1.0)
    } else {
        match b.get() {
            Some((x, y, xx, yy)) => (x, xx, y, yy),
            None => (1.0, -1.0, 1.0, -1.0),
        }
    }
}

/// `BRepAdaptor_Curve::Initialize(E)` + `BndLib_Add3dCurve::Add(BC, Tol, B)`
/// (`BRepBndLib.cxx:126-129`, `:148-150`).
fn add_edge_curve(e: &Edge, b: &mut BndBox) {
    // `BRep_Tool::IsGeometric(edge)` (`BRep_Tool.cxx:238-262`): a 3d curve or a
    // curve on surface. PARK: a pcurve-only edge is skipped (see module doc).
    let Some(curve) = BRepTool::edge_curve_world(e) else {
        return;
    };
    let (u1, u2) = BRepTool::edge_parameters(e);
    let tol = BRepTool::edge_tolerance(e);
    b.add_box(&box_curve(curve.as_ref(), u1, u2, tol));
}

/// `TopExp_Explorer(S, TopAbs_EDGE, TopAbs_FACE)` (`BRepBndLib.cxx:142`):
/// the edges of `S` that are not sub-shapes of any face of `S`.
fn edges_not_in_faces(s: &TopoShape) -> Vec<Edge> {
    let mut in_faces: HashSet<usize> = HashSet::new();
    for f in faces_of(s) {
        for e in edges_of(&f.0) {
            in_faces.insert(Arc::as_ptr(&e.0.tshape) as usize);
        }
    }
    edges_of(s)
        .into_iter()
        .filter(|e| !in_faces.contains(&(Arc::as_ptr(&e.0.tshape) as usize)))
        .collect()
}

/// `TopExp_Explorer(S, TopAbs_VERTEX, TopAbs_EDGE)` (`BRepBndLib.cxx:208`):
/// the vertices of `S` that are not sub-shapes of any edge of `S`.
fn vertices_not_in_edges(s: &TopoShape) -> Vec<Vertex> {
    let mut in_edges: HashSet<usize> = HashSet::new();
    for e in edges_of(s) {
        for v in vertices_of(&e.0) {
            in_edges.insert(Arc::as_ptr(&v.0.tshape) as usize);
        }
    }
    vertices_of(s)
        .into_iter()
        .filter(|v| !in_edges.contains(&(Arc::as_ptr(&v.0.tshape) as usize)))
        .collect()
}

/// `BRepBndLib::Add(S, B, useTriangulation)` (`BRepBndLib.cxx:81-215`).
pub fn add(shape: &TopoShape, b: &mut BndBox, use_triangulation: bool) {
    // PARK: `BRep_Tool::Triangulation` / `BRep_Tool::Polygon3D` are always null
    // in this port, so both triangulation arms are unreachable. The flag is
    // accepted for call parity with `BRepBndLib::Add(S, B, false)`.
    let _ = use_triangulation;

    // Add the faces.
    for f in faces_of(shape) {
        let Some(gs) = BRepTool::face_surface_world(&f) else {
            // `GS.IsNull()` and no triangulation: the face contributes nothing.
            continue;
        };
        if surface_kind(gs.as_ref()) != SurfaceKind::Plane {
            let (u1, u2, v1, v2) = uv_bounds(&f);
            let tol = BRepTool::face_tolerance(&f);
            add_surface(&gs, u1, u2, v1, v2, tol, b);
        } else {
            // Plane: work directly on the 3d curves of the boundary edges.
            let face_edges = edges_of(&f.0);
            if face_edges.is_empty() {
                let (u1, u2, v1, v2) = uv_bounds(&f);
                let tol = BRepTool::face_tolerance(&f);
                add_surface(&gs, u1, u2, v1, v2, tol, b);
            } else {
                for e in &face_edges {
                    add_edge_curve(e, b);
                }
                let ft = BRepTool::face_tolerance(&f);
                b.enlarge(ft);
            }
        }
    }

    // Add the edges not in faces.
    for e in edges_not_in_faces(shape) {
        add_edge_curve(&e, b);
    }

    // Add the vertices not in edges (`BRepBndLib.cxx:206-213`).
    for v in vertices_not_in_edges(shape) {
        b.add_point(&BRepTool::vertex_point_world(&v));
        b.enlarge(BRepTool::vertex_tolerance(&v));
    }
}

/// `BRepBndLib::Add(S, B, false)` + `Bnd_Box::Gap` view, mirroring the
/// `shape_bbox` entry point used by the mesh pipeline.
pub fn shape_bnd_box(shape: &TopoShape) -> BndBox {
    let mut b = BndBox::new();
    add(shape, &mut b, false);
    b
}
