//! `UpdateTolerance` (`BRepOffset_MakeOffset.cxx:4196-4302`) together with its
//! file-local helper `ComputeMaxDist` (`:4165-4192`).
//!
//! `UpdateTolerance` raises the tolerance of the result edges that the offset
//! pipeline did not re-derive (`:4240-4245` skips the caps and the faces of the
//! initial shape): on a planar face every edge gets at least the largest
//! distance between its 3D curve and the plane, sampled at `NCONTROL = 23`
//! parameters (`:4252-4263`), and every edge not seen before gets at least
//! `BRepCheck_Edge::Tolerance()` (`:4265-4275`). An updated edge passes its new
//! tolerance to both of its vertices (`:4278-4298`).
//!
//! UNPORTED / port notes:
//!
//! * `TopoDS_Shape::Locked(false)` (`:4267`, `:4285`) - the port has no locked
//!   flag on `TopoDS_TShape`, so the two calls are dropped.
//! * `BRepCheck_Vertex::Tolerance()` (`BRepCheck_Vertex.cxx:343-384`) walks the
//!   point representations of `BRep_TVertex`. The port models them
//!   ([`crate::brep_point_representation`]) but
//!   [`vertex_check_tolerance`] does not walk them, so the answer is
//!   `sqrt(Tol * 1.05)` with `Tol` left at the value the caller just set
//!   (`0`, `:4291`). The following `B.UpdateVertex(V, Tol)` (`:4297`)
//!   overwrites it with the edge tolerance in both implementations.
//! * `TV->ChangePoints().Clear()` (`:4293`) is not ported: the vertex keeps its
//!   point representations.
//! * `BRep_Tool::Curve(E, f, l)` (`:4258`) answers a null handle for an edge
//!   without a 3D curve; `ComputeMaxDist` would dereference it. The port skips
//!   the edge instead of dereferencing.
//! * `theInitShape.IsNull()` (`:4219`) is `Option::is_none` here.

use std::collections::HashSet;

use occt_core::gp::GpPln;
use occt_core::precision::{Precision, CONFUSION, INFINITE};
use occt_geom::curve::Curve;

use crate::abs::ShapeType;
use crate::brep_check;
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};

use crate::brep_offset::inter3d::ShapeIndexedMap;

/// `NCONTROL` (`cxx:4172`).
const NCONTROL: usize = 23;

/// `ComputeMaxDist(thePlane, theCrv, theFirst, theLast)` (`cxx:4165-4192`):
/// the largest distance from `theCrv` to `thePlane` over `NCONTROL` samples,
/// with the 5 % margin of `cxx:4191`.
pub(crate) fn compute_max_dist(
    the_plane: &GpPln,
    the_crv: &dyn Curve,
    the_first: f64,
    the_last: f64,
) -> f64 {
    let mut a_max_dist = 0.0f64;
    for i in 0..NCONTROL {
        // `cxx:4177`.
        let a_prm = ((NCONTROL - 1 - i) as f64 * the_first + i as f64 * the_last)
            / (NCONTROL - 1) as f64;
        let a_p = the_crv.d0(a_prm);
        // `cxx:4180-4184`.
        if Precision::is_infinite(a_p.x())
            || Precision::is_infinite(a_p.y())
            || Precision::is_infinite(a_p.z())
        {
            return INFINITE;
        }
        // `cxx:4185-4189`.
        let a_dist2 = the_plane.square_distance_pnt(&a_p);
        if a_dist2 > a_max_dist {
            a_max_dist = a_dist2;
        }
    }
    // `cxx:4191`.
    a_max_dist.sqrt() * 1.05
}

/// `BRepCheck_Vertex::Tolerance()` (`BRepCheck_Vertex.cxx:343-384`) for a
/// vertex whose tolerance has just been reset to `0` (`cxx:4291`).
///
/// `Tol` starts at `BRep_Tool::Tolerance(V)^2`, every point representation of
/// `BRep_TVertex::Points()` raises it, and the answer is `sqrt(Tol * 1.05)`
/// (`cxx:384`). The port models the representations
/// ([`crate::brep_point_representation`]) but this helper does not walk them,
/// so only the base term survives. The value is immediately overwritten by
/// `B.UpdateVertex(V, Tol)` with the edge tolerance (`:4297`) in both
/// implementations, so it does not reach the result.
fn vertex_check_tolerance(v: &Vertex) -> f64 {
    let tol = BRepTool::vertex_tolerance(v).max(CONFUSION);
    (tol * tol * 1.05).sqrt()
}

/// `UpdateTolerance(S, Faces, theInitShape)` (`cxx:4196-4302`).
///
/// `init_shape` is `myShape` for the thickening path and a null shape
/// otherwise (`:1041-1049`).
pub fn update_tolerance(s: &TopoShape, faces: &ShapeIndexedMap, init_shape: Option<&TopoShape>) {
    let reg = GeometryRegistry::global();

    // `cxx:4208-4216`: the edges of the caps are not modified.
    let mut view: HashSet<ShapeMapKey> = HashSet::new();
    for j in 1..=faces.extent() {
        let Some(f) = faces.value(j) else {
            continue;
        };
        for e in explorer_shapes(f, ShapeType::Edge) {
            view.insert(shape_map_key(&e));
        }
    }

    // `cxx:4218-4234`: neither are the faces, edges and vertices of the initial
    // shape.
    let mut a_map_init_f: HashSet<ShapeMapKey> = HashSet::new();
    if let Some(init_shape) = init_shape {
        for f in explorer_shapes(init_shape, ShapeType::Face) {
            a_map_init_f.insert(shape_map_key(&f));
            for e in explorer_shapes(&f, ShapeType::Edge) {
                view.insert(shape_map_key(&e));
                for v in explorer_shapes(&e, ShapeType::Vertex) {
                    view.insert(shape_map_key(&v));
                }
            }
        }
    }

    // `cxx:4237-4299`.
    for f_shape in explorer_shapes(s, ShapeType::Face) {
        if faces.contains(&f_shape) || a_map_init_f.contains(&shape_map_key(&f_shape)) {
            continue;
        }
        let Some(face) = Face::wrap(f_shape) else {
            continue;
        };
        // `cxx:4246`: `BRepAdaptor_Surface aBAS(TopoDS::Face(F), false)`.
        let Some(surf) = BRepTool::face_surface(&face) else {
            continue;
        };
        let is_plane = surface_kind(surf.as_ref()) == GeomAbsSurfaceType::Plane;

        for e_shape in explorer_shapes(&face.0, ShapeType::Edge) {
            let Some(e) = Edge::wrap(e_shape) else {
                continue;
            };
            let mut is_updated = false;
            let a_curr_tol = BRepTool::edge_tolerance(&e);

            // `cxx:4252-4263`: an edge without a pcurve on the plane does not
            // reach `EdgeCorrector`, so the plane distance is checked here.
            if is_plane {
                if let Some(a_crv) = BRepTool::edge_curve(&e) {
                    let (a_first, a_last) = BRepTool::edge_parameters(&e);
                    let a_max_dist = compute_max_dist(
                        &surf.gp_pln().expect("plane surface"),
                        a_crv.as_ref(),
                        a_first,
                        a_last,
                    );
                    if a_max_dist > a_curr_tol {
                        reg.set_edge_tolerance(&e.0, a_max_dist);
                        is_updated = true;
                    }
                }
            }

            // `cxx:4265-4276`.
            if view.insert(shape_map_key(&e)) {
                let tol = brep_check::edge_tolerance(&e);
                if tol > a_curr_tol {
                    reg.set_edge_tolerance(&e.0, tol);
                    is_updated = true;
                }
            }

            // `cxx:4278-4298`: update the vertices.
            if is_updated {
                let tol = BRepTool::edge_tolerance(&e);
                let (v1, v2) = edge_vertices(&e);
                for v in [v1, v2].into_iter().flatten() {
                    if view.insert(shape_map_key(&v.0)) {
                        v.set_tolerance(0.0);
                        let tol_v = vertex_check_tolerance(&v);
                        v.set_tolerance(tol_v);
                    }
                    v.set_tolerance(tol);
                }
            }
        }
    }
}
