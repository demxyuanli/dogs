//! `BRepOffset_MakeOffset::IsPlanar` (`BRepOffset_MakeOffset.cxx:5447-5552`)
//! and the file-local `RemoveSeamAndDegeneratedEdges` (`:5556-5633`) it uses.
//!
//! `IsPlanar` decides whether the whole shape can take the planar-intersection
//! short cut (`myInter && !isArc && myIsPlanar` in `analyzeProgress`,
//! `cxx:5421-5422`). Every face that is not a `GeomAbs_Plane` but whose surface
//! `GeomLib_IsPlanarSurface` accepts at `Precision::Confusion()` is replaced by
//! a `Geom_Plane` face carrying the same wires and the same tolerance, and the
//! substitution is remembered in `myFacePlanfaceMap` (and in `myFaces` when the
//! face was already registered there).
//!
//! Port guards (no OCCT branch changed):
//!
//! * `BRep_Tool::Surface(aF)` answers a null handle for a face without a
//!   surface; `GeomLib_IsPlanarSurface(null, ...)` would fault in OCCT. The
//!   port answers `aIsNonPlanarFound = true` for that face.
//! * `aPlanarityChecker.Plan()` raises `StdFail_NotDone` when `IsPlanar()` is
//!   false; the port only calls it inside the `IsPlanar()` branch, so the
//!   `Err` arm is unreachable (it is still handled to keep the `Result`).
//! * `aBAS.D1(...)` on a null surface would fault; the port treats a missing
//!   `D1` as non planar.
//! * `TopExp::Vertices(E, V1, V2, true)` in `RemoveSeamAndDegeneratedEdges`
//!   (`:5607`, `:5616`) is the oriented pair, [`BRepOffsetTool::edge_vertices`].
//! * a vertex-less edge (a degenerate without vertices) makes OCCT compare two
//!   null `TopoDS_Vertex` handles, which `IsSame` reports as equal; the port's
//!   [`same_vertex`] reproduces that on `Option<Vertex>`.

use std::sync::Arc;

use occt_core::precision::Precision;
use occt_geom::geom_lib_is_planar_surface::IsPlanarSurface;
use occt_geom::GeomPlane;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::brep_offset::BRepOffsetTool;
use crate::builder::TopoBuilder;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{children_keep_orientation, explorer_shapes, is_same, shape_map_key};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::BRepOffsetMakeOffset;

/// `aBAS.GetType() == GeomAbs_Plane` (`cxx:5457`).
fn is_plane_face(a_face: &Face) -> bool {
    SurfaceAdaptor::of_face(a_face).get_type() == GeomAbsSurfaceType::Plane
}

/// The midpoint of one parameter range, choosing `+-1` off an infinite end
/// exactly as `cxx:5474-5507` does.
fn mid_parameter(a1: f64, a2: f64) -> f64 {
    let is_inf1 = Precision::is_infinite(a1);
    let is_inf2 = Precision::is_infinite(a2);
    if !is_inf1 && !is_inf2 {
        (a1 + a2) / 2.0
    } else if is_inf1 && !is_inf2 {
        a2 - 1.0
    } else if !is_inf1 && is_inf2 {
        a1 + 1.0
    } else {
        0.0
    }
}

/// `IsPlanar()` (`cxx:5447-5552`).
pub fn is_planar(this: &mut BRepOffsetMakeOffset) -> bool {
    let mut a_is_non_planar_found = false;
    let a_builder = TopoBuilder::new();

    let Some(shape) = this.my_shape.clone() else {
        return true;
    };

    for a_f in explorer_shapes(&shape, ShapeType::Face) {
        if is_plane_face(&Face(a_f.clone())) {
            continue;
        }

        if !this.my_is_linearization_allowed {
            a_is_non_planar_found = true;
            continue;
        }

        // define the tolerance
        let a_tol_for_face = BRepTool::face_tolerance(&Face(a_f.clone()));

        let Some(a_surf) = BRepTool::face_surface(&Face(a_f.clone())) else {
            a_is_non_planar_found = true;
            continue;
        };
        let a_planarity_checker = IsPlanarSurface::new(&a_surf, Precision::CONFUSION);
        if !a_planarity_checker.is_planar() {
            a_is_non_planar_found = true;
            continue;
        }
        let Ok(mut a_pln) = a_planarity_checker.plan() else {
            a_is_non_planar_found = true;
            continue;
        };

        let (u1, u2) = a_surf.u_range();
        let (v1, v2) = a_surf.v_range();
        let um = mid_parameter(u1, u2);
        let vm = mid_parameter(v1, v2);

        let a_bas = SurfaceAdaptor::of_face(&Face(a_f.clone()));
        let Some((_, a_d1, a_d2)) = a_bas.d1(um, vm) else {
            a_is_non_planar_found = true;
            continue;
        };
        let a_norm = a_d1.crossed(&a_d2);
        let mut a_pln_norm = a_pln.position().direction();
        let pn = a_pln_norm.xyz();
        if a_norm.x() * pn.x() + a_norm.y() * pn.y() + a_norm.z() * pn.z() < 0.0 {
            a_pln_norm.reverse();
            // `aPln.SetAxis(gp_Ax1(aPln.Position().Location(), aPlnNorm))`
            // (`cxx:5518-5519`): `gp_Ax3::SetAxis` keeps the location and lets
            // `SetDirection` recompute X/Y from the old X direction, which for
            // the direct placement produced by `GeomLib_IsPlanarSurface` keeps
            // the frame direct - the same recomputation `GpAx3::set_direction`
            // performs here.
            let mut a_pos = a_pln.position();
            a_pos.set_direction(a_pln_norm);
            a_pln.set_position(&a_pos);
        }

        let a_plane = GeomPlane::new(a_pln);
        let mut a_planar_face =
            a_builder.make_face_with_tolerance(Arc::new(a_plane), &[], a_tol_for_face);

        let mut a_face_forward = a_f.clone();
        a_face_forward.set_orientation(Orientation::Forward);
        for a_wire in children_keep_orientation(&a_face_forward) {
            a_builder.add(&mut a_planar_face.0, &a_wire);
        }
        remove_seam_and_degenerated_edges(&a_planar_face, &Face(a_face_forward.clone()));

        this.my_face_planface_map
            .insert(shape_map_key(&a_f), a_planar_face.0.clone());
        if this.my_faces.contains(&a_f) {
            this.my_faces.remove_key(&a_f);
            this.my_faces.add(&a_planar_face.0);
        }
    }

    !a_is_non_planar_found
}

/// `TopoDS_Vertex::IsSame` on the nullable pairs `TopExp::Vertices` can take.
fn same_vertex(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => is_same(&x.0, &y.0),
        (None, None) => true,
        _ => false,
    }
}

/// `RemoveSeamAndDegeneratedEdges(theFace, theOldFace)` (`cxx:5556-5633`):
/// once a face became planar its seams and degenerated edges are useless, and
/// the wires are rebuilt from the remaining edges.
///
/// `theFace` is a freshly built face whose TShape is still unique to this call,
/// so mutating it through a local handle (as OCCT does at `:5594-5630`) is what
/// the port does too.
fn remove_seam_and_degenerated_edges(the_face: &Face, the_old_face: &Face) {
    let a_builder = TopoBuilder::new();

    let mut a_face = the_face.clone();
    a_face.0.set_orientation(Orientation::Forward);

    let mut a_is_deg_or_seam_found = false;
    let mut a_e_seq: Vec<TopoShape> = Vec::new();
    for an_edge in explorer_shapes(&a_face.0, ShapeType::Edge) {
        let an_edge = Edge(an_edge);
        if BRepTool::is_degenerated(&an_edge) || is_really_closed(&an_edge, the_old_face) {
            a_is_deg_or_seam_found = true;
        } else {
            a_e_seq.push(an_edge.0);
        }
    }

    if !a_is_deg_or_seam_found {
        return;
    }

    // Reconstruct wires
    let a_w_list = children_keep_orientation(&a_face.0);

    a_builder.make_free(&a_face.0);
    for a_wire in &a_w_list {
        a_builder.remove(&mut a_face.0, a_wire);
    }

    while !a_e_seq.is_empty() {
        let mut a_new_wire_edges: Vec<Edge> = Vec::new();
        let mut a_cur_edge = Edge(a_e_seq[0].clone());
        a_new_wire_edges.push(a_cur_edge.clone());
        a_e_seq.remove(0);
        let (a_first_vertex, mut a_cur_vertex) = BRepOffsetTool::edge_vertices(&a_cur_edge);
        while !same_vertex(&a_first_vertex, &a_cur_vertex) {
            let mut found: Option<(usize, Edge, Option<Vertex>)> = None;
            for (ind, candidate) in a_e_seq.iter().enumerate() {
                let candidate = Edge(candidate.clone());
                let (a_v1, a_v2) = BRepOffsetTool::edge_vertices(&candidate);
                if same_vertex(&a_v1, &a_cur_vertex) {
                    found = Some((ind, candidate, a_v2));
                    break;
                }
            }
            // error occurred: wire is not closed
            let Some((ind, edge, a_v2)) = found else {
                break;
            };

            a_cur_edge = edge;
            a_new_wire_edges.push(a_cur_edge.clone());
            a_e_seq.remove(ind);
            a_cur_vertex = a_v2;
        }

        let a_new_wire = a_builder.make_wire(&a_new_wire_edges);
        a_builder.add_wire(&mut a_face, &a_new_wire);
    }
}
