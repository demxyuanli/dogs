//! `BRep_Builder::UpdateVertex` - the three overloads that write a vertex's
//! parameter into an edge or a face.
//!
//! OCCT has one overload per carrier (`BRep_Builder.cxx`):
//!
//! * `UpdateVertex(V, Par, E, Tol)` (`cxx:1220-1312`) - every
//!   `BRep_GCurve` of the edge gets `First(Par)` / `Last(Par)`, or, for an
//!   `INTERNAL` / `EXTERNAL` vertex, an `UpdatePoints` entry.
//! * `UpdateVertex(V, Par, E, S, L, Tol)` (`cxx:1314-1408`) - the same, but
//!   only the representation whose surface is `S`.
//! * `UpdateVertex(Ve, U, V, F, Tol)` (`cxx:1420-1441`) - a
//!   `BRep_PointOnSurface` entry on the vertex, the producer that
//!   `BRep_Tool::Parameters(V, F)` reads back.
//!
//! The port keeps one entry point per overload. They are the only places that
//! write [`crate::tshape::VertexGeomCore::points`], which means the
//! representation list stays in lockstep with the parameter writes that used to
//! live in unrelated modules.
//!
//! `TopLoc_Location` is unported, so the `L` of the OCCT signatures has no
//! counterpart: `IsCurveOnSurface(S, L)` reduces to the surface handle identity
//! and the representation entries carry no location of their own.

use std::sync::Arc;

use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{children_keep_orientation, is_same};

/// `BRep_Builder::UpdateVertex(V, Par, E, Tol)` (`BRep_Builder.cxx:1220-1312`).
///
/// `cxx:1241-1270` resolves the vertex's orientation inside the edge
/// (`TopoDS_Iterator itv(E.Oriented(TopAbs_FORWARD))`, degenerate-edge fallback
/// to `V.Orientation()`), then `cxx:1272-1304` walks the edge's
/// representations: `FORWARD` writes `First(Par)`, `REVERSED` writes
/// `Last(Par)`, and anything else calls `UpdatePoints` against the
/// representation's 3D curve (`cxx:1293`) or p-curve plus surface
/// (`cxx:1299`). The tail (`cxx:1307-1311`) marks the vertex modified for the
/// non-oriented arms, raises the vertex tolerance and marks the edge modified.
///
/// `Standard_DomainError("BRep_Builder::Infinite parameter")`
/// (`cxx:1225-1228`) is the `Err` of the port.
///
/// Deviations: the port writes the edge-keyed parameter cache
/// (`GeometryRegistry::set_vertex_parameter_on_edge`) and moves the vertex's 3D
/// point to `C(Par)` in the non-oriented arm. OCCT stores the same parameter in
/// the representation list and leaves `BRep_TVertex::myPoint` untouched; the
/// cache is what `BRep_Tool::Parameter(V, E)` reads first in this port, and the
/// point write is the pre-existing behaviour of both call sites.
pub fn update_vertex_on_edge(
    v: &TopoShape,
    par: f64,
    e: &TopoShape,
    tol: f64,
) -> Result<(), String> {
    // `cxx:1225-1228`.
    if occt_core::precision::Precision::is_infinite(par) {
        return Err("BRep_Builder::Infinite parameter".into());
    }
    let reg = GeometryRegistry::global();
    let ori = vertex_orientation_in_edge(v, e);
    let (first, last) = reg.edge_parameters(e);
    match ori {
        Orientation::Forward => reg.set_ranges_all(e, par, last),
        Orientation::Reversed => reg.set_ranges_all(e, first, par),
        _ => {
            // `cxx:1272-1304`: every representation of the edge.
            let core = v.tshape.write().expect("poisoned TShape lock");
            let mut core = core;
            let vertex_core = core.vertex_core_mut();
            if let Some(c) = reg.edge_curve(e) {
                vertex_core.update_point_on_curve(par, &c);
            }
            for (surface, pc, _, _) in reg.edge_pcurve_reps(e) {
                vertex_core.update_point_on_curve_on_surface(par, &pc, &surface);
            }
            drop(core);
            reg.set_vertex_parameter_on_edge(e, v, par);
            if let Some(c) = reg.edge_curve(e) {
                reg.set_vertex(v, crate::tgeometry::VertexGeom {
                    point: c.d0(par),
                    tolerance: reg.vertex_tolerance(v),
                });
            }
        }
    }
    // `cxx:1307-1309`: `TV->Modified(true)` only for the non-oriented arms.
    if !matches!(ori, Orientation::Forward | Orientation::Reversed) {
        v.tshape.write().expect("poisoned TShape lock").set_modified(true);
    }
    // `cxx:1310`: `TV->UpdateTolerance(Tol)` - raise only.
    reg.raise_vertex_tolerance(v, tol);
    // `cxx:1311`: `TE->Modified(true)`.
    e.tshape.write().expect("poisoned TShape lock").set_modified(true);
    Ok(())
}

/// `BRep_Builder::UpdateVertex(V, Par, E, S, L, Tol)`
/// (`BRep_Builder.cxx:1314-1408`).
///
/// The representation walk stops at the first `BRep_GCurve` whose
/// `IsCurveOnSurface(S, L)` holds (`cxx:1385-1404`): `FORWARD` / `REVERSED`
/// rewrite that one representation's `First` / `Last`, anything else appends a
/// `BRep_PointOnCurveOnSurface` through `UpdatePoints` (`cxx:1392-1396`). A walk
/// that reaches the end throws `Standard_DomainError("BRep_Builder:: no
/// pcurve")` (`cxx:1400-1403`), the `Err` of the port.
///
/// In the port the p-curve representations are keyed by the surface handle
/// (`GeometryRegistry::edge_pcurves`), so `IsCurveOnSurface(S, L)` is that
/// handle identity, and the range write goes through
/// `set_pcurve_range_on_surface`.
pub fn update_vertex_on_edge_on_surface(
    v: &TopoShape,
    par: f64,
    e: &TopoShape,
    surface: &Arc<dyn Surface>,
    tol: f64,
) -> Result<(), String> {
    // `cxx:1327-1330`.
    if occt_core::precision::Precision::is_infinite(par) {
        return Err("BRep_Builder::Infinite parameter".into());
    }
    let reg = GeometryRegistry::global();
    let ori = vertex_orientation_in_edge(v, e);
    let mut matched = false;
    for (rep_surface, pc, first, last) in reg.edge_pcurve_reps(e) {
        if !Arc::ptr_eq(&rep_surface, surface) {
            continue;
        }
        match ori {
            Orientation::Forward => reg.set_pcurve_range_on_surface(e, surface, par, last),
            Orientation::Reversed => reg.set_pcurve_range_on_surface(e, surface, first, par),
            _ => {
                let mut ts = v.tshape.write().expect("poisoned TShape lock");
                ts.vertex_core_mut().update_point_on_curve_on_surface(par, &pc, surface);
                drop(ts);
                v.tshape.write().expect("poisoned TShape lock").set_modified(true);
            }
        }
        matched = true;
        break;
    }
    // `cxx:1400-1403`.
    if !matched {
        return Err("BRep_Builder:: no pcurve".into());
    }
    // `cxx:1405-1407`.
    reg.raise_vertex_tolerance(v, tol);
    e.tshape.write().expect("poisoned TShape lock").set_modified(true);
    Ok(())
}

/// `BRep_Builder::UpdateVertex(Ve, U, V, F, Tol)` (`BRep_Builder.cxx:1420-1441`).
///
/// One `UpdatePoints(lpr, U, V, S, L)` call (`cxx:1437`), i.e. a
/// `BRep_PointOnSurface` entry carrying the two parameters, followed by
/// `TV->UpdateTolerance(Tol)` and `TV->Modified(true)` (`cxx:1439-1440`).
pub fn update_vertex_on_face(
    v: &TopoShape,
    u: f64,
    vv: f64,
    face: &TopoShape,
    tol: f64,
) -> Result<(), String> {
    let Some(surface) = BRepTool::face_surface(&crate::shape::Face(face.clone())) else {
        return Err("BRep_Builder::UpdateVertex: face has no surface".into());
    };
    {
        let mut ts = v.tshape.write().expect("poisoned TShape lock");
        ts.vertex_core_mut().update_point_on_surface(u, vv, &surface);
    }
    let reg = GeometryRegistry::global();
    reg.raise_vertex_tolerance(v, tol);
    v.tshape.write().expect("poisoned TShape lock").set_modified(true);
    Ok(())
}

/// `cxx:1241-1270`: the orientation the vertex carries inside the edge.
///
/// `TopoDS_Iterator itv(E.Oriented(TopAbs_FORWARD))` is the direct-children
/// walk, so the stored orientations are the ones compared; the degenerate-edge
/// fallback (`cxx:1251-1254`) only fires when the edge stores no vertex child.
fn vertex_orientation_in_edge(v: &TopoShape, e: &TopoShape) -> Orientation {
    let verts = children_keep_orientation(e);
    let mut ori = Orientation::Internal;
    let degenerate =
        e.shape_type() == ShapeType::Edge && BRepTool::is_degenerated(&crate::shape::Edge(e.clone()));
    if verts.is_empty() && degenerate {
        ori = v.orientation();
    }
    for vcur in &verts {
        if is_same(vcur, v) {
            ori = vcur.orientation();
            if ori == v.orientation() {
                break;
            }
        }
    }
    ori
}
