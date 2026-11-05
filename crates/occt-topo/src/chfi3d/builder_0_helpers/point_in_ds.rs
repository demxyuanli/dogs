//! The `ChFiDS_CommonPoint` / `SurfData` bookkeeping statics of
//! `ChFi3d_Builder_0.cxx`: point tolerance, point index and the section index
//! of a vertex along a stripe.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx` (OCCT
//! 8.0.0):
//!
//! - `ChFi3d_SetPointTolerance(DStr, box, IP)` (`cxx:319-331`),
//! - `ChFi3d_IndexPointInDS(P1, DStr)` (`cxx:2325-2338`),
//! - `ChFi3d_IndexOfSurfData(V1, CD, sens)` (`cxx:3446-3491`),
//! - `ChFi3d_EdgeFromV1(V1, CD, sens)` (`cxx:3495-3536`).
//!
//! Translation notes:
//! - `BRep_Builder B; B.UpdateVertex(P1.Vertex(), P1.Point(), P1.Tolerance())`
//!   (`cxx:2332-2334`) is the port's `Vertex::set_point` plus the
//!   max-tolerance rule of `BRep_TVertex::UpdateTolerance`, the same pair
//!   `crate::shhealing::wire_self_inter` uses for that call,
//! - `DStr.ChangePoint(IP)` answers `myEmptyPoint` for an unbound index, so a
//!   missing point is a no-op rather than an error,
//! - `CD->SetOfSurfData().IsNull()` (`cxx:3477`) is the port's empty sequence,
//! - `Standard_ConstructionError` (`cxx:3488`) is the `Err`.

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use crate::brep_tool::BRepTool;
use crate::chfi_ds::{CommonPoint, StripeHandle};
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::is_same;
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};
use crate::topopebrep::ds::{DataStructure, Point};

/// `ChFi3d_SetPointTolerance(DStr, box, IP)` (`cxx:319-331`).
///
/// The tolerance of the DS point `IP` becomes the diagonal of `box` times
/// `1.5`; a void box has no `Get` and is a no-op (OCCT reads uninitialized
/// coordinates there).
pub fn ch_fi3d_set_point_tolerance(ds: &mut DataStructure, bx: &BndBox, ip: i32) {
    let Some((a, b, c, d, e, f)) = bx.get() else {
        // `box.Get(...)` on a void box: OCCT would read the unset corner
        // values; there is no diagonal to compute.
        return;
    };
    // `d -= a; e -= b; f -= c; d *= d; e *= e; f *= f;
    // vtol = sqrt(d + e + f) * 1.5;` (`cxx:322-329`).
    let dx = d - a;
    let dy = e - b;
    let dz = f - c;
    let vtol = (dx * dx + dy * dy + dz * dz).sqrt() * 1.5;
    // `DStr.ChangePoint(IP).Tolerance(vtol);` (`cxx:330`).
    if let Some(point) = ds.change_point(ip) {
        point.set_tolerance(vtol);
    }
}

/// `ChFi3d_IndexPointInDS(P1, DStr)` (`cxx:2325-2338`).
///
/// A common point on a vertex updates that vertex (point and tolerance) and
/// answers the DS index of the shape; otherwise a new DS point is added.
pub fn ch_fi3d_index_point_in_ds(p1: &CommonPoint, ds: &mut DataStructure) -> i32 {
    // `if (P1.IsVertex()) { BRep_Builder B; B.UpdateVertex(P1.Vertex(),
    // P1.Point(), P1.Tolerance()); return DStr.AddShape(P1.Vertex()); }`
    // (`cxx:2327-2335`).
    if let Some(vertex) = p1.vertex() {
        update_vertex_point(&Vertex(vertex.clone()), *p1.point(), p1.tolerance());
        return ds.add_shape(vertex);
    }
    // `return DStr.AddPoint(TopOpeBRepDS_Point(P1.Point(), P1.Tolerance()));`
    // (`cxx:2337`).
    ds.add_point(&Point::with_tolerance(*p1.point(), p1.tolerance()))
}

/// `BRep_Builder::UpdateVertex(V, P, Tol)` (`BRep_Builder.cxx:1203-1212`):
/// `BRep_TVertex::Pnt` then `UpdateTolerance`, which only raises. Locations are
/// identity throughout this port, so the `V.Location().Inverted()` transform
/// of the C++ call is a no-op.
fn update_vertex_point(v: &Vertex, p: GpPnt, tol: f64) {
    v.set_point(p);
    v.set_tolerance(BRepTool::vertex_tolerance(v).max(tol));
}

/// `ChFi3d_IndexOfSurfData(V1, CD, sens)` (`cxx:3446-3491`).
///
/// Answers the 1-based index of the section of `stripe` that starts (or ends,
/// `sens = -1`) on the vertex `v1`.
pub fn ch_fi3d_index_of_surf_data(
    v1: &TopoShape,
    stripe: &StripeHandle,
    sens: &mut i32,
) -> Result<i32, String> {
    let (spine, nb_sections) = {
        let guard = stripe
            .read()
            .map_err(|_| "ChFi3d_IndexOfSurfData: poisoned stripe".to_string())?;
        let spine = guard
            .spine()
            .cloned()
            .ok_or_else(|| "ChFi3d_IndexOfSurfData: the stripe has no spine".to_string())?;
        (spine, guard.set_of_surf_data().len() as i32)
    };

    // `sens = 1;` (`cxx:3450`).
    *sens = 1;

    // `const TopoDS_Edge& E = spine->Edges(1); if (E.Orientation() ==
    // TopABS_REVERSED) Vref = TopExp::LastVertex(E); else
    // Vref = TopExp::FirstVertex(E);` (`cxx:3452-3460`).
    let e = {
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_IndexOfSurfData: poisoned spine".to_string())?;
        guard
            .edges(1)
            .cloned()
            .ok_or_else(|| "ChFi3d_IndexOfSurfData: the spine has no first edge".to_string())?
    };

    // `if (Vref.IsSame(V1)) Index = 1;` (`cxx:3461-3464`).
    if let Some(vref) = first_vertex(&e) {
        if is_same(&vref, v1) {
            return Ok(1);
        }
    }

    // `const TopoDS_Edge& E1 = spine->Edges(spine->NbEdges()); if
    // (E1.Orientation() == TopAbs_REVERSED) Vref = TopExp::FirstVertex(E1);
    // else Vref = TopExp::LastVertex(E1); sens = -1;` (`cxx:3465-3476`).
    let e1 = {
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_IndexOfSurfData: poisoned spine".to_string())?;
        guard
            .edges(guard.nb_edges())
            .cloned()
            .ok_or_else(|| "ChFi3d_IndexOfSurfData: the spine has no last edge".to_string())?
    };
    *sens = -1;

    // `if (CD->SetOfSurfData().IsNull()) return 0;` (`cxx:3477-3480`).
    if nb_sections == 0 {
        return Ok(0);
    }
    // `else if (Vref.IsSame(V1)) Index = CD->SetOfSurfData()->Length();`
    // (`cxx:3481-3484`).
    if let Some(vref) = last_vertex(&e1) {
        if is_same(&vref, v1) {
            return Ok(nb_sections);
        }
    }
    // `else throw Standard_ConstructionError("ChFi3d_IndexOfSurfData() -
    // wrong construction parameters");` (`cxx:3485-3489`).
    Err("ChFi3d_IndexOfSurfData() - wrong construction parameters".to_string())
}

/// `ChFi3d_EdgeFromV1(V1, CD, sens)` (`cxx:3495-3536`): the first (or last,
/// `sens = -1`) edge of the stripe's spine whose stored vertex is `v1`.
///
/// `TopExp::FirstVertex` / `LastVertex` already fold in the edge orientation
/// (`tools2d::edge_vertices_cumulative`), so the two `Orientation() ==
/// TopAbs_REVERSED` swaps of the OCCT body (`cxx:3503-3509`, `:3519-3525`) are
/// the plain calls here. `Standard_ConstructionError` (`cxx:3534`) is the
/// `Err`, and `sens` answers `1` on the first edge, `-1` on the last one
/// (`cxx:3500`, `:3528`) — it is assigned before the check, as in OCCT.
pub fn ch_fi3d_edge_from_v1(
    v1: &TopoShape,
    stripe: &StripeHandle,
    sens: &mut i32,
) -> Result<TopoShape, String> {
    let spine = {
        let guard = stripe
            .read()
            .map_err(|_| "ChFi3d_EdgeFromV1: poisoned stripe".to_string())?;
        guard
            .spine()
            .cloned()
            .ok_or_else(|| "ChFi3d_EdgeFromV1: the stripe has no spine".to_string())?
    };
    // `sens = 1;` (`cxx:3500`).
    *sens = 1;
    // `const TopoDS_Edge& E = spine->Edges(1);` (`cxx:3502`).
    let e = {
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_EdgeFromV1: poisoned spine".to_string())?;
        guard
            .edges(1)
            .cloned()
            .ok_or_else(|| "ChFi3d_EdgeFromV1: the spine has no first edge".to_string())?
    };
    // `if (Vref.IsSame(V1)) return E;` (`cxx:3510-3513`).
    if let Some(vref) = first_vertex(&e) {
        if is_same(&vref, v1) {
            return Ok(e);
        }
    }
    // `const TopoDS_Edge& E1 = spine->Edges(spine->NbEdges()); sens = -1;`
    // (`cxx:3515-3528`).
    let e1 = {
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_EdgeFromV1: poisoned spine".to_string())?;
        guard
            .edges(guard.nb_edges())
            .cloned()
            .ok_or_else(|| "ChFi3d_EdgeFromV1: the spine has no last edge".to_string())?
    };
    *sens = -1;
    // `if (Vref.IsSame(V1)) return E1;` (`cxx:3529-3532`).
    if let Some(vref) = last_vertex(&e1) {
        if is_same(&vref, v1) {
            return Ok(e1);
        }
    }
    // `else throw Standard_ConstructionError("ChFi3d_IndexOfSurfData() -
    // wrong construction parameters");` (`cxx:3533-3535`), the same message
    // OCCT reuses here.
    Err("ChFi3d_IndexOfSurfData() - wrong construction parameters".to_string())
}
