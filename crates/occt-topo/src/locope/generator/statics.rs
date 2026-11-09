//! File-local statics of `LocOpe_Generator.cxx`: the three `ToFuse`
//! overloads (`.cxx:1232-1378`) and `NewParameter` (`.cxx:1382-1447`).
//!
//! `ToFuse` answers whether two pieces of geometry are the same *analytic*
//! object - a plane coplanar with another plane, a line coaxial with another
//! line - which is what tells the generator that two faces (or two edges) came
//! from the same base entity and can be merged. Only planes and lines ever
//! answer `true` (`.cxx:1271-1287`, `.cxx:1340-1350`); every other pair of
//! dynamic types stops at the type comparison.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::precision::{epsilon, ANGULAR, CONFUSION};
use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::locope::shape_maps::ShapeSet;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::explorer_shapes;
use crate::abs::ShapeType;

/// `BRep_Tool::Curve(E, L, First, Last)` with the location applied
/// (`BRep_Tool.cxx`), then the `Geom_TrimmedCurve` peel of the call sites
/// (`.cxx:1299-1333`).
fn edge_curve_stored(e: &TopoShape) -> Option<Arc<dyn Curve>> {
    let edge = Edge::wrap(e.clone())?;
    let c = BRepTool::edge_curve_world(&edge)?;
    Some(if c.is_geom_trimmed() { c.trim_basis().unwrap_or(c) } else { c })
}

/// `static bool ToFuse(const TopoDS_Face&, const TopoDS_Face&)`
/// (`.cxx:1232-1287`).
pub(super) fn to_fuse_faces(f1: &TopoShape, f2: &TopoShape) -> bool {
    if f1.is_null() || f2.is_null() {
        return false;
    }
    let (Some(f1w), Some(f2w)) = (Face::wrap(f1.clone()), Face::wrap(f2.clone())) else {
        return false;
    };
    let (Some(s1), Some(s2)) = (BRepTool::face_surface(&f1w), BRepTool::face_surface(&f2w)) else {
        return false;
    };
    // `cxx:1254-1263`: a rectangular trim is replaced by its basis surface.
    let s1 = match s1.rectangular_trimmed_basis() {
        Some(b) => b,
        None => s1,
    };
    let s2 = match s2.rectangular_trimmed_basis() {
        Some(b) => b,
        None => s2,
    };
    // `cxx:1266-1268`: unequal dynamic types never fuse. Below only
    // `Geom_Plane` answers `true` (`cxx:1271-1287`), so a non-plane pair stops
    // here and a plane paired with anything else stops at the `gp_pln` miss.
    let (Some(mut pl1), Some(mut pl2)) = (s1.gp_pln(), s2.gp_pln()) else {
        return false;
    };
    // `cxx:1278-1279`: `pl1.Transform(loc1.Transformation())`, where `loc1` is
    // the face location `BRep_Tool::Surface(F1, loc1)` handed back.
    if !f1.location().is_identity() {
        pl1.transform(&f1.location().transformation());
    }
    if !f2.location().is_identity() {
        pl2.transform(&f2.location().transformation());
    }
    pl1.position()
        .is_coplanar(&pl2.position(), CONFUSION, ANGULAR)
}

/// `static bool ToFuse(const TopoDS_Edge&, const TopoDS_Edge&)`
/// (`.cxx:1292-1352`).
pub(super) fn to_fuse_edges(e1: &TopoShape, e2: &TopoShape) -> bool {
    if e1.is_null() || e2.is_null() {
        return false;
    }
    let (Some(c1), Some(c2)) = (edge_curve_stored(e1), edge_curve_stored(e2)) else {
        return false;
    };
    // `cxx:1335-1338`: unequal dynamic types never fuse. Only `Geom_Line`
    // answers `true` (`cxx:1340-1350`), so a miss below covers the rest.
    let (Some(li1), Some(li2)) = (c1.gp_line(), c2.gp_line()) else {
        return false;
    };
    li1.position().is_coaxial_tol(li2.position(), ANGULAR, CONFUSION)
}

/// `static bool ToFuse(const TopoDS_Edge&, const TopoDS_Face&,
/// const TopoDS_Vertex&, const NCollection_Map<...>&)` (`.cxx:1357-1378`).
///
/// True when every *other* edge of `f` that shares `v` is already scheduled for
/// removal, i.e. nothing but `e` still hangs on `v`.
pub(super) fn to_fuse_edge_face_vertex(
    e: &TopoShape,
    f: &TopoShape,
    v: &TopoShape,
    to_remove: &ShapeSet,
) -> bool {
    for edge in explorer_shapes(f, ShapeType::Edge) {
        if crate::topo_tools::is_same(&edge, e) {
            continue;
        }
        let (vf, vl) = super::edge_vertex_shapes(&edge);
        let touches = |vtx: &Option<TopoShape>| vtx.as_ref().is_some_and(|x| crate::topo_tools::is_same(x, v));
        if (touches(&vf) || touches(&vl)) && !to_remove.contains(&edge) {
            return false;
        }
    }
    true
}

/// `static double NewParameter(const TopoDS_Edge&, const TopoDS_Vertex&,
/// const TopoDS_Edge&, const TopoDS_Vertex&)` (`.cxx:1382-1447`).
///
/// The parameter of `new_vtx` on `edg`, re-expressed on the period of a circle
/// so that it lands next to the parameter the matching vertex already carries
/// on `new_edg` (`.cxx:1416-1442`).
pub(super) fn new_parameter(
    edg: &TopoShape,
    vtx: &TopoShape,
    new_edg: &TopoShape,
    new_vtx: &TopoShape,
) -> f64 {
    let Some(new_v) = Vertex::wrap(new_vtx.clone()) else {
        return 0.0;
    };
    let p = BRepTool::vertex_point(&new_v);
    let Some(c) = edge_curve_stored(edg) else {
        return 0.0;
    };
    if let Some(lin) = c.gp_line() {
        // `cxx:1407-1410`: `ElCLib::Parameter(Line, P)`.
        return clib::line_parameter(&lin, &p);
    }
    let Some(circ) = c.gp_circ() else {
        // `cxx:1446`: no other curve type is handled.
        return 0.0;
    };
    // `cxx:1411-1414`.
    let mut prm = clib::circle_parameter(&circ.position(), &p);
    // `cxx:1416-1418`: `Vtx` comes from an exploration of `Edg` oriented
    // FORWARD, so its own orientation is the reverse of the edge's view.
    let orient = vtx.orientation().reversed();
    if matches!(orient, crate::abs::Orientation::Forward | crate::abs::Orientation::Reversed) {
        // `cxx:1419-1430`.
        let forward = new_edg.oriented(crate::abs::Orientation::Forward);
        let mut found = None;
        for v in explorer_shapes(&forward, ShapeType::Vertex) {
            if v.orientation() == orient {
                found = Some(v);
                break;
            }
        }
        if let Some(found) = found {
            let prmmax = BRepTool::parameter_of_vertex(new_edg, &found);
            if (prmmax - prm).abs() <= epsilon(2.0 * std::f64::consts::PI) {
                if orient == crate::abs::Orientation::Reversed {
                    prm -= 2.0 * std::f64::consts::PI;
                } else {
                    prm += 2.0 * std::f64::consts::PI;
                }
            }
        }
    }
    prm
}
