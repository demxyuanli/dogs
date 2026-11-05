//! `FindInternals` (`BRepFill_AdvancedEvolved.cxx:1630-1655`),
//! `RemoveInternalWires` (`cxx:1657-1679`), `ProcessVertex`
//! (`cxx:1682-1827`) and `ReduceVertexTolerance` (`cxx:1831-1850`).
//!
//! `ReduceVertexTolerance` walks every vertex of a face, recomputes the largest
//! distance between the vertex point and the points the incident edges and
//! pcurves answer at that vertex, and — when that distance is not larger than
//! the stored vertex tolerance — shrinks the tolerance to the recomputed value
//! plus a `0.01%` margin.
//!
//! Translation notes:
//! - `aLE` / `aLF` come from `TopExp::MapShapesAndUniqueAncestors`
//!   ([`map_shapes_and_unique_ancestors`]); `ProcessVertex` is independent per
//!   vertex, so the map *iteration order* (the reference's
//!   `IndexedDataMap`) is not observable.
//! - `TE->Curves()` splits into the `BRep_Curve3D` arm
//!   ([`BRepTool::edge_curve_world`]) and the `IsCurveOnSurface` arm
//!   ([`GeometryRegistry::edge_pcurve_reps`]). OCCT compares the points in the
//!   vertex's own location frame (`L = (Eloc * loc).Predivided(aV.Location())`,
//!   `cxx:1734-1736`); the same distance is obtained in world coordinates
//!   because both sides are moved by one rigid transformation, and the port has
//!   no per-representation location (they are folded into
//!   `edge_curve_world` / `face_surface_world`).
//! - `BRep_Tool::Parameter(aVx, aE, aS, L)` (`cxx:1766`) has no exact port: the
//!   port does not store `BRep_PointOnCurveOnSurface` representations. The
//!   closest equivalent is [`BRepTool::parameter_of_vertex_on_face`], which
//!   answers the pcurve parameter for the matching face and falls back to
//!   `Parameter(V, E)` exactly like the reference does when no point
//!   representation matches the surface.
//! - `BRep_Builder::Remove(aF, W)` and `Face::Free` are
//!   [`TopoBuilder::remove`] / [`TopoShape::set_free`].

use std::sync::Arc;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    explorer_shapes, is_same, map_shapes_and_unique_ancestors, shape_map_key,
};

use super::AdvancedEvolved;

/// `FindInternals(theS, theLInt)` (`cxx:1630-1655`): the INTERNAL children of
/// `the_s`, or the children having an INTERNAL child.
fn find_internals(the_s: &TopoShape, the_l_int: &mut Vec<TopoShape>) {
    for a_ss in ShapeIterator::of_shape(the_s) {
        if a_ss.orientation() == Orientation::Internal {
            the_l_int.push(a_ss);
        } else {
            for it_ss in ShapeIterator::of_shape(&a_ss) {
                if it_ss.orientation() == Orientation::Internal {
                    the_l_int.push(a_ss.clone());
                    break;
                }
            }
        }
    }
}

/// `RemoveInternalWires(theShape)` (`cxx:1657-1679`): removes the internal
/// wires of every face of `the_shape` in place.
fn remove_internal_wires(the_shape: &TopoShape) {
    let builder = TopoBuilder::new();
    for face_shape in explorer_shapes(the_shape, ShapeType::Face) {
        let mut a_f = face_shape;
        let mut a_lw_to_remove: Vec<TopoShape> = Vec::new();
        find_internals(&a_f, &mut a_lw_to_remove);
        if !a_lw_to_remove.is_empty() {
            a_f.set_free(true);
            for w in &a_lw_to_remove {
                builder.remove(&mut a_f, w);
            }
            a_f.set_free(false);
        }
    }
}

/// `ProcessVertex(aV, aLE, aLF)` (`cxx:1682-1827`).
fn process_vertex(a_v: &TopoShape, a_le: &[TopoShape], a_lf: &[TopoShape]) {
    let mut a_tol_max2 = -1.0e6;

    // `cxx:1699-1701`: `aPV3D = TV->Pnt()`, `aTol = BRep_Tool::Tolerance(aV)`.
    let a_pv3d = BRepTool::vertex_point_world(&Vertex(a_v.clone()));
    let a_tol = BRepTool::vertex_tolerance(&Vertex(a_v.clone()));

    // `cxx:1703-1783`.
    for e_shape in a_le {
        let a_e = Edge(e_shape.clone());

        // `aVExp.Init(aE, TopAbs_VERTEX)`.
        for a_vx in explorer_shapes(&a_e.0, ShapeType::Vertex) {
            if !is_same(&a_vx, a_v) {
                continue;
            }

            // `cxx:1720-1725`.
            let an_or_v = a_vx.orientation();
            if an_or_v != Orientation::Forward && an_or_v != Orientation::Reversed {
                continue;
            }

            // `cxx:1736-1756`: the `BRep_Curve3D` representation.
            if let Some(a_c3d) = BRepTool::edge_curve_world(&a_e) {
                let a_param = BRepTool::parameter_of_vertex(&a_e.0, &a_vx);
                let a_pc3d = a_c3d.d0(a_param);
                let a_d2 = a_pv3d.square_distance(&a_pc3d);
                if a_d2 > a_tol_max2 {
                    a_tol_max2 = a_d2;
                }
            }

            // `cxx:1758-1779`: the `BRep_CurveOnSurface` representations.
            for (a_s, a_c2d, _first, _last) in
                GeometryRegistry::global().edge_pcurve_reps(&a_e.0)
            {
                // `cxx:1766`: `BRep_Tool::Parameter(aVx, aE, aS, L)`.
                let a_face = a_lf.iter().find(|f| {
                    BRepTool::face_surface(&Face((*f).clone()))
                        .map(|s| Arc::ptr_eq(&s, &a_s))
                        .unwrap_or(false)
                });
                let a_param = match &a_face {
                    Some(f) => {
                        BRepTool::parameter_of_vertex_on_face(&a_vx, &a_e, &Face((*f).clone()))
                    }
                    None => BRepTool::parameter_of_vertex(&a_e.0, &a_vx),
                };

                // `cxx:1768-1774`: pcurve point, surface point, distance.
                let a_pc2d = a_c2d.d0(a_param);
                let a_s_world = match &a_face {
                    Some(f) => BRepTool::face_surface_world(&Face((*f).clone())),
                    None => Some(a_s.clone()),
                };
                let Some(a_s_world) = a_s_world else {
                    continue;
                };
                let a_pc3d = a_s_world.d0(a_pc2d.x(), a_pc2d.y());
                let a_d2 = a_pv3d.square_distance(&a_pc3d);
                if a_d2 > a_tol_max2 {
                    a_tol_max2 = a_d2;
                }
            }
        }
    }

    // `cxx:1786-1789`.
    if a_tol_max2 < 0.0 {
        return;
    }

    a_tol_max2 = a_tol_max2.sqrt();
    if a_tol_max2 > a_tol {
        return;
    }

    // `cxx:1795-1818`.
    for e in a_le {
        let a_tol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
        if a_tol_max2 < a_tol_e {
            a_tol_max2 = a_tol_e;
        }
    }
    for f in a_lf {
        let a_tol_e = BRepTool::face_tolerance(&Face(f.clone()));
        if a_tol_max2 < a_tol_e {
            a_tol_max2 = a_tol_e;
        }
    }

    // `cxx:1820-1823`.
    if a_tol_max2 > a_tol {
        return;
    }

    // `cxx:1825-1827`: "update tolerance with a small margin".
    GeometryRegistry::global()
        .set_vertex_tolerance(a_v, a_tol_max2 + a_tol_max2 * 0.0001);
}

impl AdvancedEvolved {
    /// `ReduceVertexTolerance(aS)` (`cxx:1831-1850`).
    pub(super) fn reduce_vertex_tolerance(&self, a_s: &TopoShape) {
        // `cxx:1835-1840`.
        let a_ve_map = map_shapes_and_unique_ancestors(
            a_s,
            ShapeType::Vertex,
            ShapeType::Edge,
            false,
        );
        let a_vf_map = map_shapes_and_unique_ancestors(
            a_s,
            ShapeType::Vertex,
            ShapeType::Face,
            false,
        );

        // `cxx:1841-1849`.
        for a_v in explorer_shapes(a_s, ShapeType::Vertex) {
            let key = shape_map_key(&a_v);
            let a_le = a_ve_map.get(&key).cloned().unwrap_or_default();
            let a_lf = a_vf_map.get(&key).cloned().unwrap_or_default();
            process_vertex(&a_v, &a_le, &a_lf);
        }
    }

    /// `RemoveInternalWires(myResult)` (`cxx:828`).
    pub(super) fn remove_internal_wires(&mut self) {
        if let Some(result) = &self.result {
            remove_internal_wires(result);
        }
    }
}
