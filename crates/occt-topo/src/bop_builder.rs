//! Exact boolean operations for B-Rep solids.
//!
//! [`boolean`] is `BOPAlgo_BOP` via [`crate::bop_builder2`]. Empty-face
//! operands go to the same engine, whose `CheckData` skips them with
//! `BOPAlgo_AlertEmptyShape` (`BOPAlgo_BOP.cxx:162-167`, `:203-209`) — the
//! previous voxel/mesh fallback has no OCCT counterpart (audit A5/T-41).
//! Disjoint bounding boxes short-circuit. There is no planar 2-D arrangement on
//! this path (`boolean_planar_legacy` remains only as a leftover module).

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::{GeomPlane, Surface};

use crate::abs::ShapeType;
use crate::brep_extrema::is_inside;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::inttools::{edge_edge_intersections, edge_face_intersections};
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::shell_check::{shell_invariants, shell_is_closed};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};

pub use crate::bop_builder_core::{BoolOp, BooleanResult};
pub use crate::bop_builder_planar::boolean_planar_legacy;
pub(crate) use crate::bop_builder_core::{
    disjoint_result, empty_result, single_shape_result, validate,
};
pub(crate) use crate::bop_builder_planar::*;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Exact boolean on two solids (`BRepAlgoAPI_*` / `BOPAlgo_BOP`).
///
/// Green path: [`crate::bop_builder2::builder_bop_with_fuzzy`]. A BOPAlgo
/// error is returned as-is; the planar 2-D arrangement is not a fallback.
pub fn boolean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let fa = faces_of(a);
    let fb = faces_of(b);

    if fa.is_empty() || fb.is_empty() {
        // OCCT has no mesh/voxel boolean: an empty argument is skipped with
        // `BOPAlgo_AlertEmptyShape` and the operation proceeds
        // (`BOPAlgo_BOP.cxx:162-167`, `:203-209`), which is what the engine does.
        let shape = crate::bop_builder2::builder_bop_with_fuzzy(
            std::slice::from_ref(a),
            std::slice::from_ref(b),
            crate::bop_builder_dispatch::to_bool_op2(op),
            tol,
        )?;
        return Ok(single_shape_result(&shape));
    }

    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }

    let mut r = boolean_via_bopalgo(a, b, op, tol)?;
    validate(a, b, op, &mut r);
    Ok(r)
}

/// `BOPAlgo_BOP` wrapper: Fuse/Cut/Common through [`crate::bop_builder2`].
fn boolean_via_bopalgo(
    a: &TopoShape,
    b: &TopoShape,
    op: BoolOp,
    tol: f64,
) -> Result<BooleanResult, String> {
    let op2 = match op {
        BoolOp::Fuse => crate::bop_builder2::BoolOp2::Fuse,
        BoolOp::Cut => crate::bop_builder2::BoolOp2::Cut,
        BoolOp::Common => crate::bop_builder2::BoolOp2::Common,
    };
    let shape = crate::bop_builder2::builder_bop_with_fuzzy(
        &[a.clone()],
        &[b.clone()],
        op2,
        tol,
    )?;
    let shape = unwrap_bopalgo_shape(shape);
    let mut r = single_shape_result(&shape);
    if let Some(sh) = r.shells.first() {
        if !shell_is_closed(sh) {
            r.warnings.push("result shell is not closed".into());
        }
    }
    Ok(r)
}

/// A BOPAlgo result is a compound of argument images. Unwrap a single solid.
fn unwrap_bopalgo_shape(shape: TopoShape) -> TopoShape {
    if !shape.is_compound() {
        return shape;
    }
    let solids = shapes_of(&shape, ShapeType::Solid);
    if solids.len() == 1 {
        solids[0].clone()
    } else {
        shape
    }
}

pub use crate::bop_builder_dispatch::*;
pub use crate::bop_builder_repair::*;
pub use crate::bop_builder_report::*;
pub use crate::bop_builder_heal::*;
pub use crate::bop_builder_splitapi::*;

#[cfg(test)]
#[path = "bop_builder_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "bop_builder_tests_api.rs"]
mod tests_api;

#[cfg(test)]
#[path = "bop_builder_tests_heal.rs"]
mod tests_heal;
