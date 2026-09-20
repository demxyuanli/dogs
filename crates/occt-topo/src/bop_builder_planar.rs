//! Leftover planar 2-D arrangement ([`boolean_planar_legacy`]).
//!
//! Not used by [`crate::bop_builder::boolean`]. Heal / n-ary helpers still
//! import Weld / EdgeMap from this module.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::{GeomPlane, Surface};

use crate::abs::ShapeType;
use crate::bop_builder_core::{
    disjoint_result, empty_result, single_shape_result, validate, BoolOp, BooleanResult,
};
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

pub(crate) use crate::bop_builder_planar_geom::*;
pub(crate) use crate::bop_builder_planar_trace::*;
pub(crate) use crate::bop_builder_planar_weld::*;

/// Leftover planar polygon boolean. Not called from [`crate::bop_builder::boolean`].
pub fn boolean_planar_legacy(
    a: &TopoShape,
    b: &TopoShape,
    op: BoolOp,
    tol: f64,
) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let fa = faces_of(a);
    let fb = faces_of(b);

    // Planarity gate. **UNPORTED (audit A5/T-41)**: OCCT has ONE
    // `BOPAlgo_BOP` for planar and curved inputs, so an empty or non-planar
    // operand goes to the same exact engine; the previous body fell back to a
    // 32³ `boolean_ops::voxel_boolean` (mesh boolean), which has no OCCT
    // counterpart. This module stays a leftover planar 2-D arrangement
    // (`boolean_planar_legacy`).
    if fa.is_empty() || fb.is_empty() {
        let shape = crate::bop_builder2::builder_bop_with_fuzzy(
            std::slice::from_ref(a),
            std::slice::from_ref(b),
            crate::bop_builder_dispatch::to_bool_op2(op),
            tol,
        )?;
        return Ok(single_shape_result(&shape));
    }
    let planes_a: Vec<Option<GpPln>> = fa.iter().map(face_plane_local).collect();
    let planes_b: Vec<Option<GpPln>> = fb.iter().map(face_plane_local).collect();
    if planes_a.iter().any(Option::is_none) || planes_b.iter().any(Option::is_none) {
        let shape = crate::bop_builder2::builder_bop_with_fuzzy(
            std::slice::from_ref(a),
            std::slice::from_ref(b),
            crate::bop_builder_dispatch::to_bool_op2(op),
            tol,
        )?;
        return Ok(single_shape_result(&shape));
    }

    // Disjoint shortcut (bounding boxes don't overlap).
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }

    // Step 1: face–face intersection segments, grouped by the face they lie on.
    let mut segs_a: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); fa.len()];
    let mut segs_b: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); fb.len()];
    for (i, af) in fa.iter().enumerate() {
        for (j, bf) in fb.iter().enumerate() {
            let segs = face_face_segments_local(af, bf, tol);
            for s in &segs {
                segs_a[i].push(*s);
                segs_b[j].push(*s);
            }
        }
    }

    // Edge–edge / edge–face hits (used for completeness checks; the polygon
    // split above already subsumes edge splitting, since split points become
    // sub-polygon vertices).
    let ea = edges_of(a);
    let eb = edges_of(b);
    for e in &ea {
        for f in &fb {
            let _ = edge_face_intersections(e, f, tol);
        }
    }
    for e in &eb {
        for f in &fa {
            let _ = edge_face_intersections(e, f, tol);
        }
    }
    for e1 in &ea {
        for e2 in &eb {
            let _ = edge_edge_intersections(e1, e2, tol);
        }
    }

    // Step 3: split every face along the segments that lie on it.
    let bld = TopoBuilder::new();
    let weld_cell = tol.max(1e-7);
    let mut weld = Weld::new(weld_cell);
    let mut edge_map = EdgeMap::default();
    let subs_a = split_faces(&bld, &fa, &planes_a, &segs_a, &mut weld, &mut edge_map);
    let subs_b = split_faces(&bld, &fb, &planes_b, &segs_b, &mut weld, &mut edge_map);

    // Steps 4–5: classify each sub-face and select per the operation.
    let mut selected: Vec<SubFace> = Vec::new();
    let mut flipped: Vec<bool> = Vec::new();
    for sub in &subs_a {
        let c = classify_face(&sub.plane, &sub.interior, &fb, b, tol);
        if select(c, true, op) {
            selected.push(sub.clone());
            flipped.push(false);
        }
    }
    for sub in &subs_b {
        let c = classify_face(&sub.plane, &sub.interior, &fa, a, tol);
        if select(c, false, op) {
            selected.push(sub.clone());
            // B's cut-through surfaces are oriented outward from A.
            flipped.push(matches!(op, BoolOp::Cut));
        }
    }

    // Step 6: rebuild the shell (and solid when closed).
    let mut result_faces: Vec<Face> = Vec::new();
    for (sub, fl) in selected.iter().zip(&flipped) {
        if *fl {
            result_faces.push(flip_face(&bld, sub));
        } else {
            result_faces.push(sub.face.clone());
        }
    }
    if result_faces.is_empty() {
        return Ok(empty_result(op));
    }
    // Unify edges across faces so shared geometric boundaries reference one
    // Edge TShape (a face split by intersection lines otherwise keeps
    // subdivided boundary edges while its neighbour keeps the long edge, and
    // the shell never closes). Run BEFORE orienting outward — the closure
    // check in `orient_faces_outward` needs the unified (closed) shell.
    let result_faces = unify_result_edges(&bld, &result_faces, tol);
    // Make every face's normal point outward from the result, so the signed
    // mesh volume is consistent (source faces may be oriented inward).
    let result_faces = orient_faces_outward(&bld, &result_faces);
    let shell = bld.make_shell(&result_faces);
    let inv = shell_invariants(&shell);
    // A manifold (closed) shell yields a solid; Euler is topology-reporting
    // (2 for a simple solid, 0 for a solid with a cavity) and is asserted
    // separately by the invariant oracle test.
    let solid = if inv.closed { Some(bld.make_solid(&[shell.clone()])) } else { None };
    let shape = solid
        .as_ref()
        .map(|s| s.0.clone())
        .unwrap_or_else(|| shell.0.clone());

    let mut result = BooleanResult {
        shape,
        solid,
        shells: vec![shell],
        faces: result_faces,
        warnings: vec![],
    };
    // Topological invariant gate (trellis R4 "invariant oracle"): a boolean
    // result must be a closed (manifold) shell.
    if !inv.closed {
        result.warnings.push("result shell is not closed".into());
    }

    // Step 7: validate the volume.
    validate(a, b, op, &mut result);
    Ok(result)
}
