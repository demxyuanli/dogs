//! Geometric construction of one edge from a mergeable chain.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:1665-1750`,
//! `:2164-2499`, and `:2864-2905`.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpAx2, GpDir, GpLin, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::convert_bspl::{curve_to_bspline_curve, CompCurveToBSplineCurve};
use occt_geom::geom_convert_concat::{c0_bspline_to_c1_bspline_curve, concat_c1};
use occt_geom::{Curve, GeomCircle, GeomLine, GeomTrimmedCurve};
use occt_geom2d::{curve::Curve2d, Geom2dLine};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex};
use crate::tgeometry::{GeometryRegistry, VertexGeom};
use crate::topo_tools::{common_vertex, edge_vertices, is_same};

use super::context::UnifyContext;
use super::edge_chain::{find_closest_points, oriented_vertices, AncestorMap, SubSequenceOfEdges};

fn safe_vertex(
    vertex: Vertex,
    context: &mut UnifyContext,
    safe_input_mode: bool,
) -> Result<Vertex, String> {
    if !safe_input_mode {
        return Ok(vertex);
    }
    if context.is_recorded(&vertex.0) {
        return Vertex::wrap(context.apply(&vertex.0))
            .ok_or_else(|| "MergeSubSeq: vertex replacement is not a vertex".into());
    }
    let copy = Vertex(vertex.0.empty_copied());
    GeometryRegistry::global().set_vertex(
        &copy.0,
        VertexGeom {
            point: BRepTool::vertex_point(&vertex),
            tolerance: BRepTool::vertex_tolerance(&vertex),
        },
    );
    context.replace(&vertex.0, &copy.0)?;
    Ok(copy)
}

fn attach_vertices(edge: &mut Edge, first: &Vertex, last: &Vertex, first_p: f64, last_p: f64) {
    let builder = TopoBuilder::new();
    builder.add_edge_vertices(edge, first, last);
    let registry = GeometryRegistry::global();
    registry.set_vertex_parameter_on_edge(&edge.0, &first.0, first_p);
    registry.set_vertex_parameter_on_edge(&edge.0, &last.0, last_p);
}

/// `GlueEdgesWith3DCurves` (`cxx:1665-1750`).
fn glue_edges_with_3d_curves(
    chain: &[Edge],
    first_vertex: &Vertex,
    last_vertex: &Vertex,
) -> Result<Edge, String> {
    let mut previous_vertex = first_vertex.clone();
    let mut previous_edge = chain
        .first()
        .cloned()
        .ok_or_else(|| "GlueEdgesWith3DCurves: empty chain".to_string())?;
    let mut max_tolerance: f64 = 0.0;
    let mut splines = Vec::with_capacity(chain.len());
    let mut tolerances = Vec::with_capacity(chain.len().saturating_sub(1));
    for (index, edge) in chain.iter().enumerate() {
        let (first, last) = edge_vertices(edge);
        let first = first.ok_or("GlueEdgesWith3DCurves: null first vertex")?;
        let last = last.ok_or("GlueEdgesWith3DCurves: null last vertex")?;
        let reverse = !is_same(&first.0, &previous_vertex.0);
        max_tolerance = max_tolerance
            .max(BRepTool::vertex_tolerance(&first))
            .max(BRepTool::vertex_tolerance(&last));
        if index > 0 {
            let common = common_vertex(&previous_edge, edge)
                .and_then(Vertex::wrap)
                .ok_or("GlueEdgesWith3DCurves: edges have no common vertex")?;
            tolerances.push(BRepTool::vertex_tolerance(&common));
        }
        let curve = BRepTool::edge_curve_adaptor_world(edge)
            .ok_or("GlueEdgesWith3DCurves: edge has no 3D curve")?;
        let (first_parameter, last_parameter) = BRepTool::edge_parameters(edge);
        let trimmed = GeomTrimmedCurve::new(curve, first_parameter, last_parameter);
        let spline = curve_to_bspline_curve(&trimmed, ParameterisationType::TgtThetaOver2)
            .map_err(|error| format!("GlueEdgesWith3DCurves: CurveToBSplineCurve ({error:?})"))?;
        let mut spline =
            c0_bspline_to_c1_bspline_curve(&spline, CONFUSION, occt_core::precision::ANGULAR)?;
        if reverse {
            spline.reverse();
        }
        splines.push(spline);
        previous_vertex = if reverse { first } else { last };
        previous_edge = edge.clone();
    }
    let mut closed = false;
    let (_, mut concatenated) = concat_c1(&mut splines, &tolerances, &mut closed, CONFUSION)?;
    let first = concatenated
        .drain(..1)
        .next()
        .ok_or("GlueEdgesWith3DCurves: concatenation returned no curve")?;
    let mut compound =
        CompCurveToBSplineCurve::from_bspline(first, ParameterisationType::TgtThetaOver2);
    for spline in &concatenated {
        if !compound
            .add(spline, max_tolerance, true, true, 0)
            .map_err(|error| format!("GlueEdgesWith3DCurves: concatenate ({error:?})"))?
        {
            return Err("GlueEdgesWith3DCurves: disconnected concatenated curves".into());
        }
    }
    let result = compound
        .into_curve()
        .ok_or("GlueEdgesWith3DCurves: missing result curve")?;
    let first_parameter = result.first_parameter();
    let last_parameter = result.last_parameter();
    let mut edge = TopoBuilder::new().make_edge(Arc::new(result), first_parameter, last_parameter);
    attach_vertices(
        &mut edge,
        first_vertex,
        last_vertex,
        first_parameter,
        last_parameter,
    );
    let registry = GeometryRegistry::global();
    registry.set_same_range(&edge.0, false);
    registry.set_same_parameter(&edge.0, false);
    Ok(edge)
}

fn merge_degenerated(chain: &[Edge], vertex_faces: &AncestorMap) -> Result<Option<Edge>, String> {
    let first_edge = chain.first().unwrap();
    let last_edge = chain.last().unwrap();
    let Some(points) = find_closest_points(first_edge, last_edge, vertex_faces) else {
        return Ok(None);
    };
    let first_index = 1 - points.first_index;
    let second_index = 1 - points.second_index;
    let mut start = points.first_points[first_index];
    let mut end = points.second_points[second_index];
    if (points.first_orientation == Orientation::Forward && first_index == 1)
        || (points.first_orientation == Orientation::Reversed && first_index == 0)
    {
        std::mem::swap(&mut start, &mut end);
    }
    let vector = occt_core::gp::GpVec2d::new(end.x() - start.x(), end.y() - start.y());
    let direction = match occt_core::gp::GpDir2d::from_vec2d(&vector) {
        Ok(direction) => direction,
        Err(_) => return Ok(None),
    };
    let line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(start, direction));
    let vertex = oriented_vertices(first_edge)
        .0
        .ok_or("MergeSubSeq: degenerated edge has no vertex")?;
    let point = BRepTool::vertex_point_world(&vertex);
    let mut edge = TopoBuilder::new().make_degenerated_edge(&point);
    let distance = start.distance(&end);
    GeometryRegistry::global().set_edge_range(&edge.0, 0.0, distance);
    GeometryRegistry::global().set_edge_pcurve(
        &edge.0,
        GeometryRegistry::shape_key(&points.face.0),
        line,
    );
    let mut first_vertex = vertex.clone();
    first_vertex.0.set_orientation(Orientation::Forward);
    let mut last_vertex = vertex;
    last_vertex.0.set_orientation(Orientation::Reversed);
    TopoBuilder::new().add(&mut edge.0, &first_vertex.0);
    TopoBuilder::new().add(&mut edge.0, &last_vertex.0);
    Ok(Some(edge))
}

fn merge_lines(
    chain: &[Edge],
    context: &mut UnifyContext,
    safe_input_mode: bool,
    edge_faces: &AncestorMap,
    face_plane_map: &super::same_domain::FacePlaneMap,
    face_new_face: &std::collections::HashMap<crate::topo_tools::ShapeMapKey, crate::shape::Face>,
) -> Result<Edge, String> {
    let (first, _) = oriented_vertices(chain.first().unwrap());
    let (_, last) = oriented_vertices(chain.last().unwrap());
    let first = safe_vertex(
        first.ok_or("MergeSubSeq: line has no first vertex")?,
        context,
        safe_input_mode,
    )?;
    let last = safe_vertex(
        last.ok_or("MergeSubSeq: line has no last vertex")?,
        context,
        safe_input_mode,
    )?;
    let first_point = BRepTool::vertex_point_world(&first);
    let last_point = BRepTool::vertex_point_world(&last);
    let vector = GpVec::from_pnts(&first_point, &last_point);
    let direction = GpDir::from_vec(&vector).map_err(|_| "MergeSubSeq: zero-length line union")?;
    let distance = first_point.distance(&last_point);
    let curve: Arc<dyn Curve> =
        Arc::new(GeomLine::new(GpLin::from_pnt_dir(first_point, direction)));
    let mut edge = TopoBuilder::new().make_edge(curve, 0.0, distance);
    attach_vertices(&mut edge, &first, &last, 0.0, distance);
    super::union_pcurves::union_pcurves(chain, &edge, edge_faces, face_plane_map, face_new_face)?;
    Ok(edge)
}

fn merge_circles(
    chain: &[Edge],
    context: &mut UnifyContext,
    safe_input_mode: bool,
    edge_faces: &AncestorMap,
    face_plane_map: &super::same_domain::FacePlaneMap,
    face_new_face: &std::collections::HashMap<crate::topo_tools::ShapeMapKey, crate::shape::Face>,
) -> Result<Option<Edge>, String> {
    let first_edge = chain.first().unwrap();
    let curve =
        BRepTool::edge_curve_adaptor_world(first_edge).ok_or("MergeSubSeq: circle has no curve")?;
    let circle = curve
        .gp_circ()
        .ok_or("MergeSubSeq: curve is not a circle")?;
    let (first, _) = oriented_vertices(first_edge);
    let (_, last) = oriented_vertices(chain.last().unwrap());
    let mut first = first.ok_or("MergeSubSeq: circle has no first vertex")?;
    let mut last = last.ok_or("MergeSubSeq: circle has no last vertex")?;
    let mut closed = is_same(&first.0, &last.0);
    if !closed {
        let first_point = BRepTool::vertex_point_world(&first);
        let last_point = BRepTool::vertex_point_world(&last);
        let tolerance = BRepTool::vertex_tolerance(&first).max(BRepTool::vertex_tolerance(&last));
        if first_point.square_distance(&last_point) < tolerance * tolerance {
            closed = true;
            last = first.clone();
            last.0.reverse();
        }
    }
    let mut edge;
    if closed {
        edge = TopoBuilder::new().make_edge(
            Arc::new(GeomCircle::new(circle)),
            0.0,
            2.0 * std::f64::consts::PI,
        );
        attach_vertices(&mut edge, &first, &last, 0.0, 2.0 * std::f64::consts::PI);
        edge.0.set_orientation(first_edge.0.orientation());
    } else {
        first = safe_vertex(first, context, safe_input_mode)?;
        last = safe_vertex(last, context, safe_input_mode)?;
        let first_point = BRepTool::vertex_point_world(&first);
        let last_point = BRepTool::vertex_point_world(&last);
        let origin = circle.location();
        let first_direction = GpDir::from_vec(&GpVec::from_pnts(&origin, &first_point))
            .map_err(|_| "MergeSubSeq: circle first point is at center")?;
        let mut normal = *circle.axis().direction();
        let (first_parameter, last_parameter) = BRepTool::edge_parameters(first_edge);
        let tangent_parameter = if first_edge.0.orientation() == Orientation::Forward {
            (first_parameter + last_parameter) * 0.5
        } else {
            (first_parameter + last_parameter) * 0.5
        };
        let (_, mut tangent) = curve.d1(tangent_parameter);
        if first_edge.0.orientation() == Orientation::Reversed {
            tangent.reverse();
        }
        let first_vector = GpVec::new(
            first_direction.x(),
            first_direction.y(),
            first_direction.z(),
        );
        let normal_vector = GpVec::new(normal.x(), normal.y(), normal.z());
        if first_vector.crossed(&tangent).dot(&normal_vector) < 0.0 {
            normal.reverse();
        }
        let axis = GpAx2::new(origin, normal, first_direction)
            .map_err(|_| "MergeSubSeq: invalid circle frame")?;
        let new_circle = occt_core::gp::GpCirc::new(axis, circle.radius());
        let last_direction = GpDir::from_vec(&GpVec::from_pnts(&origin, &last_point))
            .map_err(|_| "MergeSubSeq: circle last point is at center")?;
        let mut last_parameter = first_direction.angle_with_ref(&last_direction, &normal);
        if last_parameter < 0.0 {
            last_parameter += 2.0 * std::f64::consts::PI;
        }
        edge = TopoBuilder::new().make_edge(
            Arc::new(GeomCircle::new(new_circle)),
            0.0,
            last_parameter,
        );
        attach_vertices(&mut edge, &first, &last, 0.0, last_parameter);
    }
    super::union_pcurves::union_pcurves(chain, &edge, edge_faces, face_plane_map, face_new_face)?;
    Ok(Some(edge))
}

/// `MergeSubSeq` (`cxx:2164-2499`).
pub(super) fn merge_subsequence(
    chain: &[Edge],
    vertex_faces: &AncestorMap,
    context: &mut UnifyContext,
    safe_input_mode: bool,
    concat_bsplines: bool,
    angular_tolerance: f64,
    edge_faces: &AncestorMap,
    face_plane_map: &super::same_domain::FacePlaneMap,
    face_new_face: &std::collections::HashMap<crate::topo_tools::ShapeMapKey, crate::shape::Face>,
) -> Result<Option<Edge>, String> {
    if chain.len() < 2 {
        return Ok(None);
    }
    let all_degenerated = chain.iter().all(BRepTool::is_degenerated);
    if all_degenerated {
        return merge_degenerated(chain, vertex_faces);
    }
    if chain.iter().any(BRepTool::is_degenerated) {
        return Ok(None);
    }

    let mut lines_possible = true;
    let mut circles_possible = true;
    for pair in chain.windows(2) {
        let first = BRepTool::edge_curve_adaptor_world(&pair[0])
            .ok_or("MergeSubSeq: edge has no 3D curve")?;
        let second = BRepTool::edge_curve_adaptor_world(&pair[1])
            .ok_or("MergeSubSeq: edge has no 3D curve")?;
        let first_range = BRepTool::edge_parameters(&pair[0]);
        let second_range = BRepTool::edge_parameters(&pair[1]);
        match (
            super::edge_chain::linear_direction(first.as_ref(), first_range.0, first_range.1),
            super::edge_chain::linear_direction(second.as_ref(), second_range.0, second_range.1),
        ) {
            (Some(first_direction), Some(second_direction))
                if first_direction.is_parallel_tol(&second_direction, angular_tolerance) => {}
            _ => lines_possible = false,
        }
        match (first.gp_circ(), second.gp_circ()) {
            (Some(first_circle), Some(second_circle))
                if first_circle.location().distance(&second_circle.location()) <= CONFUSION => {}
            _ => circles_possible = false,
        }
    }
    if lines_possible && circles_possible {
        return Ok(None);
    }
    if lines_possible {
        return merge_lines(
            chain,
            context,
            safe_input_mode,
            edge_faces,
            face_plane_map,
            face_new_face,
        )
        .map(Some);
    }
    if circles_possible {
        return merge_circles(
            chain,
            context,
            safe_input_mode,
            edge_faces,
            face_plane_map,
            face_new_face,
        );
    }
    if concat_bsplines
        && chain.iter().all(|edge| {
            BRepTool::edge_curve_adaptor_world(edge).is_some_and(|curve| {
                curve.bspline_poles().is_some() || curve.bezier_poles().is_some()
            })
        })
    {
        let (first, _) = oriented_vertices(chain.first().unwrap());
        let (_, last) = oriented_vertices(chain.last().unwrap());
        let first = first.ok_or("MergeSubSeq: concatenated curve has no first vertex")?;
        let last = last.ok_or("MergeSubSeq: concatenated curve has no last vertex")?;
        let edge = glue_edges_with_3d_curves(chain, &first, &last)?;
        return Ok(Some(edge));
    }
    Ok(None)
}

/// `MergeEdges` result construction and `MergeSeq` context merge
/// (`cxx:2864-2905`).
pub(super) fn merge_sequences(
    subsequences: &mut [SubSequenceOfEdges],
    vertex_faces: &AncestorMap,
    context: &mut UnifyContext,
    safe_input_mode: bool,
    concat_bsplines: bool,
    angular_tolerance: f64,
    edge_faces: &AncestorMap,
    face_plane_map: &super::same_domain::FacePlaneMap,
    face_new_face: &std::collections::HashMap<crate::topo_tools::ShapeMapKey, crate::shape::Face>,
) -> Result<bool, String> {
    for subsequence in subsequences.iter_mut() {
        if subsequence.edges.len() < 2 {
            continue;
        }
        subsequence.union_edge = merge_subsequence(
            &subsequence.edges,
            vertex_faces,
            context,
            safe_input_mode,
            concat_bsplines,
            angular_tolerance,
            edge_faces,
            face_plane_map,
            face_new_face,
        )?;
    }
    for subsequence in subsequences {
        let Some(union) = &subsequence.union_edge else {
            continue;
        };
        let parts: Vec<_> = subsequence
            .edges
            .iter()
            .map(|edge| edge.0.clone())
            .collect();
        context.merge(&parts, &union.0)?;
    }
    Ok(true)
}
