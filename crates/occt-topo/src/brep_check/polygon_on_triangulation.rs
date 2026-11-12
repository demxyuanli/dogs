//! `BRepCheck_Edge::CheckPolygonOnTriangulation`.
//!
//! Source: `BRepCheck_Edge.cxx:711-839`.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_geom::curve::Curve;

use crate::brep_tool::BRepTool;
use crate::shape::Edge;
use crate::tshape::PolygonOnTriangulationRepr;

use super::Status;

const NCONTROL: usize = 23;

/// Checks every `BRep_PolygonOnTriangulation` representation against the
/// edge's 3D curve. `Err` is the exception path caught by
/// `BRepCheck_Analyzer`, which records `CheckFail`.
pub(super) fn check_polygon_on_triangulation(edge: &Edge) -> Result<Status, ()> {
    let representations = BRepTool::edge_polygons_on_triangulation(edge);
    let Some(curve) = BRepTool::edge_curve(edge) else {
        return Ok(Status::NoError);
    };
    if representations.is_empty() {
        return Ok(Status::NoError);
    }

    let curve: Arc<dyn Curve> = if edge.0.location().is_identity() {
        curve
    } else {
        Arc::from(curve.transformed(&edge.0.location().transformation()))
    };
    let (first, last) = BRepTool::edge_parameters(edge);

    for representation in &representations {
        if check_representation(edge, curve.as_ref(), first, last, representation)? {
            return Ok(Status::InvalidPolygonOnTriangulation);
        }
    }
    Ok(Status::NoError)
}

fn check_representation(
    edge: &Edge,
    curve: &dyn Curve,
    first: f64,
    last: f64,
    representation: &PolygonOnTriangulationRepr,
) -> Result<bool, ()> {
    // A closed representation is checked through PolygonOnTriangulation2,
    // exactly as `cxx:766-769`, independently of the edge orientation.
    let polygon = representation
        .polygon2
        .as_ref()
        .unwrap_or(&representation.polygon);
    let triangulation = &representation.triangulation;
    let tolerance = polygon.deflection() + BRepTool::edge_tolerance(edge);
    let tolerance2 = tolerance * tolerance;

    if let Some(parameters) = &polygon.parameters {
        for (i, &parameter) in parameters.iter().enumerate() {
            let node_index = *polygon.nodes.get(i).ok_or(())?;
            let point = triangulation_node(triangulation, node_index)?;
            if curve.value(parameter).square_distance(&point) > tolerance2 {
                return Ok(true);
            }
        }
        return Ok(false);
    }

    let mut bounds = BndBox::new();
    for &node_index in &polygon.nodes {
        bounds.add_point(&triangulation_node(triangulation, node_index)?);
    }
    bounds.enlarge(tolerance);

    let step = (last - first) / NCONTROL as f64;
    let mut parameter = first;
    for _ in 1..NCONTROL {
        if bounds.is_out(&curve.value(parameter)) {
            return Ok(true);
        }
        parameter += step;
    }
    Ok(bounds.is_out(&curve.value(last)))
}

fn triangulation_node(
    triangulation: &occt_core::poly::Triangulation,
    one_based_index: usize,
) -> Result<occt_core::gp::GpPnt, ()> {
    let index = one_based_index.checked_sub(1).ok_or(())?;
    triangulation.nodes.get(index).copied().ok_or(())
}
