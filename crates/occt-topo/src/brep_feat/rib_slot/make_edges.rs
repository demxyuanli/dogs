//! The `BRepLib_MakeEdge` constructions `SlidingProfile` /
//! `NoSlidingProfile` use while they assemble the profile wire
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:1418-2670`).
//!
//! Three reference constructors appear there:
//!
//! * `BRepLib_MakeEdge(P1, P2)` (`BRepLib_MakeEdge.cxx:202-215`) -> [`edge_from_points`].
//! * `BRepLib_MakeEdge(C, P1, P2)` (`cxx:403-408` -> `cxx:522-541`) ->
//!   [`edge_from_curve_points`].
//! * `BRepLib_MakeEdge(C, V1, V2)` (`cxx:412-417` -> `cxx:598+`) is reached as
//!   `WW.Add(E)` with an existing edge, which [`super::MakeWire`] handles.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::curve::Curve;

use crate::brep_lib_make_edge;
use crate::builder::TopoBuilder;
use crate::shape::Edge;

/// `BRepLib_MakeEdge(P1, P2)` (`BRepLib_MakeEdge.cxx:202-215`).
///
/// `None` is `BRepLib_LineThroughIdenticPoints` (`cxx:205-209`); the caller
/// then adds nothing to the wire, as `BRepLib_MakeWire::Add` does for a null
/// edge.
pub(super) fn edge_from_points(p1: &GpPnt, p2: &GpPnt) -> Option<Edge> {
    let (curve, build) = brep_lib_make_edge::init_points(p1, p2)?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(curve, build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}

/// `BRepLib_MakeEdge(C, P1, P2)` (`BRepLib_MakeEdge.cxx:403-408`): the two
/// points are projected onto `C` (`Init(C, P1, P2)`, `cxx:522-541`).
///
/// `None` is `BRepLib_PointProjectionFailed` (`cxx:555` / `:566`).
pub(super) fn edge_from_curve_points(
    curve: &Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
) -> Option<Edge> {
    let build = brep_lib_make_edge::init_curve_points(curve, p1, p2)?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}
