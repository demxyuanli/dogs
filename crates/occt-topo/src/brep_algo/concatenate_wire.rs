//! `BRepAlgo::ConcatenateWire` (`BRepAlgo.cxx:245-400`).

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::curve::Curve;
use occt_geom::geom_convert_concat::{c0_bspline_to_c1_bspline_curve, concat_c1, concat_g1};
use occt_geom::geom_lprop::continuity;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex, Wire};
use crate::topo_tools::edge_vertices;

/// `BRepLib_MakeEdge(C)` (`BRepLib_MakeEdge.cxx:604-610`): `Init(C)` over the
/// curve's own bounds, i.e. `Init(C, null, null, p1, p2)`.
fn make_edge_of_curve(curve: &Arc<dyn Curve>) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d(curve, curve.first_parameter(), curve.last_parameter())?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(Arc::clone(curve), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}

/// `BRepAlgo::ConcatenateWire(W, Option, TolAngular)` (`cxx:245-400`).
///
/// Walks `W` with `BRepTools_WireExplorer`, turns each edge's `Geom_Curve` into
/// a B-spline over the edge range (`cxx:286-290`), fixes its C0 knots
/// (`cxx:290`, `GeomConvert::C0BSplineToC1BSplineCurve` with
/// `Precision::Confusion()`), reverses the REVERSED edges of a junction
/// (`cxx:293-296`) and then merges with `GeomConvert::ConcatG1` for
/// `GeomAbs_G1` and `GeomConvert::ConcatC1` otherwise (`cxx:354-370`).
///
/// `option` is `GeomAbs_G1` for the G1 arm; every other value takes the C1 arm,
/// exactly as the `if (Option == GeomAbs_G1) ... else ...` of `cxx:354` does.
pub fn concatenate_wire(w: &Wire, option: Shape, tol_angular: f64) -> Result<Wire, String> {
    let mut wexp = WireExplorer::of_wire(w);

    // `cxx:265-270`: curve count.
    let mut nb_curve = 0usize;
    while wexp.more() {
        nb_curve += 1;
        wexp.next();
    }

    let mut make_result = MakeWire::new();

    if nb_curve > 1 {
        let mut tab: Vec<GeomBSplineCurve> = Vec::with_capacity(nb_curve);
        let mut tab_tol_vertex = vec![0.0f64; nb_curve.saturating_sub(1)];
        let mut first0 = 0.0f64;
        let mut last = 0.0f64;
        let mut v_first: Option<Vertex> = None;
        let mut v_last: Option<Vertex> = None;

        wexp.init(w, None);
        for index in 0..nb_curve {
            let edge = wexp
                .current_edge()
                .ok_or_else(|| "BRepAlgo ConcatenateWire: no edge".to_string())?;
            // `cxx:285-290`: `Geom_TrimmedCurve(aCurve, First, Last)` then
            // `CurveToBSplineCurve` + `Transform(L.Transformation())` then the
            // C0 -> C1 knot fix.
            let curve = BRepTool::edge_curve(&edge)
                .ok_or_else(|| "BRepAlgo ConcatenateWire: null curve".to_string())?;
            let (first, edge_last) = BRepTool::edge_parameters(&edge);
            let trimmed = GeomTrimmedCurve::new(curve, first, edge_last);
            let mut bs = curve_to_bspline_curve(&trimmed, ParameterisationType::TgtThetaOver2)
                .map_err(|e| format!("BRepAlgo ConcatenateWire: CurveToBSplineCurve ({e:?})"))?;
            bs.transform(&edge.0.location().transformation());
            let bs = c0_bspline_to_c1_bspline_curve(&bs, CONFUSION, ANGULAR)?;

            if index >= 1 {
                // `cxx:291-299`: the junction tolerances, and the REVERSED
                // edge flipped so that consecutive curves run the same way.
                let mut bs = bs;
                if edge.0.orientation() == Orientation::Reversed {
                    bs.reverse();
                }
                let (vf, vl) = edge_vertices(&edge);
                let tolleft = vl
                    .as_ref()
                    .map(BRepTool::vertex_tolerance)
                    .ok_or_else(|| "BRepAlgo ConcatenateWire: null last vertex".to_string())?;
                let tolright = vf
                    .as_ref()
                    .map(BRepTool::vertex_tolerance)
                    .ok_or_else(|| "BRepAlgo ConcatenateWire: null first vertex".to_string())?;
                tab_tol_vertex[index - 1] = tolleft.max(tolright);
                tab.push(bs);
            } else {
                // `cxx:302-315`: the first edge's features; its own orientation
                // decides `Vfirst` and whether the curve is flipped.
                first0 = first;
                let mut bs = bs;
                let (vf, vl) = edge_vertices(&edge);
                if edge.0.orientation() == Orientation::Reversed {
                    v_first = vl;
                    bs.reverse();
                } else {
                    v_first = vf;
                }
                tab.push(bs);
            }

            if index == nb_curve - 1 {
                // `cxx:317-326`: the last edge's features.
                let (vf, vl) = edge_vertices(&edge);
                v_last = if edge.0.orientation() == Orientation::Reversed {
                    vf
                } else {
                    vl
                };
            }

            last = edge_last;
            wexp.next();
        }

        let v_first = v_first
            .ok_or_else(|| "BRepAlgo ConcatenateWire: null first vertex".to_string())?;
        let v_last =
            v_last.ok_or_else(|| "BRepAlgo ConcatenateWire: null last vertex".to_string())?;

        // `cxx:329-337`: the closing tolerance.
        let toler = BRepTool::vertex_tolerance(&v_first).max(BRepTool::vertex_tolerance(&v_last));

        // `cxx:339-351`: the C0 / G1 closing test.
        let mut closed_flag = false;
        let mut closed_tolerance = 0.0f64;
        let pfirst = BRepTool::vertex_point_world(&v_first);
        let plast = BRepTool::vertex_point_world(&v_last);
        if pfirst.distance(&plast) <= toler {
            let cont = continuity(
                &tab[nb_curve - 1],
                &tab[0],
                last,
                first0,
                true,
                true,
                toler,
                tol_angular,
            )?;
            if cont >= Shape::G1 {
                closed_tolerance = toler;
                closed_flag = true;
            }
        }

        let concatenated = if option == Shape::G1 {
            // `cxx:356-360`.
            concat_g1(&mut tab, &tab_tol_vertex, &mut closed_flag, closed_tolerance)?
        } else {
            // `cxx:364-369`.
            let (_indices, curves) =
                concat_c1(&mut tab, &tab_tol_vertex, &mut closed_flag, closed_tolerance)?;
            curves
        };

        // `cxx:372-377`: one edge per concatenated curve.
        for bs in &concatenated {
            let arc: Arc<dyn Curve> = Arc::new(bs.clone());
            let edge = make_edge_of_curve(&arc)
                .ok_or_else(|| "BRepAlgo ConcatenateWire: BRepLib_MakeEdge not done".to_string())?;
            let _ = make_result.add(&edge);
        }
    } else {
        // `cxx:380-398`: a wire of at most one edge keeps that edge's curve.
        wexp.init(w, None);
        let edge = wexp
            .current_edge()
            .ok_or_else(|| "BRepAlgo ConcatenateWire: empty wire".to_string())?;
        let curve = BRepTool::edge_curve(&edge)
            .ok_or_else(|| "BRepAlgo ConcatenateWire: null curve".to_string())?;
        let (first, edge_last) = BRepTool::edge_parameters(&edge);
        let trimmed = GeomTrimmedCurve::new(curve, first, edge_last);
        let mut bs = curve_to_bspline_curve(&trimmed, ParameterisationType::TgtThetaOver2)
            .map_err(|e| format!("BRepAlgo ConcatenateWire: CurveToBSplineCurve ({e:?})"))?;
        bs.transform(&edge.0.location().transformation());
        let mut bs = c0_bspline_to_c1_bspline_curve(&bs, CONFUSION, ANGULAR)?;
        if edge.0.orientation() == Orientation::Reversed {
            bs.reverse();
        }
        let arc: Arc<dyn Curve> = Arc::new(bs);
        let edge = make_edge_of_curve(&arc)
            .ok_or_else(|| "BRepAlgo ConcatenateWire: BRepLib_MakeEdge not done".to_string())?;
        let _ = make_result.add(&edge);
    }

    Ok(make_result.wire())
}
