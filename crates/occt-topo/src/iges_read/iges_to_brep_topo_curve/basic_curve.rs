use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_core::toploc::TopLocLocation;

use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_basic_curve::{
    convert_location, edge_from_span, entity_location, CurveSpan,
};
use super::super::iges_to_brep_curve_and_surface::Session;

use super::CONTINUITY;

/// `TransferTopoBasicCurve` (`IGESToBRep_TopoCurve.cxx:950-1123`). One edge,
/// except that a C0 B-spline is replaced by its C1 pieces (`cxx:1013-1038`);
/// the result is then that single edge (`cxx:1108-1110`) or the wire of the
/// pieces (`cxx:1099-1104`).
///
/// `TheCurves` holds exactly the one `TransferBasicCurve` result: the
/// `GetModeApprox()` arm (`cxx:979-991`) is off because `myModeApprox` is set
/// false by the constructor and reset false by `Init`
/// (`IGESToBRep_CurveAndSurface.cxx:64`, `:110`), so `ApproxBSplineCurve` is
/// unreachable and unported.
pub(super) fn transfer_topo_basic_curve(span: &CurveSpan) -> Result<TopoShape, String> {
    // `cxx:1015-1018`: `IsKind(STANDARD_TYPE(Geom_BSplineCurve))` &&
    // `Continuity() < GeomAbs_C1` && `GetContinuity() >= 1`.
    let mut pieces = None;
    if CONTINUITY >= 1
        && occt_geom::convert_bspl::is_bspline(span.curve.as_ref())
        && span.curve.continuity() < crate::shape_upgrade::GEOMABS_C1
    {
        // `cxx:1021-1022`:
        // `ShapeAlgo::AlgoContainer()->C0BSplineToSequenceOfC1BSplineCurve`.
        pieces = crate::shape_upgrade::c0_bspline_to_sequence_of_c1_bspline_curve(
            span.curve.as_ref(),
        );
    }
    let Some(pieces) = pieces else {
        // `cxx:1039-1053`: `BRepBuilderAPI_MakeEdge ME(mycurve, a, b)`.
        return edge_from_span(span);
    };

    // `cxx:1023-1037`: `BRepBuilderAPI_MakeEdge ME(seqBS->Value(i))` ->
    // `BRepLib_MakeEdge::Init(C)` -> `Init(C, FirstParameter, LastParameter)`
    // (`BRepLib_MakeEdge.cxx:505-508`).
    let mut edges = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let pspan = CurveSpan {
            first: piece.first_parameter(),
            last: piece.last_parameter(),
            curve: piece,
        };
        edges.push(Edge(edge_from_span(&pspan)?));
    }
    if edges.is_empty() {
        return Err(
            "IGES_1156 TopoBasicCurve produced no curve: IGESToBRep_TopoCurve.cxx:1000-1008".into(),
        );
    }
    if edges.len() == 1 {
        // `cxx:1108-1110`: `myedge = TopoDS::Edge(TopoDS_Iterator(myshape).Value())`.
        return Ok(edges.pop().expect("one edge").0);
    }
    // `cxx:1059-1064`: `ShapeFix_Wire sfw; sfw.Load(sewd); sfw.FixConnected();
    // myshape = sewd->Wire()`. `myClosedMode` is true
    // (`ShapeFix_Wire.cxx:117-126`, `:167-172`), so every junction including
    // the last-to-first one is visited, and the default `prec = -1.` resolves
    // to `MaxTolerance() = Precision::Confusion()` (`ShapeFix_Wire.cxx:1485`,
    // `ShapeFix_Root.cxx:24-28`).
    let mut wire = TopoBuilder::new().make_wire(&edges);
    let _ = crate::shhealing::fix_connected_closed(&mut wire, CONFUSION, CONFUSION, true);
    // `cxx:1099-1104`: more than one edge -> the wire (plus the `IGES_1120`
    // warning).
    Ok(wire.0)
}

pub(super) fn transfer_point(session: &Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `TransferPoint` (`IGESToBRep_TopoCurve.cxx:129-171`).
    // `GetModeTransfer` is true (`IGESToBRep_Actor.cxx:167`), so the point is
    // `Value()`, then scaled, and a directory 124 is `Move`d (`cxx:157-168`).
    let xyz = match &session.model.entity(index)?.data {
        EntityData::Point { xyz } => *xyz,
        _ => return Err("IGES 116 has no point parameters".into()),
    };
    let unit = session.model.unit_factor;
    let point = GpPnt::new(xyz[0] * unit, xyz[1] * unit, xyz[2] * unit);
    let mut shape: TopoShape = TopoBuilder::new().make_vertex(point, CONFUSION).into();
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, unit) {
            shape.move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    }
    Ok(shape)
}
