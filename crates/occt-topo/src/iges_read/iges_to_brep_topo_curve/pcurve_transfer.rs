use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::builder::TopoBuilder;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_curve_and_surface::Session;

use super::pcurve_edge::{attach_pcurve, mapped_pcurve};
use super::CONTINUITY;

pub(in crate::iges_read) fn transfer_pcurve(
    from: &Edge,
    to: &Edge,
    face: &TopoShape,
) -> bool {
    // `IGESToBRep::TransferPCurve` (`IGESToBRep.cxx:400-476`) onto an edge
    // that does not already carry a pcurve on `face`.
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let pcs = reg.edge_pcurves(&from.0, key);
    let Some(pc) = pcs.into_iter().next() else {
        return false;
    };
    let (first, last) = reg
        .pcurve_range(&from.0, key)
        .unwrap_or_else(|| reg.edge_parameters(&from.0));
    attach_pcurve(to, face, pc, first, last);
    true
}

pub(in crate::iges_read) fn curve_sequence(
    model: &super::super::iges_data_iges_model::IgesModel,
    index: usize,
) -> Result<Vec<usize>, String> {
    // `IGESToBRep::IGESCurveToSequenceOfIGESCurve`.
    match &model.entity(index)?.data {
        EntityData::Composite { curves } => {
            let mut out = Vec::new();
            for child in curves.clone() {
                out.extend(curve_sequence(model, child)?);
            }
            Ok(out)
        }
        _ => Ok(vec![index]),
    }
}

pub(in crate::iges_read) fn edges_from_2d(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &occt_core::gp::GpTrsf2d,
    u_fact: f64,
) -> Result<Vec<Edge>, String> {
    // `Transfer2dTopoCurve` (`cxx:910-946`). Entity 130 is
    // `Transfer2dOffsetCurve` (`cxx:941-943`), which may be a wire.
    if session.model.entity(index)?.type_number == 130 {
        let shape = super::super::offset_curve::transfer_2d_offset_curve(
            session, index, face, trans, u_fact,
        )?;
        return edges_of(&shape);
    }
    edges_from_2d_basic(session, index, face, trans, u_fact)
}

/// `Transfer2dTopoBasicCurve` (`IGESToBRep_TopoCurve.cxx:1127-1309`): the
/// single pcurve edge, except that a C0 B-spline pcurve becomes one edge per
/// C1 piece (`cxx:1232-1256`).
pub(in crate::iges_read) fn edges_from_2d_basic(
    session: &Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &occt_core::gp::GpTrsf2d,
    u_fact: f64,
) -> Result<Vec<Edge>, String> {
    let (pc, first, last) = mapped_pcurve(session, index, face, trans, u_fact)?;
    // `cxx:1232-1237`: `IsKind(STANDARD_TYPE(Geom2d_BSplineCurve))` &&
    // `Continuity() < GeomAbs_C1` && `GetContinuity() >= 1`.
    let mut pieces = None;
    if CONTINUITY >= 1
        && pc.is_bspline2d()
        && pc.continuity() < crate::shape_upgrade::GEOMABS_C1
    {
        // `cxx:1238`: `ShapeAlgo::AlgoContainer()->C0BSplineToSequenceOfC1BSplineCurve`.
        pieces = crate::shape_upgrade::c0_bspline2d_to_sequence_of_c1_bspline_curve(pc.as_ref());
    }
    let Some(pieces) = pieces else {
        // `cxx:1257-1274`: `ShapeBuild_Edge().MakeEdge(myedge, mycurve2d, face, a, b)`.
        return Ok(vec![edge_from_pcurve(pc, first, last, face)?]);
    };
    // `cxx:1239-1253`: `ShapeBuild_Edge().MakeEdge(myedge, seqBS->Value(i), face)`
    // -> `BRepLib_MakeEdge(C2d, S)` -> `Init(C2d, S, FirstParameter, LastParameter)`.
    let mut edges = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let f = piece.first_parameter();
        let l = piece.last_parameter();
        edges.push(edge_from_pcurve(piece, f, l, face)?);
    }
    Ok(edges)
}

/// The edge tail of `Transfer2dTopoBasicCurve` (`cxx:1259-1274`) for one
/// pcurve span, shared by the split and unsplit arms.
fn edge_from_pcurve(
    pc: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
    face: &TopoShape,
) -> Result<Edge, String> {
    let surf = GeometryRegistry::global()
        .face_surface(face)
        .ok_or_else(|| "IGES face has no surface".to_string())?;
    // `ShapeBuild_Edge().MakeEdge(myedge, mycurve2d, face, a, b)`
    // (`IGESToBRep_TopoCurve.cxx:1263`) -> `BRepBuilderAPI_MakeEdge(C2d, S,
    // p1, p2)` -> `BRepLib_MakeEdge::Init(C2d, S, null, null, p1, p2)`
    // (`BRepLib_MakeEdge.cxx:905-1094`). `Err` is the `Transfer2dTopoBasicCurve`
    // failure arm for a null edge (`cxx:1265-1271`).
    let Some(build) = crate::brep_lib_make_edge::init_2d(&pc, &surf, first, last) else {
        return Err(
            "IGES_1005 edge construction error: IGESToBRep_TopoCurve.cxx:1265-1271".into(),
        );
    };
    let curve3d: Arc<dyn occt_geom::Curve> =
        Arc::new(CurveOnSurface::new(pc.clone(), surf, build.first, build.last));
    let b = TopoBuilder::new();
    let mut edge = b.make_edge(curve3d, build.first, build.last);
    // `ShapeBuild_Edge::MakeEdge(edge, C2d, S, L, p1, p2)`
    // (`ShapeBuild_Edge.cxx:851-870`): the edge carries the pcurve only, so
    // `BRep_Tool::Curve` is null. `ShapeFix_Edge::FixAddCurve3d`
    // (`ShapeFix_Edge.cxx:618-638`) rebuilds it from the pcurve.
    GeometryRegistry::global().set_curve3d_null(&edge.0, true);
    crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut edge, &build);
    // `cxx:1089-1092`: a swapped non-periodic range stores the edge REVERSED.
    if build.reverse {
        edge.0.set_orientation(crate::abs::Orientation::Reversed);
    }
    attach_pcurve(&edge, face, pc, build.first, build.last);
    Ok(edge)
}

pub(in crate::iges_read) fn edges_of(shape: &TopoShape) -> Result<Vec<Edge>, String> {
    if shape.is_edge() {
        return Ok(vec![Edge(shape.clone())]);
    }
    if shape.is_wire() {
        let kids = shape.tshape.read().expect("lock").children.clone();
        let mut edges = Vec::new();
        for kid in kids {
            if kid.is_edge() {
                edges.push(Edge(kid));
            } else {
                return Err("IGES wire child is not an edge".into());
            }
        }
        return Ok(edges);
    }
    Err("IGES curve transfer did not produce an edge or a wire".into())
}
