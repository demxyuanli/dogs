//! `IGESToBRep_TopoSurface::TransferRuledSurface` (`cxx:386-695`).
//!
//! `ReparamCurve` / `reparamBSpline` are `cxx:335-382`. The C0 retransfer
//! (`cxx:440-463`) and the `shape2.Reverse()` before the type switch
//! (`cxx:518-521`) stay commented out in OCCT and are not ported.
//! `TC.SetContinuity(0)` (`cxx:403`) is the actor default (`myContinuity = 0`,
//! `IGESToBRep_CurveAndSurface.cxx:68`), so the `GetContinuity() >= 1` splits
//! in `IGESToBRep_TopoCurve.cxx:1016` are not entered.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::curve::Curve;

use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Wire};
use crate::shhealing::fix_reorder_wire_3d;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::edges_of_wire;

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::{convert_location, entity_location};
use super::iges_to_brep_curve_and_surface::{is_topo_curve, Session};
use super::iges_to_brep_topo_curve::transfer_topo_curve;

pub(in crate::iges_read) fn transfer_ruled_surface(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    let (first, second, dirflag, _dev_flag, ruled_by_parameter) = {
        let entity = session.model.entity(index)?;
        let ruled_by_parameter = entity.form_number != 0;
        match &entity.data {
            EntityData::RuledSurface {
                first,
                second,
                direction,
                dev_flag,
            } => (*first, *second, *direction, *dev_flag, ruled_by_parameter),
            _ => return Err("IGES 118 has no ruled-surface parameters".into()),
        }
    };
    if !is_topo_curve(session.model.entity(first)?.type_number) {
        return Err("XSTEP_148: IGESToBRep_TopoSurface.cxx:492-497".into());
    }
    if !is_topo_curve(session.model.entity(second)?.type_number) {
        return Err("XSTEP_149: IGESToBRep_TopoSurface.cxx:561-566".into());
    }
    let mut shape1 = transfer_topo_curve(session, first)?;
    let mut shape2 = transfer_topo_curve(session, second)?;
    let (nb1, mut wire1) = prepare_directrix(&mut shape1, false)?;
    let (nb2, mut wire2) = prepare_directrix(&mut shape2, dirflag == 1)?;
    let result = if nb1 != nb2 {
        if nb1 == 1 {
            wire1 = Some(edge_wire(&shape1)?);
        } else if nb2 == 1 {
            wire2 = Some(edge_wire(&shape2)?);
        }
        let wire1 = wire1.ok_or_else(|| "IGES_1255: IGESToBRep_TopoSurface.cxx:588".to_string())?;
        let wire2 = wire2.ok_or_else(|| "IGES_1255: IGESToBRep_TopoSurface.cxx:588".to_string())?;
        // `IsRuledByParameter` is the unused `const bool` of `HomoWires` (`cxx:428`).
        let (new1, new2) = homo_wires(&wire1, &wire2, ruled_by_parameter)?;
        crate::brep_fill::shell(&new1, &new2)?.0
    } else if nb1 != 1 {
        let wire1 = wire1.ok_or_else(|| "IGES_1255: IGESToBRep_TopoSurface.cxx:667".to_string())?;
        let wire2 = wire2.ok_or_else(|| "IGES_1255: IGESToBRep_TopoSurface.cxx:667".to_string())?;
        crate::brep_fill::shell(&wire1, &wire2)?.0
    } else {
        let mut edge1 = Edge::wrap(shape1).ok_or_else(|| {
            "TopoDS::Edge: IGESToBRep_TopoSurface.cxx:614".to_string()
        })?;
        let mut edge2 = Edge::wrap(shape2).ok_or_else(|| {
            "TopoDS::Edge: IGESToBRep_TopoSurface.cxx:614".to_string()
        })?;
        reparam_span_to_unit(&mut edge1)?;
        reparam_span_to_unit(&mut edge2)?;
        // `surf->Continuity() == GeomAbs_C0` (`cxx:659-662`) only sends IGES_1250.
        crate::brep_fill::face(&edge1, &edge2)?.0
    };
    apply_ruled_location(session, index, result)
}

/// `ReparamCurve` on an edge, or on every edge of a wire (`cxx:468-551`).
/// `dirflag == 1` reverses the edge after that, or reverses the wire through
/// `ShapeExtend_WireData` (`cxx:545-550`).
fn prepare_directrix(
    shape: &mut TopoShape,
    reverse_join: bool,
) -> Result<(i32, Option<Wire>), String> {
    if shape.is_edge() {
        let mut edge = Edge::wrap(shape.clone()).expect("edge");
        reparam_line(&mut edge)?;
        if reverse_join {
            edge.0.reverse();
        }
        *shape = edge.0;
        Ok((1, None))
    } else if shape.is_wire() {
        let wire = Wire::wrap(shape.clone()).expect("wire");
        let edges = edges_of_wire(&wire);
        for edge in &edges {
            let mut edge = edge.clone();
            reparam_line(&mut edge)?;
        }
        let wire = if reverse_join {
            reverse_wire_data(&wire)
        } else {
            wire
        };
        let count = edges_of_wire(&wire).len() as i32;
        *shape = wire.0.clone();
        Ok((count, Some(wire)))
    } else {
        Err("IGESToBRep_TopoSurface.cxx:485-489".into())
    }
}

fn edge_wire(shape: &TopoShape) -> Result<Wire, String> {
    let edge = Edge::wrap(shape.clone())
        .ok_or_else(|| "ShapeExtend_WireData::Add: IGESToBRep_TopoSurface.cxx:573".to_string())?;
    Ok(TopoBuilder::new().make_wire(&[edge]))
}

/// `ReparamCurve` (`cxx:363-382`): only a `Geom_Line` is rebuilt on `[0, 1]`.
fn reparam_line(edge: &mut Edge) -> Result<(), String> {
    let reg = GeometryRegistry::global();
    let curve = reg
        .edge_curve(&edge.0)
        .ok_or_else(|| "BRep_Tool::Curve".to_string())?;
    if !curve.is_line() {
        return Ok(());
    }
    let (first, last) = reg.edge_parameters(&edge.0);
    let owned = curve.clone_dyn();
    let bspline = reparam_bspline(Arc::from(owned), first, last)?;
    let Some(bspline) = bspline else {
        return Ok(());
    };
    reg.set_edge_curve3d(&edge.0, Arc::new(bspline), CONFUSION);
    reg.set_ranges_all(&edge.0, 0.0, 1.0);
    Ok(())
}

/// The one-edge span rebuild (`cxx:610-646`), including curves already on `[0, 1]`.
fn reparam_span_to_unit(edge: &mut Edge) -> Result<(), String> {
    let reg = GeometryRegistry::global();
    let curve = reg
        .edge_curve(&edge.0)
        .ok_or_else(|| "BRep_Tool::Curve".to_string())?;
    let (first, last) = reg.edge_parameters(&edge.0);
    if first.abs() <= PCONFUSION && (last - 1.0).abs() <= PCONFUSION {
        return Ok(());
    }
    let bspline = reparam_bspline(Arc::from(curve.clone_dyn()), first, last)?
        .ok_or_else(|| "GeomConvert::CurveToBSplineCurve".to_string())?;
    reg.set_edge_curve3d(&edge.0, Arc::new(bspline), CONFUSION);
    reg.set_ranges_all(&edge.0, 0.0, 1.0);
    Ok(())
}

/// `reparamBSpline` (`cxx:335-361`). A null conversion leaves the curve unchanged.
fn reparam_bspline(
    mut curve: Arc<dyn Curve>,
    first: f64,
    last: f64,
) -> Result<Option<GeomBSplineCurve>, String> {
    let mut bspline = if !is_bspline(curve.as_ref()) {
        // `Geom_TrimmedCurve` keeps the basis parameter (`BRepFill` uses the
        // same constructor). Not the `[0, 1]` remap in `GeomTrimmedCurve`.
        if curve.first_parameter() < first || curve.last_parameter() > last {
            curve = Arc::new(occt_geom::geom_fill::ParameterTrim::new(curve, first, last));
        }
        match curve_to_bspline_curve(curve.as_ref(), ParameterisationType::RationalC1) {
            Ok(curve) => curve,
            Err(_) => return Ok(None),
        }
    } else {
        let mut curve = bspline_copy(curve.as_ref())?;
        curve
            .segment(first, last, PCONFUSION)
            .map_err(|_| "Geom_BSplineCurve::Segment".to_string())?;
        curve
    };
    let (mut knots, _) = bspline.distinct_knots_and_mults();
    occt_core::bspl::knots::reparametrize(0.0, 1.0, &mut knots);
    bspline
        .set_knots(&knots)
        .map_err(|_| "Geom_BSplineCurve::SetKnots".to_string())?;
    Ok(Some(bspline))
}

fn is_bspline(curve: &dyn Curve) -> bool {
    !curve.is_geom_trimmed() && curve.bspline_poles().is_some()
}

fn bspline_copy(curve: &dyn Curve) -> Result<GeomBSplineCurve, String> {
    let poles = curve
        .bspline_poles()
        .ok_or_else(|| "Geom_BSplineCurve".to_string())?
        .to_vec();
    let knots = curve
        .bspline_knots()
        .ok_or_else(|| "Geom_BSplineCurve".to_string())?
        .to_vec();
    let degree = curve
        .nurbs_degree()
        .ok_or_else(|| "Geom_BSplineCurve".to_string())?;
    Ok(GeomBSplineCurve::from_parts(
        poles,
        curve.bspline_weights().map(|w| w.to_vec()),
        knots,
        degree,
        curve.is_periodic(),
    ))
}

/// `ShapeExtend_WireData::Reverse` (`ShapeExtend_WireData.cxx:483-506`).
fn reverse_wire_data(wire: &Wire) -> Wire {
    let mut edges = edges_of_wire(wire);
    for edge in &mut edges {
        edge.0.reverse();
    }
    edges.reverse();
    TopoBuilder::new().make_wire(&edges)
}

/// `ShapeAlgo_AlgoContainer::HomoWires` (`ShapeAlgo_AlgoContainer.cxx:423-631`).
/// The last argument is unused in OCCT.
fn homo_wires(wire_in1: &Wire, wire_in2: &Wire, _ruled_by_parameter: bool) -> Result<(Wire, Wire), String> {
    let mut wire1 = copy_wire(wire_in1);
    let mut wire2 = copy_wire(wire_in2);
    let _ = fix_reorder_wire_3d(&mut wire1);
    let _ = fix_reorder_wire_3d(&mut wire2);
    let edges1 = edges_of_wire(&wire1);
    let edges2 = edges_of_wire(&wire2);
    let n1 = edges1.len();
    let n2 = edges2.len();
    let mut length1 = 0.0;
    let mut length2 = 0.0;
    for edge in &edges1 {
        let (first, last) = GeometryRegistry::global().edge_parameters(&edge.0);
        length1 += last - first;
    }
    for edge in &edges2 {
        let (first, last) = GeometryRegistry::global().edge_parameters(&edge.0);
        length2 += last - first;
    }
    let epsilon = PCONFUSION * (length1 + length2);
    if n1 == 1 && n2 == 1 {
        return Ok((wire1, wire2));
    }
    if length1 < epsilon {
        let edge1 = edges1
            .first()
            .cloned()
            .ok_or_else(|| "ShapeAlgo_AlgoContainer::HomoWires".to_string())?;
        let repeated = vec![edge1; n2.max(1)];
        return Ok((TopoBuilder::new().make_wire(&repeated), wire2));
    }
    if length2 < epsilon {
        let edge2 = edges2
            .first()
            .cloned()
            .ok_or_else(|| "ShapeAlgo_AlgoContainer::HomoWires".to_string())?;
        let repeated = vec![edge2; n1.max(1)];
        return Ok((wire1, TopoBuilder::new().make_wire(&repeated)));
    }
    let ratio = length2 / length1;
    let mut i1 = 0usize;
    let mut i2 = 0usize;
    let mut out1 = Vec::new();
    let mut out2 = Vec::new();
    let mut created = 0i32;
    let limit = (n1 + n2 - 1) as i32;
    if edges1.is_empty() || edges2.is_empty() {
        return Err("ShapeAlgo_AlgoContainer::HomoWires".into());
    }
    let (mut curve1, mut loc1, mut first1, mut last1) = edge_span(&edges1[i1])?;
    let (mut curve2, mut loc2, mut first2, mut last2) = edge_span(&edges2[i2])?;
    let mut reverse1 = edges1[i1].0.orientation() == crate::abs::Orientation::Reversed;
    let mut reverse2 = edges2[i2].0.orientation() == crate::abs::Orientation::Reversed;
    let mut delta1 = last1 - first1;
    let mut delta2 = last2 - first2;
    while created < limit {
        let mut step1 = false;
        let mut step2 = false;
        let (piece1, piece2) = if delta1 * ratio - delta2 > epsilon {
            let (edge_a, next_first, next_last) = if !reverse1 {
                let end = first1 + delta2 / ratio;
                (span_edge(&curve1, first1, end)?, end, last1)
            } else {
                let start = last1 - delta2 / ratio;
                (span_edge(&curve1, start, last1)?, first1, start)
            };
            first1 = next_first;
            last1 = next_last;
            delta1 = last1 - first1;
            step2 = true;
            created += 1;
            (edge_a, span_edge(&curve2, first2, last2)?)
        } else if (delta1 * ratio - delta2).abs() <= epsilon {
            step1 = true;
            step2 = true;
            created += 2;
            (
                span_edge(&curve1, first1, last1)?,
                span_edge(&curve2, first2, last2)?,
            )
        } else {
            let (edge_b, next_first, next_last) = if !reverse2 {
                let end = first2 + delta1 * ratio;
                (span_edge(&curve2, first2, end)?, end, last2)
            } else {
                let start = last2 - delta1 * ratio;
                (span_edge(&curve2, start, last2)?, first2, start)
            };
            first2 = next_first;
            last2 = next_last;
            delta2 = last2 - first2;
            step1 = true;
            created += 1;
            (span_edge(&curve1, first1, last1)?, edge_b)
        };
        let mut edge_a = piece1;
        let mut edge_b = piece2;
        edge_a.0.move_location(&loc1);
        edge_b.0.move_location(&loc2);
        if reverse1 {
            edge_a.0.reverse();
        }
        if reverse2 {
            edge_b.0.reverse();
        }
        out1.push(edge_a);
        out2.push(edge_b);
        if step1 && step2 && i1 + 1 >= n1 {
            return Ok((
                TopoBuilder::new().make_wire(&out1),
                TopoBuilder::new().make_wire(&out2),
            ));
        }
        if step1 {
            i1 += 1;
            if i1 >= n1 {
                return Err("ShapeAlgo_AlgoContainer::HomoWires".into());
            }
            let span = edge_span(&edges1[i1])?;
            curve1 = span.0;
            loc1 = span.1;
            first1 = span.2;
            last1 = span.3;
            delta1 = last1 - first1;
            reverse1 = edges1[i1].0.orientation() == crate::abs::Orientation::Reversed;
        }
        if step2 {
            i2 += 1;
            if i2 >= n2 {
                return Err("ShapeAlgo_AlgoContainer::HomoWires".into());
            }
            let span = edge_span(&edges2[i2])?;
            curve2 = span.0;
            loc2 = span.1;
            first2 = span.2;
            last2 = span.3;
            delta2 = last2 - first2;
            reverse2 = edges2[i2].0.orientation() == crate::abs::Orientation::Reversed;
        }
    }
    Err("ShapeAlgo_AlgoContainer::HomoWires".into())
}

fn copy_wire(wire: &Wire) -> Wire {
    TopoBuilder::new().make_wire(&edges_of_wire(wire))
}

fn edge_span(
    edge: &Edge,
) -> Result<(Arc<dyn Curve>, TopLocLocation, f64, f64), String> {
    let reg = GeometryRegistry::global();
    let curve = reg
        .edge_curve(&edge.0)
        .ok_or_else(|| "BRep_Tool::Curve".to_string())?;
    let (first, last) = reg.edge_parameters(&edge.0);
    Ok((curve, edge.0.location().clone(), first, last))
}

fn span_edge(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Result<Edge, String> {
    let mut edge = TopoBuilder::new().make_edge(curve.clone(), first, last);
    let v1 = TopoBuilder::new().make_vertex(curve.d0(first), CONFUSION);
    let v2 = TopoBuilder::new().make_vertex(curve.d0(last), CONFUSION);
    TopoBuilder::new().add_edge_vertices(&mut edge, &v1, &v2);
    Ok(edge)
}

/// `HasTransf` / `ConvertLocation` / `Move` (`cxx:676-692`). A failed conversion
/// keeps the shape and only warns (`IGES_1035`).
fn apply_ruled_location(
    session: &Session<'_>,
    index: usize,
    mut shape: TopoShape,
) -> Result<TopoShape, String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            if tr.form() != occt_core::gp::TrsfForm::Identity {
                shape.move_location(&TopLocLocation::composed(
                    &tr,
                    TopLocLocation::identity(),
                ));
            }
        }
    }
    Ok(shape)
}
