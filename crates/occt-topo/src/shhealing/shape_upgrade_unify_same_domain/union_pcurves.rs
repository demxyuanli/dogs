//! Pcurve concatenation for a newly merged edge.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:1754-2157`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::bspl::knots::reparametrize;
use occt_core::convert::ParameterisationType;
use occt_core::elib::clib2d;
use occt_core::gp::{GpAx2d, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom2d::comp_curve_to_bspline::CompCurveToBSplineCurve2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_convert::{curve_to_bspline_curve_bspl, DEFAULT_PARAMETERISATION};
use occt_geom2d::geom2d_convert_concat::{c0_bspline_to_c1_bspline_curve, concat_c1};
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{Geom2dCircle, Geom2dLine};

use crate::abs::Orientation;
use crate::brep_offset::BRepOffsetTool;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{common_vertex, shape_map_key, ShapeMapKey};

use super::edge_chain::AncestorMap;
use super::same_domain::FacePlaneMap;

struct CurvePart {
    curve: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
    forward: bool,
}

fn unique_non_planar_faces(
    chain: &[Edge],
    edge_faces: &AncestorMap,
    face_plane_map: &FacePlaneMap,
    face_new_face: &HashMap<ShapeMapKey, Face>,
) -> Vec<Face> {
    let Some(first) = chain.first() else {
        return Vec::new();
    };
    let mut faces = Vec::new();
    let mut seen_surfaces: Vec<(usize, occt_core::toploc::TopLocLocation)> = Vec::new();
    for shape in edge_faces
        .get(&shape_map_key(&first.0))
        .into_iter()
        .flatten()
    {
        let Some(mut face) = Face::wrap(shape.clone()) else {
            continue;
        };
        if face_plane_map.contains_key(&shape_map_key(&face.0)) {
            continue;
        }
        if let Some(replacement) = face_new_face.get(&shape_map_key(&face.0)) {
            face = replacement.clone();
        }
        face.0.set_orientation(Orientation::Forward);
        let Some(surface) = BRepTool::face_surface(&face) else {
            continue;
        };
        if surface.gp_pln().is_some() {
            continue;
        }
        let key = Arc::as_ptr(&surface) as *const () as usize;
        let location = face.0.location();
        if !seen_surfaces
            .iter()
            .any(|(seen_key, seen_location)| *seen_key == key && seen_location.is_equal(&location))
        {
            seen_surfaces.push((key, location.clone()));
            faces.push(face);
        }
    }
    faces
}

fn untrimmed_copy(curve: &Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    if let Some(basis) = curve.trimmed_basis() {
        Arc::from(basis.clone_dyn())
    } else {
        Arc::from(curve.clone_dyn())
    }
}

fn same_curve_parameters(
    previous: &Arc<dyn Curve2d>,
    curve: &Arc<dyn Curve2d>,
    first: f64,
    last: f64,
) -> Option<(f64, f64)> {
    if Arc::ptr_eq(previous, curve) {
        return Some((first, last));
    }
    if let (Some(previous_line), Some(_)) = (previous.gp_lin2d(), curve.gp_lin2d()) {
        let first_point = curve.d0(first);
        let last_point = curve.d0(last);
        if previous_line.contains(&first_point, CONFUSION)
            && previous_line.contains(&last_point, CONFUSION)
        {
            let mut new_first = clib2d::parameter_lin2d(&previous_line, &first_point);
            let mut new_last = clib2d::parameter_lin2d(&previous_line, &last_point);
            if new_first > new_last {
                std::mem::swap(&mut new_first, &mut new_last);
            }
            return Some((new_first, new_last));
        }
    }
    if let (Some(previous_circle), Some(circle)) = (previous.gp_circ2d(), curve.gp_circ2d()) {
        if circle.location().distance(&previous_circle.location()) <= CONFUSION
            && (circle.radius() - previous_circle.radius()).abs() <= CONFUSION
        {
            let mut new_first = clib2d::parameter_circ2d(&previous_circle, &curve.d0(first));
            let mut new_last = clib2d::parameter_circ2d(&previous_circle, &curve.d0(last));
            if new_first > new_last {
                std::mem::swap(&mut new_first, &mut new_last);
            }
            return Some((new_first, new_last));
        }
    }
    None
}

fn collect_parts(
    chain: &[Edge],
    face: &Face,
    reverse_for_seam: bool,
    vertex_tolerances: &mut Vec<f64>,
) -> Vec<CurvePart> {
    let mut parts: Vec<CurvePart> = Vec::new();
    let mut previous_edge: Option<Edge> = None;
    for source_edge in chain {
        let mut edge = source_edge.clone();
        let forward = edge.0.orientation() != Orientation::Reversed;
        if reverse_for_seam {
            edge.0.reverse();
        }
        let Some((source, first, last)) =
            crate::boptools_2d::curve_on_surface_oriented(&edge, face, false)
        else {
            continue;
        };
        let curve: Arc<dyn Curve2d> = if (source.is_bspline2d() || source.is_bezier2d())
            && super::curve_helpers::try_make_line(source.as_ref(), first, last).is_some()
        {
            Arc::new(
                super::curve_helpers::try_make_line(source.as_ref(), first, last)
                    .expect("line recognition was checked"),
            )
        } else {
            source.clone()
        };
        if let Some(previous) = parts.last_mut() {
            if let Some((new_first, new_last)) =
                same_curve_parameters(&previous.curve, &curve, first, last)
            {
                if previous.forward {
                    previous.last = new_last;
                } else {
                    previous.first = new_first;
                }
                previous_edge = Some(edge);
                continue;
            }
            if let Some(common) = previous_edge
                .as_ref()
                .and_then(|previous_edge| common_vertex(previous_edge, &edge))
                .and_then(crate::shape::Vertex::wrap)
            {
                vertex_tolerances.push(BRepTool::vertex_tolerance(&common));
            }
        }
        parts.push(CurvePart {
            curve: untrimmed_copy(&curve),
            first,
            last,
            forward,
        });
        previous_edge = Some(edge);
    }
    parts
}

fn concatenate_parts(
    parts: &[CurvePart],
    vertex_tolerances: &[f64],
    max_tolerance: &mut f64,
) -> Result<Option<(Arc<dyn Curve2d>, f64, f64)>, String> {
    if parts.is_empty() {
        return Ok(None);
    }
    if parts.len() == 1 {
        let part = &parts[0];
        let mut curve = part.curve.clone_dyn();
        let (mut first, mut last) = (part.first, part.last);
        if !part.forward {
            let new_last = curve.reversed_parameter(first);
            let new_first = curve.reversed_parameter(last);
            curve.reverse();
            first = new_first;
            last = new_last;
        }
        return Ok(Some((Arc::from(curve), first, last)));
    }

    let mut splines = Vec::with_capacity(parts.len());
    for part in parts {
        let mut trimmed = Geom2dTrimmedCurve::new(part.curve.clone(), part.first, part.last);
        if !part.forward {
            trimmed.reverse();
        }
        let mut spline = curve_to_bspline_curve_bspl(&trimmed, DEFAULT_PARAMETERISATION)
            .ok_or("UnionPCurves: CurveToBSplineCurve failed")?;
        c0_bspline_to_c1_bspline_curve(&mut spline, CONFUSION).map_err(str::to_string)?;
        splines.push(spline);
    }
    let tolerances = &vertex_tolerances[..parts.len().saturating_sub(1)];
    for tolerance in tolerances {
        *max_tolerance = max_tolerance.max(*tolerance);
    }
    let mut closed = false;
    let (_, mut concatenated) = concat_c1(
        &mut splines,
        tolerances,
        &mut closed,
        CONFUSION,
        occt_core::precision::ANGULAR,
    )
    .map_err(str::to_string)?;
    let first = concatenated
        .drain(..1)
        .next()
        .ok_or("UnionPCurves: ConcatC1 returned no curve")?;
    let mut compound =
        CompCurveToBSplineCurve2d::from_bspline(first, ParameterisationType::TgtThetaOver2);
    for spline in &concatenated {
        if !compound
            .add(spline, *max_tolerance, true)
            .map_err(str::to_string)?
        {
            return Err("UnionPCurves: concatenated pcurves are disconnected".into());
        }
    }
    let result = compound
        .into_curve()
        .ok_or("UnionPCurves: missing concatenated pcurve")?;
    let first = result.first_parameter();
    let last = result.last_parameter();
    Ok(Some((Arc::new(result), first, last)))
}

fn reparameterized(
    curve: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
    requested_first: f64,
    requested_last: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    if let Some(line) = curve.gp_lin2d() {
        let direction = line.direction();
        let point = curve.d0(first);
        let origin = GpPnt2d::new(
            point.x() - direction.x() * requested_first,
            point.y() - direction.y() * requested_first,
        );
        return Ok(Arc::new(Geom2dLine::new(GpAx2d::new(origin, *direction))));
    }
    if let Some(circle) = curve.gp_circ2d() {
        let radius = circle.radius();
        let location = circle.location();
        let mut position = *circle.position();
        position.rotate(&location, first - requested_first);
        return Ok(Arc::new(Geom2dCircle::new(occt_core::gp::GpCirc2d::new(
            position, radius,
        ))));
    }
    let trimmed = Geom2dTrimmedCurve::new(curve, first, last);
    let mut spline = curve_to_bspline_curve_bspl(&trimmed, DEFAULT_PARAMETERISATION)
        .ok_or("UnionPCurves: CurveToBSplineCurve failed during reparameterization")?;
    let (mut knots, _) = spline.distinct_knots_and_mults();
    reparametrize(requested_first, requested_last, &mut knots);
    spline.set_knots(&knots).map_err(str::to_string)?;
    Ok(Arc::new(spline))
}

fn project_or_reparameterize(
    edge: &Edge,
    face: &Face,
    curve: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
    requested_first: f64,
    requested_last: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    let curve3d = BRepTool::edge_curve(edge).ok_or("UnionPCurves: merged edge has no 3D curve")?;
    let surface = BRepTool::face_surface_world(face).ok_or("UnionPCurves: face has no surface")?;
    if let Some(projected) = BRepOffsetTool::geom_projlib_curve2d(
        &curve3d,
        requested_first,
        requested_last,
        &surface,
        CONFUSION,
    ) {
        Ok(projected)
    } else {
        reparameterized(curve, first, last, requested_first, requested_last)
    }
}

/// `UnionPCurves` (`cxx:1754-2157`).
pub(super) fn union_pcurves(
    chain: &[Edge],
    edge: &Edge,
    edge_faces: &AncestorMap,
    face_plane_map: &FacePlaneMap,
    face_new_face: &HashMap<ShapeMapKey, Face>,
) -> Result<(), String> {
    let faces = unique_non_planar_faces(chain, edge_faces, face_plane_map, face_new_face);
    let (first3d, last3d) = BRepTool::edge_parameters(edge);
    let edge_tolerance = BRepTool::edge_tolerance(edge);
    let mut max_tolerance = edge_tolerance;
    let mut results: Vec<(Face, Arc<dyn Curve2d>, f64, f64)> = Vec::new();
    let mut seam = false;
    let mut face_index = 0;
    while face_index < faces.len() + usize::from(seam) {
        let face = faces[face_index.min(faces.len().saturating_sub(1))].clone();
        let reverse_for_seam = seam && face_index == 1;
        let mut vertex_tolerances = Vec::new();
        let parts = collect_parts(chain, &face, reverse_for_seam, &mut vertex_tolerances);
        if !seam && faces.len() == 1 {
            seam = chain.iter().any(|source| {
                let mut reversed = source.clone();
                reversed.0.reverse();
                match (
                    crate::boptools_2d::curve_on_surface_oriented(source, &face, false),
                    crate::boptools_2d::curve_on_surface_oriented(&reversed, &face, false),
                ) {
                    (Some((first, _, _)), Some((second, _, _))) => !Arc::ptr_eq(&first, &second),
                    _ => false,
                }
            });
        }
        if let Some((curve, first, last)) =
            concatenate_parts(&parts, &vertex_tolerances, &mut max_tolerance)?
        {
            results.push((face, curve, first, last));
        }
        face_index += 1;
    }

    let bad_range = results
        .iter()
        .any(|(_, _, first, last)| ((last3d - first3d) - (last - first)).abs() > max_tolerance);
    for (face, curve, first, last) in &mut results {
        if bad_range {
            *curve = project_or_reparameterize(
                edge,
                face,
                curve.clone(),
                *first,
                *last,
                first3d,
                last3d,
            )?;
            *first = first3d;
            *last = last3d;
        } else if (first3d - *first).abs() > max_tolerance || (last3d - *last).abs() > max_tolerance
        {
            *curve = reparameterized(curve.clone(), *first, *last, first3d, last3d)?;
            *first = first3d;
            *last = last3d;
        }
    }

    let registry = GeometryRegistry::global();
    if seam && results.len() >= 2 {
        let face_key = GeometryRegistry::shape_key(&results[0].0 .0);
        registry.set_edge_pcurves(
            &edge.0,
            face_key,
            vec![results[0].1.clone(), results[1].1.clone()],
        );
        registry.set_pcurve_range(&edge.0, face_key, first3d, last3d);
    } else {
        for (face, curve, _, _) in results {
            let face_key = GeometryRegistry::shape_key(&face.0);
            registry.set_edge_pcurve(&edge.0, face_key, curve);
            registry.set_pcurve_range(&edge.0, face_key, first3d, last3d);
        }
    }
    Ok(())
}
