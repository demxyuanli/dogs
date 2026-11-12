//! Periodic and closed-surface preparation for `IntUnifyFaces`.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:3540-3828`.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::{GeomConvertApproxSurface, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_tools::{is_really_closed, outer_wire};
use crate::builder::TopoBuilder;
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edges_of_wire, shape_map_key, ShapeMapKey};

use super::uv_bounds::EdgeFaceMap;

pub(super) struct PeriodicFaceState {
    pub surface: Arc<dyn Surface>,
    pub reference_face: Face,
    pub u_period: f64,
    pub v_period: f64,
}

fn copy_pcurves(
    edges: &[Edge],
    old_face: &Face,
    new_face: &Face,
    temporary_pcurves: &HashSet<ShapeMapKey>,
) {
    let registry = GeometryRegistry::global();
    let old_key = GeometryRegistry::shape_key(&old_face.0);
    let new_key = GeometryRegistry::shape_key(&new_face.0);
    for edge in edges {
        let curves = registry.edge_pcurves(&edge.0, old_key);
        let range = registry.pcurve_range(&edge.0, old_key);
        if temporary_pcurves.contains(&shape_map_key(&edge.0)) {
            registry.remove_pcurves_on_surface(&edge.0, &old_face.0);
        }
        if !curves.is_empty() {
            registry.set_edge_pcurves(&edge.0, new_key, curves);
        }
        if let Some((first, last)) = range {
            registry.set_pcurve_range(&edge.0, new_key, first, last);
        }
    }
}

fn replacement_face(
    surface: Arc<dyn Surface>,
    old_face: &Face,
    edges: &[Edge],
    temporary_pcurves: &HashSet<ShapeMapKey>,
) -> Face {
    let builder = TopoBuilder::new();
    let mut face = builder.make_face_with_tolerance(surface, &[], 0.0);
    face.0.set_orientation(Orientation::Forward);
    copy_pcurves(edges, old_face, &face, temporary_pcurves);
    face
}

fn seam_state(faces: &[Face]) -> (bool, bool, Option<Edge>) {
    let mut u_seam = false;
    let mut v_seam = false;
    let mut edge_with_two_pcurves = None;
    for face in faces {
        let Some(wire) = outer_wire(face) else {
            continue;
        };
        for edge in edges_of_wire(&wire) {
            if !BRepTool::is_closed_edge_face(&edge, face) {
                continue;
            }
            if is_really_closed(&edge, face) {
                if super::wire_helpers::is_u_iso(&edge, face).unwrap_or(false) {
                    u_seam = true;
                } else {
                    v_seam = true;
                }
            } else {
                edge_with_two_pcurves = Some(edge);
            }
        }
    }
    (u_seam, v_seam, edge_with_two_pcurves)
}

fn periodic_bspline(
    surface: &Arc<dyn Surface>,
    make_u: bool,
    make_v: bool,
) -> Result<Arc<dyn Surface>, String> {
    let mut bspline = if let Some(bspline) = surface.osculating_bspline() {
        bspline
    } else {
        GeomConvertApproxSurface::new(
            surface.as_ref(),
            1.0e-4,
            Shape::C1,
            Shape::C1,
            14,
            14,
            16,
            1,
        )
        .into_surface()
        .ok_or("IntUnifyFaces: GeomConvert_ApproxSurface produced no surface")?
    };
    if make_u {
        bspline.set_u_periodic();
    }
    if make_v {
        bspline.set_v_periodic();
    }
    Ok(Arc::new(bspline))
}

fn smooth_two_pcurve_edge(edge: &Edge, edge_faces: &EdgeFaceMap) -> bool {
    let Some(adjacent) = edge_faces.get(&shape_map_key(&edge.0)) else {
        return false;
    };
    let (Some(first), Some(last)) = (adjacent.first(), adjacent.last()) else {
        return false;
    };
    let (Some(first), Some(last)) = (Face::wrap(first.clone()), Face::wrap(last.clone())) else {
        return false;
    };
    is_tangent_faces(edge, &first, &last, Shape::G1).unwrap_or(false)
}

fn closed_direction(edge: &Edge, face: &Face) -> Option<bool> {
    let registry = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&face.0);
    let curves = registry.edge_pcurves(&edge.0, key);
    if curves.len() < 2 {
        return None;
    }
    let first = registry
        .pcurve_range(&edge.0, key)
        .unwrap_or_else(|| BRepTool::edge_parameters(edge))
        .0;
    let p1 = curves[0].d0(first);
    let p2 = curves[1].d0(first);
    Some((p1.x() - p2.x()).abs() > (p1.y() - p2.y()).abs())
}

fn complete_relocated_pcurves(
    edges: &[Edge],
    reference_face: &Face,
    relocated: &mut HashMap<ShapeMapKey, Arc<dyn Curve2d>>,
    used: &HashSet<ShapeMapKey>,
) -> Result<(), String> {
    let registry = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&reference_face.0);
    for edge in edges {
        let key = shape_map_key(&edge.0);
        if used.contains(&key) {
            continue;
        }
        let curve = registry
            .edge_pcurve(&edge.0, face_key)
            .ok_or("IntUnifyFaces: boundary edge has no pcurve")?;
        let (first, last) = registry
            .pcurve_range(&edge.0, face_key)
            .unwrap_or_else(|| BRepTool::edge_parameters(edge));
        relocated.insert(key, Arc::new(Geom2dTrimmedCurve::new(curve, first, last)));
    }
    Ok(())
}

fn install_relocated_pcurves(
    edges: &[Edge],
    old_face: &Face,
    new_face: &Face,
    relocated: &HashMap<ShapeMapKey, Arc<dyn Curve2d>>,
    temporary_pcurves: &HashSet<ShapeMapKey>,
) {
    let registry = GeometryRegistry::global();
    let old_key = GeometryRegistry::shape_key(&old_face.0);
    let new_key = GeometryRegistry::shape_key(&new_face.0);
    for edge in edges {
        let key = shape_map_key(&edge.0);
        let range = registry.pcurve_range(&edge.0, old_key);
        if temporary_pcurves.contains(&key) {
            registry.remove_pcurves_on_surface(&edge.0, &old_face.0);
        }
        if let Some(curve) = relocated.get(&key) {
            registry.set_edge_pcurve(&edge.0, new_key, curve.clone());
        }
        if let Some((first, last)) = range {
            registry.set_pcurve_range(&edge.0, new_key, first, last);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_periodic_faces(
    faces: &[Face],
    edges: &[Edge],
    mut surface: Arc<dyn Surface>,
    mut reference_face: Face,
    edge_faces: &EdgeFaceMap,
    boundary_edges: &HashSet<ShapeMapKey>,
    temporary_pcurves: &HashSet<ShapeMapKey>,
    concat_bsplines: bool,
    _angle_tolerance: f64,
) -> Result<PeriodicFaceState, String> {
    let (mut u_seam, mut v_seam, edge_with_two_pcurves) = seam_state(faces);
    let mut u_period = if surface.is_u_periodic() {
        surface.u_period()
    } else {
        0.0
    };
    let mut v_period = if surface.is_v_periodic() {
        surface.v_period()
    } else {
        0.0
    };

    if concat_bsplines && !u_seam && !v_seam {
        if let Some(edge) =
            edge_with_two_pcurves.filter(|edge| smooth_two_pcurve_edge(edge, edge_faces))
        {
            let u_closed = closed_direction(&edge, &reference_face).unwrap_or(false);
            let make_u = u_closed && u_period == 0.0;
            let make_v = !u_closed && v_period == 0.0;
            if make_u || make_v {
                surface = periodic_bspline(&surface, make_u, make_v)?;
                reference_face =
                    replacement_face(surface.clone(), &reference_face, edges, temporary_pcurves);
                u_period = if surface.is_u_periodic() {
                    surface.u_period()
                } else {
                    0.0
                };
                v_period = if surface.is_v_periodic() {
                    surface.v_period()
                } else {
                    0.0
                };
            }
        }
    }

    let coordinate_tolerance = super::uv_bounds::minimum_edge_size(edges, &reference_face)
        .map(|value| (value.0 / 10.0).max(CONFUSION))
        .unwrap_or(CONFUSION);
    let periods = [u_period, v_period];
    let mut seams = [u_seam, v_seam];
    for direction in 0..2 {
        let period = periods[direction];
        if period == 0.0 || seams[direction] {
            continue;
        }
        let coordinate_index = direction + 1;
        let Some((min_coordinate, max_coordinate, interval_count, max_face_index)) =
            super::uv_bounds::find_coord_bounds(
                faces,
                &reference_face,
                edge_faces,
                boundary_edges,
                coordinate_index,
                period,
            )
        else {
            break;
        };
        if max_coordinate - min_coordinate > period - 1.0e-5 {
            seams[direction] = true;
            continue;
        }
        if interval_count != 2 {
            continue;
        }
        let mut vertex_edges = super::wire_helpers::vertex_edge_map(edges);
        let Some((mut relocated, used)) = super::uv_bounds::relocate_pcurves_to_new_origin(
            edges,
            &faces[max_face_index],
            &reference_face,
            coordinate_tolerance,
            coordinate_index,
            period,
            &mut vertex_edges,
        ) else {
            continue;
        };
        complete_relocated_pcurves(edges, &reference_face, &mut relocated, &used)?;
        let mut new_min = f64::INFINITY;
        let mut new_max = f64::NEG_INFINITY;
        for edge in edges {
            let curve = &relocated[&shape_map_key(&edge.0)];
            super::uv_bounds::update_boundaries(
                curve.as_ref(),
                curve.first_parameter(),
                curve.last_parameter(),
                coordinate_index,
                &mut new_min,
                &mut new_max,
            );
        }
        if new_max - new_min >= period - coordinate_tolerance
            || (new_min > -CONFUSION
                && new_min < period + CONFUSION
                && new_max > -CONFUSION
                && new_max < period + CONFUSION)
        {
            continue;
        }
        let range = if direction == 0 {
            surface.u_range()
        } else {
            surface.v_range()
        };
        let remaining = period - (new_max - new_min);
        let origin = (new_min - remaining / 2.0).max(range.0);
        let new_surface: Arc<dyn Surface> = if origin == range.0 {
            surface.clone()
        } else {
            Arc::new(GeomRectangularTrimmedSurface::one_param(
                surface.clone(),
                origin,
                origin + period,
                direction == 0,
            )?)
        };
        let mut new_face =
            TopoBuilder::new().make_face_with_tolerance(new_surface.clone(), &[], 0.0);
        new_face.0.set_orientation(Orientation::Forward);
        install_relocated_pcurves(
            edges,
            &reference_face,
            &new_face,
            &relocated,
            temporary_pcurves,
        );
        surface = new_surface;
        reference_face = new_face;
    }
    u_seam = seams[0];
    v_seam = seams[1];
    let _ = (u_seam, v_seam);

    let closed_basis = surface
        .rectangular_trimmed_basis()
        .unwrap_or_else(|| surface.clone());
    if u_period == 0.0 && closed_basis.is_u_closed() {
        let (first, last) = closed_basis.u_range();
        u_period = last - first;
    }
    if v_period == 0.0 && closed_basis.is_v_closed() {
        let (first, last) = closed_basis.v_range();
        v_period = last - first;
    }
    Ok(PeriodicFaceState {
        surface,
        reference_face,
        u_period,
        v_period,
    })
}
