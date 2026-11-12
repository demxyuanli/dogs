//! Same-domain predicates from `ShapeUpgrade_UnifySameDomain.cxx:1356-1645`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpCylinder, GpPln, GpVec};
use occt_core::precision::ANGULAR;
use occt_geom::geom_lib_is_planar_surface::IsPlanarSurface;
use occt_geom::{GeomPlane, Surface};

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::geom_int::TopolTool;
use crate::intpatch::ImpImpIntersection;
use crate::shape::{Edge, Face};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

pub(super) type FacePlaneMap = HashMap<ShapeMapKey, Arc<dyn Surface>>;

/// `getCylinder` (`ShapeUpgrade_UnifySameDomain.cxx:1356-1434`).
pub(super) fn get_cylinder(surface: &Arc<dyn Surface>) -> Option<GpCylinder> {
    if let Some(cylinder) = surface.gp_cylinder() {
        return Some(cylinder);
    }

    if surface.is_surface_of_revolution() {
        let axis = surface.revolution_axis()?;
        let basis = untrimmed_basis(surface.revolution_basis_curve()?);
        let line = basis.gp_line()?;
        if !line.direction().is_parallel_tol(axis.direction(), ANGULAR) {
            return None;
        }
        let location = *axis.location();
        let position = GpAx3::from_point_dir(location, *axis.direction()).ok()?;
        return GpCylinder::new(position, line.distance(&location)).ok();
    }

    if surface.is_surface_of_linear_extrusion() {
        let direction = surface.extrusion_direction()?;
        let basis = untrimmed_basis(surface.extrusion_basis_curve()?);
        let circle = basis.gp_circ()?;
        if !circle
            .position()
            .direction()
            .is_parallel_tol(&direction, ANGULAR)
        {
            return None;
        }
        let position = GpAx3::from_point_dir(circle.location(), direction).ok()?;
        return GpCylinder::new(position, circle.radius()).ok();
    }

    None
}

fn untrimmed_basis(mut curve: Arc<dyn occt_geom::Curve>) -> Arc<dyn occt_geom::Curve> {
    while curve.is_geom_trimmed() {
        let Some(basis) = curve.trim_basis() else {
            break;
        };
        curve = basis;
    }
    curve
}

/// `GetNormalToSurface` (`ShapeUpgrade_UnifySameDomain.cxx:1447-1517`).
pub(super) fn normal_to_surface_on_edge(
    face: &Face,
    edge: &Edge,
    parameter: f64,
) -> Option<occt_core::gp::GpDir> {
    let mut forward_face = face.clone();
    forward_face.0.set_orientation(Orientation::Forward);

    let pcurve = if BRepTool::is_closed_edge_face(edge, face) {
        let edge_in_face = explorer_shapes(&forward_face.0, crate::abs::ShapeType::Edge)
            .into_iter()
            .find(|candidate| is_same(candidate, &edge.0))
            .and_then(Edge::wrap)?;
        curve_on_surface_oriented(&edge_in_face, &forward_face, false)?.0
    } else {
        curve_on_surface_oriented(edge, face, false)?.0
    };

    let uv = pcurve.d0(parameter);
    let surface = BRepTool::face_surface(face)?;
    let (_, du, dv) = surface.d1(uv.x(), uv.y());
    let mut normal = du.crossed(&dv);
    if normal.magnitude() < occt_core::precision::CONFUSION {
        return None;
    }
    if face.0.orientation() == Orientation::Reversed {
        normal.reverse();
    }
    normal.transform(&face.0.location().transformation());
    occt_core::gp::GpDir::from_vec(&normal).ok()
}

/// `IsSameDomain` (`ShapeUpgrade_UnifySameDomain.cxx:1518-1645`).
pub(super) fn is_same_domain(
    face: &Face,
    checked_face: &Face,
    linear_tolerance: f64,
    angular_tolerance: f64,
    face_plane_map: &mut FacePlaneMap,
) -> bool {
    let (Some(local1), Some(local2)) = (
        BRepTool::face_surface(face),
        BRepTool::face_surface(checked_face),
    ) else {
        return false;
    };
    if Arc::ptr_eq(&local1, &local2) && face.0.location().is_equal(checked_face.0.location()) {
        return true;
    }

    let (Some(world1), Some(world2)) = (
        BRepTool::face_surface_world(face),
        BRepTool::face_surface_world(checked_face),
    ) else {
        return false;
    };
    let surface1 = super::same_surface::clear_rectangular_trim(&world1);
    let surface2 = super::same_surface::clear_rectangular_trim(&world2);

    let planar1 = IsPlanarSurface::new(&surface1, linear_tolerance);
    if planar1.is_planar() {
        let planar2 = IsPlanarSurface::new(&surface2, linear_tolerance);
        if planar2.is_planar() {
            let (Ok(plane1), Ok(plane2)) = (planar1.plan(), planar2.plan()) else {
                return false;
            };
            if planes_coincide(&plane1, &plane2, linear_tolerance, angular_tolerance) {
                bind_common_plane(face, checked_face, plane1, face_plane_map);
                return true;
            }
        }
    }

    if is_elementary(surface1.as_ref()) && is_elementary(surface2.as_ref()) {
        let domain1 = TopolTool::from_surface(surface1.as_ref());
        let domain2 = TopolTool::from_surface(surface2.as_ref());
        let mut intersection = ImpImpIntersection::new();
        intersection.perform(
            surface1.as_ref(),
            &domain1,
            surface2.as_ref(),
            &domain2,
            linear_tolerance,
            linear_tolerance,
        );
        return intersection.is_done() && !intersection.is_empty() && intersection.tangent_faces();
    }

    if is_cylindrical_or_swept(surface1.as_ref()) && is_cylindrical_or_swept(surface2.as_ref()) {
        let (Some(cylinder1), Some(cylinder2)) = (get_cylinder(&surface1), get_cylinder(&surface2))
        else {
            return false;
        };
        if (cylinder1.radius() - cylinder2.radius()).abs() >= linear_tolerance {
            return false;
        }
        let direction1 = cylinder1.position().direction();
        let direction2 = cylinder2.position().direction();
        if !direction1.is_parallel_tol(&direction2, ANGULAR) {
            return false;
        }
        let locations = GpVec::from_pnts(&cylinder1.location(), &cylinder2.location());
        return locations.square_magnitude() < linear_tolerance * linear_tolerance
            || locations.is_parallel_ang(&GpVec::from_xyz(direction1.xyz()), ANGULAR);
    }

    false
}

fn planes_coincide(
    first: &GpPln,
    second: &GpPln,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> bool {
    first
        .position()
        .direction()
        .is_parallel_tol(&second.position().direction(), angular_tolerance)
        && first.distance_pnt(&second.location()) < linear_tolerance
}

fn bind_common_plane(
    face: &Face,
    checked_face: &Face,
    plane: GpPln,
    face_plane_map: &mut FacePlaneMap,
) {
    let first_key = shape_map_key(&face.0);
    let second_key = shape_map_key(&checked_face.0);
    let common = face_plane_map
        .get(&first_key)
        .or_else(|| face_plane_map.get(&second_key))
        .cloned()
        .unwrap_or_else(|| Arc::new(GeomPlane::new(plane)));
    face_plane_map.insert(first_key, common.clone());
    face_plane_map.insert(second_key, common);
}

fn is_elementary(surface: &dyn Surface) -> bool {
    matches!(
        classify_surface(surface),
        SurfaceKind::Plane
            | SurfaceKind::Cylinder
            | SurfaceKind::Cone
            | SurfaceKind::Sphere
            | SurfaceKind::Torus
    )
}

fn is_cylindrical_or_swept(surface: &dyn Surface) -> bool {
    surface.gp_cylinder().is_some()
        || surface.is_surface_of_revolution()
        || surface.is_surface_of_linear_extrusion()
}
