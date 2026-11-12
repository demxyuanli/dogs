//! Pcurve transfer from `ShapeUpgrade_UnifySameDomain.cxx:1057-1300`.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpAx3, GpDir2d, GpPnt2d, GpTrsf2d, GpVec, GpVec2d};
use occt_core::precision::{APPROXIMATION, PCONFUSION};
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_offset::BRepOffsetTool;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

/// `TransformPCurves` (`ShapeUpgrade_UnifySameDomain.cxx:1057-1269`).
pub(super) fn transform_pcurves(
    reference_face: &Face,
    face: &Face,
    temporary_pcurves: &mut HashSet<ShapeMapKey>,
) -> Result<(), String> {
    let reference_surface = BRepTool::face_surface_world(reference_face)
        .ok_or("TransformPCurves: reference face has no surface")?;
    let face_surface =
        BRepTool::face_surface_world(face).ok_or("TransformPCurves: face has no surface")?;
    let reference_surface = super::same_surface::clear_rectangular_trim(&reference_surface);
    let face_surface = super::same_surface::clear_rectangular_trim(&face_surface);

    let (transformation, to_project) = match (
        elementary_position(reference_surface.as_ref()),
        elementary_position(face_surface.as_ref()),
    ) {
        (Some(reference_axis), Some(face_axis)) => {
            (TransformSpec::between(&reference_axis, &face_axis), false)
        }
        _ if super::same_surface::same_surface(&reference_surface, &face_surface) => {
            (TransformSpec::default(), false)
        }
        _ => (TransformSpec::default(), true),
    };

    let registry = GeometryRegistry::global();
    let reference_key = GeometryRegistry::shape_key(&reference_face.0);
    let face_key = GeometryRegistry::shape_key(&face.0);
    let mut seen = HashSet::new();

    for edge_shape in explorer_shapes(&face.0, ShapeType::Edge) {
        let edge = Edge(edge_shape);
        if BRepTool::is_degenerated(&edge) && transformation.modifies() {
            continue;
        }
        if !seen.insert(shape_map_key(&edge.0)) {
            continue;
        }

        let mut forward_reference = reference_face.clone();
        forward_reference.0.set_orientation(Orientation::Forward);
        let mut forward_edge = edge.clone();
        forward_edge.0.set_orientation(Orientation::Forward);
        let mut reversed_edge = forward_edge.clone();
        reversed_edge.0.reverse();

        let existing_first =
            curve_on_surface_oriented(&forward_edge, &forward_reference, false).map(|v| v.0);
        let existing_second =
            curve_on_surface_oriented(&reversed_edge, &forward_reference, false).map(|v| v.0);
        if let (Some(first), Some(second)) = (&existing_first, &existing_second) {
            if !Arc::ptr_eq(first, second) {
                continue;
            }
        }

        let source_curves = registry.edge_pcurves(&edge.0, face_key);
        if source_curves.is_empty() {
            return Err("TransformPCurves: source edge has no stored pcurve; \
                 OCCT lines 1182-1186 require one"
                .into());
        }
        let curve_count = source_curves.len().min(2);
        let mut transformed = Vec::with_capacity(curve_count);
        for source in source_curves.iter().take(curve_count) {
            if to_project {
                let curve3d = registry
                    .edge_curve(&edge.0)
                    .ok_or("TransformPCurves: source edge has no 3D curve")?;
                let (first, last) = registry.edge_parameters(&edge.0);
                let tolerance = BRepTool::edge_tolerance(&edge).min(APPROXIMATION);
                let projected = BRepOffsetTool::geom_projlib_curve2d(
                    &curve3d,
                    first,
                    last,
                    &reference_surface,
                    tolerance,
                )
                .ok_or("TransformPCurves: GeomProjLib::Curve2d failed")?;
                transformed.push(projected);
            } else {
                let mut copied = source.clone_dyn();
                transformation.apply(copied.as_mut())?;
                transformed.push(Arc::<dyn Curve2d>::from(copied));
            }
        }

        if transformed.len() == 1 {
            let closed = reference_surface.is_u_closed() || reference_surface.is_v_closed();
            if existing_second.is_none() || !closed {
                registry.set_edge_pcurve(&edge.0, reference_key, transformed[0].clone());
                temporary_pcurves.insert(shape_map_key(&edge.0));
            } else if let Some(existing) = existing_first {
                let (first, _) = registry.edge_parameters(&edge.0);
                let old_point = existing.d0(first);
                let new_point = transformed[0].d0(first);
                let (umin, umax) = reference_surface.u_range();
                let (vmin, vmax) = reference_surface.v_range();
                let u_period = if reference_surface.is_u_closed() {
                    umax - umin
                } else {
                    0.0
                };
                let v_period = if reference_surface.is_v_closed() {
                    vmax - vmin
                } else {
                    0.0
                };
                let separated = (u_period != 0.0
                    && (old_point.x() - new_point.x()).abs() > u_period / 2.0)
                    || (v_period != 0.0 && (old_point.y() - new_point.y()).abs() > v_period / 2.0);
                if separated {
                    let pair = if edge.0.orientation() == Orientation::Forward {
                        vec![transformed[0].clone(), existing]
                    } else {
                        vec![existing, transformed[0].clone()]
                    };
                    registry.set_edge_pcurves(&edge.0, reference_key, pair);
                    temporary_pcurves.insert(shape_map_key(&edge.0));
                }
            }
        } else {
            registry.set_edge_pcurves(&edge.0, reference_key, transformed);
            temporary_pcurves.insert(shape_map_key(&edge.0));
        }

        let range = registry
            .pcurve_range(&edge.0, face_key)
            .unwrap_or_else(|| registry.edge_parameters(&edge.0));
        registry.set_pcurve_range(&edge.0, reference_key, range.0, range.1);
    }
    Ok(())
}

/// `AddPCurves` (`ShapeUpgrade_UnifySameDomain.cxx:1270-1300`).
pub(super) fn add_pcurves(
    faces: &[Face],
    reference_face: &Face,
    temporary_pcurves: &mut HashSet<ShapeMapKey>,
) -> Result<(), String> {
    let reference_surface = BRepTool::face_surface_world(reference_face)
        .ok_or("AddPCurves: reference face has no surface")?;
    if reference_surface.gp_pln().is_some() {
        return Ok(());
    }
    for face in faces {
        let mut forward = face.clone();
        forward.0.set_orientation(Orientation::Forward);
        if is_same(&forward.0, &reference_face.0) {
            continue;
        }
        transform_pcurves(reference_face, &forward, temporary_pcurves)?;
    }
    Ok(())
}

#[derive(Clone, Copy, Default)]
struct TransformSpec {
    translation: f64,
    angle: f64,
    reverse_x: bool,
    reverse_y: bool,
}

impl TransformSpec {
    fn between(reference: &GpAx3, face: &GpAx3) -> Self {
        let origin = reference.location();
        let face_origin = face.location();
        let face_direction = face.direction();
        let delta = GpVec::from_pnts(&face_origin, &origin);
        let parameter = delta.x() * face_direction.x()
            + delta.y() * face_direction.y()
            + delta.z() * face_direction.z();
        let translation = if parameter.abs() > PCONFUSION {
            -parameter
        } else {
            0.0
        };

        let reference_cross = reference
            .x_direction()
            .crossed(reference.y_direction())
            .ok();
        let face_cross = face.x_direction().crossed(face.y_direction()).ok();
        let reverse_x = matches!(
            (reference_cross, face_cross),
            (Some(a), Some(b)) if a.dot(&b) < 0.0
        );
        let reverse_y = face.direction().dot(&reference.direction()) < 0.0;

        let angle = if !reverse_x && !reverse_y {
            let mut reference_direction = reference.direction();
            if !reference.is_direct() {
                reference_direction.reverse();
            }
            reference
                .x_direction()
                .angle_with_ref(face.x_direction(), &reference_direction)
        } else {
            reference.x_direction().angle(face.x_direction())
        };

        Self {
            translation,
            angle: if angle.abs() > PCONFUSION { angle } else { 0.0 },
            reverse_x,
            reverse_y,
        }
    }

    fn modifies(self) -> bool {
        self.translation != 0.0 || self.angle != 0.0 || self.reverse_x || self.reverse_y
    }

    fn apply(self, curve: &mut dyn Curve2d) -> Result<(), String> {
        if self.translation != 0.0 {
            transform_curve(curve, translation(0.0, self.translation));
        }
        if self.reverse_y {
            transform_curve(curve, mirror_axis(1.0, 0.0)?);
        }
        if self.reverse_x {
            transform_curve(curve, mirror_axis(0.0, 1.0)?);
            transform_curve(curve, translation(2.0 * std::f64::consts::PI, 0.0));
        }
        if self.angle != 0.0 {
            transform_curve(curve, translation(self.angle, 0.0));
        }
        Ok(())
    }
}

fn elementary_position(surface: &dyn Surface) -> Option<GpAx3> {
    surface
        .gp_pln()
        .map(|v| v.position())
        .or_else(|| surface.gp_cylinder().map(|v| v.position()))
        .or_else(|| surface.gp_cone().map(|v| v.position()))
        .or_else(|| surface.gp_sphere().map(|v| *v.position()))
        .or_else(|| surface.gp_torus().map(|v| *v.position()))
}

fn translation(x: f64, y: f64) -> GpTrsf2d {
    let mut transform = GpTrsf2d::identity();
    transform.set_translation_vec(&GpVec2d::new(x, y));
    transform
}

fn mirror_axis(x: f64, y: f64) -> Result<GpTrsf2d, String> {
    let direction = GpDir2d::new(x, y).map_err(str::to_string)?;
    let axis = GpAx2d::new(GpPnt2d::new(0.0, 0.0), direction);
    let mut transform = GpTrsf2d::identity();
    transform.set_mirror_ax2d(&axis);
    Ok(transform)
}

fn transform_curve(curve: &mut dyn Curve2d, transform: GpTrsf2d) {
    curve.transform(&transform);
}
