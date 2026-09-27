//! Shape-level geometric transforms — apply a `GpTrsf` to the *registered*
//! geometry (vertex points, edge curves, face surfaces) of an entire shape
//! subtree.
//!
//! This differs from `crate::transform::move_shape`, which only relabels
//! `TopoShape::location` and leaves the side-table geometry in its local
//! frame. `transform_shape` rewrites the geometry itself, so the world
//! coordinates change (what `BRep_Tool::Transform` does in OCCT).
//!
//! Source: `BRep_Tool::Transform`, `TopoDS_Shape::Move` (geometry-in-place).

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpPnt, GpTrsf, GpVec};
use occt_core::toploc::TopLocLocation;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::tgeometry::{EdgeGeom, FaceGeom, GeometryRegistry, VertexGeom};

/// Apply `t` to every vertex point, edge curve and face surface in the shape's
/// subtree, then reset the shape's location to identity (the transform is now
/// baked into the geometry). Mirrors `BRep_Tool::Transform`.
pub fn transform_shape(shape: &mut TopoShape, t: &GpTrsf) -> Result<(), String> {
    let reg = GeometryRegistry::global();
    let mut stack = vec![shape.clone()];
    let mut seen = std::collections::HashSet::new();
    while let Some(s) = stack.pop() {
        if !seen.insert(std::sync::Arc::as_ptr(&s.tshape)) {
            continue;
        }
        match s.shape_type() {
            ShapeType::Vertex => {
                if let Some(g) = reg.vertex_geom(&s) {
                    let p = g.point.transformed(t);
                    reg.set_vertex(&s, VertexGeom { point: p, tolerance: g.tolerance });
                }
            }
            ShapeType::Edge => {
                if let Some(g) = reg.edge_geom(&s) {
                    let curve = Arc::from(g.curve.transformed(t));
                    reg.set_edge(&s, EdgeGeom {
                        curve,
                        first: g.first,
                        last: g.last,
                        tolerance: g.tolerance,
                        same_parameter: g.same_parameter,
                        same_range: g.same_range,
                        degenerated: g.degenerated,
                        pcurves: g.pcurves.clone(),
                        pcurve_ranges: g.pcurve_ranges.clone(),
                    });
                }
            }
            ShapeType::Face => {
                if let Some(g) = reg.face_geom(&s) {
                    let surface = Arc::from(g.surface.transformed(t));
                    reg.set_face(&s, FaceGeom {
                        surface,
                        tolerance: g.tolerance,
                        natural_restriction: g.natural_restriction,
                    });
                }
            }
            _ => {}
        }
        for k in s.tshape.read().unwrap().children.clone() {
            stack.push(k);
        }
    }
    shape.location = TopLocLocation::identity();
    Ok(())
}

/// Translate a shape's geometry by `v` (in place).
pub fn translate_shape(shape: &mut TopoShape, v: &GpVec) -> Result<(), String> {
    let mut t = GpTrsf::identity();
    t.set_translation_vec(v);
    transform_shape(shape, &t)
}

/// Rotate a shape's geometry about an axis (in place).
pub fn rotate_shape(shape: &mut TopoShape, axis: &GpAx1, angle: f64) -> Result<(), String> {
    let mut t = GpTrsf::identity();
    t.set_rotation_ax1(axis, angle).map_err(|e| e.to_string())?;
    transform_shape(shape, &t)
}

/// Scale a shape's geometry about a center (in place).
pub fn scale_shape(shape: &mut TopoShape, center: &GpPnt, s: f64) -> Result<(), String> {
    let mut t = GpTrsf::identity();
    t.set_scale(center, s).map_err(|e| e.to_string())?;
    transform_shape(shape, &t)
}

/// Deep-copy a shape subtree with fresh `TShape` handles and transformed
/// geometry. The source geometry is read, transformed, and re-registered on
/// new shapes, so the two shapes are fully independent. Source:
/// `TopoDS::Transformed`, `BRep_Tool::Copy`.
pub fn transformed_copy(shape: &TopoShape, t: &GpTrsf) -> Result<TopoShape, String> {
    let reg = GeometryRegistry::global();
    // Map source TShape ptr → copied TopoShape.
    let mut map: std::collections::HashMap<usize, TopoShape> = std::collections::HashMap::new();
    // Source/copy pairs, for the pcurve re-keying pass below.
    let mut pairs: Vec<(TopoShape, TopoShape)> = Vec::new();
    let mut stack = vec![shape.clone()];
    while let Some(s) = stack.pop() {
        let key = std::sync::Arc::as_ptr(&s.tshape) as usize;
        if map.contains_key(&key) {
            continue;
        }
        let copy = TopoShape::new(s.shape_type());
        match s.shape_type() {
            ShapeType::Vertex => {
                if let Some(g) = reg.vertex_geom(&s) {
                    let p = g.point.transformed(t);
                    reg.set_vertex(&copy, VertexGeom { point: p, tolerance: g.tolerance });
                }
            }
            ShapeType::Edge => {
                if let Some(g) = reg.edge_geom(&s) {
                    reg.set_edge(&copy, EdgeGeom {
                        curve: Arc::from(g.curve.transformed(t)),
                        first: g.first,
                        last: g.last,
                        tolerance: g.tolerance,
                        same_parameter: g.same_parameter,
                        same_range: g.same_range,
                        degenerated: g.degenerated,
                        pcurves: g.pcurves.clone(),
                        pcurve_ranges: g.pcurve_ranges.clone(),
                    });
                }
            }
            ShapeType::Face => {
                if let Some(g) = reg.face_geom(&s) {
                    reg.set_face(&copy, FaceGeom {
                        surface: Arc::from(g.surface.transformed(t)),
                        tolerance: g.tolerance,
                        natural_restriction: g.natural_restriction,
                    });
                }
            }
            _ => {}
        }
        map.insert(key, copy.clone());
        pairs.push((s.clone(), copy.clone()));
        for k in s.tshape.read().unwrap().children.clone() {
            stack.push(k);
        }
    }
    // Rebuild children: attach each copied child to its copied parent.
    let mut stack = vec![(shape.clone(), map[&(std::sync::Arc::as_ptr(&shape.tshape) as usize)].clone())];
    let mut seen = std::collections::HashSet::new();
    while let Some((src, dst)) = stack.pop() {
        if !seen.insert(std::sync::Arc::as_ptr(&src.tshape) as usize) {
            continue;
        }
        for k in src.tshape.read().unwrap().children.clone() {
            let ck = std::sync::Arc::as_ptr(&k.tshape) as usize;
            if let Some(copy) = map.get(&ck) {
                // A copy is structural: it keeps the child's stored orientation
                // and location. `BRepTools_ShapeSet::Add` +
                // `BRepTools_ShapeSet::Read` and `BRepBuilderAPI_Copy` never
                // touch `TopoDS_Shape::myOrientation`, so a vertex stored
                // `REVERSED` on its edge stays `REVERSED` in the copy; dropping
                // it would make `TopExp::FirstVertex`/`LastVertex` (which read
                // that flag, `TopExp.cxx:182-210`) return nothing.
                let mut child_shape = copy.clone();
                child_shape.set_orientation(k.orientation());
                child_shape.set_location(k.location());
                dst.tshape.write().unwrap().add_child(child_shape);
                stack.push((k, copy.clone()));
            }
        }
    }
    // `EdgeGeom::pcurves` / `pcurve_ranges` are keyed by the surface data
    // pointer (`GeometryRegistry::repr_key`), so re-key them onto the copied
    // faces' surfaces. `BRepBuilderAPI_Copy` / `BRepTools_ShapeSet` copy a
    // `BRep_GCurve` together with the face it lies on, so a copied edge must
    // still answer `BRep_Tool::CurveOnSurface(E, copiedFace)` — otherwise the
    // analytic passes (`BRepGProp`) see pcurve-less faces and report a zero
    // volume for every pattern/transform copy.
    let mut surface_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for (src, dst) in &pairs {
        if src.shape_type() != ShapeType::Face {
            continue;
        }
        if let (Some(a), Some(b)) = (reg.face_surface(src), reg.face_surface(dst)) {
            surface_map.insert(
                Arc::as_ptr(&a) as *const () as usize,
                Arc::as_ptr(&b) as *const () as usize,
            );
        }
    }
    for (src, dst) in &pairs {
        if src.shape_type() != ShapeType::Edge {
            continue;
        }
        let Some(g) = reg.edge_geom(src) else {
            continue;
        };
        let pcurves: std::collections::HashMap<usize, Vec<Arc<dyn occt_geom2d::Curve2d>>> = g
            .pcurves
            .iter()
            .map(|(k, v)| (surface_map.get(k).copied().unwrap_or(*k), v.clone()))
            .collect();
        let ranges: std::collections::HashMap<usize, (f64, f64)> = g
            .pcurve_ranges
            .iter()
            .map(|(k, v)| (surface_map.get(k).copied().unwrap_or(*k), *v))
            .collect();
        if pcurves.is_empty() {
            continue;
        }
        if let Some(mut updated) = reg.edge_geom(dst) {
            updated.pcurves = pcurves;
            updated.pcurve_ranges = ranges;
            reg.set_edge(dst, updated);
        }
    }
    Ok(map[&(std::sync::Arc::as_ptr(&shape.tshape) as usize)].clone())
}

/// Translate-copy helper.
pub fn translated_copy(shape: &TopoShape, v: &GpVec) -> Result<TopoShape, String> {
    let mut t = GpTrsf::identity();
    t.set_translation_vec(v);
    transformed_copy(shape, &t)
}

/// Whether a shape subtree carries any registered geometry at all.
pub fn has_geometry(shape: &TopoShape) -> bool {
    let reg = GeometryRegistry::global();
    let mut stack = vec![shape.clone()];
    let mut seen = std::collections::HashSet::new();
    while let Some(s) = stack.pop() {
        if !seen.insert(std::sync::Arc::as_ptr(&s.tshape)) {
            continue;
        }
        match s.shape_type() {
            ShapeType::Vertex => {
                if reg.vertex_geom(&s).is_some() {
                    return true;
                }
            }
            ShapeType::Edge => {
                if reg.edge_geom(&s).is_some() {
                    return true;
                }
            }
            ShapeType::Face => {
                if reg.face_geom(&s).is_some() {
                    return true;
                }
            }
            _ => {}
        }
        for k in s.tshape.read().unwrap().children.clone() {
            stack.push(k);
        }
    }
    false
}

/// Recursively free every registry entry for the subtree (`BRep_Tool::Free`).
pub fn clear_shape_geometry(shape: &TopoShape) {
    let reg = GeometryRegistry::global();
    let mut stack = vec![shape.clone()];
    let mut seen = std::collections::HashSet::new();
    while let Some(s) = stack.pop() {
        if !seen.insert(std::sync::Arc::as_ptr(&s.tshape)) {
            continue;
        }
        reg.clear_shape(&s);
        for k in s.tshape.read().unwrap().children.clone() {
            stack.push(k);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::{GpAx1, GpDir};

    fn approx(a: f64, b: f64) -> bool { (a - b).abs() < 1e-9 * b.abs().max(1.0) }

    #[test]
    fn translate_moves_geometry() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut solid = b.solid.0;
        translate_shape(&mut solid, &GpVec::new(5.0, 0.0, 0.0)).unwrap();
        // No vertex remains at the origin (all moved +5 X).
        let verts = crate::topo_tools_full::vertices_of(&solid);
        assert!(!verts.iter().any(|v| {
            let p = BRepTool::vertex_point(v);
            approx(p.x(), 0.0) && approx(p.y(), 0.0) && approx(p.z(), 0.0)
        }), "origin vertex moved away");
        // Some vertex is now at (5,0,0).
        assert!(verts.iter().any(|v| {
            let p = BRepTool::vertex_point(v);
            approx(p.x(), 5.0) && approx(p.y(), 0.0) && approx(p.z(), 0.0)
        }), "origin vertex now at (5,0,0)");
        // A face surface's sample point moved with the geometry.
        let faces = crate::topo_tools_full::faces_of(&solid);
        assert!(faces.iter().any(|f| {
            let s = BRepTool::face_surface(f).unwrap();
            approx(s.d0(0.0, 0.0).x(), 5.0) || approx(s.d0(0.0, 0.0).x(), 6.0)
        }), "a face surface now passes through x≈5..6");
    }

    #[test]
    fn translated_copy_is_independent() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let copy = translated_copy(&b.solid.0, &GpVec::new(2.0, 0.0, 0.0)).unwrap();
        // Original unchanged, copy moved.
        let orig_verts = crate::topo_tools_full::vertices_of(&b.solid.0);
        let copy_verts = crate::topo_tools_full::vertices_of(&copy);
        assert_eq!(orig_verts.len(), copy_verts.len());
        let has_2x = copy_verts.iter().any(|v| approx(BRepTool::vertex_point(v).x(), 2.0));
        assert!(has_2x, "copy has a vertex at x=2");
        assert!(!crate::topo_tools_full::vertices_of(&b.solid.0)
            .iter().any(|v| approx(BRepTool::vertex_point(v).x(), 2.0)));
    }

    #[test]
    fn rotate_about_z_spins_corners() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut solid = b.solid.0;
        let z = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        rotate_shape(&mut solid, &z, std::f64::consts::FRAC_PI_2).unwrap();
        // Corner (1,0,0) → (0,1,0) after +90° about Z.
        let verts = crate::topo_tools_full::vertices_of(&solid);
        let rotated = verts.iter().any(|v| {
            let p = BRepTool::vertex_point(v);
            approx(p.x(), 0.0) && approx(p.y(), 1.0) && approx(p.z(), 0.0)
        });
        assert!(rotated, "corner (1,0,0) rotated to (0,1,0)");
    }

    #[test]
    fn has_geometry_and_clear() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(has_geometry(&b.solid.0));
        clear_shape_geometry(&b.solid.0);
        assert!(!has_geometry(&b.solid.0));
    }
}
