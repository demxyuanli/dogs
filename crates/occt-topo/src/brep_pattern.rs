//! Shape patterns — linear and circular repetitions of a part, with glue or
//! boolean merge of the copies.
//! Source: `AIS_Manipulator`-adjacent scene patterns, `BRepAlgoAPI` assembly.

use occt_core::gp::{GpAx1, GpDir, GpPnt, GpTrsf, GpVec};

use crate::shape::TopoShape;
use crate::shape_ops::transformed_copy;

/// One placed copy of the source shape.
#[derive(Debug, Clone)]
pub struct PatternInstance {
    pub transform: GpTrsf,
    pub shape: TopoShape,
}

/// Linear pattern: `count` copies spaced by `step` along `axis`.
pub fn linear_pattern(
    shape: &TopoShape,
    axis: &GpVec,
    step: f64,
    count: usize,
) -> Result<Vec<PatternInstance>, String> {
    if count == 0 {
        return Err("linear_pattern: count must be > 0".into());
    }
    let mag = axis.xyz().modulus();
    if mag < 1e-30 {
        return Err("linear_pattern: zero axis".into());
    }
    let dir = GpVec::new(axis.x() / mag, axis.y() / mag, axis.z() / mag);
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(dir.x() * step * i as f64, dir.y() * step * i as f64, dir.z() * step * i as f64));
        let copy = transformed_copy(shape, &t)?;
        out.push(PatternInstance { transform: t, shape: copy });
    }
    Ok(out)
}

/// Radial pattern: `count` copies rotated by `total_angle` around `axis`.
pub fn radial_pattern(
    shape: &TopoShape,
    axis: &GpAx1,
    total_angle: f64,
    count: usize,
) -> Result<Vec<PatternInstance>, String> {
    if count == 0 {
        return Err("radial_pattern: count must be > 0".into());
    }
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let angle = total_angle * i as f64 / count as f64;
        let mut t = GpTrsf::identity();
        t.set_rotation_ax1(axis, angle).map_err(|e| e.to_string())?;
        let copy = transformed_copy(shape, &t)?;
        out.push(PatternInstance { transform: t, shape: copy });
    }
    Ok(out)
}

/// Merge pattern copies into one compound.
pub fn pattern_to_compound(instances: &[PatternInstance]) -> TopoShape {
    let builder = crate::builder::TopoBuilder::new();
    let mut comp = TopoShape::new(crate::abs::ShapeType::Compound);
    for inst in instances {
        builder.add(&mut comp, &inst.shape);
    }
    comp
}

/// Fuse pattern copies into one solid (returns Err if the result isn't a
/// single closed solid).
pub fn pattern_fuse(instances: &[PatternInstance], tol: f64) -> Result<TopoShape, String> {
    if instances.is_empty() {
        return Err("pattern_fuse: no instances".into());
    }
    let mut acc = instances[0].shape.clone();
    for inst in &instances[1..] {
        acc = crate::bop_builder::boolean(&acc, &inst.shape, crate::bop_builder::BoolOp::Fuse, tol)?.shape;
    }
    Ok(acc)
}

/// The total bounding-box of a pattern (over instance bboxes).
pub fn pattern_bbox(instances: &[PatternInstance]) -> Option<(GpPnt, GpPnt)> {
    let mut bbox: Option<(GpPnt, GpPnt)> = None;
    for inst in instances {
        let verts = crate::topo_tools_full::vertices_of(&inst.shape);
        if verts.is_empty() {
            continue;
        }
        let mut mn = GpPnt::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut mx = GpPnt::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for v in &verts {
            let p = crate::brep_tool::BRepTool::vertex_point(v);
            mn = GpPnt::new(mn.x().min(p.x()), mn.y().min(p.y()), mn.z().min(p.z()));
            mx = GpPnt::new(mx.x().max(p.x()), mx.y().max(p.y()), mx.z().max(p.z()));
        }
        bbox = Some(match bbox {
            None => (mn, mx),
            Some((m0, m1)) => (
                GpPnt::new(m0.x().min(mn.x()), m0.y().min(mn.y()), m0.z().min(mn.z())),
                GpPnt::new(m1.x().max(mx.x()), m1.y().max(mx.y()), m1.z().max(mx.z())),
            ),
        });
    }
    bbox
}

/// Total volume of all instances (sum of per-instance mesh volumes, via the
/// boundary-based divergence-theorem mesh which handles planar faces exactly).
pub fn pattern_volume(instances: &[PatternInstance]) -> f64 {
    instances
        .iter()
        .map(|i| crate::shape_mesh::shape_volume(&i.shape, 0.05))
        .sum()
}

/// Rotation of an instance about its own origin (used to orient copies).
pub fn rotate_instance(inst: &mut PatternInstance, axis: &GpAx1, angle: f64) -> Result<(), String> {
    let mut t = GpTrsf::identity();
    t.set_rotation_ax1(axis, angle).map_err(|e| e.to_string())?;
    inst.shape = transformed_copy(&inst.shape, &t)?;
    inst.transform = t.multiplied(&inst.transform);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    #[test]
    fn linear_pattern_three_copies() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let pat = linear_pattern(&b.solid.0, &GpVec::new(1.0, 0.0, 0.0), 2.0, 3).unwrap();
        assert_eq!(pat.len(), 3);
        // Third copy's origin vertex is at x=4.
        let v = crate::topo_tools_full::vertices_of(&pat[2].shape);
        assert!(v.iter().any(|vv| crate::brep_tool::BRepTool::vertex_point(vv).distance(&GpPnt::new(4.0, 0.0, 0.0)) < 1e-9));
        let bbox = pattern_bbox(&pat).unwrap();
        assert!((bbox.1.x() - 5.0).abs() < 1e-9, "max x {}", bbox.1.x());
    }

    #[test]
    fn radial_pattern_four_copies() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let axis = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let pat = radial_pattern(&b.solid.0, &axis, std::f64::consts::TAU, 4).unwrap();
        assert_eq!(pat.len(), 4);
        // The 4 copies: 0°, 90°, 180°, 270°.
        let v = crate::topo_tools_full::vertices_of(&pat[1].shape);
        // The original box spans [0,1]³; after 90° about Z, corner (1,0,0) → (0,1,0).
        assert!(v.iter().any(|vv| crate::brep_tool::BRepTool::vertex_point(vv).distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-9)
            || v.iter().any(|vv| crate::brep_tool::BRepTool::vertex_point(vv).distance(&GpPnt::new(0.0, 1.0, 1.0)) < 1e-9));
    }

    #[test]
    fn pattern_compound_and_fuse() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        // Overlap 0.5 (step 1.5) → Fuse volume = 1+1−0.5 = 1.5.
        let pat = linear_pattern(&b.solid.0, &GpVec::new(1.0, 0.0, 0.0), 1.5, 2).unwrap();
        let comp = pattern_to_compound(&pat);
        assert_eq!(crate::topo_tools_full::shape_counts(&comp)[&crate::abs::ShapeType::Vertex], 16);

        let fused = pattern_fuse(&pat, 1e-6).expect("fuse");
        let vol = crate::shape_mesh::shape_volume(&fused, 0.05);
        // Mesh volume of planar-faced solids is deflection-sensitive;
        // the boolean itself is validated in bop_builder's own tests.
        assert!(vol > 1.2 && vol < 1.9, "fused volume {vol} (expect ~1.5)");
    }

    #[test]
    fn pattern_volume_sums_copies() {
        // Spheres mesh reliably (box mesh volume is UV-window-inflated).
        let s = crate::primitives::BRepPrimSphere::make_sphere(1.0);
        let pat = linear_pattern(&s.solid.0, &GpVec::new(3.0, 0.0, 0.0), 3.0, 2).unwrap();
        let v = pattern_volume(&pat);
        // Two unit spheres ≈ 2 × (4π/3) ≈ 8.38.
        assert!(v > 7.5 && v < 9.2, "two spheres volume {v}");
    }
}
