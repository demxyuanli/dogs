//! Phase 5 module: feature — feature modeling (`BRepFeat_MakePrism` simplified).
//!
//! Ports the common boss/hole/protrusion/pocket feature operations. Each
//! operation builds a tool solid (a swept prism or a cylinder), then combines
//! it with the target solid with a boolean (Fuse for protrusions, Cut for
//! pockets). The exact BRep boolean (`BRepAlgoAPI_Fuse`/`Cut`) is not yet
//! available in this port, so the operation falls back to the voxel boolean
//! (`boolean_ops::voxel_boolean`) and records a warning. `FeatureResult`
//! carries the result solid, its `TopoShape`, the collected warnings and the
//! result volume (computed from the boolean mesh, so it stays accurate even
//! though re-tessellating the converted solid can overshoot planar faces).
//!
//! ponytail: the cylindrical tool for `boss`/`hole` is built as an exact
//! triangle mesh converted to BRep, rather than `BRepPrim_Cylinder`, because
//! the analytic `GeomCylinder` lateral surface cannot be meshed correctly by
//! `face_to_triangles` (`d1` reports zero partials, so `face_uv_bounds` falls
//! back to a [0,1]² window). A mesh-built tool has only planar faces, which
//! mesh exactly. Replace when the analytic cylinder meshing is fixed.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::poly::triangulation::Triangle;

use crate::abs::ShapeType;
use crate::boolean_ops::{voxel_boolean, BoolOp};
use crate::mesh::ShapeMesh;
use crate::mesh_to_brep::shape_mesh_to_brep;
use crate::shape::{Solid, TopoShape};

/// Result of a feature operation.
#[derive(Debug, Clone)]
pub struct FeatureResult {
    pub solid: Solid,
    pub shape: TopoShape,
    pub warnings: Vec<String>,
    /// Volume of the result, computed from the boolean mesh (divergence
    /// theorem) before it was converted to BRep.
    pub volume: f64,
}

/// Extrude a planar sketch polygon by `dir × height` into a prism and fuse it
/// onto `solid`. The prism starts at the sketch's plane and extends in `dir`.
pub fn protrusion(solid: &Solid, sketch: &[GpPnt], dir: &GpVec, height: f64, tol: f64) -> Result<FeatureResult, String> {
    if sketch.len() < 3 {
        return Err("protrusion: sketch needs at least 3 points".into());
    }
    if height <= 0.0 {
        return Err("protrusion: height must be positive".into());
    }
    let prism = crate::sweep::prism_from_polygon(sketch, dir, height);
    boolean_feature(solid, &prism.solid.0, BoolOp::Union, tol)
}

/// Extrude a planar sketch polygon by `dir × depth` into a prism and cut it
/// out of `solid` (`solid` minus the prism).
pub fn pocket(solid: &Solid, sketch: &[GpPnt], dir: &GpVec, depth: f64, tol: f64) -> Result<FeatureResult, String> {
    if sketch.len() < 3 {
        return Err("pocket: sketch needs at least 3 points".into());
    }
    if depth <= 0.0 {
        return Err("pocket: depth must be positive".into());
    }
    let prism = crate::sweep::prism_from_polygon(sketch, dir, depth);
    boolean_feature(solid, &prism.solid.0, BoolOp::Subtraction, tol)
}

/// Add a cylindrical boss fused onto `solid`: a cylinder of `radius` and
/// `height` whose base circle is centered at `center` (axis +Z), fused with a
/// Fuse boolean.
pub fn boss(solid: &Solid, center: &GpPnt, radius: f64, height: f64, tol: f64) -> Result<FeatureResult, String> {
    if radius <= 0.0 || height <= 0.0 {
        return Err("boss: radius and height must be positive".into());
    }
    let tool = cylinder_tool(center, radius, height)?;
    boolean_feature(solid, &tool, BoolOp::Union, tol)
}

/// Cut a cylindrical hole out of `solid`: a cylinder of `radius` and `depth`
/// whose base circle is centered at `center` (axis +Z), removed with a Cut
/// boolean.
pub fn hole(solid: &Solid, center: &GpPnt, radius: f64, depth: f64, tol: f64) -> Result<FeatureResult, String> {
    if radius <= 0.0 || depth <= 0.0 {
        return Err("hole: radius and depth must be positive".into());
    }
    let tool = cylinder_tool(center, radius, depth)?;
    boolean_feature(solid, &tool, BoolOp::Subtraction, tol)
}

/// `|vol(after) − vol(before)|`, the material added or removed by the feature.
pub fn feature_volume_before_after(before: &Solid, after: &FeatureResult) -> f64 {
    let v0 = crate::shape_mesh::shape_volume(&before.0, 0.1);
    (after.volume - v0).abs()
}

// ---------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------

/// Run a voxel boolean between `solid` and `tool`, converting the resulting
/// mesh to a BRep solid and recording the voxel fallback as a warning.
///
/// **UNPORTED (audit A5/A12)**: OCCT has no mesh/voxel boolean — a feature runs
/// `BRepFeat_MakePrism`/`BRepFeat_MakeRevol` (or `LocOpe_*`) on top of
/// `BOPAlgo_BOP` with a real tool solid. The port's feature subsystem is
/// mesh-based here; replacing it belongs to A12/T-48, which also owns the
/// analytic `brepfeat` path.
fn boolean_feature(solid: &Solid, tool: &TopoShape, op: BoolOp, tol: f64) -> Result<FeatureResult, String> {
    let resolution = resolution_for(solid, tool, tol);
    let mesh = voxel_boolean(&solid.0, tool, resolution, op)
        .map_err(|e| format!("feature boolean ({op:?}): {e}"))?;

    let mut vol = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t.n0].coord;
        let b = mesh.vertices[t.n1].coord;
        let c = mesh.vertices[t.n2].coord;
        vol += a.dot_cross(&b, &c);
    }
    let volume = (vol / 6.0).abs();

    let brep = shape_mesh_to_brep(&mesh);
    let solid_out = match brep.solid {
        Some(s) => s,
        // A voxel staircase can be non-manifold in places; still expose it as
        // a solid so the feature result is usable for further processing.
        None => {
            let b = crate::builder::TopoBuilder::new();
            b.make_solid(&[brep.shell])
        }
    };
    let warnings = vec![format!(
        "voxel boolean fallback used (resolution {resolution}); exact BRep boolean unavailable"
    )];
    Ok(FeatureResult {
        solid: solid_out.clone(),
        shape: solid_out.0,
        warnings,
        volume,
    })
}

/// Build a closed BRep cylinder solid (axis +Z) whose base circle is centered
/// at `center`, from an exact triangle mesh.
fn cylinder_tool(center: &GpPnt, radius: f64, height: f64) -> Result<TopoShape, String> {
    const SLICES: usize = 24;
    let mesh = z_cylinder_mesh(radius, height, SLICES);
    let vertices: Vec<GpPnt> = mesh
        .vertices
        .iter()
        .map(|p| GpPnt::new(p.x() + center.x(), p.y() + center.y(), p.z() + center.z()))
        .collect();
    let mesh = ShapeMesh { vertices, triangles: mesh.triangles, source_shape: ShapeType::Solid };
    let brep = shape_mesh_to_brep(&mesh);
    brep.solid
        .map(|s| s.0)
        .ok_or("feature: cylinder tool mesh is not a closed solid".into())
}

/// Exact cylinder mesh along +Z: `radius`, height `height`, `slices` azimuthal
/// segments. Winding mirrors `crate::mesh::mesh_cylinder` (orientation is
/// irrelevant to the even-odd classifier used by the voxel boolean).
fn z_cylinder_mesh(radius: f64, height: f64, slices: usize) -> ShapeMesh {
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..=slices {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
        let x = radius * theta.cos();
        let y = radius * theta.sin();
        vertices.push(GpPnt::new(x, y, 0.0));
        vertices.push(GpPnt::new(x, y, height));
    }
    for i in 0..slices {
        let a = 2 * i;
        let b = 2 * i + 1;
        let c = 2 * i + 2;
        let d = 2 * i + 3;
        triangles.push(Triangle::new(a, b, c));
        triangles.push(Triangle::new(b, d, c));
    }
    let top_idx = vertices.len();
    let bot_idx = vertices.len() + 1;
    vertices.push(GpPnt::new(0.0, 0.0, height));
    vertices.push(GpPnt::new(0.0, 0.0, 0.0));
    for i in 0..slices {
        let tb = 2 * i + 1;
        let tt = 2 * i + 3;
        let bb = 2 * i;
        let bt = 2 * i + 2;
        triangles.push(Triangle::new(top_idx, tt, tb));
        triangles.push(Triangle::new(bot_idx, bb, bt));
    }
    ShapeMesh { vertices, triangles, source_shape: ShapeType::Solid }
}

/// Choose the voxel resolution from the combined bounding-box span so a cell
/// is about `tol` wide (clamped for speed and robustness).
fn resolution_for(solid: &Solid, tool: &TopoShape, tol: f64) -> usize {
    let mut bb = crate::bbox_from_geometry::shape_bbox(&solid.0);
    bb.add_box(&crate::bbox_from_geometry::shape_bbox(tool));
    let (x0, x1, y0, y1, z0, z1) = bb.get().unwrap_or((0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
    let span = (x1 - x0).max(y1 - y0).max(z1 - z0).max(1e-9);
    let tol = if tol.is_finite() && tol > 0.0 { tol } else { span / 40.0 };
    let r = (span / tol.max(span / 64.0)).ceil() as usize;
    r.clamp(16, 64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    #[test]
    fn boss_on_box_increases_volume() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (r, h) = (0.25, 1.0);
        let after = boss(&b.solid, &GpPnt::new(0.5, 0.5, 1.0), r, h, 0.05).expect("boss ok");
        let expect = std::f64::consts::PI * r * r * h;
        let diff = feature_volume_before_after(&b.solid, &after);
        assert!(
            (diff - expect).abs() < 0.15 * expect,
            "boss added {diff}, expected ≈{expect}"
        );
        assert!(
            after.warnings.iter().any(|w| w.contains("voxel boolean fallback")),
            "expected a fallback warning, got {:?}",
            after.warnings
        );
    }

    #[test]
    fn hole_on_box_decreases_volume() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (r, h) = (0.25, 1.0);
        let after = hole(&b.solid, &GpPnt::new(0.5, 0.5, 0.0), r, h, 0.05).expect("hole ok");
        let expect = std::f64::consts::PI * r * r * h;
        let diff = feature_volume_before_after(&b.solid, &after);
        assert!(
            (diff - expect).abs() < 0.15 * expect,
            "hole removed {diff}, expected ≈{expect}"
        );
    }

    #[test]
    fn protrusion_square_sketch_fuses_closed() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let sketch = [
            GpPnt::new(0.25, 0.25, 1.0),
            GpPnt::new(0.75, 0.25, 1.0),
            GpPnt::new(0.75, 0.75, 1.0),
            GpPnt::new(0.25, 0.75, 1.0),
        ];
        let after = protrusion(&b.solid, &sketch, &GpVec::new(0.0, 0.0, 1.0), 1.0, 0.05)
            .expect("protrusion ok");
        // 0.5×0.5×1 block fused on top of the unit box.
        let diff = feature_volume_before_after(&b.solid, &after);
        assert!(
            (diff - 0.25).abs() < 0.1,
            "protrusion added {diff}, expected ≈0.25"
        );
    }

    #[test]
    fn pocket_removes_material() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let sketch = [
            GpPnt::new(0.25, 0.25, 1.0),
            GpPnt::new(0.75, 0.25, 1.0),
            GpPnt::new(0.75, 0.75, 1.0),
            GpPnt::new(0.25, 0.75, 1.0),
        ];
        let after = pocket(&b.solid, &sketch, &GpVec::new(0.0, 0.0, -1.0), 0.5, 0.05)
            .expect("pocket ok");
        // 0.5×0.5×0.5 block removed from the top.
        let diff = feature_volume_before_after(&b.solid, &after);
        assert!(
            (diff - 0.125).abs() < 0.1,
            "pocket removed {diff}, expected ≈0.125"
        );
    }
}
