
use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::shape::TopoShape;

use super::*;

/// Even-odd point-in-solid test.
///
/// Casts a ray in +X and counts crossings against a coarse surface mesh of the
/// faces (grid-triangulated). Points strictly outside the bounding box are
/// rejected first. If the shape cannot be meshed (no faces with surfaces), the
/// bounding-box containment is returned as a fallback. Best suited to closed
/// solids; results on open surfaces are not meaningful.
/// Whether `p` lies inside the solid `shape` (audit A2 / task T-38).
///
/// Solids go through `BRepClass3d_SolidClassifier` (the port's
/// [`SolidClassifier::classify`], `brep_class3d/mod.rs:463`, itself a port of
/// `BRepClass3d_SClassifier::Perform`); `In` is the answer, matching the
/// strictly-interior semantics the previous body had.
///
/// For a shape that is **not** a solid (`BRepClass3d_SolidClassifier` is only
/// defined for solids/shells with a closed volume) this keeps the port-internal
/// approximation: a 7x7 UV mesh per face plus a jittered +X ray parity count.
/// UNPORTED there: OCCT classifies a shell through
/// `BRepClass3d_SolidExplorer` on the shell's own faces
/// (`BRepClass3d_SolidExplorer.cxx:80-120`), which this port has not wired for
/// non-solid shapes yet.
pub fn is_inside(shape: &TopoShape, p: &GpPnt) -> bool {
    match shape.shape_type() {
        ShapeType::Solid | ShapeType::CompSolid => matches!(
            crate::brep_class3d::SolidClassifier::classify(
                shape,
                p,
                occt_core::precision::CONFUSION
            ),
            crate::fclass2d::FaceState::In
        ),
        _ => is_inside_mesh(shape, p),
    }
}

/// Port-internal mesh + ray-parity point classification, kept only for
/// non-solid shapes (see [`is_inside`]): a 7x7 UV mesh of every face, then the
/// parity of a slightly jittered +X ray. UNPORTED for that case.
pub(super) fn is_inside_mesh(shape: &TopoShape, p: &GpPnt) -> bool {
    let bbox = shape_bbox(shape);
    if let Some(b) = &bbox {
        if b.is_out(p) {
            return false;
        }
    }

    let meshes = mesh_faces(shape, 7, 7);
    if meshes.is_empty() {
        // Unmeshed shape: fall back to bounding-box containment.
        return match bbox {
            Some(b) => !b.is_out(p),
            None => false,
        };
    }

    // Jitter the ray origin by a tiny amount perpendicular to the +X ray so it
    // does not land exactly on a grid cell diagonal (a degenerate ray that can
    // be double-counted or skipped by both triangles of a cell).
    let jittered = GpPnt::new(p.x(), p.y() + 1e-7, p.z() + 1e-7);
    let mut hits: Vec<f64> = Vec::new();
    for (verts, tris) in &meshes {
        for (i, j, k) in tris {
            if let Some(t) = ray_plus_x_hits(&verts[*i], &verts[*j], &verts[*k], &jittered) {
                hits.push(t);
            }
        }
    }
    // A ray passing exactly through a shared triangle edge reports two hits at
    // the same parameter; deduplicate coincident hits so the parity stays
    // correct for closed meshes.
    hits.sort_by(f64::total_cmp);
    let mut unique = 0usize;
    let mut prev: Option<f64> = None;
    for t in hits {
        if prev.map_or(true, |q| (t - q).abs() > 1e-9) {
            unique += 1;
            prev = Some(t);
        }
    }
    unique % 2 == 1
}

/// Minimum separation between two axis-aligned bounding boxes; 0 if they
/// overlap.
pub fn shape_bbox_distance(a: &BndBox, b: &BndBox) -> f64 {
    let (a0, a1, ay0, ay1, az0, az1) = match a.get() {
        Some(v) => v,
        None => return f64::INFINITY,
    };
    let (b0, b1, by0, by1, bz0, bz1) = match b.get() {
        Some(v) => v,
        None => return f64::INFINITY,
    };
    let gx = axis_gap(a0, a1, b0, b1);
    let gy = axis_gap(ay0, ay1, by0, by1);
    let gz = axis_gap(az0, az1, bz0, bz1);
    (gx * gx + gy * gy + gz * gz).sqrt()
}

pub(super) fn axis_gap(amin: f64, amax: f64, bmin: f64, bmax: f64) -> f64 {
    if amax < bmin {
        return bmin - amax;
    }
    if bmax < amin {
        return amin - bmax;
    }
    0.0
}

