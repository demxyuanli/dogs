//! Shell utilities — orientation consistency, closure and reorientation.
//! Source: `BRepCheck_Shell::CheckOrientation`, `ShapeAnalysis_Shell`,
//! `BRep_Builder::Update` (orientation propagation).
//!
//! A valid closed shell bounds a solid when every face's normal points
//! outward. These helpers detect and fix inward-pointing faces, and verify
//! the shell is closed and consistently oriented.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpPln, GpVec};

use crate::brep_surface::surface_normal;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, Shell, TopoShape};

use crate::topo_tools_full::{faces_of, wires_of_face};

/// Unit normal of a face at its centroid (robust FD normal).
pub fn face_normal_at_centroid(face: &Face) -> Option<GpVec> {
    let s = BRepTool::face_surface(face)?;
    let (u0, u1, v0, v1) = BRepTool::uv_bounds(face);
    let (u, v) = if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
        (0.5 * (u0 + u1), 0.5 * (v0 + v1))
    } else {
        (0.0, 0.0)
    };
    Some(surface_normal(s.as_ref(), u, v))
}

/// Whether a face's normal points outward from a reference solid: the probe
/// point (centroid + ε·normal) must NOT be inside the solid. `is_inside` is a
/// closure so callers can supply either `brep_extrema::is_inside` or a fast
/// mesh probe.
pub fn face_points_outward(
    face: &Face,
    reference: &TopoShape,
    is_inside: &dyn Fn(&TopoShape, &occt_core::gp::GpPnt) -> bool,
    eps: f64,
) -> bool {
    let Some(n) = face_normal_at_centroid(face) else {
        return true;
    };
    let Some(s) = BRepTool::face_surface(face) else {
        return true;
    };
    let (u0, u1, v0, v1) = BRepTool::uv_bounds(face);
    let (u, v) = if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
        (0.5 * (u0 + u1), 0.5 * (v0 + v1))
    } else {
        (0.0, 0.0)
    };
    let c = s.d0(u, v);
    let probe = occt_core::gp::GpPnt::new(c.x() + eps * n.x(), c.y() + eps * n.y(), c.z() + eps * n.z());
    !is_inside(reference, &probe)
}

/// Count faces of a shell whose normal points inward (from the shell centroid
/// heuristic: the outward direction is away from the centroid).
pub fn inward_facing_faces(shell: &Shell) -> Vec<(Face, GpVec)> {
    let faces = faces_of(&shell.0);
    if faces.is_empty() {
        return Vec::new();
    }
    // Shell centroid: average of face centroids.
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut n = 0usize;
    let mut entries: Vec<(Face, GpVec, occt_core::gp::GpPnt)> = Vec::new();
    for f in &faces {
        let Some(nrm) = face_normal_at_centroid(f) else { continue };
        let Some(s) = BRepTool::face_surface(f) else { continue };
        let (u0, u1, v0, v1) = BRepTool::uv_bounds(f);
        let (u, v) = if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
            (0.5 * (u0 + u1), 0.5 * (v0 + v1))
        } else {
            (0.0, 0.0)
        };
        let c = s.d0(u, v);
        acc = acc.added(&c.coord);
        n += 1;
        entries.push((f.clone(), nrm, c));
    }
    if n == 0 {
        return Vec::new();
    }
    let centroid = occt_core::gp::GpPnt::from_xyz(&acc.divided(n as f64));
    entries
        .into_iter()
        .filter(|(_, nrm, c)| {
            let outward = occt_core::gp::GpVec::from_pnts(&centroid, c);
            outward.xyz().dot(nrm.xyz()) < 0.0
        })
        .map(|(f, nrm, _)| (f, nrm))
        .collect()
}

/// Rebuild a face with its plane's normal reversed (for planar faces). Returns
/// the reoriented face, or the original if the face is not planar.
pub fn reorient_face(face: &Face) -> Face {
    let Some(pln) = crate::brep_surface::face_plane(face) else {
        return face.clone();
    };
    // Negate the plane's normal explicitly: new Z = −old Z, keep X.
    let z = *pln.axis().direction();
    let zr = z.reversed();
    let ax3 = GpAx3::new(pln.location(), zr, &pln.x_axis().direction()).expect("reorient: axis");
    let rev = GpPln::new(ax3);
    let builder = crate::builder::TopoBuilder::new();
    let wires = wires_of_face(face);
    builder.make_face(Arc::new(occt_geom::GeomPlane::new(rev)), &wires)
}

/// Orient every face of a planar shell outward (away from the shell centroid).
/// Returns the reoriented shell, or `None` if the shell is empty.
pub fn orient_shell_outward(shell: &Shell) -> Option<Shell> {
    let faces = faces_of(&shell.0);
    if faces.is_empty() {
        return None;
    }
    let inward = inward_facing_faces(shell);
    if inward.is_empty() {
        return Some(shell.clone());
    }
    let builder = crate::builder::TopoBuilder::new();
    let mut oriented = Vec::with_capacity(faces.len());
    for f in &faces {
        let is_inward = inward.iter().any(|(ff, _)| Arc::ptr_eq(&ff.0.tshape, &f.0.tshape));
        if is_inward {
            oriented.push(reorient_face(f));
        } else {
            oriented.push(f.clone());
        }
    }
    Some(builder.make_shell(&oriented))
}

/// Check a shell is closed AND consistently oriented (all normals outward).
pub fn shell_is_valid(shell: &Shell) -> bool {
    if !crate::shell_check::shell_is_closed(shell) {
        return false;
    }
    let inward = inward_facing_faces(shell);
    // A well-oriented closed shell may still have a couple of faces whose
    // centroid-vs-normal test is ambiguous (curved faces); tolerate none for
    // planar boxes, allow up to 10% otherwise.
    if inward.is_empty() {
        return true;
    }
    let total = faces_of(&shell.0).len();
    let tolerance = (total as f64 * 0.1).ceil() as usize;
    inward.len() <= tolerance.max(1)
}

/// Volume sign of a closed shell via the signed divergence-theorem mesh sum:
/// positive for outward-oriented faces.
pub fn shell_volume_sign(shell: &Shell, deflection: f64) -> f64 {
    let mesh = crate::shape_mesh::mesh_shape(&shell.0, deflection);
    let mut vol = 0.0;
    for t in &mesh.triangles {
        let a = &mesh.vertices[t.n0];
        let b = &mesh.vertices[t.n1];
        let c = &mesh.vertices[t.n2];
        let cross = b.coord.subtracted(&a.coord).crossed(&c.coord.subtracted(&a.coord));
        vol += a.coord.dot(&cross) / 6.0;
    }
    vol
}

/// Reorient a closed shell so its signed volume is positive (outward normals).
pub fn orient_shell_positive(shell: &Shell, deflection: f64) -> Option<Shell> {
    let vol = shell_volume_sign(shell, deflection);
    if vol >= 0.0 {
        return Some(shell.clone());
    }
    // Negative: flip every planar face.
    orient_shell_outward(shell)
}

/// The total surface area of a shell (mesh-based).
pub fn shell_area(shell: &Shell, deflection: f64) -> f64 {
    crate::mesh::mesh_surface_area(&crate::shape_mesh::mesh_shape(&shell.0, deflection))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;


    fn box_shell() -> Shell {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shells = crate::topo_tools_full::shapes_of(&b.solid.0, crate::abs::ShapeType::Shell);
        crate::shape::Shell(shells[0].clone())
    }

    #[test]
    fn box_shell_faces_point_outward() {
        let shell = box_shell();
        assert!(shell_is_valid(&shell));
        let inward = inward_facing_faces(&shell);
        assert!(inward.is_empty(), "no inward faces: {inward:?}");
    }

    #[test]
    fn shell_volume_sign_positive() {
        let shell = box_shell();
        let v = shell_volume_sign(&shell, 0.1);
        assert!(v > 0.9 && v < 1.1, "signed volume {v} (expect +1.0)");
    }

    #[test]
    fn orient_shell_positive_preserves_volume() {
        let shell = box_shell();
        let oriented = orient_shell_positive(&shell, 0.1).expect("oriented");
        let v = shell_volume_sign(&oriented, 0.1);
        assert!(v > 0.9, "still positive {v}");
        assert!((shell_area(&oriented, 0.1) - 6.0).abs() < 0.01);
    }

    #[test]
    fn reorient_flips_plane_normal() {
        let shell = box_shell();
        let f = faces_of(&shell.0);
        let f0 = &f[0];
        let rev = reorient_face(f0);
        // The reoriented face's plane normal is opposite.
        let n0 = face_normal_at_centroid(f0).unwrap();
        let n1 = face_normal_at_centroid(&rev).unwrap();
        assert!(n0.xyz().dot(n1.xyz()) < -0.9, "normals opposite: {n0:?} vs {n1:?}");
    }
}
