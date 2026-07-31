//! Mass properties of a topological shape.
//!
//! Approximate port of `BRepGProp` (TKBRep): computes surface area, volume,
//! centroid and density-weighted mass from the shape's geometry. The shape is
//! first tessellated (via the local face-mesh fallback in [`crate::brep_extrema`]
//! until `shape_mesh::mesh_shape` lands), then:
//!
//! * area — sum of triangle areas ([`crate::mesh::mesh_surface_area`]);
//! * volume — signed tetrahedron sum about the origin, `|Σ (1/6)·p0·(p1×p2)|`
//!   (reusing [`occt_core::gprop::inertia::compute_inertia`]);
//! * centroid — volume-weighted average of tetrahedron centroids.
//!
//! Results are exact for planar faces and approximate for curved ones (the
//! error shrinks as `deflection` decreases).

use occt_core::gp::GpPnt;
use occt_core::gprop::centroid_of_points;
use occt_core::gprop::inertia::compute_inertia;
use occt_core::poly::triangulation::Triangle;

use crate::brep_extrema::mesh_faces;
use crate::mesh::{mesh_surface_area, ShapeMesh};
use crate::shape::TopoShape;

/// Mass properties of a shape (or point set).
#[derive(Debug, Clone)]
pub struct MassProperties {
    pub volume: f64,
    pub area: f64,
    pub center: GpPnt,
    pub centroid_ok: bool,
}

/// Full mass properties of a closed solid: surface area, volume and centroid.
///
/// `centroid_ok` is `false` when the shape has no measurable volume (open
/// surface, empty, or degenerate mesh).
pub fn shape_mass_properties(shape: &TopoShape, deflection: f64) -> MassProperties {
    let mesh = mesh_shape(shape, deflection);
    let area = mesh_surface_area(&mesh);
    let tris: Vec<(usize, usize, usize)> = mesh
        .triangles
        .iter()
        .map(|t| (t.n0, t.n1, t.n2))
        .collect();
    let props = compute_inertia(&mesh.vertices, &tris, 1.0);
    let volume = props.mass.abs();
    let centroid_ok = volume > 1e-7;
    MassProperties {
        volume,
        area,
        center: props.center,
        centroid_ok,
    }
}

/// Surface area of a shape, from its tessellation.
pub fn surface_area(shape: &TopoShape, deflection: f64) -> f64 {
    mesh_surface_area(&mesh_shape(shape, deflection))
}

/// Volume of a closed solid, from its tessellation.
pub fn volume(shape: &TopoShape, deflection: f64) -> f64 {
    shape_mass_properties(shape, deflection).volume
}

/// Volume centroid of a closed solid; `None` when the volume is ~0.
pub fn centroid(shape: &TopoShape, deflection: f64) -> Option<GpPnt> {
    let m = shape_mass_properties(shape, deflection);
    if m.centroid_ok {
        Some(m.center)
    } else {
        None
    }
}

/// Mass of a solid with uniform density (`volume × density`).
pub fn solid_density_mass(shape: &TopoShape, density: f64, deflection: f64) -> f64 {
    volume(shape, deflection) * density
}

/// Mass properties of a bare point cloud: mean point, zero volume/area.
pub fn point_mass_properties(points: &[GpPnt]) -> MassProperties {
    MassProperties {
        volume: 0.0,
        area: 0.0,
        center: centroid_of_points(points),
        centroid_ok: false,
    }
}

/// Tessellate a shape into a triangle mesh.
///
/// Uses the local face-grid fallback from `brep_extrema`. When
/// `shape_mesh::mesh_shape` is completed by the concurrent agent, this body
/// can be swapped to delegate to it.
fn mesh_shape(shape: &TopoShape, deflection: f64) -> ShapeMesh {
    let n = ((1.0 / deflection.max(1e-4)).ceil() as usize).clamp(1, 32);
    let mut vertices: Vec<GpPnt> = Vec::new();
    let mut triangles: Vec<Triangle> = Vec::new();
    for (fv, ft) in mesh_faces(shape, n, n) {
        let base = vertices.len();
        vertices.extend(fv);
        for (a, b, c) in ft {
            triangles.push(Triangle::new(base + a, base + b, base + c));
        }
    }
    ShapeMesh {
        vertices,
        triangles,
        source_shape: shape.shape_type(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box;

    #[test]
    fn box_mass_properties() {
        let b = test_box::unit_box();
        let m = shape_mass_properties(&b.solid.0, 0.05);
        assert!((m.area - 6.0).abs() < 0.05, "area {}", m.area);
        assert!((m.volume - 1.0).abs() < 0.05, "volume {}", m.volume);
        assert!(m.centroid_ok);
        assert!((m.center.x() - 0.5).abs() < 0.05, "cx {}", m.center.x());
        assert!((m.center.y() - 0.5).abs() < 0.05, "cy {}", m.center.y());
        assert!((m.center.z() - 0.5).abs() < 0.05, "cz {}", m.center.z());
    }

    #[test]
    fn surface_area_and_volume() {
        let b = test_box::unit_box();
        assert!((surface_area(&b.solid.0, 0.05) - 6.0).abs() < 0.05);
        assert!((volume(&b.solid.0, 0.05) - 1.0).abs() < 0.05);
        let c = centroid(&b.solid.0, 0.05).expect("centroid");
        assert!(c.distance(&GpPnt::new(0.5, 0.5, 0.5)) < 0.05);
    }

    #[test]
    fn density_mass() {
        let b = test_box::unit_box();
        let m = solid_density_mass(&b.solid.0, 2.0, 0.05);
        assert!((m - 2.0).abs() < 0.1, "mass {m}");
    }

    #[test]
    fn point_cloud_center() {
        let pts = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(2.0, 0.0, 0.0)];
        let m = point_mass_properties(&pts);
        assert!((m.center.x() - 1.0).abs() < 1e-12);
        assert_eq!(m.volume, 0.0);
        assert!(!m.centroid_ok);
    }
}
