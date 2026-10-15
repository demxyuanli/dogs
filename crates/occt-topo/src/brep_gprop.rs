//! Mass properties of a topological shape.
//!
//! `BRepGProp` (TKTopAlgo) entry points: surface area, volume, centroid and
//! density-weighted mass. The integration itself is the exact `BRepGProp_Gauss`
//! port in [`crate::brep_gprop_full`] (`SurfaceProperties` -> `BRepGProp_Sinert`
//! / `BRepGProp_Gauss::Sinert`; `VolumeProperties` -> `BRepGProp_Vinert`), so
//! area/volume/centroid no longer depend on a deflection or a face grid.
//!
//! The one fallback is the `NoSurf` arm of
//! `BRepGProp::VolumeProperties` (`BRepGProp.cxx:344-360`): a face whose
//! analytic integration cannot be built is taken from its triangulation through
//! `BRepGProp_MeshProps`, here over the `BRepMesh` tessellation of
//! [`crate::shape_mesh::mesh_shape`]. OCCT orders the same three-way test
//! (`NoSurf`/`NoTri`/surface) per face; the port applies it per shape, so an
//! integrator failure falls back to the triangulation for the whole shape.

use occt_core::gp::GpPnt;
use occt_core::gprop::centroid_of_points;
use occt_core::poly::triangulation::{Triangle, Triangulation};

use crate::brep_gprop_full::{mesh_props, surface_properties, volume_properties};
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
/// surface, empty, or degenerate geometry).
pub fn shape_mass_properties(shape: &TopoShape, deflection: f64) -> MassProperties {
    // `BRepGProp::SurfaceProperties` / `VolumeProperties` over the analytic
    // boundary: exact, no deflection.
    if let (Ok(v), Ok((s, _))) = (volume_properties(shape), surface_properties(shape)) {
        let volume = v.mass().abs();
        return MassProperties {
            volume,
            area: s.mass().abs(),
            center: v.center(),
            centroid_ok: volume > 1e-7,
        };
    }
    mesh_mass_properties(shape, deflection)
}

/// `BRepGProp_MeshProps` arm (`BRepGProp.cxx:359`): the properties of a shape
/// whose faces carry no integrable surface come from its triangulation.
fn mesh_mass_properties(shape: &TopoShape, deflection: f64) -> MassProperties {
    let m = crate::shape_mesh::mesh_shape(shape, deflection);
    let tri = Triangulation::new(m.vertices.clone(), m.triangles.iter().copied().collect::<Vec<Triangle>>());
    let origin = GpPnt::zero();
    let v = mesh_props(&tri, &origin, false, true);
    let s = mesh_props(&tri, &origin, false, false);
    let volume = v.mass().abs();
    MassProperties {
        volume,
        area: s.mass().abs(),
        center: v.center(),
        centroid_ok: volume > 1e-7,
    }
}

/// Surface area of a shape, from the exact boundary integration
/// (`BRepGProp::SurfaceProperties`); `deflection` only reaches the
/// triangulation fallback.
pub fn surface_area(shape: &TopoShape, deflection: f64) -> f64 {
    shape_mass_properties(shape, deflection).area
}

/// Volume of a closed solid, from the exact boundary integration
/// (`BRepGProp::VolumeProperties`); `deflection` only reaches the
/// triangulation fallback.
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
