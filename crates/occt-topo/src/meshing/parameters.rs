//! Meshing parameters — port of `IMeshTools_Parameters` + `IMeshTools_MeshAlgoType`.
//!
//! `MeshParameters` carries the full set of deflection/quality switches the
//! BRepMesh pipeline consumes. Defaults mirror OCCT's `IMeshTools_Parameters`
//! constructor exactly.

/// Built-in 2D Delaunay triangulation algorithms.
///
/// Maps to `IMeshTools_MeshAlgoType`. `Delaunay` is the classic Watson
/// triangulator and the effective default (`IMeshTools_MeshAlgoType_DEFAULT`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshAlgoType {
    /// No algorithm selected / unknown factory.
    Unknown,
    /// Watson 2D Delaunay triangulation (`BRepMesh_MeshAlgoFactory`).
    Delaunay,
    /// Divide-and-Fix PreOscar Delaunay.
    DFPreOscar,
    /// Divide-and-Fix Presta Delaunay.
    DFPresta,
    /// Delabella Delaunay triangulation (`BRepMesh_DelabellaMeshAlgoFactory`).
    Delabella,
}

/// Structure storing meshing parameters. Source: `IMeshTools_Parameters.hxx`.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshParameters {
    /// 2D Delaunay triangulation algorithm factory to use.
    pub mesh_algo: MeshAlgoType,
    /// Angular deflection used to tessellate the boundary edges.
    pub angle: f64,
    /// Linear deflection used to tessellate the boundary edges.
    pub deflection: f64,
    /// Angular deflection used to tessellate the face interior.
    pub angle_interior: f64,
    /// Linear deflection used to tessellate the face interior.
    pub deflection_interior: f64,
    /// Minimum size limiting triangle edge length to avoid amplification
    /// on distorted curves and surfaces.
    pub min_size: f64,
    /// Switches on/off multi-thread computation.
    pub in_parallel: bool,
    /// Relative deflection: per-edge deflection is `deflection * edge size`.
    pub relative: bool,
    /// Take internal face vertices into account in triangulation.
    pub internal_vertices_mode: bool,
    /// Check deviation of triangulation and interior of the face.
    pub control_surface_deflection: bool,
    /// Apply the `ControlSurfaceDeflection` check to all surface types,
    /// including analytical ones.
    pub enable_control_surface_deflection_all_surfaces: bool,
    /// Clean temporary data model when the algorithm finishes.
    pub clean_model: bool,
    /// Locally adjust min size depending on edge size (disabled by default).
    pub adjust_min_size: bool,
    /// Use shape tolerances for computing face deflection (disabled by default).
    pub force_face_deflection: bool,
    /// Allow decreasing the quality of an existing generated mesh.
    pub allow_quality_decrease: bool,
}

impl Default for MeshParameters {
    /// OCCT `IMeshTools_Parameters` default constructor values.
    fn default() -> Self {
        Self {
            mesh_algo: MeshAlgoType::Delaunay,
            angle: 0.5,
            deflection: 0.001,
            angle_interior: -1.0,
            deflection_interior: -1.0,
            min_size: -1.0,
            in_parallel: false,
            relative: false,
            internal_vertices_mode: true,
            control_surface_deflection: true,
            enable_control_surface_deflection_all_surfaces: false,
            clean_model: true,
            adjust_min_size: false,
            force_face_deflection: false,
            allow_quality_decrease: false,
        }
    }
}

impl MeshParameters {
    /// Factor used to compute the default value of `min_size` from deflection.
    /// Source: `IMeshTools_Parameters::RelMinSize()`.
    pub fn rel_min_size() -> f64 {
        0.1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parameters_defaults_match_occt() {
        let p = MeshParameters::default();
        assert_eq!(p.mesh_algo, MeshAlgoType::Delaunay);
        assert_eq!(p.angle, 0.5);
        assert_eq!(p.deflection, 0.001);
        assert_eq!(p.angle_interior, -1.0);
        assert_eq!(p.deflection_interior, -1.0);
        assert_eq!(p.min_size, -1.0);
        assert!(!p.in_parallel);
        assert!(!p.relative);
        assert!(p.internal_vertices_mode);
        assert!(p.control_surface_deflection);
        assert!(!p.enable_control_surface_deflection_all_surfaces);
        assert!(p.clean_model);
        assert!(!p.adjust_min_size);
        assert!(!p.force_face_deflection);
        assert!(!p.allow_quality_decrease);
    }

    #[test]
    fn rel_min_size_factor() {
        assert_eq!(MeshParameters::rel_min_size(), 0.1);
    }

    #[test]
    fn parameters_are_mutable() {
        let mut p = MeshParameters::default();
        p.mesh_algo = MeshAlgoType::Delabella;
        p.deflection = 0.01;
        p.relative = true;
        assert_eq!(p.mesh_algo, MeshAlgoType::Delabella);
        assert_eq!(p.deflection, 0.01);
        assert!(p.relative);
    }
}
