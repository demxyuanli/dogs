use super::super::prelude::*;

/// Analytic surface type used to pick a range splitter. Source:
/// `GeomAbs_SurfaceType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    BezierSurface,
    BSplineSurface,
    OffsetSurface,
    OtherSurface,
}

/// Classifies a `dyn Surface` into a [`MeshSurfaceType`].
///
/// Exact transcription of `GeomAdaptor_Surface::Load`
/// (`GeomAdaptor_Surface.cxx:422-513`), which compares **`DynamicType`**
/// (exact class, not `IsKind`) in a fixed order: rectangular-trimmed (recurse
/// on the basis), Plane, Cylinder, Cone, Sphere, Torus, SurfaceOfRevolution,
/// SurfaceOfLinearExtrusion, **Bezier** (`cxx:480`, before BSpline),
/// BSpline, Offset, else `GeomAbs_OtherSurface`.
///
/// The previous body guessed Plane/Cylinder/Cone/Sphere/Torus from the
/// periodicity flags plus a sampled constant-radius test (`is_cylinder_like`),
/// which turned any V-bounded cylinder or cone into a `Sphere` and swapped the
/// whole internal node grid (audit A24).
pub fn mesh_surface_type(s: &dyn Surface) -> MeshSurfaceType {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return mesh_surface_type(basis.as_ref());
    }
    if s.gp_pln().is_some() {
        return MeshSurfaceType::Plane;
    }
    if s.gp_cylinder().is_some() {
        return MeshSurfaceType::Cylinder;
    }
    if s.gp_cone().is_some() {
        return MeshSurfaceType::Cone;
    }
    if s.gp_sphere().is_some() {
        return MeshSurfaceType::Sphere;
    }
    if s.gp_torus().is_some() {
        return MeshSurfaceType::Torus;
    }
    if s.is_surface_of_revolution() {
        return MeshSurfaceType::SurfaceOfRevolution;
    }
    if s.is_surface_of_linear_extrusion() {
        return MeshSurfaceType::SurfaceOfExtrusion;
    }
    if s.is_bezier_surface() {
        return MeshSurfaceType::BezierSurface;
    }
    if s.is_bspline_surface() {
        return MeshSurfaceType::BSplineSurface;
    }
    if s.is_offset_surface() {
        return MeshSurfaceType::OffsetSurface;
    }
    MeshSurfaceType::OtherSurface
}

/// Radius of a cylinder. `BRepAdaptor_Surface::Cylinder().Radius()` when the
/// face is a `Geom_CylindricalSurface`; otherwise the chord `(0,0)`--`(π,0)`.
pub(in crate::meshing::range_splitter) fn cylinder_radius(s: &dyn Surface) -> f64 {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return cylinder_radius(basis.as_ref());
    }
    if let Some(c) = s.gp_cylinder() {
        return c.radius();
    }
    s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5
}

/// Radius of a sphere. `BRepAdaptor_Surface::Sphere().Radius()` when the face
/// is a `Geom_SphericalSurface`; otherwise the chord `(0,0)`--`(π,0)`.
pub(in crate::meshing::range_splitter) fn sphere_radius(s: &dyn Surface) -> f64 {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return sphere_radius(basis.as_ref());
    }
    if let Some(sph) = s.gp_sphere() {
        return sph.radius();
    }
    s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5
}

/// `(major, minor)` radii. `BRepAdaptor_Surface::Torus()` when the face is a
/// `Geom_ToroidalSurface`; otherwise the same chord measure as before.
pub(in crate::meshing::range_splitter) fn torus_radii(s: &dyn Surface) -> (f64, f64) {
    if let Some(t) = s.gp_torus() {
        return (t.major_radius(), t.minor_radius());
    }
    let minor = s.d0(0.0, 0.0).distance(&s.d0(0.0, PI)) * 0.5;
    let major_plus_minor = s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5;
    (major_plus_minor - minor, minor)
}

/// Radius of the U-circle of a cone at parameter `v` measured from the surface.
pub(in crate::meshing::range_splitter) fn cone_radius_at(s: &dyn Surface, v: f64) -> f64 {
    s.d0(0.0, v).distance(&s.d0(PI, v)) * 0.5
}

/// Angular step for a circle of `radius` under linear (`lin_deflection`) and
/// angular (`ang_deflection`) deflections, with an optional `min_length` floor.
/// Source: `GCPnts_TangentialDeflection::ArcAngularStep`.
pub fn arc_angular_step(
    radius: f64,
    lin_deflection: f64,
    ang_deflection: f64,
    min_length: f64,
) -> f64 {
    let mut du = 0.0;
    let mut min_size_ang = 0.0;
    if radius > CONFUSION {
        du = (1.0 - lin_deflection / radius).max(0.0);
        if min_length > CONFUSION {
            min_size_ang = (min_length / radius).min(PI * 0.5);
        }
    }
    du = 2.0 * du.acos();
    du.min(ang_deflection).max(min_size_ang)
}

/// Parametric step along U or V that corresponds to a 3D displacement of `tol`.
/// BSpline uses `Geom_BSplineSurface::Resolution` (`cxx:2197-2221`,
/// `GeomAdaptor_Surface.cxx:1877-1880`). Other types keep the FD estimate.
pub(in crate::meshing::range_splitter) fn param_resolution(
    s: &dyn Surface,
    tol: f64,
    is_u: bool,
) -> f64 {
    if let Some((ru, rv)) = s.uv_resolution(tol) {
        return if is_u { ru } else { rv };
    }
    let (r0, r1) = if is_u { s.u_range() } else { s.v_range() };
    let (o0, o1) = if is_u { s.v_range() } else { s.u_range() };
    let mid = |a: f64, b: f64| if a.is_finite() && b.is_finite() { 0.5 * (a + b) } else { 0.0 };
    let u = mid(r0, r1);
    let v = mid(o0, o1);
    let h = 1e-4;
    let p0 = s.d0(u, v);
    let p1 = if is_u { s.d0(u + h, v) } else { s.d0(u, v + h) };
    let dist = p0.distance(&p1);
    if dist > 1e-12 {
        tol * h / dist
    } else {
        tol
    }
}
