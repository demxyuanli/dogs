//! Surface classification and p-curve helpers for faces.
//! Source: `GeomAdaptor_Surface.hxx`, `BRep_Tool::Surface`,
//! `GeomAPI_ProjectPointOnSurf.hxx`.
//!
//! OCCT classifies a surface with its `GeomAdaptor_Surface::Load` type tag
//! (`GeomAdaptor_Surface.cxx:422-513`), i.e. the exact `DynamicType` of the
//! surface object; the port's equivalent is the `Surface` type queries
//! (`gp_pln` / `gp_sphere` / `gp_cylinder` / `gp_cone` / `gp_torus`, which also
//! answer through a `Geom_RectangularTrimmedSurface` basis, the documented
//! `GeomAdaptor_Surface.cxx:423-425` recursion). [`classify_surface`],
//! [`is_plane_type`], [`face_is_plane_type`], [`face_plane`] and
//! [`sphere_center`] all use those queries; the 8x8 normal-sampling
//! `is_planar` / `face_is_planar` pair was removed with finding 14 in
//! `specs/_audit/topo-rest.md`.

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_geom::Surface;


use crate::shape::Face;
use crate::tgeometry::GeometryRegistry;

pub use occt_geom::surface_kind::SurfaceKind;

/// Unit normal of a surface at (u, v), robust to surfaces whose `d1` returns
/// zero vectors (some ported surfaces only implement `d0`).
pub fn surface_normal(s: &dyn Surface, u: f64, v: f64) -> GpVec {
    let (_, du, dv) = s.d1(u, v);
    let n = du.xyz().crossed(dv.xyz());
    if n.square_modulus() > 1e-30 {
        let m = n.modulus();
        return GpVec::new(n.x / m, n.y / m, n.z / m);
    }
    // Fallback: central finite differences of d0.
    let eps = 1e-6;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { eps };
    let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { eps };
    let p0 = s.d0(u, v);
    let duv = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
    let dvv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
    let n = duv.xyz().crossed(dvv.xyz());
    let m = n.modulus();
    if m > 1e-30 {
        GpVec::new(n.x / m, n.y / m, n.z / m)
    } else {
        GpVec::zero()
    }
}

/// Classify a surface by its **exact analytic type**.
///
/// Mirrors `GeomAdaptor_Surface::Load` (`GeomAdaptor_Surface.cxx:422-513`),
/// which compares the surface's exact class in a fixed order
/// (rectangular-trimmed → recurse on the basis; Plane, Cylinder, Cone, Sphere,
/// Torus, …; anything else is `GeomAbs_OtherSurface`). The port exposes the same
/// information through the `Surface` type queries, so no sampling is needed.
///
/// The previous body sampled an 8×8 grid and compared normals / center
/// distances with `1e-6` / `1e-4` thresholds, and could only ever report
/// `Plane`, `Sphere` or `Other` (audit A30).
pub fn classify_surface(s: &dyn Surface) -> SurfaceKind {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return classify_surface(basis.as_ref());
    }
    if s.gp_pln().is_some() {
        return SurfaceKind::Plane;
    }
    if s.gp_cylinder().is_some() {
        return SurfaceKind::Cylinder;
    }
    if s.gp_cone().is_some() {
        return SurfaceKind::Cone;
    }
    if s.gp_sphere().is_some() {
        return SurfaceKind::Sphere;
    }
    if s.gp_torus().is_some() {
        return SurfaceKind::Torus;
    }
    SurfaceKind::Other
}

/// `GeomAdaptor_Surface::GetType() == GeomAbs_Plane`, i.e. [`classify_surface`]
/// says `Plane`.
///
/// This is the predicate OCCT uses wherever plane behaviour is selected:
/// `BOPAlgo_PaveFiller_6.cxx:380` (`aBAS1.GetType() == GeomAbs_Plane` gating
/// `CheckPlanes`), `_6.cxx:404`, `_5.cxx:393` and `_5.cxx:1038` (the
/// `bLinePlane` / `bUseAddTol` tests), `IntTools_FaceFace.cxx:395` (before
/// `PerformPlanes`) and `Extrema_ExtPS.cxx:279` (the analytic
/// `Extrema_ExtPElS` projection). It is the exact `GeomAdaptor_Surface::Load`
/// type (`GeomAdaptor_Surface.cxx:422-513`), so a trimmed `Geom_Plane` answers
/// true and a geometrically flat B-spline patch (`GeomAbs_BSplineSurface`)
/// answers false. The previous predicate, [`is_planar`], sampled 8x8 normals
/// and could not tell those two apart (audit `topo-rest.md` finding 14).
pub fn is_plane_type(s: &dyn Surface) -> bool {
    classify_surface(s) == SurfaceKind::Plane
}

/// [`is_plane_type`] for a face, via the geometry registry.
pub fn face_is_plane_type(face: &Face) -> bool {
    match GeometryRegistry::global().face_surface(&face.0) {
        Some(s) => is_plane_type(s.as_ref()),
        None => false,
    }
}

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// If the surface is spherical, return its center.
///
/// `GeomAdaptor_Surface::Sphere()` (`GeomAdaptor_Surface.cxx:1997-2006`), i.e.
/// the `Geom_SphericalSurface` downcast; `None` for every other type. The
/// previous body solved the circumcenter of three sampled points, which is
/// exact for a sphere but also answers for a cylinder, a cone or a B-spline
/// patch, so callers saw a "sphere" where OCCT sees `GeomAbs_OtherSurface`
/// (audit `topo-rest.md` finding 14).
pub fn sphere_center(s: &dyn Surface) -> Option<GpPnt> {
    s.gp_sphere().map(|sp| sp.pos.location())
}

/// The plane underlying a planar face, if it is a `GeomAbs_Plane`.
///
/// The classification is [`is_plane_type`] - `GeomAdaptor_Surface::Plane()`
/// (`GeomAdaptor_Surface.cxx:1863-1870`), the `Geom_Plane` downcast that
/// `BRepGProp_Face.cxx:43-46` and `BRepOffset_*` gate on - which accepts a
/// `Geom_Plane` and rejects a geometrically flat B-spline patch. The frame is
/// read from the surface's own parameterization (`d0` / `d1` at the natural
/// window origin) so it stays consistent with `d0`, the frame every other UV
/// routine in the port uses; on a real `Geom_Plane` that is exactly its
/// placement.
pub fn face_plane(face: &Face) -> Option<GpPln> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    if !is_plane_type(surf.as_ref()) {
        return None;
    }
    let (u0, _u1, v0, _v1) = sample_bounds(surf.as_ref());
    let o = surf.d0(u0, v0);
    let n = surface_normal(surf.as_ref(), u0, v0);
    let d = GpDir::from_vec(&n).ok()?;
    // Pick an in-plane X direction perpendicular to the normal (a fixed axis
    // may be parallel to the normal, which GpAx3 rejects).
    let z_axis = occt_core::gp::GpDir::new(0.0, 0.0, 1.0).unwrap();
    let x_dir = if d.is_normal(&z_axis) {
        z_axis
    } else {
        occt_core::gp::GpDir::new(1.0, 0.0, 0.0).unwrap()
    };
    Some(GpPln::new(GpAx3::new(o, d, &x_dir).ok()?))
}

/// Surface normal at a face's natural (u, v) parameter point.
pub fn face_normal(face: &Face, u: f64, v: f64) -> Option<GpVec> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    Some(surface_normal(surf.as_ref(), u, v))
}

/// Centroid of a face surface sampled on a `nu × nv` grid.
pub fn face_centroid(face: &Face, nu: usize, nv: usize) -> Option<GpPnt> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    let (u0, u1, v0, v1) = sample_bounds(surf.as_ref());
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut count = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let v = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            acc = acc.added(&surf.d0(u, v).coord);
            count += 1.0;
        }
    }
    if count == 0.0 {
        None
    } else {
        Some(GpPnt::from_xyz(&acc.divided(count)))
    }
}


/// Face UV bounds (`BRep_Tool::UVBounds(F, UMin, UMax, VMin, VMax)`,
/// `BRep_Tool.cxx`), i.e. the box of the face's pcurves, not the surface's
/// natural range:
///
///   BRep_Tool::UVBounds -> BRepTools::UVBounds(F, B) -> AddUVBounds
/// (`BRepTools.cxx:64-75`, `:126-160`); a face with no usable pcurve box falls
/// back to the surface bounds (`BRepTools.cxx:141-153`). This is also what
/// `BRepAdaptor_Surface(F)` restricts to (`BRepAdaptor_Surface.cxx:72-75`,
/// `BRepTools::UVBounds`), which is the adaptor `IntTools_FaceFace::Perform`
/// loads into its `GeomAdaptor_Surface` (`IntTools_FaceFace.cxx`).
pub fn face_uv_bounds(face: &Face) -> (f64, f64, f64, f64) {
    crate::brep_tools::uv_bounds(face)
}
