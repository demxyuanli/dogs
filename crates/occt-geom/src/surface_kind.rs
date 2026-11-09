//! Coarse surface classification shared by the geometry and topology layers.
//!
//! `Adaptor3d_Surface::GetType()` (`GeomAdaptor_Surface.cxx`) is the single
//! source of truth for "is this surface analytic, and which analytic family".
//! Both `pcurve_full` (projection dispatch) and the topology registry (debug
//! surface descriptions) need it, so it lives here rather than in either of
//! them.

use std::sync::Arc;

use crate::surface::Surface;

/// Coarse classification of a surface's analytic type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

/// Classify a surface by its concrete geometry type.
///
/// This mirrors `Adaptor3d_Surface::GetType()`: `ProjLib_ProjectedCurve.cxx`
/// (line 375) selects its analytic projection from `GeomAbs_SurfaceType`
/// (`case GeomAbs_Plane:` ... `case GeomAbs_Torus:`), and
/// `ShapeFix_Edge::FixAddPCurve` (`cxx:487`) tests `surf->IsKind(Geom_Plane)`.
/// Both are *type* tests, so a surface that is only geometrically a cylinder /
/// cone / torus (a `Geom_BSplineSurface` patch) must stay `Other` and be
/// projected by the general path; treating it as analytic would feed the
/// sampled analytic frame a parameterization it does not have.
pub fn classify_surface_kind(s: &dyn Surface) -> SurfaceKind {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return classify_surface_kind(basis.as_ref());
    }
    if s.gp_pln().is_some() {
        SurfaceKind::Plane
    } else if s.gp_cylinder().is_some() {
        SurfaceKind::Cylinder
    } else if s.gp_cone().is_some() {
        SurfaceKind::Cone
    } else if s.gp_sphere().is_some() {
        SurfaceKind::Sphere
    } else if s.gp_torus().is_some() {
        SurfaceKind::Torus
    } else {
        SurfaceKind::Other
    }
}

/// Scratch type kept so callers that need the concrete basis can reuse the
/// trimmed-surface unwrapping without re-implementing the recursion.
pub fn trimmed_basis(s: &dyn Surface) -> Option<Arc<dyn Surface>> {
    s.rectangular_trimmed_basis()
}

/// `GeomAdaptor_Surface::GetType() == GeomAbs_Plane`, i.e. [`classify_surface_kind`]
/// says `Plane`.
///
/// This is the predicate OCCT uses wherever plane behaviour is selected:
/// `BOPAlgo_PaveFiller_6.cxx:380` (`aBAS1.GetType() == GeomAbs_Plane` gating
/// `CheckPlanes`), `_6.cxx:404`, `_5.cxx:393` and `_5.cxx:1038` (the
/// `bLinePlane` / `bUseAddTol` tests), `IntTools_FaceFace.cxx:395` (before
/// `PerformPlanes`) and `Extrema_ExtPS.cxx:279`.
pub fn is_plane_type(s: &dyn Surface) -> bool {
    classify_surface_kind(s) == SurfaceKind::Plane
}
