//! The `ChFi3d_Builder_C1.cxx` file statics.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx`
//! (OCCT 8.0.0). The corner bodies themselves (`PerformOneCorner`,
//! `PerformMoreThreeCorner`, `PerformIntersectionAtEnd`, `PerformMoreSurfdata`,
//! `IntersectMoreCorner`) are still UNPORTED - they consume
//! `ChFi3d_FilBuilder`'s edge-contact surfaces and `BRepBlend`. This module
//! carries the helpers those bodies and the two ported member functions
//! ([`crate::chfi3d::builder::ChFi3dBuilder::find_face`] /
//! [`crate::chfi3d::builder::ChFi3dBuilder::more_surfdata`]) call:
//!
//! | OCCT | here |
//! |---|---|
//! | `recadre` (`:118-138`) | `update::recadre` |
//! | `Update(fb, pcfb, surf, fi, cp, p2dbout, isfirst, pared, wop, tol)` (`:147-181`) | [`update::update_pcurve_on_pcurve`] |
//! | `Update(fb, ct, fi, cp, p2dbout, isfirst, wop)` (`:191-278`) | [`update::update_surface_curve`] |
//! | `IntersUpdateOnSame(...)` (`:292-340`) | [`update::inter_update_on_same`] |
//! | `Update(face, edonface, surf, fi, cp, isfirst)` (`:349-417`) | [`update::update_curve_on_surface`] |
//! | `ChFi3d_ExtendSurface` (`:421-449`) | [`extend_surface::ch_fi3d_extend_surface`] |
//! | `ComputeCurve2d` (`:456-483`) | [`compute_curve2d::compute_curve2d`] |
//! | `ChFi3d_Recale` (`:487-551`) | [`recale::ch_fi3d_recale`] |
//! | `ChFi3d_SelectStripe` (`:558-588`) | [`select_stripe::ch_fi3d_select_stripe`] |
//! | `cherche_face` (`:1660-1694`) | [`face_queries::cherche_face`] |
//! | `cherche_edge1` (`:1701-1730`) | [`face_queries::cherche_edge1`] |
//! | `containV` (`:1737-1758`) | [`face_queries::contain_v`] |
//! | `containE` (`:1764-1785`) | [`face_queries::contain_e`] |
//! | `IsShrink` (`:1792-1825`) | [`face_queries::is_shrink`] |
//!
//! UNPORTED, with the OCCT range of each body:
//!
//! - `PerformOneCorner` (`:611-1827`),
//! - `PerformMoreThreeCorner` (`:1237-1652`),
//! - `PerformIntersectionAtEnd` (`:1828-3770`),
//! - `PerformMoreSurfdata` (`:3771-4531`),
//! - `IntersectMoreCorner` (`:4708-...`),
//! - `ComputeCurve2d` (`:456-483`): depends on `BRepAlgo_NormalProjection`
//!   (the whole `BRepAlgo_NormalProjection.cxx` + `BRepAlgo_Approx`), which is
//!   not ported; [`compute_curve2d::compute_curve2d`] documents the gap and
//!   answers `None`.
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//!
//! - `occ::handle<Adaptor3d_Surface>` / `occ::handle<Geom_Surface>` /
//!   `occ::handle<Geom_Curve>` / `occ::handle<Geom2d_Curve>` become
//!   `Arc<dyn Surface>` / `Arc<dyn Curve>` / `Arc<dyn Curve2d>`; a null handle
//!   that OCCT would dereference becomes `None` / `false` (a null
//!   `Geom2dAdaptor_Curve` has no point to evaluate),
//! - `Geom2dAdaptor_Curve(pc, f, l)` - a pcurve restricted to a parameter
//!   window - is [`occt_geom2d::trimmed::Geom2dTrimmedCurve`] with the same
//!   `[f, l]` window: same `FirstParameter` / `LastParameter` for the extrema
//!   search and the same evaluation inside the window,
//! - `Extrema_LocateExtCC` is [`occt_geom::extrema_cc::ExtremaLocateExtCC`],
//!   `IntCurveSurface_HInter::Perform` is
//!   [`crate::intcurvesurface::perform_curve_surface`], `ChFi3d::IsTangentFaces`
//!   is [`crate::chfi3d::tangent_faces::is_tangent_faces`],
//! - the `OCCT_DEBUG` blocks of `cherche_face` / `containV` / `containE` and
//!   the `std::cout` prints are not translated.

mod compute_curve2d;
mod extend_surface;
mod face_queries;
mod recale;
mod select_stripe;
mod update;

pub use compute_curve2d::compute_curve2d;
pub use extend_surface::ch_fi3d_extend_surface;
pub use face_queries::{cherche_edge1, cherche_face, contain_e, contain_v, is_shrink};
pub use recale::ch_fi3d_recale;
pub use select_stripe::{ch_fi3d_select_stripe, StripeIterator};
pub use update::{
    inter_update_on_same, update_curve_on_surface, update_pcurve_on_pcurve, update_surface_curve,
};
