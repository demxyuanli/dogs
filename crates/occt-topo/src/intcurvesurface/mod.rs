//! Exact curve–surface intersection.
//!
//! Port of `IntCurveSurface_HInter` (TKGeomAlgo) at the level the precise
//! boolean needs: given a parametric curve and a surface, compute every point
//! (and, for the coincidence degeneracy, segment) where the curve meets the
//! surface.
//!
//! Dispatch mirrors the OCCT `PerformBounds` → `PerformConicSurf` /
//! `InternalPerform` split:
//!
//! - **Analytic path** for conic curves against quadric surfaces, mirroring
//!   `PerformConicSurf{Line,Circle,Ellipse,Parabola,Hyperbola}`
//!   (`IntCurveSurface_Inter.pxx:523-966`): every conic/quadric arm hands the
//!   conic and the surface quadric to `IntAna_IntConicQuad`
//!   (`IntAna_IntConicQuad.cxx`), except the line/torus arm, which runs
//!   `IntAna_IntLinTorus`. Tangency / coincidence produces a zero-length
//!   tangent point or an `On` segment.
//! - **Polyhedron path** (OCCT's general arm) for non-quadric surfaces
//!   (torus, B-spline, offset, ...): `DecomposeSurfaceIntervals` →
//!   `SamplePars` polygon + `Adaptor3d_HSurfaceTool::NbSamplesU/V` polyhedron →
//!   `Intf_InterferencePolygonPolyhedron` → `SectionPointToParameters` start
//!   points → `IntCurveSurface_TheExactHInter` (`math_FunctionSetRoot`). See
//!   `polyhedron_curve_surface` and the `polygon_utils` / `interference` /
//!   `cs_function` / `exact_inter` / `sorted_points` modules.
//! - **Quadric arm** for a non-conic curve against a plane / cylinder / cone /
//!   sphere (`InternalPerformCurveQuadric`): `math_FunctionAllRoots` on the
//!   signed distance `Q(w)` over each C1 interval of the curve
//!   (`IntCurveSurface_TheQuadCurvExactHInter` +
//!   `IntSurf_Quadric::Distance/Gradient/ValAndGrad`). See `quadric_exact`.
//!
//! Trait objects (`Arc<dyn Curve>` / `Arc<dyn Surface>`) cannot be downcast, so
//! curve dispatch classifies curves by geometric invariants — the established
//! pattern in this port (see `brep_surface::classify_surface`). Surface dispatch
//! instead uses the faithful `Adaptor3d_Surface::GetType()` classifier
//! (`geom_bnd_lib_surface3d::surface_kind`, mirroring
//! `GeomBndLib_Surface::initFromSurface`), so a B-spline patch that is only
//! geometrically a quadric takes the polyhedron arm, exactly as OCCT does.
mod prelude {

    pub(crate) use std::cmp::Ordering;

    pub(crate) use occt_core::elib::clib::in_period;
    pub(crate) use occt_core::elib::slib;
    pub(crate) use occt_core::gp::{GpAx3, GpDir, GpLin, GpPnt, GpTorus, GpVec};
    pub(crate) use occt_geom::intana::{line_torus_intersect, IntAnaIntConicQuad, IntAnaQuadric};
    pub(crate) use occt_geom::{Curve, Surface};

    pub(crate) use crate::brep_surface::surface_normal;

}


mod types;
mod solvers;
mod polygon;
mod polyhedron;
mod section_point_params;
mod polygon_utils;
mod interference;
mod cs_function;
mod exact_inter;
mod sorted_points;
mod quadric_exact;
pub use types::*;
pub use solvers::*;
pub use polygon::{ThePolygonOfHInter, ThePolygonToolOfHInter};
pub use polyhedron::{ThePolyhedronOfHInter, ThePolyhedronToolOfHInter};
pub use section_point_params::section_point_to_parameters;
pub use polygon_utils::{
    nb_samples, sample_pars, surface_nb_samples_u, surface_nb_samples_u_range,
    surface_nb_samples_v, surface_nb_samples_v_range, unique_knot_count, MY_MIN_PNTS,
};
pub use interference::TheInterferenceOfHInter;
pub use cs_function::TheCSFunctionOfHInter;
pub use exact_inter::{TheExactHInter, THE_TOLTANGENCY};
pub use sorted_points::{
    clamp_uv_parameters, collect_interference_points, compute_append_point, compute_transitions,
    decompose_surface_intervals, process_sorted_points, sort_start_points, SortedStartPoints,
    UVBounds, THE_TOLERANCE_ANGULAIRE,
};
pub use quadric_exact::{
    int_surf_quadric_of, TheQuadCurvExactHInter, TheQuadCurvFuncOfTheQuadCurvExactHInter,
    EPSDIST, EPSNUL, EPSX,
};

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
