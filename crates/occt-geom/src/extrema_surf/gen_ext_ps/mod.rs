//! Extrema_GenExtPS - general point/surface extrema engine.
//!
//! Source: Extrema_GenExtPS.cxx (1195 lines) + Extrema_GenExtPS.hxx (181),
//! src/ModelingData/TKGeomBase/Extrema/. Line references below are to the
//! .cxx unless another file is named. The helper classes it needs live here too:
//!
//! * [ExtremaPOnSurf] - Extrema_POnSurf.hxx:26-71.
//! * [ExtremaPOnSurfParams] - Extrema_POnSurfParams.hxx:29-88.
//! * [ExtremaElementType] - Extrema_ElementType.hxx:20-26.
//! * [ExtremaFuncPsNorm] - the myF functional, Extrema_FuncPSNorm.hxx /
//!   .cxx (194 lines); this is the math_FunctionSetWithDerivatives the Newton
//!   refinement runs on.
//!
//! The control flow is ported branch for branch:
//! Initialize (cxx:275-319), GetGridPoints + fillParams (cxx:321-454),
//! ComputeEdgeParameters (cxx:460-526), BuildGrid (cxx:528-753), LengthOfIso /
//! CorrectNbSamples (cxx:755-854), BuildTree (cxx:856-931), FindSolution
//! (cxx:933-952), SetFlag / SetAlgo (cxx:954-966), Perform (cxx:968-1151) and
//! the accessors (cxx:1155-1193).
//!
//! Two implementation devices are **not** control-flow differences:
//!
//! * GeomGridEval_Surface::EvaluateGrid (cxx:569-572, cxx:914-917) is an
//!   evaluator with span caching; the grid it produces is exactly
//!   S.Value(U_i, V_j) for every pair, so the port calls
//!   [crate::surface::Surface::d0] node by node.
//! * The NCollection_UBTree built by BuildTree (cxx:907-930) is an acceleration
//!   structure. Bnd_SphereUBTreeSelectorMin/Max (cxx:77-169) selects the grid
//!   sphere minimizing / maximizing Bnd_Sphere::Distance(P); the port keeps
//!   OCCT's mySphereArray (cxx:911) and performs the equivalent scan in
//!   [ExtremaGenExtPs::perform] instead of maintaining the tree. Bnd_Sphere
//!   itself is occt_core::bnd::BndSphere.
//!
//! Extrema_ExtPExtS / Extrema_ExtPRevS are ported in
//! extrusion_point_extrema.rs / revolution_point_extrema.rs and dispatched from
//! point_surface_extrema.rs (Extrema_ExtPS.cxx:292-343); they drive this engine
//! whenever their basis curve is not analytically computable.

mod engine;
mod func_ps_norm;
mod helpers;
mod point_params;

pub use self::engine::ExtremaGenExtPs;
pub use self::func_ps_norm::ExtremaFuncPsNorm;
pub use self::point_params::{ExtremaElementType, ExtremaPOnSurf, ExtremaPOnSurfParams};
