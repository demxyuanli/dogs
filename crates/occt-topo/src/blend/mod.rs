//! `Blend` — the blending package of `TKFillet`.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/` (OCCT 8.0.0).
//!
//! Ported so far:
//! - [`point`]: `Blend_Point` (`Blend_Point.hxx:31-308`, `.lxx:19-200`,
//!   `.cxx:21-483`), the point of a blending line,
//! - [`status`]: `Blend_Status` / `Blend_DecrochStatus`
//!   (`Blend_Status.hxx:19-32`, `Blend_DecrochStatus.hxx:19-28`),
//! - [`app_function`]: `Blend_AppFunction` (`hxx:39-192`, `.cxx:20-23`), the
//!   deferred function an approximation driver consumes,
//! - [`function`]: `Blend_Function` (`hxx:44-121`, `.cxx:18-58`), two surfaces
//!   plus a guide line,
//! - [`cs_function`]: `Blend_CSFunction` (`hxx:44-199`, `.cxx:18-52`), a surface
//!   plus a curve plus a guide line,
//! - [`surf_rst_function`]: `Blend_SurfRstFunction` (`hxx:38-223`, `.cxx:19-35`),
//!   a surface plus a pcurve plus a guide line,
//! - [`rst_rst_function`]: `Blend_RstRstFunction` (`hxx:42-233`, `.cxx:19-35`),
//!   two pcurves plus a guide line,
//! - [`func_inv`]: `Blend_FuncInv` (`hxx:37-94`, `.cxx:18-21`),
//! - [`inv_functions`]: `Blend_SurfPointFuncInv` / `Blend_SurfCurvFuncInv` /
//!   `Blend_CurvPointFuncInv` (`hxx:31-88`, `.cxx:21-25`),
//! - [`blend_func`]: the `BlendFunc` package helpers `GetShape` /
//!   `GetMinimalWeights` / `NextShape` and `BlendFunc_SectionShape`
//!   (`ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc.cxx:40-152`,
//!   `BlendFunc_SectionShape.hxx:20-26`).
//!
//! Convention shared by every trait here: OCCT writes the body of a redefined
//! base method inside an intermediate class, and Rust sub-traits cannot supply
//! that body, so the implementor writes it. Each trait's doc comment names the
//! OCCT `.cxx` line of the required body. See [`app_function`].
//!
//! UNPORTED, with their OCCT location (`Blend/` unless noted):
//! - `Blend_CurvPointRadInv` (`Blend_CurvPointRadInv.cxx`) — the concrete
//!   inverse for a constant radius between a curve and a point,
//! - `Blend_RstRstEval` (`Blend_RstRstEval.cxx`) — the evaluator that switches
//!   between `BlendFunc_*` and `BRepBlend_*` solutions,
//! - `Blend_AppSurface` (`Blend_AppSurface.cxx`) — the approximation surface
//!   driver (its `AppBlend_Approx` base is ported separately),
//! - the `BlendFunc_*` solvers of the `BlendFunc` package
//!   (`ModelingAlgorithms/TKFillet/BlendFunc/`), which implement the traits
//!   above; the package helpers `GetShape` / `GetMinimalWeights` /
//!   `NextShape` are ported in [`blend_func`], `ComputeNormal` /
//!   `ComputeDNormal` still wait for `CSLib`.

pub mod app_function;
pub mod blend_func;
pub mod cs_function;
pub mod func_inv;
pub mod function;
pub mod inv_functions;
pub mod point;
pub mod rst_rst_function;
pub mod status;
pub mod surf_rst_function;

pub use app_function::BlendAppFunction;
pub use blend_func::{get_minimal_weights, get_shape, next_shape, BlendFuncSectionShape};
pub use cs_function::BlendCsFunction;
pub use func_inv::BlendFuncInv;
pub use function::BlendFunction;
pub use inv_functions::{BlendCurvPointFuncInv, BlendSurfCurvFuncInv, BlendSurfPointFuncInv};
pub use point::BlendPoint;
pub use rst_rst_function::BlendRstRstFunction;
pub use status::{BlendDecrochStatus, BlendStatus};
pub use surf_rst_function::BlendSurfRstFunction;
