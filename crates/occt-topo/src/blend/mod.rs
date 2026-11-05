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
//!   `GetMinimalWeights` / `NextShape`, `ComputeNormal` / `ComputeDNormal` and
//!   `BlendFunc_SectionShape` (`ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc.cxx:40-280`,
//!   `BlendFunc_SectionShape.hxx:20-26`),
//! - [`const_rad`]: `BlendFunc_ConstRad` (`hxx:49-190`, `.cxx:35-1374`) and
//!   `BlendFunc_ConstRadInv` (`hxx:27-79`, `.cxx:25-646`), the constant-radius
//!   solvers between two surfaces,
//! - [`evol_rad`]: `BlendFunc_EvolRad` (`hxx:36-282`, `.cxx:100-2036`) and
//!   `BlendFunc_EvolRadInv` (`hxx:25-77`, `.cxx:27-470`), the same two systems
//!   with the radius following a `Law_Function` along the guide line,
//! - [`ruled`]: `BlendFunc_Ruled` (`hxx:34-207`, `.cxx:31-744`) and
//!   `BlendFunc_RuledInv` (`hxx:25-73`, `.cxx:21-508`), the ruled segment
//!   between the two supports (no radius).
//! - [`chamfer`]: `BlendFunc_GenChamfer` (`hxx:31-149`, `.cxx:29-300`),
//!   `BlendFunc_ConstThroat` (`hxx:27-110`, `.cxx:27-281`),
//!   `BlendFunc_ConstThroatWithPenetration` (`hxx:28-100`, `.cxx:27-185`),
//!   `BlendFunc_GenChamfInv` (`hxx:31-107`, `.cxx:23-122`),
//!   `BlendFunc_ConstThroatInv` (`hxx:26-76`, `.cxx:22-137`) and
//!   `BlendFunc_ConstThroatWithPenetrationInv` (`hxx:24-49`, `.cxx:22-199`),
//!   the symmetric / right-angled chamfer solvers,
//! - [`tensor`]: `BlendFunc_Tensor` (`hxx:29-70`, `.lxx:17-24`, `.cxx:21-53`),
//!   the "gradient of gradient" cube and its `Multiply`.
//!
//! Convention shared by every trait here: OCCT writes the body of a redefined
//! base method inside an intermediate class, and Rust sub-traits cannot supply
//! that body, so the implementor writes it. Each trait's doc comment names the
//! OCCT `.cxx` line of the required body. See [`app_function`].
//!
//! No computational `UNPORTED` class is left in the `Blend` package: OCCT
//! 8.0.0's `ModelingAlgorithms/TKFillet/Blend/` holds eleven `.cxx` and ten of
//! them are translated (`Blend_AppFunction.cxx:1-20` → [`app_function`],
//! `Blend_CSFunction.cxx:1-46` → [`cs_function`], `Blend_Function.cxx:1-49` →
//! [`function`], `Blend_Point.cxx:1-464` → [`point`], `Blend_FuncInv.cxx:1-19`
//! + `Blend_CurvPointFuncInv.cxx:1-21` + `Blend_SurfPointFuncInv.cxx:1-21` +
//! `Blend_SurfCurvFuncInv.cxx:1-22` → [`func_inv`] plus the `nb_variables`
//! bodies (`2` for `Blend_CurvPointFuncInv.cxx:21-24`, `3` for
//! `Blend_SurfPointFuncInv.cxx:21-24` and `Blend_SurfCurvFuncInv.cxx:22-25`,
//! the only body each of the three `*FuncInv` subclasses adds) written by the
//! concrete `BRepBlend_*Inv` ports, `Blend_RstRstFunction.cxx:1-31` →
//! [`rst_rst_function`], `Blend_SurfRstFunction.cxx:1-31` →
//! [`surf_rst_function`]).
//!
//! UNPORTED — `Blend_Debug.cxx:1-50`: the two Draw-only globals
//! `Blend_traceDRAWSECT` (`cxx:23-33`) and `Blend_contextNOTESTDEFL`
//! (`cxx:39-49`) with their setters and getters. Both default to `false`; the
//! first only gates `#ifdef OCCT_DEBUG` section dumps
//! (`BRepBlend_SurfRstLineBuilder.cxx:296-945`, `BRepBlend_Walking.cxx:272`),
//! the second is read only inside `#ifdef OCCT_DEBUG`
//! (`BRepBlend_Walking.cxx:2342-2345`, `testdefl = !Blend_GetcontextNOTESTDEFL()`),
//! so outside an `OCCT_DEBUG` build the deflection test always runs. The port
//! keeps that always-on path and does not carry the globals.
//!
//! `Blend_CurvPointRadInv`, `Blend_RstRstEval` and `Blend_AppSurface` are
//! absent from the 8.0.0 tree — they were removed before 8.0.0 and their roles
//! now sit in `BRepBlend_CurvPointRadInv.cxx:1-105`
//! (`brep_blend/.../curv_point_rad_inv.rs`), `BRepBlend_RstRstConstRad.cxx:1-799`
//! / `BRepBlend_RstRstEvolRad.cxx:1-896` and `BRepBlend_AppSurface.cxx:1-124`
//! (`brep_blend/app_surface.rs`, over the `AppBlend_Approx` ported in
//! `occt_geom::geom_fill::app_surf`).
//! The `BlendFunc_*` package is complete: `Chamfer` / `ChamfInv`,
//! `ChAsym` / `ChAsymInv`, `Corde`, `CSConstRad` and `CSCircular` live in
//! [`chamfer`] / [`cs_const_rad`] / [`cs_circular`], `GenChamfer` /
//! `GenChamfInv` and the four throat classes in [`chamfer`],
//! `BlendFunc_Tensor` in [`tensor`], and the package helpers `GetShape` /
//! `GetMinimalWeights` / `NextShape` / `ComputeNormal` / `ComputeDNormal`
//! in [`blend_func`].

pub mod app_function;
pub mod blend_func;
pub mod chamfer;
pub mod const_rad;
pub mod cs_circular;
pub mod cs_const_rad;
pub mod cs_function;
pub mod evol_rad;
pub mod func_inv;
pub mod function;
pub mod inv_functions;
pub mod point;
pub mod rst_rst_function;
pub mod ruled;
pub mod status;
pub mod surf_rst_function;
pub mod tensor;

pub use app_function::BlendAppFunction;
pub use blend_func::{get_minimal_weights, get_shape, next_shape, BlendFuncSectionShape};
pub use chamfer::{
    ChAsym, ChAsymInv, ChamfInv, Chamfer, ConstThroat, ConstThroatInv,
    ConstThroatWithPenetration, ConstThroatWithPenetrationInv, Corde, GenChamfInv, GenChamfer,
};
pub use const_rad::{ConstRad, ConstRadInv};
pub use cs_circular::CsCircular;
pub use cs_const_rad::CsConstRad;
pub use cs_function::BlendCsFunction;
pub use evol_rad::{EvolRad, EvolRadInv};
pub use func_inv::BlendFuncInv;
pub use function::BlendFunction;
pub use inv_functions::{BlendCurvPointFuncInv, BlendSurfCurvFuncInv, BlendSurfPointFuncInv};
pub use point::BlendPoint;
pub use rst_rst_function::BlendRstRstFunction;
pub use ruled::{Ruled, RuledInv};
pub use status::{BlendDecrochStatus, BlendStatus};
pub use surf_rst_function::BlendSurfRstFunction;
pub use tensor::Tensor;
