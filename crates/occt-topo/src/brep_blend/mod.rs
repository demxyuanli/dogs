//! `BRepBlend` — the value layer of the `TKFillet` blending package.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/` (OCCT 8.0.0).
//!
//! Ported so far:
//! - [`point_on_rst`]: `BRepBlend_PointOnRst` (`BRepBlend_PointOnRst.hxx:34-75`,
//!   `.lxx:16-35`, `.cxx:22-47`),
//! - [`extremity`]: `BRepBlend_Extremity` (`BRepBlend_Extremity.hxx:42-155`,
//!   `.lxx:16-87`, `.cxx:26-142`),
//! - [`line`]: `BRepBlend_Line` (`BRepBlend_Line.hxx:37-186`, `.lxx:16-130`,
//!   `.cxx:25-58`),
//! - [`blend_tool`]: `BRepBlend_BlendTool` (`hxx:35-90`, `.cxx:17-142`,
//!   `.lxx:17-58`), the geometric queries of the walking builders,
//! - [`curve_tool`]: the `NbSamples` body shared by `BRepBlend_HCurveTool`
//!   (`hxx:41-121`, `.cxx:25-57`) and `BRepBlend_HCurve2dTool`
//!   (`hxx:43-125`, `.cxx:25-57`); their other members are branch-free
//!   forwarders to the curve adaptor and are therefore the `Curve` / `Curve2d`
//!   trait calls themselves.
//! - [`inv`]: the five inverse functions the walking builders hand to
//!   `math_FunctionSetRoot`, `BRepBlend_SurfPointConstRadInv`
//!   (`hxx:32-84`, `.cxx:17-276`), `BRepBlend_SurfPointEvolRadInv`
//!   (`.cxx:17-269`), `BRepBlend_SurfCurvConstRadInv` (`.cxx:17-284`),
//!   `BRepBlend_SurfCurvEvolRadInv` (`.cxx:17-303`) and
//!   `BRepBlend_CurvPointRadInv` (`.cxx:17-130`).
//! - [`walking`]: `BRepBlend_Walking` (`hxx:41-186`, `.cxx:95-2769`), the
//!   surface/surface walking driver: `Perform`, `PerformFirstSection`,
//!   `TestArret`, `Recadre` / `ArcToRecadre`, `InternalPerform`, `Transition`,
//!   `MakeExtremity` / `MakeSingularExtremity` / `CorrectExtremityOnOneRst`,
//!   `Continu` / `Complete`, `SetDomainsToRecadre` and `AddSingularPoint`.
//! - [`cs_walking`]: `BRepBlend_CSWalking` (`hxx:46-117`, `.cxx:88-1118`), the
//!   surface/curve walking driver: `Perform`, `Complete`, `TestArret`,
//!   `CheckDeflectionOnSurf` / `CheckDeflectionOnCurv`, `InternalPerform`,
//!   `Transition` and `MakeExtremity`.
//! - [`surf_rst_line_builder`]: `BRepBlend_SurfRstLineBuilder`
//!   (`hxx:64-176`, `.lxx:19-50`, `.cxx:137-1777`), the surface / restriction
//!   walking driver: `Perform`, `PerformFirstSection`, `ArcToRecadre`, the
//!   three `Recadre` overloads, `InternalPerform`, `Transition`,
//!   `MakeExtremity`, `TestArret`, `CheckDeflectionOnSurf` /
//!   `CheckDeflectionOnRst`, `CheckInside` and `Complete`.
//! - [`rst_rst_line_builder`]: `BRepBlend_RstRstLineBuilder`
//!   (`hxx:64-205`, `.lxx:19-61`, `.cxx:151-1964`), the restriction /
//!   restriction walking driver: `Perform`, `PerformFirstSection`, the four
//!   `Recadre` overloads, `InternalPerform`, `Transition`, `MakeExtremity`,
//!   `TestArret`, `CheckDeflectionOnRst1` / `CheckDeflectionOnRst2`,
//!   `CheckInside` and `Complete`.
//! - [`surf_rst_const_rad`]: `BRepBlend_SurfRstConstRad`
//!   (`hxx:41-244`, `.cxx:51-829`), the constant-radius blending function on a
//!   surface / restriction pair.
//! - [`surf_rst_evol_rad`]: `BRepBlend_SurfRstEvolRad`
//!   (`hxx:41-248`, `.cxx:113-1166`), its evolutive-radius form driven by a
//!   `Law_Function`.
//! - [`rst_rst_const_rad`]: `BRepBlend_RstRstConstRad`
//!   (`hxx:41-265`, `.cxx:53-981`), the constant-radius blending function on a
//!   restriction pair.
//! - [`rst_rst_evol_rad`]: `BRepBlend_RstRstEvolRad`
//!   (`hxx:41-268`, `.cxx:117-1092`), its evolutive-radius form.
//!
//! UNPORTED, with their OCCT location:
//! - the approximation half, which needs the `AppBlend`/`AppSurf` driver that is
//!   still open in `occt_geom::app_def`: `BRepBlend_AppFunc`,
//!   `BRepBlend_AppFuncRoot`, `BRepBlend_AppFuncRst`, `BRepBlend_AppFuncRstRst`,
//!   `BRepBlend_AppSurf`, `BRepBlend_AppSurface` (`BRepBlend_AppSurf.cxx`,
//!   `BRepBlend_AppSurface.cxx`).

pub mod blend_tool;
pub mod cs_walking;
pub mod curve_tool;
pub mod extremity;
pub mod inv;
pub mod line;
pub mod point_on_rst;
pub mod rst_rst_const_rad;
pub mod rst_rst_evol_rad;
pub mod rst_rst_line_builder;
pub mod surf_rst_const_rad;
pub mod surf_rst_evol_rad;
pub mod surf_rst_line_builder;
pub mod walking;

pub use blend_tool::BlendTool;
pub use cs_walking::CsWalking;
pub use curve_tool::{nb_samples, nb_samples2d};
pub use extremity::Extremity;
pub use inv::{
    CurvPointRadInv, SurfCurvConstRadInv, SurfCurvEvolRadInv, SurfPointConstRadInv,
    SurfPointEvolRadInv,
};
pub use line::Line;
pub use point_on_rst::PointOnRst;
pub use rst_rst_const_rad::SurfRstRstConstRad;
pub use rst_rst_evol_rad::SurfRstRstEvolRad;
pub use rst_rst_line_builder::RstRstLineBuilder;
pub use surf_rst_const_rad::SurfRstConstRad;
pub use surf_rst_evol_rad::SurfRstEvolRad;
pub use surf_rst_line_builder::SurfRstLineBuilder;
