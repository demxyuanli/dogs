//! `BRepApprox` — the `BRep`-face instantiation of the `ApproxInt` templates.
//!
//! Source package: `ModelingAlgorithms/TKTopAlgo/BRepApprox` (OCCT 8.0.0).
//! Every `.cxx` there is a `#define ThePSurface BRepAdaptor_Surface … #include
//! <ApproxInt_….gxx>` block, so the package body lives in the generic bodies of
//! [`crate::approx_int`] and this module holds only the concrete pieces:
//!
//! | OCCT | lines | here |
//! |---|---|---|
//! | `BRepApprox_ApproxLine.cxx:24-94` | 94 | [`approx_line`] |
//! | `BRepApprox_SurfaceTool.cxx:19-40`, `.lxx:22-251` | 291 | [`surface_tool`] |
//! | `BRepApprox_TheMultiLineOfApprox_0.cxx:16-29` | 29 | [`multi_line`] |
//! | `BRepApprox_ThePrmPrmSvSurfacesOfApprox_0.cxx:19-47` | 47 | [`prm_prm_sv_surfaces`] |
//! | `BRepApprox_Approx_0.cxx:16-98` | 98 | [`approx`] |
//! | `BRepApprox_TheMultiLineToolOfApprox.hxx:24-141` | 141 | `crate::approx_int::multi_line_tool` (the `LineTool` impl on the multiline) |
//! | `BRepApprox_TheComputeLineOfApprox_0.cxx:99` / `…TheComputeLineBezierOfApprox_0.cxx:69` | 2 | the `Approx_BSplComputeLine.gxx` / `Approx_ComputeLine.gxx` instantiations, ported generically as `occt_geom::approx::bspl_compute_line` / `compute_line` over `LineTool` |
//! | `BRepApprox_TheImpPrmSvSurfacesOfApprox_0.cxx:19-44` | 44 | UNPORTED: needs `ApproxInt_ImpPrmSvSurfaces` |
//! | `BRepApprox_MyBSplGradient*` / `...ParLeastSquare*` / `...ResConstraint*` / `...ParFunction*` / `...Gradient_BFGS*` / `...MyGradientbis*` | 1,232 | instantiation-only wrappers of the `AppDef` / `AppParCurves` templates; their bodies are `occt_geom::app_par_curves` (`LeastSquare`, `ResConstraint`, `Gradient`, `BspGradient`, `ParFunction`, `BspFunction`) generalised over `LineTool` |
//!
//! The `TheInt2S` / `TheFunctionOfTheInt2S` wrappers of the package instantiate
//! `IntPatch` / `AppParCurves` templates over the same line, in the same way;
//! `ThePrmPrmSvSurfaces` reuses [`crate::intwalk::TheInt2S`] instead of its own
//! `TheInt2S` (`prm_prm_sv_surfaces.rs` header explains why).

pub mod approx;
pub mod approx_line;
pub mod imp_prm_initial_guess;
pub mod imp_prm_sv_surfaces;
pub mod imp_prm_tangent_processing;
pub mod multi_line;
pub mod prm_prm_sv_surfaces;
pub mod surface_tool;
pub mod zer_imp_func;

pub use approx::BRepApproxApprox;
pub use approx_line::ApproxLine;
pub use imp_prm_sv_surfaces::ImpPrmSvSurfaces;
pub use multi_line::BRepMultiLine;
pub use prm_prm_sv_surfaces::PrmPrmSvSurfaces;
pub use surface_tool::SurfaceTool;
pub use zer_imp_func::ZerImpFunc;
