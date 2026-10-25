//! `ApproxInt` — the surface-pair evaluation layer the `BRepApprox` fit drives.
//!
//! Source package: `ModelingAlgorithms/TKGeomAlgo/ApproxInt` (OCCT 8.0.0):
//!
//! | OCCT | lines | here |
//! |---|---|---|
//! | `ApproxInt_SvSurfaces.hxx:22-101`, `.cxx:5-19` | 101 | [`sv_surfaces`] |
//! | `ApproxInt_MultiLine.gxx:30-778` | 778 | [`multi_line`] |
//! | `ApproxInt_MultiLineTool.lxx:22-186` | 186 | [`multi_line_tool`] |
//! | `ApproxInt_PrmPrmSvSurfaces.gxx:28-373` | 373 | `crate::brep_approx::prm_prm_sv_surfaces` |
//! | `ApproxInt_KnotTools.cxx:1-848`, `.hxx:38-132` | 898 | [`knot_tools`] |
//! | `ApproxInt_Approx.gxx:1-751` | 751 | [`approx`], [`approx_data`] |
//! | `ApproxInt_ImpPrmSvSurfaces.gxx:1-1020` | 1,020 | `crate::brep_approx::imp_prm_sv_surfaces`, `imp_prm_tangent_processing`, `imp_prm_initial_guess` |
//!
//! OCCT writes these as `.gxx` templates instantiated per surface kind
//! (`#define ThePSurface BRepAdaptor_Surface` … `#include <ApproxInt_….gxx>`).
//! Rust generics carry the same parameterisation, so the generic bodies are
//! translated once here and instantiated by [`crate::brep_approx`]. The one
//! exception is `ApproxInt_ImpPrmSvSurfaces`: its `TheISurface` is always
//! `IntSurf_Quadric` and its only instantiation in the port is the `BRepApprox`
//! one, so it lives beside `BRepApprox` in `crate::brep_approx` (its
//! `.gxx` body is not generic over the surface classes beyond those two).
//!
//! UNPORTED:
//! - `ApproxInt_KnotTools::DefineParType` (`cxx:669-848`): builds a
//!   `GeomInt_TheMultiLineOfWLApprox` over an `IntPatch_WLine` (`cxx:686-699`);
//!   it is not on the `ApproxInt_Approx` path (`gxx:610` calls only
//!   `BuildKnots`).

pub mod approx;
pub mod approx_data;
pub mod knot_tools;
pub mod multi_line;
pub mod multi_line_tool;
pub mod sv_surfaces;

pub use approx::ApproxIntApprox;
pub use approx_data::ApproxData;
pub use multi_line::{LineSource, MultiLine};
pub use sv_surfaces::{SolverFlag, SvSurfaces};
