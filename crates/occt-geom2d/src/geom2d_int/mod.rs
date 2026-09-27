//! Port of the 2D intersection packages of TKGeomAlgo used by
//! `ShapeFix_ComposeShell::SplitByLine`: the `Geom2dInt` curve adaptor and the
//! `Geom2dInt_GInter` dispatch with its Line/Line conic-conic kernel.
//!
//! `IntRes2d`, `IntCurve_IConicTool`, `IntImpParGen` and the `math` root
//! finders are already ported in `occt-core` (`intres2d`, `intcurve`,
//! `intimpargen`, `math_function_all_roots`) and are reused from there.

mod conic_conic;
mod conic_curve;
mod curve_locator;
mod curve_sampling;
mod curve_tool;
mod gen_locate_ext_pc;
mod gfunc_ext_pc;
mod ginter;
mod int_conic_curve;
mod int_conic_conic;
mod proj_p_cur;

pub use conic_conic::*;
pub use conic_curve::*;
pub use curve_locator::*;
pub use curve_sampling::*;
pub use curve_tool::*;
pub use gen_locate_ext_pc::*;
pub use ginter::*;
pub use int_conic_curve::*;
pub use gfunc_ext_pc::*;
pub use proj_p_cur::*;
pub use int_conic_conic::*;
