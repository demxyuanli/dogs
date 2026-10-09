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
mod dist_between_pcurves;
mod exact_intersection_point;
mod gen_locate_ext_pc;
mod gfunc_ext_pc;
mod ginter;
mod int_conic_curve;
mod int_conic_conic;
mod int_conic_conic_ana_bounds;
mod int_conic_conic_circle_circle;
mod int_conic_conic_circle_conics;
mod int_conic_conic_ellips_hypr;
mod int_conic_conic_line_circle;
mod int_conic_conic_line_ellipse;
mod int_conic_conic_line_parab_hypr;
mod int_conic_conic_parab;
mod int_conic_conic_tool;
mod int_poly_poly_gen;
mod polygon2d;
mod proj_p_cur;

pub use conic_conic::*;
pub use conic_curve::*;
pub use curve_locator::*;
pub use curve_sampling::*;
pub use curve_tool::*;
pub use dist_between_pcurves::*;
pub use exact_intersection_point::*;
pub use gen_locate_ext_pc::*;
pub use ginter::*;
pub use int_conic_curve::*;
pub use gfunc_ext_pc::*;
pub use int_poly_poly_gen::*;
pub use polygon2d::*;
pub use proj_p_cur::*;
pub use int_conic_conic::*;
