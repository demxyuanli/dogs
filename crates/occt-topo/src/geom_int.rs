//! GeomInt — intersection curves between two `Geom` surfaces.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomInt/` (OCCT 8.0.0):
//! `GeomInt.cxx`, `GeomInt_IntSS.cxx` / `_1.cxx`, `GeomInt_LineTool.cxx`,
//! `GeomInt_LineConstructor.cxx`, `GeomInt_WLApprox.hxx`,
//! `GeomInt_ParameterAndOrientation.cxx`.
//!
//! `GeomInt_WLApprox` generated ApproxInt templates are not instantiated; the
//! public Perform path falls back to `MakeBSpline` when approximation is not
//! done, matching `GeomInt_IntSS::MakeCurve` when `theapp3d.IsDone()` is false.

#[path = "geom_int_types.rs"]
mod types;
#[path = "geom_int_param_ori.rs"]
mod param_ori;
#[path = "geom_int_topol.rs"]
mod topol;
#[path = "geom_int_quadric.rs"]
mod quadric;
#[path = "geom_int_line_tool.rs"]
mod line_tool;
#[path = "geom_int_line_ctor.rs"]
mod line_ctor;
#[path = "geom_int_line_ctor_gline.rs"]
mod line_ctor_gline;
#[path = "geom_int_line_ctor_rline.rs"]
mod line_ctor_rline;
#[path = "geom_int_intss.rs"]
mod intss;
#[path = "geom_int_intss_make.rs"]
mod intss_make;
#[path = "geom_int_intss_pcurve.rs"]
mod intss_pcurve;
#[path = "geom_int_intss_bspline.rs"]
mod intss_bspline;
#[path = "geom_int_wl_approx.rs"]
mod wl_approx;

pub use types::{ALine, GLine, GLineKind, GeomIntLine, IntPatchIType, RLine};
pub use param_ori::ParameterAndOrientation;
pub use topol::{RestrictionArc, TopolTool};
pub use line_tool::{first_parameter, last_parameter, nb_vertex, vertex};
pub use line_ctor::LineConstructor;
pub use intss::{IntSS, IntSSLine};
pub use intss_bspline::{make_bspline, make_bspline2d};
pub use intss_pcurve::{
    adjust_u_periodic_curve2d, build_pcurves, treat_rline, trim_iline_on_surf_boundaries,
};
pub use wl_approx::{MultiBSpCurve, ParametrizationType, WlApprox};
pub(crate) use quadric::surface_parameters;

use crate::int_tools_wline::adjust_periodic as wline_adjust_periodic;

/// `GeomInt::AdjustPeriodic` (`GeomInt.cxx:21`). Default `theEps` is 0.
pub fn adjust_periodic(
    the_par: f64,
    the_par_min: f64,
    the_par_max: f64,
    the_period: f64,
    the_eps: f64,
) -> (f64, f64) {
    wline_adjust_periodic(the_par, the_par_min, the_par_max, the_period, the_eps)
}
