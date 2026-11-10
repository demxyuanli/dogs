//! `Extrema_ExtPC` / `Extrema_GGExtPC::Perform` + `IntervalPerform`.
//!
//! Source: `Extrema_GGExtPC.hxx` (Perform ~132-528, IntervalPerform ~590-614,
//! AddSol ~617-632), `Extrema_GGenExtPC.hxx` (Perform ~160-173),
//! `Extrema_GFuncExtPC.hxx` (Value/Values/GetStateNumber),
//! `Extrema_CurveTool::DeflCurvIntervals` (`Extrema_CurveTool.cxx:41-91`).
//!
//! Used by `ShapeAnalysis_Curve::ProjectAct` (`ShapeAnalysis_Curve.cxx:277`).
//! Singular DN / MaxDerivOrder: `GFuncExtPC::Value` (`hxx:167-241`) and
//! `Values` three-point DF (`hxx:297-338`).

use occt_core::bspl::knots::unique_knots_mults;
use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::gcpnts::{perform_tangential_curve, CurveSecondDeriv};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::{CONFUSION, PCONFUSION, Precision};

use crate::curve::Curve;

const EPS_DL_DT: f64 = 1.0e-3;
const MAX_DEFL: f64 = 1.0e3;
const MIN_DEFL: f64 = 1.0e-3;
const DEFL_LENGTH_SAMPLES: i32 = 23;

const TOL_F: f64 = 1.0e-10;
const MIN_TOL: f64 = 1.0e-20;
const TOL_FACTOR: f64 = 1.0e-12;
const MIN_STEP: f64 = 1.0e-7;
const MAX_ORDER: i32 = 3;
const MAX_SAMPLE: i32 = 17;

/// One `Extrema_ExtPC` solution (`SquareDistance` / `IsMin` / `Point`).

mod api;
mod defl;
mod function;
mod general;
mod gfunc;
mod perform;
mod state;

pub use api::*;
pub(crate) use defl::*;
pub use function::*;
pub(crate) use general::*;
pub(crate) use perform::*;
pub(crate) use state::*;
