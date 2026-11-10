//! `Approx_SameParameter`.
//! Source: `Approx_SameParameter.cxx` (`Build` 318-539, `BuildInitialDistribution`
//! 547-579, `IncreaseInitialNbSamples` 588-648, `CheckSameParameter` 653-766,
//! `ComputeTangents` 770-809, `Interpolate` 813-853, `Check` 192-266,
//! Evaluator 64-102, `ComputeTolReached` 153-188, `IncreaseNbPoles` 857-991).
//!
//! `CheckSameParameter` uses `Extrema_LocateExtPC` then `Extrema_ExtPC`
//! whole-space search (`cxx:700-758`).

use std::sync::Arc;

use occt_core::bspl::banded_interp::{interpolate_contact, knot_sequence};
use occt_core::bspl::eval::{eval_curve, eval_curve_d1};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};

use crate::adv_approx::ApproxAFunction1dPair;
use crate::curve::Curve;
use crate::surface::Surface;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

const NB_SAMPLES: usize = 22;
const MAX_ARRAY: usize = 1000;

/// Result of `Approx_SameParameter`.
pub struct ApproxSameParameter {
    pub done: bool,
    pub same_parameter: bool,
    pub tol_reached: f64,
    /// Rebuilt 2d curve when `IsDone && !IsSameParameter` (`cxx:1648-1654`).
    pub curve2d: Option<Arc<dyn Curve2d>>,
}

mod approx;
mod checks;
mod cons;
mod distance;
mod perform;
mod resolution;

pub(crate) use approx::*;
pub(crate) use checks::*;
pub(crate) use cons::*;
pub(crate) use distance::*;
pub use resolution::*;
