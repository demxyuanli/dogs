//! `ShapeConstruct_ProjectCurveOnSurface::getLine`
//! (`ShapeConstruct_ProjectCurveOnSurface.cxx:902-1102`): the first branch of
//! `approxPCurve` (`cxx:1119`), taken before `isAnIsoparametric` (`cxx:1170`)
//! and before the interpolation/approximation arms. It decides whether the
//! sampled projection is a straight line in `(u, v)` within `myPreci` and, if
//! so, returns an exact `Geom2d_Line` (or a 2-pole degree-1 B-spline when the
//! 2D chord length differs from the parameter span).
//!
//! Helpers ported with it:
//! - `fixPeriodicityTroubles` (`cxx:286-412`) with
//!   `ShapeAnalysis::AdjustByPeriod` / `AdjustToPeriod`
//!   (`ShapeAnalysis.cxx:48-69`). `AdjustByPeriod` has a single definition in
//!   [`crate::shhealing::adjust_by_period`]; `AdjustToPeriod` is the thin
//!   wrapper below.
//! - the `myCache` endpoint cache (`cxx:66-70` `CachePoint`, filled at
//!   `cxx:1122-1142` and `cxx:1892-1903`) and the
//!   `SurfaceProjectorWithCache` B-spline corner cache
//!   (`cxx:78-217`), both consumed by `getLine`'s `projectPoint` lambda
//!   (`cxx:944-965`). A hit replaces the global projection with
//!   `ShapeAnalysis_Surface::NextValueOfUV` seeded from the cached point and,
//!   for the first probe (`theIndex == 0`), hands that point to
//!   `fixPeriodicityTroubles` as `theSavedPoint` / `theSavedParam`
//!   (`cxx:957-961`). Without it `theSavedPoint` stays `-1` and the fix anchors
//!   on the middle of the period (`cxx:295-301`), which mis-unwraps a pcurve
//!   that starts on the seam.
//!
//! PARKED (reported, not fudged):
//! - nothing from the `isoParam` chain: `isAnIsoparametric`
//!   (`surface_projector::is_an_isoparametric`, `cxx:2461-2796`), the `isoParam` arms of the
//!   sample loop (`cxx:1317-1375`), the `p1OnIso` / `p2OnIso` endpoint
//!   overrides (`cxx:1379-1385`) and the `if (!isoPar2d3d)` guard of
//!   `projectDegeneratedPoints` (`cxx:1467`) are all ported. A model with an
//!   isoparametric boundary edge therefore takes the iso arm instead of the
//!   `myCache` arm below. The walking-Newton projection of the sample loop
//!   itself (`cxx:1432`) stays ported, in
//!   `surface_projector::project_curve_on_surface_perform`.

use super::surface_projector::{bspline_from_samples, next_value_of_uv, value_of_uv_with_gap};
use crate::shhealing::adjust_by_period;
use occt_core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};

mod cache;
mod interpolate;
mod line_pcurve;
mod out;
mod periodic;

pub use cache::*;
pub(crate) use interpolate::*;
pub(crate) use line_pcurve::*;
pub use out::*;
pub use periodic::*;
