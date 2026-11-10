//! `GeomConvert::ConcatC1`, `GeomConvert::C0BSplineToC1BSplineCurve` and
//! `GeomConvert::C0BSplineToArrayOfC1BSplineCurve`, with the file-static
//! helpers they need.
//!
//! Source: `ModelingData/TKGeomBase/GeomConvert/GeomConvert.cxx` (OCCT 8.0.0):
//! * `GeomConvert_law_evaluator` `:463-492`;
//! * `MultNumandDenom` `:496-583`;
//! * `Pretreatment` `:587-606`;
//! * `NeedToBeTreated` `:612-632`;
//! * `Need2DegRepara` `:638-652`;
//! * `Indexmin` `:658-670`;
//! * `ReorderArrayOfG1Curves` `:676-728`;
//! * `GeomConvert_reparameterise_evaluator` `:734-760`;
//! * `ConcatG1` `:764-997`;
//! * `ConcatC1` 7-argument `:1022-1310` and the 6-argument overload
//!   `:1004-1019`;
//! * `C0BSplineToC1BSplineCurve` `:1313-1338`;
//! * `C0BSplineToArrayOfC1BSplineCurve` 3-argument `:1344-1354` and
//!   4-argument `:1356-1425`.
//!
//! The 2D twin of this file lives in `occt-geom2d::geom2d_convert_concat`. The
//! two `.cxx` files are not identical — the differences that matter are called
//! out inline where the port follows the 3D source:
//! * `C0BSplineToArrayOfC1BSplineCurve` counts the C1 spans from
//!   `FirstUKnotIndex()` (`cxx:1366`), the 2D version from `FirstUKnotIndex()+1`;
//! * the closed flag of that function is `BS->IsClosed() &&
//!   V1.IsParallel(V2, AngularTolerance)` (`cxx:1408-1413`), the 2D version
//!   tests the endpoint distance itself;
//! * `ConcatC1`'s closed branch calls `C.Add(Curve1, local_tolerance(j-1))`
//!   (`cxx:1232-1233`, i.e. `After = false`), the 2D version passes `true`;
//! * `Pretreatment` has no pole-count guard (`cxx:593-605`).
//!
//! Divergences, all in the failure channel only:
//! * `Standard_ConstructionError` / `Standard_DimensionError` and the `throw`
//!   of `Replacer` become `Err`.
//! * The two "uninitialised read" windows of the `.cxx` are bounded:
//!   `Pretreatment` skips a curve with fewer than two poles (OCCT reads
//!   `Weight(2)` past the end), and the `Segment` knot walk tests the stored
//!   multiplicity array bounds before the loop test that OCCT evaluates first
//!   (`cxx:1392`).

use crate::bspline_curve::GeomBSplineCurve;
use crate::convert_bspl::CompCurveToBSplineCurve;
use occt_core::convert::ParameterisationType;
use occt_geom2d::curve::Curve2d;
use crate::curve::Curve;
use crate::geom_lprop::{continuity, continuity_default};
use crate::hermit;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::function_ops::{
    flat_knot_sequence, function_multiply, function_reparameterise, merge_bspline_knots,
};
use occt_core::bspl::knots::reparametrize;
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

mod c0_to_c1;
mod closed_repara;
mod concat_c1;
mod concat_g1;
mod helpers;
mod mult_num_and_denom;

pub use c0_to_c1::*;
pub(crate) use closed_repara::*;
pub use concat_c1::*;
pub use concat_g1::*;
pub(crate) use helpers::*;
pub use mult_num_and_denom::*;
