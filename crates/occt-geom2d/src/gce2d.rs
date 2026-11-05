//! `GCE2d_MakeSegment` (`GCE2d_MakeSegment.cxx`): the 2D segment constructor
//! the `BRepMAT2d_Explorer` gap filler (`BRepMAT2d_Explorer.cxx:585`) and
//! `Bisector_BisecAna`'s collapsed-segment branch (`Bisector_BisecAna.cxx:96`)
//! both call.
//!
//! Only the `(P1, P2)` constructor is ported: it is the sole overload the two
//! consumers use.

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d};

use crate::curve::Curve2d;
use crate::line::Geom2dLine;
use crate::trimmed::Geom2dTrimmedCurve;

/// `GCE2d_MakeSegment(P1, P2).Value()` (`GCE2d_MakeSegment.cxx`): the line
/// through `P1` towards `P2`, trimmed from `0` to `P1.Distance(P2)`.
///
/// `gp_Dir2d(P2 - P1)` raises `Standard_ConstructionError` for a zero vector in
/// OCCT; the port falls back to the default direction (`unwrap_or_default`),
/// the same deviation the `Bisector_BisecAna` caller already documents.
pub fn make_segment(p1: &GpPnt2d, p2: &GpPnt2d) -> Geom2dTrimmedCurve {
    let dir = GpDir2d::new(p2.x() - p1.x(), p2.y() - p1.y()).unwrap_or_default();
    let line = Geom2dLine::new(GpAx2d::new(*p1, dir));
    let basis: Arc<dyn Curve2d> = Arc::from(line.clone_dyn());
    Geom2dTrimmedCurve::new(basis, 0.0, p1.distance(p2))
}
