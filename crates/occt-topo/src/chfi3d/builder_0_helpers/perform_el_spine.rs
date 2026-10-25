//! `ChFi3d_PerformElSpine` — **UNPORTED**.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4797-5352`
//! (OCCT 8.0.0).
//!
//! The function rebuilds one `ChFiDS_ElSpine`'s 3D curve from the edges of the
//! spine: it trims the first / last edge to the `ElSpine` range, converts every
//! piece to a B-spline, smooths non-`C2` pieces with `ChFi3d_ApproxByC2`,
//! concatenates them with `GeomConvert_CompCurveToBSplineCurve`, reparametrizes
//! the knots onto `[Wrefdeb, Wreffin]` and then builds the tangent extensions at
//! both extremities.
//!
//! Blockers that keep it out of this batch (each is a `.cxx` unit of its own):
//! - `ChFi3d_ApproxByC2` (`ChFi3d_Builder_0.cxx:4655-4730`) and `CurveCleaner`
//!   (`:4700-4745` region) are not ported,
//! - `ChFi3d_IsSmooth` is ported ([`crate::chfi3d::builder_helpers::is_smooth`]),
//!   but the two above are what the body calls.
//!
//! `PerformSetOfKPart` calls it once per `ElSpine` of the spine
//! (`ChFi3d_Builder_2.cxx:3266-3278`), and once more per offset `ElSpine` when
//! `Spine->Mode() == ChFiDS_ConstThroatWithPenetrationChamfer`. Since the KPart
//! data itself is still [`crate::chfi3d::chfi_k_part`]'s inert placeholder, no
//! `ElSpine` reaches this function on the current baseline; the `Err` is there
//! so that a later caller cannot silently keep an `ElSpine` without its curve.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;

use crate::chfi_ds::{ElSpine, SpineHandle};

/// `ChFi3d_PerformElSpine(HES, Spine, Continuity, Tol, IsOffset = false)`
/// (`cxx:4797-4802`) — **UNPORTED**, answers `Err`.
pub fn ch_fi3d_perform_el_spine(
    _hes: &mut Arc<ElSpine>,
    _spine: &SpineHandle,
    _continuity: Shape,
    _tol: f64,
    _is_offset: bool,
) -> Result<(), String> {
    Err("ChFi3d_PerformElSpine: UNPORTED (ChFi3d_Builder_0.cxx:4797-5352)".into())
}
