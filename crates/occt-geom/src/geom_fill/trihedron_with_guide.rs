//! `GeomFill_TrihedronWithGuide` — the abstract law of a trihedron carried by
//! a guide curve, together with the `myGuide` / `myTrimG` /
//! `myCurPointOnGuide` bookkeeping the OCCT base class owns.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_TrihedronWithGuide.hxx:27-44` and
//! `GeomFill_TrihedronWithGuide.cxx:22-33` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Guide()` (`cxx:22-25`) and `Origine` (`hxx:35`, pure virtual) are the
//!   trait methods; `CurrentPointOnGuide` (`cxx:30-33`) is a default method
//!   reading [`TrihedronWithGuideState::cur_point_on_guide`].
//! - `myTrimG` is written by the derived `SetInterval` through a shared
//!   handle, so the state lives behind a `Mutex` and the trait takes `&self`
//!   (the same convention as [`crate::geom_fill::trihedron_law`]).

use std::sync::{Arc, Mutex};

use occt_core::gp::GpPnt;

use crate::curve::Curve;
use crate::geom_fill::trihedron_law::TrihedronLaw;


/// `myGuide` / `myTrimG` / `myCurPointOnGuide` (`hxx:41-43`).
#[derive(Clone, Default)]
pub struct TrihedronWithGuideState {
    /// `myGuide` (`hxx:41`).
    pub guide: Option<Arc<dyn Curve>>,
    /// `myTrimG` (`hxx:42`).
    pub trim_g: Option<Arc<dyn Curve>>,
    /// `myCurPointOnGuide` (`hxx:43`), written by `D0` / `D1` / `D2`.
    pub cur_point_on_guide: GpPnt,
}

/// `GeomFill_TrihedronWithGuide` (`GeomFill_TrihedronWithGuide.hxx:27-44`).
pub trait TrihedronWithGuide: TrihedronLaw {
    /// The shared `myGuide` / `myTrimG` / `myCurPointOnGuide`.
    fn guide_state(&self) -> &Mutex<TrihedronWithGuideState>;

    /// `GeomFill_TrihedronWithGuide::Guide` (`cxx:22-25`).
    fn guide(&self) -> Option<Arc<dyn Curve>> {
        self.guide_state()
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .guide
            .clone()
    }

    /// `myTrimG` (`hxx:42`).
    fn trim_g(&self) -> Option<Arc<dyn Curve>> {
        self.guide_state()
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .trim_g
            .clone()
    }

    /// `GeomFill_TrihedronWithGuide::CurrentPointOnGuide` (`cxx:30-33`): the
    /// point on the guide found by the last `D0` / `D1` / `D2`.
    fn current_point_on_guide(&self) -> GpPnt {
        self.guide_state()
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .cur_point_on_guide
    }

    /// `GeomFill_TrihedronWithGuide::Origine` (`hxx:35`, pure virtual).
    fn origine(&self, param1: f64, param2: f64);

    /// `GeomFill_TrihedronWithGuide::Copy` (`hxx:38`) returning a guide-typed
    /// handle. `GeomFill_LocationGuide::Copy` (`cxx:519-521`) down-casts the
    /// `Copy` result back to `GeomFill_TrihedronWithGuide`; a
    /// `Arc<dyn TrihedronLaw>` cannot be down-cast in Rust, so the guide laws
    /// expose the same construction through this method.
    fn copy_with_guide(&self) -> Arc<dyn TrihedronWithGuide>;
}
