//! `BRepBlend_PointOnRst` — a point of a blending line that sits on a
//! restriction (an arc) of a surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_PointOnRst.hxx:34-75`,
//! `BRepBlend_PointOnRst.lxx:16-35`, `BRepBlend_PointOnRst.cxx:22-47`
//! (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<Adaptor2d_Curve2d>` becomes `Option<Arc<dyn Curve2d>>`: the
//!   empty constructor (`cxx:22-25`) leaves the handle null, and `Arc()` hands
//!   it back without a guard in OCCT, so the port keeps the null case
//!   observable instead of raising,
//! - `IntSurf_Transition` becomes [`IntSurfTransition`], the crate's value
//!   type.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::intsurf_transition::IntSurfTransition;

/// `class BRepBlend_PointOnRst` (`hxx:34-75`).
#[derive(Clone)]
pub struct PointOnRst {
    /// `arc` (`hxx:70`)
    arc: Option<Arc<dyn Curve2d>>,
    /// `traline` (`hxx:71`)
    traline: IntSurfTransition,
    /// `traarc` (`hxx:72`)
    traarc: IntSurfTransition,
    /// `prm` (`hxx:73`)
    prm: f64,
}

impl Default for PointOnRst {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for PointOnRst {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PointOnRst")
            .field("has_arc", &self.arc.is_some())
            .field("traline", &self.traline)
            .field("traarc", &self.traarc)
            .field("prm", &self.prm)
            .finish()
    }
}

impl PointOnRst {
    /// `BRepBlend_PointOnRst()` (`cxx:22-25`): only `prm` is assigned; the arc
    /// stays null and the two transitions keep their default (UNDECIDED).
    pub fn new() -> Self {
        Self {
            arc: None,
            traline: IntSurfTransition::new(),
            traarc: IntSurfTransition::new(),
            prm: 0.0,
        }
    }

    /// `BRepBlend_PointOnRst(A, Param, TLine, TArc)` (`cxx:27-38`).
    pub fn with_arc(
        arc: Arc<dyn Curve2d>,
        param: f64,
        traline: IntSurfTransition,
        traarc: IntSurfTransition,
    ) -> Self {
        Self {
            arc: Some(arc),
            traline,
            traarc,
            prm: param,
        }
    }

    /// `SetArc(A, Param, TLine, TArc)` (`cxx:40-47`).
    pub fn set_arc(
        &mut self,
        arc: Arc<dyn Curve2d>,
        param: f64,
        traline: IntSurfTransition,
        traarc: IntSurfTransition,
    ) {
        self.arc = Some(arc);
        self.prm = param;
        self.traline = traline;
        self.traarc = traarc;
    }

    /// `Arc()` (`lxx:18-21`).
    pub fn arc(&self) -> Option<&Arc<dyn Curve2d>> {
        self.arc.as_ref()
    }

    /// `TransitionOnLine()` (`lxx:23-26`).
    pub fn transition_on_line(&self) -> &IntSurfTransition {
        &self.traline
    }

    /// `TransitionOnArc()` (`lxx:28-31`).
    pub fn transition_on_arc(&self) -> &IntSurfTransition {
        &self.traarc
    }

    /// `ParameterOnArc()` (`lxx:33-35`).
    pub fn parameter_on_arc(&self) -> f64 {
        self.prm
    }
}
