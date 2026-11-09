//! `LocOpe_RevolutionForm` (`LocOpe_RevolutionForm.hxx:25-90`,
//! `.cxx:35-...`).
//!
//! The feature front end of a revolved prism. `IntPerf`
//! (`LocOpe_RevolutionForm.cxx:72-191`) is a verbatim copy of
//! `LocOpe_Revol::IntPerf` (`LocOpe_Revol.cxx:94-212`) without the two curve
//! queries. The port therefore drives one [`super::revol::Revol`] body instead
//! of duplicating that control flow, and keeps the two points the reference
//! stores but never reads (`myPnt1`, `myPnt2`).
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | `Perform` (`cxx:47-65`) | this file |
//! | `IntPerf` (`cxx:72-191`) | [`super::revol::Revol::perform`] / [`super::revol::Revol::perform_with_dec`] |
//! | `Shape` / `FirstShape` / `LastShape` / `Shapes` | accessors |
//!
//! Deviations: `myRes` is `Option<TopoShape>` and `Shape()` answers `None`
//! where OCCT raises `StdFail_NotDone`; `Shapes(S)` answers `None` where the
//! reference's `myMap(S)` raises on an unbound key.

use occt_core::gp::{GpAx1, GpPnt};

use crate::locope::revol::Revol;
use crate::shape::TopoShape;

/// `LocOpe_RevolutionForm`.
#[derive(Debug, Clone, Default)]
pub struct RevolutionForm {
    /// The stored points (`hxx`): written by `Perform`, never read.
    my_pnt1: GpPnt,
    /// See [`Self::my_pnt1`].
    my_pnt2: GpPnt,
    /// The shared `IntPerf` body (`LocOpe_Revol.cxx:94-212`).
    revol: Revol,
    /// `bool myDone`
    my_done: bool,
}

impl RevolutionForm {
    /// `LocOpe_RevolutionForm()`.
    pub fn new() -> Self {
        Self {
            my_pnt1: GpPnt::new(0.0, 0.0, 0.0),
            my_pnt2: GpPnt::new(0.0, 0.0, 0.0),
            ..Default::default()
        }
    }

    /// `Perform(Base, Axis, Angle)` (`cxx:47-65`).
    pub fn perform(&mut self, base: &TopoShape, axis: &GpAx1, angle: f64) {
        self.revol.perform(base, axis, angle);
        self.my_done = self.revol.is_done();
    }

    /// `Perform(Base, Axis, Angle, angledec)`.
    pub fn perform_with_dec(&mut self, base: &TopoShape, axis: &GpAx1, angle: f64, dec: f64) {
        self.revol.perform_with_dec(base, axis, angle, dec);
        self.my_done = self.revol.is_done();
    }

    /// `Shape()`; `None` is the `StdFail_NotDone` raise.
    pub fn shape(&self) -> Option<&TopoShape> {
        if self.my_done {
            self.revol.shape()
        } else {
            None
        }
    }

    /// `FirstShape()`.
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.revol.first_shape()
    }

    /// `LastShape()`.
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.revol.last_shape()
    }

    /// `Shapes(S)`.
    pub fn shapes(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        self.revol.shapes(s)
    }

    /// The two stored points; kept so the port does not drop reference state.
    pub fn points(&self) -> (&GpPnt, &GpPnt) {
        (&self.my_pnt1, &self.my_pnt2)
    }
}
