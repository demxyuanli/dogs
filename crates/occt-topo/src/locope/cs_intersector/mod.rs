//! `LocOpe_CSIntersector` (`LocOpe_CSIntersector.hxx:35-151`, `.cxx:55-313`,
//! `.lxx:23-48`).
//!
//! Intersects a *set* of axes, circles or curves with every face of a shape and
//! keeps one sorted `LocOpe_PntFace` sequence per element of the set
//! (`cxx:68-177`). The face-by-face work goes through
//! [`crate::int_curves_face::FaceIntersector`] and the harness of `AddPoints`
//! (`cxx:517-589`); the four `Localize*` scans are the shared ones in
//! [`super::localize`], here with the caller's `Tol`.
//!
//! `myPoints` is a C array of `NCollection_Sequence` (`hxx:150`) whose length is
//! `myNbelem`; the port keeps it as a `Vec<Vec<PntFace>>` indexed `[I - 1]`.

use std::sync::Arc;

use occt_core::gp::{GpCirc, GpLin};
use occt_geom::Curve;

use crate::shape::TopoShape;

use super::localize::{self, Localized};
use super::pnt_face::PntFace;

mod add_points;
mod perform;

pub struct CSIntersector {
    pub(super) my_done: bool,
    pub(super) my_shape: TopoShape,
    pub(super) my_points: Vec<Vec<PntFace>>,
}

impl Default for CSIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl CSIntersector {
    /// `LocOpe_CSIntersector()` (`.lxx:23-29`): `myDone` false, no points.
    pub fn new() -> Self {
        Self {
            my_done: false,
            my_shape: TopoShape::new(crate::abs::ShapeType::Compound),
            my_points: Vec::new(),
        }
    }

    /// `LocOpe_CSIntersector(S)` (`.lxx:33-40`): stores `S` only, no `Perform`.
    pub fn with_shape(s: &TopoShape) -> Self {
        Self {
            my_done: false,
            my_shape: s.clone(),
            my_points: Vec::new(),
        }
    }

    /// `Init(S)` (`.cxx:57-65`): `myDone` false, drops the sequences.
    pub fn init(&mut self, s: &TopoShape) {
        self.my_done = false;
        self.my_shape = s.clone();
        self.my_points.clear();
    }

    /// `Perform(NCollection_Sequence<gp_Lin>)` (`.cxx:69-98`).
    pub fn perform_lines(&mut self, lines: &[GpLin]) -> Result<(), String> {
        perform::perform_lines(self, lines)
    }

    /// `Perform(NCollection_Sequence<gp_Circ>)` (`.cxx:102-137`).
    pub fn perform_circles(&mut self, circles: &[GpCirc]) -> Result<(), String> {
        perform::perform_circles(self, circles)
    }

    /// `Perform(NCollection_Sequence<Handle(Geom_Curve)>)` (`.cxx:140-177`).
    /// A null handle in the sequence is skipped (`cxx:162-165`), so the slice
    /// holds `Option`s.
    pub fn perform_curves(
        &mut self,
        curves: &[Option<Arc<dyn Curve>>],
    ) -> Result<(), String> {
        perform::perform_curves(self, curves)
    }

    /// `IsDone()` (`.lxx:43-46`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `NbPoints(I)` (`.cxx:182-192`), `I` 1-based.
    pub fn nb_points(&self, i: usize) -> Result<usize, String> {
        self.sequence(i).map(|s| s.len())
    }

    /// `Point(I, Index)` (`.cxx:196-208`), both 1-based.
    pub fn point(&self, i: usize, index: usize) -> Result<&PntFace, String> {
        let seq = self.sequence(i)?;
        if index == 0 || index > seq.len() {
            return Err("Standard_OutOfRange".to_string());
        }
        Ok(&seq[index - 1])
    }

    /// `LocalizeAfter(I, From, Tol, Or, IndFrom, IndTo)` (`.cxx:214-241`).
    pub fn localize_after_param(
        &self,
        i: usize,
        from: f64,
        tol: f64,
    ) -> Result<Option<Localized>, String> {
        let seq = self.sequence(i)?;
        Ok(localize::loc_after_param(seq, from, tol))
    }

    /// `LocalizeBefore(I, From, Tol, Or, IndFrom, IndTo)` (`.cxx:245-265`).
    pub fn localize_before_param(
        &self,
        i: usize,
        from: f64,
        tol: f64,
    ) -> Result<Option<Localized>, String> {
        let seq = self.sequence(i)?;
        Ok(localize::loc_before_param(seq, from, tol))
    }

    /// `LocalizeAfter(I, FromInd, Tol, Or, IndFrom, IndTo)` (`.cxx:269-292`).
    pub fn localize_after_index(
        &self,
        i: usize,
        from_ind: usize,
        tol: f64,
    ) -> Result<Option<Localized>, String> {
        let seq = self.sequence(i)?;
        Ok(localize::loc_after_index(seq, from_ind, tol))
    }

    /// `LocalizeBefore(I, FromInd, Tol, Or, IndFrom, IndTo)` (`.cxx:296-316`).
    pub fn localize_before_index(
        &self,
        i: usize,
        from_ind: usize,
        tol: f64,
    ) -> Result<Option<Localized>, String> {
        let seq = self.sequence(i)?;
        Ok(localize::loc_before_index(seq, from_ind, tol))
    }

    /// The `myPoints[I - 1]` slot with the `StdFail_NotDone` /
    /// `Standard_OutOfRange` guards every accessor opens with.
    fn sequence(&self, i: usize) -> Result<&[PntFace], String> {
        if !self.my_done {
            return Err("StdFail_NotDone".to_string());
        }
        if i == 0 || i > self.my_points.len() {
            return Err("Standard_OutOfRange".to_string());
        }
        Ok(&self.my_points[i - 1])
    }
}
