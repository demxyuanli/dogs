//! `GeomFill_PipeError` — the status a sweep law reports.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_PipeError.hxx:20-26` (OCCT 8.0.0).

/// `GeomFill_PipeError` (`GeomFill_PipeError.hxx:20-26`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomFillPipeError {
    /// `GeomFill_PipeOk`.
    PipeOk,
    /// `GeomFill_PipeNotOk`.
    PipeNotOk,
    /// `GeomFill_PlaneNotIntersectGuide`.
    PlaneNotIntersectGuide,
    /// `GeomFill_ImpossibleContact`.
    ImpossibleContact,
}
