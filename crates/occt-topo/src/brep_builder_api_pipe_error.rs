//! `BRepBuilderAPI_PipeError` — the status a pipe shell construction reports.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_PipeError.hxx:21-28`
//! (OCCT 8.0.0). A plain enumerator with no `.cxx`.
//!
//! `BRepOffsetAPI_MakePipeShell::GetStatus` (`BRepOffsetAPI_MakePipeShell.cxx:150-171`)
//! translates `GeomFill_PipeError` into this type, which is how the enum is used
//! by the port: [`PipeError::from_geom_fill_pipe_error`].

use occt_geom::geom_fill::GeomFillPipeError;

/// `BRepBuilderAPI_PipeError` (`BRepBuilderAPI_PipeError.hxx:21-28`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipeError {
    /// `BRepBuilderAPI_PipeDone`.
    Done,
    /// `BRepBuilderAPI_PipeNotDone`.
    NotDone,
    /// `BRepBuilderAPI_PlaneNotIntersectGuide`.
    PlaneNotIntersectGuide,
    /// `BRepBuilderAPI_ImpossibleContact`.
    ImpossibleContact,
}

impl PipeError {
    /// The `switch` of `BRepOffsetAPI_MakePipeShell::GetStatus`
    /// (`BRepOffsetAPI_MakePipeShell.cxx:152-170`): `GeomFill_PipeOk` maps to
    /// `PipeDone`, the two guide/contact failures keep their own values and
    /// every other `GeomFill_PipeError` (including `GeomFill_PipeNotOk`) falls
    /// into the `default:` arm, i.e. `PipeNotDone`.
    pub fn from_geom_fill_pipe_error(stat: GeomFillPipeError) -> Self {
        match stat {
            GeomFillPipeError::PipeOk => PipeError::Done,
            GeomFillPipeError::PlaneNotIntersectGuide => PipeError::PlaneNotIntersectGuide,
            GeomFillPipeError::ImpossibleContact => PipeError::ImpossibleContact,
            _ => PipeError::NotDone,
        }
    }
}
