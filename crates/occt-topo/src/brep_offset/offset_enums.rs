//! The three `BRepOffset` enumerations.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKOffset/BRepOffset/`):
//! `BRepOffset_Mode.hxx:27-32`, `BRepOffset_Error.hxx:20-34`,
//! `BRepOffset_Status.hxx:25-31`.

/// `BRepOffset_Mode` (`BRepOffset_Mode.hxx:27-32`): what the offset is applied
/// to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// `BRepOffset_Skin`: offset along the surface of a solid.
    Skin,
    /// `BRepOffset_Pipe`: offset of a curve (a pre-surface).
    Pipe,
    /// `BRepOffset_RectoVerso`: offset of a shell along both sides.
    RectoVerso,
}

/// `BRepOffset_Error` (`BRepOffset_Error.hxx:20-34`): why an offset failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// `BRepOffset_NoError`
    NoError,
    /// `BRepOffset_UnknownError`
    UnknownError,
    /// `BRepOffset_BadNormalsOnGeometry`
    BadNormalsOnGeometry,
    /// `BRepOffset_C0Geometry`
    C0Geometry,
    /// `BRepOffset_NullOffset`
    NullOffset,
    /// `BRepOffset_NotConnectedShell`
    NotConnectedShell,
    /// `BRepOffset_CannotTrimEdges`: exception while trimming edges.
    CannotTrimEdges,
    /// `BRepOffset_CannotFuseVertices`: exception while fusing vertices.
    CannotFuseVertices,
    /// `BRepOffset_CannotExtentEdge`: exception while extending edges.
    CannotExtentEdge,
    /// `BRepOffset_UserBreak`: user break.
    UserBreak,
    /// `BRepOffset_MixedConnectivity`: faces are partially C0 and tangent
    /// along the same edge.
    MixedConnectivity,
}

/// `BRepOffset_Status` (`BRepOffset_Status.hxx:25-31`): the status of one
/// offset face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// `BRepOffset_Good`
    Good,
    /// `BRepOffset_Reversed`: e.g. offset larger than a cylinder radius.
    Reversed,
    /// `BRepOffset_Degenerated`: e.g. offset equal to a cylinder radius.
    Degenerated,
    /// `BRepOffset_Unknown`: e.g. for a Bezier surface.
    Unknown,
}
