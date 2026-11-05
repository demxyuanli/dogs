//! The three status enumerations of `FilletSurf`.
//!
//! Sources (OCCT 8.0.0):
//! - `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_StatusDone.hxx:22-28`
//! - `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_ErrorTypeStatus.hxx:22-30`
//! - `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_StatusType.hxx:22-28`
//!
//! The three plain `enum`s are kept one per OCCT header, in the OCCT
//! declaration order, because the API answers them by name.

/// `FilletSurf_StatusDone` (`FilletSurf_StatusDone.hxx:22-28`): the status of
/// the computation of the fillet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusDone {
    /// `FilletSurf_IsOk`: no problem during the computation.
    IsOk,
    /// `FilletSurf_IsNotOk`: no result is produced.
    IsNotOk,
    /// `FilletSurf_IsPartial`: the result is partial.
    IsPartial,
}

/// `FilletSurf_ErrorTypeStatus` (`FilletSurf_ErrorTypeStatus.hxx:22-30`): the
/// error kind when [`StatusDone::IsNotOk`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorTypeStatus {
    /// `FilletSurf_EmptyList`: the edge list is empty.
    EmptyList,
    /// `FilletSurf_EdgeNotG1`: the edges are not G1.
    EdgeNotG1,
    /// `FilletSurf_FacesNotG1`: two connected faces on the same support are
    /// not G1.
    FacesNotG1,
    /// `FilletSurf_EdgeNotOnShape`: the edge is not on the shape.
    EdgeNotOnShape,
    /// `FilletSurf_NotSharpEdge`: the edge is not sharp.
    NotSharpEdge,
    /// `FilletSurf_PbFilletCompute`: problem during the computation of the
    /// fillet.
    PbFilletCompute,
}

/// `FilletSurf_StatusType` (`FilletSurf_StatusType.hxx:22-28`): how many
/// extremities of a section landed on the edge of the corresponding support
/// face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusType {
    /// `FilletSurf_TwoExtremityOnEdge`.
    TwoExtremityOnEdge,
    /// `FilletSurf_OneExtremityOnEdge`.
    OneExtremityOnEdge,
    /// `FilletSurf_NoExtremityOnEdge`.
    NoExtremityOnEdge,
}
