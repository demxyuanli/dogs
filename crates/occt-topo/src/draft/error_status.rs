//! `Draft_ErrorStatus` (`Draft_ErrorStatus.hxx:20-28`): why `Perform` failed.

/// `enum Draft_ErrorStatus` (`Draft_ErrorStatus.hxx:20-28`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DraftErrorStatus {
    /// `Draft_NoError`
    NoError,
    /// `Draft_FaceRecomputation`: the new surface of a face could not be built.
    FaceRecomputation,
    /// `Draft_EdgeRecomputation`: the new 3D curve or pcurve failed.
    EdgeRecomputation,
    /// `Draft_VertexRecomputation`: the new point or parameter failed.
    VertexRecomputation,
}
