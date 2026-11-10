//! Shared BRep_Tool tolerance clamp.

/// The four lines every `BRep_Tool::Tolerance` overload runs on the value it
/// reads off the `TShape` (`BRep_Tool.cxx:140-149`, `:884-893`, `:1323-1332`):
///
/// ```text
/// constexpr double pMin = Precision::Confusion();
/// if (p > pMin) return p; else return pMin;
/// ```
///
/// Note this is `BRep_Tool` behaviour only: the underlying accessors
/// (`BRep_TVertex::Tolerance`, `BRep_TEdge::Tolerance`, `BRep_TFace::Tolerance`)
/// keep answering the stored `RealEpsilon()`, which is what the persistence
/// layer (`BinTools_ShapeSet`, `XmlMDF`) writes out.
pub(super) fn clamp_to_confusion(p: f64) -> f64 {
    let p_min = occt_core::precision::CONFUSION;
    if p > p_min {
        p
    } else {
        p_min
    }
}
