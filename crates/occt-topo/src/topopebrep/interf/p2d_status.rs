//! `TopOpeBRep_P2Dstatus` — status of a `TopOpeBRep_Point2d`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_P2Dstatus.hxx:20-27`.

/// `enum TopOpeBRep_P2Dstatus` (`P2Dstatus.hxx:22-28`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum P2DStatus {
    /// `TopOpeBRep_P2DUNK`
    #[default]
    Unknown,
    /// `TopOpeBRep_P2DINT`
    Intersection,
    /// `TopOpeBRep_P2DSGF` — segment first point
    SegmentFirst,
    /// `TopOpeBRep_P2DSGL` — segment last point
    SegmentLast,
    /// `TopOpeBRep_P2DNEW`
    New,
}
