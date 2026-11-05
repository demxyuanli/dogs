//! `FairCurve_AnalysisCode` (`FairCurve_AnalysisCode.hxx:34-40`).

/// `FairCurve_AnalysisCode` (`FairCurve_AnalysisCode.hxx:34-40`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FairCurveAnalysisCode {
    /// `FairCurve_OK`: the computation is successfully completed.
    Ok,
    /// `FairCurve_NotConverged`: the algorithm does not converge.
    NotConverged,
    /// `FairCurve_InfiniteSliding`: the sliding is infinite, the computation
    /// stops.
    InfiniteSliding,
    /// `FairCurve_NullHeight`: no matter is left at one of the ends of the
    /// curve, the computation stops.
    NullHeight,
}
