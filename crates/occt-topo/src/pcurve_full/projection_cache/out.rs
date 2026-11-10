//! `GetLineOut`: the two out-parameters of `getLine`.

use super::super::prelude::*;

/// The two out-parameters of `getLine` (`cxx:905-907` `theIsRecompute` /
/// `theIsFromCache`) plus the endpoint projections it writes into `thePoints2d`
/// before it can fail (`cxx:1020-1021`), which the caller's sample loop reuses
/// verbatim when `!isRecompute` (`cxx:1390-1399`).
#[derive(Clone, Copy)]
pub struct GetLineOut {
    pub is_recompute: bool,
    pub is_from_cache: bool,
    pub first2d: GpPnt2d,
    pub last2d: GpPnt2d,
    /// `theIndex == 0` probe found in the caches, i.e. `theSavedPoint >= 0`.
    pub saved_point: GpPnt2d,
    /// The four probe projections, in `aP[0..4]` order, after
    /// `fixPeriodicityTroubles`.
    pub probes2d: [GpPnt2d; 4],
    /// `SurfaceProjectorWithCache::Gap()` left behind by the last probe
    /// projected, i.e. `mySurf->Gap()` after `getLine` returns. It is the
    /// `gap = mySurf->Gap();` seed the caller's sample loop reads at
    /// `ShapeConstruct_ProjectCurveOnSurface.cxx:1393`.
    pub last_gap: f64,
}

impl Default for GetLineOut {
    fn default() -> Self {
        Self {
            is_recompute: false,
            is_from_cache: false,
            first2d: GpPnt2d::zero(),
            last2d: GpPnt2d::zero(),
            saved_point: GpPnt2d::zero(),
            probes2d: [GpPnt2d::zero(); 4],
            last_gap: 0.0,
        }
    }
}
