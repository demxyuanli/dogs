//! `ShapeConstruct_ProjectCurveOnSurface::getLine`
//! (`ShapeConstruct_ProjectCurveOnSurface.cxx:902-1102`): the first branch of
//! `approxPCurve` (`cxx:1119`), taken before `isAnIsoparametric` (`cxx:1170`)
//! and before the interpolation/approximation arms. It decides whether the
//! sampled projection is a straight line in `(u, v)` within `myPreci` and, if
//! so, returns an exact `Geom2d_Line` (or a 2-pole degree-1 B-spline when the
//! 2D chord length differs from the parameter span).
//!
//! Helpers ported with it:
//! - `fixPeriodicityTroubles` (`cxx:286-412`) with
//!   `ShapeAnalysis::AdjustByPeriod` / `AdjustToPeriod`
//!   (`ShapeAnalysis.cxx:48-69`). `AdjustByPeriod` has a single definition in
//!   [`crate::shhealing::adjust_by_period`]; `AdjustToPeriod` is the thin
//!   wrapper below.
//! - the `myCache` endpoint cache (`cxx:66-70` `CachePoint`, filled at
//!   `cxx:1122-1142` and `cxx:1892-1903`) and the
//!   `SurfaceProjectorWithCache` B-spline corner cache
//!   (`cxx:78-217`), both consumed by `getLine`'s `projectPoint` lambda
//!   (`cxx:944-965`). A hit replaces the global projection with
//!   `ShapeAnalysis_Surface::NextValueOfUV` seeded from the cached point and,
//!   for the first probe (`theIndex == 0`), hands that point to
//!   `fixPeriodicityTroubles` as `theSavedPoint` / `theSavedParam`
//!   (`cxx:957-961`). Without it `theSavedPoint` stays `-1` and the fix anchors
//!   on the middle of the period (`cxx:295-301`), which mis-unwraps a pcurve
//!   that starts on the seam.
//!
//! PARKED (reported, not fudged):
//! - nothing from the `isoParam` chain: `isAnIsoparametric`
//!   (`p01::is_an_isoparametric`, `cxx:2461-2796`), the `isoParam` arms of the
//!   sample loop (`cxx:1317-1375`), the `p1OnIso` / `p2OnIso` endpoint
//!   overrides (`cxx:1379-1385`) and the `if (!isoPar2d3d)` guard of
//!   `projectDegeneratedPoints` (`cxx:1467`) are all ported. A model with an
//!   isoparametric boundary edge therefore takes the iso arm instead of the
//!   `myCache` arm below. The walking-Newton projection of the sample loop
//!   itself (`cxx:1432`) stays ported, in
//!   `p01::project_curve_on_surface_perform`.

use super::prelude::*;
use super::p01::{bspline_from_samples, next_value_of_uv, value_of_uv_with_gap};
use crate::shhealing::adjust_by_period;
use occt_core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};

// ---------------------------------------------------------------------------
// `ShapeConvert_ProjectCurveOnSurface` projection caches
// ---------------------------------------------------------------------------

/// `ShapeConstruct_ProjectCurveOnSurface::CachePoint` (`cxx:66-70`): a 3D
/// sample of the pcurve just built and its 2D projection.
#[derive(Clone, Copy)]
pub struct CachePoint {
    pub first: GpPnt,
    pub second: GpPnt2d,
}

/// `ShapeConstruct_ProjectCurveOnSurface::myCache` (`cxx:66`): the two endpoint
/// anchors of the previously projected pcurve on the same surface, so the next
/// pcurve can be projected along the same parameter branch. Cleared by
/// `SetSurface` (`cxx:676-683`), i.e. one cache lives as long as the
/// `ShapeFix_Edge` keeps its surface.
#[derive(Clone, Default)]
pub struct ProjectorCache {
    /// `myCache(0)` first, `myCache(1)` second; empty before the first pcurve.
    entries: Vec<CachePoint>,
}

impl ProjectorCache {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// `myCache.Length()`.
    fn len(&self) -> usize {
        self.entries.len()
    }

    fn get(&self, j: usize) -> CachePoint {
        self.entries[j]
    }

    /// `cxx:1122-1141` / `cxx:1892-1903`: store the two endpoints of the pcurve
    /// just built. `aChangeCycle` (`cxx:1123-1129`, `cxx:1298-1305`) is set when
    /// the stored anchor `myCache(0)` is closer to the LAST point of this curve
    /// than to its first one, i.e. the curve was sampled in the opposite
    /// direction, and then the roles of `myCache(0)` / `myCache(1)` swap.
    pub(super) fn store(&mut self, points: &[GpPnt], points2d: &[GpPnt2d], change_cycle: bool) {
        let n = points.len();
        if n == 0 || points2d.len() < n {
            return;
        }
        let (a, b) = if change_cycle {
            (
                (points[0], points2d[0]),
                (points[n - 1], points2d[n - 1]),
            )
        } else {
            (
                (points[n - 1], points2d[n - 1]),
                (points[0], points2d[0]),
            )
        };
        self.entries = vec![
            CachePoint {
                first: a.0,
                second: a.1,
            },
            CachePoint {
                first: b.0,
                second: b.1,
            },
        ];
    }

    /// `cxx:1123-1129` / `cxx:1298-1305`: is the incoming cache anchored on the
    /// far end of the curve about to be projected?
    pub(super) fn change_cycle(&self, first: &GpPnt, last: &GpPnt) -> bool {
        if self.entries.is_empty() {
            return false;
        }
        let anchor = self.entries[0].first;
        anchor.distance(first) > anchor.distance(last) && anchor.distance(last) < CONFUSION
    }

    /// `approxPCurve` (`cxx:1403-1421`): the `isRecompute` endpoint arm's
    /// `myCache` scan, first entry whose `SquareDistance` is below `tol2`. The
    /// C++ breaks on the first hit, so a linear scan in cache order is exact.
    pub(super) fn find(&self, p: &GpPnt, tol2: f64) -> Option<CachePoint> {
        self.entries
            .iter()
            .find(|e| e.first.square_distance(p) < tol2)
            .copied()
    }
}

/// `SurfaceProjectorWithCache` (`cxx:78-217`): the B-spline corner cache. When
/// the underlying surface is a clamped `Geom_BSplineSurface`, its four corner
/// poles are known `(u, v)` points, so a sample that lands on a corner is
/// projected from that exact UV with `NextValueOfUV` instead of the global
/// search (`cxx:106-138`).
#[derive(Default)]
struct SurfaceProjectorWithCache {
    corners3d: Vec<GpPnt>,
    corners2d: Vec<GpPnt2d>,
}

impl SurfaceProjectorWithCache {
    /// `SurfaceProjectorWithCache::SurfaceProjectorWithCache` (`cxx:84-98`) +
    /// `buildCornerCache` (`cxx:148-200`). `down_cast<Geom_BSplineSurface>`
    /// (`cxx:92-93`) is a dynamic type test, so a `Geom_RectangularTrimmedSurface`
    /// wrapper gets no cache.
    fn new(s: &dyn Surface) -> Self {
        let mut cache = Self::default();
        if !s.is_bspline_surface() {
            return cache;
        }
        let Some(bs) = s.osculating_bspline() else {
            return cache;
        };
        let nb_u = bs.nb_poles_u();
        let nb_v = bs.nb_poles_v();
        if nb_u == 0 || nb_v == 0 {
            return cache;
        }
        let (u_mults, v_mults) = (
            occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(&bs.knots_u),
            occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(&bs.knots_v),
        );
        if u_mults.1.is_empty() || v_mults.1.is_empty() {
            return cache;
        }
        // `cxx:156-170`.
        let is_u_first_clamped = u_mults.1[0] >= bs.deg_u as i32 + 1;
        let is_u_last_clamped = u_mults.1[u_mults.1.len() - 1] >= bs.deg_u as i32 + 1;
        let is_v_first_clamped = v_mults.1[0] >= bs.deg_v as i32 + 1;
        let is_v_last_clamped = v_mults.1[v_mults.1.len() - 1] >= bs.deg_v as i32 + 1;
        let u_first = bs.knots_u[0];
        let u_last = bs.knots_u[bs.knots_u.len() - 1];
        let v_first = bs.knots_v[0];
        let v_last = bs.knots_v[bs.knots_v.len() - 1];
        // `cxx:172-199`: the four corners, in the order OCCT appends them.
        if is_u_first_clamped && is_v_first_clamped {
            cache.corners3d.push(bs.poles[0][0]);
            cache.corners2d.push(GpPnt2d::new(u_first, v_first));
        }
        if is_u_last_clamped && is_v_first_clamped {
            cache.corners3d.push(bs.poles[nb_u - 1][0]);
            cache.corners2d.push(GpPnt2d::new(u_last, v_first));
        }
        if is_u_first_clamped && is_v_last_clamped {
            cache.corners3d.push(bs.poles[0][nb_v - 1]);
            cache.corners2d.push(GpPnt2d::new(u_first, v_last));
        }
        if is_u_last_clamped && is_v_last_clamped {
            cache.corners3d.push(bs.poles[nb_u - 1][nb_v - 1]);
            cache.corners2d.push(GpPnt2d::new(u_last, v_last));
        }
        cache
    }

    /// `findInCornerCache` (`cxx:205-216`).
    fn find(&self, p: &GpPnt, tol2: f64) -> Option<GpPnt2d> {
        for (c3, c2) in self.corners3d.iter().zip(self.corners2d.iter()) {
            if c3.square_distance(p) < tol2 {
                return Some(*c2);
            }
        }
        None
    }

    /// `SurfaceProjectorWithCache::ValueOfUV` (`cxx:106-115`).
    fn value_of_uv(&self, s: &dyn Surface, p: &GpPnt, tol: f64, tol_sq: f64) -> (GpPnt2d, f64) {
        if let Some(uv) = self.find(p, tol_sq) {
            return next_value_of_uv(s, &uv, p, tol, tol);
        }
        value_of_uv_with_gap(s, p, tol)
    }

    /// `SurfaceProjectorWithCache::NextValueOfUV` (`cxx:125-138`).
    fn next_value_of_uv(
        &self,
        s: &dyn Surface,
        hint: &GpPnt2d,
        p: &GpPnt,
        tol: f64,
        tol_sq: f64,
        step: f64,
    ) -> (GpPnt2d, f64) {
        if let Some(uv) = self.find(p, tol_sq) {
            return next_value_of_uv(s, &uv, p, tol, tol);
        }
        next_value_of_uv(s, hint, p, tol, step)
    }
}

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

/// `ShapeAnalysis::AdjustToPeriod` (`ShapeAnalysis.cxx:66-69`).
pub(super) fn adjust_to_period(val: f64, val_min: f64, val_max: f64) -> f64 {
    adjust_by_period(val, 0.5 * (val_min + val_max), val_max - val_min)
}

/// Coordinate accessor for `fixPeriodicityTroubles`: OCCT's `theIdx` is 1 for
/// U (the `x` coordinate) and 2 for V (`cxx:285`).
fn iso_coord(p: &GpPnt2d, idx: usize) -> f64 {
    if idx == 0 {
        p.x()
    } else {
        p.y()
    }
}

fn set_iso_coord(p: &mut GpPnt2d, idx: usize, v: f64) {
    if idx == 0 {
        p.set_x(v);
    } else {
        p.set_y(v);
    }
}

/// `fixPeriodicityTroubles` (`cxx:286-412`): fold the four probe points onto one
/// period window around `theSavedParam` and undo a period jump between them.
/// Returns `true` when a jump was found and fixed (`theIsRecompute`).
fn fix_periodicity_troubles(
    pnt: &mut [GpPnt2d],
    idx: usize,
    period: f64,
    saved_point: i32,
    saved_param: f64,
) -> bool {
    let mut a_min_param = 0.0;
    let mut a_max_param = period;
    // `cxx:295-315`: with no cached point the reference is the middle of the
    // period and the walk starts at index 0.
    let (a_saved_param, a_saved_point) = if saved_point < 0 {
        (0.5 * period, 0usize)
    } else {
        while a_min_param > saved_param {
            a_min_param -= period;
            a_max_param -= period;
        }
        while a_max_param < saved_param {
            a_min_param += period;
            a_max_param += period;
        }
        (saved_param, saved_point as usize)
    };

    // `cxx:317-324`: an iso line sits on the window border and is pinned there.
    let mut a_fix_iso_param = a_min_param;
    let mut is_iso_line = false;
    if a_max_param - a_saved_param < PCONFUSION || a_saved_param - a_min_param < PCONFUSION {
        a_fix_iso_param = a_saved_param;
        is_iso_line = true;
    }

    for i in 0..4 {
        let mut a_param = iso_coord(&pnt[i], idx);
        a_param += adjust_to_period(a_param, a_min_param, a_max_param);
        if is_iso_line {
            if a_max_param - a_param < PCONFUSION || a_param - a_min_param < PCONFUSION {
                a_param = a_fix_iso_param;
            }
        } else {
            if a_max_param - a_param < PCONFUSION {
                a_param = a_max_param;
            }
            if a_param - a_min_param < PCONFUSION {
                a_param = a_min_param;
            }
        }
        set_iso_coord(&mut pnt[i], idx, a_param);
    }

    // `cxx:355-371`: a monotonicity break is a period jump.
    let mut is_jump = false;
    let mut a_prev_diff = 0.0;
    let mut a_sum_diff = 1.0;
    for i in 0..3 {
        let a_diff = iso_coord(&pnt[i + 1], idx) - iso_coord(&pnt[i], idx);
        if a_diff < -PCONFUSION {
            a_sum_diff *= -1.0;
        }
        if a_diff * a_prev_diff < -PCONFUSION {
            is_jump = true;
        }
        a_prev_diff = a_diff;
    }
    if !is_jump {
        return false;
    }

    let sp = a_saved_point.min(3);
    if a_sum_diff > 0.0 {
        for i in (1..=sp).rev() {
            if iso_coord(&pnt[i], idx) > iso_coord(&pnt[i - 1], idx) {
                let v = iso_coord(&pnt[i - 1], idx) + period;
                set_iso_coord(&mut pnt[i - 1], idx, v);
            }
        }
        for i in sp..3 {
            if iso_coord(&pnt[i], idx) < iso_coord(&pnt[i + 1], idx) {
                let v = iso_coord(&pnt[i + 1], idx) - period;
                set_iso_coord(&mut pnt[i + 1], idx, v);
            }
        }
    } else {
        for i in (1..=sp).rev() {
            if iso_coord(&pnt[i], idx) < iso_coord(&pnt[i - 1], idx) {
                let v = iso_coord(&pnt[i - 1], idx) - period;
                set_iso_coord(&mut pnt[i - 1], idx, v);
            }
        }
        for i in sp..3 {
            if iso_coord(&pnt[i], idx) > iso_coord(&pnt[i + 1], idx) {
                let v = iso_coord(&pnt[i + 1], idx) + period;
                set_iso_coord(&mut pnt[i + 1], idx, v);
            }
        }
    }
    true
}

/// `MaxKnotMult(internal mults)` over the flat knot vector: the multiplicity of
/// the most repeated *interior* knot (`BSplCLib::MaxKnotMult`,
/// `Geom_BSplineCurve_1.cxx:54-55`).
fn max_internal_mult(flat: &[f64]) -> i32 {
    let (knots, mults) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(flat);
    if knots.len() <= 2 {
        return 0;
    }
    mults[1..knots.len() - 1].iter().copied().max().unwrap_or(0)
}

/// `Geom_Surface::IsCNu(1) && Geom_Surface::IsCNv(1)` (`cxx:1035`). The
/// `Geom_Surface` base answers `Standard_True`; `Geom_BSplineSurface` forwards
/// to `Geom_BSplineCurve::IsCN(1)` (`Geom_BSplineCurve_1.cxx:34-59`), true iff
/// the smoothness is at least `C1`, i.e. `degree - MaxKnotMult(internal) >= 1`.
/// A `Geom_BezierSurface` is polynomial on its span, so `N = 1` holds.
fn is_cn_1(s: &dyn Surface) -> bool {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return is_cn_1(basis.as_ref());
    }
    if let Some(basis) = s.offset_basis_surface() {
        return is_cn_1(basis.as_ref());
    }
    if !s.is_bspline_surface() {
        return true;
    }
    let Some(b) = s.osculating_bspline() else {
        return true;
    };
    let cn = |flat: &[f64], deg: usize| (deg as i32 - max_internal_mult(flat)) >= 1;
    cn(&b.knots_u, b.deg_u) && cn(&b.knots_v, b.deg_v)
}

/// The projection of `p` (with its residue, OCCT's `SurfaceProjector::Gap`).
/// `SurfaceProjectorWithCache::ValueOfUV` (`cxx:106-115`) is the corner-cache
/// lookup followed by `mySurf->ValueOfUV(thePoint, theTol)`, i.e.
/// [`value_of_uv`] with the working tolerance as `preci`, and `Gap()` (`cxx:141`)
/// is `mySurf->Gap()`; see [`super::p01::value_of_uv_with_gap`] for how the
/// residue is recovered.
///
/// `ShapeConstruct_ProjectCurveOnSurface::getLine` (`cxx:902-1102`): the
/// straight-pcurve shortcut. `points` / `params` are `thePoints` / `theParams`
/// from `generateCurvePoints` (`cxx:750-752`), `tol` is `myPreci`
/// (`cxx:1119` passes it through), `cache` is `myCache` and `out` receives
/// `theIsRecompute` / `theIsFromCache` and the endpoint projections.
/// `None` means "not a straight pcurve" and the caller falls through to
/// `isAnIsoparametric` / the approximation arms.
pub(super) fn get_line(
    s: &dyn Surface,
    points: &[GpPnt],
    params: &[f64],
    tol: f64,
    cache: &ProjectorCache,
    out: &mut GetLineOut,
) -> Option<Arc<dyn Curve2d>> {
    let a_nb = points.len();
    if a_nb < 2 || params.len() != a_nb {
        return None;
    }
    // `cxx:911-916`: probes 1, 2, nb-1 and nb.
    let a_p = [points[0], points[1], points[a_nb - 2], points[a_nb - 1]];

    let mut a_tol2 = tol * tol;
    let mut a_tol_working = tol;
    let is_periodic_u = s.is_u_periodic();
    let is_periodic_v = s.is_v_periodic();

    // `cxx:925-936`: protection against bad tolerance shapes.
    if a_tol2 > 1.0 {
        a_tol_working = CONFUSION;
        a_tol2 = a_tol_working * a_tol_working;
    }
    if a_tol2 < SQUARE_CONFUSION {
        a_tol2 = SQUARE_CONFUSION;
    }
    let an_old_tol2 = a_tol2;

    let mut a_p2d = [GpPnt2d::zero(); 4];

    // `cxx:941`: `SurfaceProjectorWithCache aProjector(mySurf)`.
    let projector = SurfaceProjectorWithCache::new(s);

    // `cxx:944-965` `projectPoint`: an existing endpoint cache hit first
    // (`cxx:946-961`), then the full projection (`cxx:964`).
    let mut saved_point_num: i32 = -1;
    let mut saved_point = GpPnt2d::zero();
    for i in [0usize, 3] {
        let hit = (0..cache.len()).find(|&j| cache.get(j).first.square_distance(&a_p[i]) < a_tol2);
        let a_gap;
        match hit {
            Some(j) => {
                let a_cache_pnt = cache.get(j);
                let (uv, g) = projector.next_value_of_uv(
                    s,
                    &a_cache_pnt.second,
                    &a_p[i],
                    a_tol_working,
                    a_tol2,
                    a_tol_working,
                );
                a_p2d[i] = uv;
                a_gap = g;
                saved_point_num = i as i32;
                saved_point = a_cache_pnt.second;
                if i == 0 {
                    out.is_from_cache = true;
                }
            }
            None => {
                let (uv, g) = projector.value_of_uv(s, &a_p[i], a_tol_working, a_tol2);
                a_p2d[i] = uv;
                a_gap = g;
            }
        }
        // `cxx:972-978`: `const double aDist = aProjector.Gap();` then
        // `aCurDist = aDist * aDist` grows the running tolerance. `Gap()` is
        // `mySurf->Gap()` (`ShapeConstruct_ProjectCurveOnSurface.cxx:141`),
        // which `ValueOfUV` closes with `if (myGap <= 0) myGap =
        // P3D.Distance(SurfAdapt.Value(S, T));` (`ShapeAnalysis_Surface.cxx:1502-1513`),
        // so `aDist` is the projection residue and `aTol2` only grows for a
        // poor projection or the conic-extrusion arm (`cxx:1312`), which leaves
        // `myGap = -1.` (`cxx:1249`). `getLine` also leaves the value in
        // `SurfaceProjectorWithCache::Gap()` for the caller's sample loop.
        out.last_gap = a_gap;
        let a_cur_dist = a_gap * a_gap;
        if a_tol2 < a_cur_dist {
            a_tol2 = a_cur_dist;
        }
    }
    // `cxx:1020-1021`: `thePoints2d` is written before any failure return.
    out.first2d = a_p2d[0];
    out.last2d = a_p2d[3];

    // `cxx:982-1013`: on a periodic surface the second and last-but-one points
    // seed the period-jump fix, which is anchored on the cached point when the
    // first probe came from a cache (`cxx:957-961`).
    let mut is_recompute = false;
    if is_periodic_u || is_periodic_v {
        for i in [1usize, 2] {
            let hit =
                (0..cache.len()).find(|&j| cache.get(j).first.square_distance(&a_p[i]) < a_tol2);
            let a_gap;
            match hit {
                Some(j) => {
                    let a_cache_pnt = cache.get(j);
                    let (uv, g) = projector.next_value_of_uv(
                        s,
                        &a_cache_pnt.second,
                        &a_p[i],
                        a_tol_working,
                        a_tol2,
                        a_tol_working,
                    );
                    a_p2d[i] = uv;
                    a_gap = g;
                    saved_point_num = i as i32;
                    saved_point = a_cache_pnt.second;
                }
                None => {
                    let (uv, g) = projector.value_of_uv(s, &a_p[i], a_tol_working, a_tol2);
                    a_p2d[i] = uv;
                    a_gap = g;
                }
            }
            // `cxx:987`: `aProjector.Gap()`, the value `getLine` leaves behind.
            out.last_gap = a_gap;
            let a_cur_dist = a_gap * a_gap;
            if a_tol2 < a_cur_dist {
                a_tol2 = a_cur_dist;
            }
        }
        if is_periodic_u {
            is_recompute =
                fix_periodicity_troubles(&mut a_p2d, 0, s.u_period(), saved_point_num, saved_point.x());
        }
        if is_periodic_v {
            is_recompute |=
                fix_periodicity_troubles(&mut a_p2d, 1, s.v_period(), saved_point_num, saved_point.y());
        }
    }
    out.is_recompute = is_recompute;
    out.probes2d = a_p2d;

    // `cxx:1015-1018`: a spherical surface with a periodicity jump is given up
    // (the caller's `fixPeriodicityTroubles` comment: `rln S4135 sphere is not
    // considered as V-closed anymore`).
    if is_recompute && s.gp_sphere().is_some() {
        return None;
    }

    a_tol2 = an_old_tol2;

    let d_par = params[a_nb - 1] - params[0];
    if d_par.abs() < PCONFUSION {
        return None;
    }

    // `cxx:1032-1034`: the candidate pcurve is the straight 2D chord, linear in
    // the 3D parameter.
    let a_vec0 = GpVec2d::new(a_p2d[3].x() - a_p2d[0].x(), a_p2d[3].y() - a_p2d[0].y());
    let a_vec = GpVec2d::new(a_vec0.x() / d_par, a_vec0.y() / d_par);
    let candidate = |t: f64| {
        let dt = t - params[0];
        GpPnt2d::new(a_p2d[0].x() + a_vec.x() * dt, a_p2d[0].y() + a_vec.y() * dt)
    };

    // `cxx:1036-1058`: on a C1 surface every sample must lie on the normal line
    // of its candidate `(u, v)`.
    let mut is_normal_check = is_cn_1(s);
    if is_normal_check {
        for i in 0..a_nb {
            let cu = candidate(params[i]);
            let (a_cur_p, a_du, a_dv) = s.d1(cu.x(), cu.y());
            let a_normal_vec = a_du.crossed(&a_dv);
            if a_normal_vec.square_magnitude() < SQUARE_CONFUSION {
                is_normal_check = false;
                break;
            }
            // `gp_Lin(aCurP, gp_Dir(aNormalVec)).Distance(thePoints(i))`.
            let a_dist = a_normal_vec
                .crossed(&GpVec::from_pnts(&a_cur_p, &points[i]))
                .magnitude()
                / a_normal_vec.magnitude();
            if a_dist > a_tol_working {
                return None;
            }
        }
    }

    // `cxx:1060-1076`: otherwise the surface distance must stay constant.
    if !is_normal_check {
        let a_first_point_dist = s
            .value(a_p2d[0].x(), a_p2d[0].y())
            .square_distance(&points[0]);
        a_tol2 = a_tol2.max(a_tol2 * 2.0 * a_first_point_dist);
        for i in 1..a_nb - 1 {
            let cu = candidate(params[i]);
            let a_cur_p = s.d0(cu.x(), cu.y());
            let a_dist1 = a_cur_p.square_distance(&points[i]);
            if (a_first_point_dist - a_dist1).abs() > a_tol2 {
                return None;
            }
        }
    }

    // `cxx:1079-1086`: an exact `Geom2d_Line` when the 2D chord is unit speed.
    let a_l_length = a_vec0.magnitude();
    if (a_l_length - d_par).abs() <= PCONFUSION {
        let a_dir_l = GpVec2d::new(a_vec0.x() / a_l_length, a_vec0.y() / a_l_length);
        let dir2d = GpDir2d::from_vec2d(&a_dir_l).ok()?;
        let loc = GpPnt2d::new(
            a_p2d[0].x() - params[0] * a_dir_l.x(),
            a_p2d[0].y() - params[0] * a_dir_l.y(),
        );
        return Some(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir2d))));
    }

    // `cxx:1088-1101`: a two-pole degree-1 B-spline for the rest.
    let bs = bspline_from_samples(&[a_p2d[0], a_p2d[3]], params[0], params[a_nb - 1], 1).ok()?;
    Some(Arc::new(bs))
}

// ---------------------------------------------------------------------------
// `interpolatePCurve` (`ShapeConstruct_ProjectCurveOnSurface.cxx:2152-2213`) +
// `Geom2dAPI_Interpolate` non-periodic (`Geom2dAPI_Interpolate.cxx:615-810`)
// ---------------------------------------------------------------------------

/// `ShapeConstruct_ProjectCurveOnSurface::checkPoints2d` (`cxx:2373-2460`):
/// drop every point that coincides with its last kept predecessor
/// (`SquareDistance < gp::Resolution()`, `cxx:2396`) and shrink the tolerance
/// to `0.9 * min distance` (`cxx:2419-2422`). Returns the kept points / params
/// and the updated tolerance.
fn check_points2d(points: &[GpPnt2d], params: &[f64], preci: f64) -> (Vec<GpPnt2d>, Vec<f64>, f64) {
    let n = points.len();
    if n == 0 || params.len() != n {
        return (Vec::new(), Vec::new(), preci);
    }
    // `cxx:2382-2390`: `tmpParam(i)` flags the entries that survive.
    let mut keep = vec![true; n];
    let mut nb_pnt_dropped = 0usize;
    let mut last_valid = 0usize;
    let mut dist_min2 = f64::MAX;
    let mut prev = points[0];
    for i in 1..n {
        let cur_dist2 = prev.square_distance(&points[i]);
        if cur_dist2 < f64::MIN_POSITIVE {
            // `cxx:2398-2408`: drop this point, or the previous kept point when
            // this is the last one.
            nb_pnt_dropped += 1;
            if i + 1 == n {
                keep[last_valid] = false;
            } else {
                keep[i] = false;
            }
        } else {
            if cur_dist2 < dist_min2 {
                dist_min2 = cur_dist2;
            }
            last_valid = i;
            prev = points[i];
        }
    }
    let preci = if dist_min2 < f64::MAX {
        // `cxx:2419-2422`
        0.9 * dist_min2.sqrt()
    } else {
        preci
    };
    if nb_pnt_dropped == 0 {
        return (points.to_vec(), params.to_vec(), preci);
    }
    let mut pts: Vec<GpPnt2d> = points.to_vec();
    let new_last = n - nb_pnt_dropped;
    if new_last < 2 {
        // `cxx:2430-2438`: a minimal-length pcurve, the last point shifted by
        // `(preci, preci)`.
        keep[0] = true;
        keep[n - 1] = true;
        let mut p = pts[n - 1];
        p.set_x(p.x() + preci);
        p.set_y(p.y() + preci);
        pts[n - 1] = p;
    }
    let mut out_p = Vec::new();
    let mut out_t = Vec::new();
    for i in 0..n {
        if keep[i] {
            out_p.push(pts[i]);
            out_t.push(params[i]);
        }
    }
    (out_p, out_t, preci)
}

/// `Geom2dAPI_Interpolate.cxx:156-214` `BuildTangents`, with every
/// `TangentFlag` false as the two-array constructor leaves them
/// (`Geom2dAPI_Interpolate.cxx:299-308`): the end tangents are the first
/// derivatives of the Lagrange polynomial through the first / last
/// `degree + 1` points, evaluated at the end parameter (`PLib::EvalLagrange`,
/// `cxx:184-190` / `cxx:202-208`). `None` marks the
/// `Standard_ConstructionError` thrown for fewer than three points
/// (`cxx:170-175`).
fn build_tangents(points: &[GpPnt2d], params: &[f64]) -> Option<(GpVec2d, GpVec2d)> {
    let n = points.len();
    if n < 3 || params.len() != n {
        return None;
    }
    // `cxx:169-178`: degree 2 for three points, 3 otherwise.
    let degree = if n == 3 { 2 } else { 3 };
    let mut vals = vec![[0.0f64; 2]; degree + 1];
    let mut prms = vec![0.0f64; degree + 1];
    for i in 0..=degree {
        vals[i] = [points[i].x(), points[i].y()];
        prms[i] = params[i];
    }
    let first = occt_core::bspl::plib_eval::eval_lagrange(params[0], 1, degree, &vals, &prms).ok()?;
    // `cxx:196-214`: the last `degree + 1` points, evaluated at `params(n)`.
    let base = n - 1 - degree;
    for i in 0..=degree {
        vals[i] = [points[base + i].x(), points[base + i].y()];
        prms[i] = params[base + i];
    }
    let last =
        occt_core::bspl::plib_eval::eval_lagrange(params[n - 1], 1, degree, &vals, &prms).ok()?;
    Some((GpVec2d::new(first[1][0], first[1][1]), GpVec2d::new(last[1][0], last[1][1])))
}

/// `CheckPoints` (`Geom2dAPI_Interpolate.cxx:33-44`): the `Geom2dAPI_Interpolate`
/// constructor throws `Standard_ConstructionError` unless every consecutive
/// pair is at least `tol` apart.
fn check_points(points: &[GpPnt2d], tol: f64) -> bool {
    let tol2 = tol * tol;
    points.windows(2).all(|w| w[0].square_distance(&w[1]) >= tol2)
}

/// `CheckParameters` (`Geom2dAPI_Interpolate.cxx:70-80`): consecutive
/// parameters must grow by at least `RealSmall()`.
fn check_parameters(params: &[f64]) -> bool {
    params.windows(2).all(|w| w[1] - w[0] >= f64::MIN_POSITIVE)
}

/// `Geom2dAPI_Interpolate::PerformNonPeriodic`
/// (`Geom2dAPI_Interpolate.cxx:615-810`) for `myTangentRequest == false`: a
/// single-span non-periodic degree-3 (degree 2 for three points, degree 1 for
/// two) C2 interpolation through the points, solved by `BSplCLib::Interpolate`
/// with contact orders. `None` marks OCCT's `inversion_problem`
/// (`cxx:802`) and the throw paths of the calls it makes.
fn interpolate_non_periodic(points: &[GpPnt2d], params: &[f64]) -> Option<Geom2dBSplineCurve> {
    let n = points.len();
    if n < 2 || params.len() != n {
        return None;
    }
    // `cxx:623-648`
    let (degree, num_poles) = if n == 2 {
        (1usize, n)
    } else if n == 3 {
        (2usize, 3usize)
    } else {
        (3usize, n + 2)
    };
    // `cxx:650-668`: flat knots and contact orders. `Geom2d_BSplineCurve` builds
    // the same flat knots from `knots` + `mults` (`cxx:680`, `:698`, `:806`).
    let mut flatknots = vec![0.0f64; num_poles + degree + 1];
    let mut contact = vec![0i32; num_poles];

    for ii in 1..=degree + 1 {
        flatknots[ii - 1] = params[0];
        flatknots[ii + num_poles - 1] = params[n - 1];
    }

    let mut flat_poles = vec![0.0f64; num_poles * 2];
    match degree {
        1 => {
            // `cxx:675-683`
            for ii in 0..num_poles {
                flat_poles[2 * ii] = points[ii].x();
                flat_poles[2 * ii + 1] = points[ii].y();
            }
        }
        2 => {
            // `cxx:684-701`: `knots = (params(1), params(3))`, `mults = (3, 3)`.
            for ii in 0..num_poles {
                flat_poles[2 * ii] = points[ii].x();
                flat_poles[2 * ii + 1] = points[ii].y();
            }
            occt_core::bspl::banded_interp::interpolate_contact(
                degree,
                &flatknots,
                params,
                &contact,
                &mut flat_poles,
                2,
            )
            .ok()?;
        }
        _ => {
            // `cxx:703-712`: the end tangents.
            let (tan_first, tan_last) = build_tangents(points, params)?;
            // `cxx:717-726`
            contact[1] = 1;
            let mut parameters = vec![0.0f64; num_poles];
            parameters[0] = params[0];
            parameters[1] = params[0];
            flat_poles[0] = points[0].x();
            flat_poles[1] = points[0].y();
            flat_poles[2] = tan_first.x();
            flat_poles[3] = tan_first.y();
            // `cxx:764-770`
            for ii in 0..n {
                parameters[ii + 1] = params[ii];
            }
            for ii in 1..n - 1 {
                flat_poles[2 * (ii + 1)] = points[ii].x();
                flat_poles[2 * (ii + 1) + 1] = points[ii].y();
            }
            for ii in 0..n {
                flatknots[degree + ii] = params[ii];
            }
            // `cxx:785-795`
            let last = num_poles - 1;
            flat_poles[2 * (last - 1)] = tan_last.x();
            flat_poles[2 * (last - 1) + 1] = tan_last.y();
            contact[last - 1] = 1;
            parameters[last] = params[n - 1];
            parameters[last - 1] = params[n - 1];
            flat_poles[2 * last] = points[n - 1].x();
            flat_poles[2 * last + 1] = points[n - 1].y();
            // `cxx:797-806`
            occt_core::bspl::banded_interp::interpolate_contact(
                degree,
                &flatknots,
                &parameters,
                &contact,
                &mut flat_poles,
                2,
            )
            .ok()?;
        }
    }

    let xs: Vec<f64> = (0..num_poles).map(|i| flat_poles[2 * i]).collect();
    let ys: Vec<f64> = (0..num_poles).map(|i| flat_poles[2 * i + 1]).collect();
    Geom2dBSplineCurve::new(xs, ys, flatknots, degree).ok()
}

/// `ShapeConstruct_ProjectCurveOnSurface::interpolatePCurve`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2152-2213`): the `Perform`
/// fallback (`cxx:771`) reached for every pcurve that `approxPCurve` did not
/// settle analytically. `preci` is `myPreci` (`cxx:2158`).
pub(super) fn interpolate_pcurve(
    points: &[GpPnt2d],
    params: &[f64],
    preci: f64,
) -> Option<Arc<dyn Curve2d>> {
    if points.is_empty() || points.len() != params.len() {
        return None;
    }
    // `cxx:2158`: `theTolerance2d = myPreci / (100 * theNbPnt)`.
    let tol2d = preci / (100.0 * points.len() as f64);
    // `cxx:2177-2195`: drop coincident points and shrink the tolerance.
    let (pts, prm, a_preci) = check_points2d(points, params, tol2d);
    if pts.len() < 2 || prm.len() != pts.len() {
        return None;
    }
    // `cxx:2197-2202`: `Geom2dAPI_Interpolate(aPnts2d, aParams, false, aPreci)`
    // then `Perform()`. The constructor's `CheckPoints` / `CheckParameters`
    // throws (`cxx:294-353`) surface as the `catch` at `cxx:2204`.
    if !check_points(&pts, a_preci) || !check_parameters(&prm) {
        return None;
    }
    let curve = interpolate_non_periodic(&pts, &prm)?;
    Some(Arc::new(curve))
}
