//! `CachePoint`, `ProjectorCache` and `SurfaceProjectorWithCache` (`cxx:66-217`).

use super::super::prelude::*;
use super::*;

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
    /// `ShapeExtend_DONE4` from the last `approxPCurve` (`cxx:1815`, `cxx:1886`).
    pub done4: bool,
    /// `AdjustOverDegenMode() = false` (`ShapeFix_Wire.cxx:852`): pin `ul`/`vl`
    /// instead of `uf`/`vf`. Default false keeps the constructor mode 1.
    pub adjust_over_degen_off: bool,
}

impl ProjectorCache {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            done4: false,
            adjust_over_degen_off: false,
        }
    }

    /// `myCache.Length()`.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// `myCache.IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub(crate) fn get(&self, j: usize) -> CachePoint {
        self.entries[j]
    }

    /// `cxx:1122-1141` / `cxx:1892-1903`: store the two endpoints of the pcurve
    /// just built. `aChangeCycle` (`cxx:1123-1129`, `cxx:1298-1305`) is set when
    /// the stored anchor `myCache(0)` is closer to the LAST point of this curve
    /// than to its first one, i.e. the curve was sampled in the opposite
    /// direction, and then the roles of `myCache(0)` / `myCache(1)` swap.
    pub(crate) fn store(&mut self, points: &[GpPnt], points2d: &[GpPnt2d], change_cycle: bool) {
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
    pub(crate) fn change_cycle(&self, first: &GpPnt, last: &GpPnt) -> bool {
        if self.entries.is_empty() {
            return false;
        }
        let anchor = self.entries[0].first;
        anchor.distance(first) > anchor.distance(last) && anchor.distance(last) < CONFUSION
    }

    /// `approxPCurve` (`cxx:1403-1421`): the `isRecompute` endpoint arm's
    /// `myCache` scan, first entry whose `SquareDistance` is below `tol2`. The
    /// C++ breaks on the first hit, so a linear scan in cache order is exact.
    pub(crate) fn find(&self, p: &GpPnt, tol2: f64) -> Option<CachePoint> {
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
pub(crate) struct SurfaceProjectorWithCache {
    corners3d: Vec<GpPnt>,
    corners2d: Vec<GpPnt2d>,
}

impl SurfaceProjectorWithCache {
    /// `SurfaceProjectorWithCache::SurfaceProjectorWithCache` (`cxx:84-98`) +
    /// `buildCornerCache` (`cxx:148-200`). `down_cast<Geom_BSplineSurface>`
    /// (`cxx:92-93`) is a dynamic type test, so a `Geom_RectangularTrimmedSurface`
    /// wrapper gets no cache.
    pub(crate) fn new(s: &dyn Surface) -> Self {
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
            occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(bs.knots_u()),
            occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(bs.knots_v()),
        );
        if u_mults.1.is_empty() || v_mults.1.is_empty() {
            return cache;
        }
        // `cxx:156-170`.
        let is_u_first_clamped = u_mults.1[0] >= bs.deg_u() as i32 + 1;
        let is_u_last_clamped = u_mults.1[u_mults.1.len() - 1] >= bs.deg_u() as i32 + 1;
        let is_v_first_clamped = v_mults.1[0] >= bs.deg_v() as i32 + 1;
        let is_v_last_clamped = v_mults.1[v_mults.1.len() - 1] >= bs.deg_v() as i32 + 1;
        let u_first = bs.knots_u()[0];
        let u_last = bs.knots_u()[bs.knots_u().len() - 1];
        let v_first = bs.knots_v()[0];
        let v_last = bs.knots_v()[bs.knots_v().len() - 1];
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
    pub(crate) fn value_of_uv(&self, s: &dyn Surface, p: &GpPnt, tol: f64, tol_sq: f64) -> (GpPnt2d, f64) {
        if let Some(uv) = self.find(p, tol_sq) {
            return next_value_of_uv(s, &uv, p, tol, tol);
        }
        value_of_uv_with_gap(s, p, tol)
    }

    /// `SurfaceProjectorWithCache::NextValueOfUV` (`cxx:125-138`).
    pub(crate) fn next_value_of_uv(
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
