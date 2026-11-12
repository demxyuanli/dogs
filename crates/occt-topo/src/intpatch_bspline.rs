//! B-spline patch wrapper used as a general `Surface`.
//! Source: `Geom_BSplineSurface` helpers previously in `intpatch.rs`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_geom::Surface;
use occt_math::spline_surface::{bspline_surface_d0, interpolate_grid, BSplineSurface};

/// A thin wrapper over `occt_math`'s tensor-product B-spline surface, exposing
/// it as an `occt_geom::Surface` on the unit square `[0, 1]²`.
/// Source: `Geom_BSplineSurface`.
#[derive(Debug, Clone)]
pub struct BsplinePatch {
    /// The underlying B-spline surface (poles, knots, degrees).
    pub bs: BSplineSurface,
}

impl BsplinePatch {
    /// Wrap an existing `BSplineSurface`.
    pub fn new(bs: BSplineSurface) -> Self {
        Self { bs }
    }

    /// Evaluate the patch at `(u, v)`.
    pub fn d0(&self, u: f64, v: f64) -> GpPnt {
        let a = bspline_surface_d0(&self.bs, u, v);
        GpPnt::new(a[0], a[1], a[2])
    }

    /// Convert into a boxed `Surface` handle.
    pub fn to_surface(self) -> Arc<dyn Surface> {
        Arc::new(BsplineSurfaceWrapper { patch: self })
    }
}

/// `occt_geom::Surface` implementation backed by a [`BsplinePatch`].
#[derive(Debug, Clone)]
pub struct BsplineSurfaceWrapper {
    pub patch: BsplinePatch,
}

impl BsplineSurfaceWrapper {
    pub fn new(bs: BSplineSurface) -> Self {
        Self { patch: BsplinePatch::new(bs) }
    }

    /// Evaluate the wrapped surface at `(u, v)`.
    pub fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.patch.d0(u, v)
    }
}

impl BsplineSurfaceWrapper {
    /// Flatten `poles[i][j]` to the row-major layout `BSplSLib` expects.
    fn flat_poles(&self) -> Vec<GpPnt> {
        let mut out = Vec::with_capacity(
            self.patch.bs.poles.len()
                * self.patch.bs.poles.first().map(|r| r.len()).unwrap_or(0),
        );
        for row in &self.patch.bs.poles {
            for pole in row {
                out.push(GpPnt::new(pole[0], pole[1], pole[2]));
            }
        }
        out
    }
}

impl Surface for BsplineSurfaceWrapper {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        let a = bspline_surface_d0(&self.patch.bs, u, v);
        GpPnt::new(a[0], a[1], a[2])
    }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        // `BSplSLib::D1` non-rational arm on the wrapper pole grid.
        let poles = self.flat_poles();
        let n_u = self.patch.bs.poles.len();
        let n_v = self.patch.bs.poles.first().map(|r| r.len()).unwrap_or(0);
        occt_core::bspl::surface::eval_surface_d1(
            &poles,
            n_u,
            n_v,
            &self.patch.bs.knots_u,
            &self.patch.bs.knots_v,
            self.patch.bs.deg_u,
            self.patch.bs.deg_v,
            u,
            v,
        )
    }
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        // `BSplSLib::D2` non-rational arm (`BSplSLib.cxx:1063-1247`).
        let poles = self.flat_poles();
        let n_u = self.patch.bs.poles.len();
        let n_v = self.patch.bs.poles.first().map(|r| r.len()).unwrap_or(0);
        occt_core::bspl::eval::eval_surface_d2(
            &poles,
            n_u,
            n_v,
            &self.patch.bs.knots_u,
            &self.patch.bs.knots_v,
            self.patch.bs.deg_u,
            self.patch.bs.deg_v,
            u,
            v,
        )
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 1.0) }
    fn v_range(&self) -> (f64, f64) { (0.0, 1.0) }
    fn continuity(&self) -> u8 {
        (self.patch.bs.deg_u.min(self.patch.bs.deg_v).saturating_sub(1)).min(3) as u8
    }
    fn transform(&mut self, t: &GpTrsf) {
        for row in &mut self.patch.bs.poles {
            for pole in row {
                let p = GpPnt::new(pole[0], pole[1], pole[2]).transformed(t);
                pole[0] = p.x();
                pole[1] = p.y();
                pole[2] = p.z();
            }
        }
    }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}

/// Fit a B-spline patch through a `points[i][j]` grid (`i` = u rows, `j` = v
/// columns) with degrees `deg_u × deg_v`. Uses clamped uniform knots and
/// two-pass cubic collocation interpolation (from `occt_math`), so the surface
/// passes through every grid point at `u = i/(nu−1)`, `v = j/(nv−1)`.
pub fn fit_bspline_grid(points: &[Vec<GpPnt>], deg_u: usize, deg_v: usize) -> Result<BsplinePatch, String> {
    let nu = points.len();
    let nv = points.first().map(|r| r.len()).unwrap_or(0);
    if nu == 0 || nv == 0 {
        return Err("fit_bspline_grid: empty grid".into());
    }
    let rows: Vec<Vec<[f64; 3]>> = points
        .iter()
        .map(|row| row.iter().map(|p| [p.x(), p.y(), p.z()]).collect())
        .collect();
    let refs: Vec<&[[f64; 3]]> = rows.iter().map(|r| r.as_slice()).collect();
    let bs = interpolate_grid(&refs, nu, nv, deg_u, deg_v)?;
    Ok(BsplinePatch::new(bs))
}

/// Build a `Surface` from a grid of points, suitable as a general (B-spline)
/// face surface. See [`fit_bspline_grid`].
pub fn make_bspline_surface_from_grid(points: &[Vec<GpPnt>], deg_u: usize, deg_v: usize) -> Result<Arc<dyn Surface>, String> {
    Ok(fit_bspline_grid(points, deg_u, deg_v)?.to_surface())
}
