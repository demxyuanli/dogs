//! Surface sampling to grids and triangulations.
//! Source: `Poly_Triangulation.hxx`, `GCPnts` surface point generators.

use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::poly::Triangulation;
use occt_core::poly::triangulation::Triangle;

/// Rectangular sample grid over a surface.
pub struct UVGrid {
    pub us: Vec<f64>,
    pub vs: Vec<f64>,
    /// `points[iu][iv]` is the surface point at `(us[iu], vs[iv])`.
    pub points: Vec<Vec<GpPnt>>,
}

/// Sample `nu x nv` points of `s`. Non-finite parameter ranges are clamped to
/// `[-1, 1]` so that unbounded surfaces (planes, cylinders) remain sampleable.
pub fn surface_to_grid(s: &dyn Surface, nu: usize, nv: usize) -> UVGrid {
    let (u0, u1) = finite_range(s.u_range());
    let (v0, v1) = finite_range(s.v_range());
    let nu = nu.max(2);
    let nv = nv.max(2);
    let us: Vec<f64> = (0..nu)
        .map(|i| u0 + (u1 - u0) * i as f64 / (nu - 1) as f64)
        .collect();
    let vs: Vec<f64> = (0..nv)
        .map(|j| v0 + (v1 - v0) * j as f64 / (nv - 1) as f64)
        .collect();
    let points = us
        .iter()
        .map(|&u| vs.iter().map(|&v| s.d0(u, v)).collect())
        .collect();
    UVGrid { us, vs, points }
}

/// Flatten a grid row-major into a node array plus two triangles per cell.
/// The split diagonal is chosen so triangle windings stay consistent with the
/// cell's natural normal `(u+1,v)-(u,v) x (u,v+1)-(u,v)`.
pub fn grid_to_triangles(grid: &UVGrid) -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    let nu = grid.us.len();
    let nv = grid.vs.len();
    let mut nodes = Vec::with_capacity(nu * nv);
    for iu in 0..nu {
        for iv in 0..nv {
            nodes.push(grid.points[iu][iv]);
        }
    }
    let mut tris = Vec::with_capacity(2 * (nu - 1) * (nv - 1));
    for iu in 0..nu - 1 {
        for iv in 0..nv - 1 {
            let i00 = iu * nv + iv;
            let i10 = (iu + 1) * nv + iv;
            let i01 = iu * nv + (iv + 1);
            let i11 = (iu + 1) * nv + (iv + 1);
            let p00 = grid.points[iu][iv];
            let p10 = grid.points[iu + 1][iv];
            let p01 = grid.points[iu][iv + 1];
            let n = p10
                .coord
                .subtracted(&p00.coord)
                .crossed(&p01.coord.subtracted(&p00.coord));
            let (sx, sy, sz) = (n.x.abs(), n.y.abs(), n.z.abs());
            let flip = if sx >= sy && sx >= sz {
                n.x < 0.0
            } else if sy >= sz {
                n.y < 0.0
            } else {
                n.z < 0.0
            };
            if flip {
                tris.push((i00, i10, i01));
                tris.push((i10, i11, i01));
            } else {
                tris.push((i00, i10, i11));
                tris.push((i00, i11, i01));
            }
        }
    }
    (nodes, tris)
}

/// Build a [`Triangulation`] from a uniform `nu x nv` surface grid.
pub fn surface_to_triangulation(s: &dyn Surface, nu: usize, nv: usize) -> Triangulation {
    let grid = surface_to_grid(s, nu, nv);
    let (nodes, tris) = grid_to_triangles(&grid);
    let triangles = tris
        .into_iter()
        .map(|(a, b, c)| Triangle::new(a, b, c))
        .collect();
    Triangulation::new(nodes, triangles)
}

/// Numeric surface patch area via the midpoint rule on `|du x dv|`.
///
/// Partial derivatives are obtained by central finite differences of `d0`
/// rather than `d1`, because several surface implementations (e.g.
/// `GeomSphere`) currently return zero derivative vectors from `d1`.
pub fn surface_patch_area_numeric(s: &dyn Surface, nu: usize, nv: usize) -> f64 {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
        return f64::INFINITY;
    }
    let (nu, nv) = (nu.max(2), nv.max(2));
    let hu = (u1 - u0) / nu as f64;
    let hv = (v1 - v0) / nv as f64;
    let span = (u1 - u0).min(v1 - v0);
    let h = 1e-6 * span.max(1e-12);
    let mut area = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (i as f64 + 0.5) * hu;
            let v = v0 + (j as f64 + 0.5) * hv;
            let du = fd_u(s, u, v, h);
            let dv = fd_v(s, u, v, h);
            area += du.coord.crossed(&dv.coord).modulus() * hu * hv;
        }
    }
    area
}

fn finite_range((a, b): (f64, f64)) -> (f64, f64) {
    if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

fn fd_u(s: &dyn Surface, u: f64, v: f64, h: f64) -> GpVec {
    let p1 = s.d0(u + h, v);
    let p2 = s.d0(u - h, v);
    GpVec::new(
        (p1.x() - p2.x()) / (2.0 * h),
        (p1.y() - p2.y()) / (2.0 * h),
        (p1.z() - p2.z()) / (2.0 * h),
    )
}

fn fd_v(s: &dyn Surface, u: f64, v: f64, h: f64) -> GpVec {
    let p1 = s.d0(u, v + h);
    let p2 = s.d0(u, v - h);
    GpVec::new(
        (p1.x() - p2.x()) / (2.0 * h),
        (p1.y() - p2.y()) / (2.0 * h),
        (p1.z() - p2.z()) / (2.0 * h),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::GeomSphere;
    use occt_core::gp::{GpPln, GpSphere as GpSphereCore, GpTrsf, GpVec};

    /// A finite rectangular patch of the XY plane, so area = (u1-u0)*(v1-v0).
    #[derive(Clone)]
    struct PlanePatch {
        pl: GpPln,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    }

    impl Surface for PlanePatch {
        fn d0(&self, u: f64, v: f64) -> GpPnt {
            let xd = self.pl.pos.x_direction().xyz();
            let yd = self.pl.pos.y_direction().xyz();
            let loc = self.pl.location().coord;
            GpPnt::from_xyz(&loc.added(&xd.multiplied(u)).added(&yd.multiplied(v)))
        }
        fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
            (
                self.d0(u, v),
                GpVec::from_xyz(self.pl.pos.x_direction().xyz()),
                GpVec::from_xyz(self.pl.pos.y_direction().xyz()),
            )
        }
        fn u_range(&self) -> (f64, f64) {
            (self.u0, self.u1)
        }
        fn v_range(&self) -> (f64, f64) {
            (self.v0, self.v1)
        }
        fn continuity(&self) -> u8 {
            3
        }
        fn transform(&mut self, _t: &GpTrsf) {}
        fn clone_dyn(&self) -> Box<dyn Surface> {
            Box::new(self.clone())
        }
    }

    #[test]
    fn plane_area_matches_params() {
        let patch = PlanePatch {
            pl: GpPln::default(),
            u0: 0.0,
            u1: 2.0,
            v0: 0.0,
            v1: 3.0,
        };
        let area = surface_patch_area_numeric(&patch, 16, 16);
        assert!((area - 6.0).abs() < 1e-6, "area={area}");
    }

    #[test]
    fn sphere_area_within_five_percent() {
        let s = GeomSphere::new(GpSphereCore::new(occt_core::gp::GpAx3::standard(), 2.0).unwrap());
        let area = surface_patch_area_numeric(&s, 32, 32);
        let expected = 4.0 * std::f64::consts::PI * 4.0; // 4*pi*r^2, r=2
        assert!(
            (area - expected).abs() / expected < 0.05,
            "area={area}, expected={expected}"
        );
    }

    #[test]
    fn grid_and_triangulation() {
        let patch = PlanePatch {
            pl: GpPln::default(),
            u0: 0.0,
            u1: 1.0,
            v0: 0.0,
            v1: 1.0,
        };
        let grid = surface_to_grid(&patch, 4, 3);
        assert_eq!(grid.us.len(), 4);
        assert_eq!(grid.vs.len(), 3);
        assert_eq!(grid.points.len(), 4);
        assert_eq!(grid.points[0].len(), 3);
        assert!((grid.points[3][2].distance(&GpPnt::new(1., 1., 0.))).abs() < 1e-12);

        let (nodes, tris) = grid_to_triangles(&grid);
        assert_eq!(nodes.len(), 12);
        assert_eq!(tris.len(), 2 * 3 * 2); // (4-1)*(3-1) cells * 2

        let tri = surface_to_triangulation(&patch, 4, 3);
        assert_eq!(tri.nb_nodes(), 12);
        assert_eq!(tri.nb_triangles(), 12);
    }
}
