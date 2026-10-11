//! `Geom_BezierSurface::ExchangeUV` and `Increase`.
//!
//! `Increase` packs the grid with `BSplSLib::SetPoles` (`BSplSLib.cxx:1916-2017`),
//! calls the ported `BSplCLib::IncreaseDegree`, then unpacks with `GetPoles`
//! (`BSplSLib.cxx:2021-2122`).

use super::{bezier_flat_knots, GeomBezierSurface};
use occt_core::bspl::bspline_surface_resolution;
use occt_core::bspl::increase_degree::{increase_degree_count_knots, increase_degree_dim};
use occt_core::gp::GpPnt;

impl GeomBezierSurface {
    /// `Geom_BezierSurface::Resolution` (`Geom_BezierSurface.cxx:1991-2037`).
    /// Calls `BSplSLib::Resolution` with the Bezier knot sequence. The cxx
    /// caches `myUMaxDerivInv` with tolerance 1 then multiplies; the ported
    /// `bspline_surface_resolution` returns the product for `r3d` directly.
    pub fn resolution(&self, r3d: f64) -> (f64, f64) {
        let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); self.n_v]; self.n_u];
        let mut weights = None;
        if self.u_rational || self.v_rational {
            if let Some(src) = self.weights.as_ref() {
                let mut grid = vec![vec![0.0; self.n_v]; self.n_u];
                for i in 0..self.n_u {
                    for j in 0..self.n_v {
                        poles[i][j] = self.poles[i * self.n_v + j];
                        grid[i][j] = src[i * self.n_v + j];
                    }
                }
                weights = Some(grid);
            }
        } else {
            for i in 0..self.n_u {
                for j in 0..self.n_v {
                    poles[i][j] = self.poles[i * self.n_v + j];
                }
            }
        }
        let deg_u = self.u_degree() as i32;
        let deg_v = self.v_degree() as i32;
        bspline_surface_resolution(
            &poles,
            weights.as_deref(),
            &bezier_flat_knots(deg_u),
            &bezier_flat_knots(deg_v),
            deg_u,
            deg_v,
            self.u_rational,
            self.v_rational,
            r3d,
        )
    }
}

impl GeomBezierSurface {
    /// `Geom_BezierSurface::ExchangeUV` (`Geom_BezierSurface.cxx:490-528`).
    /// `myMaxDerivInvOk` (`:527`) is not stored.
    pub fn exchange_uv(&mut self) {
        let old_nu = self.n_u;
        let old_nv = self.n_v;
        let new_nu = old_nv;
        let new_nv = old_nu;
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); new_nu * new_nv];
        let rational = self.u_rational || self.v_rational;
        let mut weights = rational.then(|| vec![0.0; new_nu * new_nv]);
        for ii in 0..new_nu {
            for jj in 0..new_nv {
                let src = jj * old_nv + ii;
                let dst = ii * new_nv + jj;
                poles[dst] = self.poles[src];
                if let (Some(dst_w), Some(src_w)) = (weights.as_mut(), self.weights.as_ref()) {
                    dst_w[dst] = src_w[src];
                }
            }
        }
        self.poles = poles;
        self.n_u = new_nu;
        self.n_v = new_nv;
        self.weights = weights;
        std::mem::swap(&mut self.u_rational, &mut self.v_rational);
    }

    /// `Geom_BezierSurface::Increase` (`Geom_BezierSurface.cxx:532-648`).
    /// `myMaxDerivInvOk` (`:647`) is not stored.
    pub fn increase(&mut self, u_deg: i32, v_deg: i32) -> Result<(), &'static str> {
        let old_u = self.u_degree() as i32;
        let old_v = self.v_degree() as i32;
        let max = Self::max_degree();
        if u_deg < old_u || u_deg > max || v_deg < old_v || v_deg > max {
            return Err("Geom_BezierSurface::Increase");
        }
        if u_deg == old_u && v_deg == old_v {
            return Ok(());
        }
        if u_deg > old_u {
            self.increase_direction(true, old_u, u_deg);
        }
        if v_deg > old_v {
            self.increase_direction(false, old_v, v_deg);
        }
        if self.u_rational || self.v_rational {
            self.store_rational_flags();
        } else {
            self.weights = None;
        }
        Ok(())
    }

    /// One direction of `BSplSLib::IncreaseDegree` (`BSplSLib.cxx:2267-2327`).
    fn increase_direction(&mut self, u_direction: bool, old_deg: i32, new_deg: i32) {
        let rational = self.u_rational || self.v_rational;
        let comp = if rational { 4 } else { 3 };
        let nu = self.n_u;
        let nv = self.n_v;
        let mut flat = vec![0.0; comp * nu * nv];
        pack_poles(&self.poles, self.weights.as_deref(), nu, nv, u_direction, &mut flat);
        let span = if u_direction { nv } else { nu };
        let dim = (comp * span) as i32;
        let knots = [0.0, 1.0];
        let mults = [old_deg + 1, old_deg + 1];
        let n_knots = increase_degree_count_knots(old_deg, new_deg, false, &mults).max(0) as usize;
        let new_along = (new_deg as usize) + 1;
        let (new_nu, new_nv) = if u_direction {
            (new_along, nv)
        } else {
            (nu, new_along)
        };
        let mut new_flat = vec![0.0; comp * new_nu * new_nv];
        let mut new_knots = vec![0.0; n_knots];
        let mut new_mults = vec![0i32; n_knots];
        increase_degree_dim(
            old_deg,
            new_deg,
            false,
            dim,
            &flat,
            &knots,
            &mults,
            &mut new_flat,
            &mut new_knots,
            &mut new_mults,
        );
        let (poles, weights) = unpack_poles(&new_flat, new_nu, new_nv, u_direction, rational);
        self.poles = poles;
        self.n_u = new_nu;
        self.n_v = new_nv;
        if rational {
            self.weights = Some(weights);
        }
    }
}

fn pack_poles(
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    nu: usize,
    nv: usize,
    u_direction: bool,
    out: &mut [f64],
) {
    let comp = if weights.is_some() { 4 } else { 3 };
    let mut l = 0usize;
    if u_direction {
        for i in 0..nu {
            for j in 0..nv {
                write_pole(poles, weights, i * nv + j, out, &mut l, comp);
            }
        }
    } else {
        for j in 0..nv {
            for i in 0..nu {
                write_pole(poles, weights, i * nv + j, out, &mut l, comp);
            }
        }
    }
}

fn write_pole(
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    index: usize,
    out: &mut [f64],
    l: &mut usize,
    comp: usize,
) {
    let p = poles[index];
    if let Some(weights) = weights {
        let w = weights[index];
        out[*l] = p.x() * w;
        out[*l + 1] = p.y() * w;
        out[*l + 2] = p.z() * w;
        out[*l + 3] = w;
    } else {
        out[*l] = p.x();
        out[*l + 1] = p.y();
        out[*l + 2] = p.z();
    }
    *l += comp;
}

fn unpack_poles(
    flat: &[f64],
    nu: usize,
    nv: usize,
    u_direction: bool,
    rational: bool,
) -> (Vec<GpPnt>, Vec<f64>) {
    let comp = if rational { 4 } else { 3 };
    let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); nu * nv];
    let mut weights = vec![1.0; nu * nv];
    let mut l = 0usize;
    let mut store = |row: usize, col: usize| {
        let index = row * nv + col;
        if rational {
            let w = flat.get(l + 3).copied().unwrap_or(0.0);
            weights[index] = w;
            poles[index] = GpPnt::new(
                flat.get(l).copied().unwrap_or(0.0) / w,
                flat.get(l + 1).copied().unwrap_or(0.0) / w,
                flat.get(l + 2).copied().unwrap_or(0.0) / w,
            );
        } else {
            poles[index] = GpPnt::new(
                flat.get(l).copied().unwrap_or(0.0),
                flat.get(l + 1).copied().unwrap_or(0.0),
                flat.get(l + 2).copied().unwrap_or(0.0),
            );
        }
        l += comp;
    };
    if u_direction {
        for i in 0..nu {
            for j in 0..nv {
                store(i, j);
            }
        }
    } else {
        for j in 0..nv {
            for i in 0..nu {
                store(i, j);
            }
        }
    }
    (poles, weights)
}
