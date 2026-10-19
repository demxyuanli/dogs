//! `Geom_BSplineSurface::RemoveUKnot` / `RemoveVKnot`
//! (`Geom_BSplineSurface_1.cxx:2027-2220`) through `BSplSLib::RemoveKnot`
//! (`BSplSLib.cxx:2195-2263`).

use occt_core::gp::GpPnt;

use crate::bspline_surface::GeomBSplineSurface;

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::RemoveUKnot` (`Geom_BSplineSurface_1.cxx:2027-2108`).
    pub fn remove_u_knot(&mut self, index: i32, m: i32, tolerance: f64) -> Result<bool, &'static str> {
        self.remove_knot_dir(true, index, m, tolerance)
    }

    /// `Geom_BSplineSurface::RemoveVKnot` (`Geom_BSplineSurface_1.cxx:2112-2193`).
    pub fn remove_v_knot(&mut self, index: i32, m: i32, tolerance: f64) -> Result<bool, &'static str> {
        self.remove_knot_dir(false, index, m, tolerance)
    }

    fn remove_knot_dir(
        &mut self,
        u_direction: bool,
        index: i32,
        m: i32,
        tolerance: f64,
    ) -> Result<bool, &'static str> {
        if m < 0 {
            return Ok(true);
        }
        let (knots, mults, periodic, degree) = if u_direction {
            let (k, mlt) = self.distinct_knots_and_mults_u();
            (k, mlt, self.u_periodic, self.deg_u as i32)
        } else {
            let (k, mlt) = self.distinct_knots_and_mults_v();
            (k, mlt, self.v_periodic, self.deg_v as i32)
        };
        let (i1, i2) = knot_index_range(periodic, degree, &mults);
        if !periodic && (index <= i1 || index >= i2) {
            return Err("Geom_BSplineSurface::RemoveUKnot: invalid Index");
        } else if periodic && (index < i1 || index > i2) {
            return Err("Geom_BSplineSurface::RemoveUKnot: invalid Index for periodic case");
        }
        let idx = (index - 1) as usize;
        if idx >= mults.len() {
            return Err("Geom_BSplineSurface::RemoveUKnot: invalid Index");
        }
        let step = mults[idx] - m;
        if step <= 0 {
            return Ok(true);
        }

        let nb_u = self.poles.len();
        let nb_v = self.poles.first().map_or(0, |row| row.len());
        let rational = self.weights.is_some();
        let coord = if rational { 4 } else { 3 };
        let (along, across) = if u_direction { (nb_u, nb_v) } else { (nb_v, nb_u) };
        let along_new = along as i32 - step;
        if along_new < 0 {
            return Err("Geom_BSplineSurface::RemoveUKnot: invalid Index");
        }
        let dim = coord * across as i32;
        let flat = flatten_poles(self, u_direction, rational);
        let nb_new_knots = knots.len() as i32 - if m == 0 { 1 } else { 0 };
        if nb_new_knots < 0 {
            return Err("Geom_BSplineSurface::RemoveUKnot: invalid Index");
        }
        let mut new_flat = vec![0.0; (dim * along_new) as usize];
        let mut new_knots = vec![0.0; nb_new_knots as usize];
        let mut new_mults = vec![0i32; nb_new_knots as usize];
        let ok = occt_core::bspl::remove_knot::remove_knot_dim(
            index,
            m,
            degree,
            periodic,
            dim,
            &flat,
            &knots,
            &mults,
            &mut new_flat,
            &mut new_knots,
            &mut new_mults,
            tolerance,
        );
        if !ok {
            return Ok(false);
        }
        let (poles, weights) = unpack_poles(&new_flat, u_direction, rational, along_new as usize, across);
        // `myMaxDerivInvOk = false` before `updateUKnots()` / `updateVKnots()`
        // (`Geom_BSplineSurface_1.cxx:2105-2106`, `:2190-2191`).
        self.res_latch.invalidate();
        self.poles = poles;
        self.weights = weights;
        if u_direction {
            self.knots_u = flat_knots(periodic, degree, &new_knots, &new_mults);
            self.update_u_knots();
        } else {
            self.knots_v = flat_knots(periodic, degree, &new_knots, &new_mults);
            self.update_v_knots();
        }
        Ok(true)
    }
}

fn knot_index_range(periodic: bool, degree: i32, mults: &[i32]) -> (i32, i32) {
    if periodic {
        (1, mults.len() as i32)
    } else {
        (
            occt_core::bspl::locate::first_u_knot_index(degree, mults).max(1),
            occt_core::bspl::locate::last_u_knot_index(degree, mults).max(1),
        )
    }
}

fn flat_knots(periodic: bool, degree: i32, knots: &[f64], mults: &[i32]) -> Vec<f64> {
    if periodic {
        occt_core::bspl::knots::knot_sequence_periodic(knots, mults, degree)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(knots, mults, degree)
    }
}

/// `BSplSLib::SetPoles` (`BSplSLib.cxx:1916-2017`).
fn flatten_poles(surf: &GeomBSplineSurface, u_direction: bool, rational: bool) -> Vec<f64> {
    let nb_u = surf.poles.len();
    let nb_v = surf.poles.first().map_or(0, |row| row.len());
    let mut flat = Vec::with_capacity(nb_u * nb_v * if rational { 4 } else { 3 });
    let push = |flat: &mut Vec<f64>, u: usize, v: usize| {
        let p = surf.poles[u][v];
        if rational {
            let w = surf.weights.as_ref().map(|g| g[u][v]).unwrap_or(1.0);
            flat.push(p.x() * w);
            flat.push(p.y() * w);
            flat.push(p.z() * w);
            flat.push(w);
        } else {
            flat.push(p.x());
            flat.push(p.y());
            flat.push(p.z());
        }
    };
    if u_direction {
        for u in 0..nb_u {
            for v in 0..nb_v {
                push(&mut flat, u, v);
            }
        }
    } else {
        for v in 0..nb_v {
            for u in 0..nb_u {
                push(&mut flat, u, v);
            }
        }
    }
    flat
}

/// `BSplSLib::GetPoles` (`BSplSLib.cxx:2021-2118`).
fn unpack_poles(
    flat: &[f64],
    u_direction: bool,
    rational: bool,
    along: usize,
    across: usize,
) -> (Vec<Vec<GpPnt>>, Option<Vec<Vec<f64>>>) {
    let (nb_u, nb_v) = if u_direction {
        (along, across)
    } else {
        (across, along)
    };
    let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); nb_v]; nb_u];
    let mut weights = if rational {
        Some(vec![vec![1.0; nb_v]; nb_u])
    } else {
        None
    };
    let stride = if rational { 4 } else { 3 };
    let mut l = 0usize;
    let take = |flat: &[f64], l: &mut usize, rational: bool| -> (GpPnt, f64) {
        if rational {
            let w = flat.get(*l + 3).copied().unwrap_or(1.0);
            let inv = if w == 0.0 { 0.0 } else { 1.0 / w };
            let p = GpPnt::new(
                flat.get(*l).copied().unwrap_or(0.0) * inv,
                flat.get(*l + 1).copied().unwrap_or(0.0) * inv,
                flat.get(*l + 2).copied().unwrap_or(0.0) * inv,
            );
            *l += 4;
            (p, w)
        } else {
            let p = GpPnt::new(
                flat.get(*l).copied().unwrap_or(0.0),
                flat.get(*l + 1).copied().unwrap_or(0.0),
                flat.get(*l + 2).copied().unwrap_or(0.0),
            );
            *l += stride;
            (p, 1.0)
        }
    };
    if u_direction {
        for u in 0..nb_u {
            for v in 0..nb_v {
                let (p, w) = take(flat, &mut l, rational);
                poles[u][v] = p;
                if let Some(grid) = weights.as_mut() {
                    grid[u][v] = w;
                }
            }
        }
    } else {
        for v in 0..nb_v {
            for u in 0..nb_u {
                let (p, w) = take(flat, &mut l, rational);
                poles[u][v] = p;
                if let Some(grid) = weights.as_mut() {
                    grid[u][v] = w;
                }
            }
        }
    }
    (poles, weights)
}
