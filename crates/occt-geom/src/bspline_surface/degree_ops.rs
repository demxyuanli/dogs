//! `Geom_BSplineSurface::IncreaseDegree(UDegree, VDegree)` — the degree
//! elevation `GeomFill_GordonBuilder::unifySurfaces` calls on each of the three
//! Gordon surfaces.
//!
//! Source: `ModelingData/TKG3d/Geom/Geom_BSplineSurface.cxx:365-520`
//! (OCCT 8.0.0), which forwards to `BSplSLib::IncreaseDegree`
//! (`BSplSLib.cxx:2267-2326`): the pole grid is flattened in memory order by
//! `BSplSLib::SetPoles` and handed to `BSplCLib::IncreaseDegree` with
//! `dim = (3|4) * perpendicular_count`, then unpacked by `BSplSLib::GetPoles`.
//!
//! Translation notes:
//! - [`GeomBSplineSurface::pack_dir`] / [`GeomBSplineSurface::unpack_dir`] are
//!   that `SetPoles` / `GetPoles` pair (`BSplSLib.cxx:1916-2121`).
//! - `dim` is `(3|4) * RowLength` for U and `(3|4) * ColLength` for V
//!   (`BSplSLib.cxx:2299-2306`). `RowLength()` is the column count, i.e. the V
//!   count, and `ColLength()` the U count (`NCollection_Array2.hxx:198-201`
//!   with `NbUPoles = ColLength()`), so a U elevation multiplies by the V
//!   count and a V elevation by the U count.
//! - `increase_degree_dim` takes the same `BSplCLib` dimension-major layout as
//!   `BSplSLib::IncreaseDegree` passes on, so no re-packing happens between the
//!   two, exactly as in OCCT.
//! - `Standard_ConstructionError` for a degree that is lower than the current
//!   one or above `Geom_BSplineSurface::MaxDegree()` (`cxx:369-372`,
//!   `cxx:434-437`) becomes `Err`.

use super::*;
use crate::bspline_curve::MAX_DEGREE;
use occt_core::bspl::increase_degree as bspl_raise;
use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::gp::GpPnt;

impl GeomBSplineSurface {
    /// `BSplSLib::SetPoles` (`BSplSLib.cxx:1916-2013`): the pole grid in
    /// memory order, homogeneous when the surface is rational.
    fn pack_dir(&self) -> Vec<f64> {
        let rational = self.weights.is_some();
        let dim = if rational { 4 } else { 3 };
        let n_u = self.poles.len();
        let n_v = self.poles.first().map_or(0, |r| r.len());
        let mut flat = Vec::with_capacity(n_u * n_v * dim);
        for u in 0..n_u {
            for v in 0..n_v {
                let p = self.poles[u][v];
                let w = match &self.weights {
                    Some(w) => w[u][v],
                    None => 1.0,
                };
                flat.push(p.x() * w);
                flat.push(p.y() * w);
                flat.push(p.z() * w);
                if rational {
                    flat.push(w);
                }
            }
        }
        flat
    }

    /// `BSplSLib::GetPoles` (`BSplSLib.cxx:2021-2121`): the inverse of
    /// [`Self::pack_dir`], dividing the homogeneous coordinates by the weight.
    fn unpack_dir(&mut self, flat: &[f64], n_u: usize, n_v: usize) {
        let rational = self.weights.is_some();
        let dim = if rational { 4 } else { 3 };
        let mut poles = vec![vec![GpPnt::zero(); n_v]; n_u];
        let mut weights = if rational {
            Some(vec![vec![1.0f64; n_v]; n_u])
        } else {
            None
        };
        for u in 0..n_u {
            for v in 0..n_v {
                let base = (u * n_v + v) * dim;
                let w = if rational { flat[base + 3] } else { 1.0 };
                poles[u][v] = GpPnt::new(flat[base] / w, flat[base + 1] / w, flat[base + 2] / w);
                if let Some(ws) = weights.as_mut() {
                    ws[u][v] = w;
                }
            }
        }
        self.poles = poles;
        self.weights = weights;
    }

    /// `Geom_BSplineSurface::IncreaseDegree(UDegree, VDegree)`
    /// (`Geom_BSplineSurface.cxx:365-520`). A degree equal to the current one
    /// is a no-op for that direction.
    pub fn increase_degree(&mut self, u_degree: usize, v_degree: usize) -> Result<(), String> {
        if u_degree != self.deg_u {
            if u_degree < self.deg_u || u_degree > MAX_DEGREE {
                return Err("Geom_BSplineSurface::IncreaseDegree: bad U degree value".into());
            }
            let n_u = self.poles.len();
            let n_v = self.poles.first().map_or(0, |r| r.len());
            let from_k1 = first_u_knot_index(self.deg_u as i32, &self.umults);
            let to_k2 = last_u_knot_index(self.deg_u as i32, &self.umults);
            let step = (u_degree - self.deg_u) as i32;
            let nb_poles = (n_u as i32 + step * (to_k2 - from_k1)).max(0) as usize;
            let nb_knots = bspl_raise::increase_degree_count_knots(
                self.deg_u as i32,
                u_degree as i32,
                self.u_periodic,
                &self.umults,
            )
            .max(0) as usize;

            let flat = self.pack_dir();
            let dim = if self.weights.is_some() { 4 } else { 3 } * n_v;
            let mut new_poles = vec![0.0f64; nb_poles * dim];
            let mut new_knots = vec![0.0f64; nb_knots];
            let mut new_mults = vec![0i32; nb_knots];
            bspl_raise::increase_degree_dim(
                self.deg_u as i32,
                u_degree as i32,
                self.u_periodic,
                dim as i32,
                &flat,
                &self.uknots,
                &self.umults,
                &mut new_poles,
                &mut new_knots,
                &mut new_mults,
            );

            self.unpack_dir(&new_poles, nb_poles, n_v);
            self.deg_u = u_degree;
            self.set_uknots(flat_knots_of(
                &new_knots,
                &new_mults,
                u_degree,
                self.u_periodic,
            ));
        }

        if v_degree != self.deg_v {
            if v_degree < self.deg_v || v_degree > MAX_DEGREE {
                return Err("Geom_BSplineSurface::IncreaseDegree: bad V degree value".into());
            }
            let n_u = self.poles.len();
            let n_v = self.poles.first().map_or(0, |r| r.len());
            let from_k1 = first_u_knot_index(self.deg_v as i32, &self.vmults);
            let to_k2 = last_u_knot_index(self.deg_v as i32, &self.vmults);
            let step = (v_degree - self.deg_v) as i32;
            let nb_poles = (n_v as i32 + step * (to_k2 - from_k1)).max(0) as usize;
            let nb_knots = bspl_raise::increase_degree_count_knots(
                self.deg_v as i32,
                v_degree as i32,
                self.v_periodic,
                &self.vmults,
            )
            .max(0) as usize;

            let flat = self.pack_dir();
            let dim = if self.weights.is_some() { 4 } else { 3 } * n_u;
            let mut new_poles = vec![0.0f64; nb_poles * dim];
            let mut new_knots = vec![0.0f64; nb_knots];
            let mut new_mults = vec![0i32; nb_knots];
            bspl_raise::increase_degree_dim(
                self.deg_v as i32,
                v_degree as i32,
                self.v_periodic,
                dim as i32,
                &flat,
                &self.vknots,
                &self.vmults,
                &mut new_poles,
                &mut new_knots,
                &mut new_mults,
            );

            self.unpack_dir(&new_poles, n_u, nb_poles);
            self.deg_v = v_degree;
            self.set_vknots(flat_knots_of(
                &new_knots,
                &new_mults,
                v_degree,
                self.v_periodic,
            ));
        }
        Ok(())
    }
}

/// The flat knot vector rebuilt from the distinct knots and multiplicities
/// `updateUKnots` / `updateVKnots` expand (`Geom_BSplineSurface.cxx:1148-1198`).
fn flat_knots_of(distinct: &[f64], mults: &[i32], degree: usize, periodic: bool) -> Vec<f64> {
    if periodic {
        occt_core::bspl::knots::knot_sequence_periodic(distinct, mults, degree as i32)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(distinct, mults, degree as i32)
    }
}
