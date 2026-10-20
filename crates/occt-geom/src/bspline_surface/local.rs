//! `Geom_BSplineSurface::LocateParameter` and the local (`[FromK1, ToK2]`)
//! derivative evaluators `LocalD1` / `LocalD2`.

use super::*;
use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpVec};

impl GeomBSplineSurface {
    pub(super) fn locate_param(u: f64, parametric_tolerance: f64, knots: &[f64], periodic: bool) -> (i32, i32) {
        if knots.is_empty() {
            return (0, 1);
        }
        let mut new_u = u;
        if periodic {
            let uf = knots[0];
            let ul = knots[knots.len() - 1];
            new_u = occt_core::bspl::locate::in_period(u, uf, ul);
        }
        let tol = parametric_tolerance.abs();
        let u_first = knots[0];
        let u_last = knots[knots.len() - 1];
        let upper = knots.len() as i32;
        if (new_u - u_first).abs() <= tol {
            return (1, 1);
        }
        if (new_u - u_last).abs() <= tol {
            return (upper, upper);
        }
        if new_u < u_first {
            return (0, 1);
        }
        if new_u > u_last {
            return (upper, upper + 1);
        }
        let mut i1 = occt_core::bspl::locate::hunt_occt(knots, new_u);
        i1 = i1.max(1).min(upper);
        while i1 + 1 <= upper
            && (knots[(i1 as usize)] - new_u).abs() <= tol
        {
            i1 += 1;
        }
        if (knots[(i1 as usize) - 1] - new_u).abs() <= tol {
            (i1, i1)
        } else {
            (i1, i1 + 1)
        }
    }

    /// `Geom_BSplineSurface::LocalD1` (`Geom_BSplineSurface_1.cxx:372-413`).
    pub fn local_d1(
        &self,
        u: f64,
        v: f64,
        from_uk1: i32,
        to_uk2: i32,
        from_vk1: i32,
        to_vk2: i32,
    ) -> (GpPnt, GpVec, GpVec) {
        if from_uk1 == to_uk2 || from_vk1 == to_vk2 {
            return self.d1(u, v);
        }
        // The flat-window arm indexes the flat knot array; on a periodic
        // representation the extension knots make those indices span a period,
        // so use the surface's own periodic-aware D1 (BSplSLib::D1 with the
        // periodic pole wrap).
        if self.has_periodic_direction() {
            return self.d1(u, v);
        }
        let (u_flat, uu) = self.local_flat_index(u, true, from_uk1, to_uk2);
        let (v_flat, vv) = self.local_flat_index(v, false, from_vk1, to_vk2);
        let rat = self.weights.is_some();
        let p = self.d0(uu, vv);
        let du = occt_core::bspl::prepare_eval::surface_dn(
            uu,
            vv,
            1,
            0,
            u_flat,
            v_flat,
            &self.poles,
            self.weights.as_deref(),
            &self.knots_u,
            &self.knots_v,
            None,
            None,
            self.deg_u as i32,
            self.deg_v as i32,
            rat,
            rat,
            self.u_periodic,
            self.v_periodic,
        );
        let dv = occt_core::bspl::prepare_eval::surface_dn(
            uu,
            vv,
            0,
            1,
            u_flat,
            v_flat,
            &self.poles,
            self.weights.as_deref(),
            &self.knots_u,
            &self.knots_v,
            None,
            None,
            self.deg_u as i32,
            self.deg_v as i32,
            rat,
            rat,
            self.u_periodic,
            self.v_periodic,
        );
        (p, du, dv)
    }

    /// `Geom_BSplineSurface::LocalD2` (`Geom_BSplineSurface_1.cxx:417-464`).
    pub fn local_d2(
        &self,
        u: f64,
        v: f64,
        from_uk1: i32,
        to_uk2: i32,
        from_vk1: i32,
        to_vk2: i32,
    ) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        if from_uk1 == to_uk2 || from_vk1 == to_vk2 {
            return self.d2(u, v);
        }
        // See local_d1: periodic surfaces take the periodic-aware D2.
        if self.has_periodic_direction() {
            return self.d2(u, v);
        }
        let (u_flat, uu) = self.local_flat_index(u, true, from_uk1, to_uk2);
        let (v_flat, vv) = self.local_flat_index(v, false, from_vk1, to_vk2);
        let rat = self.weights.is_some();
        let dn = |nu: i32, nv: i32| {
            occt_core::bspl::prepare_eval::surface_dn(
                uu,
                vv,
                nu,
                nv,
                u_flat,
                v_flat,
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                None,
                None,
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            )
        };
        let p = self.d0(uu, vv);
        let du = dn(1, 0);
        let dv = dn(0, 1);
        let duu = if self.deg_u >= 2 {
            dn(2, 0)
        } else {
            GpVec::zero()
        };
        let dvv = if self.deg_v >= 2 {
            dn(0, 2)
        } else {
            GpVec::zero()
        };
        let duv = dn(1, 1);
        (p, du, dv, duu, dvv, duv)
    }

    /// `LocateParameter` on flat knots over `[FromK1, ToK2]` then `FlatIndex`
    /// (`Geom_BSplineSurface_1.cxx:388-392`).
    fn local_flat_index(&self, t: f64, is_u: bool, from_k1: i32, to_k2: i32) -> (i32, f64) {
        let (flat, mults, deg, periodic) = if is_u {
            let (_, m) = Self::unique_knots_mults(&self.knots_u);
            (&self.knots_u[..], m, self.deg_u as i32, self.u_periodic)
        } else {
            let (_, m) = Self::unique_knots_mults(&self.knots_v);
            (&self.knots_v[..], m, self.deg_v as i32, self.v_periodic)
        };
        let (uf, ul) = if periodic && flat.len() > 2 * deg as usize {
            (flat[deg as usize], flat[flat.len() - 1 - deg as usize])
        } else {
            (0.0, 1.0)
        };
        let (idx, new_t) = occt_core::bspl::locate::locate_parameter_range(
            flat, t, periodic, from_k1, to_k2, uf, ul,
        );
        let flat_idx = occt_core::bspl::locate::flat_index(deg, idx, &mults, periodic);
        (flat_idx, new_t)
    }
}
