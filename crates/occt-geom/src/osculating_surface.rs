//! `Geom_OsculatingSurface`. Source: `Geom_OsculatingSurface.cxx`.

use std::sync::Arc;

use occt_core::bspl::build_cache::build_cache;
use occt_core::bspl::convert_grid_poly::grid_polynomial_to_poles;
use occt_core::bspl::locate::hunt_occt;
use occt_core::bspl::plib::{u_trimming, v_trimming};
use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;

use crate::bspline_surface::GeomBSplineSurface;
use crate::surface::Surface;

/// `Geom_OsculatingSurface` (`Geom_OsculatingSurface.pxx`).
#[derive(Clone)]
pub struct OsculatingSurface {
    basis: Arc<dyn Surface>,
    tol: f64,
    oscul1: Vec<GeomBSplineSurface>,
    oscul2: Vec<GeomBSplineSurface>,
    kdeg: Vec<i32>,
    along: [bool; 4],
}

impl OsculatingSurface {
    /// `Geom_OsculatingSurface::Init` (`cxx:105-391`).
    pub fn new(basis: Arc<dyn Surface>, tol: f64) -> Self {
        let mut s = Self {
            basis: basis.clone(),
            tol,
            oscul1: Vec::new(),
            oscul2: Vec::new(),
            kdeg: Vec::new(),
            along: [false; 4],
        };
        s.init();
        s
    }

    fn clear_flags(&mut self) {
        self.along = [false; 4];
    }

    fn is_along_u(&self) -> bool {
        self.along[0] || self.along[1]
    }

    fn is_along_v(&self) -> bool {
        self.along[2] || self.along[3]
    }

    /// `isQPunctual` (`cxx:784-832`).
    fn is_q_punctual(surf: &dyn Surface, param: f64, iso_v: bool, tol_min: f64, tol_max: f64) -> bool {
        let (u1, u2) = surf.u_range();
        let (v1, v2) = surf.v_range();
        if iso_v {
            let step = (u2 - u1) / 10.0;
            let mut d1_max = 0.0f64;
            let mut t = u1;
            while t <= u2 {
                let (_, d1u, _) = surf.d1(t, param);
                d1_max = d1_max.max(d1u.magnitude());
                t += step;
            }
            !(d1_max > tol_max || d1_max < tol_min)
        } else {
            let step = (v2 - v1) / 10.0;
            let mut d1_max = 0.0f64;
            let mut t = v1;
            while t <= v2 {
                let (_, _, d1v) = surf.d1(param, t);
                d1_max = d1_max.max(d1v.magnitude());
                t += step;
            }
            !(d1_max > tol_max || d1_max < tol_min)
        }
    }

    fn init(&mut self) {
        self.clear_flags();
        self.oscul1.clear();
        self.oscul2.clear();
        self.kdeg.clear();
        let Some(init_surf) = self.basis.osculating_bspline() else {
            self.clear_flags();
            return;
        };
        let (u1, u2) = init_surf.u_range();
        let (v1, v2) = init_surf.v_range();
        let tol_min = 0.0;
        self.along[0] = Self::is_q_punctual(&init_surf, v1, true, tol_min, self.tol);
        self.along[1] = Self::is_q_punctual(&init_surf, v2, true, tol_min, self.tol);
        self.along[2] = Self::is_q_punctual(&init_surf, u1, false, tol_min, self.tol);
        self.along[3] = Self::is_q_punctual(&init_surf, u2, false, tol_min, self.tol);
        if !(self.along[0] || self.along[1] || self.along[2] || self.along[3]) {
            return;
        }
        if self.is_along_u() && self.is_along_v() {
            self.clear_flags();
            return;
        }
        if !((self.is_along_u() && init_surf.deg_v > 1) || (self.is_along_v() && init_surf.deg_u > 1))
        {
            self.clear_flags();
            return;
        }
        let (u_knots, _) = GeomBSplineSurface::unique_knots_mults(&init_surf.knots_u);
        let (v_knots, _) = GeomBSplineSurface::unique_knots_mults(&init_surf.knots_v);
        let nb_uk = u_knots.len() as i32;
        let nb_vk = v_knots.len() as i32;
        let mut oscul_ok = true;
        if self.along[0] || self.along[1] {
            for i in 1..nb_uk {
                if self.along[0] {
                    let mut s = init_surf.clone();
                    let mut is_qp = true;
                    let mut uk = i;
                    let mut vk = 1i32;
                    let mut last = None;
                    while is_qp {
                        match self.build_osculating(v1, uk, vk, &s) {
                            Some(l) => {
                                is_qp = Self::is_q_punctual(&l, v1, true, 0.0, self.tol);
                                uk = 1;
                                vk = 1;
                                s = l.clone();
                                last = Some(l);
                            }
                            None => {
                                oscul_ok = false;
                                break;
                            }
                        }
                    }
                    if oscul_ok {
                        if let Some(l) = last {
                            self.oscul1.push(l);
                        }
                    } else {
                        self.clear_flags();
                        return;
                    }
                    if self.along[1] && oscul_ok {
                        let mut s = init_surf.clone();
                        let mut is_qp = true;
                        let mut uk = i;
                        let mut vk = nb_vk - 1;
                        let mut k = 0i32;
                        let mut last = None;
                        while is_qp {
                            match self.build_osculating(v2, uk, vk, &s) {
                                Some(l) => {
                                    k += 1;
                                    is_qp = Self::is_q_punctual(&l, v2, true, 0.0, self.tol);
                                    uk = 1;
                                    vk = 1;
                                    s = l.clone();
                                    last = Some(l);
                                }
                                None => {
                                    oscul_ok = false;
                                    break;
                                }
                            }
                        }
                        if oscul_ok {
                            if let Some(l) = last {
                                self.oscul2.push(l);
                                self.kdeg.push(k);
                            }
                        }
                    }
                } else {
                    let mut s = init_surf.clone();
                    let mut is_qp = true;
                    let mut uk = i;
                    let mut vk = nb_vk - 1;
                    let mut k = 0i32;
                    let mut last = None;
                    while is_qp {
                        match self.build_osculating(v2, uk, vk, &s) {
                            Some(l) => {
                                k += 1;
                                is_qp = Self::is_q_punctual(&l, v2, true, 0.0, self.tol);
                                uk = 1;
                                vk = 1;
                                s = l.clone();
                                last = Some(l);
                            }
                            None => {
                                oscul_ok = false;
                                break;
                            }
                        }
                    }
                    if oscul_ok {
                        if let Some(l) = last {
                            self.oscul2.push(l);
                            self.kdeg.push(k);
                        }
                    } else {
                        self.clear_flags();
                        return;
                    }
                }
            }
        }
        if self.along[2] || self.along[3] {
            for i in 1..nb_vk {
                if self.along[2] {
                    let mut s = init_surf.clone();
                    let mut is_qp = true;
                    let mut uk = 1i32;
                    let mut vk = i;
                    let mut last = None;
                    while is_qp {
                        match self.build_osculating(u1, uk, vk, &s) {
                            Some(l) => {
                                is_qp = Self::is_q_punctual(&l, u1, false, 0.0, self.tol);
                                uk = 1;
                                vk = 1;
                                s = l.clone();
                                last = Some(l);
                            }
                            None => {
                                oscul_ok = false;
                                break;
                            }
                        }
                    }
                    if oscul_ok {
                        if let Some(l) = last {
                            self.oscul1.push(l);
                        }
                    } else {
                        self.clear_flags();
                        return;
                    }
                    if self.along[3] && oscul_ok {
                        let mut s = init_surf.clone();
                        let mut is_qp = true;
                        let mut uk = nb_uk - 1;
                        let mut vk = i;
                        let mut k = 0i32;
                        let mut last = None;
                        while is_qp {
                            match self.build_osculating(u2, uk, vk, &s) {
                                Some(l) => {
                                    k += 1;
                                    is_qp = Self::is_q_punctual(&l, u2, false, 0.0, self.tol);
                                    uk = 1;
                                    vk = 1;
                                    s = l.clone();
                                    last = Some(l);
                                }
                                None => {
                                    oscul_ok = false;
                                    break;
                                }
                            }
                        }
                        if oscul_ok {
                            if let Some(l) = last {
                                self.oscul2.push(l);
                                self.kdeg.push(k);
                            }
                        }
                    }
                } else {
                    let mut s = init_surf.clone();
                    let mut is_qp = true;
                    let mut uk = nb_uk - 1;
                    let mut vk = i;
                    let mut k = 0i32;
                    let mut last = None;
                    while is_qp {
                        match self.build_osculating(u2, uk, vk, &s) {
                            Some(l) => {
                                k += 1;
                                is_qp = Self::is_q_punctual(&l, u2, false, 0.0, self.tol);
                                uk = 1;
                                vk = 1;
                                s = l.clone();
                                last = Some(l);
                            }
                            None => {
                                oscul_ok = false;
                                break;
                            }
                        }
                    }
                    if oscul_ok {
                        if let Some(l) = last {
                            self.oscul2.push(l);
                            self.kdeg.push(k);
                        }
                    } else {
                        self.clear_flags();
                        return;
                    }
                }
            }
        }
    }

    /// `buildOsculatingSurface` (`cxx:532-779`).
    fn build_osculating(
        &self,
        param: f64,
        su_knot: i32,
        sv_knot: i32,
        bs: &GeomBSplineSurface,
    ) -> Option<GeomBSplineSurface> {
        let udeg = bs.deg_u as i32;
        let vdeg = bs.deg_v as i32;
        if (self.is_along_u() && vdeg <= 1) || (self.is_along_v() && udeg <= 1) {
            return None;
        }
        let (u_knots, _) = GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
        let (v_knots, _) = GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
        if su_knot < 1 || sv_knot < 1 {
            return None;
        }
        let ui = (su_knot - 1) as usize;
        let vi = (sv_knot - 1) as usize;
        if ui + 1 >= u_knots.len() || vi + 1 >= v_knots.len() {
            return None;
        }
        let mut osc_u = 0i32;
        let mut osc_v = 0i32;
        if self.is_along_u() {
            osc_u = udeg + 1;
            osc_v = vdeg;
        }
        if self.is_along_v() {
            osc_u = udeg;
            osc_v = vdeg + 1;
        }
        if osc_u * osc_v * 3 == 0 {
            return None;
        }
        let ucache = u_knots[ui];
        let vcache = v_knots[vi];
        let uspan = u_knots[ui + 1] - u_knots[ui];
        let vspan = v_knots[vi + 1] - v_knots[vi];
        let is_v_neg = param > vcache + vspan / 2.0;
        let is_u_neg = param > ucache + uspan / 2.0;
        let mut ucache_p = ucache;
        let mut vcache_p = vcache;
        if self.is_along_u() && param > vcache + vspan / 2.0 {
            vcache_p = vcache + vspan;
        }
        if self.is_along_v() && param > ucache + uspan / 2.0 {
            ucache_p = ucache + uspan;
        }
        let cache = build_cache(
            ucache_p,
            vcache_p,
            uspan,
            vspan,
            false,
            false,
            udeg,
            vdeg,
            0,
            0,
            &bs.knots_u,
            &bs.knots_v,
            &bs.poles,
            None,
        );
        let mut osc = vec![vec![GpPnt::new(0.0, 0.0, 0.0); osc_v.max(0) as usize]; osc_u.max(0) as usize];
        if self.is_along_u() {
            if udeg > vdeg {
                for n in 1..=udeg + 1 {
                    for m in 1..=vdeg {
                        osc[(n - 1) as usize][(m - 1) as usize] =
                            cache_at(&cache, n, m + 1);
                    }
                }
            } else {
                for n in 1..=udeg + 1 {
                    for m in 1..=vdeg {
                        osc[(n - 1) as usize][(m - 1) as usize] =
                            cache_at(&cache, m + 1, n);
                    }
                }
            }
            if is_v_neg {
                v_trimming(-1.0, 0.0, &mut osc);
            }
        }
        if self.is_along_v() {
            if udeg > vdeg {
                for n in 1..=udeg {
                    for m in 1..=vdeg + 1 {
                        osc[(n - 1) as usize][(m - 1) as usize] =
                            cache_at(&cache, n + 1, m);
                    }
                }
            } else {
                for n in 1..=udeg {
                    for m in 1..=vdeg + 1 {
                        osc[(n - 1) as usize][(m - 1) as usize] =
                            cache_at(&cache, m, n + 1);
                    }
                }
            }
            if is_u_neg {
                u_trimming(-1.0, 0.0, &mut osc);
            }
        }
        let mut coeffs = Vec::with_capacity((osc_u * osc_v * 3) as usize);
        if self.is_along_u() {
            for n in 1..=udeg + 1 {
                for m in 1..=vdeg {
                    let p = osc[(n - 1) as usize][(m - 1) as usize];
                    coeffs.push(p.x());
                    coeffs.push(p.y());
                    coeffs.push(p.z());
                }
            }
        }
        if self.is_along_v() {
            coeffs.clear();
            for n in 1..=udeg {
                for m in 1..=vdeg + 1 {
                    let p = osc[(n - 1) as usize][(m - 1) as usize];
                    coeffs.push(p.x());
                    coeffs.push(p.y());
                    coeffs.push(p.z());
                }
            }
        }
        let mut max_u = udeg;
        let mut max_v = vdeg;
        if self.is_along_u() {
            max_v -= 1;
        }
        if self.is_along_v() {
            max_u -= 1;
        }
        let data = grid_polynomial_to_poles(
            -1,
            -1,
            max_u,
            max_v,
            osc_u,
            osc_v,
            &coeffs,
            &[0.0, 1.0],
            &[0.0, 1.0],
            &[u_knots[ui], u_knots[ui + 1]],
            &[v_knots[vi], v_knots[vi + 1]],
        )?;
        GeomBSplineSurface::from_poles_knots_mults(
            data.poles,
            data.u_knots,
            data.v_knots,
            data.u_mults,
            data.v_mults,
            data.u_degree.max(0) as usize,
            data.v_degree.max(0) as usize,
        )
        .ok()
    }

    /// `UOsculatingSurface` (`cxx:395-461`).
    pub fn u_osculating(&self, u: f64, v: f64) -> (bool, bool, Option<GeomBSplineSurface>) {
        if !(self.along[0] || self.along[1]) {
            return (false, false, None);
        }
        let Some(bs) = self.basis.osculating_bspline() else {
            return (false, false, None);
        };
        let (u_knots, _) = GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
        let (v_knots, _) = GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
        let nb_uk = u_knots.len() as i32;
        let nb_vk = v_knots.len() as i32;
        let mut nu = hunt_occt(&u_knots, u);
        let nv = hunt_occt(&v_knots, v);
        if nu < 1 {
            nu = 1;
        }
        if nu >= nb_uk {
            nu = nb_uk - 1;
        }
        let mut skip_second = false;
        if nb_vk == 2 && nv == 1 && v_knots.len() >= 2 {
            if v_knots[nb_vk as usize - 1] - v > v - v_knots[0] {
                skip_second = true;
            }
        }
        let mut opposite = false;
        let mut along = false;
        let mut out = None;
        if self.along[0] && nv == 1 {
            out = self.oscul1.get((nu - 1) as usize).cloned();
            along = out.is_some();
        }
        if self.along[1] && nv == nb_vk - 1 && !skip_second {
            if self.kdeg.get((nu - 1) as usize).copied().unwrap_or(0) % 2 != 0 {
                opposite = true;
            }
            out = self.oscul2.get((nu - 1) as usize).cloned();
            along = out.is_some();
        }
        (along, opposite, out)
    }

    /// `VOsculatingSurface` (`cxx:465-528`).
    pub fn v_osculating(&self, u: f64, v: f64) -> (bool, bool, Option<GeomBSplineSurface>) {
        if !(self.along[2] || self.along[3]) {
            return (false, false, None);
        }
        let Some(bs) = self.basis.osculating_bspline() else {
            return (false, false, None);
        };
        let (u_knots, _) = GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
        let (v_knots, _) = GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
        let nb_uk = u_knots.len() as i32;
        let nb_vk = v_knots.len() as i32;
        let nu = hunt_occt(&u_knots, u);
        let mut nv = hunt_occt(&v_knots, v);
        if nv < 1 {
            nv = 1;
        }
        if nv >= nb_vk {
            nv = nb_vk - 1;
        }
        let mut skip_second = false;
        if nb_uk == 2 && nu == 1 && u_knots.len() >= 2 {
            if u_knots[nb_uk as usize - 1] - u > u - u_knots[0] {
                skip_second = true;
            }
        }
        let mut opposite = false;
        let mut along = false;
        let mut out = None;
        if self.along[2] && nu == 1 {
            out = self.oscul1.get((nv - 1) as usize).cloned();
            along = out.is_some();
        }
        if self.along[3] && nu == nb_uk - 1 && !skip_second {
            if self.kdeg.get((nv - 1) as usize).copied().unwrap_or(0) % 2 != 0 {
                opposite = true;
            }
            out = self.oscul2.get((nv - 1) as usize).cloned();
            along = out.is_some();
        }
        (along, opposite, out)
    }
}

fn cache_at(cache: &[Vec<GpPnt>], i: i32, j: i32) -> GpPnt {
    cache
        .get((i - 1) as usize)
        .and_then(|r| r.get((j - 1) as usize))
        .copied()
        .unwrap_or(GpPnt::new(0.0, 0.0, 0.0))
}

/// Tolerance used by `Geom_OffsetSurface` (`cxx:265`).
pub const OSCULATING_TOL: f64 = CONFUSION;
