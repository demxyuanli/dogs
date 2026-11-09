//! The stop tests of `IntWalk_IWalking`.
//!
//! - `TestArretPassage` for an open line (`IntWalk_IWalking.gxx:584-766`).
//! - `TestArretPassage` for a closed line (`gxx:768-953`).
//! - `TestArretAjout` (`gxx:955-1021`).
//! - `TestArretCadre` (`gxx:1176-1388`).
//! - `TestPassedSolutionWithNegativeState` (`gxx:1392-1988`).
//!
//! OCCT returns the `Arrive` flag and writes `Irang` through a reference; the
//! port answers `Option<i32>`: `None` is `Arrive == false`, `Some(Irang)`
//! carries the index.

use occt_core::gp::GpVec2d;
use occt_math::{MathMatrix, MathVector};

use crate::int_tools_wline::PntOn2S;

use super::types::{cut_vector_by_tolerances, math_vec2, uv_on, WalkingData};
use super::{ContapIWLineOfTheIWalking, ContapTheIWalking};

/// Re-run `Values(UV, F, D)` so that `Root`, `Point`, `Direction3d` and
/// `Direction2d` are those of `uv` (the `sp.Values(...)` calls of the tests).
pub(super) fn refresh_function(
    sp: &mut super::super::surf_function::ContapSurfFunction,
    uv: &[f64; 2],
) {
    use occt_math::function_set_root::MathFunctionSetWithDerivatives;

    let mut f = MathVector::new(1, 1);
    let mut d = MathMatrix::new(1, 1, 1, 2);
    let _ = sp.values(&math_vec2(uv[0], uv[1]), &mut f, &mut d);
}

impl ContapTheIWalking {
    /// `TestArretPassage(Umult, Vmult, sp, UV, Irang)` (`gxx:584-766`), the
    /// open-line test.
    pub(super) fn test_arret_passage_open(
        &mut self,
        umult: &[f64],
        vmult: &[f64],
        sp: &mut super::super::surf_function::ContapSurfFunction,
        uv: &mut [f64; 2],
    ) -> Option<i32> {
        let tolu = self.vec_tolerance[0];
        let tolv = self.vec_tolerance[1];
        let tolu2 = 10.0 * self.vec_tolerance[0];
        let tolv2 = 10.0 * self.vec_tolerance[1];

        let mut arrive = false;
        let mut irang = 0i32;

        let (up, vp) = uv_on(&self.previous_point, self.reversed);

        // `gxx:625-663`: crossing test on the interior start points.
        for i in 1..self.wd2.len() {
            if self.wd2[i].etat > 0 {
                let utest = self.wd2[i].ustart;
                let vtest = self.wd2[i].vstart;

                let du = uv[0] - utest;
                let dv = uv[1] - vtest;
                let dup = up - utest;
                let dvp = vp - vtest;

                if (du.abs() < tolu2 && dv.abs() < tolv2)
                    || (dup.abs() < tolu2 && dvp.abs() < tolv2)
                {
                    self.wd2[i].etat = -self.wd2[i].etat;
                } else {
                    let d_du = uv[0] - up;
                    let d_dv = uv[1] - vp;
                    let ddd = d_du * d_du + d_dv * d_dv;
                    let dd1 = du * du + dv * dv;
                    if dd1 <= ddd {
                        let dd2 = dup * dup + dvp * dvp;
                        if dd2 <= ddd && (du * dup + dv * dvp * tolu / tolv) <= 0.0 {
                            self.wd2[i].etat = -self.wd2[i].etat;
                        }
                    }
                }
            }
        }

        // `gxx:666-765`: stop test on the points given at input.
        let mut i_candidates: Vec<usize> = Vec::new();
        let mut sq_dist_candidates: Vec<f64> = Vec::new();

        let mut l = 1;
        while l <= 2 && !arrive {
            for i in 1..self.wd1.len() {
                let is_to_check = if l == 1 {
                    self.wd1[i].etat > 0
                } else {
                    self.wd1[i].etat < 0
                };
                if !is_to_check {
                    continue;
                }

                let utest = self.wd1[i].ustart;
                let vtest = self.wd1[i].vstart;
                let dup = up - utest;
                let dvp = vp - vtest;
                if dup.abs() < tolu && dvp.abs() < tolv {
                    continue;
                }

                let uv1m_utest = uv[0] - utest;
                let uv2m_vtest = uv[1] - vtest;
                if (dup * uv1m_utest + dvp * uv2m_vtest) < 0.0
                    || (uv1m_utest.abs() < tolu && uv2m_vtest.abs() < tolv)
                {
                    i_candidates.push(i);
                    sq_dist_candidates.push(dup * dup + dvp * dvp);
                } else if self.nb_multiplicities[i] > 0 && i_candidates.is_empty() {
                    let mut n = 0i32;
                    for k in 1..i {
                        n += self.nb_multiplicities[k];
                    }
                    let last = n + self.nb_multiplicities[i];
                    let mut j = n + 1;
                    while j <= last {
                        let (umj, vmj) = (umult[(j - 1) as usize], vmult[(j - 1) as usize]);
                        if ((up - umj) * (uv[0] - umj) + (vp - vmj) * (uv[1] - vmj) < 0.0)
                            || (uv[0] - umj).abs() < tolu && (uv[1] - vmj).abs() < tolv
                        {
                            irang = i as i32;
                            arrive = true;
                            uv[0] = utest;
                            uv[1] = vtest;
                            break;
                        }
                        j += 1;
                    }
                }
                if arrive {
                    refresh_function(sp, uv);
                    break;
                }
            }

            if !i_candidates.is_empty() {
                let mut min_sq_dist = f64::MAX;
                for (ind, &sq) in sq_dist_candidates.iter().enumerate() {
                    if sq < min_sq_dist {
                        min_sq_dist = sq;
                        irang = i_candidates[ind] as i32;
                    }
                }
                arrive = true;
                uv[0] = self.wd1[irang as usize].ustart;
                uv[1] = self.wd1[irang as usize].vstart;
            }
            l += 1;
        }

        if arrive {
            Some(irang)
        } else {
            None
        }
    }

    /// `TestArretPassage(Umult, Vmult, UV, Index, Irang)` (`gxx:768-953`), the
    /// closed-line test.
    pub(super) fn test_arret_passage_closed(
        &mut self,
        umult: &[f64],
        vmult: &[f64],
        uv: &[f64; 2],
        index: i32,
    ) -> (bool, i32) {
        let mut arrive = false;
        let mut irang = 0i32;

        let mut tolu = self.vec_tolerance[0];
        let mut tolv = self.vec_tolerance[1];

        let (mut up, mut vp) = uv_on(&self.previous_point, self.reversed);

        let mut uv1 = uv[0];
        let mut uv2 = uv[1];

        // `gxx:809-828`: normalise the range so that a 1.0 span is the
        // reference length.
        let is_highly_anisotropic = tolu.max(tolv) > 1000.0 * tolu.min(tolv);
        let deltau = if self.s_range_u.is_void() {
            self.umax - self.um
        } else if is_highly_anisotropic {
            self.s_range_u.delta()
        } else {
            self.s_range_u.delta().max(1.0)
        };
        let deltav = if self.s_range_v.is_void() {
            self.vmax - self.vm
        } else if is_highly_anisotropic {
            self.s_range_v.delta()
        } else {
            self.s_range_v.delta().max(1.0)
        };

        up /= deltau;
        uv1 /= deltau;
        vp /= deltav;
        uv2 /= deltav;

        tolu /= deltau;
        tolv /= deltav;

        let tolu2 = tolu + tolu;
        let tolv2 = tolv + tolv;

        let d_previous_current = (up - uv1) * (up - uv1) + (vp - uv2) * (vp - uv2);
        for k in 1..self.wd2.len() {
            if self.wd2[k].etat <= 0 {
                continue;
            }
            let utest = self.wd2[k].ustart / deltau;
            let vtest = self.wd2[k].vstart / deltav;

            let uv1m_utest = uv1 - utest;
            let uv2m_vtest = uv2 - vtest;
            if (uv1m_utest < tolu2 && uv1m_utest > -tolu2)
                && (uv2m_vtest < tolv2 && uv2m_vtest > -tolv2)
            {
                if index != k as i32 {
                    self.wd2[k].etat = -self.wd2[k].etat;
                } else {
                    arrive = true;
                }
            } else {
                let upm_utest = up - utest;
                let vpm_vtest = vp - vtest;
                let d_previous_start = upm_utest * upm_utest + vpm_vtest * vpm_vtest;
                let d_current_start = uv1m_utest * uv1m_utest + uv2m_vtest * uv2m_vtest;

                let scal = upm_utest * uv1m_utest + vpm_vtest * uv2m_vtest;
                if upm_utest.abs() < tolu && vpm_vtest.abs() < tolv {
                    if index != k as i32 {
                        self.wd2[k].etat = -self.wd2[k].etat;
                    }
                } else if scal < 0.0 && (d_previous_start + d_current_start < d_previous_current) {
                    if index == k as i32 {
                        arrive = true;
                    } else {
                        self.wd2[k].etat = -self.wd2[k].etat;
                    }
                } else if k as i32 != index {
                    if d_previous_start < d_previous_current * 0.25 {
                        self.wd2[k].etat = -self.wd2[k].etat;
                    } else if d_current_start < d_previous_current * 0.25 {
                        self.wd2[k].etat = -self.wd2[k].etat;
                    } else {
                        let umid_utest = 0.5 * (uv1 + up) - utest;
                        let vmid_vtest = 0.5 * (uv2 + vp) - vtest;
                        let d_middle_start = umid_utest * umid_utest + vmid_vtest * vmid_vtest;
                        if d_middle_start < d_previous_current * 0.5 {
                            self.wd2[k].etat = -self.wd2[k].etat;
                        }
                    }
                }
            }
        }

        // `gxx:919-951`: crossing test on the passing start points.
        irang = 0;
        for i in 1..self.wd1.len() {
            if self.wd1[i].etat > 0 && self.wd1[i].etat < 11 {
                let utest = self.wd1[i].ustart / deltau;
                let vtest = self.wd1[i].vstart / deltav;

                if ((up - utest) * (uv1 - utest) + (vp - vtest) * (uv2 - vtest) < 0.0)
                    || (uv1 - utest).abs() < tolu && (uv2 - vtest).abs() < tolv
                {
                    irang = i as i32;
                } else if self.nb_multiplicities[i] > 0 {
                    let mut n = 0i32;
                    for k in 1..i {
                        n += self.nb_multiplicities[k];
                    }
                    let mut crossed = false;
                    let last = n + self.nb_multiplicities[i];
                    let mut j = n + 1;
                    while j <= last {
                        let umj = umult[(j - 1) as usize] / deltau;
                        let vmj = vmult[(j - 1) as usize] / deltav;
                        if ((up - umj) * (uv1 - umj) + (vp - vmj) * (uv2 - vmj) < 0.0)
                            || (uv1 - umj).abs() < tolu && (uv2 - vmj).abs() < tolv
                        {
                            irang = i as i32;
                            crossed = true;
                            break;
                        }
                        j += 1;
                    }
                    let _ = crossed;
                }
            }
        }

        (arrive, irang)
    }

    /// `TestArretAjout(sp, UV, Irang, Psol)` (`gxx:955-1021`).
    pub(super) fn test_arret_ajout(
        &mut self,
        sp: &mut super::super::surf_function::ContapSurfFunction,
        uv: &mut [f64; 2],
        psol: &mut PntOn2S,
    ) -> Option<i32> {
        let (up, vp) = uv_on(&self.previous_point, self.reversed);
        let mut arrive = false;
        let mut irang = 0i32;

        let seq: Vec<i32> = self.seq_ajout.clone();
        for irang_v in seq {
            irang = irang_v;
            if (irang.unsigned_abs() as usize) <= self.lines.len() {
                let line = &self.lines[irang.unsigned_abs() as usize - 1];
                let p = if irang > 0 {
                    *line.value(line.nb_points())
                } else {
                    *line.value(1)
                };
                *psol = p;
                let (u1, v1) = uv_on(&p, self.reversed);
                if ((up - u1) * (uv[0] - u1) + (vp - v1) * (uv[1] - v1)) < 0.0
                    || ((uv[0] - u1).abs() < self.vec_tolerance[0]
                        && (uv[1] - v1).abs() < self.vec_tolerance[1])
                {
                    arrive = true;
                    uv[0] = u1;
                    uv[1] = v1;
                    refresh_function(sp, uv);
                    break;
                }
            }
        }

        if arrive {
            Some(irang)
        } else {
            None
        }
    }

    /// `TestArretCadre(Umult, Vmult, Line, sp, UV, Irang)` (`gxx:1176-1388`).
    ///
    /// `Line` is cut when a stop point is found on it. `Some(Irang)` is the
    /// `Arrive` case; a negative `Irang` means the point must be added to the
    /// line without a `Cut`.
    pub(super) fn test_arret_cadre(
        &mut self,
        umult: &[f64],
        vmult: &[f64],
        line: &mut ContapIWLineOfTheIWalking,
        sp: &mut super::super::surf_function::ContapSurfFunction,
        uv: &mut [f64; 2],
    ) -> Option<i32> {
        let mut irang = 0i32;
        let mut found = false;

        for i in 1..self.wd1.len() {
            if self.wd1[i].etat >= 0 {
                continue;
            }
            let mut n = 0i32;
            if self.nb_multiplicities[i] > 0 {
                for k in 1..i {
                    n += self.nb_multiplicities[k];
                }
            }

            let (mut up, mut vp) = uv_on(line.value(1), self.reversed);
            let mut nbp = line.nb_points();
            let mut j = 2;
            while j <= nbp {
                let (uc, vc) = uv_on(line.value(j), self.reversed);

                let mut a_vec1 = GpVec2d::new(up - self.wd1[i].ustart, vp - self.wd1[i].vstart);
                let mut a_vec2 = GpVec2d::new(uc - self.wd1[i].ustart, vc - self.wd1[i].vstart);
                cut_vector_by_tolerances(&mut a_vec1, self.vec_tolerance);
                cut_vector_by_tolerances(&mut a_vec2, self.vec_tolerance);

                let scal = a_vec1.dot(&a_vec2);

                if scal < 0.0 {
                    line.cut(j);
                    nbp = line.nb_points();
                    irang = i as i32;
                    uv[0] = self.wd1[irang as usize].ustart;
                    uv[1] = self.wd1[irang as usize].vstart;
                    found = true;
                } else if (uc - self.wd1[i].ustart).abs() < self.vec_tolerance[0]
                    && (vc - self.wd1[i].vstart).abs() < self.vec_tolerance[1]
                {
                    line.cut(j);
                    nbp = line.nb_points();
                    irang = i as i32;
                    uv[0] = self.wd1[irang as usize].ustart;
                    uv[1] = self.wd1[irang as usize].vstart;
                    found = true;
                } else if self.nb_multiplicities[i] > 0 {
                    let mut k = n + 1;
                    while k <= n + self.nb_multiplicities[i] {
                        let (umk, vmk) = (umult[(k - 1) as usize], vmult[(k - 1) as usize]);
                        let mut a_vec1 = GpVec2d::new(up - umk, vp - vmk);
                        let mut a_vec2 = GpVec2d::new(uc - umk, vc - vmk);
                        cut_vector_by_tolerances(&mut a_vec1, self.vec_tolerance);
                        cut_vector_by_tolerances(&mut a_vec2, self.vec_tolerance);

                        let scal = a_vec1.dot(&a_vec2);
                        if scal < 0.0 {
                            line.cut(j);
                            nbp = line.nb_points();
                            irang = i as i32;
                            uv[0] = self.wd1[irang as usize].ustart;
                            uv[1] = self.wd1[irang as usize].vstart;
                            found = true;
                            break;
                        } else if (uc - umk).abs() < self.vec_tolerance[0]
                            && (vc - vmk).abs() < self.vec_tolerance[1]
                        {
                            line.cut(j);
                            nbp = line.nb_points();
                            irang = i as i32;
                            uv[0] = self.wd1[irang as usize].ustart;
                            uv[1] = self.wd1[irang as usize].vstart;
                            found = true;
                            break;
                        }
                        k += 1;
                    }
                }

                if found {
                    refresh_function(sp, uv);
                    let nbp_line = line.nb_points();
                    let (_, indextg) = line.tangent_vector();
                    if indextg > nbp_line && j > 3 && j <= nbp_line + 1 {
                        // `gxx:1297-1307`: `Dir3d` is a local copy, so its
                        // `Reverse()` never reaches the line; OCCT stores
                        // `previousd3d` in both branches.
                        line.set_tangent_vector(self.previous_d3d, j - 1);
                    }
                    return Some(irang);
                }
                up = uc;
                vp = vc;
                j += 1;
            }

            // `gxx:1322-1385`: last point of the line against the current one.
            let mut a_vec1 = GpVec2d::new(up - self.wd1[i].ustart, vp - self.wd1[i].vstart);
            let mut a_vec2 = GpVec2d::new(uv[0] - self.wd1[i].ustart, uv[1] - self.wd1[i].vstart);
            cut_vector_by_tolerances(&mut a_vec1, self.vec_tolerance);
            cut_vector_by_tolerances(&mut a_vec2, self.vec_tolerance);

            let scal = a_vec1.dot(&a_vec2);
            if scal < 0.0 {
                irang = i as i32;
                uv[0] = self.wd1[irang as usize].ustart;
                uv[1] = self.wd1[irang as usize].vstart;
                found = true;
            } else if (uv[0] - self.wd1[i].ustart).abs() < self.vec_tolerance[0]
                && (uv[1] - self.wd1[i].vstart).abs() < self.vec_tolerance[1]
            {
                irang = i as i32;
                uv[0] = self.wd1[irang as usize].ustart;
                uv[1] = self.wd1[irang as usize].vstart;
                found = true;
            } else if self.nb_multiplicities[i] > 0 {
                let mut j = n + 1;
                while j <= n + self.nb_multiplicities[i] {
                    let (umj, vmj) = (umult[(j - 1) as usize], vmult[(j - 1) as usize]);
                    let mut a_vec1 = GpVec2d::new(up - umj, vp - vmj);
                    let mut a_vec2 = GpVec2d::new(uv[0] - umj, uv[1] - vmj);
                    cut_vector_by_tolerances(&mut a_vec1, self.vec_tolerance);
                    cut_vector_by_tolerances(&mut a_vec2, self.vec_tolerance);

                    let scal = a_vec1.dot(&a_vec2);
                    if scal < 0.0 {
                        irang = i as i32;
                        uv[0] = self.wd1[irang as usize].ustart;
                        uv[1] = self.wd1[irang as usize].vstart;
                        found = true;
                        break;
                    } else if (uv[0] - umj).abs() < self.vec_tolerance[0]
                        && (uv[1] - vmj).abs() < self.vec_tolerance[1]
                    {
                        irang = i as i32;
                        uv[0] = self.wd1[irang as usize].ustart;
                        uv[1] = self.wd1[irang as usize].vstart;
                        found = true;
                        break;
                    }
                    j += 1;
                }
            }
            if found {
                irang = -irang; // `gxx:1382`
                refresh_function(sp, uv);
                return Some(irang);
            }
        }
        None
    }
}

/// `TestPassedSolutionWithNegativeState(wd, Umult, Vmult, prevUp, prevVp,
/// nbMultiplicities, tolerance, sp, UV, Irang)` (`gxx:1912-1988`).
pub(super) fn test_passed_solution_with_negative_state(
    wd: &[WalkingData],
    umult: &[f64],
    vmult: &[f64],
    prev_up: f64,
    prev_vp: f64,
    nb_multiplicities: &[i32],
    tolerance: [f64; 2],
    sp: &mut super::super::surf_function::ContapSurfFunction,
    uv: &mut [f64; 2],
) -> Option<i32> {
    let mut arrive = false;
    let mut irang = 0i32;
    let tolu = tolerance[0];
    let tolv = tolerance[1];

    for i in 1..wd.len() {
        if wd[i].etat >= -11 {
            continue;
        }
        let utest = wd[i].ustart;
        let vtest = wd[i].vstart;
        let dup = prev_up - utest;
        let dvp = prev_vp - vtest;
        if dup.abs() < tolu && dvp.abs() < tolv {
            continue;
        }

        let uv1m_utest = uv[0] - utest;
        let uv2m_vtest = uv[1] - vtest;
        if (dup * uv1m_utest + dvp * uv2m_vtest) < 0.0
            || (uv1m_utest.abs() < tolu && uv2m_vtest.abs() < tolv)
        {
            irang = i as i32;
            arrive = true;
            uv[0] = utest;
            uv[1] = vtest;
        } else if nb_multiplicities[i] > 0 {
            let mut n = 0i32;
            for k in 1..i {
                n += nb_multiplicities[k];
            }
            let mut j = n + 1;
            while j <= n + nb_multiplicities[i] {
                let (umj, vmj) = (umult[(j - 1) as usize], vmult[(j - 1) as usize]);
                if ((prev_up - umj) * (uv[0] - umj) + (prev_vp - vmj) * (uv[1] - vmj) < 0.0)
                    || (uv[0] - umj).abs() < tolu && (uv[1] - vmj).abs() < tolv
                {
                    irang = i as i32;
                    arrive = true;
                    uv[0] = utest;
                    uv[1] = vtest;
                    break;
                }
                j += 1;
            }
        }
        if arrive {
            refresh_function(sp, uv);
            break;
        }
    }

    if arrive {
        Some(irang)
    } else {
        None
    }
}
