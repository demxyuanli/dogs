//! `ProjLib_CompProjectedCurve::Init` and its sequence post-processing.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx:657-1227`
//! (OCCT 8.0.0).
//!
//! `Init` marches along the curve parameter `t` with an adaptive step
//! (`WalkStep = clamp(0.1 * |V1| / |V2|, MinStep, MaxStep)`), starting a new
//! continuous part wherever `ProjLib_PrjResolve` stops being solvable. Each
//! part's first point is refined onto the surface border (`ExactBound`, with
//! `DichExactBound` as the fallback) and crossing a split point found by
//! [`super::splits::build_curve_splits`] starts a new part too.
//!
//! Translation notes:
//! - `Nend` is only assigned inside the `myMaxDist > 0 && type !=
//!   SurfaceOfExtrusion` block (`cxx:676-692`), so with `myMaxDist <= 0` the
//!   extrema-based initial point search is skipped entirely (`cxx:744`). The
//!   port keeps that structure verbatim.
//! - `Extrema_ExtCS::IsParallel()` (`cxx:744`) is UNPORTED (the port's
//!   curve-surface extrema do not carry the flag, see
//!   `extrema_surf::curve_surface_extrema_all`); the port takes the
//!   non-parallel arm.
//! - `mySequence` indices are 1-based in OCCT; `[i - 1]` here.

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::SurfaceType;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::approx_same_parameter::{u_resolution, v_resolution};
use crate::extrema_surf::curve_surface_extrema_all;
use crate::projlib::projected_curve::surface_type;
use crate::projlib::prj_resolve::PrjResolve;

use super::bounds::{dich_exact_bound, exact_bound, initial_point};
use super::jets::{d1, d2_curv_on_surf};
use super::splits::build_curve_splits;
use super::{CompProjectedCurve, FUNC_TOL};

/// `GlobalMinStep` (`cxx:700`): "sufficiently small to provide solving from
/// an initial point, and sufficiently large to avoid too close solutions".
const GLOBAL_MIN_STEP: f64 = 1.0e-4;

impl CompProjectedCurve {
    /// `Init()` (`cxx:657-1227`).
    pub(crate) fn init(&mut self) {
        self.tab_int = None;
        let mut a_splits: Vec<f64> = Vec::new();
        let mut a_split_idx = 0usize;
        let mut is_splits_computed = false;

        let (Some(surface), Some(curve)) = (self.surface.clone(), self.curve.clone()) else {
            return;
        };
        let s = surface.as_ref();
        let c = curve.as_ref();

        // `Extrema_ExtCS CExt(*myCurve, *mySurface, aTolExt, aTolExt)`
        // (`cxx:669`). `CExt.IsDone() && CExt.NbExt()` is a non-empty result.
        let cext = curve_surface_extrema_all(c, s);
        let mut nend = 0usize;
        if !cext.is_empty() {
            // Search for the minimum solution. Avoid usage of extrema result
            // that can be wrong for extrusion (`cxx:673-693`).
            if self.max_dist > 0.0 && surface_type(s) != SurfaceType::SurfaceOfExtrusion {
                let mut min_val2 = cext[0].distance * cext[0].distance;
                nend = cext.len();
                for e in cext.iter().skip(1) {
                    let d2 = e.distance * e.distance;
                    if d2 < min_val2 {
                        min_val2 = d2;
                    }
                }
                if min_val2 > self.max_dist * self.max_dist {
                    return; // No near solution -> exit.
                }
            }
        }

        let first_u = c.first_parameter();
        let last_u = c.last_parameter();
        let min_step = 0.01 * (last_u - first_u);
        let max_step = 0.1 * (last_u - first_u);
        let search_step = 10.0 * min_step;
        let mut step = search_step;

        let a_low_border = (s.u_range().0, s.v_range().0);
        let a_upp_border = (s.u_range().1, s.v_range().1);
        let a_tol = (self.tol_u, self.tol_v);
        let mut a_prj_ps = PrjResolve::new(c, s, 1).expect("myFix in 1..=3");

        let mut t = first_u;
        let mut prev_deb = 0.0;
        let mut same_deb = false;
        let mut from_last_u = false;
        let mut triple = GpPnt::zero();
        let mut prev_triple;

        // Basic loop (`cxx:720`).
        while t <= last_u {
            // Search for the beginning of a new continuous part, and avoid
            // infinite computation in some difficult cases (`cxx:723-728`).
            let mut new_part = false;
            if t > first_u && (t - prev_deb).abs() <= PCONFUSION {
                same_deb = true;
            }
            while t <= last_u && !new_part && !from_last_u && !same_deb {
                prev_deb = t;
                if t == last_u {
                    from_last_u = true;
                }
                let mut initpoint = false;
                let mut u = 0.0f64;
                let mut v = 0.0f64;
                let mut cpoint = GpPnt::zero();

                // Search an initial point in the list of Extrema Curve-Surface
                // (`cxx:744-765`). `Extrema_ExtCS::IsParallel` is UNPORTED.
                if nend != 0 {
                    for e in cext.iter().take(nend) {
                        let par_t = e.u1;
                        let par_u = e.u2;
                        let par_v = e.v2.unwrap_or(0.0);
                        a_prj_ps.perform(
                            par_t,
                            par_u,
                            par_v,
                            a_tol,
                            a_low_border,
                            a_upp_border,
                            FUNC_TOL,
                            true,
                        );
                        if a_prj_ps.is_done()
                            && par_t > first_u.max(t - step + PCONFUSION)
                            && par_t <= t
                        {
                            t = par_t;
                            u = par_u;
                            v = par_v;
                            cpoint = e.p1;
                            initpoint = true;
                            break;
                        }
                    }
                }
                if !initpoint {
                    cpoint = c.d0(t);
                    // PConfusion - use geometric tolerances in extrema /
                    // optimization (`cxx:772-779`).
                    if let Some((uu, vv)) = initial_point(
                        &cpoint, t, c, s, self.tol_u, self.tol_v, self.max_dist,
                    ) {
                        u = uu;
                        v = vv;
                        initpoint = true;
                    }
                }
                if initpoint {
                    // When U or V lie on a surface joint, in some cases we
                    // cannot use them as the initial point for aPrjPS, so we
                    // switch them (`cxx:782-830`).
                    if (s.is_u_periodic()
                        && (a_upp_border.0 - a_low_border.0 - s.u_period()).abs() < CONFUSION)
                        || (s.is_v_periodic()
                            && (a_upp_border.1 - a_low_border.1 - s.v_period()).abs() < CONFUSION)
                    {
                        if (u - a_low_border.0).abs() < u_resolution(s, PCONFUSION)
                            && s.is_u_periodic()
                        {
                            let d = d1(t, u, v, c, s).expect("d1");
                            if d.x() < 0.0 {
                                u = a_upp_border.0;
                            }
                        } else if (u - a_upp_border.0).abs() < u_resolution(s, PCONFUSION)
                            && s.is_u_periodic()
                        {
                            let d = d1(t, u, v, c, s).expect("d1");
                            if d.x() > 0.0 {
                                u = a_low_border.0;
                            }
                        }

                        if (v - a_low_border.1).abs() < v_resolution(s, PCONFUSION)
                            && s.is_v_periodic()
                        {
                            let d = d1(t, u, v, c, s).expect("d1");
                            if d.y() < 0.0 {
                                v = a_upp_border.1;
                            }
                        } else if (v - a_upp_border.1).abs() <= v_resolution(s, PCONFUSION)
                            && s.is_v_periodic()
                        {
                            let d = d1(t, u, v, c, s).expect("d1");
                            if d.y() > 0.0 {
                                v = a_low_border.1;
                            }
                        }
                    }

                    if self.max_dist > 0.0 {
                        // Stop if the distance between the projection and the
                        // corresponding curve point is greater than MaxDist
                        // (`cxx:832-847`).
                        let p_on_s = s.d0(u, v);
                        if cpoint.distance(&p_on_s) > self.max_dist {
                            self.sequence.clear();
                            self.nb_curves = 0;
                            return;
                        }
                    }

                    triple = GpPnt::new(t, u, v);
                    if t != first_u {
                        // Search for an exact boundary point (`cxx:849-864`).
                        let mut tol = self.tol_u.min(self.tol_v);
                        let ad = d1(triple.x(), triple.y(), triple.z(), c, s).expect("d1");
                        tol /= ad.x().abs().max(ad.y().abs());
                        if !exact_bound(
                            &mut triple,
                            t - step,
                            tol,
                            self.tol_u,
                            self.tol_v,
                            c,
                            s,
                        ) {
                            dich_exact_bound(
                                &mut triple,
                                t - step,
                                tol,
                                self.tol_u,
                                self.tol_v,
                                c,
                                s,
                            );
                        }
                    }
                    new_part = true;
                } else {
                    if t == last_u {
                        break;
                    }
                    t += step;
                    if t > last_u {
                        step = step + last_u - t;
                        t = last_u;
                    }
                }
            }
            if !new_part {
                break;
            }

            // We have found a new continuous part (`cxx:886-890`).
            self.sequence.push(vec![triple]);
            self.nb_curves += 1;
            prev_triple = triple;

            if triple.x() == last_u {
                break;
            }

            // Computation of WalkStep (`cxx:896-911`).
            let (d1v, d2v) =
                d2_curv_on_surf(triple.x(), triple.y(), triple.z(), c, s).expect("d2CurvOnSurf");
            let magn_d1 = d1v.magnitude();
            let magn_d2 = d2v.magnitude();
            let walk_step = if magn_d2 < CONFUSION {
                max_step
            } else {
                max_step.min(min_step.max(0.1 * magn_d1 / magn_d2))
            };
            step = walk_step;

            t = triple.x() + step;
            if t > last_u {
                t = last_u;
            }
            let mut prev_step = step;

            // Here we are trying to prolong the continuous part
            // (`cxx:922-1081`).
            while t <= last_u && new_part {
                let mut u0 = triple.y() + (step / prev_step) * (triple.y() - prev_triple.y());
                let mut v0 = triple.z() + (step / prev_step) * (triple.z() - prev_triple.z());
                u0 = u0.max(a_low_border.0).min(a_upp_border.0);
                v0 = v0.max(a_low_border.1).min(a_upp_border.1);

                a_prj_ps.perform(
                    t,
                    u0,
                    v0,
                    a_tol,
                    a_low_border,
                    a_upp_border,
                    FUNC_TOL,
                    true,
                );
                if !a_prj_ps.is_done() {
                    if step <= GLOBAL_MIN_STEP {
                        // Search for an exact boundary point (`cxx:939-973`).
                        let mut tol = self.tol_u.min(self.tol_v);
                        let d = d1(triple.x(), triple.y(), triple.z(), c, s).expect("d1");
                        tol /= d.x().abs().max(d.y().abs());
                        if !exact_bound(&mut triple, t, tol, self.tol_u, self.tol_v, c, s) {
                            dich_exact_bound(&mut triple, t, tol, self.tol_u, self.tol_v, c, s);
                        }
                        let part_idx = self.nb_curves as usize - 1;
                        let last_x = self.sequence[part_idx][self.sequence[part_idx].len() - 1].x();
                        if triple.x() - last_x > 1.0e-10 {
                            self.sequence[part_idx].push(triple);
                        }
                        if (last_u - triple.x()) < tol {
                            t = last_u + 1.0;
                            break;
                        }
                        step = search_step;
                        t = triple.x() + step;
                        if t > (last_u - min_step / 2.0) {
                            step = step + last_u - t;
                            t = last_u;
                        }
                        new_part = false;
                    } else {
                        // Decrease the step (`cxx:976-988`).
                        let save_step = step;
                        step /= 2.0;
                        t = triple.x() + step;
                        if t > (last_u - min_step / 4.0) {
                            step = step + last_u - t;
                            if (step - save_step).abs() <= PCONFUSION {
                                step = GLOBAL_MIN_STEP; // to avoid looping
                            }
                            t = last_u;
                        }
                    }
                } else {
                    // Go further (`cxx:991-1080`).
                    prev_triple = triple;
                    prev_step = step;
                    let sol = a_prj_ps.solution().expect("IsDone");
                    triple = GpPnt::new(t, sol.0, sol.1);

                    // Check for possible local traps.
                    self.update_triple_by_trap_criteria(&mut triple);

                    // Protection from the case when the whole curve lies on a
                    // seam (`cxx:1001-1023`).
                    if !is_splits_computed {
                        let mut is_u_possible = false;
                        if s.is_u_periodic()
                            && (triple.y() - a_low_border.0).abs() > PCONFUSION
                            && (triple.y() - a_upp_border.0).abs() > PCONFUSION
                        {
                            is_u_possible = true;
                        }
                        let mut is_v_possible = false;
                        if s.is_v_periodic()
                            && (triple.z() - a_low_border.1).abs() > PCONFUSION
                            && (triple.z() - a_upp_border.1).abs() > PCONFUSION
                        {
                            is_v_possible = true;
                        }
                        if is_u_possible || is_v_possible {
                            // When the point is well conditioned.
                            build_curve_splits(c, &surface, self.tol_u, self.tol_v, &mut a_splits);
                            is_splits_computed = true;
                        }
                    }

                    let part_idx = self.nb_curves as usize - 1;
                    let last_x = self.sequence[part_idx][self.sequence[part_idx].len() - 1].x();
                    if triple.x() - last_x > 1.0e-10 {
                        self.sequence[part_idx].push(triple);
                    }
                    if t == last_u {
                        t = last_u + 1.0;
                        break;
                    }

                    // Computation of WalkStep (`cxx:1038-1053`).
                    let (d1v, d2v) =
                        d2_curv_on_surf(triple.x(), triple.y(), triple.z(), c, s)
                            .expect("d2CurvOnSurf");
                    let magn_d1 = d1v.magnitude();
                    let magn_d2 = d2v.magnitude();
                    step = if magn_d2 < CONFUSION {
                        max_step
                    } else {
                        max_step.min(min_step.max(0.1 * magn_d1 / magn_d2))
                    };

                    t += step;
                    if t > (last_u - min_step / 2.0) {
                        step = step + last_u - t;
                        t = last_u;
                    }

                    // We assume at least one point of cache inside of a split
                    // (`cxx:1059-1079`).
                    let a_size = a_splits.len();
                    let mut an_idx = a_split_idx;
                    while an_idx < a_size {
                        let a_param = a_splits[an_idx];
                        if (a_param - triple.x()).abs() < PCONFUSION {
                            // The current point is equal to a split point.
                            new_part = false;
                            // Move the split index to avoid checking the whole
                            // list.
                            a_split_idx += 1;
                            break;
                        } else if a_param < t + PCONFUSION {
                            // The next point crosses the split point.
                            t = a_param;
                            step = t - prev_triple.x();
                        }
                        an_idx += 1;
                    }
                }
            }
        }

        self.post_process();
    }

    /// Sequence post-proceeding (`cxx:1086-1226`).
    fn post_process(&mut self) {
        // 1. Removing poor parts (`cxx:1089-1106`).
        let nb_part = self.nb_curves;
        let mut ipart = 0usize;
        let mut i = 1;
        while i <= nb_part {
            if self.sequence[ipart].len() < 2 {
                self.sequence.remove(ipart);
                self.nb_curves -= 1;
            } else {
                ipart += 1;
            }
            i += 1;
        }

        if self.nb_curves == 0 {
            return;
        }

        // 2. Removing common parts of bounds (`cxx:1110-1118`).
        for i in 1..self.nb_curves {
            let last_x = {
                let prev = &self.sequence[i as usize - 1];
                prev[prev.len() - 1].x()
            };
            if last_x >= self.sequence[i as usize][0].x() {
                self.sequence[i as usize][0].set_x(last_x + 1.0e-12);
            }
        }

        // 3. Computation of the maximum distance from each part of the curve
        // to the surface (`cxx:1122-1139`).
        self.max_distance = vec![0.0; self.nb_curves as usize];
        let (Some(surface), Some(curve)) = (self.surface.clone(), self.curve.clone()) else {
            return;
        };
        let s = surface.as_ref();
        let c = curve.as_ref();
        for i in 1..=self.nb_curves as usize {
            for tri in &self.sequence[i - 1] {
                let p_on_c = c.d0(tri.x());
                let p_on_s = s.d0(tri.y(), tri.z());
                let distance = p_on_c.distance(&p_on_s);
                if self.max_distance[i - 1] < distance {
                    self.max_distance[i - 1] = distance;
                }
            }
        }

        // 4. Check the projection to be a single point (`cxx:1143-1170`).
        self.sngl_pnts = vec![true; self.nb_curves as usize];
        let a_max_tol = if self.tol_u < self.tol_v {
            self.tol_v
        } else {
            self.tol_u
        };
        for i in 1..=self.nb_curves as usize {
            let part = &self.sequence[i - 1];
            let n = part.len() as f64;
            let mut ave_u = 0.0;
            let mut ave_v = 0.0;
            for tri in part {
                ave_u += tri.y();
                ave_v += tri.z();
            }
            ave_u /= n;
            ave_v /= n;
            for tri in part {
                let du = tri.y() - ave_u;
                let dv = tri.z() - ave_v;
                if (du * du + dv * dv).sqrt() > a_max_tol {
                    self.sngl_pnts[i - 1] = false;
                    break;
                }
            }
        }

        // 5. Check the projection to be an isoparametric curve of the surface
        // (`cxx:1174-1225`).
        self.ui_so = vec![true; self.nb_curves as usize];
        self.v_iso = vec![true; self.nb_curves as usize];
        for i in 1..=self.nb_curves as usize {
            let is_sngl = self.sngl_pnts[i - 1];
            let len = self.sequence[i - 1].len();
            if is_sngl || len <= 2 {
                self.ui_so[i - 1] = false;
                self.v_iso[i - 1] = false;
                continue;
            }
            let part = &self.sequence[i - 1];
            let n = len as f64;
            let mut ave_u = 0.0;
            let mut ave_v = 0.0;
            for tri in part {
                ave_u += tri.y();
                ave_v += tri.z();
            }
            ave_u /= n;
            ave_v /= n;

            // is i-part U-isoparametric ?
            let mut is_u = true;
            for tri in part {
                if (tri.y() - ave_u).abs() > self.tol_u {
                    is_u = false;
                    break;
                }
            }
            self.ui_so[i - 1] = is_u;

            // is i-part V-isoparametric ?
            let mut is_v = true;
            for tri in part {
                if (tri.z() - ave_v).abs() > self.tol_v {
                    is_v = false;
                    break;
                }
            }
            self.v_iso[i - 1] = is_v;
        }
    }
}
