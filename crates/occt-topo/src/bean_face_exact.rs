//! Exact BeanFace paths: HInter, Extrema, near-boundary completion.
use occt_core::gp::GpVec;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::extrema_surf::curve_surface_extrema_all;

use crate::bean_face::{golden_1d, BeanFaceIntersector};
use crate::bean_face_analytic::adjust_periodic;
use crate::bean_face_kind::{
    classify_fast_curve, plane_distance, plane_geometry, FastCurveKind,
};
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::intcurvesurface::perform_curve_surface;

const PI: f64 = std::f64::consts::PI;

impl BeanFaceIntersector {
    /// Exact curve–surface intersection via `IntCurveSurface_HInter`, then
    /// range expansion around every intersection point / segment. Port of
    /// `ComputeAroundExactIntersection`.
    pub(crate) fn compute_around_exact_intersection(&mut self) {
        let uv = self.finite_uv_bounds();
        let res = perform_curve_surface(
            self.curve(),
            self.surface(),
            (self.first_parameter, self.last_parameter),
            uv,
        );
        let hr = match res {
            Ok(h) => h,
            Err(_) => return,
        };

        // With more than one point, tighten the criteria to avoid merging
        // distinct crossings into one range.
        if hr.nb_points() > 1 {
            self.criteria = 3.0 * CONFUSION;
            self.curve_resolution = self.resolution(self.criteria);
        }

        for i in 0..hr.nb_points() {
            let p = hr.point(i);
            let mut w = p.param();
            if w < self.first_parameter || w > self.last_parameter {
                continue;
            }
            // Refine to the exact crossing (HInter's general path is only
            // accurate to ~1e-6, coarser than the result criteria).
            w = self.refine_crossing_param(w);
            let (su, sv) = {
                let pnt = self.curve_d0(w);
                let (a, b, _) = self.closest_params_dist(&pnt);
                (a, b)
            };
            let mut u = su;
            let mut v = sv;
            let u_not_valid = self.umin > u || u > self.umax;
            let v_not_valid = self.vmin > v || v > self.vmax;
            let mut solution_is_valid = !u_not_valid && !v_not_valid;

            if u_not_valid || v_not_valid {
                let mut b_u_corrected = true;
                if u_not_valid {
                    b_u_corrected = false;
                    solution_is_valid = false;
                    if self.surface().is_u_periodic() {
                        u = adjust_periodic(u, self.umin, self.umax, 2.0 * PI);
                        solution_is_valid = true;
                        b_u_corrected = true;
                    }
                }
                if b_u_corrected && v_not_valid {
                    solution_is_valid = false;
                    if self.surface().is_v_periodic() {
                        v = adjust_periodic(v, self.vmin, self.vmax, 2.0 * PI);
                        solution_is_valid = true;
                    }
                }
            }

            if !solution_is_valid {
                continue;
            }

            let n = self.range_manager.len();
            self.compute_range_from_start_point(false, w, u, v);
            self.compute_range_from_start_point(true, w, u, v);
            if n == self.range_manager.len() {
                self.set_empty_result_range(w);
            } else {
                self.min_sq_distance = 0.0;
            }
        }

        for i in 0..hr.nb_segments() {
            let seg = hr.segment(i);
            let p1 = seg.first_point();
            let p2 = seg.second_point();
            let first_param = if p1.param() < self.first_parameter {
                self.first_parameter
            } else {
                p1.param()
            };
            let last_param = if self.last_parameter < p2.param() {
                self.last_parameter
            } else {
                p2.param()
            };
            self.range_manager.insert_range(first_param, last_param, 2);
            self.compute_range_from_start_point(false, p1.param(), p1.u(), p1.v());
            self.compute_range_from_start_point(true, p2.param(), p2.u(), p2.v());
            self.min_sq_distance = 0.0;
        }
    }

    /// Complete result ranges whose start/end boundaries are near the surface
    /// but were missed by the discrete intersection points. Port of
    /// `ComputeNearRangeBoundaries`.
    pub(crate) fn compute_near_range_boundaries(&mut self) {
        let mut u = self.umin;
        let mut v = self.vmin;

        let n = self.range_manager.len();
        for i in 0..n {
            if self.range_manager.flag(i) > 0 {
                continue;
            }
            if i > 0 && self.range_manager.flag(i - 1) > 0 {
                continue;
            }
            let r = self.range_manager.range(i);
            if self.distance_with_uv(r.first, &mut u, &mut v) < self.criteria {
                let old_len = self.range_manager.len();
                if i > 0 {
                    self.compute_range_from_start_point_idx(false, r.first, u, v, i - 1);
                }
                let idx = i + (self.range_manager.len() - old_len);
                if idx < self.range_manager.len() {
                    self.compute_range_from_start_point_idx(true, r.first, u, v, idx);
                }
                if old_len == self.range_manager.len() {
                    self.set_empty_result_range(r.first);
                }
            }
        }

        if self.range_manager.is_empty() {
            return;
        }
        let last_idx = self.range_manager.len() - 1;
        if self.range_manager.flag(last_idx) == 0 {
            let r = self.range_manager.range(last_idx);
            if self.distance_with_uv(r.last, &mut u, &mut v) < self.criteria {
                let old_len = self.range_manager.len();
                self.compute_range_from_start_point_idx(false, r.last, u, v, last_idx);
                if old_len == self.range_manager.len() {
                    self.set_empty_result_range(r.last);
                }
            }
        }
    }

    /// Refine an HInter crossing parameter to the local minimum of the exact
    /// surface distance. The HInter general path converges with the coarse
    /// projector (`~1e-6` distance error), which can exceed the tightened
    /// result criteria (`3·Confusion`); refining makes the walk's starting
    /// point lie genuinely on the surface.
    pub(crate) fn refine_crossing_param(&self, w: f64) -> f64 {
        let window = (10.0 * self.curve_resolution).max(1e-5);
        let f = |u: f64| self.closest_params_dist(&self.curve_d0(u)).2;
        golden_1d(&f, w - window, w + window, 1e-10).0
    }

    /// Use curve–surface extrema to find near-surface spans that HInter's
    /// discrete points missed (tangencies / parallel spans). Port of
    /// `ComputeUsingExtremum` including the `Extrema_ExtCS::IsParallel` branch.
    pub(crate) fn compute_using_extremum(&mut self) {
        let tol = PCONFUSION;
        let mut i = 0usize;
        while i < self.range_manager.len() {
            if self.range_manager.flag(i) > 0 {
                i += 1;
                continue;
            }
            let r = self.range_manager.range(i);
            let anarg1 = r.first;
            let anarg2 = r.last;

            if anarg2 - anarg1 < PCONFUSION {
                if (i > 0 && self.range_manager.flag(i - 1) == 2)
                    || (i + 1 < self.range_manager.len() && self.range_manager.flag(i + 1) == 2)
                {
                    self.range_manager.set_flag(i, 1);
                    i += 1;
                    continue;
                }
            }

            let old_len = self.range_manager.len();

            if let Some(a_sq_dist) = self.extcs_is_parallel(anarg1, anarg2) {
                self.min_sq_distance = self.min_sq_distance.min(a_sq_dist);
                if a_sq_dist < self.criteria * self.criteria {
                    let mut u1 = 0.0;
                    let mut v1 = 0.0;
                    let mut u2 = 0.0;
                    let mut v2 = 0.0;
                    let adistance1 = self.distance_with_uv(anarg1, &mut u1, &mut v1);
                    let adistance2 = self.distance_with_uv(anarg2, &mut u2, &mut v2);
                    let valid1 = adistance1 < self.criteria;
                    let valid2 = adistance2 < self.criteria;
                    if valid1 && valid2 {
                        let _ = self.range_manager.insert_range(anarg1, anarg2, 2);
                    } else if valid1 {
                        self.compute_range_from_start_point(true, anarg1, u1, v1);
                    } else if valid2 {
                        self.compute_range_from_start_point(false, anarg2, u2, v2);
                    } else {
                        let mut a = anarg1;
                        let mut b = anarg2;
                        let mut da = adistance1;
                        let mut db = adistance2;
                        let mut found = false;
                        let mut asolution = a;
                        let mut u_sol = u1;
                        let mut v_sol = v1;
                        while (b - a) > self.curve_resolution && !found {
                            asolution = (a + b) * 0.5;
                            let adist = self.distance_with_uv(asolution, &mut u_sol, &mut v_sol);
                            if adist < self.criteria {
                                found = true;
                            } else if da < db {
                                b = asolution;
                                db = adist;
                            } else {
                                a = asolution;
                                da = adist;
                            }
                        }
                        if found {
                            self.compute_range_from_start_point(false, asolution, u_sol, v_sol);
                            self.compute_range_from_start_point(true, asolution, u_sol, v_sol);
                        } else {
                            self.range_manager.set_flag(i, 1);
                        }
                    }
                } else {
                    self.range_manager.set_flag(i, 1);
                }
                let diff = self.range_manager.len() - old_len;
                if diff > 0 {
                    i += diff;
                } else {
                    i += 1;
                }
                continue;
            }

            let mut solution_found = false;
            let exts = curve_surface_extrema_all(self.curve(), self.surface());
            let mut candidates: Vec<(f64, f64, f64)> = Vec::new();
            for e in &exts {
                self.min_sq_distance = self.min_sq_distance.min(e.distance * e.distance);
                if e.distance * e.distance >= self.criteria * self.criteria {
                    continue;
                }
                if e.u1 < anarg1 - tol || e.u1 > anarg2 + tol {
                    continue;
                }
                let u = e.u2;
                let v = e.v2.unwrap_or(0.0);
                if u < self.umin || u > self.umax || v < self.vmin || v > self.vmax {
                    continue;
                }
                candidates.push((e.u1, u, v));
            }

            if candidates.is_empty() {
                let e = occt_geom::extrema::curve_surface_extrema(self.curve(), self.surface(), 16);
                self.min_sq_distance = self.min_sq_distance.min(e.distance * e.distance);
                if e.distance * e.distance < self.criteria * self.criteria
                    && e.u1 >= anarg1 - tol
                    && e.u1 <= anarg2 + tol
                {
                    let u = e.u2;
                    let v = e.v2.unwrap_or(0.0);
                    if u >= self.umin && u <= self.umax && v >= self.vmin && v <= self.vmax {
                        let n = self.range_manager.len();
                        self.compute_range_from_start_point(false, e.u1, u, v);
                        self.compute_range_from_start_point(true, e.u1, u, v);
                        solution_found = true;
                        if n == self.range_manager.len() {
                            self.set_empty_result_range(e.u1);
                        }
                    }
                }
            } else {
                for (t, u, v) in candidates {
                    let n = self.range_manager.len();
                    self.compute_range_from_start_point(false, t, u, v);
                    self.compute_range_from_start_point(true, t, u, v);
                    solution_found = true;
                    if n == self.range_manager.len() {
                        self.set_empty_result_range(t);
                    }
                }
            }

            if !solution_found {
                self.range_manager.set_flag(i, 1);
            }
            let diff = self.range_manager.len() - old_len;
            if diff > 0 {
                i += diff;
            } else {
                i += 1;
            }
        }
    }

    /// `Extrema_ExtCS::IsParallel` for the analytic line/plane pair (direction
    /// perpendicular to the plane normal). Returns the constant square distance.
    pub(crate) fn extcs_is_parallel(&self, anarg1: f64, anarg2: f64) -> Option<f64> {
        if classify_fast_curve(self.curve()) != FastCurveKind::Line {
            return None;
        }
        if classify_surface(self.surface()) != SurfaceKind::Plane {
            return None;
        }
        let (ploc, _px, _py, pn) = plane_geometry(self.surface())?;
        let p0 = self.curve_d0(anarg1);
        let p1 = self.curve_d0(anarg2);
        let dir = GpVec::from_pnts(&p0, &p1);
        let mag = dir.magnitude();
        if mag <= 1e-30 {
            return None;
        }
        let nd = dir.dot(&GpVec::from_xyz(pn.xyz())) / mag;
        if nd.abs() > ANGULAR {
            return None;
        }
        let d = plane_distance(&ploc, &pn, &p0);
        Some(d * d)
    }

    /// Finite `(u0, v0, u1, v1)` bounds for the HInter call; unbounded
    /// directions are clamped to a broad window (they only arise for analytic
    /// surfaces handled elsewhere).
    pub(crate) fn finite_uv_bounds(&self) -> (f64, f64, f64, f64) {
        let clamp = |a: f64, b: f64| {
            if a.is_finite() && b.is_finite() {
                (a, b)
            } else {
                (-1e4, 1e4)
            }
        };
        let (u0, u1) = clamp(self.umin, self.umax);
        let (v0, v1) = clamp(self.vmin, self.vmax);
        (u0, u1, v0, v1)
    }
}
