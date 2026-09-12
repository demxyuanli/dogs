//! Analytic BeanFace paths: line/plane, FastComputeAnalytic, TestComputeCoinside.
use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::ANGULAR;

use crate::bean_face::BeanFaceIntersector;
use crate::bean_face_kind::{
    circle_geometry, classify_fast_curve, cylinder_geometry, ellipse_geometry, line_geometry,
    plane_distance, plane_geometry, sphere_geometry, FastCurveKind,
};
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::inttools_data::IntRange;
use crate::meshing::range_splitter::{classify_surface as classify_surface_mesh, SurfaceType};

const PI: f64 = std::f64::consts::PI;

impl BeanFaceIntersector {
    /// Line × plane intersection: substitute the line into the plane equation,
    /// then emit either a single root range (expanded by the tolerance-derived
    /// parameter width) or the whole range when the line lies in the plane.
    /// Port of `ComputeLinePlane`.
    pub(crate) fn compute_line_plane(&mut self) {
        let tol_ang = 1e-9;
        self.is_done = true;

        let (ploc, _px, _py, pn) = match plane_geometry(self.surface()) {
            Some(g) => g,
            None => return,
        };
        let (orig, ld) = match line_geometry(self.curve()) {
            Some(g) => g,
            None => return,
        };
        let nrm = GpVec::from_xyz(pn.xyz());
        let (a, b, c) = (nrm.x(), nrm.y(), nrm.z());
        let dcoef = -nrm.dot(&GpVec::from_pnts(&GpPnt::zero(), &ploc));
        let (al, bl, cl) = (ld.x(), ld.y(), ld.z());
        let direc = a * al + b * bl + c * cl;
        let dis = a * orig.x() + b * orig.y() + c * orig.z() + dcoef;

        let (mut parallel, mut inplane) = (false, false);
        if direc.abs() < tol_ang {
            parallel = true;
            inplane = dis.abs() < self.criteria;
        } else {
            let p1 = self.curve_d0(self.first_parameter);
            let p2 = self.curve_d0(self.last_parameter);
            let mut d1 = a * p1.x() + b * p1.y() + c * p1.z() + dcoef;
            if d1 < 0.0 {
                d1 = -d1;
            }
            let mut d2 = a * p2.x() + b * p2.y() + c * p2.z() + dcoef;
            if d2 < 0.0 {
                d2 = -d2;
            }
            if d1 <= self.criteria && d2 <= self.criteria {
                inplane = true;
            }
        }

        if inplane {
            // OCCT `ComputeLinePlane`: the whole bean range is the result.
            // Face restriction is applied later by `IntTools_EdgeFace::IsProjectable`.
            self.results.push(IntRange::new_unchecked(
                self.first_parameter,
                self.last_parameter,
            ));
            return;
        }
        if parallel {
            return;
        }

        let t = -dis / direc;
        if t < self.first_parameter || t > self.last_parameter {
            return;
        }
        let pint = orig.translated_vec(&GpVec::from_xyz(ld.xyz()).multiplied_scalar(t));
        let (u, v) = plane_uv_of_point(&ploc, &ld, &pint, &_px, &_py);
        if self.umin > u || u > self.umax || self.vmin > v || v > self.vmax {
            return;
        }

        // Parameter half-width from the tolerances and the incidence angle.
        let angle = (PI * 0.5 - ld.angle(&pn)).abs();
        let a_dt = compute_int_range(self.bean_tolerance, self.face_tolerance, angle);
        let t1 = self.first_parameter.max(t - a_dt);
        let t2 = self.last_parameter.min(t + a_dt);
        self.results.push(IntRange::new_unchecked(t1, t2));
    }

    /// Fast analytic coincidence / no-intersection checks for conic curves
    /// against quadric surfaces. Returns `true` when a decisive verdict was
    /// reached; otherwise computation continues. Port of `FastComputeAnalytic`.
    pub(crate) fn fast_compute_analytic(&mut self) -> bool {
        let ck = classify_fast_curve(self.curve());
        if ck == FastCurveKind::Other {
            return false;
        }
        let sk = classify_surface(self.surface());
        let mut is_coincide = false;
        let mut has_intersection = true;

        match sk {
            SurfaceKind::Plane => {
                let (ploc, _px, _py, pn) = match plane_geometry(self.surface()) {
                    Some(g) => g,
                    None => return false,
                };
                let (adir, aloc) = match ck {
                    FastCurveKind::Circle => {
                        let (c, _r, n) = match circle_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        (n, c)
                    }
                    FastCurveKind::Ellipse => {
                        let (c, n) = match ellipse_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        (n, c)
                    }
                    _ => return false,
                };
                let angle = adir.angle(&pn);
                if angle > ANGULAR {
                    return false;
                }
                has_intersection = false;
                let dist = plane_distance(&ploc, &pn, &aloc);
                is_coincide = dist < self.criteria;
            }
            SurfaceKind::Sphere => {
                let (sc, sr) = match sphere_geometry(self.surface()) {
                    Some(g) => g,
                    None => return false,
                };
                if ck == FastCurveKind::Line {
                    let (lloc, ldir) = match line_geometry(self.curve()) {
                        Some(g) => g,
                        None => return false,
                    };
                    let lin = GpLin::from_pnt_dir(lloc, ldir);
                    let dist = lin.distance(&sc) - sr;
                    has_intersection = dist < self.criteria;
                } else {
                    return false;
                }
            }
            SurfaceKind::Cylinder | SurfaceKind::Other => {
                // `brep_surface::classify_surface` never returns `Cylinder`
                // directly; re-confirm via the mesh classifier.
                if classify_surface_mesh(self.surface()) != SurfaceType::Cylinder {
                    return false;
                }
                let (axis, radius) = match cylinder_geometry(self.surface()) {
                    Some(g) => g,
                    None => return false,
                };
                match ck {
                    FastCurveKind::Line => {
                        let (lloc, ldir) = match line_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        if !ldir.is_parallel(axis.direction()) {
                            return false;
                        }
                        has_intersection = false;
                        let lin = GpLin::from_pnt_dir(lloc, ldir);
                        let dist = (lin.distance(axis.location()) - radius).abs();
                        is_coincide = dist < self.criteria;
                    }
                    FastCurveKind::Circle => {
                        let (cloc, cr, cn) = match circle_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        let angle = axis.direction().angle(&cn);
                        if angle > ANGULAR {
                            return false;
                        }
                        let axis_lin = GpLin::from_pnt_dir(*axis.location(), *axis.direction());
                        let dist_loc = axis_lin.distance(&cloc);
                        let dist = dist_loc + (cr - radius).abs();
                        is_coincide = dist < self.criteria;
                        if !is_coincide {
                            has_intersection = (dist_loc - (cr + radius)) < self.criteria
                                && ((cr - radius).abs() - dist_loc) < self.criteria;
                        }
                    }
                    _ => return false,
                }
            }
            SurfaceKind::Cone | SurfaceKind::Torus => return false,
        }

        if is_coincide {
            self.results.push(IntRange::new_unchecked(self.first_parameter, self.last_parameter));
        }
        is_coincide || !has_intersection
    }

    /// Scan the whole bean for coincidence with the surface: sample 23 points,
    /// expand ranges from each. Port of `TestComputeCoinside`.
    pub(crate) fn test_compute_coinside(&mut self) -> bool {
        let cfp = self.first_parameter;
        let clp = self.last_parameter;
        let nb_seg = 23;
        let cdp = (clp - cfp) / nb_seg as f64;

        let mut u = 0.0;
        let mut v = 0.0;
        if self.distance_with_uv(cfp, &mut u, &mut v) > self.criteria {
            return false;
        }
        self.compute_range_from_start_point(true, cfp, u, v);

        let found = self.range_manager.get_index(clp, false);
        if found >= 0 && self.range_manager.flag(found as usize) == 2 {
            return true;
        }
        if self.distance_with_uv(clp, &mut u, &mut v) > self.criteria {
            return false;
        }
        self.compute_range_from_start_point(false, clp, u, v);

        for i in 1..nb_seg {
            let par = cfp + i as f64 * cdp;
            if self.distance_with_uv(par, &mut u, &mut v) > self.criteria {
                return false;
            }
            let n = self.range_manager.len();
            self.compute_range_from_start_point(false, par, u, v);
            self.compute_range_from_start_point(true, par, u, v);
            if n == self.range_manager.len() {
                self.set_empty_result_range(par);
            }
        }
        true
    }
}

/// Wrap `value` into `[min, max]` by whole periods (positive modulo).
pub(crate) fn adjust_periodic(value: f64, min: f64, max: f64, period: f64) -> f64 {
    let p = period.abs();
    if p <= 1e-30 || !min.is_finite() || !max.is_finite() {
        return value;
    }
    let mut v = value;
    while v < min {
        v += p;
    }
    while v > max {
        v -= p;
    }
    v
}

/// `(u, v)` parameters of `p` in a plane's natural frame.
pub(crate) fn plane_uv_of_point(ploc: &GpPnt, _pn: &GpDir, p: &GpPnt, px: &GpDir, py: &GpDir) -> (f64, f64) {
    let d = GpVec::from_pnts(ploc, p);
    (d.dot(&GpVec::from_xyz(px.xyz())), d.dot(&GpVec::from_xyz(py.xyz())))
}

/// Port of `IntTools_Tools::ComputeIntRange`: the parameter half-width that
/// covers the tolerance band around a crossing at incidence `angle`.
pub(crate) fn compute_int_range(tol1: f64, tol2: f64, angle: f64) -> f64 {
    if (PI * 0.5 - angle).abs() < ANGULAR {
        return tol2;
    }
    let an_angle = if angle > PI * 0.5 { PI - angle } else { angle };
    let a1 = tol1 * (PI * 0.5 - an_angle).tan();
    let a2 = tol2 / an_angle.sin();
    a1 + a2
}
