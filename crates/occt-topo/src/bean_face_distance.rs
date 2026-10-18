//! `IntTools_BeanFaceIntersector::Distance` with isoline fallback.
//!
//! Source: `IntTools_BeanFaceIntersector.cxx:397-560`. When the point-on-surface
//! projector fails, OCCT projects onto the four boundary isos
//! (`UIso`/`VIso` + `GeomAPI_ProjectPointOnCurve`). `dyn Surface` has no iso
//! handle, so the iso is sampled and refined in its parameter.

use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::bean_face::BeanFaceIntersector;

const ISO_SAMPLES: usize = 32;

impl BeanFaceIntersector {
    /// Distance from the curve point at `arg` to the surface; closest surface
    /// parameters (clamped to the window) are written to `u`/`v`.
    /// Port of `Distance(theArg, theUParameter, theVParameter)`.
    pub(crate) fn distance_with_uv(&self, arg: f64, u: &mut f64, v: &mut f64) -> f64 {
        let p = self.curve_d0(arg);
        *u = self.umin;
        *v = self.vmin;
        let (su, sv, d) = self.closest_params_dist(&p);
        let projection_found = d.is_finite();
        let mut distance_val = d;
        let mut u_val = su;
        let mut v_val = sv;

        if !projection_found {
            distance_val = f64::MAX;
            for i in 0..4 {
                let an_iso = if i == 0 {
                    self.umin
                } else if i == 1 {
                    self.umax
                } else if i == 2 {
                    self.vmin
                } else {
                    self.vmax
                };
                let min_val = if i < 2 { self.vmin } else { self.umin };
                let max_val = if i < 2 { self.vmax } else { self.umax };
                if !min_val.is_finite() || !max_val.is_finite() || !an_iso.is_finite() {
                    continue;
                }
                let mid = (min_val + max_val) * 0.5;
                let point_min = if i < 2 {
                    self.surface_d0(an_iso, min_val)
                } else {
                    self.surface_d0(min_val, an_iso)
                };
                let point_max = if i < 2 {
                    self.surface_d0(an_iso, max_val)
                } else {
                    self.surface_d0(max_val, an_iso)
                };
                let point_mid = if i < 2 {
                    self.surface_d0(an_iso, mid)
                } else {
                    self.surface_d0(mid, an_iso)
                };

                let mut use_min_max = true;
                let mut compute_iso = true;
                if point_min.distance(&point_max) <= self.criteria
                    && point_min.distance(&point_mid) <= self.criteria
                    && point_max.distance(&point_mid) <= self.criteria
                {
                    compute_iso = false;
                }

                if compute_iso {
                    if let Some((t_iso, dist_iso)) =
                        project_point_on_iso(self.surface(), &p, i < 2, an_iso, min_val, max_val)
                    {
                        use_min_max = false;
                        if distance_val > dist_iso {
                            u_val = if i <= 1 { an_iso } else { t_iso };
                            v_val = if i >= 2 { an_iso } else { t_iso };
                            distance_val = dist_iso;
                        }
                    }
                }

                if use_min_max {
                    let dmin = p.distance(&point_min);
                    if dmin < distance_val {
                        u_val = if i <= 1 { an_iso } else { min_val };
                        v_val = if i >= 2 { an_iso } else { min_val };
                        distance_val = dmin;
                    }
                    let dmax = p.distance(&point_max);
                    if dmax < distance_val {
                        u_val = if i <= 1 { an_iso } else { max_val };
                        v_val = if i >= 2 { an_iso } else { max_val };
                        distance_val = dmax;
                    }
                }
            }
        }

        *u = u_val.clamp(self.umin.min(self.umax), self.umin.max(self.umax));
        *v = v_val.clamp(self.vmin.min(self.vmax), self.vmin.max(self.vmax));
        if self.umin > *u {
            *u = self.umin;
        }
        if self.umax < *u {
            *u = self.umax;
        }
        if self.vmin > *v {
            *v = self.vmin;
        }
        if self.vmax < *v {
            *v = self.vmax;
        }
        distance_val
    }

    /// Distance from the curve point at `arg` to the surface.
    /// Port of `Distance(theArg)`.
    pub(crate) fn distance(&self, arg: f64) -> f64 {
        let p = self.curve_d0(arg);
        let (_u, _v, d) = self.closest_params_dist(&p);
        if d.is_finite() {
            return d;
        }
        let mut distance_val = f64::MAX;
        for i in 0..4 {
            let an_iso = if i == 0 {
                self.umin
            } else if i == 1 {
                self.umax
            } else if i == 2 {
                self.vmin
            } else {
                self.vmax
            };
            let min_val = if i < 2 { self.vmin } else { self.umin };
            let max = if i < 2 { self.vmax } else { self.umax };
            if !min_val.is_finite() || !max.is_finite() || !an_iso.is_finite() {
                continue;
            }
            let mid = (min_val + max) * 0.5;
            let point_min = if i < 2 {
                self.surface_d0(an_iso, min_val)
            } else {
                self.surface_d0(min_val, an_iso)
            };
            let point_max = if i < 2 {
                self.surface_d0(an_iso, max)
            } else {
                self.surface_d0(max, an_iso)
            };
            let point_mid = if i < 2 {
                self.surface_d0(an_iso, mid)
            } else {
                self.surface_d0(mid, an_iso)
            };
            let mut use_min_max = true;
            let mut compute_iso = true;
            if point_min.distance(&point_max) <= self.criteria
                && point_min.distance(&point_mid) <= self.criteria
                && point_max.distance(&point_mid) <= self.criteria
            {
                compute_iso = false;
            }
            if compute_iso {
                if let Some((_t, dist_iso)) =
                    project_point_on_iso(self.surface(), &p, i < 2, an_iso, min_val, max)
                {
                    use_min_max = false;
                    if distance_val > dist_iso {
                        distance_val = dist_iso;
                    }
                }
            }
            if use_min_max {
                distance_val = distance_val.min(p.distance(&point_min));
                distance_val = distance_val.min(p.distance(&point_max));
            }
        }
        distance_val
    }
}

/// Closest parameter and distance of `p` on a U-iso (`iso_u`) or V-iso.
fn project_point_on_iso(
    surf: &dyn Surface,
    p: &GpPnt,
    iso_u: bool,
    iso_param: f64,
    tmin: f64,
    tmax: f64,
) -> Option<(f64, f64)> {
    if !tmin.is_finite() || !tmax.is_finite() || (tmax - tmin).abs() <= f64::EPSILON {
        return None;
    }
    let eval = |t: f64| {
        if iso_u {
            surf.d0(iso_param, t)
        } else {
            surf.d0(t, iso_param)
        }
    };
    let mut best_t = tmin;
    let mut best_d = f64::MAX;
    for i in 0..=ISO_SAMPLES {
        let t = tmin + (tmax - tmin) * i as f64 / ISO_SAMPLES as f64;
        let d = eval(t).distance(p);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    let span = (tmax - tmin).abs();
    let lo = (best_t - span / ISO_SAMPLES as f64).max(tmin.min(tmax));
    let hi = (best_t + span / ISO_SAMPLES as f64).min(tmin.max(tmax));
    let f = |t: f64| eval(t).distance(p);
    let (t, d) = crate::bean_face::golden_1d(&f, lo, hi, 1e-12);
    if d.is_finite() {
        Some((t, d))
    } else {
        None
    }
}
