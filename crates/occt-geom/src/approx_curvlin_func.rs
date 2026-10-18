//! `Approx_CurvlinFunc` case 2 (`Approx_CurvlinFunc.cxx:139-793`).
//!
//! Case 2 is the `Curve2d + Surface` constructor used by
//! `Approx_CurvilinearParameter` (`Approx_CurvilinearParameter.cxx:285-320`).
//! Cases 1 and 3 are not required by `BRepLib::SameParameter` and stay
//! UNPORTED.

use crate::curve::Curve;
use crate::gcpnts::{abscissa_point_from_guess, curve_length_between};
use crate::surface::Surface;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use std::sync::Arc;

fn cubic(x: f64, xi: &[f64; 4], yi: &[f64; 4]) -> f64 {
    let i1 = (yi[0] - yi[1]) / (xi[0] - xi[1]);
    let i2 = (yi[1] - yi[2]) / (xi[1] - xi[2]);
    let i3 = (yi[2] - yi[3]) / (xi[2] - xi[3]);
    let i21 = (i1 - i2) / (xi[0] - xi[2]);
    let i22 = (i2 - i3) / (xi[1] - xi[3]);
    let i31 = (i21 - i22) / (xi[0] - xi[3]);
    yi[0] + (x - xi[0]) * (i1 + (x - xi[1]) * (i21 + (x - xi[2]) * i31))
}

fn find_four_points(
    n_interval: usize,
    si: &[f64],
    ui: &[f64],
    prev_s: f64,
    prev_u: f64,
) -> Option<([f64; 4], [f64; 4])> {
    let nb_int = si.len().saturating_sub(1);
    if nb_int < 3 {
        return None;
    }
    let mut n_interval = n_interval;
    if n_interval < 1 {
        n_interval = 1;
    } else if n_interval > nb_int - 2 {
        n_interval = nb_int - 2;
    }
    let mut xi = [0.0; 4];
    let mut yi = [0.0; 4];
    for i in 0..4 {
        xi[i] = si[n_interval - 1 + i];
        yi[i] = ui[n_interval - 1 + i];
    }
    for i in 0..3 {
        if xi[i] < prev_s && prev_s < xi[i + 1] {
            for j in 0..i {
                xi[j] = xi[j + 1];
                yi[j] = yi[j + 1];
            }
            xi[i] = prev_s;
            yi[i] = prev_u;
            break;
        }
    }
    Some((xi, yi))
}

/// `Approx_CurvlinFunc` for case 2 (`cxx:139-152`).
pub struct CurvlinFunc {
    c2d1: Arc<dyn Curve2d>,
    surf1: Arc<dyn Surface>,
    /// `Adaptor3d_CurveOnSurface` standing in for length / parameter queries.
    cos1: Arc<dyn Curve>,
    first_s: f64,
    last_s: f64,
    first_u1: f64,
    last_u1: f64,
    tol_len: f64,
    length: f64,
    si_1: Vec<f64>,
    ui_1: Vec<f64>,
    prev_s: f64,
    prev_u: f64,
}

impl CurvlinFunc {
    /// `Approx_CurvlinFunc(C2D, Surf, Tol)` (`cxx:139-152`).
    ///
    /// `cos` must be the matching `Adaptor3d_CurveOnSurface(C2D, Surf)`.
    pub fn from_curve_on_surface(
        c2d: Arc<dyn Curve2d>,
        surf: Arc<dyn Surface>,
        cos: Arc<dyn Curve>,
        tol: f64,
    ) -> Self {
        let first_u1 = cos.first_parameter();
        let last_u1 = cos.last_parameter();
        let mut s = Self {
            c2d1: c2d,
            surf1: surf,
            cos1: cos,
            first_s: 0.0,
            last_s: 1.0,
            first_u1,
            last_u1,
            tol_len: tol,
            length: 0.0,
            si_1: Vec::new(),
            ui_1: Vec::new(),
            prev_s: 0.0,
            prev_u: first_u1,
        };
        s.init_tables();
        s.length = curve_length_between(s.cos1.as_ref(), s.first_u1, s.last_u1, Some(s.tol_len));
        s
    }

    /// `Init(Adaptor3d_Curve&, Si, Ui)` (`cxx:215-263`).
    fn init_tables(&mut self) {
        let first_u = self.cos1.first_parameter();
        let last_u = self.cos1.last_parameter();
        let nb_int = 10usize;
        let intervals = self.cos1.parameter_intervals(Shape::C3 as u8);
        let disc: Vec<f64> = if intervals.len() >= 2 {
            intervals
        } else {
            vec![first_u, last_u]
        };
        let nb_int_c3 = disc.len().saturating_sub(1).max(1);
        let mut ui = Vec::with_capacity(nb_int_c3 * nb_int + 1);
        let mut si = Vec::with_capacity(nb_int_c3 * nb_int + 1);
        ui.push(first_u);
        si.push(0.0);
        for j in 0..nb_int_c3 {
            let step = (disc[j + 1] - disc[j]) / nb_int as f64;
            for _k in 0..nb_int {
                let prev_u = *ui.last().unwrap();
                let next_u = prev_u + step;
                let prev_s = *si.last().unwrap();
                let seg = curve_length_between(self.cos1.as_ref(), prev_u, next_u, Some(self.tol_len));
                ui.push(next_u);
                si.push(prev_s + seg);
            }
        }
        let len = *si.last().unwrap_or(&1.0);
        if len > 0.0 {
            for v in si.iter_mut() {
                *v /= len;
            }
        }
        self.ui_1 = ui;
        self.si_1 = si;
        self.prev_s = self.first_s;
        self.prev_u = first_u;
    }

    pub fn first_parameter(&self) -> f64 {
        self.first_s
    }

    pub fn last_parameter(&self) -> f64 {
        self.last_s
    }

    pub fn length(&self) -> f64 {
        self.length
    }

    /// `NbIntervals` case 2 (`cxx:288-291`).
    pub fn nb_intervals(&self, continuity: Shape) -> i32 {
        self.cos1.nb_intervals(continuity as u8)
    }

    /// `Intervals` case 2 (`cxx:324-327`) then remap each knot by
    /// `GetSParameter` (`cxx:351-354`).
    pub fn intervals(&self, continuity: Shape) -> Vec<f64> {
        let raw = self.cos1.parameter_intervals(continuity as u8);
        if raw.len() < 2 || self.length <= CONFUSION {
            return vec![self.first_s, self.last_s];
        }
        let origin = self.cos1.first_parameter();
        raw.iter()
            .map(|&u| {
                let arc = curve_length_between(self.cos1.as_ref(), origin, u, Some(self.tol_len));
                self.first_s + arc / self.length
            })
            .collect()
    }

    /// `Trim` case 2 (`cxx:396-412`). Rebuilds the Cos through the caller's
    /// factory after updating the 2d / surface holders.
    pub fn trim<F>(&mut self, first: f64, last: f64, tol: f64, rebuild_cos: F)
    where
        F: Fn(Arc<dyn Curve2d>, Arc<dyn Surface>) -> Arc<dyn Curve>,
    {
        if first < 0.0 || last > 1.0 {
            return;
        }
        if (last - first) < tol {
            return;
        }
        // Restore full Cos window then re-trim to [FirstU, LastU] in U.
        self.cos1 = rebuild_cos(self.c2d1.clone(), self.surf1.clone());
        let first_u = self.u_parameter(first);
        let last_u = self.u_parameter(last);
        // Trim pcurve to the U window and rebuild Cos.
        let trimmed: Arc<dyn Curve2d> = Arc::new(occt_geom2d::trimmed::Geom2dTrimmedCurve::new(
            Arc::from(self.c2d1.clone_dyn()),
            first_u,
            last_u,
        ));
        self.c2d1 = trimmed.clone();
        self.cos1 = rebuild_cos(trimmed, self.surf1.clone());
        self.first_u1 = self.cos1.first_parameter();
        self.last_u1 = self.cos1.last_parameter();
        self.first_s = first;
        self.last_s = last;
        self.init_tables();
        self.length = curve_length_between(self.cos1.as_ref(), self.first_u1, self.last_u1, Some(self.tol_len));
    }

    /// `GetUParameter` case 2 (`cxx:491-572`).
    pub fn u_parameter(&mut self, s: f64) -> f64 {
        if !(0.0..=1.0).contains(&s) {
            return self.first_u1;
        }
        let init_u = &self.ui_1;
        let init_s = &self.si_1;
        let length = self.length;
        let nb_int = init_u.len().saturating_sub(1);
        if nb_int == 0 {
            return self.first_u1;
        }
        let n_interval = if s == 1.0 {
            nb_int - 1
        } else {
            let mut i = 0;
            while i < nb_int {
                if init_s[i] <= s && s < init_s[i + 1] {
                    break;
                }
                i += 1;
            }
            i.min(nb_int - 1)
        };
        if (s - init_s[n_interval]).abs() < CONFUSION {
            return init_u[n_interval];
        }
        if n_interval + 1 < init_s.len() && (s - init_s[n_interval + 1]).abs() < CONFUSION {
            return init_u[n_interval + 1];
        }
        let base = init_u[n_interval];
        let delta_s = (s - init_s[n_interval]) * length;
        let Some((xi, yi)) =
            find_four_points(n_interval, init_s, init_u, self.prev_s, self.prev_u)
        else {
            return base;
        };
        let u_guess = cubic(s, &xi, &yi);
        let u = abscissa_point_from_guess(
            self.cos1.as_ref(),
            delta_s,
            base,
            u_guess,
            self.tol_len,
        )
        .unwrap_or(u_guess);
        self.prev_s = s;
        self.prev_u = u;
        u
    }

    /// `EvalCase2` / `EvalCurOnSur` (`cxx:639-792`).
    pub fn eval_case2(&mut self, s: f64, order: i32, result: &mut [f64]) -> bool {
        if result.len() < 5 {
            return false;
        }
        let u = self.u_parameter(s);
        let length = self.length;
        match order {
            0 => {
                let c2d = self.c2d1.d0(u);
                let c = self.surf1.d0(c2d.x(), c2d.y());
                result[0] = c2d.x();
                result[1] = c2d.y();
                result[2] = c.x();
                result[3] = c.y();
                result[4] = c.z();
                true
            }
            1 => {
                let (c2d, dc2d) = self.c2d1.d1(u);
                let dv_du = dc2d.x();
                let dw_du = dc2d.y();
                let (c, ds_dv, ds_dw) = self.surf1.d1(c2d.x(), c2d.y());
                let _ = c;
                let dc_du = ds_dv.multiplied_scalar(dv_du).added(&ds_dw.multiplied_scalar(dw_du));
                let mag = dc_du.magnitude();
                if mag <= CONFUSION {
                    return false;
                }
                let du_ds = length / mag;
                result[0] = dv_du * du_ds;
                result[1] = dw_du * du_ds;
                let dc_ds = dc_du.multiplied_scalar(du_ds);
                result[2] = dc_ds.x();
                result[3] = dc_ds.y();
                result[4] = dc_ds.z();
                true
            }
            2 => {
                let (c2d, dc2d, d2c2d) = self.c2d1.d2(u);
                let dv_du = dc2d.x();
                let dw_du = dc2d.y();
                let d2v_du2 = d2c2d.x();
                let d2w_du2 = d2c2d.y();
                let (c, ds_dv, ds_dw, d2s_dv2, d2s_dw2, d2s_dvdw) =
                    self.surf1.d2(c2d.x(), c2d.y());
                let _ = c;
                let dc_du = ds_dv.multiplied_scalar(dv_du).added(&ds_dw.multiplied_scalar(dw_du));
                let d2c_du2 = d2s_dv2
                    .multiplied_scalar(dv_du)
                    .added(&d2s_dvdw.multiplied_scalar(dw_du))
                    .multiplied_scalar(dv_du)
                    .added(&ds_dv.multiplied_scalar(d2v_du2))
                    .added(
                        &d2s_dvdw
                            .multiplied_scalar(dv_du)
                            .added(&d2s_dw2.multiplied_scalar(dw_du))
                            .multiplied_scalar(dw_du),
                    )
                    .added(&ds_dw.multiplied_scalar(d2w_du2));
                let mag = dc_du.magnitude();
                if mag <= CONFUSION {
                    return false;
                }
                let du_ds = length / mag;
                let d2u_ds2 =
                    -length * dc_du.dot(&d2c_du2) * du_ds / (mag * mag * mag);
                let dv_ds = dv_du * du_ds;
                let dw_ds = dw_du * du_ds;
                let d2v_ds2 = d2v_du2 * du_ds * du_ds + dv_du * d2u_ds2;
                let d2w_ds2 = d2w_du2 * du_ds * du_ds + dw_du * d2u_ds2;
                let d2c_ds2 = d2s_dv2
                    .multiplied_scalar(dv_ds)
                    .added(&d2s_dvdw.multiplied_scalar(dw_ds))
                    .multiplied_scalar(dv_ds)
                    .added(&ds_dv.multiplied_scalar(d2v_ds2))
                    .added(
                        &d2s_dw2
                            .multiplied_scalar(dw_ds)
                            .added(&d2s_dvdw.multiplied_scalar(dv_ds))
                            .multiplied_scalar(dw_ds),
                    )
                    .added(&ds_dw.multiplied_scalar(d2w_ds2));
                result[0] = d2v_ds2;
                result[1] = d2w_ds2;
                result[2] = d2c_ds2.x();
                result[3] = d2c_ds2.y();
                result[4] = d2c_ds2.z();
                true
            }
            _ => {
                result[..5].fill(0.0);
                false
            }
        }
    }
}
