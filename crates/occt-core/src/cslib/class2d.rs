//! 2D point-in-polygon classifier. Source: `CSLib_Class2d.hxx/.cxx`.
//!
//! The polygon is normalized to `[0,1] x [0,1]` and classified with a
//! horizontal ray-cast. `SiDans` also probes the four tolerance corners so a
//! near-boundary point is reported as uncertain (ON).

use crate::gp::GpPnt2d;
use crate::precision::PCONFUSION;

const MIN_RANGE: f64 = 1e-10;

/// Classification result (`CSLib_Class2d::Result`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class2dResult {
    /// Point is strictly inside the polygon.
    Inside,
    /// Point is strictly outside the polygon.
    Outside,
    /// Point is on the boundary or classification is uncertain.
    Uncertain,
}

impl Class2dResult {
    /// Integer encoding used by `IntTools_FClass2d`: Inside=1, Outside=-1, Uncertain=0.
    pub fn as_i32(self) -> i32 {
        match self {
            Class2dResult::Inside => 1,
            Class2dResult::Outside => -1,
            Class2dResult::Uncertain => 0,
        }
    }
}

/// Low-level 2D classifier. Source: `CSLib_Class2d`.
#[derive(Debug, Clone)]
pub struct Class2d {
    pnts_x: Vec<f64>,
    pnts_y: Vec<f64>,
    tol_u: f64,
    tol_v: f64,
    points_count: usize,
    u_min: f64,
    v_min: f64,
    u_max: f64,
    v_max: f64,
}

impl Class2d {
    /// Empty classifier (`CSLib_Class2d()`).
    pub fn new_empty() -> Self {
        Self {
            pnts_x: Vec::new(),
            pnts_y: Vec::new(),
            tol_u: 0.0,
            tol_v: 0.0,
            points_count: 0,
            u_min: 0.0,
            v_min: 0.0,
            u_max: 0.0,
            v_max: 0.0,
        }
    }

    /// Construct from polygon vertices. The polygon is closed internally.
    pub fn new(
        pnts: &[GpPnt2d],
        tol_u: f64,
        tol_v: f64,
        u_min: f64,
        v_min: f64,
        u_max: f64,
        v_max: f64,
    ) -> Self {
        let mut c = Self::new_empty();
        c.init(pnts, tol_u, tol_v, u_min, v_min, u_max, v_max);
        c
    }

    fn init(
        &mut self,
        pnts: &[GpPnt2d],
        tol_u: f64,
        tol_v: f64,
        u_min: f64,
        v_min: f64,
        u_max: f64,
        v_max: f64,
    ) {
        self.u_min = u_min;
        self.v_min = v_min;
        self.u_max = u_max;
        self.v_max = v_max;
        if u_max <= u_min || v_max <= v_min || pnts.len() < 3 {
            self.points_count = 0;
            return;
        }
        self.points_count = pnts.len();
        self.tol_u = tol_u;
        self.tol_v = tol_v;
        let du = u_max - u_min;
        let dv = v_max - v_min;
        self.pnts_x = vec![0.0; self.points_count + 1];
        self.pnts_y = vec![0.0; self.points_count + 1];
        for i in 0..self.points_count {
            self.pnts_x[i] = transform_to_normalized(pnts[i].x(), u_min, du);
            self.pnts_y[i] = transform_to_normalized(pnts[i].y(), v_min, dv);
        }
        self.pnts_x[self.points_count] = self.pnts_x[0];
        self.pnts_y[self.points_count] = self.pnts_y[0];
        if du > MIN_RANGE {
            self.tol_u /= du;
        }
        if dv > MIN_RANGE {
            self.tol_v /= dv;
        }
    }

    /// Classify a point (`CSLib_Class2d::SiDans`).
    pub fn si_dans(&self, point: &GpPnt2d) -> Class2dResult {
        if self.points_count == 0 {
            return Class2dResult::Uncertain;
        }
        let mut x = point.x();
        let mut y = point.y();
        let a_tol_u = self.tol_u * (self.u_max - self.u_min);
        let a_tol_v = self.tol_v * (self.v_max - self.v_min);
        if x < (self.u_min - a_tol_u)
            || x > (self.u_max + a_tol_u)
            || y < (self.v_min - a_tol_v)
            || y > (self.v_max + a_tol_v)
        {
            return Class2dResult::Outside;
        }
        x = transform_to_normalized(x, self.u_min, self.u_max - self.u_min);
        y = transform_to_normalized(y, self.v_min, self.v_max - self.v_min);
        let result = self.internal_si_dans_ou_on(x, y);
        if result == Class2dResult::Uncertain {
            return Class2dResult::Uncertain;
        }
        if self.tol_u > 0.0 || self.tol_v > 0.0 {
            let is_inside = result == Class2dResult::Inside;
            if is_inside != self.internal_si_dans(x - self.tol_u, y - self.tol_v)
                || is_inside != self.internal_si_dans(x + self.tol_u, y - self.tol_v)
                || is_inside != self.internal_si_dans(x - self.tol_u, y + self.tol_v)
                || is_inside != self.internal_si_dans(x + self.tol_u, y + self.tol_v)
            {
                return Class2dResult::Uncertain;
            }
        }
        result
    }

    /// Classify with an explicit ON tolerance (`CSLib_Class2d::SiDans_OnMode`).
    pub fn si_dans_on_mode(&self, point: &GpPnt2d, tol: f64) -> Class2dResult {
        if self.points_count == 0 {
            return Class2dResult::Uncertain;
        }
        let mut x = point.x();
        let mut y = point.y();
        if x < (self.u_min - tol)
            || x > (self.u_max + tol)
            || y < (self.v_min - tol)
            || y > (self.v_max + tol)
        {
            return Class2dResult::Outside;
        }
        x = transform_to_normalized(x, self.u_min, self.u_max - self.u_min);
        y = transform_to_normalized(y, self.v_min, self.v_max - self.v_min);
        let result = self.internal_si_dans_ou_on(x, y);
        if tol > 0.0 {
            let is_inside = result == Class2dResult::Inside;
            if is_inside != self.internal_si_dans(x - tol, y - tol)
                || is_inside != self.internal_si_dans(x + tol, y - tol)
                || is_inside != self.internal_si_dans(x - tol, y + tol)
                || is_inside != self.internal_si_dans(x + tol, y + tol)
            {
                return Class2dResult::Uncertain;
            }
        }
        result
    }

    fn internal_si_dans(&self, px: f64, py: f64) -> bool {
        let mut nb_crossings = 0i32;
        let mut prev_dx = self.pnts_x[0] - px;
        let mut prev_dy = self.pnts_y[0] - py;
        let mut prev_y_neg = prev_dy < 0.0;
        for next_idx in 1..=self.points_count {
            let curr_dx = self.pnts_x[next_idx] - px;
            let curr_dy = self.pnts_y[next_idx] - py;
            let curr_y_neg = curr_dy < 0.0;
            if curr_y_neg != prev_y_neg {
                if prev_dx > 0.0 && curr_dx > 0.0 {
                    nb_crossings += 1;
                } else if prev_dx > 0.0 || curr_dx > 0.0 {
                    let x_intersect = prev_dx - prev_dy * (curr_dx - prev_dx) / (curr_dy - prev_dy);
                    if x_intersect > 0.0 {
                        nb_crossings += 1;
                    }
                }
                prev_y_neg = curr_y_neg;
            }
            prev_dx = curr_dx;
            prev_dy = curr_dy;
        }
        (nb_crossings & 1) != 0
    }

    fn internal_si_dans_ou_on(&self, px: f64, py: f64) -> Class2dResult {
        let mut nb_crossings = 0i32;
        let mut prev_dx = self.pnts_x[0] - px;
        let mut prev_dy = self.pnts_y[0] - py;
        let mut prev_y_neg = prev_dy < 0.0;
        for next_idx in 1..=self.points_count {
            let prev_idx = next_idx - 1;
            let curr_dx = self.pnts_x[next_idx] - px;
            let curr_dy = self.pnts_y[next_idx] - py;
            if curr_dx < self.tol_u
                && curr_dx > -self.tol_u
                && curr_dy < self.tol_v
                && curr_dy > -self.tol_v
            {
                return Class2dResult::Uncertain;
            }
            let edge_dx = self.pnts_x[next_idx] - self.pnts_x[prev_idx];
            if (self.pnts_x[prev_idx] - px) * curr_dx < 0.0 && edge_dx.abs() > PCONFUSION {
                let interp_y = self.pnts_y[next_idx]
                    - (self.pnts_y[next_idx] - self.pnts_y[prev_idx]) / edge_dx * curr_dx;
                let delta_y = interp_y - py;
                if delta_y >= -self.tol_v && delta_y <= self.tol_v {
                    return Class2dResult::Uncertain;
                }
            }
            let curr_y_neg = curr_dy < 0.0;
            if curr_y_neg != prev_y_neg {
                if prev_dx > 0.0 && curr_dx > 0.0 {
                    nb_crossings += 1;
                } else if prev_dx > 0.0 || curr_dx > 0.0 {
                    let x_intersect = prev_dx - prev_dy * (curr_dx - prev_dx) / (curr_dy - prev_dy);
                    if x_intersect > 0.0 {
                        nb_crossings += 1;
                    }
                }
                prev_y_neg = curr_y_neg;
            }
            prev_dx = curr_dx;
            prev_dy = curr_dy;
        }
        if (nb_crossings & 1) != 0 {
            Class2dResult::Inside
        } else {
            Class2dResult::Outside
        }
    }
}

fn transform_to_normalized(u: f64, u_min: f64, u_range: f64) -> f64 {
    if u_range > MIN_RANGE {
        (u - u_min) / u_range
    } else {
        u
    }
}
