//! `CSLib_Class2d` — point-in-polygon classification with an on-boundary band.
//!
//! Source: `FoundationClasses/TKMath/CSLib/CSLib_Class2d.cxx:45-342` (OCCT
//! 8.0.0), declared `CSLib_Class2d.hxx:41-207`.
//!
//! The constructor normalises the polygon and its tolerances to `[0, 1]` on
//! both axes, closes the polygon by repeating the first point (`cxx:85-87`) and
//! keeps `myPointsCount` as the number of input points, so the closed ring has
//! `myPointsCount + 1` entries.
//!
//! Translation notes:
//! - the three OCCT constructors (`NCollection_Array1` / `Sequence` /
//!   `DynamicArray`) all forward to the same `init` template, so the port has
//!   one `new` taking a slice,
//! - `Result` becomes [`Class2dResult`]; `SiDans` / `SiDans_OnMode` keep their
//!   two `internalSiDans*` helpers as private methods.

use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;

/// `CSLib_Class2d::Result` (`CSLib_Class2d.hxx:42-48`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class2dResult {
    /// `Result_Uncertain` — on the boundary or undecidable.
    Uncertain,
    /// `Result_Inside`.
    Inside,
    /// `Result_Outside`.
    Outside,
}

/// `transformToNormalized` (`cxx:30-37`): `THE_MIN_RANGE = 1e-10`.
const THE_MIN_RANGE: f64 = 1e-10;

fn transform_to_normalized(u: f64, u_min: f64, u_range: f64) -> f64 {
    if u_range > THE_MIN_RANGE {
        (u - u_min) / u_range
    } else {
        u
    }
}

/// `CSLib_Class2d`.
#[derive(Debug, Clone)]
pub struct Class2d {
    /// `myPnts2dX` — `myPointsCount + 1` normalised X coordinates.
    pnts2d_x: Vec<f64>,
    /// `myPnts2dY` — `myPointsCount + 1` normalised Y coordinates.
    pnts2d_y: Vec<f64>,
    /// `myPointsCount`
    points_count: usize,
    /// `myTolU` — normalised.
    tol_u: f64,
    /// `myTolV` — normalised.
    tol_v: f64,
    /// `myUMin`
    u_min: f64,
    /// `myVMin`
    v_min: f64,
    /// `myUMax`
    u_max: f64,
    /// `myVMax`
    v_max: f64,
}

impl Default for Class2d {
    fn default() -> Self {
        Self {
            pnts2d_x: Vec::new(),
            pnts2d_y: Vec::new(),
            points_count: 0,
            tol_u: 0.0,
            tol_v: 0.0,
            u_min: 0.0,
            v_min: 0.0,
            u_max: 0.0,
            v_max: 0.0,
        }
    }
}

impl Class2d {
    /// `init(thePnts2d, theTolU, theTolV, theUMin, theVMin, theUMax, theVMax)`
    /// (`cxx:45-96`), shared by the three constructors.
    pub fn new(
        points: &[GpPnt2d],
        tol_u: f64,
        tol_v: f64,
        u_min: f64,
        v_min: f64,
        u_max: f64,
        v_max: f64,
    ) -> Self {
        let mut this = Self::default();
        this.u_min = u_min;
        this.v_min = v_min;
        this.u_max = u_max;
        this.v_max = v_max;

        // `if (theUMax <= theUMin || theVMax <= theVMin || thePnts2d.Length() < 3)`
        // (`cxx:58-63`).
        if u_max <= u_min || v_max <= v_min || points.len() < 3 {
            this.points_count = 0;
            return this;
        }

        this.points_count = points.len();
        this.tol_u = tol_u;
        this.tol_v = tol_v;

        let a_du = u_max - u_min;
        let a_dv = v_max - v_min;

        this.pnts2d_x = Vec::with_capacity(this.points_count + 1);
        this.pnts2d_y = Vec::with_capacity(this.points_count + 1);
        for p in points {
            this.pnts2d_x
                .push(transform_to_normalized(p.x(), u_min, a_du));
            this.pnts2d_y
                .push(transform_to_normalized(p.y(), v_min, a_dv));
        }
        // `myPnts2dX.ChangeLast() = myPnts2dX.First();` (`cxx:86-87`).
        this.pnts2d_x.push(this.pnts2d_x[0]);
        this.pnts2d_y.push(this.pnts2d_y[0]);

        // `cxx:89-97`.
        if a_du > THE_MIN_RANGE {
            this.tol_u /= a_du;
        }
        if a_dv > THE_MIN_RANGE {
            this.tol_v /= a_dv;
        }
        this
    }

    /// `SiDans(thePoint)` (`cxx:142-187`).
    pub fn si_dans(&self, point: &GpPnt2d) -> Class2dResult {
        if self.points_count == 0 {
            return Class2dResult::Uncertain;
        }

        let mut a_x = point.x();
        let mut a_y = point.y();

        // `cxx:150-161`: the tolerances back in the original space.
        let a_tol_u = self.tol_u * (self.u_max - self.u_min);
        let a_tol_v = self.tol_v * (self.v_max - self.v_min);

        if a_x < (self.u_min - a_tol_u)
            || a_x > (self.u_max + a_tol_u)
            || a_y < (self.v_min - a_tol_v)
            || a_y > (self.v_max + a_tol_v)
        {
            return Class2dResult::Outside;
        }

        a_x = transform_to_normalized(a_x, self.u_min, self.u_max - self.u_min);
        a_y = transform_to_normalized(a_y, self.v_min, self.v_max - self.v_min);

        let a_result = self.internal_si_dans_ou_on(a_x, a_y);
        if a_result == Class2dResult::Uncertain {
            return Class2dResult::Uncertain;
        }

        // `cxx:175-185`: the four corners of the ±tol box.
        if self.tol_u > 0.0 || self.tol_v > 0.0 {
            let is_inside = a_result == Class2dResult::Inside;
            if is_inside != self.internal_si_dans(a_x - self.tol_u, a_y - self.tol_v)
                || is_inside != self.internal_si_dans(a_x + self.tol_u, a_y - self.tol_v)
                || is_inside != self.internal_si_dans(a_x - self.tol_u, a_y + self.tol_v)
                || is_inside != self.internal_si_dans(a_x + self.tol_u, a_y + self.tol_v)
            {
                return Class2dResult::Uncertain;
            }
        }

        a_result
    }

    /// `SiDans_OnMode(thePoint, theTol)` (`cxx:192-230`).
    pub fn si_dans_on_mode(&self, point: &GpPnt2d, tol: f64) -> Class2dResult {
        if self.points_count == 0 {
            return Class2dResult::Uncertain;
        }

        let mut a_x = point.x();
        let mut a_y = point.y();

        // `cxx:202-208`.
        if a_x < (self.u_min - tol)
            || a_x > (self.u_max + tol)
            || a_y < (self.v_min - tol)
            || a_y > (self.v_max + tol)
        {
            return Class2dResult::Outside;
        }

        a_x = transform_to_normalized(a_x, self.u_min, self.u_max - self.u_min);
        a_y = transform_to_normalized(a_y, self.v_min, self.v_max - self.v_min);

        let a_result = self.internal_si_dans_ou_on(a_x, a_y);

        // `cxx:218-228`.
        if tol > 0.0 {
            let is_inside = a_result == Class2dResult::Inside;
            if is_inside != self.internal_si_dans(a_x - tol, a_y - tol)
                || is_inside != self.internal_si_dans(a_x + tol, a_y - tol)
                || is_inside != self.internal_si_dans(a_x - tol, a_y + tol)
                || is_inside != self.internal_si_dans(a_x + tol, a_y + tol)
            {
                return Class2dResult::Uncertain;
            }
        }

        a_result
    }

    /// `internalSiDans(thePx, thePy)` (`cxx:235-276`): the ray-casting count
    /// without the on-boundary detection.
    fn internal_si_dans(&self, px: f64, py: f64) -> bool {
        let mut nb_crossings = 0;

        let mut prev_dx = self.pnts2d_x[0] - px;
        let mut prev_dy = self.pnts2d_y[0] - py;
        let mut prev_y_is_negative = prev_dy < 0.0;

        for next_idx in 1..=self.points_count {
            let curr_dx = self.pnts2d_x[next_idx] - px;
            let curr_dy = self.pnts2d_y[next_idx] - py;
            let curr_y_is_negative = curr_dy < 0.0;

            if curr_y_is_negative != prev_y_is_negative {
                if prev_dx > 0.0 && curr_dx > 0.0 {
                    nb_crossings += 1;
                } else if prev_dx > 0.0 || curr_dx > 0.0 {
                    let x_intersect = prev_dx - prev_dy * (curr_dx - prev_dx) / (curr_dy - prev_dy);
                    if x_intersect > 0.0 {
                        nb_crossings += 1;
                    }
                }
                prev_y_is_negative = curr_y_is_negative;
            }

            prev_dx = curr_dx;
            prev_dy = curr_dy;
        }

        (nb_crossings & 1) != 0
    }

    /// `internalSiDansOuOn(thePx, thePy)` (`cxx:279-342`).
    fn internal_si_dans_ou_on(&self, px: f64, py: f64) -> Class2dResult {
        let mut nb_crossings = 0;

        let mut prev_dx = self.pnts2d_x[0] - px;
        let mut prev_dy = self.pnts2d_y[0] - py;
        let mut prev_y_is_negative = prev_dy < 0.0;

        for next_idx in 1..=self.points_count {
            let prev_idx = next_idx - 1;
            let curr_dx = self.pnts2d_x[next_idx] - px;
            let curr_dy = self.pnts2d_y[next_idx] - py;

            // `cxx:296-300`.
            if curr_dx < self.tol_u && curr_dx > -self.tol_u && curr_dy < self.tol_v
                && curr_dy > -self.tol_v
            {
                return Class2dResult::Uncertain;
            }

            // `cxx:302-316`.
            let edge_dx = self.pnts2d_x[next_idx] - self.pnts2d_x[prev_idx];
            if (self.pnts2d_x[prev_idx] - px) * curr_dx < 0.0 && edge_dx.abs() > PCONFUSION {
                let interp_y = self.pnts2d_y[next_idx]
                    - (self.pnts2d_y[next_idx] - self.pnts2d_y[prev_idx]) / edge_dx * curr_dx;
                let delta_y = interp_y - py;
                if delta_y >= -self.tol_v && delta_y <= self.tol_v {
                    return Class2dResult::Uncertain;
                }
            }

            let curr_y_is_negative = curr_dy < 0.0;
            if curr_y_is_negative != prev_y_is_negative {
                if prev_dx > 0.0 && curr_dx > 0.0 {
                    nb_crossings += 1;
                } else if prev_dx > 0.0 || curr_dx > 0.0 {
                    let x_intersect = prev_dx - prev_dy * (curr_dx - prev_dx) / (curr_dy - prev_dy);
                    if x_intersect > 0.0 {
                        nb_crossings += 1;
                    }
                }
                prev_y_is_negative = curr_y_is_negative;
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
