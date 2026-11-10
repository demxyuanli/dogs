//! `ExtPcState` and `ElCLib::InPeriod`.

use super::*;

pub(crate) struct ExtPcState {
    pub(crate) sols: Vec<ExtPcSolution>,
    pub(crate) done: bool,
    pub(crate) dist1: f64,
    pub(crate) dist2: f64,
    pub(crate) pf: GpPnt,
    pub(crate) pl: GpPnt,
    pub(crate) tol_u: f64,
    pub(crate) u_inf: f64,
    pub(crate) u_sup: f64,
    pub(crate) sample: i32,
    pub(crate) int_u_inf: f64,
    pub(crate) int_u_sup: f64,
}

impl ExtPcState {
    pub(crate) fn new(curve: &dyn Curve, u_inf: f64, u_sup: f64) -> Self {
        Self {
            sols: Vec::new(),
            done: false,
            dist1: f64::MAX,
            dist2: f64::MAX,
            pf: GpPnt::zero(),
            pl: GpPnt::zero(),
            tol_u: curve.resolution(CONFUSION),
            u_inf,
            u_sup,
            sample: 17,
            int_u_inf: u_inf,
            int_u_sup: u_sup,
        }
    }

    /// `Extrema_GGExtPC::AddSol` (`hxx:617-632`).
    pub(crate) fn add_sol(&mut self, u: f64, p: GpPnt, sq_dist: f64, is_min: bool) {
        for s in &self.sols {
            if (s.u - u).abs() <= self.tol_u {
                return;
            }
        }
        self.sols.push(ExtPcSolution {
            u,
            point: p,
            sq_dist,
            is_min,
        });
    }

    /// `IntervalPerform` (`hxx:590-614`) via `math_FunctionRoots` on `GFuncExtPC`.
    pub(crate) fn interval_perform(&mut self, curve: &dyn Curve, point: &GpPnt) {
        let mut f = ExtPcf::new(curve, point, self.int_u_inf, self.int_u_sup);
        f.sub_interval_initialize(self.int_u_inf, self.int_u_sup);
        f.clear_sols();
        let roots = FunctionRoots::new(
            &mut f,
            self.int_u_inf,
            self.int_u_sup,
            self.sample,
            self.tol_u,
            TOL_F,
            TOL_F,
            0.0,
        );
        self.done = roots.is_done() && !roots.is_all_null();
        if !self.done {
            return;
        }
        let period = if curve.is_periodic() {
            curve.period()
        } else {
            0.0
        };
        for s in f.sols {
            let mut u = s.u;
            if period.abs() > 0.0 {
                u = elclib_in_period(u, self.u_inf, self.u_inf + period);
            }
            if u >= self.u_inf - self.tol_u && u <= self.u_sup + self.tol_u {
                self.add_sol(u, s.point, s.sq_dist, s.is_min);
            }
        }
    }
}

/// `ElCLib::InPeriod` (`ElCLib.cxx`).
pub(crate) fn elclib_in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    let period = u_last - u_first;
    if period <= 0.0 {
        return u;
    }
    let mut x = u;
    while x < u_first {
        x += period;
    }
    while x > u_last || (x - u_last).abs() <= PCONFUSION {
        if (x - u_last).abs() <= PCONFUSION {
            return u_last;
        }
        x -= period;
        if x < u_first {
            break;
        }
    }
    if x < u_first {
        x += period;
    }
    x
}

