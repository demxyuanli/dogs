//! `GFuncExtPC` construction and `MaxDerivOrder` / `SearchOfTolerance`.

use super::*;

/// `GFuncExtPC` Initialize GetType arm (`hxx:87-99`, `115-128`).
pub(crate) fn gfunc_max_deriv_order(curve: &dyn Curve) -> i32 {
    if curve.bspline_poles().is_some()
        || curve.bspline_knots().is_some()
        || curve.bezier_poles().is_some()
    {
        return MAX_ORDER;
    }
    // GeomAbs_OffsetCurve / GeomAbs_OtherCurve (no elementary line/circle).
    if curve.is_line() || curve.gp_circ().is_some() || curve.circle_radius().is_some() {
        0
    } else {
        MAX_ORDER
    }
}

/// `GFuncExtPC::SearchOfTolerance` (`hxx:428-457`) on curve First/Last.
pub(crate) fn gfunc_search_of_tolerance(curve: &dyn Curve) -> f64 {
    const N_POINT: i32 = 10;
    let u_inf = curve.first_parameter();
    let u_sup = curve.last_parameter();
    let step = (u_sup - u_inf) / (N_POINT as f64);
    let mut num = 0;
    let mut max_val = -Precision::INFINITE;
    loop {
        let mut u = u_inf + (num as f64) * step;
        if u > u_sup {
            u = u_sup;
        }
        let (_p, v_der) = curve.d1(u);
        if !(Precision::is_infinite(v_der.x()) || Precision::is_infinite(v_der.y())) {
            let vm = v_der.magnitude();
            if vm > max_val {
                max_val = vm;
            }
        }
        num += 1;
        if num >= N_POINT + 1 {
            break;
        }
    }
    (max_val * TOL_FACTOR).max(MIN_TOL)
}

impl<'a> ExtPcf<'a> {
    pub(crate) fn new(curve: &'a dyn Curve, point: &GpPnt, u_inf: f64, u_sup: f64) -> Self {
        let max_deriv_order = gfunc_max_deriv_order(curve);
        let tol = if max_deriv_order != 0 {
            gfunc_search_of_tolerance(curve)
        } else {
            MIN_TOL
        };
        Self {
            curve,
            point: *point,
            u: 0.0,
            pc: GpPnt::zero(),
            d1f: 0.0,
            u_inf,
            u_sup,
            max_deriv_order,
            tol,
            sols: Vec::new(),
        }
    }

    pub(crate) fn clear_sols(&mut self) {
        self.sols.clear();
    }

    pub(crate) fn sub_interval_initialize(&mut self, u_first: f64, u_last: f64) {
        self.u_inf = u_first;
        self.u_sup = u_last;
    }

    /// Subinterval span for singular delta (`hxx:170-176`, `301-307`).
    pub(crate) fn singular_du(&self) -> f64 {
        if !self.u_sup.is_finite() || !self.u_inf.is_finite() {
            0.0
        } else {
            self.u_sup - self.u_inf
        }
    }

    /// `GFuncExtPC::Value` (`hxx:146-253`) including MaxDerivOrder singular D1.
    pub(crate) fn value_at(&mut self, u_par: f64, f: &mut f64) -> bool {
        self.u = u_par;
        let (pc, mut d1c) = self.curve.d1(u_par);
        self.pc = pc;
        if Precision::is_infinite(d1c.x()) || Precision::is_infinite(d1c.y()) {
            *f = f64::INFINITY;
            return false;
        }
        let mut ndu = d1c.magnitude();

        if self.max_deriv_order != 0 && ndu <= self.tol {
            // Singular case (`hxx:167-241`).
            const DIVISION_FACTOR: f64 = 1.0e-3;
            let delta = (self.singular_du() * DIVISION_FACTOR).max(MIN_STEP);

            let mut n = 1;
            let mut v = GpVec::zero();
            let mut is_derive_found = false;
            while !is_derive_found && n < self.max_deriv_order {
                n += 1;
                v = self.curve.eval_dn(u_par, n);
                ndu = v.magnitude();
                is_derive_found = ndu > self.tol;
            }

            if is_derive_found {
                let u = if self.u - self.u_inf < delta {
                    self.u + delta
                } else {
                    self.u - delta
                };
                let p1 = self.curve.d0(self.u.min(u));
                let p2 = self.curve.d0(self.u.max(u));
                let v1 = GpVec::from_pnts(&p1, &p2);
                if v.dot(&v1) < 0.0 {
                    d1c = v.reversed();
                } else {
                    d1c = v;
                }
            } else {
                // Three-point Taylor (`hxx:213-238`); Ptemp default origin.
                let (p1, p2, p3, grown) = if self.u - self.u_inf < 2.0 * delta {
                    (
                        self.curve.d0(self.u),
                        self.curve.d0(self.u + delta),
                        self.curve.d0(self.u + 2.0 * delta),
                        true,
                    )
                } else {
                    (
                        self.curve.d0(self.u - 2.0 * delta),
                        self.curve.d0(self.u - delta),
                        self.curve.d0(self.u),
                        false,
                    )
                };
                let v1 = GpVec::new(p1.x(), p1.y(), p1.z());
                let v2 = GpVec::new(p2.x(), p2.y(), p2.z());
                let v3 = GpVec::new(p3.x(), p3.y(), p3.z());
                d1c = if grown {
                    v1.multiplied_scalar(-3.0)
                        .added(&v2.multiplied_scalar(4.0))
                        .subtracted(&v3)
                } else {
                    v1.subtracted(&v2.multiplied_scalar(4.0))
                        .added(&v3.multiplied_scalar(3.0))
                };
            }
            ndu = d1c.magnitude();
        }

        if ndu <= MIN_TOL {
            return false;
        }
        let ppc = GpVec::from_pnts(&self.point, &self.pc);
        *f = ppc.dot(&d1c) / ndu;
        true
    }

    /// `GFuncExtPC::Values` (`hxx:274-352`) including singular DF (`hxx:297-338`).
    pub(crate) fn values_at(&mut self, u_par: f64, f: &mut f64, df: &mut f64) -> bool {
        let pc_old = self.pc;
        let p_old = self.point;
        if !self.value_at(u_par, f) {
            return false;
        }
        self.u = u_par;
        self.pc = pc_old;
        self.point = p_old;
        let (pc, d1c, d2c) = self.curve.d2(u_par);
        self.pc = pc;
        let ndu = d1c.magnitude();
        if ndu <= self.tol {
            // Singular DF via three F samples (`hxx:297-338`).
            const DIVISION_FACTOR: f64 = 0.01;
            let delta = (self.singular_du() * DIVISION_FACTOR).max(MIN_STEP);
            let f1;
            let f2;
            let f3;
            if self.u - self.u_inf < 2.0 * delta {
                f1 = *f;
                let u2 = self.u + delta;
                let u3 = self.u + delta * 2.0;
                let mut f2v = 0.0;
                let mut f3v = 0.0;
                if !(self.value_at(u2, &mut f2v) && self.value_at(u3, &mut f3v)) {
                    return false;
                }
                f2 = f2v;
                f3 = f3v;
                *df = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * delta);
            } else {
                f3 = *f;
                let u1 = self.u - delta * 2.0;
                let u2 = self.u - delta;
                let mut f1v = 0.0;
                let mut f2v = 0.0;
                if !(self.value_at(u2, &mut f2v) && self.value_at(u1, &mut f1v)) {
                    return false;
                }
                f1 = f1v;
                f2 = f2v;
                *df = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * delta);
            }
            self.u = u_par;
            self.pc = pc_old;
            self.point = p_old;
        } else {
            let ppc = GpVec::from_pnts(&self.point, &self.pc);
            *df = ndu + (ppc.dot(&d2c) / ndu) - *f * (d1c.dot(&d2c)) / (ndu * ndu);
        }
        self.d1f = *df;
        true
    }
}

