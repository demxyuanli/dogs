//! Trimmed 2D curve. Source: `Geom2d_TrimmedCurve.cxx`.
use std::sync::Arc;
use occt_core::gp::{GpPnt2d, GpVec2d, GpTrsf2d};
use occt_core::precision::{PCONFUSION, Precision};
use crate::curve::Curve2d;

/// `ElCLib::AdjustPeriodic` (`ElCLib.cxx:115`).
fn adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    if Precision::is_infinite(u_first) || Precision::is_infinite(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    if a_period < f64::EPSILON {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}

#[derive(Clone)]
pub struct Geom2dTrimmedCurve {
    basis: Arc<dyn Curve2d>,
    u_trim1: f64,
    u_trim2: f64,
}

impl Geom2dTrimmedCurve {
    pub fn new(curve: Arc<dyn Curve2d>, u1: f64, u2: f64) -> Self {
        Self::new_sense(curve, u1, u2, true, true)
    }

    /// `Geom2d_TrimmedCurve(C, U1, U2, Sense, theAdjustPeriodic)`.
    pub fn new_sense(
        curve: Arc<dyn Curve2d>,
        u1: f64,
        u2: f64,
        sense: bool,
        adjust_periodic: bool,
    ) -> Self {
        let basis = if let Some(b) = curve.trimmed_basis() {
            Arc::from(b.clone_dyn())
        } else {
            curve
        };
        let mut s = Self {
            basis,
            u_trim1: u1,
            u_trim2: u2,
        };
        s.set_trim(u1, u2, sense, adjust_periodic);
        s
    }

    pub fn basis_curve(&self) -> &Arc<dyn Curve2d> {
        &self.basis
    }

    /// `Geom2d_TrimmedCurve::SetTrim` (`cxx:98-154`).
    fn set_trim(&mut self, u1: f64, u2: f64, sense: bool, adjust_periodic_flag: bool) {
        if u1 == u2 {
            return;
        }
        let udeb = self.basis.first_parameter();
        let ufin = self.basis.last_parameter();
        let mut same_sense = true;
        if self.basis.is_periodic() {
            same_sense = sense;
            self.u_trim1 = u1;
            self.u_trim2 = u2;
            if adjust_periodic_flag {
                let preci = ((self.u_trim2 - self.u_trim1).abs() / 2.0).min(PCONFUSION);
                adjust_periodic(udeb, ufin, preci, &mut self.u_trim1, &mut self.u_trim2);
            }
        } else if u1 < u2 {
            same_sense = sense;
            self.u_trim1 = u1;
            self.u_trim2 = u2;
        } else {
            same_sense = !sense;
            self.u_trim1 = u2;
            self.u_trim2 = u1;
        }
        if !same_sense {
            self.reverse();
        }
    }
}

impl Curve2d for Geom2dTrimmedCurve {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.basis.d0(u)
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        self.basis.d1(u)
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        self.basis.d2(u)
    }
    /// `Geom2d_TrimmedCurve::EvalD3` (`Geom2d_TrimmedCurve.cxx:273-276`): direct
    /// delegation to the basis, like `d0`/`d1`/`d2` above.
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        self.basis.d3(u)
    }
    /// `Geom2d_TrimmedCurve::EvalDN` (`Geom2d_TrimmedCurve.cxx:280-283`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        self.basis.eval_dn(u, n)
    }
    fn first_parameter(&self) -> f64 {
        self.u_trim1
    }
    fn last_parameter(&self) -> f64 {
        self.u_trim2
    }
    fn is_periodic(&self) -> bool {
        if self.basis.is_periodic() {
            let period = self.basis.period();
            let length = self.u_trim2 - self.u_trim1;
            if length > PCONFUSION && period > 0.0 {
                let rem = length - period * (length / period).round();
                if rem.abs() <= PCONFUSION {
                    return true;
                }
            }
        }
        false
    }
    fn period(&self) -> f64 {
        self.basis.period()
    }
    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }
    fn transform(&mut self, t: &GpTrsf2d) {
        let mut b = self.basis.clone_dyn();
        b.transform(t);
        self.basis = Arc::from(b);
    }
    fn reverse(&mut self) {
        let mut b = self.basis.clone_dyn();
        b.reverse();
        self.basis = Arc::from(b);
        std::mem::swap(&mut self.u_trim1, &mut self.u_trim2);
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
    /// `Geom2d_TrimmedCurve::TransformedParameter` (`Geom2d_TrimmedCurve.cxx:297-300`).
    fn transformed_parameter(&self, u: f64, t: &GpTrsf2d) -> f64 {
        self.basis.transformed_parameter(u, t)
    }
    fn is_line(&self) -> bool {
        self.basis.is_line()
    }
    fn gp_lin2d(&self) -> Option<occt_core::gp::GpLin2d> {
        self.basis.gp_lin2d()
    }
    fn gp_circ2d(&self) -> Option<occt_core::gp::GpCirc2d> {
        self.basis.gp_circ2d()
    }
    /// `Geom2dAdaptor_Curve::load` (`cxx:96-120`) unwraps a
    /// `Geom2d_TrimmedCurve` and keeps the **basis**, so the adaptor's
    /// `GetType()` / `Ellipse()` / `Parabola()` / `Hyperbola()` queries resolve
    /// to the basis.
    fn gp_elips2d(&self) -> Option<occt_core::gp::GpElips2d> {
        self.basis.gp_elips2d()
    }
    fn gp_parab2d(&self) -> Option<occt_core::gp::GpParab2d> {
        self.basis.gp_parab2d()
    }
    fn gp_hypr2d(&self) -> Option<occt_core::gp::GpHypr2d> {
        self.basis.gp_hypr2d()
    }
    fn trimmed_basis(&self) -> Option<&dyn Curve2d> {
        Some(self.basis.as_ref())
    }
}
