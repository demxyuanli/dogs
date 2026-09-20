//! Trimmed 3D curve. Source: `Geom_TrimmedCurve.hxx`
use std::sync::Arc;
use occt_core::gp::{GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;

#[derive(Clone)]
pub struct GeomTrimmedCurve {
    basis: Arc<dyn Curve>,
    first: f64,
    last: f64,
    /// Flat knots remapped into the trimmed `[0, 1]` parameter (BSpline only).
    remapped_knots: Option<Vec<f64>>,
    nurbs_deg: Option<usize>,
}

impl GeomTrimmedCurve {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        let a = first.min(last);
        let b = first.max(last);
        let den = b - a;
        let remapped_knots = curve.bspline_knots().map(|knots| {
            if den.abs() <= 0.0 {
                knots.to_vec()
            } else {
                knots
                    .iter()
                    .map(|&u| (u - a) / den)
                    .collect()
            }
        });
        let nurbs_deg = curve.nurbs_degree();
        Self {
            basis: curve,
            first: a,
            last: b,
            remapped_knots,
            nurbs_deg,
        }
    }
    pub fn basis_curve(&self) -> &Arc<dyn Curve> { &self.basis }
}

impl Curve for GeomTrimmedCurve {
    fn d0(&self, u: f64) -> GpPnt { self.basis.d0(self.first + u * (self.last - self.first)) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let t = self.first + u * (self.last - self.first);
        let (p, d) = self.basis.d1(t);
        let s = self.last - self.first;
        (p, GpVec::new(d.x()*s, d.y()*s, d.z()*s))
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let t = self.first + u * (self.last - self.first);
        let (p, d1, d2) = self.basis.d2(t);
        let s = self.last - self.first;
        (p, GpVec::new(d1.x()*s,d1.y()*s,d1.z()*s), GpVec::new(d2.x()*s*s,d2.y()*s*s,d2.z()*s*s))
    }
    /// `Geom_TrimmedCurve::EvalD3` (`Geom_TrimmedCurve.cxx:231-236`) delegates to
    /// the basis; this view reparameterises to `[0, 1]`, so the chain rule of
    /// `d1`/`d2` above applies with `s = last - first`.
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let t = self.first + u * (self.last - self.first);
        let (p, d1, d2, d3) = self.basis.d3(t);
        let s = self.last - self.first;
        let s2 = s * s;
        let s3 = s2 * s;
        (
            p,
            GpVec::new(d1.x()*s, d1.y()*s, d1.z()*s),
            GpVec::new(d2.x()*s2, d2.y()*s2, d2.z()*s2),
            GpVec::new(d3.x()*s3, d3.y()*s3, d3.z()*s3),
        )
    }
    /// `Geom_TrimmedCurve::EvalDN` (`Geom_TrimmedCurve.cxx:240-243`) delegates to
    /// the basis; same `s^N` chain rule as `d1`/`d2`/`d3`.
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        let t = self.first + u * (self.last - self.first);
        let s = self.last - self.first;
        self.basis.eval_dn(t, n).multiplied_scalar(s.powi(n))
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 1.0 }
    fn continuity(&self) -> u8 { self.basis.continuity() }
    fn circle_radius(&self) -> Option<f64> { self.basis.circle_radius() }
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> { self.basis.gp_circ() }
    /// `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps a `Geom_TrimmedCurve`
    /// and keeps the **basis**, so the adaptor's `GetType()` / `Line()` /
    /// `Ellipse()` / `Hyperbola()` / `Parabola()` queries resolve to the basis.
    fn gp_line(&self) -> Option<occt_core::gp::GpLin> { self.basis.gp_line() }
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> { self.basis.gp_ellipse() }
    fn gp_hyperbola(&self) -> Option<occt_core::gp::GpHypr> { self.basis.gp_hyperbola() }
    fn gp_parabola(&self) -> Option<occt_core::gp::GpParab> { self.basis.gp_parabola() }
    fn is_geom_trimmed(&self) -> bool { true }
    fn trimmed_basis_range(&self) -> Option<(f64, f64)> {
        // `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps nested trims.
        if let Some((bf, bl)) = self.basis.trimmed_basis_range() {
            Some((
                bf + self.first * (bl - bf),
                bf + self.last * (bl - bf),
            ))
        } else {
            Some((self.first, self.last))
        }
    }
    fn untrimmed_basis(&self) -> Option<(std::sync::Arc<dyn Curve>, f64, f64)> {
        // `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps to the basis curve
        // while the adaptor range stays in the caller's parameter space. Our
        // `Geom_TrimmedCurve` remaps the trim onto `[0, 1]`, so the reported
        // range is expressed in the basis curve's own parameters.
        match self.basis.untrimmed_basis() {
            Some((innermost, ib1, ib2)) => {
                let den = ib2 - ib1;
                Some((
                    innermost,
                    ib1 + self.first * den,
                    ib1 + self.last * den,
                ))
            }
            None => Some((self.basis.clone(), self.first, self.last)),
        }
    }
    fn is_line(&self) -> bool { self.basis.is_line() }
    fn nurbs_degree(&self) -> Option<usize> { self.nurbs_deg.or_else(|| self.basis.nurbs_degree()) }
    fn bspline_knots(&self) -> Option<&[f64]> { self.remapped_knots.as_deref() }
    fn bspline_poles(&self) -> Option<&[GpPnt]> { self.basis.bspline_poles() }
    fn bspline_weights(&self) -> Option<&[f64]> { self.basis.bspline_weights() }
    fn bezier_poles(&self) -> Option<&[GpPnt]> { self.basis.bezier_poles() }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        let den = self.last - self.first;
        if den.abs() <= 0.0 {
            return vec![0.0, 1.0];
        }
        self.basis
            .parameter_intervals(continuity)
            .into_iter()
            .map(|u| (u - self.first) / den)
            .collect()
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }
    fn resolution(&self, r3d: f64) -> f64 {
        let den = (self.last - self.first).abs();
        if den <= 0.0 {
            return self.basis.resolution(r3d);
        }
        self.basis.resolution(r3d) / den
    }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn reverse(&mut self) { std::mem::swap(&mut self.first, &mut self.last); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}

/// `Geom_TrimmedCurve` evaluated in the BASIS curve's parameters, i.e. the
/// faithful `Geom_Curve` view, as opposed to [`GeomTrimmedCurve`] above whose
/// `[0, 1]` remap matches `GeomAdaptor_Curve` consumers in this port.
///
/// Source: `Geom_TrimmedCurve.cxx:212-215` (`EvalD0` / `EvalD1` / `EvalD2`
/// forward to `basisCurve` with the same parameter), `:255-267`
/// (`FirstParameter` / `LastParameter` return the trimming parameters
/// `uTrim1` / `uTrim2`), `:155-167` (`IsClosed`) and `:169-180` (`IsPeriodic`,
/// true only when the trim covers a full period).
#[derive(Clone)]
pub struct GeomTrimmedCurveBasis {
    basis: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl GeomTrimmedCurveBasis {
    /// `Geom_TrimmedCurve(theBasisCurve, theU1, theU2)`; `theU1` / `theU2` are
    /// basis parameters (`Geom_TrimmedCurve.cxx:60-80`).
    pub fn new(basis: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { basis, first, last }
    }

    pub fn basis_curve(&self) -> &Arc<dyn Curve> { &self.basis }
}

impl Curve for GeomTrimmedCurveBasis {
    fn d0(&self, u: f64) -> GpPnt { self.basis.d0(u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { self.basis.d1(u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { self.basis.d2(u) }
    /// `Geom_TrimmedCurve::EvalD3` (`Geom_TrimmedCurve.cxx:231-236`): direct
    /// delegation — this view is already expressed in basis parameters.
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) { self.basis.d3(u) }
    /// `Geom_TrimmedCurve::EvalDN` (`Geom_TrimmedCurve.cxx:240-243`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec { self.basis.eval_dn(u, n) }
    fn first_parameter(&self) -> f64 { self.first }
    fn last_parameter(&self) -> f64 { self.last }
    /// `Geom_TrimmedCurve::IsPeriodic` (`cxx:169-180`): the basis must be
    /// periodic and the trim length a multiple of the period.
    fn is_periodic(&self) -> bool {
        if !self.basis.is_periodic() {
            return false;
        }
        let period = self.basis.period();
        if period <= 0.0 {
            return false;
        }
        let length = self.last - self.first;
        if length <= occt_core::precision::PCONFUSION {
            return false;
        }
        let rem = length - period * (length / period).round();
        rem.abs() <= occt_core::precision::PCONFUSION
    }
    fn period(&self) -> f64 { self.basis.period() }
    fn continuity(&self) -> u8 { self.basis.continuity() }
    fn circle_radius(&self) -> Option<f64> { self.basis.circle_radius() }
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> { self.basis.gp_circ() }
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> { self.basis.gp_ellipse() }
    fn is_geom_trimmed(&self) -> bool { true }
    fn trimmed_basis_range(&self) -> Option<(f64, f64)> {
        // `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps nested trims. This
        // view already reports basis parameters, so no remap is applied.
        if let Some((bf, bl)) = self.basis.trimmed_basis_range() {
            Some((bf, bl))
        } else {
            Some((self.first, self.last))
        }
    }
    fn untrimmed_basis(&self) -> Option<(Arc<dyn Curve>, f64, f64)> {
        match self.basis.untrimmed_basis() {
            Some((innermost, ib1, ib2)) => Some((innermost, ib1, ib2)),
            None => Some((self.basis.clone(), self.first, self.last)),
        }
    }
    fn is_line(&self) -> bool { self.basis.is_line() }
    fn nurbs_degree(&self) -> Option<usize> { self.basis.nurbs_degree() }
    fn bspline_knots(&self) -> Option<&[f64]> { self.basis.bspline_knots() }
    fn bspline_poles(&self) -> Option<&[GpPnt]> { self.basis.bspline_poles() }
    fn bspline_weights(&self) -> Option<&[f64]> { self.basis.bspline_weights() }
    fn bezier_poles(&self) -> Option<&[GpPnt]> { self.basis.bezier_poles() }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.parameter_intervals(continuity)
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }
    fn resolution(&self, r3d: f64) -> f64 { self.basis.resolution(r3d) }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn reverse(&mut self) { std::mem::swap(&mut self.first, &mut self.last); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
