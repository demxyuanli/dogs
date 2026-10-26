//! `Adaptor2d_Curve2d::Trim` as `Geom2dAdaptor_Curve` reports it.
//!
//! Source: `ModelingData/TKG2d/Geom2dAdaptor/Geom2dAdaptor_Curve.cxx`
//! (`Trim` `:577-584`, the `(Curve, UFirst, ULast)` constructor `:250-258`,
//! `NbIntervals` `:409-486`, `Intervals` `:490-...`, `Continuity`
//! `:...`).
//!
//! `Geom2dAdaptor_Curve::Trim(First, Last, Tol)` is
//! `new Geom2dAdaptor_Curve(myCurve, First, Last)`: the stored curve is kept
//! **unchanged** and only `myFirst` / `myLast` move to the requested range, so
//! evaluation stays at the basis parameter (`Value(U)` forwards to the curve)
//! while the span queries are reported over `[First, Last]`. This wrapper is
//! exactly that: same evaluation, same type table, new range.
//!
//! Deviations, both outside the queries `Approx_CurveOnSurface` reads:
//! - `Continuity` is the B-spline local continuity over the new range when the
//!   wrapped curve is a B-spline, and the wrapped curve's answer otherwise
//!   (OCCT's `Geom2dAdaptor_Curve::Continuity` also maps the offset-curve arm
//!   to the basis, which the wrapped curve already reports).
//! - `NbIntervals` / `Intervals` follow `Geom2dAdaptor_Curve.cxx:409-486`:
//!   B-spline knots restricted to `[myFirst, myLast]`, an offset curve recursed
//!   into its basis with the continuity shifted by one (C0->C1, C1->C2, C2->C3,
//!   else CN), and a single span for every other type.

use std::sync::Arc;

use occt_core::bspl::adaptor_intervals;
use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::PCONFUSION;

use crate::curve::Curve2d;
use crate::Geom2dBSplineCurve;

/// `Geom2dAdaptor_Curve(Curve, First, Last)` — the result of `Trim`.
#[derive(Clone)]
pub struct TrimmedAdaptor2d {
    curve: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
}

impl TrimmedAdaptor2d {
    /// `Geom2dAdaptor_Curve(C, UFirst, ULast)` (`cxx:250-258`): `load` keeps a
    /// `Geom2d_TrimmedCurve`'s basis in the type table and stores the two
    /// parameters as the range.
    pub fn new(curve: Arc<dyn Curve2d>, first: f64, last: f64) -> Self {
        Self { curve, first, last }
    }

    pub fn curve(&self) -> &Arc<dyn Curve2d> {
        &self.curve
    }

    fn shift_continuity(continuity: u8) -> Option<u8> {
        // `cxx:456-474`: C0->C1, C1->C2, C2->C3, CN->CN; G1/G2 raise
        // `Standard_DomainError` (`cxx:459-461`).
        match continuity {
            0 => Some(2),
            2 => Some(4),
            4 => Some(5),
            6 => Some(6),
            _ => None,
        }
    }
}

impl Curve2d for TrimmedAdaptor2d {
    // `Geom2dAdaptor_Curve::Value` / `D0` / `D1` / `D2` evaluate the stored
    // curve at the requested parameter: the range only bounds the queries.
    fn d0(&self, u: f64) -> GpPnt2d {
        self.curve.d0(u)
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        self.curve.d1(u)
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        self.curve.d2(u)
    }
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        self.curve.d3(u)
    }
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        self.curve.eval_dn(u, n)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }
    fn period(&self) -> f64 {
        self.curve.period()
    }

    /// `Geom2dAdaptor_Curve::Continuity` (`cxx:...`): the local continuity of a
    /// B-spline over `[myFirst, myLast]`, the stored curve's answer otherwise.
    fn continuity(&self) -> u8 {
        if let (Some(knots), Some(degree)) = (self.curve.bspline_knots2d(), self.curve.bspline_degree())
        {
            return occt_core::bspl::local_continuity(
                knots,
                degree,
                self.curve.is_periodic(),
                self.first,
                self.last,
            );
        }
        self.curve.continuity()
    }

    fn transform(&mut self, t: &GpTrsf2d) {
        let mut c = self.curve.clone_dyn();
        c.transform(t);
        self.curve = Arc::from(c);
    }

    fn reverse(&mut self) {
        let mut c = self.curve.clone_dyn();
        c.reverse();
        self.curve = Arc::from(c);
    }

    // Type table: every query resolves on the stored curve, whose own trimmed
    // unwrapping (`Geom2dAdaptor_Curve::load`, `cxx:279-346`) is already applied.
    fn poles2d(&self) -> Option<Vec<GpPnt2d>> {
        self.curve.poles2d()
    }
    fn set_poles2d(&mut self, poles: &[GpPnt2d]) {
        let mut c = self.curve.clone_dyn();
        c.set_poles2d(poles);
        self.curve = Arc::from(c);
    }
    fn is_bspline2d(&self) -> bool {
        self.curve.is_bspline2d()
    }
    fn is_geom2d_line(&self) -> bool {
        self.curve.is_geom2d_line()
    }
    fn is_geom2d_circle(&self) -> bool {
        self.curve.is_geom2d_circle()
    }
    fn is_bezier2d(&self) -> bool {
        self.curve.is_bezier2d()
    }
    fn is_line(&self) -> bool {
        self.curve.is_line()
    }
    fn gp_lin2d(&self) -> Option<GpLin2d> {
        self.curve.gp_lin2d()
    }
    fn gp_circ2d(&self) -> Option<GpCirc2d> {
        self.curve.gp_circ2d()
    }
    fn gp_elips2d(&self) -> Option<GpElips2d> {
        self.curve.gp_elips2d()
    }
    fn gp_parab2d(&self) -> Option<GpParab2d> {
        self.curve.gp_parab2d()
    }
    fn gp_hypr2d(&self) -> Option<GpHypr2d> {
        self.curve.gp_hypr2d()
    }
    fn trimmed_basis(&self) -> Option<&dyn Curve2d> {
        self.curve.trimmed_basis()
    }
    fn offset_basis(&self) -> Option<&dyn Curve2d> {
        self.curve.offset_basis()
    }
    fn offset_value(&self) -> Option<f64> {
        self.curve.offset_value()
    }
    fn bspline_weights2d(&self) -> Option<&[f64]> {
        self.curve.bspline_weights2d()
    }
    fn bezier_weights2d(&self) -> Option<&[f64]> {
        self.curve.bezier_weights2d()
    }
    fn bspline_distinct_knots_mults(&self) -> Option<(Vec<f64>, Vec<i32>)> {
        self.curve.bspline_distinct_knots_mults()
    }
    fn bezier_nb_poles(&self) -> Option<usize> {
        self.curve.bezier_nb_poles()
    }
    fn is_rational(&self) -> bool {
        self.curve.is_rational()
    }
    fn bspline_nb_knots(&self) -> Option<usize> {
        self.curve.bspline_nb_knots()
    }
    fn bspline_degree(&self) -> Option<usize> {
        self.curve.bspline_degree()
    }
    fn bspline_poles2d(&self) -> Option<(&[f64], &[f64])> {
        self.curve.bspline_poles2d()
    }
    fn bspline_copy2d(&self) -> Option<Geom2dBSplineCurve> {
        self.curve.bspline_copy2d()
    }
    fn bspline_knots2d(&self) -> Option<&[f64]> {
        self.curve.bspline_knots2d()
    }

    /// `Geom2dAdaptor_Curve::Resolution` (`cxx:1186-1219`) switches on the type
    /// only, never on the range, so the stored curve's answer is the same.
    fn resolution(&self, ruv: f64) -> f64 {
        self.curve.resolution(ruv)
    }

    /// `Geom2dAdaptor_Curve::NbIntervals` / `Intervals`
    /// (`cxx:409-486`, `:490-...`).
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        if let (Some(knots), Some(degree)) =
            (self.curve.bspline_knots2d(), self.curve.bspline_degree())
        {
            // `cxx:411-450`: `BSplCLib::Intervals(Knots, Mults, Degree,
            // IsPeriodic, aCont, myFirst, myLast, anEps, ...)`.
            let eps = self.curve.resolution(occt_core::precision::CONFUSION).min(PCONFUSION);
            return adaptor_intervals(
                knots,
                degree,
                self.curve.is_periodic(),
                continuity,
                self.first,
                self.last,
                eps,
            );
        }
        if let Some(basis) = self.curve.offset_basis() {
            // `cxx:453-480`: recurse into the basis with the continuity shifted
            // by one and the same range. `G1` / `G2` raise
            // `Standard_DomainError` (`cxx:459-461`); the port reports a single
            // span there, which is what OCCT's caller would get for C0.
            if let Some(base_s) = Self::shift_continuity(continuity) {
                let inner = TrimmedAdaptor2d::new(Arc::from(basis.clone_dyn()), self.first, self.last);
                return inner.parameter_intervals(base_s);
            }
            return vec![self.first, self.last];
        }
        // `cxx:482-485`: every other type is a single interval.
        vec![self.first, self.last]
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
}
