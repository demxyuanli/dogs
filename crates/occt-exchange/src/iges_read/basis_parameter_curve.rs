//! Trim that keeps the basis parameter, unlike `GeomTrimmedCurveBasis`.
//!
//! `Geom_TrimmedCurve` stores basis parameters (`Geom_TrimmedCurve.cxx:255-267`).
//! `GeomTrimmedCurveBasis::untrimmed_basis` is read by `box_curve` as a `[0, 1]`
//! fraction (`geom_bnd_lib_curve3d.rs:538-540`), which stretches a line of length
//! L out to L*L. This view reports the same parameters and does not advertise
//! that remap. `trim_basis` still unwraps it for `ShapeFix_Wire.cxx:1686-1697`.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCirc, GpLin, GpPnt, GpTrsf, GpVec};
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;

use super::iges_to_brep_basic_curve::CurveSpan;

/// Rotate a transferred generatrix. Rigid motion keeps the edge parameters
/// (`Geom_Curve::Transform` / `TransformedParameter` for isometries).
pub(in crate::iges_read) fn rotate_span(
    span: &CurveSpan,
    axis: &GpAx1,
    angle: f64,
) -> Result<Arc<dyn Curve>, String> {
    let mut rot = GpTrsf::identity();
    rot.set_rotation_ax1(axis, angle)
        .map_err(|e| e.to_string())?;
    let (basis, first, last) = if let Some((b, f, l)) = span.curve.untrimmed_basis() {
        (b, f, l)
    } else {
        (span.curve.clone(), span.first, span.last)
    };
    let moved: Arc<dyn Curve> = Arc::from(basis.transformed(&rot));
    let cf = moved.first_parameter();
    let cl = moved.last_parameter();
    if first != cf || last != cl {
        Ok(Arc::new(BasisParameterCurve::new(moved, first, last)))
    } else {
        Ok(moved)
    }
}

/// Trim evaluated in the basis parameter, without `untrimmed_basis`.
struct BasisParameterCurve {
    basis: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl BasisParameterCurve {
    fn new(basis: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { basis, first, last }
    }
}

impl Curve for BasisParameterCurve {
    fn d0(&self, u: f64) -> GpPnt {
        self.basis.d0(u)
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.basis.d1(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.basis.d2(u)
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }
    /// `Geom_TrimmedCurve::IsPeriodic` (`Geom_TrimmedCurve.cxx:172-185`): a
    /// periodic basis counts only when the trim length is a whole period.
    fn is_periodic(&self) -> bool {
        if !self.basis.is_periodic() {
            return false;
        }
        let period = self.basis.period();
        if period <= 0.0 {
            return false;
        }
        let length = self.last - self.first;
        if length <= PCONFUSION {
            return false;
        }
        let rem = length - period * (length / period).round();
        rem.abs() <= PCONFUSION
    }
    fn circle_radius(&self) -> Option<f64> {
        self.basis.circle_radius()
    }
    fn trim_basis(&self) -> Option<Arc<dyn Curve>> {
        Some(self.basis.clone())
    }
    /// Forward `Geom_OffsetCurve` so `ShapeFix_Wire::FixShifted` (`cxx:1686-1697`)
    /// can unwrap an offset generatrix after this trim view.
    fn offset_curve(&self) -> Option<(Arc<dyn Curve>, f64)> {
        self.basis.offset_curve()
    }
    fn period(&self) -> f64 {
        self.basis.period()
    }
    /// `GeomAdaptor_Curve::load` unwraps `Geom_TrimmedCurve` (`cxx:252-254`)
    /// and `Resolution` (`cxx:1116-1148`) uses that basis type. The stored
    /// parameters are still the basis parameters, so the value is not rescaled.
    fn resolution(&self, r3d: f64) -> f64 {
        self.basis.resolution(r3d)
    }
    fn is_line(&self) -> bool {
        self.basis.is_line()
    }
    fn gp_line(&self) -> Option<GpLin> {
        self.basis.gp_line()
    }
    fn gp_circ(&self) -> Option<GpCirc> {
        self.basis.gp_circ()
    }
    fn transform(&mut self, t: &GpTrsf) {
        let mut curve = self.basis.clone_dyn();
        curve.transform(t);
        let scale = t.scale_factor().abs();
        if self.basis.is_line() && (scale - 1.0).abs() > 0.0 {
            self.first *= scale;
            self.last *= scale;
        }
        self.basis = Arc::from(curve);
    }
    fn reverse(&mut self) {
        std::mem::swap(&mut self.first, &mut self.last);
    }
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(Self {
            basis: self.basis.clone(),
            first: self.first,
            last: self.last,
        })
    }
}
