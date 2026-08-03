//! IntTools_Curve data class — the face/face intersection-curve container.
//!
//! Ports `IntTools_Curve` (TKBO, `ModelingAlgorithms/TKBO/IntTools/IntTools_Curve.hxx`,
//! 216 lines). The class stores one 3D curve plus two optional 2D pcurves (one
//! per intersecting face) and two tolerance values. It is the payload Face/Face
//! intersection produces: the 3D curve is the intersection curve, the 2D
//! pcurves map the curve's parameter onto each face's UV space, and the
//! tolerances record the geometric validity of the result.
//!
//! Unlike the light `inttools_data::IntCurve` (a classification-only record),
//! this module carries the actual curve geometry, so it can evaluate points
//! ([`IntCurve::d0`]), classify the curve family ([`IntCurve::curve_type`]) and
//! report the parameter bounds ([`IntCurve::has_bounds`] / [`IntCurve::bounds`]).
//!
//! ponytail: `curve_type` classifies Line and Circle from geometric invariants;
//! Ellipse / Parabola / Hyperbola collapse to BSpline / Other (the `CurveKind`
//! doc folds Bezier into BSpline anyway). A root-exact conic classifier can
//! replace the sampling checks if a caller needs the fine-grained families.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::{Curve, GeomLine};
use occt_geom2d::curve::Curve2d;

use crate::inttools_data::CurveKind;

/// A face/face (or edge/face) intersection curve container.
///
/// Mirrors `IntTools_Curve`. `curve` is always present; the two pcurves are
/// optional because a valid intersection may predate pcurve computation.
#[derive(Clone)]
pub struct IntCurve {
    /// The 3D intersection curve.
    pub curve: Arc<dyn Curve>,
    /// Pcurve of `curve` on the first face (when known).
    pub pcurve1: Option<Arc<dyn Curve2d>>,
    /// Pcurve of `curve` on the second face (when known).
    pub pcurve2: Option<Arc<dyn Curve2d>>,
    /// Valid tolerance for the 3D curve (max deviation of the 3D curve from the
    /// 2D curves, or from the surfaces when no 2D curves exist).
    pub tolerance: f64,
    /// Tangential tolerance — max distance from the 3D curve to the end of the
    /// tangential zone between the faces.
    pub tangential_tolerance: f64,
}

/// Sentinel curve for the empty constructor. OCCT leaves the 3D handle null in
/// `IntTools_Curve()`; here a degenerate X-axis line stands in until
/// [`IntCurve::set_curves`] supplies real geometry.
fn empty_curve() -> Arc<dyn Curve> {
    Arc::new(GeomLine::from_pnt_dir(GpPnt::zero(), occt_core::gp::GpDir::default_dir()))
}

impl IntCurve {
    /// Empty constructor (OCCT `IntTools_Curve()`). Fill in the curves with
    /// [`set_curves`](Self::set_curves) before evaluating.
    pub fn new() -> Self {
        Self {
            curve: empty_curve(),
            pcurve1: None,
            pcurve2: None,
            tolerance: 0.0,
            tangential_tolerance: 0.0,
        }
    }

    /// Full constructor (OCCT `IntTools_Curve(curve, pc1, pc2, tol, tanTol)`).
    pub fn with_curves(
        curve: Arc<dyn Curve>,
        pcurve1: Option<Arc<dyn Curve2d>>,
        pcurve2: Option<Arc<dyn Curve2d>>,
        tolerance: f64,
        tangential_tolerance: f64,
    ) -> Self {
        Self { curve, pcurve1, pcurve2, tolerance, tangential_tolerance }
    }

    /// Sets the three curves at once (OCCT `SetCurves`).
    pub fn set_curves(
        &mut self,
        curve: Arc<dyn Curve>,
        pcurve1: Option<Arc<dyn Curve2d>>,
        pcurve2: Option<Arc<dyn Curve2d>>,
    ) {
        self.curve = curve;
        self.pcurve1 = pcurve1;
        self.pcurve2 = pcurve2;
    }

    /// Sets the 3D curve (OCCT `SetCurve`).
    pub fn set_curve(&mut self, curve: Arc<dyn Curve>) {
        self.curve = curve;
    }

    /// Sets the first 2D pcurve (OCCT `SetFirstCurve2d`).
    pub fn set_first_curve2d(&mut self, pc: Option<Arc<dyn Curve2d>>) {
        self.pcurve1 = pc;
    }

    /// Sets the second 2D pcurve (OCCT `SetSecondCurve2d`).
    pub fn set_second_curve2d(&mut self, pc: Option<Arc<dyn Curve2d>>) {
        self.pcurve2 = pc;
    }

    /// Sets the curve tolerance (OCCT `SetTolerance`).
    pub fn set_tolerance(&mut self, t: f64) {
        self.tolerance = t;
    }

    /// Sets the tangential tolerance (OCCT `SetTangentialTolerance`).
    pub fn set_tangential_tolerance(&mut self, t: f64) {
        self.tangential_tolerance = t;
    }

    /// Returns the 3D curve (OCCT `Curve()`).
    pub fn curve(&self) -> &Arc<dyn Curve> {
        &self.curve
    }

    /// Returns the first 2D pcurve (OCCT `FirstCurve2d()`).
    pub fn first_curve2d(&self) -> Option<&Arc<dyn Curve2d>> {
        self.pcurve1.as_ref()
    }

    /// Returns the second 2D pcurve (OCCT `SecondCurve2d()`).
    pub fn second_curve2d(&self) -> Option<&Arc<dyn Curve2d>> {
        self.pcurve2.as_ref()
    }

    /// Returns the curve tolerance (OCCT `Tolerance()`).
    pub fn get_tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns the tangential tolerance (OCCT `TangentialTolerance()`).
    pub fn get_tangential_tolerance(&self) -> f64 {
        self.tangential_tolerance
    }

    /// Whether the 3D curve is a bounded curve (OCCT `HasBounds()`).
    ///
    /// A curve is bounded when both its parameter extremes are finite; this is
    /// the Rust-port analogue of OCCT's `Geom_BoundedCurve` down-cast.
    pub fn has_bounds(&self) -> bool {
        let (a, b) = (self.curve.first_parameter(), self.curve.last_parameter());
        a.is_finite() && b.is_finite()
    }

    /// The boundary parameters `(first, last)` of the 3D curve, or `None` when
    /// the curve is unbounded (OCCT `Bounds`).
    pub fn bounds(&self) -> Option<(f64, f64)> {
        if self.has_bounds() {
            Some((self.curve.first_parameter(), self.curve.last_parameter()))
        } else {
            None
        }
    }

    /// Computes the 3D point at `param` (OCCT `D0`).
    ///
    /// `Err` when the curve is bounded and `param` lies outside `[first, last]`,
    /// or when `param` is non-finite. Unbounded curves accept any finite
    /// parameter.
    pub fn d0(&self, param: f64) -> Result<GpPnt, String> {
        if !param.is_finite() {
            return Err(format!("IntCurve::d0: non-finite parameter {param}"));
        }
        let (a, b) = (self.curve.first_parameter(), self.curve.last_parameter());
        if a.is_finite() && b.is_finite() && (param < a || param > b) {
            return Err(format!("IntCurve::d0: parameter {param} outside bounds [{a}, {b}]"));
        }
        Ok(self.curve.d0(param))
    }

    /// The analytic family of the 3D curve (OCCT `Type()`, mapped to
    /// [`CurveKind`]).
    ///
    /// Classified from geometric invariants by sampling: a constant-tangent
    /// curve is a [`CurveKind::Line`]; a coplanar, equidistant-sample curve is a
    /// [`CurveKind::Circle`]; any other bounded curve is [`CurveKind::BSpline`]
    /// (Bezier folded in, per the `CurveKind` doc) and any unbounded non-line is
    /// [`CurveKind::Other`].
    pub fn curve_type(&self) -> CurveKind {
        let c = &*self.curve;
        if is_line_like(c) {
            return CurveKind::Line;
        }
        if is_circle_like(c) {
            return CurveKind::Circle;
        }
        let (a, b) = (c.first_parameter(), c.last_parameter());
        if a.is_finite() && b.is_finite() {
            CurveKind::BSpline
        } else {
            CurveKind::Other
        }
    }
}

impl Default for IntCurve {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Geometric classification helpers (sampling-based)
// ---------------------------------------------------------------------------

/// Whether the curve has a constant tangent direction (a real line).
fn is_line_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let u0 = if a.is_finite() { a } else { 0.0 };
    let (_, d) = c.d1(u0);
    let m = d.magnitude();
    if m < 1e-12 {
        return false;
    }
    if a.is_finite() && b.is_finite() && b > a {
        for i in 1..=6 {
            let (_, di) = c.d1(a + (b - a) * i as f64 / 6.0);
            if di.cross_magnitude(&d) > 1e-6 * m * di.magnitude().max(m) {
                return false;
            }
        }
    }
    true
}

/// Whether every sample of the bounded curve is coplanar and equidistant from
/// a single center (a real circle).
fn is_circle_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite() && (b - a).abs() > 1e-15) {
        return false;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let Some((p0, p1, p2)) = first_three_spanning(&pts) else { return false };
    let Some(center) = circumcenter(&p0, &p1, &p2) else { return false };
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return false;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let nm = nrm.magnitude();
    if nm <= 1e-30 {
        return false;
    }
    let nrm = nrm.divided(nm);
    let scale = radius.max(1.0).max((b - a).abs());
    let tol = 1e-6 * scale;
    for p in pts {
        let v = GpVec::from_pnts(&center, &p);
        if v.dot(&nrm).abs() > tol {
            return false;
        }
        if (v.magnitude() - radius).abs() > tol {
            return false;
        }
    }
    true
}

/// First three samples that span a plane (non-collinear), if any.
fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
    let p0 = pts[0];
    let mut i1 = None;
    for (i, p) in pts.iter().enumerate().skip(1) {
        if GpVec::from_pnts(&p0, p).magnitude() > 1e-9 {
            i1 = Some(i);
            break;
        }
    }
    let i1 = i1?;
    let p1 = pts[i1];
    let d0 = GpVec::from_pnts(&p0, &p1);
    for p in pts.iter().skip(i1 + 1) {
        if GpVec::from_pnts(&p0, p).cross_magnitude(&d0) > 1e-9 * d0.magnitude().max(1e-9) {
            return Some((p0, p1, *p));
        }
    }
    None
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear points (perpendicular-bisector system).
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.crossed(&d2);
    if n.magnitude() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.xyz().x, n.xyz().y, n.xyz().z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n.xyz())];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpDir2d, GpPnt2d};
    use occt_geom::GeomCircle;
    use occt_geom2d::Geom2dLine;

    /// A pcurve for the test — a 2D line along X at the given offset.
    fn pc_x(off: f64) -> Arc<dyn Curve2d> {
        Arc::new(Geom2dLine::from_pnt_dir(
            GpPnt2d::new(off, 0.0),
            GpDir2d::new(1.0, 0.0).unwrap(),
        ))
    }

    #[test]
    fn intcurve_set_curves_line_d0_and_type() {
        let lin = GeomLine::from_pnt_dir(
            GpPnt::new(1.0, 2.0, 3.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        );
        let mut ic = IntCurve::new();
        ic.set_curves(Arc::new(lin), Some(pc_x(0.0)), None);
        ic.set_tolerance(1e-7);
        ic.set_tangential_tolerance(2e-7);

        assert_eq!(ic.get_tolerance(), 1e-7);
        assert_eq!(ic.get_tangential_tolerance(), 2e-7);
        assert!(ic.first_curve2d().is_some(), "first pcurve stored");
        assert!(ic.second_curve2d().is_none(), "second pcurve not set");
        assert_eq!(ic.curve_type(), CurveKind::Line);

        // An unbounded line has no bounds.
        assert!(!ic.has_bounds());
        assert!(ic.bounds().is_none());

        // D0 along the X line through (1,2,3).
        let p = ic.d0(0.5).unwrap();
        assert!((p.x() - 1.5).abs() < 1e-12, "x={}", p.x());
        assert!((p.y() - 2.0).abs() < 1e-12, "y={}", p.y());
        assert!((p.z() - 3.0).abs() < 1e-12, "z={}", p.z());
        // Unbounded → any finite parameter is accepted.
        let p = ic.d0(-2.0).unwrap();
        assert!((p.x() + 1.0).abs() < 1e-12, "x={}", p.x());
    }

    #[test]
    fn intcurve_bounded_circle_has_bounds() {
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let mut ic = IntCurve::new();
        ic.set_curve(Arc::new(circ));

        assert!(ic.has_bounds(), "a circle is a bounded curve");
        let (a, b) = ic.bounds().unwrap();
        assert!((a - 0.0).abs() < 1e-12, "first={a}");
        assert!((b - 2.0 * PI).abs() < 1e-12, "last={b}");
        assert_eq!(ic.curve_type(), CurveKind::Circle);

        let p = ic.d0(0.0).unwrap();
        assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-12, "p={p:?}");
        // Out-of-bounds and non-finite parameters error.
        assert!(ic.d0(3.0 * PI).is_err(), "beyond last parameter");
        assert!(ic.d0(-1.0).is_err(), "before first parameter");
        assert!(ic.d0(f64::NAN).is_err(), "NaN rejected");
    }

    #[test]
    fn intcurve_full_constructor_and_accessors() {
        let lin = GeomLine::from_pnt_dir(GpPnt::zero(), GpDir::default_dir());
        let ic = IntCurve::with_curves(Arc::new(lin), Some(pc_x(1.0)), Some(pc_x(2.0)), 1e-6, 1e-5);
        assert_eq!(ic.get_tolerance(), 1e-6);
        assert_eq!(ic.get_tangential_tolerance(), 1e-5);
        assert!(ic.curve().first_parameter().is_infinite(), "unbounded line");
        assert!(ic.first_curve2d().is_some());
        assert!(ic.second_curve2d().is_some());

        // Default → sentinel line, then per-field setters.
        let mut d = IntCurve::default();
        assert_eq!(d.curve_type(), CurveKind::Line, "sentinel is a line");
        assert!(!d.has_bounds());
        d.set_first_curve2d(Some(pc_x(0.5)));
        d.set_second_curve2d(None);
        assert!(d.first_curve2d().is_some());
        assert!(d.second_curve2d().is_none());
        d.set_tolerance(0.25);
        assert_eq!(d.get_tolerance(), 0.25);
    }

    #[test]
    fn intcurve_all_curve_kinds_classify() {
        // A bounded quadratic Bezier is neither line nor circle → BSpline
        // (Bezier folded into BSpline per the CurveKind doc).
        let bs = Arc::new(
            occt_geom::GeomBSplineCurve::new(
                vec![
                    GpPnt::new(0.0, 0.0, 0.0),
                    GpPnt::new(1.0, 1.0, 0.0),
                    GpPnt::new(2.0, 0.0, 0.0),
                ],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                2,
            )
            .unwrap(),
        );
        let mut ic = IntCurve::new();
        ic.set_curve(bs);
        assert_eq!(ic.curve_type(), CurveKind::BSpline, "quadratic Bezier is not line/circle");
        assert!(ic.has_bounds());
        assert!(ic.d0(0.5).is_ok());
    }
}
