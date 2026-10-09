//! `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter` and
//! `IntCurveSurface_TheQuadCurvExactHInter`.
//!
//! Source: `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter.hxx/.cxx`
//! (instantiation of `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter`)
//! and `IntCurveSurface_TheQuadCurvExactHInter.hxx/.cxx` +
//! `IntCurveSurface_QuadricCurveExactInterUtils.pxx`.
//!
//! This is the `InternalPerformCurveQuadric` arm of the OCCT
//! `IntCurveSurface_HInter` dispatch: a (generally non-conic) curve against a
//! plane / cylinder / cone / sphere. `math_FunctionAllRoots` finds every root
//! of the signed distance `Q(w)` on each C1 interval of the curve.

use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_all_roots::FunctionAllRoots;
use occt_core::math_function_sample::FunctionSample;
use occt_geom::{Curve, Surface};

use crate::int_surf_quadric::IntSurfQuadric;

use super::polygon_utils::nb_samples;

/// `EPSX` (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:29`).
pub const EPSX: f64 = 0.00000000000001;
/// `EPSDIST` (`...pxx:30`).
pub const EPSDIST: f64 = 0.00000001;
/// `EPSNUL` (`...pxx:31`).
pub const EPSNUL: f64 = 0.00000001;

/// The port's `Adaptor3d_Surface::GetType()` dispatch + `SetValue`
/// (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:66-87`):
/// Plane / Cylinder / Cone / Sphere build an [`IntSurfQuadric`]; every other
/// surface type (including Torus) hits `default: return`.
pub fn int_surf_quadric_of(surface: &dyn Surface) -> Option<IntSurfQuadric> {
    if let Some(p) = surface.gp_pln() {
        return Some(IntSurfQuadric::from_plane(&p));
    }
    if let Some(c) = surface.gp_cylinder() {
        return Some(IntSurfQuadric::from_cylinder(&c));
    }
    if let Some(c) = surface.gp_cone() {
        return Some(IntSurfQuadric::from_cone(&c));
    }
    if let Some(s) = surface.gp_sphere() {
        return Some(IntSurfQuadric::from_sphere(&s));
    }
    None
}

/// `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter`
/// (`...hxx:27-56`, `.cxx:24-55`): `math_FunctionWithDerivative` returning the
/// signed distance `Q(w)` and `T(w) . grad Q`.
pub struct TheQuadCurvFuncOfTheQuadCurvExactHInter<'a> {
    quadric: IntSurfQuadric,
    curve: &'a dyn Curve,
}

impl<'a> TheQuadCurvFuncOfTheQuadCurvExactHInter<'a> {
    /// `...OfTheQuadCurvExactHInter(Q, C)` (`.cxx:24-30`).
    pub fn new(quadric: IntSurfQuadric, curve: &'a dyn Curve) -> Self {
        Self { quadric, curve }
    }
}

impl MathFunction for TheQuadCurvFuncOfTheQuadCurvExactHInter<'_> {
    /// `Value(Param, F)` (`.cxx:32-36`):
    /// `F = myQuadric.Distance(CurveTool::Value(myCurve, Param))`.
    fn value(&mut self, param: f64, f: &mut f64) -> bool {
        *f = self.quadric.distance(&self.curve.d0(param));
        true
    }
}

impl MathFunctionWithDerivative for TheQuadCurvFuncOfTheQuadCurvExactHInter<'_> {
    /// `Derivative(Param, D)` (`.cxx:38-45`):
    /// `D = T.Dot(myQuadric.Gradient(P))`.
    fn derivative(&mut self, param: f64, d: &mut f64) -> bool {
        let (p, t) = self.curve.d1(param);
        *d = t.dot(&self.quadric.gradient(&p));
        true
    }

    /// `Values(Param, F, D)` (`.cxx:47-56`): `ValAndGrad(P, F, Grad)` then
    /// `D = T.Dot(Grad)`.
    fn values(&mut self, param: f64, f: &mut f64, d: &mut f64) -> bool {
        let (p, t) = self.curve.d1(param);
        let (dist, grad) = self.quadric.val_and_grad(&p);
        *f = dist;
        *d = t.dot(&grad);
        true
    }
}

/// `IntCurveSurface_TheQuadCurvExactHInter` (`...hxx:28-53`, `.cxx:28-73`).
pub struct TheQuadCurvExactHInter {
    nbpnts: i32,
    pnts: Vec<f64>,
    nbintv: i32,
    intv: Vec<f64>,
}

impl TheQuadCurvExactHInter {
    /// `IntCurveSurface_TheQuadCurvExactHInter(S, C)` (`.cxx:28-43`), which
    /// calls `IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection`
    /// (`.pxx:51-135`).
    pub fn new(surface: &dyn Surface, curve: &dyn Curve) -> Self {
        let mut r = Self {
            nbpnts: -1,
            pnts: Vec::new(),
            nbintv: -1,
            intv: Vec::new(),
        };

        let quadric = match int_surf_quadric_of(surface) {
            Some(q) => q,
            // `default: { return; }` (`pxx:82-84`) leaves `theNbPnts` at -1.
            None => return r,
        };

        let nb_intervals = curve.nb_intervals(2); // GeomAbs_C1
        let intervals = curve.parameter_intervals(2);

        let mut completed = true;
        for ii in 0..nb_intervals.max(0) as usize {
            let u1 = intervals[ii];
            let u2 = intervals[ii + 1];

            let sample = FunctionSample::new(u1, u2, nb_samples(curve, u1, u2) as i32);
            let mut func = TheQuadCurvFuncOfTheQuadCurvExactHInter::new(quadric.clone(), curve);
            let roots = FunctionAllRoots::new(&mut func, &sample, EPSX, EPSDIST, EPSNUL);

            if roots.is_done() {
                for i in 1..=roots.nb_points() {
                    r.pnts.push(roots.get_point(i));
                }
                for i in 1..=roots.nb_intervals() {
                    let (a, b) = roots.get_interval(i);
                    r.intv.push(a);
                    r.intv.push(b);
                }
            } else {
                // `else { break; }` (`pxx:121-123`): `ii` stays `<= aNbIntervals`
                // so the trailing `if (ii > aNbIntervals)` is false.
                completed = false;
                break;
            }
        }

        if completed {
            r.nbpnts = r.pnts.len() as i32;
            r.nbintv = (r.intv.len() / 2) as i32;
        }
        r
    }

    /// `IsDone()` (`.cxx:45-48`): `nbpnts != -1`.
    pub fn is_done(&self) -> bool {
        self.nbpnts != -1
    }

    /// `NbRoots()` (`.cxx:52-55`).
    pub fn nb_roots(&self) -> i32 {
        self.nbpnts
    }

    /// `NbIntervals()` (`.cxx:57-61`).
    pub fn nb_intervals(&self) -> i32 {
        self.nbintv
    }

    /// `Root(Index)` (`.cxx:65-68`), 1-based as in OCCT.
    pub fn root(&self, index: i32) -> f64 {
        self.pnts[(index - 1) as usize]
    }

    /// `Intervals(Index, a, b)` (`.cxx:70-78`), 1-based as in OCCT.
    pub fn intervals(&self, index: i32) -> (f64, f64) {
        let index2 = index + index - 1;
        (self.intv[(index2 - 1) as usize], self.intv[index2 as usize])
    }
}
