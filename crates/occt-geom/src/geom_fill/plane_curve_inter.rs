//! `IntCurveSurface_HInter` restricted to a plane surface: the parameter/point
//! pairs where a 3D curve crosses a plane.
//!
//! Source: `IntCurveSurface_HInter::Perform(C, S)` with
//! `S = GeomAdaptor_Surface(Geom_Plane)` (`IntCurveSurface_HInter.cxx:55-120`)
//! drives the exact quadric/curve path
//! `IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection`
//! (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:50-136`) on a
//! `GeomAbs_Plane` quadric, whose scalar function is
//! `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter` (`...cxx:25-56`)
//! and whose root finder is `math_FunctionAllRoots` over a
//! `math_FunctionSample` of the interval.
//!
//! Two GeomFill classes need this walk: `GeomFill_GuideTrihedronPlan::Init`
//! (`GeomFill_GuideTrihedronPlan.cxx:138-175`) and
//! `GeomFill_SectionPlacement::Perform` (`GeomFill_SectionPlacement.cxx:467-489`),
//! which is why it lives here instead of inside either of them.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_all_roots::FunctionAllRoots;
use occt_core::math_function_sample::FunctionSample;
use occt_core::precision::PCONFUSION;

use crate::curve::Curve;

/// `EPSX` (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:29`).
const EPSX: f64 = 1.0e-14;
/// `EPSDIST` (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:30`).
const EPSDIST: f64 = 1.0e-8;
/// `EPSNUL` (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:31`).
const EPSNUL: f64 = 1.0e-8;

/// `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter` (`...cxx:25-56`)
/// for the plane case of `IntSurf_Quadric::Distance`
/// (`IntSurf_Quadric.cxx:173-179`) and `Gradient`.
struct QuadCurvFunc<'a> {
    plane_pnt: &'a GpPnt,
    normal: &'a GpVec,
    curve: &'a dyn Curve,
}

impl MathFunction for QuadCurvFunc<'_> {
    /// `Value` (`...cxx:32-36`): `myQuadric.Distance(Value(Param))`.
    fn value(&mut self, param: f64, f: &mut f64) -> bool {
        let p = self.curve.d0(param);
        *f = plane_distance(self.plane_pnt, self.normal, &p);
        true
    }
}

impl MathFunctionWithDerivative for QuadCurvFunc<'_> {
    /// `Derivative` (`...cxx:38-45`): `T . myQuadric.Gradient(P)`.
    fn derivative(&mut self, param: f64, d: &mut f64) -> bool {
        let (_, t) = self.curve.d1(param);
        *d = t.dot(self.normal);
        true
    }

    /// `Values` (`...cxx:47-56`).
    fn values(&mut self, param: f64, f: &mut f64, d: &mut f64) -> bool {
        let (p, t) = self.curve.d1(param);
        *f = plane_distance(self.plane_pnt, self.normal, &p);
        *d = t.dot(self.normal);
        true
    }
}

/// `IntSurf_Quadric::Distance` for `GeomAbs_Plane`: `prm1*X + prm2*Y + prm3*Z
/// + prm4` (`IntSurf_Quadric.cxx:177-178`), i.e. `normal . (P - myPnt)`.
fn plane_distance(plane_pnt: &GpPnt, normal: &GpVec, p: &GpPnt) -> f64 {
    GpVec::from_pnts(plane_pnt, p).dot(normal)
}

/// The `W` / point pairs of the intersections of `curve` with the plane
/// through `plane_pnt` normal to `normal`, in parameter order.
///
/// `IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection`
/// (`...pxx:50-136`) for `GeomAbs_Plane`: per C1 interval, a
/// `math_FunctionSample` of `NbSamples` points and a `math_FunctionAllRoots`
/// on the quadric function.
///
/// `pub` because `BRepFill_PipeShell::Add` (`BRepFill_PipeShell.cxx:505-529`)
/// runs the same `IntCurveSurface_HInter::Perform(C, GeomAdaptor_Surface)`
/// against the plane of the profile face.
pub fn plane_intersections(
    curve: &dyn Curve,
    plane_pnt: &GpPnt,
    normal: &GpVec,
) -> Vec<(f64, GpPnt)> {
    // `cxx:135` intersects `myTrimG`; an `Adaptor3d_Curve` trim only exposes
    // the parameter window, while `GeomTrimmedCurveBasis` forwards the basis
    // intervals, so the window is clipped here.
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    let nb_intervals = curve.nb_intervals(Shape::C1 as u8);
    if nb_intervals <= 0 {
        return Vec::new();
    }
    let disc = curve.parameter_intervals(Shape::C1 as u8);

    let mut points = Vec::new();
    for i in 0..nb_intervals as usize {
        if i + 1 >= disc.len() {
            break;
        }
        let u1 = disc[i].max(first.min(last));
        let u2 = disc[i + 1].min(first.max(last));
        if u2 - u1 <= PCONFUSION {
            continue;
        }
        let sample = FunctionSample::new(u1, u2, nb_samples(curve, u1, u2) as i32);
        let mut func = QuadCurvFunc {
            plane_pnt,
            normal,
            curve,
        };
        let roots = FunctionAllRoots::new(&mut func, &sample, EPSX, EPSDIST, EPSNUL);
        if !roots.is_done() {
            break;
        }
        for k in 1..=roots.nb_points() {
            let w = roots.point(k);
            points.push((w, curve.d0(w)));
        }
    }
    points
}

/// `IntCurveSurface_TheHCurveTool::NbSamples` (`IntCurveSurface_TheHCurveTool.cxx:31-70`),
/// also ported as `occt_topo::intcurvesurface::curve_nb_samples`, which this
/// crate cannot reach (it would invert the `occt-geom` / `occt-topo` layering).
fn nb_samples(curve: &dyn Curve, u0: f64, u1: f64) -> usize {
    const NBS_OTHER: f64 = 10.0;
    let mut nbs = NBS_OTHER;

    if curve.is_line() {
        nbs = 2.0;
    } else if let Some(poles) = curve.bezier_poles() {
        nbs = 3.0 + poles.len() as f64;
    } else if let Some(degree) = curve.nurbs_degree() {
        let nb_knots = curve.bspline_knots().map_or(0, unique_knot_count) as f64;
        nbs = nb_knots * degree as f64 * (curve.last_parameter() - curve.first_parameter());
        let range = u1 - u0;
        if range.abs() > PCONFUSION {
            nbs /= range;
        }
        if nbs < 2.0 {
            nbs = 2.0;
        }
    }

    if nbs > 50.0 {
        nbs = 50.0;
    }
    nbs as usize
}

/// `NbKnots` of a flat knot vector.
fn unique_knot_count(knots: &[f64]) -> usize {
    let mut n = 0usize;
    let mut prev = f64::NAN;
    for &k in knots {
        if !(prev == k) {
            n += 1;
            prev = k;
        }
    }
    n
}
