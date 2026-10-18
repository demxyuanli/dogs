//! `Adaptor3d_CurveOnSurface::NbIntervals` / `Intervals`
//! (`Adaptor3d_CurveOnSurface.cxx:1010-1129`) and `Adaptor3d_InterFunc`
//! (`Adaptor3d_InterFunc.cxx:24-75`).

use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::PCONFUSION;
use occt_geom::Surface;
use occt_geom2d::Curve2d;

/// `Adaptor3d_InterFunc` (`Adaptor3d_InterFunc.cxx:24-75`):
/// `Fix == 1` solves `C.X() - FixVal = 0`, `Fix == 2` solves `C.Y() - FixVal = 0`.
struct InterFunc<'a> {
    curve: &'a dyn Curve2d,
    fix_val: f64,
    fix: i32,
}

impl MathFunction for InterFunc<'_> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let c = self.curve.d0(x);
        *f = if self.fix == 1 {
            c.x() - self.fix_val
        } else {
            c.y() - self.fix_val
        };
        true
    }
}

impl MathFunctionWithDerivative for InterFunc<'_> {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values(x, &mut f, d)
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (c, dc) = self.curve.d1(x);
        if self.fix == 1 {
            *f = c.x() - self.fix_val;
            *d = dc.x();
        } else {
            *f = c.y() - self.fix_val;
            *d = dc.y();
        }
        true
    }
}

/// `AddIntervals` (`Adaptor3d_CurveOnSurface.cxx:1010-1041`).
fn add_intervals(parameters: &mut Vec<f64>, roots: &FunctionRoots, tol: f64) {
    if !roots.is_done() || roots.is_all_null() {
        return;
    }
    let nsol = roots.nb_solutions();
    for i in 1..=nsol {
        let param = roots.value(i);
        if parameters.is_empty() || param - parameters[0] < tol {
            continue;
        }
        let mut inserted = false;
        for j in 1..parameters.len() {
            let delta = parameters[j] - param;
            if delta > tol {
                parameters.insert(j, param);
                inserted = true;
                break;
            } else if delta >= -tol {
                inserted = true;
                break;
            }
        }
        let _ = inserted;
    }
}

/// Sorted break points for `Adaptor3d_CurveOnSurface::Intervals`
/// (`cxx:1045-1129`): start from the 2D curve's intervals, then add roots where
/// the pcurve crosses each interior U/V surface discontinuity.
pub fn curve_on_surface_intervals(
    pcurve: &dyn Curve2d,
    surface: &dyn Surface,
    continuity: u8,
) -> Vec<f64> {
    let nu = surface.nb_u_intervals(continuity);
    let nv = surface.nb_v_intervals(continuity);
    let tdeb = pcurve.first_parameter();
    let tfin = pcurve.last_parameter();
    let mut intervals = pcurve.parameter_intervals(continuity);
    if intervals.is_empty() {
        intervals.push(tdeb);
        intervals.push(tfin);
    }

    let tol = PCONFUSION / 10.0;
    let nb_sample = 20;

    if nu > 1 {
        let tab_u = surface.u_intervals(continuity);
        // Interior breaks: `iu = 2..=nu` on a 1-based array of length `nu+1`
        // (`cxx:1084-1090`) — skip the first and last knots.
        for iu in 1..tab_u.len().saturating_sub(1) {
            let u = tab_u[iu];
            let mut func = InterFunc {
                curve: pcurve,
                fix_val: u,
                fix: 1,
            };
            let resol = FunctionRoots::new(&mut func, tdeb, tfin, nb_sample, tol, tol, tol, 0.0);
            add_intervals(&mut intervals, &resol, tol);
        }
    }
    if nv > 1 {
        let tab_v = surface.v_intervals(continuity);
        for iv in 1..tab_v.len().saturating_sub(1) {
            let v = tab_v[iv];
            let mut func = InterFunc {
                curve: pcurve,
                fix_val: v,
                fix: 2,
            };
            let resol = FunctionRoots::new(&mut func, tdeb, tfin, nb_sample, tol, tol, tol, 0.0);
            add_intervals(&mut intervals, &resol, tol);
        }
    }

    // `cxx:1107-1110`: a single-value sequence gets a duplicate end.
    if intervals.len() == 1 {
        let v = intervals[0];
        intervals.push(v);
    }
    intervals
}
