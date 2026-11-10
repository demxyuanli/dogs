//! `interpolatePCurve` (`cxx:2152-2213`) and the non-periodic `Geom2dAPI_Interpolate` helpers.

use super::super::prelude::*;

// ---------------------------------------------------------------------------
// `interpolatePCurve` (`ShapeConstruct_ProjectCurveOnSurface.cxx:2152-2213`) +
// `Geom2dAPI_Interpolate` non-periodic (`Geom2dAPI_Interpolate.cxx:615-810`)
// ---------------------------------------------------------------------------

/// `ShapeConstruct_ProjectCurveOnSurface::checkPoints2d` (`cxx:2373-2460`):
/// drop every point that coincides with its last kept predecessor
/// (`SquareDistance < gp::Resolution()`, `cxx:2396`) and shrink the tolerance
/// to `0.9 * min distance` (`cxx:2419-2422`). Returns the kept points / params
/// and the updated tolerance.
fn check_points2d(points: &[GpPnt2d], params: &[f64], preci: f64) -> (Vec<GpPnt2d>, Vec<f64>, f64) {
    let n = points.len();
    if n == 0 || params.len() != n {
        return (Vec::new(), Vec::new(), preci);
    }
    // `cxx:2382-2390`: `tmpParam(i)` flags the entries that survive.
    let mut keep = vec![true; n];
    let mut nb_pnt_dropped = 0usize;
    let mut last_valid = 0usize;
    let mut dist_min2 = f64::MAX;
    let mut prev = points[0];
    for i in 1..n {
        let cur_dist2 = prev.square_distance(&points[i]);
        if cur_dist2 < f64::MIN_POSITIVE {
            // `cxx:2398-2408`: drop this point, or the previous kept point when
            // this is the last one.
            nb_pnt_dropped += 1;
            if i + 1 == n {
                keep[last_valid] = false;
            } else {
                keep[i] = false;
            }
        } else {
            if cur_dist2 < dist_min2 {
                dist_min2 = cur_dist2;
            }
            last_valid = i;
            prev = points[i];
        }
    }
    let preci = if dist_min2 < f64::MAX {
        // `cxx:2419-2422`
        0.9 * dist_min2.sqrt()
    } else {
        preci
    };
    if nb_pnt_dropped == 0 {
        return (points.to_vec(), params.to_vec(), preci);
    }
    let mut pts: Vec<GpPnt2d> = points.to_vec();
    let new_last = n - nb_pnt_dropped;
    if new_last < 2 {
        // `cxx:2430-2438`: a minimal-length pcurve, the last point shifted by
        // `(preci, preci)`.
        keep[0] = true;
        keep[n - 1] = true;
        let mut p = pts[n - 1];
        p.set_x(p.x() + preci);
        p.set_y(p.y() + preci);
        pts[n - 1] = p;
    }
    let mut out_p = Vec::new();
    let mut out_t = Vec::new();
    for i in 0..n {
        if keep[i] {
            out_p.push(pts[i]);
            out_t.push(params[i]);
        }
    }
    (out_p, out_t, preci)
}

/// `Geom2dAPI_Interpolate.cxx:156-214` `BuildTangents`, with every
/// `TangentFlag` false as the two-array constructor leaves them
/// (`Geom2dAPI_Interpolate.cxx:299-308`): the end tangents are the first
/// derivatives of the Lagrange polynomial through the first / last
/// `degree + 1` points, evaluated at the end parameter (`PLib::EvalLagrange`,
/// `cxx:184-190` / `cxx:202-208`). `None` marks the
/// `Standard_ConstructionError` thrown for fewer than three points
/// (`cxx:170-175`).
fn build_tangents(points: &[GpPnt2d], params: &[f64]) -> Option<(GpVec2d, GpVec2d)> {
    let n = points.len();
    if n < 3 || params.len() != n {
        return None;
    }
    // `cxx:169-178`: degree 2 for three points, 3 otherwise.
    let degree = if n == 3 { 2 } else { 3 };
    let mut vals = vec![[0.0f64; 2]; degree + 1];
    let mut prms = vec![0.0f64; degree + 1];
    for i in 0..=degree {
        vals[i] = [points[i].x(), points[i].y()];
        prms[i] = params[i];
    }
    let first = occt_core::bspl::plib_eval::eval_lagrange(params[0], 1, degree, &vals, &prms).ok()?;
    // `cxx:196-214`: the last `degree + 1` points, evaluated at `params(n)`.
    let base = n - 1 - degree;
    for i in 0..=degree {
        vals[i] = [points[base + i].x(), points[base + i].y()];
        prms[i] = params[base + i];
    }
    let last =
        occt_core::bspl::plib_eval::eval_lagrange(params[n - 1], 1, degree, &vals, &prms).ok()?;
    Some((GpVec2d::new(first[1][0], first[1][1]), GpVec2d::new(last[1][0], last[1][1])))
}

/// `CheckPoints` (`Geom2dAPI_Interpolate.cxx:33-44`): the `Geom2dAPI_Interpolate`
/// constructor throws `Standard_ConstructionError` unless every consecutive
/// pair is at least `tol` apart.
fn check_points(points: &[GpPnt2d], tol: f64) -> bool {
    let tol2 = tol * tol;
    points.windows(2).all(|w| w[0].square_distance(&w[1]) >= tol2)
}

/// `CheckParameters` (`Geom2dAPI_Interpolate.cxx:70-80`): consecutive
/// parameters must grow by at least `RealSmall()`.
fn check_parameters(params: &[f64]) -> bool {
    params.windows(2).all(|w| w[1] - w[0] >= f64::MIN_POSITIVE)
}

/// `Geom2dAPI_Interpolate::PerformNonPeriodic`
/// (`Geom2dAPI_Interpolate.cxx:615-810`) for `myTangentRequest == false`: a
/// single-span non-periodic degree-3 (degree 2 for three points, degree 1 for
/// two) C2 interpolation through the points, solved by `BSplCLib::Interpolate`
/// with contact orders. `None` marks OCCT's `inversion_problem`
/// (`cxx:802`) and the throw paths of the calls it makes.
fn interpolate_non_periodic(points: &[GpPnt2d], params: &[f64]) -> Option<Geom2dBSplineCurve> {
    let n = points.len();
    if n < 2 || params.len() != n {
        return None;
    }
    // `cxx:623-648`
    let (degree, num_poles) = if n == 2 {
        (1usize, n)
    } else if n == 3 {
        (2usize, 3usize)
    } else {
        (3usize, n + 2)
    };
    // `cxx:650-668`: flat knots and contact orders. `Geom2d_BSplineCurve` builds
    // the same flat knots from `knots` + `mults` (`cxx:680`, `:698`, `:806`).
    let mut flatknots = vec![0.0f64; num_poles + degree + 1];
    let mut contact = vec![0i32; num_poles];

    for ii in 1..=degree + 1 {
        flatknots[ii - 1] = params[0];
        flatknots[ii + num_poles - 1] = params[n - 1];
    }

    let mut flat_poles = vec![0.0f64; num_poles * 2];
    match degree {
        1 => {
            // `cxx:675-683`
            for ii in 0..num_poles {
                flat_poles[2 * ii] = points[ii].x();
                flat_poles[2 * ii + 1] = points[ii].y();
            }
        }
        2 => {
            // `cxx:684-701`: `knots = (params(1), params(3))`, `mults = (3, 3)`.
            for ii in 0..num_poles {
                flat_poles[2 * ii] = points[ii].x();
                flat_poles[2 * ii + 1] = points[ii].y();
            }
            occt_core::bspl::banded_interp::interpolate_contact(
                degree,
                &flatknots,
                params,
                &contact,
                &mut flat_poles,
                2,
            )
            .ok()?;
        }
        _ => {
            // `cxx:703-712`: the end tangents.
            let (tan_first, tan_last) = build_tangents(points, params)?;
            // `cxx:717-726`
            contact[1] = 1;
            let mut parameters = vec![0.0f64; num_poles];
            parameters[0] = params[0];
            parameters[1] = params[0];
            flat_poles[0] = points[0].x();
            flat_poles[1] = points[0].y();
            flat_poles[2] = tan_first.x();
            flat_poles[3] = tan_first.y();
            // `cxx:764-770`
            for ii in 0..n {
                parameters[ii + 1] = params[ii];
            }
            for ii in 1..n - 1 {
                flat_poles[2 * (ii + 1)] = points[ii].x();
                flat_poles[2 * (ii + 1) + 1] = points[ii].y();
            }
            for ii in 0..n {
                flatknots[degree + ii] = params[ii];
            }
            // `cxx:785-795`
            let last = num_poles - 1;
            flat_poles[2 * (last - 1)] = tan_last.x();
            flat_poles[2 * (last - 1) + 1] = tan_last.y();
            contact[last - 1] = 1;
            parameters[last] = params[n - 1];
            parameters[last - 1] = params[n - 1];
            flat_poles[2 * last] = points[n - 1].x();
            flat_poles[2 * last + 1] = points[n - 1].y();
            // `cxx:797-806`
            occt_core::bspl::banded_interp::interpolate_contact(
                degree,
                &flatknots,
                &parameters,
                &contact,
                &mut flat_poles,
                2,
            )
            .ok()?;
        }
    }

    let xs: Vec<f64> = (0..num_poles).map(|i| flat_poles[2 * i]).collect();
    let ys: Vec<f64> = (0..num_poles).map(|i| flat_poles[2 * i + 1]).collect();
    Geom2dBSplineCurve::new(xs, ys, flatknots, degree).ok()
}

/// `ShapeConstruct_ProjectCurveOnSurface::interpolatePCurve`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2152-2213`): the `Perform`
/// fallback (`cxx:771`) reached for every pcurve that `approxPCurve` did not
/// settle analytically. `preci` is `myPreci` (`cxx:2158`).
pub(crate) fn interpolate_pcurve(
    points: &[GpPnt2d],
    params: &[f64],
    preci: f64,
) -> Option<Arc<dyn Curve2d>> {
    if points.is_empty() || points.len() != params.len() {
        return None;
    }
    // `cxx:2158`: `theTolerance2d = myPreci / (100 * theNbPnt)`.
    let tol2d = preci / (100.0 * points.len() as f64);
    // `cxx:2177-2195`: drop coincident points and shrink the tolerance.
    let (pts, prm, preci_val) = check_points2d(points, params, tol2d);
    if pts.len() < 2 || prm.len() != pts.len() {
        return None;
    }
    // `cxx:2197-2202`: `Geom2dAPI_Interpolate(aPnts2d, aParams, false, aPreci)`
    // then `Perform()`. The constructor's `CheckPoints` / `CheckParameters`
    // throws (`cxx:294-353`) surface as the `catch` at `cxx:2204`.
    if !check_points(&pts, preci_val) || !check_parameters(&prm) {
        return None;
    }
    let curve = interpolate_non_periodic(&pts, &prm)?;
    Some(Arc::new(curve))
}
