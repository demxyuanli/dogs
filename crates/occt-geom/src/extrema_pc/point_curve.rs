use super::prelude::*;
use super::*;


/// Deduplicate by point proximity and sort ascending by distance.
pub(super) fn dedupe_sort(v: Vec<ExtremaPair>) -> Vec<ExtremaPair> {
    let mut out: Vec<ExtremaPair> = Vec::new();
    for e in v {
        if let Some(o) = out.iter_mut().find(|o| o.p2.distance(&e.p2) < 1e-6) {
            if e.distance < o.distance {
                *o = e;
            }
        } else {
            out.push(e);
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

// ---------------------------------------------------------------------------
// Public dispatch.
// ---------------------------------------------------------------------------

/// `Extrema_ExtPElC` arms selected by `Extrema_GGExtPC`'s curve-type switch
/// (`Extrema_GGExtPC.hxx:390-405`, ported in `general_extrema_pc`): Line / Circle / Ellipse /
/// Hyperbola / Parabola go to `Extrema_ExtPElC`.
///
/// Each entry carries OCCT's `myIsMin` flag: line is always a minimum
/// (`Extrema_ExtPElC.cxx:77`), for a circle `Usol[0]` (near point) is the
/// minimum and `Usol[1] = Usol[0] + PI` the maximum (`cxx:177-188`); for an
/// ellipse the flag compares against `C(Us + 0.1)` (`cxx:277`) and for a
/// hyperbola/parabola against `C(Us + 1)` (`cxx:381`, `:473`).
///
/// Returns `None` when the curve is not an elementary type (`Extrema_GGExtPC`
/// then takes its `default:` arm in `general_extrema_pc`).
pub(super) fn ext_pelc_all(
    c: &dyn Curve,
    p: &GpPnt,
    uinf: f64,
    usup: f64,
) -> Option<Vec<(ExtremaPair, bool)>> {
    if let Some(l) = c.gp_line() {
        // `Extrema_GGExtPC` reaches `Extrema_ExtPElC` with the curve's **own**
        // `gp_Lin` (`theCurve.Line()`, `Extrema_GGExtPC.hxx:390-405`), so the
        // returned parameter lives in the curve's parameter space. Building the
        // analytic line at `d0(uinf)` instead made every parameter off by
        // `-uinf` (measured on `data/occ/OffsetPlaneHoleEdge.step`: the
        // two-`Project` edge range came out `(2, 11)` where the true projection
        // is `(0, 10)`, i.e. shifted by the window lower bound).
        return Some(line_all(&l, p, uinf, usup).into_iter().map(|e| (e, true)).collect());
    }
    if c.is_line() {
        // Line-like curve without its own `gp_Lin`: rebuild it from the curve
        // origin `d0(0)`, which is the line's `Location` under the port's
        // arc-length parameterisation (`GeomLine::d0(u) = Loc + u * Dir`).
        let loc = c.d0(0.0);
        let dir = match GpDir::from_vec(&c.d1(0.0).1) {
            Ok(d) => d,
            Err(_) => return Some(Vec::new()),
        };
        let l = GpLin::from_pnt_dir(loc, dir);
        return Some(line_all(&l, p, uinf, usup).into_iter().map(|e| (e, true)).collect());
    }
    if let Some(gc) = c.gp_circ() {
        return Some(
            circle_all(&gc, p, uinf, usup)
                .into_iter()
                .enumerate()
                .map(|(i, e)| (e, i == 0))
                .collect(),
        );
    }
    if let Some(e) = c.gp_ellipse() {
        // `Extrema_ExtPElC.cxx:270-279`: `myIsMin = sqDist(Us) < |P − C(Us+0.1)|²`.
        return Some(
            ellipse_all(&e, p, uinf, usup)
                .into_iter()
                .map(|pair| {
                    let is_min = pair.distance * pair.distance
                        < p.square_distance(&clib::ellipse_value(&e, pair.u1 + 0.1));
                    (pair, is_min)
                })
                .collect(),
        );
    }
    if let Some(h) = c.gp_hyperbola() {
        // `Extrema_ExtPElC.cxx:377-384`: `myIsMin` at parameter step `+1`.
        return Some(
            hyperbola_all(&h, p, uinf, usup)
                .into_iter()
                .map(|pair| {
                    let is_min = pair.distance * pair.distance
                        < p.square_distance(&clib::hyperbola_value(&h, pair.u1 + 1.0));
                    (pair, is_min)
                })
                .collect(),
        );
    }
    if let Some(pa) = c.gp_parabola() {
        // `Extrema_ExtPElC.cxx:469-476`: `myIsMin` at parameter step `+1`.
        return Some(
            parabola_all(&pa, p, uinf, usup)
                .into_iter()
                .map(|pair| {
                    let is_min = pair.distance * pair.distance
                        < p.square_distance(&clib::parabola_value(&pa, pair.u1 + 1.0));
                    (pair, is_min)
                })
                .collect(),
        );
    }
    None
}

/// All local extrema (minima and maxima) of the point–curve distance,
/// deduplicated and sorted by distance.
///
/// `Extrema_ExtPC` is the `Extrema_GGExtPC` type alias in OCCT 8.0.0
/// (`Extrema_ExtPC.hxx:31-38`): its type switch sends the elementary curves to
/// `Extrema_ExtPElC` and BSpline / Bezier / OtherCurve to the `default:` arm
/// (`hxx:390-502`). Both arms are ported (`poly_roots` analytic arms, `general_extrema_pc` exact
/// bracketing + Newton refinement); the previous grid + 60-step Newton sampler
/// (audit A7) is gone.
///
/// Note: like OCCT, the range ends are **not** extrema by themselves —
/// `Extrema_ExtPElC` keeps only solutions inside `[Uinf, Usup]` (`cxx:180`) and
/// the `default:` arm adds an end only when the query point coincides with it
/// (`postprocess_ends`, `GGExtPC.hxx:474-502`). Callers that need endpoint
/// handling (e.g. `ShapeAnalysis_Curve::Project`, `cxx:161-182`) do it
/// themselves.
pub fn point_curve_extrema_all(c: &dyn Curve, p: &GpPnt) -> Vec<ExtremaPair> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let pairs: Vec<ExtremaPair> = match ext_pelc_all(c, p, a, b) {
        Some(v) => v.into_iter().map(|(e, _)| e).collect(),
        None => extrema_ext_pc_range(c, p, a, b)
            .into_iter()
            .map(|s| ExtremaPair {
                p1: *p,
                p2: s.point,
                distance: s.sq_dist.sqrt(),
                u1: s.u,
                v1: None,
                u2: s.u,
                v2: None,
            })
            .collect(),
    };
    dedupe_sort(pairs)
}

/// Minimum distance from `p` to `c` (with the closest point and parameter).
///
/// `None` means OCCT's `IsDone() == false`: `Extrema_GGExtPC::Perform` leaves the
/// solution list empty and callers must test `IsDone` before using the results
/// (`Extrema_GGExtPC.hxx:531`, `:545-550`). The previous body fabricated a
/// golden-section "refine" result instead (audit A15), which made a failure
/// indistinguishable from a real extremum.
pub fn point_curve_extrema(c: &dyn Curve, p: &GpPnt) -> Option<ExtremaPair> {
    point_curve_extrema_all(c, p).into_iter().next()
}

/// Maximum distance from `p` to `c` (farthest local extremum); `None` as above.
pub fn point_curve_max_extrema(c: &dyn Curve, p: &GpPnt) -> Option<ExtremaPair> {
    point_curve_extrema_all(c, p).into_iter().last()
}
