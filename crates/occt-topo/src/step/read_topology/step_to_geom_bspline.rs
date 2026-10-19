//! `StepToGeom::MakeBSplineCurveCommon` (`StepToGeom.cxx:750-927`) and its 2D arm `MakeBSplineCurve2d` (`cxx:952-963`): multiplicity normalisation, periodic promotion, and the shared `BsplineDescriptor`.
use super::*;

/// The pole-type-independent part of `StepToGeom::MakeBSplineCurveCommon`
/// (`StepToGeom.cxx:776-894`), shared by the 3D and 2D arms (OCCT implements
/// both with one template, `:938-944` vs `:952-963`).
pub(super) struct BsplineDescriptor {
    uknots: Vec<f64>,
    umults: Vec<i32>,
    /// Pole window `[start, end)` left after the multiplicity trim.
    pole_range: std::ops::Range<usize>,
    should_be_periodic: bool,
}

pub(super) fn bspline_descriptor(
    degree: usize,
    n_poles: usize,
    mults: &[usize],
    knots: &[f64],
) -> Result<BsplineDescriptor, String> {
    // 1. unique knots + summed multiplicities (`cxx:784-821`).
    let mut uknots: Vec<f64> = Vec::new();
    let mut umults: Vec<i32> = Vec::new();
    for (i, &k) in knots.iter().enumerate() {
        let m = mults.get(i).copied().unwrap_or(1) as i32;
        if let Some(&last) = uknots.last() {
            if k - last <= occt_core::precision::epsilon(last.abs()) {
                *umults.last_mut().unwrap() += m;
                continue;
            }
        }
        uknots.push(k);
        umults.push(m);
    }
    if uknots.len() <= 1 {
        // `if (NbUniqueKnots <= 1) return nullptr;` (`cxx:795-798`).
        return Err("B_SPLINE_CURVE: at most one unique knot".into());
    }

    // 2. clamp multiplicities above `degree + 1` (`cxx:823-845`).
    let deg1 = degree as i32 + 1;
    let n_unique = umults.len();
    let mut first_diff = 0usize;
    let mut last_diff = 0usize;
    for i in 0..n_unique {
        if umults[i] > deg1 {
            if i == 0 {
                first_diff = (umults[i] - deg1) as usize;
            }
            if i == n_unique - 1 {
                last_diff = (umults[i] - deg1) as usize;
            }
            umults[i] = deg1;
        }
    }

    // 3. drop the extra poles (`cxx:847-871`).
    let diff = first_diff + last_diff;
    let n_unique_poles = n_poles
        .checked_sub(diff)
        .filter(|n| *n > 0)
        .ok_or_else(|| "B_SPLINE_CURVE: no poles left after multiplicity trim".to_string())?;

    // 4. periodic-looking descriptor (`cxx:873-894`).
    let summary: i32 = umults.iter().sum();
    let should_be_periodic = summary != (n_unique_poles as i32 + degree as i32 + 1)
        && umults[0] == umults[n_unique - 1]
        && (summary - umults[0]) == n_unique_poles as i32;

    Ok(BsplineDescriptor {
        uknots,
        umults,
        pole_range: first_diff..first_diff + n_unique_poles,
        should_be_periodic,
    })
}

/// `StepToGeom::MakeBSplineCurveCommon` (`StepToGeom.cxx:750-927`), the part
/// shared by the `B_SPLINE_CURVE_WITH_KNOTS` and Bezier / uniform /
/// quasi-uniform arms once the knot and multiplicity lists are known:
///
/// 1. duplicate knots are merged within `Epsilon(|lastKnot|)` and their
///    multiplicities summed (`cxx:784-821`);
/// 2. multiplicities above `degree + 1` are clamped to `degree + 1` and the
///    corresponding leading / trailing poles (and weights) are dropped
///    (`cxx:823-871`, `cxx:902-906`);
/// 3. a descriptor that "looks periodic" (`cxx:873-894`: the total multiplicity
///    is *not* `NbPoles + degree + 1`, the end multiplicities are equal and
///    `Σmults - mults(1) == NbPoles`) is built in the periodic representation,
///    i.e. `Geom_BSplineCurve(Poles, Knots, Mults, Degree, true)` whose flat
///    sequence is the periodic `BSplCLib::KnotSequence`;
/// 4. a closed curve (`ClosedCurve() && Degree() > 1 && IsClosed()`,
///    `cxx:920-926`) is forced periodic with `SetPeriodic()`.
///
/// `IsClosed()` is `StartPoint().SquareDistance(EndPoint()) <=
/// Precision::Computational()` (`Geom_BSplineCurve.cxx:146-149`).
pub(super) fn make_bspline_curve_with_knots(
    degree: usize,
    poles: Vec<GpPnt>,
    weights: Option<Vec<f64>>,
    mults: &[usize],
    knots: &[f64],
    closed: bool,
) -> Result<GeomBSplineCurve, String> {
    let d = bspline_descriptor(degree, poles.len(), mults, knots)?;
    // 5. drop the extra poles / weights (`cxx:847-871`).
    let weights = match weights {
        Some(w) if w.len() >= d.pole_range.end => Some(w[d.pole_range.clone()].to_vec()),
        Some(_) => return Err("B_SPLINE_CURVE: weight count mismatch".into()),
        None => None,
    };
    let poles = poles[d.pole_range.clone()].to_vec();
    let (uknots, umults) = (&d.uknots, &d.umults);

    let mut curve = if d.should_be_periodic {
        // `new TBSplineCurve(Poles, [Weights,] UniqueKnots, UniqueMults, Degree,
        // true)`: `CheckCurveData` (`Geom_BSplineCurve.cxx:91-94`) requires
        // `NbPoles(Degree, true, Mults) == Poles.Length()`.
        if occt_core::bspl::knots::nb_poles(degree as i32, true, umults) as usize != poles.len() {
            return Err("B_SPLINE_CURVE: periodic pole/degree mismatch".into());
        }
        let flat = occt_core::bspl::knots::knot_sequence_periodic(uknots, umults, degree as i32);
        GeomBSplineCurve::from_parts(poles, weights, flat, degree, true)
    } else {
        let flat = occt_core::bspl::banded_interp::knot_sequence(uknots, umults, degree as i32);
        match weights {
            Some(w) => {
                GeomBSplineCurve::rational(poles, w, flat, degree).map_err(|e| e.to_string())?
            }
            None => GeomBSplineCurve::new(poles, flat, degree).map_err(|e| e.to_string())?,
        }
    };

    // 6. force periodicity on closed curves (`cxx:920-926`).
    if closed && degree > 1 {
        let sp = curve.d0(curve.first_parameter());
        let ep = curve.d0(curve.last_parameter());
        if sp.square_distance(&ep) <= occt_core::precision::COMPUTATIONAL {
            curve.set_periodic();
        }
    }
    Ok(curve)
}

/// The 2D arm of `StepToGeom::MakeBSplineCurveCommon`
/// (`StepToGeom::MakeBSplineCurve2d`, `StepToGeom.cxx:952-963`): identical
/// knot/multiplicity surgery, built into `Geom2d_BSplineCurve` (non-rational:
/// this port's 2D B-spline stores `xs`/`ys` without weights).
pub(super) fn make_bspline_curve_2d_with_knots(
    degree: usize,
    xs: Vec<f64>,
    ys: Vec<f64>,
    mults: &[usize],
    knots: &[f64],
    closed: bool,
) -> Result<Geom2dBSplineCurve, String> {
    let d = bspline_descriptor(degree, xs.len(), mults, knots)?;
    let xs = xs[d.pole_range.clone()].to_vec();
    let ys = ys[d.pole_range.clone()].to_vec();
    let (uknots, umults) = (&d.uknots, &d.umults);

    let mut curve = if d.should_be_periodic {
        if occt_core::bspl::knots::nb_poles(degree as i32, true, umults) as usize != xs.len() {
            return Err("B_SPLINE_CURVE_2D: periodic pole/degree mismatch".into());
        }
        let flat = occt_core::bspl::knots::knot_sequence_periodic(uknots, umults, degree as i32);
        Geom2dBSplineCurve { xs, ys, weights: None, knots: flat, degree, periodic: true }
    } else {
        let flat = occt_core::bspl::banded_interp::knot_sequence(uknots, umults, degree as i32);
        Geom2dBSplineCurve::new(xs, ys, flat, degree).map_err(|e| e.to_string())?
    };

    if closed && degree > 1 {
        let sp = curve.d0(curve.first_parameter());
        let ep = curve.d0(curve.last_parameter());
        if sp.square_distance(&ep) <= occt_core::precision::COMPUTATIONAL {
            curve.set_periodic();
        }
    }
    Ok(curve)
}

/// The STEP `closed_curve` / `closed` flag of a curve record.
///
/// The resolver sees two layouts: a plain instance `(name, degree, poles,
/// curve_form, closed, …)` with `closed` at index 4, and a merged rational
/// complex `(name, degree, poles, weights, curve_form, closed, …)` with
/// `closed` at index 5 (`transfer.rs::merge_complex_body`). The weights slot
/// distinguishes them: the merged body puts `SELF` or a `(…)` list there, a
/// plain body a `.FORM.` marker.
pub(super) fn curve_record_closed(rec: &Record) -> bool {
    let merged = rec
        .args
        .get(3)
        .map(|s| {
            let t = s.trim();
            t == "SELF" || t.starts_with('(')
        })
        .unwrap_or(false);
    let idx = if merged { 5 } else { 4 };
    parse_logical(rec.args.get(idx).map(|s| s.as_str()), false)
}
