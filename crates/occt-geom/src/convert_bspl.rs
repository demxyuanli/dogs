//! BSpline conversion depth. Port of the `GeomConvert` gap-fillers
//! (TKGeomBase): `GeomConvert_BSplineCurveKnotSplitting`,
//! `GeomConvert_BSplineSurfaceKnotSplitting` and
//! `GeomConvert_BSplineSurfaceToBezierSurface`.
//!
//! Curve→Bezier is already covered by `crate::bspline_to_bezier` and is not
//! duplicated here; `GeomConvert_ApproxCurve` (adaptive spline approximation)
//! is covered pragmatically by `curve_approx` (polyline) plus
//! `curve_reparam::resample_bspline` (interpolation) and is deferred.
//!
//! **UNPORTED (audit A8 / task T-44)**: `GeomConvert_CompCurveToBSplineCurve` and
//! the arms of `GeomConvert::CurveToBSplineCurve` that need `Geom_BSplineCurve::Segment`,
//! `Geom_BezierCurve::Segment`, `GeomConvert_ApproxCurve` or the periodic
//! `Geom_BSplineCurve` representation (see the notes on `curve_to_bspline_curve`
//! further down); the port's earlier sampling substitute was removed in batch 66.
//!
//! `GeomConvert::CurveToBSplineCurve` itself is ported for the trimmed
//! line/conic arms and the Bezier / B-spline copy arms (batch 92).

use occt_core::bspl::bezier::boehm_insert;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::knots::{self, hunt, insert_knot, multiplicity};
use occt_core::convert::{
    circle_to_bspline_curve, circle_to_bspline_curve_range, ellipse_to_bspline_curve,
    ellipse_to_bspline_curve_range, hyperbola_to_bspline_curve, parabola_to_bspline_curve,
    ConicToBSplineCurve, ConvertError, ParameterisationType,
};
use occt_core::gp::{
    GpAx2, GpAx22d, GpAx3, GpCirc2d, GpElips2d, GpHypr2d, GpParab2d, GpPnt, GpTrsf,
};

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;

// ---------------------------------------------------------------------------
// Knot splitting (GeomConvert_BSplineCurveKnotSplitting /
// BSplineSurfaceKnotSplitting).
// ---------------------------------------------------------------------------

/// Distinct knot values and their multiplicities from a full knot vector.
fn distinct_knots(knots: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut vals: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &k in knots {
        if let Some(last) = vals.last() {
            if (k - last).abs() < 1e-15 {
                *mults.last_mut().unwrap() += 1;
                continue;
            }
        }
        vals.push(k);
        mults.push(1);
    }
    (vals, mults)
}

/// Split indices for `degree`/`mults` at the requested continuity range.
///
/// Returns 0-based indices into the distinct-knot array (the OCCT
/// `splitIndexes`), starting at the first and ending at the last distinct
/// knot. Port of `GeomConvert_BSplineCurveKnotSplitting` / the surface
/// counterpart.
fn knot_splits(degree: usize, mults: &[usize], continuity_range: usize) -> Vec<usize> {
    let first = 0usize;
    let last = mults.len() - 1;
    if continuity_range == 0 {
        return vec![first, last];
    }
    let mmax = *mults.iter().max().unwrap_or(&0);
    if degree >= mmax + continuity_range {
        return vec![first, last];
    }
    let mut out = vec![first];
    for i in 1..last {
        // Continuity at knot i is (degree − multiplicity); split where it
        // drops below `continuity_range`.
        if degree < mults[i] + continuity_range {
            out.push(i);
        }
    }
    out.push(last);
    out
}

/// Split indices of `c` where its continuity drops below C1 (the OCCT default
/// `ContinuityRange = 1`). Indices are into the curve's distinct-knot array.
pub fn knot_splitting_curve(c: &GeomBSplineCurve) -> Vec<usize> {
    knot_splitting_curve_with_continuity(c, 1)
}

/// `knot_splitting_curve` with an explicit continuity range.
pub fn knot_splitting_curve_with_continuity(c: &GeomBSplineCurve, continuity_range: usize) -> Vec<usize> {
    let (_, mults) = distinct_knots(&c.knots);
    knot_splits(c.degree, &mults, continuity_range)
}

/// Split indices of `s` in u and v (each a `Vec<usize>` into the respective
/// distinct-knot array), at the OCCT default continuity range 1.
pub fn knot_splitting_surface(s: &GeomBSplineSurface) -> (Vec<usize>, Vec<usize>) {
    knot_splitting_surface_with_continuity(s, 1, 1)
}

/// `knot_splitting_surface` with explicit u/v continuity ranges.
pub fn knot_splitting_surface_with_continuity(
    s: &GeomBSplineSurface,
    u_continuity: usize,
    v_continuity: usize,
) -> (Vec<usize>, Vec<usize>) {
    let (_, um) = distinct_knots(&s.knots_u);
    let (_, vm) = distinct_knots(&s.knots_v);
    (
        knot_splits(s.deg_u, &um, u_continuity),
        knot_splits(s.deg_v, &vm, v_continuity),
    )
}

// ---------------------------------------------------------------------------
// BSpline surface → Bezier surface patches
// (GeomConvert_BSplineSurfaceToBezierSurface).
// ---------------------------------------------------------------------------

/// A single Bezier patch of a split B-spline surface, over the parameter
/// rectangle `[u0, u1] × [v0, v1]` with a `(deg_u+1) × (deg_v+1)` pole grid.
#[derive(Debug, Clone)]
pub struct BezierSurfacePatch {
    pub poles: Vec<Vec<GpPnt>>,
    pub weights: Option<Vec<Vec<f64>>>,
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

fn insert_surface_u(s: &GeomBSplineSurface, u: f64, mult: usize) -> GeomBSplineSurface {
    let nu = s.poles.len();
    let nv = s.poles[0].len();
    let deg = s.deg_u;
    let new_nu = nu + mult;
    let mut new_poles = vec![vec![GpPnt::zero(); nv]; new_nu];
    let mut new_weights = s.weights.as_ref().map(|_| vec![vec![0.0; nv]; new_nu]);
    for j in 0..nv {
        let mut col: Vec<GpPnt> = (0..nu).map(|i| s.poles[i][j]).collect();
        let mut wcol = s.weights.as_ref().map(|w| (0..nu).map(|i| w[i][j]).collect::<Vec<_>>());
        let idx = hunt(&s.knots_u, u);
        for _ in 0..mult {
            boehm_insert(&mut col, &s.knots_u, idx, u, deg, wcol.as_mut());
        }
        for i in 0..new_nu {
            new_poles[i][j] = col[i];
            if let (Some(nw), Some(wc)) = (new_weights.as_mut(), wcol.as_ref()) {
                nw[i][j] = wc[i];
            }
        }
    }
    let mut knots_u = s.knots_u.clone();
    for _ in 0..mult {
        knots_u = insert_knot(&knots_u, u, 1);
    }
    GeomBSplineSurface {
        poles: new_poles,
        knots_u,
        knots_v: s.knots_v.clone(),
        deg_u: s.deg_u,
        deg_v: s.deg_v,
        weights: new_weights,
        u_periodic: s.u_periodic,
        v_periodic: s.v_periodic,
    }
}

fn insert_surface_v(s: &GeomBSplineSurface, v: f64, mult: usize) -> GeomBSplineSurface {
    let nu = s.poles.len();
    let nv = s.poles[0].len();
    let deg = s.deg_v;
    let new_nv = nv + mult;
    let mut new_poles = vec![vec![GpPnt::zero(); new_nv]; nu];
    let mut new_weights = s.weights.as_ref().map(|_| vec![vec![0.0; new_nv]; nu]);
    for i in 0..nu {
        let mut row: Vec<GpPnt> = s.poles[i].clone();
        let mut wrow = s.weights.as_ref().map(|w| w[i].clone());
        let idx = hunt(&s.knots_v, v);
        for _ in 0..mult {
            boehm_insert(&mut row, &s.knots_v, idx, v, deg, wrow.as_mut());
        }
        for j in 0..new_nv {
            new_poles[i][j] = row[j];
            if let (Some(nw), Some(wr)) = (new_weights.as_mut(), wrow.as_ref()) {
                nw[i][j] = wr[j];
            }
        }
    }
    let mut knots_v = s.knots_v.clone();
    for _ in 0..mult {
        knots_v = insert_knot(&knots_v, v, 1);
    }
    GeomBSplineSurface {
        poles: new_poles,
        knots_u: s.knots_u.clone(),
        knots_v,
        deg_u: s.deg_u,
        deg_v: s.deg_v,
        weights: new_weights,
        u_periodic: s.u_periodic,
        v_periodic: s.v_periodic,
    }
}

fn distinct_in_range(knots: &[f64], first: f64, last: f64) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if k < first - 1e-15 || k > last + 1e-15 {
            continue;
        }
        if out.last().map_or(true, |&b| (b - k).abs() > 1e-15) {
            out.push(k);
        }
    }
    out
}

/// Convert a clamped B-spline surface into its Bezier patches by knot
/// insertion: every interior knot is raised to multiplicity `degree`, after
/// which each span is a Bezier patch. Mirrors
/// `GeomConvert_BSplineSurfaceToBezierSurface`.
pub fn bspline_surface_to_bezier_surface(s: &GeomBSplineSurface) -> Vec<BezierSurfacePatch> {
    let mut s = s.clone();
    let du = s.deg_u;
    let dv = s.deg_v;
    let (fu, lu) = (s.knots_u[du], s.knots_u[s.knots_u.len() - 1 - du]);
    let (fv, lv) = (s.knots_v[dv], s.knots_v[s.knots_v.len() - 1 - dv]);

    let interior_u: Vec<f64> = s
        .knots_u
        .iter()
        .copied()
        .filter(|&k| k > fu + 1e-15 && k < lu - 1e-15)
        .collect::<Vec<_>>()
        .into_iter()
        .fold(Vec::new(), |mut acc, k| {
            if acc.last().map_or(true, |&l| (l - k).abs() > 1e-15) {
                acc.push(k);
            }
            acc
        });
    for u in interior_u {
        let m = multiplicity(&s.knots_u, u);
        if m < du {
            s = insert_surface_u(&s, u, du - m);
        }
    }
    let interior_v: Vec<f64> = s
        .knots_v
        .iter()
        .copied()
        .filter(|&k| k > fv + 1e-15 && k < lv - 1e-15)
        .collect::<Vec<_>>()
        .into_iter()
        .fold(Vec::new(), |mut acc, k| {
            if acc.last().map_or(true, |&l| (l - k).abs() > 1e-15) {
                acc.push(k);
            }
            acc
        });
    for v in interior_v {
        let m = multiplicity(&s.knots_v, v);
        if m < dv {
            s = insert_surface_v(&s, v, dv - m);
        }
    }

    let fu = s.knots_u[du];
    let lu = s.knots_u[s.knots_u.len() - 1 - du];
    let fv = s.knots_v[dv];
    let lv = s.knots_v[s.knots_v.len() - 1 - dv];
    let bpu = distinct_in_range(&s.knots_u, fu, lu);
    let bpv = distinct_in_range(&s.knots_v, fv, lv);
    if bpu.len() < 2 || bpv.len() < 2 {
        return Vec::new();
    }

    let mut patches = Vec::with_capacity((bpu.len() - 1) * (bpv.len() - 1));
    for iu in 0..bpu.len() - 1 {
        for iv in 0..bpv.len() - 1 {
            let a_u = iu * du;
            let a_v = iv * dv;
            let mut poles = Vec::with_capacity(du + 1);
            let mut weights = s.weights.as_ref().map(|_| Vec::with_capacity(du + 1));
            for i in a_u..=a_u + du {
                poles.push(s.poles[i][a_v..=a_v + dv].to_vec());
                if let (Some(w), Some(sw)) = (weights.as_mut(), s.weights.as_ref()) {
                    w.push(sw[i][a_v..=a_v + dv].to_vec());
                }
            }
            patches.push(BezierSurfacePatch {
                poles,
                weights,
                u0: bpu[iu],
                u1: bpu[iu + 1],
                v0: bpv[iv],
                v1: bpv[iv + 1],
            });
        }
    }
    patches
}

// ---------------------------------------------------------------------------
// CurveToBSplineCurve.
// ---------------------------------------------------------------------------

/// `Geom_BSplineCurve`'s file-static `Rational(Weights)`
/// (`Geom_BSplineCurve.cxx:97-108`): the curve is rational as soon as two
/// consecutive weights differ by more than `gp::Resolution()` (`REAL_SMALL`).
fn weights_are_rational(weights: &[f64]) -> bool {
    for w in weights.windows(2) {
        if (w[0] - w[1]).abs() > occt_core::precision::REAL_SMALL {
            return true;
        }
    }
    false
}

/// `BSplineCurveBuilder` (`GeomConvert.cxx:60-83`): lift the converter's 2D
/// poles to `z = 0`, build the B-spline from its knot/multiplicity tables, then
/// move it from the conic's own frame into the world by
/// `gp_Trsf::SetTransformation(TheConic->Position(), gp::XOY())`.
///
/// OCCT builds the **rational** constructor and lets `CheckRational` decide
/// (`Geom_BSplineCurve.cxx:205-223`); [`weights_are_rational`] reproduces that
/// decision so `IsRational()` matches (a parabola's unit weights end up
/// non-rational, a circular arc's do not).
fn bspline_curve_builder(
    conic_position: &GpAx2,
    convert: &ConicToBSplineCurve,
) -> Result<GeomBSplineCurve, ConvertError> {
    // `CheckCurveData` (`Geom_BSplineCurve.cxx:91-94`):
    // `Poles.Length() == BSplCLib::NbPoles(Degree, Periodic, Mults)`.
    if knots::nb_poles(
        convert.degree() as i32,
        convert.is_periodic(),
        convert.multiplicities(),
    ) as usize
        != convert.poles().len()
    {
        return Err(ConvertError::ConstructionError);
    }
    let poles: Vec<GpPnt> = convert
        .poles()
        .iter()
        .map(|p| GpPnt::new(p.x(), p.y(), 0.0))
        .collect();
    // `Geom_BSplineCurve`'s rational/non-rational constructors call
    // `updateKnots()`, i.e. `BSplCLib::KnotSequence(…, Periodic)`; for a
    // periodic result that is the sequence extended by one period each side.
    let flat = if convert.is_periodic() {
        knots::knot_sequence_periodic(
            convert.knots(),
            convert.multiplicities(),
            convert.degree() as i32,
        )
    } else {
        knot_sequence(
            convert.knots(),
            convert.multiplicities(),
            convert.degree() as i32,
        )
    };
    let weights = convert.weights().to_vec();
    let weights = if weights_are_rational(&weights) { Some(weights) } else { None };
    let mut curve = GeomBSplineCurve {
        poles,
        weights,
        knots: flat,
        degree: convert.degree(),
        periodic: convert.is_periodic(),
    };
    let mut trsf = GpTrsf::identity();
    trsf.set_transformation_from_to(&conic_position.to_ax3(), &GpAx3::standard());
    Curve::transform(&mut curve, &trsf);
    Ok(curve)
}

/// `GeomConvert::CurveToBSplineCurve(C, Parameterisation)` (`GeomConvert.cxx:163-430`);
/// `Parameterisation` defaults to `Convert_TgtThetaOver2`
/// (`GeomConvert.hxx:270-272`).
///
/// **UNPORTED arms** (they return [`ConvertError::Unported`]; each names the
/// OCCT control flow that is still missing in this repository):
/// - a trimmed `Geom_BezierCurve` (`:300-321`) needs `Geom_BezierCurve::Segment`;
/// - a trimmed `Geom_BSplineCurve` (`:322-339`) needs
///   `Geom_BSplineCurve::Segment` (`Geom_BSplineCurve.cxx:527-660`) and
///   `SetNotPeriodic` (`:974-...`);
/// - the `U2 - U1 >= 6` sub-arm of a trimmed circle/ellipse under
///   `Convert_RationalC1` (`:224-242`, `:262-280`) needs
///   `GeomConvert_CompCurveToBSplineCurve`;
/// - a `Geom_OffsetCurve` (`:340-354`, `:436-450`) needs `GeomConvert_ApproxCurve`.
///
/// The non-trimmed `Geom_Circle`/`Geom_Ellipse` arms (`:363-408`) end with
/// `TheCurve->SetPeriodic()` (`:378`, `:383`, `:406`); the periodic
/// `Geom_BSplineCurve` representation they need is ported
/// (`GeomBSplineCurve::set_periodic`, `Geom_BSplineCurve.cxx:777-815`) —
/// board card **R2-21**.
///
/// An unrecognised curve type throws `Standard_DomainError("No such curve")`
/// (`:355-358`, `:451-454`) — returned as [`ConvertError::DomainError`].
pub fn curve_to_bspline_curve(
    c: &dyn Curve,
    parameterisation: ParameterisationType,
) -> Result<GeomBSplineCurve, ConvertError> {
    if c.is_geom_trimmed() {
        // `Curv = Ctrim->BasisCurve()`, `U1 = FirstParameter()`,
        // `U2 = LastParameter()` (`GeomConvert.cxx:171-175`).
        let (basis, _b_first, _b_last) = c.untrimmed_basis().ok_or(ConvertError::DomainError)?;
        let mut u1 = c.first_parameter();
        let mut u2 = c.last_parameter();
        // `cxx:177-189`: for a non-periodic basis the range is clamped, so that
        // `BS->Segment` cannot raise.
        if !basis.is_periodic() {
            if u1 < basis.first_parameter() {
                u1 = basis.first_parameter();
            }
            if u2 > basis.last_parameter() {
                u2 = basis.last_parameter();
            }
        }

        if basis.gp_line().is_some() {
            // `cxx:191-206`: `Poles = {StartPoint, EndPoint}` (the trim's own
            // values), `Knots = {FirstParameter, LastParameter}`, both with
            // multiplicity 2, `Degree = 1`.
            let poles = vec![c.d0(c.first_parameter()), c.d0(c.last_parameter())];
            let knots = vec![c.first_parameter(), c.last_parameter()];
            let flat = knot_sequence(&knots, &[2, 2], 1);
            return GeomBSplineCurve::new(poles, flat, 1)
                .map_err(|_| ConvertError::ConstructionError);
        }

        if let Some(circ) = basis.gp_circ() {
            // `cxx:208-244`: the 2D conic is centred at the origin with the same
            // radius; `BSplineCurveBuilder` moves it back through
            // `TheConic->Position()`.
            let c2d = GpCirc2d::new(GpAx22d::standard(), circ.radius());
            if parameterisation != ParameterisationType::RationalC1 || (u2 - u1) < 6.0 {
                let convert = circle_to_bspline_curve_range(&c2d, u1, u2, parameterisation)?;
                return bspline_curve_builder(&circ.position(), &convert);
            }
            // `cxx:225-242`: `U2 - U1 >= 6` splits the circle in two halves and
            // stitches them with `GeomConvert_CompCurveToBSplineCurve` to avoid
            // the numerical overflow at `U2 - U1 ~ 2*PI`.
            return Err(ConvertError::Unported);
        }

        if let Some(elips) = basis.gp_ellipse() {
            // `cxx:246-282`.
            let e2d = GpElips2d::new(
                GpAx22d::standard(),
                elips.major_radius,
                elips.minor_radius,
            );
            if parameterisation != ParameterisationType::RationalC1 || (u2 - u1) < 6.0 {
                let convert = ellipse_to_bspline_curve_range(&e2d, u1, u2, parameterisation)?;
                return bspline_curve_builder(elips.position(), &convert);
            }
            return Err(ConvertError::Unported);
        }

        if let Some(hypr) = basis.gp_hyperbola() {
            // `cxx:284-290`.
            let h2d = GpHypr2d::new(GpAx22d::standard(), hypr.major_radius, hypr.minor_radius);
            let convert = hyperbola_to_bspline_curve(&h2d, u1, u2)?;
            return bspline_curve_builder(hypr.position(), &convert);
        }

        if let Some(parab) = basis.gp_parabola() {
            // `cxx:292-298`.
            let p2d = GpParab2d::new(GpAx22d::standard(), parab.focal);
            let convert = parabola_to_bspline_curve(&p2d, u1, u2)?;
            return bspline_curve_builder(parab.position(), &convert);
        }

        // `cxx:300-321` (`Geom_BezierCurve`: `CBez->Segment(U1, U2)`) and
        // `cxx:322-339` (`Geom_BSplineCurve`: `AdjustPeriodic` +
        // `SetNotPeriodic` + `Segment`).
        if basis.bezier_poles().is_some() || basis.bspline_poles().is_some() {
            return Err(ConvertError::Unported);
        }
        // `cxx:340-354`: `GeomConvert_ApproxCurve(C, 1e-4, C2, 16, 14)`.
        if basis.offset_curve().is_some() {
            return Err(ConvertError::Unported);
        }
        // `throw Standard_DomainError("No such curve")` (`cxx:355-358`).
        return Err(ConvertError::DomainError);
    }

    // Non-trimmed arm (`cxx:361-455`).
    if let Some(elips) = c.gp_ellipse() {
        // `cxx:363-385`: `Convert_EllipseToBSplineCurve(E2d, Parameterisation)`
        // followed by `TheCurve->SetPeriodic()`.
        let e2d = GpElips2d::new(GpAx22d::standard(), elips.major_radius, elips.minor_radius);
        let convert = ellipse_to_bspline_curve(&e2d, parameterisation)?;
        let mut curve = bspline_curve_builder(elips.position(), &convert)?;
        curve.set_periodic();
        return Ok(curve);
    }

    if let Some(circ) = c.gp_circ() {
        // `cxx:387-408`.
        let c2d = GpCirc2d::new(GpAx22d::standard(), circ.radius());
        let convert = circle_to_bspline_curve(&c2d, parameterisation)?;
        let mut curve = bspline_curve_builder(&circ.position(), &convert)?;
        curve.set_periodic();
        return Ok(curve);
    }

    if let Some(poles) = c.bezier_poles() {
        // `cxx:410-429`: a Bezier is exactly the clamped B-spline with the same
        // poles, knots `{0, 1}` of multiplicity `degree + 1`; a rational Bezier
        // passes `WeightsArray()` (`cxx:421-424`). This port's `GeomBezierCurve`
        // is non-rational (`bezier_curve.rs`), which matches
        // `Geom_BezierCurve::IsRational() == false` for that representation.
        let degree = poles
            .len()
            .checked_sub(1)
            .ok_or(ConvertError::ConstructionError)?;
        let mut knots = vec![0.0f64; degree + 1];
        knots.extend(std::iter::repeat(1.0).take(degree + 1));
        return GeomBSplineCurve::new(poles.to_vec(), knots, degree)
            .map_err(|_| ConvertError::ConstructionError);
    }

    if let (Some(poles), Some(knots), Some(degree)) =
        (c.bspline_poles(), c.bspline_knots(), c.nurbs_degree())
    {
        // `cxx:431-434`: `TheCurve = C->Copy()`. Trait objects cannot be
        // downcast in this port, so the copy is rebuilt from the same
        // `Geom_BSplineCurve` queries (`bspline_poles`/`bspline_weights`/
        // `bspline_knots`/`nurbs_degree`/`is_periodic`) — equal data, not a
        // re-approximation.
        return Ok(GeomBSplineCurve {
            poles: poles.to_vec(),
            weights: c.bspline_weights().map(|w| w.to_vec()),
            knots: knots.to_vec(),
            degree,
            periodic: c.is_periodic(),
        });
    }

    // `cxx:436-450`: `GeomConvert_ApproxCurve` for a `Geom_OffsetCurve`.
    if c.offset_curve().is_some() {
        return Err(ConvertError::Unported);
    }

    // `throw Standard_DomainError("No such curve")` (`cxx:451-454`).
    Err(ConvertError::DomainError)
}

// ---------------------------------------------------------------------------
// CompCurveToBSplineCurve.
// ---------------------------------------------------------------------------

// UNPORTED (audit A8 / task T-44): `GeomConvert_CompCurveToBSplineCurve`
// (`GeomConvert_CompCurveToBSplineCurve.cxx:135-215`) concatenates the **exact**
// B-spline images of its segments - each produced by
// `GeomConvert::CurveToBSplineCurve` (`GeomConvert_CurveToBSplineCurve.cxx` →
// `Convert_LineToBSplineCurve` / `Convert_CircleToBSplineCurve` /
// `Convert_EllipseToBSplineCurve` / `Convert_ParabolaToBSplineCurve` /
// `Convert_HyperbolaToBSplineCurve`, all rational and exact) - and stitches them
// by knot/pole surgery, reparametrising and inserting knots as needed.
//
// This port used to substitute a **sampled polyline fitted as a degree-1
// B-spline** (`n = clamp((b-a)/0.1, 2, 64)` points per segment), which is not
// OCCT's construction. That function had no production caller
// (`git grep comp_curve_to_bspline` reached only its own tests), so it was
// removed in batch 66 rather than left as an invented rule; porting the OCCT
// algorithm above is deferred until a consumer needs it.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface::Surface;
    use occt_core::gp::GpPnt;

    fn degree3_multiknot_curve() -> GeomBSplineCurve {
        // Cubic; interior 0.5 (mult 1) keeps C2, interior 0.75 (mult 3) is C0.
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
            GpPnt::new(3.0, 0.0, 0.0),
            GpPnt::new(4.0, 0.0, 0.0),
            GpPnt::new(5.0, 0.0, 0.0),
            GpPnt::new(6.0, 0.0, 0.0),
            GpPnt::new(7.0, 0.0, 0.0),
        ];
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.75, 0.75, 0.75, 1.0, 1.0, 1.0, 1.0];
        GeomBSplineCurve::new(poles, knots, 3).unwrap()
    }

    #[test]
    fn knot_splitting_curve_finds_c0_knot() {
        let c = degree3_multiknot_curve();
        let splits = knot_splitting_curve(&c);
        // Distinct knots: [0 (m4), 0.5 (m1), 0.75 (m3), 1 (m4)].
        // Only 0.75 (index 2) drops below C1.
        assert_eq!(splits, vec![0, 2, 3], "splits {splits:?}");
    }

    #[test]
    fn knot_splitting_curve_uniform_stays_whole() {
        let c = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0., 0., 0.),
                GpPnt::new(1., 1., 0.),
                GpPnt::new(2., 0., 0.),
                GpPnt::new(3., -1., 0.),
                GpPnt::new(4., 0., 0.),
                GpPnt::new(5., 1., 0.),
            ],
            vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0],
            2,
        )
        .unwrap();
        // All interior knots have multiplicity 1 → degree−mult = C1 → no split.
        let splits = knot_splitting_curve(&c);
        assert_eq!(splits, vec![0, 4], "splits {splits:?}");
    }

    #[test]
    fn knot_splitting_surface_u_and_v() {
        // Degree 2 in u (interior mult 1 → C1), degree 2 in v (interior mult 2 → C0).
        // u: 4 poles (7 knots), v: 5 poles (8 knots).
        let poles = vec![vec![GpPnt::new(0., 0., 0.); 5]; 4];
        let ku = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        let kv = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];
        let s = GeomBSplineSurface::new(poles, ku, kv, 2, 2).unwrap();
        let (us, vs) = knot_splitting_surface(&s);
        // u: distinct [0,0.5,1], mult 1 → no split → [0, 2].
        // v: distinct [0,0.5,1], mult 2 = degree → split at index 1 → [0,1,2].
        assert_eq!(us, vec![0, 2], "u splits {us:?}");
        assert_eq!(vs, vec![0, 1, 2], "v splits {vs:?}");
    }

    fn bilinear_surface() -> GeomBSplineSurface {
        let poles = vec![
            vec![GpPnt::new(0., 0., 0.), GpPnt::new(0., 1., 0.)],
            vec![GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.)],
        ];
        let ku = vec![0.0, 0.0, 1.0, 1.0];
        let kv = vec![0.0, 0.0, 1.0, 1.0];
        GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap()
    }

    #[test]
    fn surface_to_bezier_single_patch_bilinear() {
        let s = bilinear_surface();
        let patches = bspline_surface_to_bezier_surface(&s);
        assert_eq!(patches.len(), 1);
        let p = &patches[0];
        assert_eq!(p.poles.len(), 2);
        assert_eq!(p.poles[0].len(), 2);
        assert!((p.u0 - 0.0).abs() < 1e-12 && (p.u1 - 1.0).abs() < 1e-12);
        assert!((p.v0 - 0.0).abs() < 1e-12 && (p.v1 - 1.0).abs() < 1e-12);
        // Patch corners equal the surface corners.
        assert!(p.poles[0][0].distance(&s.d0(0.0, 0.0)) < 1e-9);
        assert!(p.poles[0][1].distance(&s.d0(0.0, 1.0)) < 1e-9);
        assert!(p.poles[1][0].distance(&s.d0(1.0, 0.0)) < 1e-9);
        assert!(p.poles[1][1].distance(&s.d0(1.0, 1.0)) < 1e-9);
    }

    #[test]
    fn surface_to_bezier_patches_match_at_shared_corners() {
        // Degree-2 surface in u with an interior knot → two u-patches; the
        // shared corner must coincide with the original surface.
        let poles = vec![
            vec![GpPnt::new(0., 0., 0.), GpPnt::new(0., 1., 0.)],
            vec![GpPnt::new(0.5, 0.5, 0.), GpPnt::new(0.5, 0.5, 0.)],
            vec![GpPnt::new(1.5, 0.5, 0.), GpPnt::new(1.5, 0.5, 0.)],
            vec![GpPnt::new(2., 0., 0.), GpPnt::new(2., 1., 0.)],
        ];
        let ku = vec![0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0];
        let kv = vec![0.0, 0.0, 1.0, 1.0];
        let s = GeomBSplineSurface::new(poles, ku, kv, 2, 1).unwrap();

        let patches = bspline_surface_to_bezier_surface(&s);
        // u: distinct [0,1,2] → 2 u-patches; v: distinct [0,1] → 1 v-patch.
        assert_eq!(patches.len(), 2, "patches {patches:?}");
        let (p0, p1) = (&patches[0], &patches[1]);
        assert!((p0.u1 - p1.u0).abs() < 1e-12, "u1 {} vs u0 {}", p0.u1, p1.u0);
        // Shared corner between the patches equals the surface at the breakpoint.
        let shared_u = p0.u1;
        for v in [0.0, 1.0] {
            let got = p0.poles[p0.poles.len() - 1][if v < 0.5 { 0 } else { 1 }];
            let expected = s.d0(shared_u, v);
            assert!(got.distance(&expected) < 1e-9, "corner v={v}: {got:?} vs {expected:?}");
            let got2 = p1.poles[0][if v < 0.5 { 0 } else { 1 }];
            assert!(got2.distance(&expected) < 1e-9);
        }
        // Every patch's four corners sit on the original surface.
        for p in &patches {
            let (nu, nv) = (p.poles.len() - 1, p.poles[0].len() - 1);
            for &(i, j) in &[(0, 0), (0, nv), (nu, 0), (nu, nv)] {
                let u = if i == 0 { p.u0 } else { p.u1 };
                let v = if j == 0 { p.v0 } else { p.v1 };
                assert!(p.poles[i][j].distance(&s.d0(u, v)) < 1e-9, "corner ({i},{j})");
            }
        }
    }
}
