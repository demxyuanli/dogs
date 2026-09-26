//! BSpline conversion depth. Port of the `GeomConvert` gap-fillers
//! (TKGeomBase): `GeomConvert_BSplineCurveKnotSplitting`,
//! `GeomConvert_BSplineSurfaceKnotSplitting` and
//! `GeomConvert_BSplineSurfaceToBezierSurface`.
//!
//! Curve→Bezier is already covered by `crate::bspline_to_bezier` and is not
//! duplicated here; `GeomConvert_ApproxCurve` (adaptive spline approximation)
//! is ported in `crate::convert_approx_curve` on top of
//! `crate::adv_approx::ApproxAFunction3d` (the `AdvApprox_ApproxAFunction`
//! port) and `AdvApprox_PrefAndRec` cutting.
//!
//! `GeomConvert_CompCurveToBSplineCurve` is ported (constructor + `Add`,
//! `GeomConvert_CompCurveToBSplineCurve.cxx:32-273`); it needs
//! `Geom_BSplineCurve::IncreaseDegree` and `Geom_BSplineCurve::RemoveKnot`
//! (both ported on `bspline_curve.rs`).
//!
//! Both `Geom_OffsetCurve` arms of `GeomConvert::CurveToBSplineCurve`
//! (`GeomConvert.cxx:340-354`, `:436-450`) call `GeomConvert_ApproxCurve(C,
//! 1e-4, C2, 16, 14)`; they are wired to `crate::convert_approx_curve`.
//!
//! `GeomConvert::CurveToBSplineCurve` itself is ported for the trimmed
//! line/conic arms (including the `RationalC1` + `U2-U1>=6` split-stitch arm)
//! and the Bezier / B-spline copy arms (batch 92).

use crate::bezier_curve::GeomBezierCurve;
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
use crate::convert_approx_curve::GeomConvertApproxCurve;
use crate::curve::Curve;
use occt_core::kernel::geomabs::Shape;

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

/// Dynamic type test of OCCT's `IsKind(STANDARD_TYPE(Geom_BSplineCurve))`
/// down_cast. The port's `Curve` trait is not `Any`, so the B-spline query is
/// the discriminator; a `Geom_TrimmedCurve` delegates that query to its basis
/// (`trimmed.rs:127`) and is excluded here exactly as the OCCT down_cast
/// excludes it.
fn is_bspline(c: &dyn Curve) -> bool {
    !c.is_geom_trimmed() && c.bspline_poles().is_some()
}

/// `Geom_BSplineCurve::Copy()` (`Geom_BSplineCurve.cxx:112-115`): a verbatim
/// copy of the `Geom_BSplineCurve` fields, recovered from the `Curve`
/// queries (`bspline_poles`/`bspline_weights`/`bspline_knots`/
/// `nurbs_degree`/`is_periodic`). Equal data, not a re-approximation.
fn bspline_copy(c: &dyn Curve) -> Result<GeomBSplineCurve, ConvertError> {
    let poles = c.bspline_poles().ok_or(ConvertError::DomainError)?;
    let knots = c.bspline_knots().ok_or(ConvertError::DomainError)?;
    let degree = c.nurbs_degree().ok_or(ConvertError::DomainError)?;
    Ok(GeomBSplineCurve {
        poles: poles.to_vec(),
        weights: c.bspline_weights().map(|w| w.to_vec()),
        knots: knots.to_vec(),
        degree,
        periodic: c.is_periodic(),
    })
}

/// `GeomConvert::CurveToBSplineCurve(C, Parameterisation)` (`GeomConvert.cxx:163-430`);
/// `Parameterisation` defaults to `Convert_TgtThetaOver2`
/// (`GeomConvert.hxx:270-272`).
///
/// The two `Geom_OffsetCurve` arms (`:340-354`, `:436-450`) approximate the
/// (possibly trimmed) curve through `GeomConvert_ApproxCurve(C, 1e-4, C2, 16,
/// 14)` and raise `Standard_ConstructionError` when there is no result.
///
/// The other arms are ported: the trimmed line/conic arms (including the
/// `RationalC1` + `U2 - U1 >= 6` split-and-stitch, `:224-242`, `:262-280`),
/// the trimmed Bezier arm (`:300-321`, its rational sub-branch is noted
/// in place) and the trimmed/non-trimmed B-spline copy + `Segment` arms
/// (`:322-339`, `:431-434`).
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
            let u_med = (u1 + u2) * 0.5;
            let convert1 = circle_to_bspline_curve_range(&c2d, u1, u_med, parameterisation)?;
            let curve1 = bspline_curve_builder(&circ.position(), &convert1)?;
            let convert2 = circle_to_bspline_curve_range(&c2d, u_med, u2, parameterisation)?;
            let curve2 = bspline_curve_builder(&circ.position(), &convert2)?;
            // `cxx:237-241`: `Add(TheCurve2, Precision::PConfusion(), true)`
            // (After = true, WithRatio = true, MinM = 0).
            let mut cctbspl = CompCurveToBSplineCurve::from_bspline(curve1, parameterisation);
            cctbspl.add(&curve2, occt_core::precision::PCONFUSION, true, true, 0)?;
            return cctbspl
                .into_curve()
                .ok_or(ConvertError::ConstructionError);
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
            // `cxx:262-280`: the same split-and-stitch as the circle arm.
            let u_med = (u1 + u2) * 0.5;
            let convert1 = ellipse_to_bspline_curve_range(&e2d, u1, u_med, parameterisation)?;
            let curve1 = bspline_curve_builder(elips.position(), &convert1)?;
            let convert2 = ellipse_to_bspline_curve_range(&e2d, u_med, u2, parameterisation)?;
            let curve2 = bspline_curve_builder(elips.position(), &convert2)?;
            // `cxx:275-279`: `Add(TheCurve2, Precision::PConfusion(), true)`.
            let mut cctbspl = CompCurveToBSplineCurve::from_bspline(curve1, parameterisation);
            cctbspl.add(&curve2, occt_core::precision::PCONFUSION, true, true, 0)?;
            return cctbspl
                .into_curve()
                .ok_or(ConvertError::ConstructionError);
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

        // `cxx:300-321` (`Geom_BezierCurve`: `CBez->Segment(U1, U2)`) needs
        // `Geom_BezierCurve::Segment` (`PLib::Trimming` + `PLib::CoefficientsPoles`)
        // — ported (`bezier_curve.rs::segment`, OCCT `Geom_BezierCurve.cxx:388-425`).
        // `cxx:300-321`: the copy is trimmed and re-wrapped as a B-spline over
        // `[0, 1]` with multiplicities `{Degree+1, Degree+1}`. UNPORTED: the
        // rational branch (`CBez->IsRational()` → `WeightsArray`,
        // `GeomConvert.cxx:313-321`) — the port's `GeomBezierCurve` carries no
        // weights.
        if let Some(poles) = basis.bezier_poles() {
            let mut bez = GeomBezierCurve::new(poles.to_vec())
                .map_err(|_| ConvertError::ConstructionError)?;
            bez.segment(u1, u2);
            let degree = bez.degree();
            let m = (degree + 1) as i32;
            let flat = knot_sequence(&[0.0, 1.0], &[m, m], degree as i32);
            return GeomBSplineCurve::new(bez.poles.clone(), flat, degree)
                .map_err(|_| ConvertError::ConstructionError);
        }
        // `cxx:322-339` (`Geom_BSplineCurve`): the basis is copied, its range is
        // folded into the period (`ElCLib::AdjustPeriodic`), a full-period trim
        // drops the periodic representation (`SetNotPeriodic`) and the copy is
        // then cut with `Segment(U1, U2)`.
        if basis.bspline_poles().is_some() {
            let mut bs = GeomBSplineCurve {
                poles: basis.bspline_poles().ok_or(ConvertError::DomainError)?.to_vec(),
                weights: basis.bspline_weights().map(|w| w.to_vec()),
                knots: basis.bspline_knots().ok_or(ConvertError::DomainError)?.to_vec(),
                degree: basis.nurbs_degree().ok_or(ConvertError::DomainError)?,
                periodic: basis.is_periodic(),
            };
            if bs.is_periodic() {
                let (uf, ul) = (bs.first_parameter(), bs.last_parameter());
                occt_core::elib::clib2d::adjust_periodic(
                    uf,
                    ul,
                    occt_core::precision::CONFUSION,
                    &mut u1,
                    &mut u2,
                );
                if (u1 - uf).abs() <= occt_core::precision::CONFUSION
                    && (u2 - ul).abs() <= occt_core::precision::CONFUSION
                {
                    bs.set_not_periodic();
                }
            }
            // `theTolerance` defaults to `Precision::PConfusion()`
            // (`Geom_BSplineCurve.hxx:322-324`); the throw is
            // `Standard_DomainError`.
            bs.segment(u1, u2, occt_core::precision::PCONFUSION)
                .map_err(|_| ConvertError::DomainError)?;
            return Ok(bs);
        }
        // `cxx:340-354`: `GeomConvert_ApproxCurve(C, 1e-4, C2, 16, 14)` on the
        // original (possibly trimmed) curve `C`; a missing result raises
        // `Standard_ConstructionError` (`:352`).
        if basis.offset_curve().is_some() {
            let appr = GeomConvertApproxCurve::new(c, 1.0e-4, Shape::C2, 16, 14);
            return match appr.curve() {
                Some(curve) => Ok(curve.clone()),
                None => Err(ConvertError::ConstructionError),
            };
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

    if is_bspline(c) {
        // `cxx:431-434`: `TheCurve = C->Copy()`.
        return bspline_copy(c);
    }

    // `cxx:436-450`: `GeomConvert_ApproxCurve(C, 1e-4, C2, 16, 14)`; a missing
    // result raises `Standard_ConstructionError` (`:448`).
    if c.offset_curve().is_some() {
        let appr = GeomConvertApproxCurve::new(c, 1.0e-4, Shape::C2, 16, 14);
        return match appr.curve() {
            Some(curve) => Ok(curve.clone()),
            None => Err(ConvertError::ConstructionError),
        };
    }

    // `throw Standard_DomainError("No such curve")` (`cxx:451-454`).
    Err(ConvertError::DomainError)
}

// ---------------------------------------------------------------------------
// CompCurveToBSplineCurve.
// ---------------------------------------------------------------------------

/// `GeomConvert_CompCurveToBSplineCurve`
/// (`GeomConvert_CompCurveToBSplineCurve.hxx:30-78`): converts and
/// concatenates several curves into one B-spline
/// (`GeomConvert_CompCurveToBSplineCurve.cxx:32-273`).
pub struct CompCurveToBSplineCurve {
    my_curve: Option<GeomBSplineCurve>,
    my_tol: f64,
    my_type: ParameterisationType,
}

impl CompCurveToBSplineCurve {
    /// `GeomConvert_CompCurveToBSplineCurve(Parameterisation)`
    /// (`cxx:32-37`): `myTol = Precision::Confusion()`.
    pub fn new(parameterisation: ParameterisationType) -> Self {
        Self {
            my_curve: None,
            my_tol: occt_core::precision::CONFUSION,
            my_type: parameterisation,
        }
    }

    /// `GeomConvert_CompCurveToBSplineCurve(BasisCurve, Parameterisation)`
    /// (`cxx:41-56`): a `Geom_BSplineCurve` basis is copied, anything else is
    /// converted with `GeomConvert::CurveToBSplineCurve(BasisCurve, myType)`.
    pub fn from_curve(
        basis: &dyn Curve,
        parameterisation: ParameterisationType,
    ) -> Result<Self, ConvertError> {
        let mut s = Self::new(parameterisation);
        s.my_curve = Some(if is_bspline(basis) {
            bspline_copy(basis)?
        } else {
            curve_to_bspline_curve(basis, parameterisation)?
        });
        Ok(s)
    }

    /// The constructor's `down_cast` branch (`cxx:47-51`) for a caller that
    /// already owns a `GeomBSplineCurve`; equivalent to `Copy()`.
    pub fn from_bspline(curve: GeomBSplineCurve, parameterisation: ParameterisationType) -> Self {
        Self {
            my_curve: Some(curve),
            my_tol: occt_core::precision::CONFUSION,
            my_type: parameterisation,
        }
    }

    /// `BSplineCurve()` (`cxx:264-267`).
    pub fn bspline_curve(&self) -> Option<&GeomBSplineCurve> {
        self.my_curve.as_ref()
    }

    /// `Clear()` (`cxx:271-274`).
    pub fn clear(&mut self) {
        self.my_curve = None;
    }

    /// `BSplineCurve()` by value.
    pub fn into_curve(self) -> Option<GeomBSplineCurve> {
        self.my_curve
    }

    /// `Add(NewCurve, Tolerance, After, WithRatio, MinM)` (`cxx:60-131`).
    ///
    /// `After`/`WithRatio`/`MinM` reproduce the OCCT defaults explicitly
    /// (`GeomConvert_CompCurveToBSplineCurve.hxx:56-60`: `false` / `true` /
    /// `0`).
    pub fn add(
        &mut self,
        new_curve: &dyn Curve,
        tolerance: f64,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> Result<bool, ConvertError> {
        // Conversion (`cxx:66-75`).
        let mut bs = if is_bspline(new_curve) {
            bspline_copy(new_curve)?
        } else {
            curve_to_bspline_curve(new_curve, self.my_type)?
        };
        if self.my_curve.is_none() {
            self.my_curve = Some(bs);
            return Ok(true);
        }
        self.my_tol = tolerance;

        // Use actual curve endpoints instead of poles for the G0 continuity
        // check (`cxx:85-94`).
        let (a_curve_start, a_curve_end) = {
            let cur = self.my_curve.as_ref().unwrap();
            (cur.d0(cur.first_parameter()), cur.d0(cur.last_parameter()))
        };
        let a_bs_start = bs.d0(bs.first_parameter());
        let a_bs_end = bs.d0(bs.last_parameter());

        let mut avant = a_curve_start.distance(&a_bs_start) < self.my_tol
            || a_curve_start.distance(&a_bs_end) < self.my_tol;
        let mut apres = a_curve_end.distance(&a_bs_start) < self.my_tol
            || a_curve_end.distance(&a_bs_end) < self.my_tol;

        // Will myCurve be (or become) closed? Resolve the ambiguity
        // (`cxx:96-107`).
        if avant && apres {
            if after {
                avant = false;
            } else {
                apres = false;
            }
        }

        if apres {
            // Append after? (`cxx:109-118`).
            if a_curve_end.distance(&a_bs_end) < self.my_tol {
                bs.reverse();
            }
            let mut first = self.my_curve.take().unwrap();
            self.add_pair(&mut first, &mut bs, true, with_ratio, min_m)?;
            Ok(true)
        } else if avant {
            // Prepend before? (`cxx:119-128`).
            if a_curve_start.distance(&a_bs_start) < self.my_tol {
                bs.reverse();
            }
            let mut first = bs;
            let mut second = self.my_curve.take().unwrap();
            self.add_pair(&mut first, &mut second, false, with_ratio, min_m)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// The private `Add(FirstCurve, SecondCurve, After, WithRatio, MinM)`
    /// (`cxx:135-260`): harmonize the degrees, reparameterize onto the common
    /// knot, concatenate the poles/weights and optionally lower the common
    /// knot's multiplicity down to `MinM`. Sets `myCurve`.
    fn add_pair(
        &mut self,
        first: &mut GeomBSplineCurve,
        second: &mut GeomBSplineCurve,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> Result<(), ConvertError> {
        // Harmonize the degrees (`cxx:141-150`).
        let deg = first.degree().max(second.degree());
        if first.degree() < deg {
            first
                .increase_degree(deg)
                .map_err(|_| ConvertError::ConstructionError)?;
        }
        if second.degree() < deg {
            second
                .increase_degree(deg)
                .map_err(|_| ConvertError::ConstructionError)?;
        }

        // Reparameterization ratio (C1 if possible) (`cxx:163-177`).
        let mut ratio = 1.0f64;
        if with_ratio {
            let l1 = first.eval_dn(first.last_parameter(), 1).magnitude();
            let l2 = second.eval_dn(second.first_parameter(), 1).magnitude();
            if l1 > occt_core::precision::CONFUSION && l2 > occt_core::precision::CONFUSION {
                ratio = l1 / l2;
            }
            if ratio < occt_core::precision::CONFUSION
                || ratio > 1.0 / occt_core::precision::CONFUSION
            {
                ratio = 1.0;
            }
        }

        let (uf, umf) = first.distinct_knots_and_mults();
        let (us, ums) = second.distinct_knots_and_mults();
        let nb_p1 = first.nb_poles();
        let nb_p2 = second.nb_poles();
        let nb_k1 = uf.len();
        let nb_k2 = us.len();
        if nb_k1 == 0 || nb_k2 == 0 || nb_p1 == 0 || nb_p2 == 0 {
            return Err(ConvertError::ConstructionError);
        }

        let (ratio1, delta1, ratio2, delta2);
        if after {
            // Do not move the first curve (`cxx:179-186`).
            ratio1 = 1.0;
            delta1 = 0.0;
            ratio2 = 1.0 / ratio;
            delta2 = ratio2 * us[0] - uf[nb_k1 - 1];
        } else {
            // Do not move the second curve (`cxx:187-194`).
            ratio1 = ratio;
            delta1 = ratio1 * uf[nb_k1 - 1] - us[0];
            ratio2 = 1.0;
            delta2 = 0.0;
        }

        // The knots (`cxx:196-229`).
        let mut noeuds = vec![0.0f64; nb_k1 + nb_k2 - 1];
        let mut mults_out = vec![0i32; nb_k1 + nb_k2 - 1];
        for ii in 1..=nb_k1 {
            noeuds[ii - 1] = ratio1 * uf[ii - 1] - delta1;
            if ii > 1 {
                let mut eps = occt_core::precision::epsilon(noeuds[ii - 2].abs());
                if eps < 5.0e-10 {
                    eps = 5.0e-10;
                }
                if noeuds[ii - 1] - noeuds[ii - 2] <= eps {
                    noeuds[ii - 1] += eps;
                }
            }
            mults_out[ii - 1] = umf[ii - 1];
        }
        mults_out[nb_k1 - 1] = first.degree() as i32;
        for ii in 2..=nb_k2 {
            let jj = nb_k1 + ii - 1;
            noeuds[jj - 1] = ratio2 * us[ii - 1] - delta2;
            let mut eps = occt_core::precision::epsilon(noeuds[jj - 2].abs());
            if eps < 5.0e-10 {
                eps = 5.0e-10;
            }
            if noeuds[jj - 1] - noeuds[jj - 2] <= eps {
                noeuds[jj - 1] += eps;
            }
            mults_out[jj - 1] = ums[ii - 1];
        }

        // The poles and weights (`cxx:231-247`).
        let mut r = first.weight(nb_p1 as i32);
        r /= second.weight(1);
        let mut poles_out: Vec<GpPnt> = Vec::with_capacity(nb_p1 + nb_p2 - 1);
        let mut weights_out: Vec<f64> = Vec::with_capacity(nb_p1 + nb_p2 - 1);
        for ii in 1..nb_p1 {
            poles_out.push(*first.pole(ii - 1));
            weights_out.push(first.weight(ii as i32));
        }
        for ii in 1..=nb_p2 {
            poles_out.push(*second.pole(ii - 1));
            weights_out.push(r * second.weight(ii as i32));
        }

        // Create the BSpline (`cxx:249-250`): the rational constructor with
        // `CheckRational` decides `myRational` (`Geom_BSplineCurve.cxx:205-223`).
        let flat = knot_sequence(&noeuds, &mults_out, deg as i32);
        let mut curve = if weights_are_rational(&weights_out) {
            GeomBSplineCurve::rational(poles_out, weights_out, flat, deg)
                .map_err(|_| ConvertError::ConstructionError)?
        } else {
            GeomBSplineCurve::new(poles_out, flat, deg)
                .map_err(|_| ConvertError::ConstructionError)?
        };

        // Optionally reduce multiplicity down to MinM (`cxx:252-259`).
        let mut ok = true;
        let mut m = mults_out[nb_k1 - 1];
        while m > min_m && ok {
            m -= 1;
            ok = curve
                .remove_knot(nb_k1 as i32, m, self.my_tol)
                .map_err(|_| ConvertError::ConstructionError)?;
        }

        self.my_curve = Some(curve);
        Ok(())
    }
}

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
