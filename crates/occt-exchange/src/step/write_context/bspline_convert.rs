use super::super::prelude::*;

// ---------------------------------------------------------------------------
// B-spline surface access
// ---------------------------------------------------------------------------

/// The surface's own B-spline data: exact poles / knots / weights.
///
/// `GeomToStep_MakeBoundedSurface.cxx:41-70`: the
/// `IsKind(STANDARD_TYPE(Geom_BSplineSurface))` arm reads the surface itself,
/// and the `IsKind(STANDARD_TYPE(Geom_BezierSurface))` arm first converts it
/// with `GeomConvert::SurfaceToBSplineSurface`. A
/// `Geom_RectangularTrimmedSurface` is peeled to its basis, which goes through
/// this same chain (`cxx:82-88` -> `GeomToStep_MakeRectangularTrimmedSurface.cxx:46-53`).
/// `Surface::osculating_bspline` covers the first two arms exactly (the
/// B-spline clone, the Bezier conversion), so the pole accessors below only
/// serve a basis reached through a wrapper.
fn exact_bspline_of(s: &dyn Surface) -> Option<GeomBSplineSurface> {
    if let Some(bs) = s.osculating_bspline() {
        return Some(bs);
    }
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return exact_bspline_of(basis.as_ref());
    }
    let poles = s.bspline_surface_poles()?.to_vec();
    let knots_u = s.bspline_surface_uknots()?.to_vec();
    let knots_v = s.bspline_surface_vknots()?.to_vec();
    let deg_u = usize::try_from(s.u_degree()).ok()?;
    let deg_v = usize::try_from(s.v_degree()).ok()?;
    match s.bspline_surface_weights() {
        Some(w) => {
            GeomBSplineSurface::rational(poles, w.to_vec(), knots_u, knots_v, deg_u, deg_v).ok()
        }
        None => GeomBSplineSurface::new(poles, knots_u, knots_v, deg_u, deg_v).ok(),
    }
}

/// The B-spline surface `GeomToStep_MakeBoundedSurface` writes, unperiodized.
///
/// `GeomToStep_MakeBoundedSurface.cxx:46-52`: a periodic B-spline is copied and
/// `SetUNotPeriodic` / `SetVNotPeriodic` is applied before the poles and knots
/// are read, so the written entity carries the open knot vector
/// (`BSplSLib::Unperiodize`, [`occt_core::bspl::unperiodize_direction`], the
/// same call the IGES writer makes for entity 128).
///
/// Sampling the surface on a 6x6 grid and re-fitting it (the previous body) has
/// no counterpart in OCCT: it cannot reproduce the source poles, and OCCT's
/// reader is given the surface's own control net.
///
/// UNPORTED: the trim of a `Geom_RectangularTrimmedSurface` is not wrapped in a
/// `RECTANGULAR_TRIMMED_SURFACE` entity (`cxx:82-88`); only the basis surface is
/// written, as `iges/emit_curve.rs` does for the same class. The trim itself is carried by
/// the `ADVANCED_FACE` bounds.
pub(super) fn bounded_bspline_surface(s: &dyn Surface) -> Option<GeomBSplineSurface> {
    let mut bs = exact_bspline_of(s)?;
    if bs.u_periodic() {
        let (knots, map) = occt_core::bspl::unperiodize_direction(
            bs.deg_u() as i32,
            bs.knots_u(),
            bs.nb_poles_u(),
        );
        let poles: Vec<Vec<GpPnt>> = map.iter().map(|k| bs.poles[*k].clone()).collect();
        bs.poles = poles;
        if let Some(w) = bs.weights.take() {
            bs.weights = Some(map.iter().map(|k| w[*k].clone()).collect());
        }
        bs.set_u_periodic_flag(false);
        bs.set_uknots(knots);
    }
    if bs.v_periodic() {
        let (knots, map) = occt_core::bspl::unperiodize_direction(
            bs.deg_v() as i32,
            bs.knots_v(),
            bs.nb_poles_v(),
        );
        for row in bs.poles.iter_mut() {
            *row = map.iter().map(|k| row[*k].clone()).collect();
        }
        if let Some(w) = bs.weights.as_mut() {
            for row in w.iter_mut() {
                *row = map.iter().map(|k| row[*k]).collect();
            }
        }
        bs.set_v_periodic_flag(false);
        bs.set_vknots(knots);
    }
    let nu = bs.poles.len();
    let nv = bs.poles.first().map_or(0, |r| r.len());
    if bs.knots_u().len() != nu + bs.deg_u() + 1 || bs.knots_v().len() != nv + bs.deg_v() + 1 {
        return None;
    }
    // `poles` and the flat knot sequences were reassigned above, so the
    // `Resolution` bound latched on `bs` no longer describes it
    // (`Geom_BSplineSurface_1.cxx:2197-2221`).
    bs.res_latch.invalidate();
    Some(bs)
}

/// `GeomConvert::CurveToBSplineCurve` (`GeomConvert.cxx:163-430`) through the
/// ported `occt-geom` entry point; a Bezier is converted this way by
/// `GeomToStep_MakeBoundedCurve.cxx:64-75` before it is written.
pub(super) fn curve_to_bspline(c: &dyn Curve) -> Option<GeomBSplineCurve> {
    occt_geom::convert_bspl::curve_to_bspline_curve(
        c,
        occt_core::convert::ParameterisationType::TgtThetaOver2,
    )
    .ok()
}

/// `BS = B->Copy(); BS->SetNotPeriodic()` (`GeomToStep_MakeBoundedCurve.cxx:44-51`).
///
/// `Geom_BSplineCurve::Copy` keeps the source periodicity, and the port's
/// constructors (`GeomBSplineCurve::new` / `rational`) never set the flag, so a
/// copy rebuilt from the curve's poles and knots has to have it restored before
/// `set_not_periodic()` (`Geom_BSplineCurve.cxx:974-1019`) can open the knot
/// vector. `periodic` is `C->IsPeriodic()` of the curve as received (`cxx:45`).
fn open_periodic_bspline(bs: &mut GeomBSplineCurve, periodic: bool) {
    if periodic {
        bs.set_periodic_flag(true);
        bs.set_not_periodic();
    }
}

/// The B-spline `GeomToStep_MakeBoundedCurve` writes for a bounded curve.
///
/// `GeomToStep_MakeBoundedCurve.cxx:42-62`: the `Geom_BSplineCurve` arm reads
/// the curve itself (through the non-periodic copy when `C->IsPeriodic()`), and
/// the `Geom_BezierCurve` arm is converted first by
/// `GeomConvert::CurveToBSplineCurve` (that conversion is `curve_to_bspline`
/// above, used by the Bezier caller).
pub(super) fn bounded_bspline_of_curve(c: &dyn Curve) -> Option<GeomBSplineCurve> {
    let poles = c.bspline_poles()?.to_vec();
    let knots = c.bspline_knots()?.to_vec();
    let degree = c.nurbs_degree()?;
    let weights = c.bspline_weights().map(|w| w.to_vec());
    let mut bs = match weights {
        Some(w) if w.len() == poles.len() => {
            GeomBSplineCurve::rational(poles, w, knots, degree).ok()?
        }
        _ => GeomBSplineCurve::new(poles, knots, degree).ok()?,
    };
    open_periodic_bspline(&mut bs, c.is_periodic());
    Some(bs)
}
