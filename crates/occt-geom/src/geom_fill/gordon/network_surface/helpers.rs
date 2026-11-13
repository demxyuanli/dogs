//! Anonymous-namespace helpers of `GeomFill_NetworkSurface.cxx` (OCCT 8.0.1).

use occt_core::bspl::banded_interp::{
    interpolate_contact, knot_sequence, schoenberg_points,
};
use occt_core::bspl::knots::reparametrize;
use occt_core::bspl::surface_d0::homogeneous_d0;
use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};

use crate::bspline_curve::{GeomBSplineCurve, MAX_DEGREE};
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;
use crate::geom_fill::profiler::Profiler;

/// Candidate knot with multiplicity collected before tolerance compaction
/// (`GeomFill_NetworkSurface.cxx:245-249`).
#[derive(Clone, Copy)]
pub(super) struct KnotEntry {
    pub parameter: f64,
    pub multiplicity: i32,
}

/// Interpolation basis used to skin one curve family
/// (`GeomFill_NetworkSurface.cxx:252-263`).
pub(super) struct SkinningBasis {
    pub degree: i32,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub flat_knots: Vec<f64>,
}

impl SkinningBasis {
    /// `SkinningBasis::Init` (`cxx:259-262`).
    pub fn init(parameters: &[f64]) -> Self {
        let (degree, knots, mults, flat_knots) = make_interpolation_basis(parameters);
        Self {
            degree,
            knots,
            mults,
            flat_knots,
        }
    }
}

/// `checkParameters` (`cxx:35-49`).
pub(super) fn check_parameters(parameters: &[f64]) -> bool {
    if parameters.len() < 2 {
        return false;
    }
    for i in 1..parameters.len() {
        if parameters[i] <= parameters[i - 1] + PCONFUSION {
            return false;
        }
    }
    true
}

/// `checkWeights` (`cxx:51-65`).
pub(super) fn check_weights(weights: &[Vec<f64>]) -> bool {
    if weights.is_empty() || weights.iter().any(|row| row.is_empty()) {
        return false;
    }
    for row in weights {
        for &w in row {
            if w <= CONFUSION {
                return false;
            }
        }
    }
    true
}

/// `reparametrizeSurface` (`cxx:74-82`).
pub(super) fn reparametrize_surface(surface: &mut GeomBSplineSurface) {
    let mut u_knots = surface.uknots().to_vec();
    let mut v_knots = surface.vknots().to_vec();
    reparametrize(0.0, 1.0, &mut u_knots);
    reparametrize(0.0, 1.0, &mut v_knots);
    let flat_u = knot_sequence(&u_knots, surface.umults(), surface.deg_u() as i32);
    let flat_v = knot_sequence(&v_knots, surface.vmults(), surface.deg_v() as i32);
    surface.set_uknots(flat_u);
    surface.set_vknots(flat_v);
}

/// `hasRationalCurve` (`cxx:84-94`).
pub(super) fn has_rational_curve(curves: &[GeomBSplineCurve]) -> bool {
    curves.iter().any(|c| c.is_rational())
}

/// `hasPeriodicCurve` (`cxx:96-106`).
pub(super) fn has_periodic_curve(curves: &[GeomBSplineCurve]) -> bool {
    curves.iter().any(|c| c.is_periodic())
}

/// `isCompatibleCurveFamily` (`cxx:108-137`).
pub(super) fn is_compatible_curve_family(curves: &[GeomBSplineCurve]) -> bool {
    let Some(reference) = curves.first() else {
        return false;
    };
    let is_rational = reference.is_rational();
    for curve in curves.iter().skip(1) {
        if curve.is_rational() != is_rational
            || curve.degree() != reference.degree()
            || curve.nb_poles() != reference.nb_poles()
            || curve.distinct_knots().len() != reference.distinct_knots().len()
        {
            return false;
        }
        let n_knots = curve.distinct_knots().len();
        for knot_idx in 0..n_knots {
            let knot = (knot_idx + 1) as i32;
            if curve.multiplicity(knot) != reference.multiplicity(knot)
                || (curve.distinct_knots()[knot_idx] - reference.distinct_knots()[knot_idx]).abs()
                    > PCONFUSION
            {
                return false;
            }
        }
    }
    true
}

/// `prepareCurveFamily` (`cxx:139-168`).
pub(super) fn prepare_curve_family(curves: &mut Vec<GeomBSplineCurve>) -> bool {
    if has_rational_curve(curves) {
        return is_compatible_curve_family(curves);
    }

    let mut profiler = Profiler::default();
    for curve in curves.iter() {
        if profiler.add_curve(curve).is_err() {
            return false;
        }
    }
    if profiler.perform(PCONFUSION).is_err() {
        return false;
    }
    let prepared = profiler.curves();
    if prepared.len() != curves.len() || prepared.iter().any(|c| c.nb_poles() == 0) {
        return false;
    }
    *curves = prepared.to_vec();
    true
}

/// `makeInterpolationBasis` (`cxx:170-207`).
pub(super) fn make_interpolation_basis(
    parameters: &[f64],
) -> (i32, Vec<f64>, Vec<i32>, Vec<f64>) {
    let nb_params = parameters.len();
    let degree = 3i32.min(nb_params as i32 - 1);
    let nb_internal_knots = nb_params.saturating_sub(degree as usize + 1);
    let mut knots = vec![0.0; nb_internal_knots + 2];
    let mut mults = vec![0i32; nb_internal_knots + 2];

    knots[0] = parameters[0];
    mults[0] = degree + 1;
    for knot_idx in 0..nb_internal_knots {
        let mut sum = 0.0;
        for param_offset in 1..=degree as usize {
            sum += parameters[knot_idx + param_offset];
        }
        knots[knot_idx + 1] = sum / degree as f64;
        mults[knot_idx + 1] = 1;
    }
    let last_knot = knots.len() - 1;
    let last_mult = mults.len() - 1;
    knots[last_knot] = parameters[parameters.len() - 1];
    mults[last_mult] = degree + 1;

    let mut flat_knots = Vec::with_capacity(degree as usize + nb_params + 1);
    for (knot_idx, &knot) in knots.iter().enumerate() {
        for _ in 0..mults[knot_idx] {
            flat_knots.push(knot);
        }
    }
    (degree, knots, mults, flat_knots)
}

/// `poleWeight` (`cxx:227-230`).
pub(super) fn pole_weight(curve: &GeomBSplineCurve, pole_idx_1based: i32) -> f64 {
    if curve.is_rational() {
        curve.weight(pole_idx_1based)
    } else {
        1.0
    }
}

/// `interpolatePoles` (`cxx:209-225`): `BSplCLib::Interpolate` with poles and weights.
pub(super) fn interpolate_poles(
    degree: i32,
    flat_knots: &[f64],
    parameters: &[f64],
    poles: &mut [GpPnt],
    weights: &mut [f64],
) -> bool {
    let n = poles.len();
    if n != parameters.len() || n != weights.len() || n == 0 {
        return false;
    }
    let mut data = vec![0.0; n * 4];
    for i in 0..n {
        data[i * 4] = poles[i].x();
        data[i * 4 + 1] = poles[i].y();
        data[i * 4 + 2] = poles[i].z();
        data[i * 4 + 3] = weights[i];
    }
    let contact = vec![0i32; n];
    if interpolate_contact(
        degree as usize,
        flat_knots,
        parameters,
        &contact,
        &mut data,
        4,
    )
    .is_err()
    {
        return false;
    }
    for i in 0..n {
        poles[i] = GpPnt::new(data[i * 4], data[i * 4 + 1], data[i * 4 + 2]);
        weights[i] = data[i * 4 + 3];
    }
    true
}

/// `sameKnotRange` (`cxx:232-242`).
pub(super) fn same_knot_range(
    knots1: &[f64],
    knots2: &[f64],
) -> Option<(f64, f64)> {
    let start = *knots1.first()?;
    let end = *knots1.last()?;
    let start2 = *knots2.first()?;
    let end2 = *knots2.last()?;
    if (start - start2).abs() <= PCONFUSION
        && (end - end2).abs() <= PCONFUSION
        && start < end - PCONFUSION
    {
        Some((start, end))
    } else {
        None
    }
}

/// `normalizedKnot` (`cxx:492-503`).
pub(super) fn normalized_knot(start: f64, end: f64, knot: f64) -> f64 {
    if (knot - start).abs() <= PCONFUSION {
        start
    } else if (knot - end).abs() <= PCONFUSION {
        end
    } else {
        knot
    }
}

/// `collectKnots` (`cxx:506-523`).
pub(super) fn collect_knots(
    degree: i32,
    start: f64,
    end: f64,
    knots: &[f64],
    mults: &[i32],
    entries: &mut Vec<KnotEntry>,
) {
    for (knot_idx, &knot_raw) in knots.iter().enumerate() {
        let knot = normalized_knot(start, end, knot_raw);
        let is_boundary = knot == start || knot == end;
        let max_mult = if is_boundary { degree + 1 } else { degree };
        let mult = mults.get(knot_idx).copied().unwrap_or(1).min(max_mult);
        entries.push(KnotEntry {
            parameter: knot,
            multiplicity: mult,
        });
    }
}

/// `compactKnots` (`cxx:526-558`).
pub(super) fn compact_knots(entries: &mut [KnotEntry]) -> (Vec<f64>, Vec<i32>) {
    if entries.is_empty() {
        return (Vec::new(), Vec::new());
    }
    entries.sort_by(|a, b| {
        a.parameter
            .partial_cmp(&b.parameter)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut knots = Vec::new();
    let mut mults = Vec::new();
    let mut current = entries[0];
    for entry in entries.iter().skip(1) {
        if (entry.parameter - current.parameter).abs() <= PCONFUSION {
            current.multiplicity = current.multiplicity.max(entry.multiplicity);
            continue;
        }
        knots.push(current.parameter);
        mults.push(current.multiplicity);
        current = *entry;
    }
    knots.push(current.parameter);
    mults.push(current.multiplicity);
    (knots, mults)
}

/// `uniteKnots` (`cxx:560-592`).
pub(super) fn unite_knots(
    degree: i32,
    knots1: &[f64],
    mults1: &[i32],
    knots2: &[f64],
    mults2: &[i32],
) -> Option<(Vec<f64>, Vec<i32>)> {
    let (start, end) = same_knot_range(knots1, knots2)?;
    let mut entries = Vec::new();
    collect_knots(degree, start, end, knots1, mults1, &mut entries);
    collect_knots(degree, start, end, knots2, mults2, &mut entries);
    let (knots, mults) = compact_knots(&mut entries);
    if knots.is_empty() || mults.is_empty() {
        None
    } else {
        Some((knots, mults))
    }
}

/// `denominatorValue` (`cxx:266-291`).
pub(super) fn denominator_value(surface: &GeomBSplineSurface, u: f64, v: f64) -> f64 {
    let (w, _) = homogeneous_d0(
        u,
        v,
        0,
        0,
        &surface.poles,
        surface.weights.as_deref(),
        surface.uknots(),
        surface.vknots(),
        Some(surface.umults()),
        Some(surface.vmults()),
        surface.deg_u() as i32,
        surface.deg_v() as i32,
        true,
        true,
        false,
        false,
    );
    w
}

/// Homogeneous numerator point from `BSplSLib::HomogeneousD0`
/// (`cxx:266-290` companion).
pub(super) fn homogeneous_numerator(surface: &GeomBSplineSurface, u: f64, v: f64) -> (f64, GpPnt) {
    homogeneous_d0(
        u,
        v,
        0,
        0,
        &surface.poles,
        surface.weights.as_deref(),
        surface.uknots(),
        surface.vknots(),
        Some(surface.umults()),
        Some(surface.vmults()),
        surface.deg_u() as i32,
        surface.deg_v() as i32,
        true,
        true,
        false,
        false,
    )
}

/// `makeTripleProductBasis` (`cxx:325-360`).
pub(super) fn make_triple_product_basis(
    knots: &[f64],
    mults: &[i32],
    degree: i32,
) -> Option<(i32, Vec<f64>, Vec<i32>, Vec<f64>)> {
    let product_degree = 3 * degree;
    if product_degree as usize > MAX_DEGREE {
        return None;
    }
    let mut product_knots = knots.to_vec();
    let mut product_mults = vec![0i32; mults.len()];
    for knot_idx in 0..knots.len() {
        product_knots[knot_idx] = knots[knot_idx];
        if knot_idx == 0 || knot_idx + 1 == knots.len() {
            product_mults[knot_idx] = product_degree + 1;
        } else {
            product_mults[knot_idx] = product_degree.min(2 * degree + mults[knot_idx]);
        }
    }
    let product_flat = knot_sequence(&product_knots, &product_mults, product_degree);
    Some((product_degree, product_knots, product_mults, product_flat))
}

/// Local port of `BSplSLib::FunctionMultiply` used by `multiplyByDenominators`
/// (`cxx:362-392`, evaluator `cxx:294-323`).
///
/// Samples the surface numerator/denominator times `W_first * W_second` at
/// Schoenberg points of the product basis, then interpolates.
pub(super) fn multiply_by_denominators(
    surface: &GeomBSplineSurface,
    first_denominator: &GeomBSplineSurface,
    second_denominator: &GeomBSplineSurface,
    u_degree: i32,
    v_degree: i32,
    u_flat_knots: &[f64],
    v_flat_knots: &[f64],
) -> Option<(Vec<Vec<GpPnt>>, Vec<Vec<f64>>)> {
    let nb_u = u_flat_knots.len() as i32 - u_degree - 1;
    let nb_v = v_flat_knots.len() as i32 - v_degree - 1;
    if nb_u < 2 || nb_v < 2 {
        return None;
    }
    let nb_u = nb_u as usize;
    let nb_v = nb_v as usize;
    let u_params = schoenberg_points(u_degree as usize, u_flat_knots, nb_u);
    let v_params = schoenberg_points(v_degree as usize, v_flat_knots, nb_v);

    let mut numerator = vec![vec![GpPnt::zero(); nb_v]; nb_u];
    let mut denominator = vec![vec![0.0; nb_v]; nb_u];
    for (iu, &u) in u_params.iter().enumerate() {
        for (iv, &v) in v_params.iter().enumerate() {
            let factor = denominator_value(first_denominator, u, v)
                * denominator_value(second_denominator, u, v);
            let (w, num) = homogeneous_numerator(surface, u, v);
            numerator[iu][iv] = GpPnt::new(num.x() * factor, num.y() * factor, num.z() * factor);
            denominator[iu][iv] = w * factor;
        }
    }

    // Interpolate homogeneous numerator (3D) and denominator (as weight channel).
    if !interpolate_surface_weighted(
        u_degree,
        v_degree,
        u_flat_knots,
        v_flat_knots,
        &u_params,
        &v_params,
        &mut numerator,
        &mut denominator,
    ) {
        return None;
    }
    Some((numerator, denominator))
}

/// `BSplSLib::Interpolate` with poles and weights (`cxx:1096-1104` call site).
/// Packs `(x,y,z,w)` as dimension-4 rows matching the non-rational
/// [`occt_core::bspl::surface_interpolate::interpolate_surface`] layout.
pub(super) fn interpolate_surface_weighted(
    u_degree: i32,
    v_degree: i32,
    u_flat: &[f64],
    v_flat: &[f64],
    u_params: &[f64],
    v_params: &[f64],
    poles: &mut [Vec<GpPnt>],
    weights: &mut [Vec<f64>],
) -> bool {
    let u_len = u_params.len();
    let v_len = v_params.len();
    if u_len == 0 || v_len == 0 || poles.len() != u_len || weights.len() != u_len {
        return false;
    }
    let dim_u = 4 * u_len;
    let mut points = vec![0.0; v_len * dim_u];
    for ii in 0..v_len {
        for jj in 0..u_len {
            let p = poles[jj][ii];
            let w = weights[jj][ii];
            let ll = jj * 4;
            points[ii * dim_u + ll] = p.x();
            points[ii * dim_u + ll + 1] = p.y();
            points[ii * dim_u + ll + 2] = p.z();
            points[ii * dim_u + ll + 3] = w;
        }
    }
    if interpolate_contact(
        v_degree as usize,
        v_flat,
        v_params,
        &vec![0i32; v_len],
        &mut points,
        dim_u,
    )
    .is_err()
    {
        return false;
    }
    let dim_v = 4 * v_len;
    let mut iso = vec![0.0; u_len * dim_v];
    for ii in 0..u_len {
        let kk = ii * 4;
        for jj in 0..v_len {
            let ll = jj * 4;
            iso[ii * dim_v + ll] = points[jj * dim_u + kk];
            iso[ii * dim_v + ll + 1] = points[jj * dim_u + kk + 1];
            iso[ii * dim_v + ll + 2] = points[jj * dim_u + kk + 2];
            iso[ii * dim_v + ll + 3] = points[jj * dim_u + kk + 3];
        }
    }
    if interpolate_contact(
        u_degree as usize,
        u_flat,
        u_params,
        &vec![0i32; u_len],
        &mut iso,
        dim_v,
    )
    .is_err()
    {
        return false;
    }
    for ii in 0..u_len {
        for jj in 0..v_len {
            let ll = jj * 4;
            poles[ii][jj] = GpPnt::new(
                iso[ii * dim_v + ll],
                iso[ii * dim_v + ll + 1],
                iso[ii * dim_v + ll + 2],
            );
            weights[ii][jj] = iso[ii * dim_v + ll + 3];
        }
    }
    true
}

/// Build a surface from poles/weights/knots (`cxx:433-440`, `:482-489`, etc.).
pub(super) fn make_surface(
    poles: Vec<Vec<GpPnt>>,
    weights: Option<Vec<Vec<f64>>>,
    u_knots: Vec<f64>,
    v_knots: Vec<f64>,
    u_mults: Vec<i32>,
    v_mults: Vec<i32>,
    deg_u: usize,
    deg_v: usize,
) -> Option<GeomBSplineSurface> {
    let mut surface =
        GeomBSplineSurface::from_poles_knots_mults(poles, u_knots, v_knots, u_mults, v_mults, deg_u, deg_v)
            .ok()?;
    if let Some(w) = weights {
        let rational = w.iter().any(|row| row.iter().any(|ww| (*ww - 1.0).abs() > 0.0));
        if rational {
            surface.weights = Some(w);
        }
    }
    Some(surface)
}

/// Guard used by weight-division after reference interpolation (`cxx:1115`).
#[inline]
pub(super) fn gp_resolution() -> f64 {
    RESOLUTION
}
