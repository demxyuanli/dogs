//! Sampling-count helpers for 2D/3D bounding-box evaluation.
//!
//! Source: `GeomBndLib_SamplingHelpers.pxx` (`ComputeNbSamples2d`,
//! `ComputeNbUSamples`, `ComputeNbVSamples`). PerformAreas uses
//! `GeomBndLib_OtherCurve2d::Box` (`N = 33`) rather than these counts; they
//! are the sample budget of `BoxOptimal` and of surface-box helpers, kept
//! here so a later `BoxOptimal` port does not invent a different `N`.

use occt_geom2d::curve::Curve2d;

/// Curve kind used only to pick the sample budget, matching `GeomAbs_CurveType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleCurveKind {
    Bezier,
    BSpline,
    Other,
}

/// `RealToInt` used by the OCCT helper: truncate toward zero after adding 0.5
/// is *not* what OCCT does; `RealToInt` is a C-style cast. Truncate toward zero.
fn real_to_int(x: f64) -> i32 {
    x as i32
}

/// `GeomBndLib_SamplingHelpers::ComputeNbSamples2d`.
///
/// * Bezier: `N = 2 * NbPoles`; if `(UMax-UMin) < 0.9` scale and floor at 5.
/// * BSpline: `N = 2 * (Degree+1) * (NbKnots-1)` then the same span scale.
/// * Other: `N = 17`.
/// * Finally `min(23, N)`.
pub fn compute_nb_samples2d(
    kind: SampleCurveKind,
    nb_poles: usize,
    degree: usize,
    nb_knots: usize,
    first: f64,
    last: f64,
    the_u_min: f64,
    the_u_max: f64,
) -> i32 {
    let mut n = match kind {
        SampleCurveKind::Bezier => {
            let mut n = 2 * nb_poles as i32;
            let du = the_u_max - the_u_min;
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleCurveKind::BSpline => {
            let mut n = 2 * (degree as i32 + 1) * (nb_knots as i32 - 1).max(0);
            let umin = first;
            let umax = last;
            let den = umax - umin;
            let du = if den.abs() > f64::EPSILON {
                (the_u_max - the_u_min) / den
            } else {
                1.0
            };
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleCurveKind::Other => 17,
    };
    n = n.min(23);
    n
}

/// Surface kind for U/V sample counts (`GeomAbs_SurfaceType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleSurfaceKind {
    Bezier,
    BSpline,
    Other,
}

/// `GeomBndLib_SamplingHelpers::ComputeNbUSamples`.
///
/// * Bezier: `N = 2 * NbUPoles`; span `< 0.9` scales and floors at 5.
/// * BSpline: `N = 2 * (UDegree+1) * (NbUKnots-1)` then the same span scale.
/// * Other: `N = 33`.
/// * Finally `min(50, N)`.
pub fn compute_nb_u_samples(
    kind: SampleSurfaceKind,
    nb_u_poles: usize,
    u_degree: usize,
    nb_u_knots: usize,
    umin: f64,
    umax: f64,
    the_u_min: f64,
    the_u_max: f64,
) -> i32 {
    let n = match kind {
        SampleSurfaceKind::Bezier => {
            let mut n = 2 * nb_u_poles as i32;
            let du = the_u_max - the_u_min;
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleSurfaceKind::BSpline => {
            let mut n = 2 * (u_degree as i32 + 1) * (nb_u_knots as i32 - 1).max(0);
            let den = umax - umin;
            let du = if den.abs() > f64::EPSILON {
                (the_u_max - the_u_min) / den
            } else {
                1.0
            };
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleSurfaceKind::Other => 33,
    };
    n.min(50)
}

/// `GeomBndLib_SamplingHelpers::ComputeNbVSamples` — V-axis copy of the U table
/// (`NbVPoles` / `VDegree` / `NbVKnots`, span over `[vmin, vmax]`).
pub fn compute_nb_v_samples(
    kind: SampleSurfaceKind,
    nb_v_poles: usize,
    v_degree: usize,
    nb_v_knots: usize,
    vmin: f64,
    vmax: f64,
    the_v_min: f64,
    the_v_max: f64,
) -> i32 {
    let n = match kind {
        SampleSurfaceKind::Bezier => {
            let mut n = 2 * nb_v_poles as i32;
            let dv = the_v_max - the_v_min;
            if dv < 0.9 {
                n = real_to_int(dv * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleSurfaceKind::BSpline => {
            let mut n = 2 * (v_degree as i32 + 1) * (nb_v_knots as i32 - 1).max(0);
            let den = vmax - vmin;
            let dv = if den.abs() > f64::EPSILON {
                (the_v_max - the_v_min) / den
            } else {
                1.0
            };
            if dv < 0.9 {
                n = real_to_int(dv * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleSurfaceKind::Other => 33,
    };
    n.min(50)
}

/// Classify a 2D curve the way `Geom2dAdaptor_Curve::GetType()` does for
/// sample budgets (`Geom2dAdaptor_Curve.cxx:285-345`): a `Geom2d_TrimmedCurve`
/// takes its basis's type (`cxx:285-288`), `GeomAbs_BezierCurve` and
/// `GeomAbs_BSplineCurve` are reported as such, and every other curve type
/// falls to `GeomAbs_OtherCurve` (`cxx:343-345`), which
/// `ComputeNbSamples2d` samples 17 times (`GeomBndLib_SamplingHelpers.pxx:103-104`).
pub fn sample_kind_of(curve: &dyn Curve2d) -> SampleCurveKind {
    if let Some(basis) = curve.trimmed_basis() {
        return sample_kind_of(basis);
    }
    if curve.bezier_nb_poles().is_some() {
        SampleCurveKind::Bezier // `cxx:316-320`
    } else if curve.bspline_nb_knots().is_some() {
        SampleCurveKind::BSpline // `cxx:322-328`
    } else {
        SampleCurveKind::Other
    }
}
pub fn compute_nb_u_samples_full(kind: SampleSurfaceKind, nb_u_poles: usize, u_degree: usize, nb_u_knots: usize) -> i32 {
    let n = match kind {
        SampleSurfaceKind::Bezier => 2 * nb_u_poles as i32,
        SampleSurfaceKind::BSpline => 2 * (u_degree as i32 + 1) * (nb_u_knots as i32 - 1).max(0),
        SampleSurfaceKind::Other => 33,
    };
    n.min(50)
}

/// `GeomBndLib_SamplingHelpers::ComputeNbVSamples` (full-range variant).
pub fn compute_nb_v_samples_full(kind: SampleSurfaceKind, nb_v_poles: usize, v_degree: usize, nb_v_knots: usize) -> i32 {
    compute_nb_u_samples_full(kind, nb_v_poles, v_degree, nb_v_knots)
}

/// `GeomBndLib_SamplingHelpers::ComputeNbSamples` (3D curve).
///
/// Bezier: `2 * NbPoles`, span `< 0.9` scales and floors at 5.
/// BSpline: `2 * (Degree+1) * (NbKnots-1)` then the same span scale.
/// Other: `N = 33`. Finally `min(500, N)`.
pub fn compute_nb_samples_3d(
    kind: SampleCurveKind,
    nb_poles: usize,
    degree: usize,
    nb_knots: usize,
    first: f64,
    last: f64,
    the_u_min: f64,
    the_u_max: f64,
) -> i32 {
    let n = match kind {
        SampleCurveKind::Bezier => {
            let mut n = 2 * nb_poles as i32;
            let du = the_u_max - the_u_min;
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleCurveKind::BSpline => {
            let mut n = 2 * (degree as i32 + 1) * (nb_knots as i32 - 1).max(0);
            let den = last - first;
            let du = if den.abs() > f64::EPSILON {
                (the_u_max - the_u_min) / den
            } else {
                1.0
            };
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SampleCurveKind::Other => 33,
    };
    n.min(500)
}

/// `GeomBndLib_SamplingHelpers::ComputeNbSamplesT` for a 2D adaptor.
pub fn compute_nb_samples_t_2d(
    kind: SampleCurveKind,
    nb_poles: usize,
    degree: usize,
    nb_knots: usize,
    first: f64,
    last: f64,
    u1: f64,
    u2: f64,
) -> i32 {
    compute_nb_samples2d(kind, nb_poles, degree, nb_knots, first, last, u1, u2)
}

/// `GeomBndLib_SamplingHelpers::ComputeNbSamplesT` for a 3D adaptor.
pub fn compute_nb_samples_t_3d(
    kind: SampleCurveKind,
    nb_poles: usize,
    degree: usize,
    nb_knots: usize,
    first: f64,
    last: f64,
    u1: f64,
    u2: f64,
) -> i32 {
    compute_nb_samples_3d(kind, nb_poles, degree, nb_knots, first, last, u1, u2)
}
