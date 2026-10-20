use occt_geom::Curve;

use super::set::ParamSet;
use super::splitter::RangeSplitter;
use super::surface_type::{mesh_surface_type, MeshSurfaceType};
use super::super::prelude::*;

/// Clamps a discrete range to a surface's geometric range, handling periodic
/// surfaces (whose discrete range may cross the seam). Source:
/// `BRepMesh_DefaultRangeSplitter::updateRange`.
pub(in crate::meshing::range_splitter) fn update_range(
    geom_first: f64,
    geom_last: f64,
    is_periodic: bool,
    d_first: &mut f64,
    d_last: &mut f64,
) {
    if *d_first < geom_first || *d_last > geom_last {
        if is_periodic {
            if (*d_last - *d_first) > (geom_last - geom_first) {
                *d_last = *d_first + (geom_last - geom_first);
            }
        } else if (*d_first < geom_last) && (*d_last > geom_first) {
            // Protection against pcurves out of the surface domain (OCCT #23675).
            if geom_first > *d_first {
                *d_first = geom_first;
            }
            if geom_last < *d_last {
                *d_last = geom_last;
            }
        }
    }
}

/// Adds `param` to the set when it lies within `range`. Returns true on insert.
pub(in crate::meshing::range_splitter) fn add_param(
    param: f64,
    range: (f64, f64),
    params: &mut ParamSet,
) -> bool {
    if param < range.0 || param > range.1 {
        return false;
    }
    params.insert(param);
    true
}

/// Initializes a parameter set from breakpoint intervals, optionally splitting
/// each interval at its midpoint. Source: `initParamsFromIntervals`.
pub(in crate::meshing::range_splitter) fn init_params_from_intervals(
    intervals: &[f64],
    range: (f64, f64),
    is_split: bool,
    params: &mut ParamSet,
) -> bool {
    let mut added = false;
    for i in 0..intervals.len() {
        if add_param(intervals[i], range, params) {
            added = true;
        }
        if is_split && i + 1 < intervals.len() {
            let mid = 0.5 * (intervals[i] + intervals[i + 1]);
            if add_param(mid, range, params) {
                added = true;
            }
        }
    }
    added
}

/// `GeomAdaptor_Curve::NbIntervals` / `Intervals` on `[range]` (`cxx:371-413`).
/// `C0` or `S <= Continuity` is one span. Otherwise keep unique knots that
/// lie inside the trimmed range (CN uses `aCont = degree`, every knot breaks).
fn adaptor_curve_intervals(curve: &dyn Curve, continuity: u8, range: (f64, f64)) -> Vec<f64> {
    // `GeomAdaptor_Curve::NbIntervals` (`cxx:371-413`): BSpline uses
    // `LocalContinuity(First, Last)` then `BSplCLib::Intervals`.
    if let (Some(knots), Some(degree)) = (curve.bspline_knots(), curve.nurbs_degree()) {
        return occt_core::bspl::adaptor_intervals(
            knots,
            degree,
            curve.is_periodic(),
            continuity,
            range.0,
            range.1,
            curve.resolution(CONFUSION).min(PCONFUSION),
        );
    }
    if continuity == 0 || (!curve.is_periodic() && continuity <= curve.continuity()) {
        return vec![range.0, range.1];
    }
    let raw = curve.parameter_intervals(continuity);
    let mut iv = vec![range.0];
    for &k in &raw {
        if k > range.0 + 1e-14 && k < range.1 - 1e-14 && iv.last().map_or(true, |p| (k - p).abs() > 1e-14)
        {
            iv.push(k);
        }
    }
    if iv.last().map_or(true, |p| (range.1 - p).abs() > 1e-14) {
        iv.push(range.1);
    }
    iv
}

/// `GeomAdaptor_Surface` Offset interval shape (`cxx:664-682`, `722-740`).
/// C0->C1, C1->C2, C2->C3; C3/CN keep CN. Mesh `initParameters` passes CN.
fn offset_interval_continuity(s: u8) -> u8 {
    match s {
        0 => 2,
        2 => 4,
        4 => 5,
        _ => 6,
    }
}

/// `GeomAdaptor_Surface::NbUIntervals` / `UIntervals` (`cxx:643-818`).
/// BSpline uses `VIso(FirstVKnot)` trimmed to the face U range; V uses `UIso`.
/// Offset forwards to the basis adaptor (`cxx:664-684`), never Offset UIso.
/// Revolution U is `[UFirst,ULast]`; V uses a BSpline generatrix (`cxx:712-720`).
fn adaptor_surface_intervals(
    surf: &dyn Surface,
    is_u: bool,
    continuity: u8,
    range: (f64, f64),
) -> Vec<f64> {
    match mesh_surface_type(surf) {
        MeshSurfaceType::OffsetSurface => {
            let base_s = offset_interval_continuity(continuity);
            if let Some(basis) = surf.offset_basis_surface() {
                return adaptor_surface_intervals(basis.as_ref(), is_u, base_s, range);
            }
        }
        MeshSurfaceType::SurfaceOfRevolution => {
            if !is_u {
                if let Some(c) = surf.revolution_basis_curve() {
                    if c.bspline_knots().is_some() {
                        return adaptor_curve_intervals(c.as_ref(), continuity, range);
                    }
                }
            }
            return vec![range.0, range.1];
        }
        MeshSurfaceType::BSplineSurface => {}
        _ => {
            return if is_u {
                surf.u_intervals(continuity)
            } else {
                surf.v_intervals(continuity)
            };
        }
    }
    let iso = if is_u {
        let v0 = surf
            .v_intervals(0)
            .first()
            .copied()
            .unwrap_or(surf.v_range().0);
        surf.v_iso_curve(v0)
    } else {
        let u0 = surf
            .u_intervals(0)
            .first()
            .copied()
            .unwrap_or(surf.u_range().0);
        surf.u_iso_curve(u0)
    };
    match iso {
        Some(c) => adaptor_curve_intervals(c.as_ref(), continuity, range),
        None => {
            if is_u {
                surf.u_intervals(continuity)
            } else {
                surf.v_intervals(continuity)
            }
        }
    }
}

/// `BRepMesh_NURBSRangeSplitter::getUndefinedInterval` (`cxx:414-449`).
///
/// The interval count and break array come from the untrimmed surface adaptor
/// (`cxx:420-421`, `cxx:438-447`): `GetSurface()` is the `BRepAdaptor_Surface`
/// built with `R = false` by `IMeshData_Face` (`IMeshData_Face.hxx`), so
/// `NbUIntervals` / `UIntervals` see the full `Geom_Surface` bounds. The
/// adjusted range `GetRangeU()` (`cxx:460-461`) is only fed to the uniform
/// fallback grid (`cxx:430`).
pub(in crate::meshing::range_splitter) fn undefined_interval(
    s: &dyn RangeSplitter,
    is_u: bool,
    continuity: u8,
    range: (f64, f64),
) -> Vec<f64> {
    let Some(surf) = s.surface() else {
        return vec![range.0, range.1];
    };
    let geom_range = if is_u { surf.u_range() } else { surf.v_range() };
    let iv = adaptor_surface_intervals(surf.as_ref(), is_u, continuity, geom_range);
    if std::env::var("NURBSGRID").is_ok() {
        eprint!("NURBSGRID iv is_u={is_u} geom=[{:.17e},{:.17e}] n={} vals=", geom_range.0, geom_range.1, iv.len());
        for (i, &v) in iv.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{v:.17e}");
        }
        eprintln!();
    }
    let mut intervals_nb = iv.len().saturating_sub(1) as i32;
    if intervals_nb == 1 {
        intervals_nb = s.undefined_interval_nb(is_u, continuity);
        if intervals_nb > 1 {
            let diff = (range.1 - range.0) / intervals_nb as f64;
            return (1..intervals_nb)
                .map(|i| range.0 + i as f64 * diff)
                .collect();
        }
    }
    if iv.is_empty() {
        vec![range.0, range.1]
    } else {
        iv
    }
}

/// `toSplitIntervals`: split when `GeomLib::NormEstim` fails at a corner.
pub(in crate::meshing::range_splitter) fn to_split_intervals(
    surface: &dyn Surface,
    u_iv: &[f64],
    v_iv: &[f64],
) -> bool {
    use occt_core::precision::CONFUSION;
    use crate::meshing::geomlib_norm::norm_estim;
    for &u in u_iv {
        if !u.is_finite() {
            continue;
        }
        for &v in v_iv {
            if !v.is_finite() {
                continue;
            }
            if norm_estim(surface, u, v, CONFUSION).0 != 0 {
                return true;
            }
        }
    }
    false
}
