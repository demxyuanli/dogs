//! Split-point search of `ProjLib_CompProjectedCurve`: `BuildCurveSplits`,
//! `SplitOnDirection` and `FindSplitPoint`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx:2243-2390`
//! (OCCT 8.0.0), plus the `SplitDS` bundle (`cxx:80-116`).
//!
//! The splits are the curve parameters where the projection crosses the
//! periodic seam of the surface, so `Init` can stop the marching walk there.
//! For a periodic direction OCCT builds the 3D image of the seam in parameter
//! space (a `Geom2d_Line` through `(FirstU, FirstV)` along the periodic axis,
//! carried onto the surface by `Adaptor3d_CurveOnSurface`), runs
//! `Extrema_ExtCC` against the curve to be projected, and keeps only the
//! extrema whose point projects back onto the seam:
//! `Extrema_ExtPS` is run on the curve point and its nearest solution is
//! required to sit on `First`/`Last` of the periodic parameter.
//!
//! Translation notes:
//! - `Extrema_ExtCC` with `SetCurve(1, seam)`, `SetCurve(2, curve)`,
//!   `SetSingleSolutionFlag(true)`, `SetRange(1, 0, Last2DParam)` and
//!   `SetRange(2, MinParam, MaxParam)` is [`GGenExtCC::new`] with those ranges
//!   plus [`GGenExtCC::set_single_solution_flag`].
//! - `std::sort(..., Comparator)` (`cxx:2272`) is a stable ascending sort.
//! - `anExtPS.IsDone()` with `NbExt() == 0` would index `Point(-1)` in OCCT
//!   (`cxx:2356`); [`ExtPs::is_done`] implies `nb_ext() > 0`, and the port
//!   skips the extremum instead of raising.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::curve::Curve;
use crate::curve_on_surface::CurveOnSurface;
use crate::extrema_cc::GGenExtCC;
use crate::extrema_surf::{ExtPs, ExtremaExtFlag};
use crate::surface::Surface;

/// `BuildCurveSplits(Curve, Surface, TolU, TolV, Splits)` (`cxx:2243-2273`).
pub(crate) fn build_curve_splits(
    curve: &dyn Curve,
    surface: &Arc<dyn Surface>,
    tol_u: f64,
    tol_v: f64,
    splits: &mut Vec<f64>,
) {
    let mut ext_ps = ExtPs::new();
    ext_ps.initialize(
        surface.as_ref(),
        surface.u_range().0,
        surface.u_range().1,
        surface.v_range().0,
        surface.v_range().1,
        tol_u,
        tol_v,
    );
    ext_ps.set_flag(ExtremaExtFlag::Min);

    if surface.is_u_periodic() {
        split_on_direction(0, curve, surface, &mut ext_ps, splits);
    }
    if surface.is_v_periodic() {
        split_on_direction(1, curve, surface, &mut ext_ps, splits);
    }
    splits.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
}

/// `SplitOnDirection(SplitDS&)` (`cxx:2276-2308`).
fn split_on_direction(
    periodic_dir: i32,
    curve: &dyn Curve,
    surface: &Arc<dyn Surface>,
    ext_ps: &mut ExtPs<'_>,
    splits: &mut Vec<f64>,
) {
    let start = GpPnt2d::new(surface.u_range().0, surface.v_range().0);
    // `gp_Dir2d(theSplitDS.myPeriodicDir, (int)!theSplitDS.myPeriodicDir)`
    // (`cxx:2286`).
    let dir = if periodic_dir == 0 {
        GpDir2d::new(1.0, 0.0).expect("unit direction")
    } else {
        GpDir2d::new(0.0, 1.0).expect("unit direction")
    };

    let per_min_param = if periodic_dir == 0 {
        surface.u_range().0
    } else {
        surface.v_range().0
    };
    let per_max_param = if periodic_dir == 0 {
        surface.u_range().1
    } else {
        surface.v_range().1
    };
    // `aLast2DParam` (`cxx:2293-2296`).
    let last_2d_param = if periodic_dir == 1 {
        surface.u_range().1 - surface.u_range().0
    } else {
        surface.v_range().1 - surface.v_range().0
    };

    // `Geom2d_Line(aStartPnt, aDir)` + `Geom2dAdaptor_Curve(aC2GC, 0,
    // aLast2DParam)` + `Adaptor3d_CurveOnSurface` (`cxx:2299-2301`).
    let a_c2gc: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(start, dir));
    let a_c_on_s = CurveOnSurface::new(a_c2gc, Arc::clone(surface));

    find_split_point(
        curve,
        &a_c_on_s,
        ext_ps,
        periodic_dir,
        per_min_param,
        per_max_param,
        last_2d_param,
        splits,
        curve.first_parameter(),
        curve.last_parameter(),
    );
}

/// `FindSplitPoint(SplitDS&, theMinParam, theMaxParam)` (`cxx:2311-2390`).
#[allow(clippy::too_many_arguments)]
fn find_split_point(
    curve: &dyn Curve,
    border: &CurveOnSurface,
    ext_ps: &mut ExtPs<'_>,
    periodic_dir: i32,
    per_min_param: f64,
    per_max_param: f64,
    last_2d_param: f64,
    splits: &mut Vec<f64>,
    min_param: f64,
    max_param: f64,
) {
    // A fresh `Extrema_ExtCC` per recursion level (`cxx:2314-2322`).
    let mut ext_cc = GGenExtCC::new(
        border,
        curve,
        (0.0, min_param),
        (last_2d_param, max_param),
    );
    ext_cc.set_single_solution_flag(true);
    let _ = ext_cc.perform();

    if ext_cc.is_done() && !ext_cc.is_parallel() {
        let nb_ext = ext_cc.nb_ext();
        for idx in 1..=nb_ext {
            let (_, _, a_param, p2) = ext_cc.points(idx);
            ext_ps.perform(&p2);
            if !ext_ps.is_done() {
                return;
            }

            // Nearest `Extrema_ExtPS` solution, to avoid false positives
            // (`cxx:2338-2355`).
            let nb_pext = ext_ps.nb_ext();
            let mut a_min_idx = 0usize;
            let mut a_min_sq_dist = f64::MAX;
            for a_p_idx in 1..=nb_pext {
                let a_curr_sq_dist = ext_ps.square_distance(a_p_idx);
                if a_curr_sq_dist < a_min_sq_dist {
                    a_min_sq_dist = a_curr_sq_dist;
                    a_min_idx = a_p_idx;
                }
            }
            if a_min_idx == 0 {
                continue;
            }

            // Check that the point is projected onto the periodic border
            // (`cxx:2357-2367`).
            let (u, v, _) = ext_ps.point(a_min_idx);
            let a_proj_param = if periodic_dir != 0 { v } else { u };
            if (a_proj_param - per_min_param).abs() < occt_core::precision::PCONFUSION
                || (a_proj_param - per_max_param).abs() < occt_core::precision::PCONFUSION
            {
                let a_cf_param = curve.first_parameter();
                let a_cl_param = curve.last_parameter();

                if a_param > a_cf_param + occt_core::precision::PCONFUSION
                    && a_param < a_cl_param - occt_core::precision::PCONFUSION
                {
                    // Add only inner points (`cxx:2368-2372`).
                    splits.push(a_param);
                }

                let a_delta_coeff = 0.01;
                let a_delta = (max_param - min_param + a_cl_param - a_cf_param) * a_delta_coeff;

                if a_param - a_delta > min_param + occt_core::precision::PCONFUSION {
                    find_split_point(
                        curve,
                        border,
                        ext_ps,
                        periodic_dir,
                        per_min_param,
                        per_max_param,
                        last_2d_param,
                        splits,
                        min_param,
                        a_param - a_delta,
                    );
                }
                if a_param + a_delta < max_param - occt_core::precision::PCONFUSION {
                    find_split_point(
                        curve,
                        border,
                        ext_ps,
                        periodic_dir,
                        per_min_param,
                        per_max_param,
                        last_2d_param,
                        splits,
                        a_param + a_delta,
                        max_param,
                    );
                }
            }
        }
    }
}
