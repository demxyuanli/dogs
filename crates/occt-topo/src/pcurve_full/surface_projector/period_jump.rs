//! Seam / degeneracy jump resolution on a closed surface - the `needResolveUJump` / `needResolveVJump` and `AdjustOverDegen` arms of `ShapeConstruct_ProjectCurveOnSurface::Perform`.

use super::*;

/// `approxPCurve` (`ShapeConstruct_ProjectCurveOnSurface.cxx:1499-1690`):
/// "Handle U-closed surfaces" / "Handle V-closed surfaces". A closed (not
/// necessarily periodic) surface's point-wise projection returns every sample
/// reduced into the parameter window, so a curve that runs once around the
/// closed direction carries a full-period jump in the middle of its `(u, v)`
/// sequence. Interpolating through that sequence produces a curve that sweeps
/// the whole window instead of running along the surface (measured on
/// ATU01038: 0.77 to 16.0 against a 0.01 projection tolerance), and the
/// mesher's UV frontier then degenerates.
///
/// `theIsFromCache` / `aSavedPoint` are produced by the sample loop
/// (`cxx:1295-1465`, ported in [`project_curve_on_surface_perform`]): in the
/// `!isRecompute` endpoint arm they are `getLine`'s `theIsFromCache` / the
/// `thePoints2d(1)` projection (`cxx:1390-1399`), and in the `isRecompute`
/// endpoint arm the `myCache` lookup sets them from the cache anchor
/// (`cxx:1401-1424`). `is_from_cache` switches off the "pull the first sample
/// into the window" and the final recentring steps (`cxx:1508`,
/// `cxx:1521-1536`, `cxx:1571`, `cxx:1599`, `cxx:1613-1628`, `cxx:1664`) so
/// that the cached parameter branch is preserved.
///
/// NOT PORTED HERE (reported): the `isoParam` arms of the sampling loop
/// (`cxx:1317-1375`), so this function's inputs are always the `else` arm of
/// `cxx:1376-1385`.
///
/// The `needResolveUJump` / `needResolveVJump` flags (`cxx:1436-1451`, latched
/// in the sampling loop of the caller) arrive as `need_resolve_u_jump` /
/// `need_resolve_v_jump` and widen the two triggers below exactly as in OCCT
/// (`cxx:1504`, `cxx:1596-1598`).
///
/// The `Handle AdjustOverDegen` block (`cxx:1746-1890`) runs at the end of this
/// function, after the V-closed seam adjustment, exactly as in OCCT.
pub(super) fn resolve_closed_surface_period_jump(
    surf: &dyn Surface,
    curve: &dyn Curve,
    singularities: &[super::super::singularities::Singularity],
    preci: f64,
    need_resolve_u_jump: bool,
    need_resolve_v_jump: bool,
    is_from_cache: bool,
    saved_point: GpPnt2d,
    adjust_off: bool,
    params: &mut Vec<f64>,
    pts3d: &mut Vec<GpPnt>,
    pts2d: &mut Vec<GpPnt2d>,
) -> bool {
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let up = ul - uf;
    let vp = vl - vf;
    // `cxx:1501-1502`
    let tol_on_u_period = occt_core::precision::CONFUSION * up;
    let tol_on_v_period = occt_core::precision::CONFUSION * vp;

    if sa_is_u_closed(surf, preci) || need_resolve_u_jump {
        // `cxx:1505-1519`: pull the first sample into `[uf, ul]` unless it came
        // from the cache, in which case that representative is kept.
        let mut first_x = pts2d[0].x();
        if !is_from_cache {
            while first_x < uf {
                first_x += up;
                pts2d[0].set_x(first_x);
            }
            while first_x > ul {
                first_x -= up;
                pts2d[0].set_x(first_x);
            }
        }
        // `cxx:1521-1536`: with a cached anchor the period window is placed
        // around it and the first sample is folded into that window.
        if surf.is_u_periodic() && is_from_cache {
            let mut min_param = uf;
            let mut max_param = ul;
            while min_param > saved_point.x() {
                min_param -= up;
                max_param -= up;
            }
            while max_param < saved_point.x() {
                min_param += up;
                max_param += up;
            }
            first_x += super::super::projection_cache::adjust_to_period(first_x, min_param, max_param);
            pts2d[0].set_x(first_x);
        }
        let mut prev_x = first_x;
        let mut min_x = first_x;
        let mut max_x = first_x;
        let mut to_adjust = false;
        // `cxx:1544-1569`
        let mut i = 1usize;
        while i < pts2d.len() {
            let raw_x = pts2d[i].x();
            let mut cur_x = raw_x;
            if (cur_x - prev_x).abs() > 0.5 * up {
                insert_additional_point_or_adjust(
                    &mut to_adjust,
                    1,
                    up,
                    tol_on_u_period,
                    raw_x,
                    &mut cur_x,
                    prev_x,
                    preci,
                    curve,
                    &mut i,
                    params,
                    pts3d,
                    pts2d,
                    surf,
                );
            }
            prev_x = cur_x;
            if min_x > cur_x {
                min_x = cur_x;
            } else if max_x < cur_x {
                max_x = cur_x;
            }
            i += 1;
        }
        // `cxx:1571-1591`: recenter the sequence on the surface window. Skipped
        // when the anchor came from the cache (`cxx:1571`).
        if !is_from_cache {
            let mid_x = 0.5 * (min_x + max_x);
            let shift_x = if mid_x > ul {
                -up
            } else if mid_x < uf {
                up
            } else {
                0.0
            };
            if shift_x != 0.0 {
                for p in pts2d.iter_mut() {
                    p.set_x(p.x() + shift_x);
                }
            }
        }
    }

    // `cxx:1596-1598`; the `Geom_SphericalSurface` arm is the extra trigger
    // OCCT adds because a sphere is closed in V without being V-periodic.
    if sa_is_v_closed(surf, preci) || need_resolve_v_jump || surf.gp_sphere().is_some() {
        // `cxx:1599-1611`
        let mut first_y = pts2d[0].y();
        if !is_from_cache {
            while first_y < vf {
                first_y += vp;
                pts2d[0].set_y(first_y);
            }
            while first_y > vl {
                first_y -= vp;
                pts2d[0].set_y(first_y);
            }
        }
        // `cxx:1613-1628`
        if surf.is_v_periodic() && is_from_cache {
            let mut min_param = vf;
            let mut max_param = vl;
            while min_param > saved_point.y() {
                min_param -= vp;
                max_param -= vp;
            }
            while max_param < saved_point.y() {
                min_param += vp;
                max_param += vp;
            }
            first_y += super::super::projection_cache::adjust_to_period(first_y, min_param, max_param);
            pts2d[0].set_y(first_y);
        }
        let mut prev_y = first_y;
        let mut min_y = first_y;
        let mut max_y = first_y;
        let mut to_adjust = false;
        // `cxx:1636-1662`
        let mut i = 1usize;
        while i < pts2d.len() {
            let raw_y = pts2d[i].y();
            let mut cur_y = raw_y;
            if (cur_y - prev_y).abs() > 0.5 * vp {
                insert_additional_point_or_adjust(
                    &mut to_adjust,
                    2,
                    vp,
                    tol_on_v_period,
                    raw_y,
                    &mut cur_y,
                    prev_y,
                    preci,
                    curve,
                    &mut i,
                    params,
                    pts3d,
                    pts2d,
                    surf,
                );
            }
            prev_y = cur_y;
            if min_y > cur_y {
                min_y = cur_y;
            } else if max_y < cur_y {
                max_y = cur_y;
            }
            i += 1;
        }
        // `cxx:1664-1685`. Skipped when the anchor came from the cache.
        if !is_from_cache {
            let mid_y = 0.5 * (min_y + max_y);
            let shift_y = if mid_y > vl {
                -vp
            } else if mid_y < vf {
                vp
            } else {
                0.0
            };
            if shift_y != 0.0 {
                for p in pts2d.iter_mut() {
                    p.set_y(p.y() + shift_y);
                }
            }
        }

        // `cxx:1687-1746`: "Handle V-closed seam adjustment". A pair that still
        // straddles the seam by more than half a period really crosses it, so
        // whichever of the two samples is nearest to a seam boundary is pinned
        // exactly onto it (`vf` or `vl`); leaving it a tolerance away from the
        // seam makes the mesher's UV frontier cut the surface instead of
        // closing along it.
        for i in 1..pts2d.len() {
            let prev_y = pts2d[i - 1].y();
            let cur_y = pts2d[i].y();
            if (cur_y - prev_y).abs() > 0.5 * vp {
                let dist_prev_vf = (prev_y - vf).abs();
                let dist_prev_vl = (prev_y - vl).abs();
                let dist_curr_vf = (cur_y - vf).abs();
                let dist_curr_vl = (cur_y - vl).abs();

                let mut prev_on_first = true;
                let mut prev_on_last = false;
                let mut curr_on_first = false;
                let mut curr_on_last = false;
                let mut min = dist_prev_vf;
                if dist_prev_vl < min {
                    min = dist_prev_vl;
                    prev_on_first = false;
                    prev_on_last = true;
                }
                if dist_curr_vf < min {
                    min = dist_curr_vf;
                    prev_on_first = false;
                    prev_on_last = false;
                    curr_on_first = true;
                }
                if dist_curr_vl < min {
                    prev_on_first = false;
                    prev_on_last = false;
                    curr_on_first = false;
                    curr_on_last = true;
                }

                if prev_on_first {
                    pts2d[i - 1].set_y(vf);
                } else if prev_on_last {
                    pts2d[i - 1].set_y(vl);
                } else if curr_on_first {
                    pts2d[i].set_y(vf);
                } else if curr_on_last {
                    pts2d[i].set_y(vl);
                }
            }
        }
    }

    // `cxx:1746-1890`: "Handle AdjustOverDegen".
    super::super::singularities::adjust_over_degenerated(
        surf,
        singularities,
        preci,
        pts3d,
        pts2d,
        adjust_off,
    )
}

