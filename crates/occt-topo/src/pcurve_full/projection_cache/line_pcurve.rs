//! `getLine` (`cxx:902-1102`) plus `MaxKnotMult` and `isCn`.

use super::super::prelude::*;
use super::*;

/// `MaxKnotMult(internal mults)` over the flat knot vector: the multiplicity of
/// the most repeated *interior* knot (`BSplCLib::MaxKnotMult`,
/// `Geom_BSplineCurve_1.cxx:54-55`).
fn max_internal_mult(flat: &[f64]) -> i32 {
    let (knots, mults) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(flat);
    if knots.len() <= 2 {
        return 0;
    }
    mults[1..knots.len() - 1].iter().copied().max().unwrap_or(0)
}

/// `Geom_Surface::IsCNu(1) && Geom_Surface::IsCNv(1)` (`cxx:1035`). The
/// `Geom_Surface` base answers `Standard_True`; `Geom_BSplineSurface` forwards
/// to `Geom_BSplineCurve::IsCN(1)` (`Geom_BSplineCurve_1.cxx:34-59`), true iff
/// the smoothness is at least `C1`, i.e. `degree - MaxKnotMult(internal) >= 1`.
/// A `Geom_BezierSurface` is polynomial on its span, so `N = 1` holds.
fn is_cn_1(s: &dyn Surface) -> bool {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return is_cn_1(basis.as_ref());
    }
    if let Some(basis) = s.offset_basis_surface() {
        return is_cn_1(basis.as_ref());
    }
    if !s.is_bspline_surface() {
        return true;
    }
    let Some(b) = s.osculating_bspline() else {
        return true;
    };
    let cn = |flat: &[f64], deg: usize| (deg as i32 - max_internal_mult(flat)) >= 1;
    cn(b.knots_u(), b.deg_u()) && cn(b.knots_v(), b.deg_v())
}

/// The projection of `p` (with its residue, OCCT's `SurfaceProjector::Gap`).
/// `SurfaceProjectorWithCache::ValueOfUV` (`cxx:106-115`) is the corner-cache
/// lookup followed by `mySurf->ValueOfUV(thePoint, theTol)`, i.e.
/// [`value_of_uv`] with the working tolerance as `preci`, and `Gap()` (`cxx:141`)
/// is `mySurf->Gap()`; see [`super::surface_projector::value_of_uv_with_gap`] for how the
/// residue is recovered.
///
/// `ShapeConstruct_ProjectCurveOnSurface::getLine` (`cxx:902-1102`): the
/// straight-pcurve shortcut. `points` / `params` are `thePoints` / `theParams`
/// from `generateCurvePoints` (`cxx:750-752`), `tol` is `myPreci`
/// (`cxx:1119` passes it through), `cache` is `myCache` and `out` receives
/// `theIsRecompute` / `theIsFromCache` and the endpoint projections.
/// `None` means "not a straight pcurve" and the caller falls through to
/// `isAnIsoparametric` / the approximation arms.
///
/// KNIFE-EDGE (reported, not fudged): on a sphere whose `V` is NOT periodic
/// (`Geom_SphericalSurface::IsVPeriodic`, `Geom_SphericalSurface.cxx:130-133`)
/// only the U half of `cxx:982-1013` runs, so the verdict of this function for
/// an equator edge is decided entirely by `fix_periodicity_troubles` on U. For
/// a probe whose projected `u` is exactly `0.` the `AdjustByPeriod`
/// (`ShapeAnalysis.cxx:48-60`) shift is `0.` (`D <= 0.5 * P` holds because
/// `0. - PI` is exactly `-UPeriod/2`), the point is pinned to `aMinParam`
/// (`cxx:348-350`) and the following probe folded to `~UPeriod` makes
/// `cxx:355-371` read `D > 0.5 * P` as a sign change -> `theIsRecompute` ->
/// the `cxx:1015-1018` sphere guard -> `nullptr` -> `interpolatePCurve`. For
/// the same probe projected as `-4.1e-13` instead of `0.` the shift is
/// `+UPeriod`, the point pins to `aMaxParam` and no sign change is found, so
/// `getLine` returns the straight pcurve and the entity stays a type 110
/// instead of the type 126 OCCT writes. `Shape-1` has three such edges, and
/// the projected `u` differs only because the edge's 3D range there is
/// `3*PI/2 + 7.2e-13 .. 2*PI + 7.2e-13` rather than OCCT's exact
/// `3*PI/2 .. 2*PI`; the `theParams` used by `generateCurvePoints`
/// (`cxx:531-575`) come straight from the edge, so the sample point lands
/// ~1e-11 off the sphere's principal plane and `atan2` (`ElSLib.cxx:1640-1643`)
/// no longer returns the exact `0.`.
pub(crate) fn line_pcurve_from_points(
    s: &dyn Surface,
    points: &[GpPnt],
    params: &[f64],
    tol: f64,
    cache: &ProjectorCache,
    out: &mut GetLineOut,
) -> Option<Arc<dyn Curve2d>> {
    let nb = points.len();
    if nb < 2 || params.len() != nb {
        return None;
    }
    // `cxx:911-916`: probes 1, 2, nb-1 and nb.
    let p = [points[0], points[1], points[nb - 2], points[nb - 1]];

    let mut tol2 = tol * tol;
    let mut tol_working = tol;
    let is_periodic_u = s.is_u_periodic();
    let is_periodic_v = s.is_v_periodic();

    // `cxx:925-936`: protection against bad tolerance shapes.
    if tol2 > 1.0 {
        tol_working = CONFUSION;
        tol2 = tol_working * tol_working;
    }
    if tol2 < SQUARE_CONFUSION {
        tol2 = SQUARE_CONFUSION;
    }
    let an_old_tol2 = tol2;

    let mut p2d = [GpPnt2d::zero(); 4];

    // `cxx:941`: `SurfaceProjectorWithCache aProjector(mySurf)`.
    let projector = SurfaceProjectorWithCache::new(s);

    // `cxx:944-965` `projectPoint`: an existing endpoint cache hit first
    // (`cxx:946-961`), then the full projection (`cxx:964`).
    let mut saved_point_num: i32 = -1;
    let mut saved_point = GpPnt2d::zero();
    for i in [0usize, 3] {
        let hit = (0..cache.len()).find(|&j| cache.get(j).first.square_distance(&p[i]) < tol2);
        let gap;
        match hit {
            Some(j) => {
                let cache_pnt = cache.get(j);
                let (uv, g) = projector.next_value_of_uv(
                    s,
                    &cache_pnt.second,
                    &p[i],
                    tol_working,
                    tol2,
                    tol_working,
                );
                p2d[i] = uv;
                gap = g;
                saved_point_num = i as i32;
                saved_point = cache_pnt.second;
                if i == 0 {
                    out.is_from_cache = true;
                }
            }
            None => {
                let (uv, g) = projector.value_of_uv(s, &p[i], tol_working, tol2);
                p2d[i] = uv;
                gap = g;
            }
        }
        // `cxx:972-978`: `const double aDist = aProjector.Gap();` then
        // `aCurDist = aDist * aDist` grows the running tolerance. `Gap()` is
        // `mySurf->Gap()` (`ShapeConstruct_ProjectCurveOnSurface.cxx:141`),
        // which `ValueOfUV` closes with `if (myGap <= 0) myGap =
        // P3D.Distance(SurfAdapt.Value(S, T));` (`ShapeAnalysis_Surface.cxx:1502-1513`),
        // so `aDist` is the projection residue and `aTol2` only grows for a
        // poor projection or the conic-extrusion arm (`cxx:1312`), which leaves
        // `myGap = -1.` (`cxx:1249`). `getLine` also leaves the value in
        // `SurfaceProjectorWithCache::Gap()` for the caller's sample loop.
        out.last_gap = gap;
        let cur_dist = gap * gap;
        if tol2 < cur_dist {
            tol2 = cur_dist;
        }
    }
    // `cxx:1020-1021`: `thePoints2d` is written before any failure return.
    out.first2d = p2d[0];
    out.last2d = p2d[3];

    // `cxx:982-1013`: on a periodic surface the second and last-but-one points
    // seed the period-jump fix, which is anchored on the cached point when the
    // first probe came from a cache (`cxx:957-961`).
    let mut is_recompute = false;
    if is_periodic_u || is_periodic_v {
        for i in [1usize, 2] {
            let hit =
                (0..cache.len()).find(|&j| cache.get(j).first.square_distance(&p[i]) < tol2);
            let gap;
            match hit {
                Some(j) => {
                    let cache_pnt = cache.get(j);
                    let (uv, g) = projector.next_value_of_uv(
                        s,
                        &cache_pnt.second,
                        &p[i],
                        tol_working,
                        tol2,
                        tol_working,
                    );
                    p2d[i] = uv;
                    gap = g;
                    saved_point_num = i as i32;
                    saved_point = cache_pnt.second;
                }
                None => {
                    let (uv, g) = projector.value_of_uv(s, &p[i], tol_working, tol2);
                    p2d[i] = uv;
                    gap = g;
                }
            }
            // `cxx:987`: `aProjector.Gap()`, the value `getLine` leaves behind.
            out.last_gap = gap;
            let cur_dist = gap * gap;
            if tol2 < cur_dist {
                tol2 = cur_dist;
            }
        }
        // `cxx:996` / `cxx:1007`: both calls ASSIGN `theIsRecompute`, so on a
        // doubly periodic surface (torus, closed B-spline) the U verdict is
        // overwritten by the V one and only the V result survives; on a surface
        // periodic in U only (cylinder, cone, sphere) the U verdict is the one
        // that reaches the spherical-surface guard below.
        if is_periodic_u {
            is_recompute =
                fix_periodicity_troubles(&mut p2d, 0, s.u_period(), saved_point_num, saved_point.x());
        }
        if is_periodic_v {
            is_recompute =
                fix_periodicity_troubles(&mut p2d, 1, s.v_period(), saved_point_num, saved_point.y());
        }
    }
    out.is_recompute = is_recompute;
    out.probes2d = p2d;

    // `cxx:1015-1018`: a spherical surface with a periodicity jump is given up
    // (the caller's `fixPeriodicityTroubles` comment: `rln S4135 sphere is not
    // considered as V-closed anymore`).
    if is_recompute && s.gp_sphere().is_some() {
        return None;
    }

    tol2 = an_old_tol2;

    let d_par = params[nb - 1] - params[0];
    if d_par.abs() < PCONFUSION {
        return None;
    }

    // `cxx:1032-1034`: the candidate pcurve is the straight 2D chord, linear in
    // the 3D parameter.
    let vec0 = GpVec2d::new(p2d[3].x() - p2d[0].x(), p2d[3].y() - p2d[0].y());
    let vec = GpVec2d::new(vec0.x() / d_par, vec0.y() / d_par);
    let candidate = |t: f64| {
        let dt = t - params[0];
        GpPnt2d::new(p2d[0].x() + vec.x() * dt, p2d[0].y() + vec.y() * dt)
    };

    // `cxx:1036-1058`: on a C1 surface every sample must lie on the normal line
    // of its candidate `(u, v)`.
    let mut is_normal_check = is_cn_1(s);
    if is_normal_check {
        for i in 0..nb {
            let cu = candidate(params[i]);
            let (cur_p, du, dv) = s.d1(cu.x(), cu.y());
            let normal_vec = du.crossed(&dv);
            if normal_vec.square_magnitude() < SQUARE_CONFUSION {
                is_normal_check = false;
                break;
            }
            // `gp_Lin(aCurP, gp_Dir(aNormalVec)).Distance(thePoints(i))`.
            let dist = normal_vec
                .crossed(&GpVec::from_pnts(&cur_p, &points[i]))
                .magnitude()
                / normal_vec.magnitude();
            if dist > tol_working {
                return None;
            }
        }
    }

    // `cxx:1060-1076`: otherwise the surface distance must stay constant.
    if !is_normal_check {
        let first_point_dist = s
            .value(p2d[0].x(), p2d[0].y())
            .square_distance(&points[0]);
        tol2 = tol2.max(tol2 * 2.0 * first_point_dist);
        for i in 1..nb - 1 {
            let cu = candidate(params[i]);
            let cur_p = s.d0(cu.x(), cu.y());
            let dist1 = cur_p.square_distance(&points[i]);
            if (first_point_dist - dist1).abs() > tol2 {
                return None;
            }
        }
    }

    // `cxx:1079-1086`: an exact `Geom2d_Line` when the 2D chord is unit speed.
    let l_length = vec0.magnitude();
    if (l_length - d_par).abs() <= PCONFUSION {
        let dir_l = GpVec2d::new(vec0.x() / l_length, vec0.y() / l_length);
        let dir2d = GpDir2d::from_vec2d(&dir_l).ok()?;
        let loc = GpPnt2d::new(
            p2d[0].x() - params[0] * dir_l.x(),
            p2d[0].y() - params[0] * dir_l.y(),
        );
        return Some(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir2d))));
    }

    // `cxx:1088-1101`: a two-pole degree-1 B-spline for the rest.
    let bs = bspline_from_samples(&[p2d[0], p2d[3]], params[0], params[nb - 1], 1).ok()?;
    Some(Arc::new(bs))
}
