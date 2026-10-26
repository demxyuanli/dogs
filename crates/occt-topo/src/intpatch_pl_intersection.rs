//! `FUN_GetUiso` / `FUN_GetViso` / `FUN_PL_Intersection`
//! (`IntPatch_Intersection.cxx:451-888`).
//!
//! `static` helpers of `IntPatch_Intersection.cxx` for the degenerate
//! "both domains infinite" arm of `ParamParamPerfom` (`:1710-1750`): when each
//! surface has at least one straight iso line and a pair of those isos is
//! parallel, the intersection of two infinite surfaces is a family of lines
//! (`IntPatch_GLine`) instead of something a walker could trace.
//!
//! The `.cxx` builds the isos (`FUN_GetViso` / `FUN_GetUiso`), projects both on
//! the plane through `puvS1` with normal `DV`, intersects the two 2D images
//! (`Geom2dAPI_InterCurveCurve`) and lifts every point back with
//! `ElCLib::To3d`. The port keeps that sequence; the projection step is folded
//! into [`crate::pcurve::project_curve_on_plane`], which is
//! `GeomProjLib::ProjectOnPlane(KeepParam=true)` followed by
//! `GeomProjLib::Curve2d` on a plane (see its doc comment), and the 2D solver
//! is [`occt_geom2d::geom2d_api::intersect_curves`]
//! (`Geom2dAPI_InterCurveCurve`).

use std::f64::consts::TAU;
use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpDir, GpPnt, GpPln, GpVec};
use occt_core::precision::Precision;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom::{Curve, Surface};
use occt_geom2d::geom2d_api::intersect_curves;

use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::pcurve::project_curve_on_plane;

/// `FUN_GetUiso(GS, T, FirstV, LastV, IsVC, IsVP, U, I)`
/// (`cxx:451-553`): the U-iso of `GS` at `U`, trimmed to `[FirstV, LastV]`
/// unless V is periodic over a full turn. `None` is the `.cxx`'s null `I`
/// (`return;` for an offset hyperbola iso outside `|4|`, `cxx:527-539`).
pub(crate) fn get_uiso(
    gs: &dyn Surface,
    t: GeomAbsSurfaceType,
    first_v: f64,
    last_v: f64,
    is_vc: bool,
    is_vp: bool,
    u: f64,
) -> Option<Arc<dyn Curve>> {
    if t != GeomAbsSurfaceType::OffsetSurface {
        let gc = gs.u_iso_curve(u)?;
        if is_vp && first_v == 0.0 && last_v == TAU {
            return Some(gc);
        }
        return Some(Arc::new(GeomTrimmedCurve::new(gc, first_v, last_v)));
    }
    // `else // OffsetSurface` (`cxx:478-552`).
    let gos = gs;
    let bs = gos.offset_basis_surface()?;
    let gcbs = bs.u_iso_curve(u)?;
    let gact = curve_kind(gcbs.as_ref());
    if is_vp
        || is_vc
        || gact == GeomAbsCurveType::BSplineCurve
        || gact == GeomAbsCurveType::BezierCurve
        || (last_v - first_v).abs() < 1.0e5
    {
        let gc = gos.u_iso_curve(u)?;
        if is_vp && first_v == 0.0 && last_v == TAU {
            return Some(gc);
        }
        return Some(Arc::new(GeomTrimmedCurve::new(gc, first_v, last_v)));
    }
    // `else // Offset Line, Parab, Hyperb` (`cxx:499-552`).
    let (vm_tr, v_m_tr) = if gact != GeomAbsCurveType::Hyperbola {
        if first_v >= 0.0 && last_v >= 0.0 {
            (first_v, if (last_v - first_v) > 1.0e4 { first_v + 1.0e4 } else { last_v })
        } else if first_v < 0.0 && last_v < 0.0 {
            (if (first_v - last_v) < -1.0e4 { last_v - 1.0e4 } else { first_v }, last_v)
        } else {
            (
                if first_v < -1.0e4 { -1.0e4 } else { first_v },
                if last_v > 1.0e4 { 1.0e4 } else { last_v },
            )
        }
    } else {
        // `else // Hyperbola` (`cxx:520-547`), clamped to `|4|`.
        if first_v >= 0.0 && last_v >= 0.0 {
            if first_v > 4.0 {
                return None;
            }
            (first_v, if last_v > 4.0 { 4.0 } else { last_v })
        } else if first_v < 0.0 && last_v < 0.0 {
            if last_v < -4.0 {
                return None;
            }
            (if first_v < -4.0 { -4.0 } else { first_v }, last_v)
        } else {
            (
                if first_v < -4.0 { -4.0 } else { first_v },
                if last_v > 4.0 { 4.0 } else { last_v },
            )
        }
    };
    // `new Geom_RectangularTrimmedSurface(gos, VmTr, VMTr, true); I = rts->UIso(U);`
    // (`cxx:548-549`). The `true` is the `.cxx`'s `UTrim` flag even though the
    // bounds come from the V range; it is kept as written.
    let rts =
        GeomRectangularTrimmedSurface::one_param(Arc::from(gos.clone_dyn()), vm_tr, v_m_tr, true)
            .ok()?;
    rts.u_iso_curve(u)
}

/// `FUN_GetViso(GS, T, FirstU, LastU, IsUC, IsUP, V, I)` (`cxx:557-657`): the
/// V mirror of [`get_uiso`], clamped to `|4|` in U for an offset hyperbola
/// (`cxx:626-646`).
pub(crate) fn get_viso(
    gs: &dyn Surface,
    t: GeomAbsSurfaceType,
    first_u: f64,
    last_u: f64,
    is_uc: bool,
    is_up: bool,
    v: f64,
) -> Option<Arc<dyn Curve>> {
    if t != GeomAbsSurfaceType::OffsetSurface {
        let gc = gs.v_iso_curve(v)?;
        if is_up && first_u == 0.0 && last_u == TAU {
            return Some(gc);
        }
        return Some(Arc::new(GeomTrimmedCurve::new(gc, first_u, last_u)));
    }
    let gos = gs;
    let bs = gos.offset_basis_surface()?;
    let gcbs = bs.v_iso_curve(v)?;
    let gact = curve_kind(gcbs.as_ref());
    if is_up
        || is_uc
        || gact == GeomAbsCurveType::BSplineCurve
        || gact == GeomAbsCurveType::BezierCurve
        || (last_u - first_u).abs() < 1.0e5
    {
        let gc = gos.v_iso_curve(v)?;
        if is_up && first_u == 0.0 && last_u == TAU {
            return Some(gc);
        }
        return Some(Arc::new(GeomTrimmedCurve::new(gc, first_u, last_u)));
    }
    let (um_tr, u_m_tr) = if gact != GeomAbsCurveType::Hyperbola {
        if first_u >= 0.0 && last_u >= 0.0 {
            (first_u, if (last_u - first_u) > 1.0e4 { first_u + 1.0e4 } else { last_u })
        } else if first_u < 0.0 && last_u < 0.0 {
            (if (first_u - last_u) < -1.0e4 { last_u - 1.0e4 } else { first_u }, last_u)
        } else {
            (
                if first_u < -1.0e4 { -1.0e4 } else { first_u },
                if last_u > 1.0e4 { 1.0e4 } else { last_u },
            )
        }
    } else {
        if first_u >= 0.0 && last_u >= 0.0 {
            if first_u > 4.0 {
                return None;
            }
            (first_u, if last_u > 4.0 { 4.0 } else { last_u })
        } else if first_u < 0.0 && last_u < 0.0 {
            if last_u < -4.0 {
                return None;
            }
            (if first_u < -4.0 { -4.0 } else { first_u }, last_u)
        } else {
            (
                if first_u < -4.0 { -4.0 } else { first_u },
                if last_u > 4.0 { 4.0 } else { last_u },
            )
        }
    };
    // `new Geom_RectangularTrimmedSurface(gos, UmTr, UMTr, true); I = rts->VIso(V);`
    // (`cxx:650-652`).
    let rts =
        GeomRectangularTrimmedSurface::one_param(Arc::from(gos.clone_dyn()), um_tr, u_m_tr, true)
            .ok()?;
    rts.v_iso_curve(v)
}

/// Result of `FUN_PL_Intersection`: `IsOk`, the point sequence `SP` and the
/// direction `DV` the `.cxx` builds each `IntPatch_GLine` from.
pub(crate) struct PlIntersection {
    pub is_ok: bool,
    pub points: Vec<GpPnt>,
    pub dir: GpVec,
}

/// `FUN_PL_Intersection(S1, T1, S2, T2, IsOk, SP, DV)` (`cxx:660-888`).
///
/// `IsOk` is `true` once two parallel straight isos were found, even when the
/// projected curves do not meet (`SP` empty, `:860-869`); the caller turns an
/// empty `SP` into `done = false` (`:1734-1738`).
pub(crate) fn pl_intersection(
    s1: &dyn Surface,
    t1: GeomAbsSurfaceType,
    s2: &dyn Surface,
    t2: GeomAbsSurfaceType,
) -> PlIntersection {
    // `IsOk = false; DV = gp_Vec(0., 0., 1.);` (`cxx:668-670`).
    let mut dir = GpVec::new(0.0, 0.0, 1.0);
    let mut is_ok = false;
    let mut points = Vec::new();

    let (u1f, u1l) = s1.u_range();
    let (v1f, v1l) = s1.v_range();
    let (u2f, u2l) = s2.u_range();
    let (v2f, v2l) = s2.v_range();
    // `MS1[0] = 0.5 * (LastU + FirstU)` etc. (`cxx:679-683`).
    let ms1 = [0.5 * (u1l + u1f), 0.5 * (v1l + v1f)];
    let ms2 = [0.5 * (u2l + u2f), 0.5 * (v2l + v2f)];

    // 1. Check: both surfaces have U(V)isos - lines (`cxx:685-762`).
    let mut iso_s1_is_line = [false, false];
    let mut iso_s2_is_line = [false, false];
    if t1 == GeomAbsSurfaceType::SurfaceOfExtrusion {
        iso_s1_is_line[0] = true;
    } else if !s1.is_v_periodic() && !s1.is_v_closed() {
        let c1 = if t1 != GeomAbsSurfaceType::OffsetSurface {
            s1.u_iso_curve(ms1[0])
        } else {
            s1.offset_basis_surface().and_then(|b| b.u_iso_curve(ms1[0]))
        };
        if let Some(c1) = c1 {
            if curve_kind(c1.as_ref()) == GeomAbsCurveType::Line {
                iso_s1_is_line[0] = true;
            }
        }
    }
    if !s1.is_u_periodic() && !s1.is_u_closed() {
        let c1 = if t1 != GeomAbsSurfaceType::OffsetSurface {
            s1.v_iso_curve(ms1[1])
        } else {
            s1.offset_basis_surface().and_then(|b| b.v_iso_curve(ms1[1]))
        };
        if let Some(c1) = c1 {
            if curve_kind(c1.as_ref()) == GeomAbsCurveType::Line {
                iso_s1_is_line[1] = true;
            }
        }
    }
    if t2 == GeomAbsSurfaceType::SurfaceOfExtrusion {
        iso_s2_is_line[0] = true;
    } else if !s2.is_v_periodic() && !s2.is_v_closed() {
        let c2 = if t2 != GeomAbsSurfaceType::OffsetSurface {
            s2.u_iso_curve(ms2[0])
        } else {
            s2.offset_basis_surface().and_then(|b| b.u_iso_curve(ms2[0]))
        };
        if let Some(c2) = c2 {
            if curve_kind(c2.as_ref()) == GeomAbsCurveType::Line {
                iso_s2_is_line[0] = true;
            }
        }
    }
    if !s2.is_u_periodic() && !s2.is_u_closed() {
        let c2 = if t2 != GeomAbsSurfaceType::OffsetSurface {
            s2.v_iso_curve(ms2[1])
        } else {
            s2.offset_basis_surface().and_then(|b| b.v_iso_curve(ms2[1]))
        };
        if let Some(c2) = c2 {
            if curve_kind(c2.as_ref()) == GeomAbsCurveType::Line {
                iso_s2_is_line[1] = true;
            }
        }
    }
    let is_both_lines =
        (iso_s1_is_line[0] || iso_s1_is_line[1]) && (iso_s2_is_line[0] || iso_s2_is_line[1]);
    if !is_both_lines {
        return PlIntersection { is_ok, points, dir };
    }

    // 2. Check: Uiso lines of both surfaces are collinear (`cxx:769-862`).
    let (puv_s1, der_s1_0, der_s1_1) = s1.d1(ms1[0], ms1[1]);
    let (_, der_s2_0, der_s2_1) = s2.d1(ms2[0], ms2[1]);
    let der_s1 = [der_s1_0, der_s1_1];
    let der_s2 = [der_s2_0, der_s2_1];
    let ang = Precision::ANGULAR;

    let (iso, c1, c2): (usize, Option<Arc<dyn Curve>>, Option<Arc<dyn Curve>>) = if iso_s1_is_line[0]
        && iso_s2_is_line[0]
        && der_s1[1].is_parallel_ang(&der_s2[1], ang)
    {
        (
            1,
            get_viso(s1, t1, u1f, u1l, s1.is_u_closed(), s1.is_u_periodic(), ms1[1]),
            get_viso(s2, t2, u2f, u2l, s2.is_u_closed(), s2.is_u_periodic(), ms2[1]),
        )
    } else if iso_s1_is_line[0] && iso_s2_is_line[1] && der_s1[1].is_parallel_ang(&der_s2[0], ang) {
        (
            1,
            get_viso(s1, t1, u1f, u1l, s1.is_u_closed(), s1.is_u_periodic(), ms1[1]),
            get_uiso(s2, t2, v2f, v2l, s2.is_v_closed(), s2.is_v_periodic(), ms2[0]),
        )
    } else if iso_s1_is_line[1] && iso_s2_is_line[0] && der_s1[0].is_parallel_ang(&der_s2[1], ang) {
        (
            0,
            get_uiso(s1, t1, v1f, v1l, s1.is_v_closed(), s1.is_v_periodic(), ms1[0]),
            get_viso(s2, t2, u2f, u2l, s2.is_u_closed(), s2.is_u_periodic(), ms2[1]),
        )
    } else if iso_s1_is_line[1] && iso_s2_is_line[1] && der_s1[0].is_parallel_ang(&der_s2[0], ang) {
        (
            0,
            get_uiso(s1, t1, v1f, v1l, s1.is_v_closed(), s1.is_v_periodic(), ms1[0]),
            get_uiso(s2, t2, v2f, v2l, s2.is_v_closed(), s2.is_v_periodic(), ms2[0]),
        )
    } else {
        // `IsOk = false; return;` (`cxx:857-861`).
        return PlIntersection { is_ok, points, dir };
    };
    is_ok = true;

    // 3. Make intersections of V(U)isos (`cxx:863-887`).
    let (Some(c1), Some(c2)) = (c1, c2) else {
        return PlIntersection { is_ok, points, dir };
    };
    dir = der_s1[iso];
    let Ok(nrm) = GpDir::from_vec(&dir) else {
        return PlIntersection { is_ok, points, dir };
    };
    // `gp_Ax2(puvS1, gp_Dir(DV))` for the lift and `gp_Pln(puvS1, gp_Dir(DV))`
    // for the projection are the same frame: the plane axis is built from the
    // `gp_Ax2`'s X direction so that `ElCLib::To3d` inverts `Curve2d` exactly.
    let ax3 = GpAx3::from_ax1(&GpAx1::new(puv_s1, nrm));
    let pln = GpPln::new(ax3.clone());
    let Ok(ax2) = GpAx2::new(puv_s1, nrm, *ax3.x_direction()) else {
        return PlIntersection { is_ok, points, dir };
    };
    let (Some(c1_prj), Some(c2_prj)) = (
        project_curve_on_plane(c1.as_ref(), &pln),
        project_curve_on_plane(c2.as_ref(), &pln),
    ) else {
        return PlIntersection { is_ok, points, dir };
    };
    // `Geom2dAPI_InterCurveCurve ICC(C1Prj2d, C2Prj2d, 1.0e-7)` (`cxx:877`).
    for hit in intersect_curves(c1_prj.as_ref(), c2_prj.as_ref(), 1.0e-7) {
        // `ElCLib::To3d(gp_Ax2(puvS1, gp_Dir(DV)), P)` (`cxx:883`).
        points.push(clib::to_3d_pnt(&ax2, &hit.point));
    }
    PlIntersection { is_ok, points, dir }
}
