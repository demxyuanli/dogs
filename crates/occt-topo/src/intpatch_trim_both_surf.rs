//! `FUN_NewFirstLast` / `FUN_TrimBothSurf`
//! (`IntPatch_Intersection.cxx:889-1062`).
//!
//! `static` helpers of `IntPatch_Intersection.cxx` used by `ParamParamPerfom`
//! when *both* domains are infinite and the two surfaces do not share a
//! collinear pair of iso lines (`:1744-1747`), and by `GeomParamPerfom`
//! (`:1930-1938`). Each surface's iso line through the mid parameter is
//! classified (`GeomAdaptor_Curve::GetType`) and an unbounded direction is cut
//! at `TV` (`1.e+8` on both call sites), because
//! `IntPatch_PrmPrmIntersection` cannot sample an infinite domain.

use std::sync::Arc;

use occt_geom::Surface;

use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;

use super::trim_inf_surf::{u_trim, v_trim};

/// `FUN_NewFirstLast(ga_ct, Fst, Lst, TrVal, NewFst, NewLst, NeedTr)`
/// (`cxx:891-961`). Returns `NeedTr`; `NewFst` / `NewLst` start as
/// `(Fst, Lst)`.
///
/// Note the `.cxx` braces: in both arms the `else` at `:918` / `:948` binds to
/// the *second* `if` (`Fst < 0. && Lst < 0.`), not to the first, so for a
/// range that is entirely non-negative the first `if`'s clamp is overwritten
/// by the `else` body. The port keeps the same shape instead of folding the
/// two tests into one three-way chain.
pub(crate) fn new_first_last(
    ga_ct: GeomAbsCurveType,
    fst: f64,
    lst: f64,
    tr_val: f64,
    new_fst: &mut f64,
    new_lst: &mut f64,
) -> bool {
    *new_fst = fst;
    *new_lst = lst;
    let mut need_tr = false;
    match ga_ct {
        GeomAbsCurveType::Line | GeomAbsCurveType::Parabola => {
            if (lst - fst).abs() > tr_val {
                if fst >= 0.0 && lst >= 0.0 {
                    *new_fst = fst;
                    *new_lst = if (fst + tr_val) < lst { fst + tr_val } else { lst };
                }
                if fst < 0.0 && lst < 0.0 {
                    *new_lst = lst;
                    *new_fst = if (lst - tr_val) > fst { lst - tr_val } else { fst };
                } else {
                    *new_fst = if fst < -tr_val { -tr_val } else { fst };
                    *new_lst = if lst > tr_val { tr_val } else { lst };
                }
                need_tr = true;
            }
        }
        GeomAbsCurveType::Hyperbola => {
            if (lst - fst).abs() > 10.0 {
                if fst >= 0.0 && lst >= 0.0 {
                    if fst > 4.0 {
                        return false;
                    }
                    *new_fst = fst;
                    *new_lst = if lst > 4.0 { 4.0 } else { lst };
                }
                if fst < 0.0 && lst < 0.0 {
                    if lst < -4.0 {
                        return false;
                    }
                    *new_lst = lst;
                    *new_fst = if fst < -4.0 { -4.0 } else { fst };
                } else {
                    *new_fst = if fst < -4.0 { -4.0 } else { fst };
                    *new_lst = if lst > 4.0 { 4.0 } else { lst };
                }
                need_tr = true;
            }
        }
        _ => {}
    }
    need_tr
}

/// `FUN_TrimBothSurf(S1, T1, S2, T2, TV, NS1, NS2)` (`cxx:966-1062`).
///
/// `NS1` / `NS2` are `None` when the surface is left untrimmed; the caller then
/// keeps passing the original `S1` / `S2` (`NS1 = S1` at `:1009`, and the
/// handles start out equal to the arguments).
pub(crate) fn trim_both_surf(
    s1: &dyn Surface,
    t1: GeomAbsSurfaceType,
    s2: &dyn Surface,
    t2: GeomAbsSurfaceType,
    tv: f64,
) -> (Option<Arc<dyn Surface>>, Option<Arc<dyn Surface>>) {
    let (u1, u2) = s1.u_range();
    let (v1, v2) = s1.v_range();
    let um1 = 0.5 * (u2 + u1);
    let vm1 = 0.5 * (v2 + v1);
    let (u1b, u2b) = s2.u_range();
    let (v1b, v2b) = s2.v_range();
    let um2 = 0.5 * (u2b + u1b);
    let vm2 = 0.5 * (v2b + v1b);

    // `if (T1 != GeomAbs_OffsetSurface) gs1->VIso(VM1) / gs1->UIso(UM1); else`
    // the basis surface's (`cxx:983-996`).
    let (viso_s1, uiso_s1) = if t1 != GeomAbsSurfaceType::OffsetSurface {
        (s1.v_iso_curve(vm1), s1.u_iso_curve(um1))
    } else {
        match s1.offset_basis_surface() {
            Some(b) => (b.v_iso_curve(vm1), b.u_iso_curve(um1)),
            None => (None, None),
        }
    };
    let (viso_s2, uiso_s2) = if t2 != GeomAbsSurfaceType::OffsetSurface {
        (s2.v_iso_curve(vm2), s2.u_iso_curve(um2))
    } else {
        match s2.offset_basis_surface() {
            Some(b) => (b.v_iso_curve(vm2), b.u_iso_curve(um2)),
            None => (None, None),
        }
    };

    let (Some(uiso_s1), Some(uiso_s2), Some(viso_s1), Some(viso_s2)) =
        (uiso_s1, uiso_s2, viso_s1, viso_s2)
    else {
        // `NS1 = S1; NS2 = S2; return;` (`cxx:1009-1012`).
        return (None, None);
    };

    // `GeomAdaptor_Curve gau1(uisoS1); ... GA_U1 = gau1.GetType();`
    // (`cxx:1013-1020`). The V-parameter range is clamped by the type of the
    // U-iso curve, and vice versa (`FUN_NewFirstLast(GA_U1, FirstV, LastV)`
    // at `:1026`).
    let ga_u1 = curve_kind(uiso_s1.as_ref());
    let ga_v1 = curve_kind(viso_s1.as_ref());
    let ga_u2 = curve_kind(uiso_s2.as_ref());
    let ga_v2 = curve_kind(viso_s2.as_ref());

    let (mut v1s1, mut v2s1, mut u1s1, mut u2s1) = (v1, v2, u1, u2);
    let (mut v1s2, mut v2s2, mut u1s2, mut u2s2) = (v1b, v2b, u1b, u2b);
    let trm_v1 = new_first_last(ga_u1, v1, v2, tv, &mut v1s1, &mut v2s1);
    let trm_u1 = new_first_last(ga_v1, u1, u2, tv, &mut u1s1, &mut u2s1);
    let trm_v2 = new_first_last(ga_u2, v1b, v2b, tv, &mut v1s2, &mut v2s2);
    let trm_u2 = new_first_last(ga_v2, u1b, u2b, tv, &mut u1s2, &mut u2s2);

    // `if (TrmV1) NS1 = S1->VTrim(...)`, then
    // `if (TrmU1) { if (TrmV1) TS = NS1; NS1 = TS->UTrim(...); }`
    // (`cxx:1031-1062`). The U trim always runs on the V-trimmed patch when
    // both fire.
    let mut ns1: Option<Arc<dyn Surface>> = None;
    if trm_v1 {
        ns1 = v_trim(s1, v1s1, v2s1);
    }
    if trm_u1 {
        let base: &dyn Surface = match ns1.as_deref() {
            Some(t) => t,
            None => s1,
        };
        ns1 = u_trim(base, u1s1, u2s1).or(ns1);
    }
    let mut ns2: Option<Arc<dyn Surface>> = None;
    if trm_v2 {
        ns2 = v_trim(s2, v1s2, v2s2);
    }
    if trm_u2 {
        let base: &dyn Surface = match ns2.as_deref() {
            Some(t) => t,
            None => s2,
        };
        ns2 = u_trim(base, u1s2, u2s2).or(ns2);
    }
    (ns1, ns2)
}
