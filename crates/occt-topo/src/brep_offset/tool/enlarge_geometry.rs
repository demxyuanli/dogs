//! File-local `EnlargeGeometry` of `BRepOffset_Tool`
//! (`BRepOffset_Tool.cxx:2923-3222`), the geometry half of `EnLargeFace`
//! (`:3437-3634`).
//!
//! OCCT takes the surface by `handle&` and replaces it in place; the port takes
//! `&mut Arc<dyn Surface>` and assigns the new surface, which is the same
//! contract (`GeomLib::ExtendSurfByLength` in `crates/occt-geom/src/
//! geom_lib_extend_surf.rs` also answers a new handle instead of mutating).
//!
//! UNPORTED arms:
//! * `cxx:3109` `Geom_BezierSurface` shares the b-spline branch. A Bezier
//!   surface has no `ExtendSurfByLength` support in
//!   `geom_lib_extend_surf` (its `bspline_target` only converts through
//!   `osculating_bspline`, which OCCT's `GeomConvert::SurfaceToBSplineSurface`
//!   equivalent does not provide for a Bezier), so the extension comes back
//!   `None` and the surface is left untrimmed-extended; the branch still
//!   reports `SurfaceChange = true` and the extended bounds, as OCCT does.
//! * `cxx:2988-3009` `Geom_OffsetSurface::SetBasisSurface` mutates the offset
//!   handle in place; the port rebuilds `GeomOffsetSurface::new(basis,
//!   distance)`. That changes the surface pointer, so pcurves keyed by the old
//!   handle become unreachable for the *input* face once `EnLargeFace`
//!   re-registers the surface; the face built by `EnLargeFace` is unaffected.

use std::sync::Arc;

use occt_core::precision::Precision;
use occt_geom::gcpnts::curve_length;
use occt_geom::geom_lib_extend_surf::extend_surf_by_length;
use occt_geom::{GeomOffsetSurface, GeomRectangularTrimmedSurface, Surface};

use super::BRepOffsetTool;

/// `TolApex` of `EnlargeGeometry` (`cxx:2944`), also used by the `MakeFace`
/// degeneracy test (`cxx:2695`).
pub(super) const TOL_APEX: f64 = 1.0e-5;

/// The four UV bounds `EnlargeGeometry` extends and answers
/// (`double& U1, U2, V1, V2`, `cxx:2925-2928`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct EnlargeBounds {
    pub u1: f64,
    pub u2: f64,
    pub v1: f64,
    pub v2: f64,
}

/// `EnlargeGeometry(S, U1, U2, V1, V2, IsV1degen, IsV2degen, uf1, uf2, vf1,
/// vf2, coeff, theGlobalEnlargeU, theGlobalEnlargeVfirst,
/// theGlobalEnlargeVlast, theLenBeforeUfirst, theLenAfterUlast,
/// theLenBeforeVfirst, theLenAfterVlast)` (`cxx:2923-3222`).
///
/// `(uf1, uf2)` / `(vf1, vf2)` are the face's own UV bounds; the `len_*`
/// arguments are `-1` when the caller did not fix a length (`cxx:3037-3039`).
/// Returns `SurfaceChange`.
pub(super) fn enlarge_geometry(
    s: &mut Arc<dyn Surface>,
    bounds: &mut EnlargeBounds,
    is_v1_degen: &mut bool,
    is_v2_degen: &mut bool,
    uf1: f64,
    uf2: f64,
    vf1: f64,
    vf2: f64,
    coeff: f64,
    global_enlarge_u: bool,
    global_enlarge_vfirst: bool,
    global_enlarge_vlast: bool,
    len_before_ufirst: f64,
    len_after_ulast: f64,
    len_before_vfirst: f64,
    len_after_vlast: f64,
) -> bool {
    // `cxx:2945`.
    let mut surface_change = false;

    if let Some(basis) = s.rectangular_trimmed_basis() {
        // `cxx:2947-2986`.
        let mut bs = basis;
        enlarge_geometry(
            &mut bs,
            bounds,
            is_v1_degen,
            is_v2_degen,
            uf1,
            uf2,
            vf1,
            vf2,
            coeff,
            global_enlarge_u,
            global_enlarge_vfirst,
            global_enlarge_vlast,
            len_before_ufirst,
            len_after_ulast,
            len_before_vfirst,
            len_after_vlast,
        );
        if !global_enlarge_vfirst {
            bounds.v1 = vf1;
        }
        if !global_enlarge_vlast {
            bounds.v2 = vf2;
        }
        if !global_enlarge_vfirst || !global_enlarge_vlast {
            // `S = new Geom_RectangularTrimmedSurface(BS, U1, U2, V1, V2)`
            // (`cxx:2980`). The constructor raises on a degenerate range; the
            // port keeps the basis in that case instead of faulting.
            *s = match GeomRectangularTrimmedSurface::uv(
                bs.clone(),
                bounds.u1,
                bounds.u2,
                bounds.v1,
                bounds.v2,
            ) {
                Ok(trimmed) => Arc::new(trimmed),
                Err(_) => bs,
            };
        } else {
            *s = bs;
        }
        surface_change = true;
    } else if s.is_offset_surface() {
        // `cxx:2988-3009`: recurse into the basis, then `SetBasisSurface`.
        if let (Some(basis), Some(distance)) = (s.offset_basis_surface(), s.offset_distance()) {
            let mut bs = basis;
            surface_change = enlarge_geometry(
                &mut bs,
                bounds,
                is_v1_degen,
                is_v2_degen,
                uf1,
                uf2,
                vf1,
                vf2,
                coeff,
                global_enlarge_u,
                global_enlarge_vfirst,
                global_enlarge_vlast,
                len_before_ufirst,
                len_after_ulast,
                len_before_vfirst,
                len_after_vlast,
            );
            *s = Arc::new(GeomOffsetSurface::new(bs, distance));
        }
    } else if s.is_surface_of_linear_extrusion() || s.is_surface_of_revolution() {
        // `cxx:3010-3108`.
        let mut enlarge_u = global_enlarge_u;
        let mut enlarge_v = true;
        let mut enlarge_u_first = enlarge_u;
        let mut enlarge_u_last = enlarge_u;
        let mut enlarge_v_first = global_enlarge_vfirst;
        let mut enlarge_v_last = global_enlarge_vlast;

        let (mut u1, mut u2) = s.u_range();
        let (mut v1, mut v2) = s.v_range();

        let mut du_first = 0.0;
        let mut du_last = 0.0;
        let mut dv_first = 0.0;
        let mut dv_last = 0.0;

        // `cxx:3022-3049`.
        if Precision::is_infinite(u1) || Precision::is_infinite(u2) {
            du_first = uf2 - uf1;
            du_last = du_first;
            u1 = uf1 - du_first;
            u2 = uf2 + du_last;
            enlarge_u = false;
        } else if s.is_u_closed() {
            enlarge_u = false;
        } else if let Some(viso) = s.v_iso_curve(vf1) {
            let du_default = curve_length(viso.as_ref()) * coeff;
            du_first = if len_before_ufirst == -1.0 {
                du_default
            } else {
                len_before_ufirst
            };
            du_last = if len_after_ulast == -1.0 {
                du_default
            } else {
                len_after_ulast
            };
            if let Some(uiso1) = s.u_iso_curve(uf1) {
                if BRepOffsetTool::gabarit(uiso1.as_ref()) <= TOL_APEX {
                    enlarge_u_first = false;
                }
            }
            if let Some(uiso2) = s.u_iso_curve(uf2) {
                if BRepOffsetTool::gabarit(uiso2.as_ref()) <= TOL_APEX {
                    enlarge_u_last = false;
                }
            }
        } else {
            // Port guard: `GeomAdaptor_Curve(viso)` on the null `VIso` of
            // `cxx:3036` faults; without an iso curve the direction cannot be
            // measured, so it is not extended.
            enlarge_u = false;
        }

        // `cxx:3051-3081`.
        if Precision::is_infinite(v1) || Precision::is_infinite(v2) {
            dv_first = vf2 - vf1;
            dv_last = dv_first;
            v1 = vf1 - dv_first;
            v2 = vf2 + dv_last;
            enlarge_v = false;
        } else if s.is_v_closed() {
            enlarge_v = false;
        } else if let Some(uiso) = s.u_iso_curve(uf1) {
            let dv_default = curve_length(uiso.as_ref()) * coeff;
            dv_first = if len_before_vfirst == -1.0 {
                dv_default
            } else {
                len_before_vfirst
            };
            dv_last = if len_after_vlast == -1.0 {
                dv_default
            } else {
                len_after_vlast
            };
            if let Some(viso1) = s.v_iso_curve(vf1) {
                if BRepOffsetTool::gabarit(viso1.as_ref()) <= TOL_APEX {
                    enlarge_v_first = false;
                    *is_v1_degen = true;
                }
            }
            if let Some(viso2) = s.v_iso_curve(vf2) {
                if BRepOffsetTool::gabarit(viso2.as_ref()) <= TOL_APEX {
                    enlarge_v_last = false;
                    *is_v2_degen = true;
                }
            }
        } else {
            // Port guard, the V mirror of the one above (`cxx:3063`).
            enlarge_v = false;
        }

        // `cxx:3082-3107`.
        let mut a_surf: Arc<dyn Surface> =
            match GeomRectangularTrimmedSurface::uv(s.clone(), u1, u2, v1, v2) {
                Ok(trimmed) => Arc::new(trimmed),
                // `new Geom_RectangularTrimmedSurface` raises on a degenerate
                // range; keep the untrimmed surface instead of faulting.
                Err(_) => s.clone(),
            };
        if enlarge_u {
            if enlarge_u_first && du_first != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, du_first, 1, true, false) {
                    a_surf = extended;
                }
            }
            if enlarge_u_last && du_last != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, du_last, 1, true, true) {
                    a_surf = extended;
                }
            }
        }
        if enlarge_v {
            if enlarge_v_first && dv_first != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, dv_first, 1, false, false) {
                    a_surf = extended;
                }
            }
            if enlarge_v_last && dv_last != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, dv_last, 1, false, true) {
                    a_surf = extended;
                }
            }
        }
        *s = a_surf;
        let (nu1, nu2) = s.u_range();
        let (nv1, nv2) = s.v_range();
        bounds.u1 = nu1;
        bounds.u2 = nu2;
        bounds.v1 = nv1;
        bounds.v2 = nv2;
        surface_change = true;
    } else if s.is_bezier_surface() || s.is_bspline_surface() {
        // `cxx:3109-3210`.
        let mut enlarge_u = global_enlarge_u;
        let mut enlarge_v = true;
        let mut enlarge_u_first = enlarge_u;
        let mut enlarge_u_last = enlarge_u;
        let mut enlarge_v_first = global_enlarge_vfirst;
        let mut enlarge_v_last = global_enlarge_vlast;
        if s.is_u_closed() {
            enlarge_u = false;
        }
        if s.is_v_closed() {
            enlarge_v = false;
        }

        let duf = uf2 - uf1;
        let dvf = vf2 - vf1;
        let (u1, u2) = s.u_range();
        let (v1, v2) = s.v_range();

        let mut du_first = 0.0;
        let mut du_last = 0.0;
        let mut dv_first = 0.0;
        let mut dv_last = 0.0;

        let uiso1 = s.u_iso_curve(u1);
        let uiso2 = s.u_iso_curve(u2);
        let viso1 = s.v_iso_curve(v1);
        let viso2 = s.v_iso_curve(v2);
        // OCCT dereferences every `UIso` / `VIso` (`cxx:3133-3139`). A null
        // handle is treated as a degenerate (point-sized) iso, which disables
        // the same directions OCCT's `<= TolApex` tests would.
        let gabarit_uiso1 = uiso1
            .as_ref()
            .map(|c| BRepOffsetTool::gabarit(c.as_ref()))
            .unwrap_or(0.0);
        let gabarit_uiso2 = uiso2
            .as_ref()
            .map(|c| BRepOffsetTool::gabarit(c.as_ref()))
            .unwrap_or(0.0);
        let gabarit_viso1 = viso1
            .as_ref()
            .map(|c| BRepOffsetTool::gabarit(c.as_ref()))
            .unwrap_or(0.0);
        let gabarit_viso2 = viso2
            .as_ref()
            .map(|c| BRepOffsetTool::gabarit(c.as_ref()))
            .unwrap_or(0.0);

        // `cxx:3140-3148`.
        if gabarit_viso1 <= TOL_APEX || gabarit_viso2 <= TOL_APEX {
            enlarge_u = false;
        }
        if gabarit_uiso1 <= TOL_APEX || gabarit_uiso2 <= TOL_APEX {
            enlarge_v = false;
        }

        // `cxx:3149-3180`.
        if enlarge_u {
            if let Some(viso1) = viso1.as_ref() {
                let du_default = curve_length(viso1.as_ref()) * coeff;
                du_first = if len_before_ufirst == -1.0 {
                    du_default
                } else {
                    len_before_ufirst
                };
                du_last = if len_after_ulast == -1.0 {
                    du_default
                } else {
                    len_after_ulast
                };
            }
            if gabarit_uiso1 <= TOL_APEX {
                enlarge_u_first = false;
            }
            if gabarit_uiso2 <= TOL_APEX {
                enlarge_u_last = false;
            }
        }
        if enlarge_v {
            if let Some(uiso1) = uiso1.as_ref() {
                let dv_default = curve_length(uiso1.as_ref()) * coeff;
                dv_first = if len_before_vfirst == -1.0 {
                    dv_default
                } else {
                    len_before_vfirst
                };
                dv_last = if len_after_vlast == -1.0 {
                    dv_default
                } else {
                    len_after_vlast
                };
            }
            if gabarit_viso1 <= TOL_APEX {
                enlarge_v_first = false;
                *is_v1_degen = true;
            }
            if gabarit_viso2 <= TOL_APEX {
                enlarge_v_last = false;
                *is_v2_degen = true;
            }
        }

        // `cxx:3182-3208`.
        let mut a_surf = s.clone();
        if enlarge_u {
            if enlarge_u_first && uf1 - u1 < duf && du_first != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, du_first, 1, true, false) {
                    a_surf = extended;
                }
            }
            if enlarge_u_last && u2 - uf2 < duf && du_last != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, du_last, 1, true, true) {
                    a_surf = extended;
                }
            }
        }
        if enlarge_v {
            if enlarge_v_first && vf1 - v1 < dvf && dv_first != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, dv_first, 1, false, false) {
                    a_surf = extended;
                }
            }
            if enlarge_v_last && v2 - vf2 < dvf && dv_last != 0.0 {
                if let Some(extended) = extend_surf_by_length(&a_surf, dv_last, 1, false, true) {
                    a_surf = extended;
                }
            }
        }
        *s = a_surf;
        let (nu1, nu2) = s.u_range();
        let (nv1, nv2) = s.v_range();
        bounds.u1 = nu1;
        bounds.u2 = nu2;
        bounds.v1 = nv1;
        bounds.v2 = nv2;
        surface_change = true;
    } else {
        // `cxx:3211-3220`: no extension beyond the surface's own bounds.
        let (uu1, uu2) = s.u_range();
        let (vv1, vv2) = s.v_range();
        bounds.u1 = bounds.u1.max(uu1);
        bounds.v1 = bounds.v1.max(vv1);
        bounds.u2 = bounds.u2.min(uu2);
        bounds.v2 = bounds.v2.min(vv2);
    }
    surface_change
}
