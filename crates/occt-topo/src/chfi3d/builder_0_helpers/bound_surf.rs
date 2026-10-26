//! `ChFi3d_BoundSurf` — the surface of a `ChFiDS_SurfData`, trimmed to a box
//! around the two interferences so an intersection can be computed on it.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4453-4505`
//! (OCCT 8.0.0). It uses `ChFi3d_Boite` (`cxx:259-284`), [`ch_fi3d_boite4`].
//!
//! Translation notes:
//! - `occ::handle<GeomAdaptor_Surface>` is
//!   [`SurfaceAdaptor`], whose `of_surface_with_range` is
//!   `GeomAdaptor_Surface::Load(S, U1, U2, V1, V2)`,
//! - `DStr.Surface(Fd1->Surf()).Surface()` is the `Arc<dyn Surface>` of the DS
//!   entry; a null handle there is the OCCT null dereference and answers
//!   `Err`,
//! - the two `IFaCo1 == 0 || IFaArc1 == 0` guards (`cxx:4466-4471`) answer the
//!   unbound surface, exactly like the source.

use std::sync::Arc;

use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom::surface::Surface;

use crate::chfi3d::builder_0_helpers::boite::ch_fi3d_boite4;
use crate::chfi_ds::SurfData;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_BoundSurf(DStr, Fd1, IFaCo1, IFaArc1)` (`cxx:4453-4505`).
///
/// `Err` carries the `DStr.Surface(Fd1->Surf()).Surface()` null-handle
/// dereference OCCT would take, and the missing p-curve on either
/// interference (OCCT null-dereferences `PCurveOnSurf()`).
///
/// OCCT also writes `Du` / `Dv` in every branch (`cxx:4486-4503`) even where
/// that branch never reads them back; those writes are part of the flow and
/// are kept verbatim.
#[allow(unused_assignments)]
pub fn ch_fi3d_bound_surf(
    ds: &DataStructure,
    fd1: &SurfData,
    i_fa_co1: i32,
    i_fa_arc1: i32,
) -> Result<SurfaceAdaptor, String> {
    // `S1.Load(DStr.Surface(Fd1->Surf()).Surface());` (`cxx:4465-4467`).
    let surface: Arc<dyn Surface> = ds
        .surface(fd1.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_BoundSurf: DStr.Surface(Surf()).Surface() is null".to_string())?;

    // `if ((IFaCo1 == 0) || (IFaArc1 == 0)) return HS1;` (`cxx:4469-4472`).
    if i_fa_co1 == 0 || i_fa_arc1 == 0 {
        return Ok(SurfaceAdaptor::of_surface(surface));
    }

    // `const ChFiDS_FaceInterference& FiCo1 = Fd1->Interference(IFaCo1);`
    // (`cxx:4474-4475`).
    let fi_co1 = fd1.interference(i_fa_co1);
    let fi_arc1 = fd1.interference(i_fa_arc1);

    // `UVf1 = FiCo1.PCurveOnSurf()->Value(FiCo1.FirstParameter());` and the
    // three following lines (`cxx:4478-4482`).
    let pc_co1 = fi_co1
        .p_curve_on_surf()
        .ok_or_else(|| "ChFi3d_BoundSurf: FiCo1 has no p-curve on surf".to_string())?;
    let pc_arc1 = fi_arc1
        .p_curve_on_surf()
        .ok_or_else(|| "ChFi3d_BoundSurf: FiArc1 has no p-curve on surf".to_string())?;
    let uv_f1 = pc_co1.value(fi_co1.first_parameter());
    let uv_l1 = pc_co1.value(fi_co1.last_parameter());
    let uv_f2 = pc_arc1.value(fi_arc1.first_parameter());
    let uv_l2 = pc_arc1.value(fi_arc1.last_parameter());

    // `ChFi3d_Boite(UVf1, UVf2, UVl1, UVl2, Du, Dv, mu, Mu, mv, Mv);`
    // (`cxx:4484`).
    let b = ch_fi3d_boite4(&uv_f1, &uv_f2, &uv_l1, &uv_l2);
    let mut du = b.du;
    let mut dv = b.dv;

    // `GeomAbs_SurfaceType styp = S1.GetType();` (`cxx:4485`).
    let mut first_u = b.u_min;
    let mut last_u = b.u_max;
    let mut first_v = b.v_min;
    let mut last_v = b.v_max;
    match SurfaceAdaptor::of_surface(Arc::clone(&surface)).get_type() {
        // `cxx:4486-4489`.
        GeomAbsSurfaceType::Cylinder => {
            let radius = surface
                .gp_cylinder()
                .ok_or_else(|| "ChFi3d_BoundSurf: GeomAbs_Cylinder without a cylinder".to_string())?
                .radius();
            dv = (0.5 * dv).max(4.0 * radius);
            du = 0.0;
            first_v = b.v_min - dv;
            last_v = b.v_max + dv;
        }
        // `cxx:4492-4497`.
        GeomAbsSurfaceType::Torus | GeomAbsSurfaceType::Cone => {
            du = (std::f64::consts::PI - 0.5 * du).min(0.1 * du);
            dv = 0.0;
            first_u = b.u_min - du;
            last_u = b.u_max + du;
        }
        // `cxx:4498-4503`.
        GeomAbsSurfaceType::Plane => {
            du = (0.5 * du).max(4.0 * dv);
            dv = 0.0;
            first_u = b.u_min - du;
            last_u = b.u_max + du;
        }
        _ => {}
    }
    Ok(SurfaceAdaptor::of_surface_with_range(
        surface, first_u, last_u, first_v, last_v,
    ))
}

/// `ChFi3d_BoundSrfImpl(SurfaceAdaptor&, Umin, Umax, Vmin, Vmax,
/// theCheckNaturalBounds)` (`cxx:157-236`), the `BRepAdaptor_Surface`
/// instantiation of the template: `ChFi3d_BaseSurface` is
/// `GeomSurfaceOriginal()` (`cxx:118-122`), `ChFi3d_ApplyBounds` is
/// `Load(GeomSurface, U1, U2, V1, V2, Trsf, TolU, TolV)` (`cxx:141-154`), and
/// the port's `SurfaceAdaptor` already carries the world surface, so both are
/// [`SurfaceAdaptor::of_surface_with_range`].
///
/// `aStepU / aScaleU` is verbatim: OCCT divides the parameter span by the
/// adaptor's *parameter* resolution. The two `aPeriodU > 0.` arms then
/// overwrite the step with `0.1 * (period - span)`.
fn ch_fi3d_bound_srf_impl(
    adaptor: &SurfaceAdaptor,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    check_natural_bounds: bool,
) -> SurfaceAdaptor {
    // `occ::handle<Geom_Surface> aSurface = ChFi3d_BaseSurface(theSurface);`
    // then the `Geom_RectangularTrimmedSurface` unwrap (`cxx:169-176`). A null
    // handle there is the OCCT null dereference; the port answers the adaptor
    // unchanged.
    let Some(base) = adaptor.surface() else {
        return adaptor.clone();
    };
    let base = base.rectangular_trimmed_basis().unwrap_or_else(|| Arc::clone(base));

    // `aSurface->Bounds(aU1, aU2, aV1, aV2);` (`cxx:178-179`).
    let (a_u1, a_u2) = base.u_range();
    let (a_v1, a_v2) = base.v_range();

    // `cxx:180-188`.
    let a_period_u = if base.is_u_periodic() { base.u_period() } else { 0.0 };
    let a_period_v = if base.is_v_periodic() { base.v_period() } else { 0.0 };

    // `cxx:190-197`.
    let mut a_step_u = umax - umin;
    let mut a_step_v = vmax - vmin;
    let a_scale_u = u_resolution(base.as_ref(), 1.0);
    let a_scale_v = v_resolution(base.as_ref(), 1.0);
    let a_step3d_u = a_step_u / a_scale_u;
    let a_step3d_v = a_step_v / a_scale_v;

    // `cxx:198-207`.
    if a_step3d_u > a_step3d_v {
        a_step_v = a_step3d_u * a_scale_v;
    }
    if a_step3d_v > a_step3d_u {
        a_step_u = a_step3d_v * a_scale_u;
    }

    // `cxx:208-215`.
    if a_period_u > 0.0 {
        a_step_u = 0.1 * (a_period_u - (umax - umin));
    }
    if a_period_v > 0.0 {
        a_step_v = 0.1 * (a_period_v - (vmax - vmin));
    }

    // `cxx:217-219`.
    let mut a_uu1 = umin - a_step_u;
    let mut a_uu2 = umax + a_step_u;
    let mut a_vv1 = vmin - a_step_v;
    let mut a_vv2 = vmax + a_step_v;

    // `cxx:220-232`.
    if check_natural_bounds {
        if !base.is_u_periodic() {
            a_uu1 = a_uu1.max(a_u1);
            a_uu2 = a_uu2.min(a_u2);
        }
        if !base.is_v_periodic() {
            a_vv1 = a_vv1.max(a_v1);
            a_vv2 = a_vv2.min(a_v2);
        }
    }

    // `ChFi3d_ApplyBounds(theSurface, aSurface, aUU1, aUU2, aVV1, aVV2);`
    // (`cxx:234`).
    SurfaceAdaptor::of_surface_with_range(base, a_uu1, a_uu2, a_vv1, a_vv2)
}

/// `ChFi3d_BoundFac(BRepAdaptor_Surface&, Umin, Umax, Vmin, Vmax,
/// checknaturalbounds = true)` (`cxx:691-698`, default `hxx:125`):
/// [`ch_fi3d_bound_srf_impl`] with `theCheckNaturalBounds = true`, which is
/// what both `PerformOneCorner` call sites (`cxx:1005`, `cxx:1471` / `:1474`)
/// take.
pub fn ch_fi3d_bound_fac(
    adaptor: &SurfaceAdaptor,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> SurfaceAdaptor {
    ch_fi3d_bound_srf_impl(adaptor, umin, umax, vmin, vmax, true)
}
