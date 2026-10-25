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
