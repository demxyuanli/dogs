//! `IntPatch_SpecialPoints` — cone apex / sphere pole / UV seam.

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::geom_int::surface_parameters;
use crate::int_tools_wline::PntOn2S;

/// Result of inserting a singular / seam point into a walking line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecPntType {
    None,
    Pole,
    SeamU,
    SeamUv,
}

/// `IntPatch_SpecialPoints::AddSingularPole`.
pub fn add_singular_pole(
    s_quad: &dyn Surface,
    s_other: &dyn Surface,
    p_ref: &PntOn2S,
    pnt: &GpPnt,
    reversed: bool,
) -> Option<PntOn2S> {
    match classify_surface(s_quad) {
        SurfaceKind::Cone => {
            let apex = s_quad.d0(0.0, 0.0);
            if pnt.square_distance(&apex) > 1.0e-10 {
                return None;
            }
            let (uq, vq) = surface_parameters(s_quad, &apex)?;
            let (uo, vo) = surface_parameters(s_other, &apex)
                .unwrap_or((p_ref.u2, p_ref.v2));
            Some(pack(apex, uq, vq, uo, vo, reversed))
        }
        SurfaceKind::Sphere => {
            let (_u, v) = surface_parameters(s_quad, pnt)?;
            let pole = std::f64::consts::FRAC_PI_2;
            if (v - pole).abs() > 1.0e-5 && (v + pole).abs() > 1.0e-5 {
                return None;
            }
            let (uo, vo) = surface_parameters(s_other, pnt)
                .unwrap_or((p_ref.u2, p_ref.v2));
            let uref = if reversed { p_ref.u2 } else { p_ref.u1 };
            Some(pack(*pnt, uref, v, uo, vo, reversed))
        }
        _ => None,
    }
}

/// `IntPatch_SpecialPoints::AddCrossUVIsoPoint`.
pub fn add_cross_uv_iso_point(
    s_quad: &dyn Surface,
    s_other: &dyn Surface,
    p_ref: &PntOn2S,
    tol3d: f64,
    reversed: bool,
) -> Option<PntOn2S> {
    if classify_surface(s_quad) != SurfaceKind::Torus {
        return None;
    }
    let (u, v) = surface_parameters(s_quad, &p_ref.p)?;
    let (u0, u1) = s_quad.u_range();
    let (v0, v1) = s_quad.v_range();
    let near_u = near_bound(u, u0, u1);
    let near_v = near_bound(v, v0, v1);
    if !near_u || !near_v {
        return None;
    }
    let p = s_quad.d0(u0, v0);
    if p.distance(&p_ref.p) > tol3d.max(CONFUSION) * 100.0 {
        return None;
    }
    let (uo, vo) = surface_parameters(s_other, &p).unwrap_or((p_ref.u2, p_ref.v2));
    Some(pack(p, u0, v0, uo, vo, reversed))
}

fn near_bound(x: f64, a: f64, b: f64) -> bool {
    (x - a).abs() <= PCONFUSION * 10.0 || (x - b).abs() <= PCONFUSION * 10.0
}

fn pack(p: GpPnt, uq: f64, vq: f64, uo: f64, vo: f64, reversed: bool) -> PntOn2S {
    if reversed {
        PntOn2S {
            p,
            u1: uo,
            v1: vo,
            u2: uq,
            v2: vq,
        }
    } else {
        PntOn2S {
            p,
            u1: uq,
            v1: vq,
            u2: uo,
            v2: vo,
        }
    }
}
