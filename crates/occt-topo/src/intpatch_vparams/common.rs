//! Shared statics of `IntPatch_WLine::ComputeVertexParameters` and
//! `IntPatch_RLine::ComputeVertexParameters`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_WLine.cxx:88-321`.

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;

use crate::int_tools_wline::{PatchPoint, PntOn2S};

/// `RecadreMemePeriode(u1, v1, u2, v2, ...)` (`IntPatch_WLine.cxx:88-145`):
/// shift the parametric coordinates by whole periods so that each one stays
/// within `0.8 * period` of the reference value.
#[allow(clippy::too_many_arguments)]
pub fn recadre_meme_periode(
    u1: &mut f64,
    v1: &mut f64,
    u2: &mut f64,
    v2: &mut f64,
    anu1: f64,
    anv1: f64,
    anu2: f64,
    anv2: f64,
    up1: f64,
    vp1: f64,
    up2: f64,
    vp2: f64,
) {
    if up1 != 0.0 {
        while anu1 - *u1 > 0.8 * up1 {
            *u1 += up1;
        }
        while *u1 - anu1 > 0.8 * up1 {
            *u1 -= up1;
        }
    }
    if up2 != 0.0 {
        while anu2 - *u2 > 0.8 * up2 {
            *u2 += up2;
        }
        while *u2 - anu2 > 0.8 * up2 {
            *u2 -= up2;
        }
    }
    if vp1 != 0.0 {
        while anv1 - *v1 > 0.8 * vp1 {
            *v1 += vp1;
        }
        while *v1 - anv1 > 0.8 * vp1 {
            *v1 -= vp1;
        }
    }
    if vp2 != 0.0 {
        while anv2 - *v2 > 0.8 * vp2 {
            *v2 += vp2;
        }
        while *v2 - anv2 > 0.8 * vp2 {
            *v2 -= vp2;
        }
    }
}

/// `RecadreMemePeriode(POn2S, RefPOn2S, ...)` (`IntPatch_WLine.cxx:147-159`).
pub fn recadre_pnt_meme_periode(
    pon2s: &mut PntOn2S,
    ref_pon2s: &PntOn2S,
    up1: f64,
    vp1: f64,
    up2: f64,
    vp2: f64,
) {
    let (mut u1, mut v1, mut u2, mut v2) = (pon2s.u1, pon2s.v1, pon2s.u2, pon2s.v2);
    recadre_meme_periode(
        &mut u1,
        &mut v1,
        &mut u2,
        &mut v2,
        ref_pon2s.u1,
        ref_pon2s.v1,
        ref_pon2s.u2,
        ref_pon2s.v2,
        up1,
        vp1,
        up2,
        vp2,
    );
    pon2s.u1 = u1;
    pon2s.v1 = v1;
    pon2s.u2 = u2;
    pon2s.v2 = v2;
}

/// `CompareVertexAndPoint(V, P, Tol)` (`IntPatch_WLine.cxx:161-166`).
pub fn compare_vertex_and_point(v: &GpPnt, p: &GpPnt, tol: f64) -> bool {
    v.square_distance(p) <= tol * tol
}

/// `SameVtxRst(vtx1, vtx2)` (`IntPatch_WLine.cxx:227-291`).
pub fn same_vtx_rst(vtx1: &PatchPoint, vtx2: &PatchPoint) -> bool {
    if vtx1.is_on_dom_s1() {
        if vtx2.is_on_dom_s1() {
            if vtx1.arc_id(true) == vtx2.arc_id(true) {
                if vtx1.param_on_arc(true) != vtx2.param_on_arc(true) {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else if vtx2.is_on_dom_s1() {
        return false;
    }
    if vtx1.is_on_dom_s2() {
        if vtx2.is_on_dom_s2() {
            if vtx1.arc_id(false) == vtx2.arc_id(false) {
                if vtx1.param_on_arc(false) != vtx2.param_on_arc(false) {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else if vtx2.is_on_dom_s2() {
        return false;
    }
    true
}

/// `CompareVerticesOnSurf(vtx1, vtx2, onFirst)` (`IntPatch_WLine.cxx:293-311`).
fn compare_vertices_on_surf(vtx1: &PatchPoint, vtx2: &PatchPoint, on_first: bool) -> bool {
    let (u1, v1, u2, v2) = if on_first {
        let (u1, v1) = vtx1.parameters_on_s1();
        let (u2, v2) = vtx2.parameters_on_s1();
        (u1, v1, u2, v2)
    } else {
        let (u1, v1) = vtx1.parameters_on_s2();
        let (u2, v2) = vtx2.parameters_on_s2();
        (u1, v1, u2, v2)
    };
    let tol_u = PCONFUSION;
    let tol_v = PCONFUSION;
    (u1 - u2).abs() <= tol_u && (v1 - v2).abs() <= tol_v
}

/// `CompareVerticesOnS1` (`IntPatch_WLine.cxx:313-316`).
pub fn compare_vertices_on_s1(vtx1: &PatchPoint, vtx2: &PatchPoint) -> bool {
    compare_vertices_on_surf(vtx1, vtx2, true)
}

/// `CompareVerticesOnS2` (`IntPatch_WLine.cxx:318-321`).
pub fn compare_vertices_on_s2(vtx1: &PatchPoint, vtx2: &PatchPoint) -> bool {
    compare_vertices_on_surf(vtx1, vtx2, false)
}
