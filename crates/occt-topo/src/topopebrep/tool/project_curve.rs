//! Point projection onto an edge, and the parameter transfer built on it.
//!
//! Sources:
//! - `FUN_tool_bounds` (`TKBool/TopOpeBRepTool/TopOpeBRepTool_PROJECT.cxx:33-38`)
//! - `FUN_tool_getindex` (`TopOpeBRepTool_PROJECT.cxx:41-75`)
//! - `FUN_tool_projPonC` (`TopOpeBRepTool_PROJECT.cxx:75-144`)
//! - `FUN_tool_projPonE` (declared in `TopOpeBRepTool_PROJECT.hxx:108-118`; its
//!   body is the `BAC` overload applied to `BRepAdaptor_Curve(E)`)
//! - `FUN_tool_value`, `FUN_tool_parE` (`TopOpeBRepTool_TOPOLOGY.cxx:313-325`,
//!   `:932-966`)
//!
//! `Extrema_ExtPC` is `crate::...::extrema_pc`'s `extrema_ext_pc_range_done`,
//! which keeps OCCT's `IsDone` alongside the solution list because the
//! free-form arm can report `false` with a non-empty list
//! (`Extrema_GGExtPC.hxx:468`).

use occt_core::gp::GpPnt;
use occt_geom::extrema_pc::{extrema_ext_pc_range_done, ExtPcSolution};

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

use super::adaptor::CurveAdaptor;

/// `FUN_tool_bounds(E, f, l)` (`TopOpeBRepTool_PROJECT.cxx:33-38`).
pub fn fun_tool_bounds(e: &TopoShape) -> (f64, f64) {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    (bc.first_parameter(), bc.last_parameter())
}

/// `FUN_tool_getindex(ponc)` (`TopOpeBRepTool_PROJECT.cxx:41-54`): index of the
/// smallest square distance, 0-based here. The C++ starts at `1` and only
/// replaces on a strict `<`, so the first solution wins ties.
pub(crate) fn fun_tool_getindex(sols: &[ExtPcSolution]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, s) in sols.iter().enumerate() {
        match best {
            None => best = Some(i),
            Some(b) if s.sq_dist < sols[b].sq_dist => best = Some(i),
            _ => {}
        }
    }
    best
}

/// `FUN_tool_projPonC(P, tole, BAC, pmin, pmax, param, dist)`
/// (`TopOpeBRepTool_PROJECT.cxx:76-114`).
///
/// When `Extrema_ExtPC` reports no solution, the two interval ends are tested
/// against `tole` and one of them is returned if it lies within tolerance.
pub fn fun_tool_proj_pon_c_range(
    p: &GpPnt,
    tole: f64,
    bc: &CurveAdaptor,
    pmin: f64,
    pmax: f64,
) -> Option<(f64, f64)> {
    let (sols, done) = match bc.curve() {
        Some(c) => extrema_ext_pc_range_done(c.as_ref(), p, pmin, pmax),
        None => (Vec::new(), false),
    };
    if done && !sols.is_empty() {
        let i = fun_tool_getindex(&sols)?;
        return Some((sols[i].u, sols[i].sq_dist.sqrt()));
    }
    for par in [pmin, pmax] {
        let pt = bc.d0(par)?;
        let d2 = pt.square_distance(p);
        if d2 < tole * tole {
            return Some((par, d2.sqrt()));
        }
    }
    None
}

/// `FUN_tool_projPonC(P, BAC, pmin, pmax, param, dist)`
/// (`TopOpeBRepTool_PROJECT.cxx:118-126`).
pub fn fun_tool_proj_pon_c_tol(
    p: &GpPnt,
    bc: &CurveAdaptor,
    pmin: f64,
    pmax: f64,
) -> Option<(f64, f64)> {
    fun_tool_proj_pon_c_range(p, bc.tolerance(), bc, pmin, pmax)
}

/// `FUN_tool_projPonC(P, BAC, param, dist)`
/// (`TopOpeBRepTool_PROJECT.cxx:132-144`).
pub fn fun_tool_proj_pon_c(p: &GpPnt, bc: &CurveAdaptor) -> Option<(f64, f64)> {
    fun_tool_proj_pon_c_tol(p, bc, bc.first_parameter(), bc.last_parameter())
}

/// `FUN_tool_projPonE(P, tole, E, param, dist)`
/// (`TopOpeBRepTool_PROJECT.cxx:260-290`): the `tole` overload; a parameter
/// within `Parametric(Confusion())` of either bound is snapped to it.
pub fn fun_tool_proj_pon_e_tol(
    p: &GpPnt,
    tole: f64,
    e: &TopoShape,
) -> Option<(f64, f64)> {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let (mut param, dist) =
        fun_tool_proj_pon_c_range(p, tole, &bc, bc.first_parameter(), bc.last_parameter())?;
    let (f, l) = fun_tool_bounds(e);
    let tolp = occt_core::precision::PCONFUSION;
    if (f - param).abs() < tolp {
        param = f;
    }
    if (l - param).abs() < tolp {
        param = l;
    }
    Some((param, dist))
}

/// `FUN_tool_projPonE(P, E, param, dist)`: the `BRepAdaptor_Curve` overload on
/// an edge.
pub fn fun_tool_proj_pon_e(p: &GpPnt, e: &TopoShape) -> Option<(f64, f64)> {
    fun_tool_proj_pon_e_tol(p, BRepTool::edge_tolerance(&Edge(e.clone())), e)
}

/// `FUN_tool_value(par, E, P)` (`TopOpeBRepTool_TOPOLOGY.cxx:313-325`): the
/// point is only evaluated when `par` is inside the closed range.
pub fn fun_tool_value(par: f64, e: &TopoShape) -> Option<GpPnt> {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let f = bc.first_parameter();
    let l = bc.last_parameter();
    if f <= par && par <= l {
        bc.d0(par)
    } else {
        None
    }
}

/// `FUN_tool_parE(E0, par0, E, par, tol)` (`TopOpeBRepTool_TOPOLOGY.cxx:932-954`).
///
/// The point at `par0` on `E0` is projected onto `E`; the projection only
/// counts when it lands within `tol` of the point.
pub fn fun_tool_par_e_tol(
    e0: &TopoShape,
    par0: f64,
    e: &TopoShape,
    tol: f64,
) -> Option<f64> {
    let p = fun_tool_value(par0, e0)?;
    let (par, dist) = fun_tool_proj_pon_e(&p, e)?;
    if dist < tol {
        Some(par)
    } else {
        None
    }
}

/// `FUN_tool_parE(E0, par0, E, par)` (`TopOpeBRepTool_TOPOLOGY.cxx:956-966`):
/// the tolerance is `BRep_Tool::Tolerance(E) * 1.e2`, kept verbatim from
/// `cxx:964` (the C++ itself flags it as provisional).
pub fn fun_tool_par_e(e0: &TopoShape, par0: f64, e: &TopoShape) -> Option<f64> {
    let tol3d = BRepTool::edge_tolerance(&Edge(e.clone())) * 1.0e2;
    fun_tool_par_e_tol(e0, par0, e, tol3d)
}
