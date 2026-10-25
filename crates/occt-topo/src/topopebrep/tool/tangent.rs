//! Tangents on edges and the "same oriented" curve test.
//!
//! Sources:
//! - `TopOpeBRepTool_TOOL::TggeomE` (`TKBool/TopOpeBRepTool/TopOpeBRepTool_TOOL.cxx:666-713`)
//! - `TopOpeBRepTool_TOOL::TgINSIDE` (`TKBool/TopOpeBRepTool/TopOpeBRepTool_TOOL.cxx:619-660`)
//! - `FUN_tool_dirC` (`TKBool/TopOpeBRepTool/TopOpeBRepTool_TOPOLOGY.cxx:1039-1046`)
//! - `FUN_tool_tggeomE` (`TopOpeBRepTool_TOPOLOGY.cxx:1049-1057`)
//! - `FUN_tool_curvesSO` (`TopOpeBRepTool_TOPOLOGY.cxx:1443-1548`)
//!
//! Every function returns `Option`/`bool` where the C++ returns `bool` plus
//! out-parameters, so the failure paths the C++ writes as `return false`
//! (`cxx:685-687`, `:701-712`) stay visible.

use occt_core::gp::{GpDir, GpVec};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

use super::adaptor::CurveAdaptor;
use super::closed::{
    check_edge_length, ori_in_sor, par_e, CLOSING, FORWARD, REVERSED,
};

/// `FUN_tool_dirC(const double par, const Handle(Geom_Curve)& C)`
/// (`TopOpeBRepTool_GEOMETRY.cxx:76-84`): `D1` then `gp_Dir(tgE)`.
pub fn fun_tool_dir_c(par: f64, curve: &dyn occt_geom::curve::Curve) -> Option<GpDir> {
    let (_, d1) = curve.d1(par);
    GpDir::from_vec(&d1).ok()
}

/// `TopOpeBRepTool_TOOL::TggeomE(par, const BRepAdaptor_Curve& BC, Tg)`
/// (`TOOL.cxx:666-694`).
///
/// `par` must be inside `[FirstParameter, LastParameter]` up to
/// `BC.Resolution(BC.Tolerance())`; the out tangent is normalised. The commented
/// `apoles` probe at `cxx:667-674` is debug-only and not translated.
pub fn tg_geom_e_curve(par: f64, bc: &CurveAdaptor) -> Option<GpVec> {
    let f = bc.first_parameter();
    let l = bc.last_parameter();
    let tolp = bc.resolution(bc.tolerance());

    let onf = (f - par).abs() < tolp;
    let onl = (l - par).abs() < tolp;
    let inbounds = f < par && par < l;
    if !inbounds && !onf && !onl {
        return None;
    }
    let (_, mut tg) = bc.d1(par)?;
    tg.normalize();
    Some(tg)
}

/// `TopOpeBRepTool_TOOL::TggeomE(par, const TopoDS_Edge& E, Tg)`
/// (`TOOL.cxx:698-713`).
pub fn tg_geom_e(par: f64, e: &TopoShape) -> Option<GpVec> {
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return None;
    }
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    if !check_edge_length(e) {
        return None;
    }
    tg_geom_e_curve(par, &bc)
}

/// `FUN_tool_dirC(const double par, const BRepAdaptor_Curve& BAC)`
/// (`TopOpeBRepTool_TOPOLOGY.cxx:1039-1046`).
pub fn fun_tool_dir_c_adaptor(par: f64, bc: &CurveAdaptor) -> Option<GpDir> {
    let (_, d1) = bc.d1(par)?;
    GpDir::from_vec(&d1).ok()
}

/// `FUN_tool_tggeomE(paronE, E)` (`TopOpeBRepTool_TOPOLOGY.cxx:1049-1057`): a
/// degenerated edge answers the zero vector instead of failing.
pub fn fun_tool_tg_geom_e(par_on_e: f64, e: &TopoShape) -> GpVec {
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return GpVec::zero();
    }
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    match fun_tool_dir_c_adaptor(par_on_e, &bc) {
        Some(d) => GpVec::new(d.x(), d.y(), d.z()),
        None => GpVec::zero(),
    }
}

/// `TopOpeBRepTool_TOOL::TgINSIDE(v, E, Tg, OvinE)` (`TOOL.cxx:619-660`).
///
/// Returns the tangent at `v` on `E` together with the `OriinSor` code of `v`
/// in the `FORWARD`-oriented edge. The tangent is reversed for a `REVERSED`
/// occurrence, and a `CLOSING` occurrence is evaluated at the edge's first
/// parameter.
pub fn tg_inside(v: &TopoShape, e: &TopoShape) -> Option<(GpVec, i32)> {
    let efor = e.oriented(Orientation::Forward);
    let ov_e = ori_in_sor(v, &efor, true);
    if ov_e == 0 {
        return None;
    }
    let iv = if ov_e == CLOSING {
        FORWARD
    } else if ov_e == FORWARD || ov_e == REVERSED {
        ov_e
    } else {
        0
    };
    let par_e = if iv == 0 {
        BRepTool::parameter_of_vertex(e, v)
    } else {
        par_e(iv, &efor)?
    };
    let mut tg = tg_geom_e(par_e, &efor)?;
    if ov_e == REVERSED {
        tg.reverse();
    }
    Some((tg, ov_e))
}

/// `FUN_tool_curvesSO(E1, p1, E2, p2, so)` (`TOPOLOGY.cxx:1443-1479`).
///
/// `Precision::Angular() * 1.e3` is the parallel/opposite tolerance; a
/// non-parallel pair fails rather than guessing.
pub fn curves_same_oriented_at(e1: &TopoShape, p1: f64, e2: &TopoShape, p2: f64) -> Option<bool> {
    let tg1 = {
        let bc = CurveAdaptor::of_edge(&Edge(e1.clone()));
        let d = fun_tool_dir_c_adaptor(p1, &bc)?;
        GpVec::new(d.x(), d.y(), d.z())
    };
    let tg2 = {
        let bc = CurveAdaptor::of_edge(&Edge(e2.clone()));
        let d = fun_tool_dir_c_adaptor(p2, &bc)?;
        GpVec::new(d.x(), d.y(), d.z())
    };
    let tola = occt_core::precision::ANGULAR * 1.0e3;
    if tg1.is_opposite_ang(&tg2, tola) {
        return Some(false);
    }
    if tg1.is_parallel_ang(&tg2, tola) {
        return Some(true);
    }
    None
}

/// `FUN_tool_curvesSO(E1, p1, E2, so)` (`TOPOLOGY.cxx:1481-1495`): projects
/// `p1` onto `E2` first.
pub fn curves_same_oriented_point(e1: &TopoShape, p1: f64, e2: &TopoShape) -> Option<bool> {
    let p2 = super::project_curve::fun_tool_par_e(e1, p1, e2)?;
    curves_same_oriented_at(e1, p1, e2, p2)
}

/// `FUN_tool_curvesSO(E1, E2, so)` (`TOPOLOGY.cxx:1498-1548`).
///
/// Two non-closed edges decide from their shared vertices; anything with a
/// closing vertex, or without a shared vertex, projects a point near the first
/// end (`x = 0.45678`) of `E1` onto `E2` and compares tangents there.
pub fn curves_same_oriented(e1: &TopoShape, e2: &TopoShape) -> Option<bool> {
    let (vf1, vl1) = crate::topo_tools::edge_vertices(&Edge(e1.clone()));
    let (vf2, vl2) = crate::topo_tools::edge_vertices(&Edge(e2.clone()));
    let same = |a: &Option<crate::shape::Vertex>, b: &Option<crate::shape::Vertex>| match (a, b) {
        (Some(a), Some(b)) => crate::topo_tools::is_same(&a.0, &b.0),
        _ => false,
    };
    let closed1 = same(&vf1, &vl1);
    let closed2 = same(&vf2, &vl2);
    let project;
    let so;
    if closed1 || closed2 {
        project = true;
        so = false;
    } else if same(&vf1, &vf2) {
        project = false;
        so = true;
    } else if same(&vl1, &vl2) {
        project = false;
        so = true;
    } else if same(&vf1, &vl2) {
        project = false;
        so = false;
    } else if same(&vl1, &vf2) {
        project = false;
        so = false;
    } else {
        project = true;
        so = false;
    }
    if !project {
        return Some(so);
    }
    let (f, l) = BRepTool::edge_parameters(&Edge(e1.clone()));
    let x = 0.45678;
    let p1 = (1.0 - x) * l + x * f;
    curves_same_oriented_point(e1, p1, e2)
}
