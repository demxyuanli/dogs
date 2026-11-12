//! `BOPTools_AlgoTools2D::AttachExistingPCurve` and its closed-seam helper
//! `UpdateClosedPCurve` (`BOPTools_AlgoTools2D_1.cxx:44-311`).
//!
//! Distinct from [`crate::boptools_2d::attach_existing_pcurve`], which is only
//! the `BRep_Builder::UpdateEdge(edge, curve2d, face, tol)` registration step.
//! This module ports the full copy path: reverse/SameRange/ComputeTolerance
//! gates, a temporary edge, `BRepLib::SameParameter`/`SameRange`, optional
//! closed-pcurve rebuild, then `BRep_Builder::Transfert` onto the new edge.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpTrsf2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::geom_lib_same_range_2d::same_range_2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_topology::brep_builder_transfert::transfert_edge;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools_face::is_split_to_reverse_edge;
use crate::boptools_2d::{curve_on_surface_range, intermediate_point};
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_lib_same_range::same_range;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::cumulated_children;
use crate::int_tools_full::IntToolsContext;
use crate::int_tools_segpln::compute_tolerance;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same};

/// `BOPTools_AlgoTools2D::AttachExistingPCurve` (`cxx:44-160`).
///
/// Copies the p-curve of `e2_old` onto `e1_new` for `face`, returning the OCCT
/// `iRet` codes (`0` success; `1`..`6` the documented failure arms).
pub fn attach_existing_pcurve_edges(
    e2_old: &Edge,
    e1_new: &Edge,
    face: &Face,
    ctx: &IntToolsContext,
) -> i32 {
    let a_f = Face(face.0.oriented(Orientation::Forward));
    let a_e1 = Edge(e1_new.0.oriented(Orientation::Forward));
    let a_e2 = Edge(e2_old.0.oriented(Orientation::Forward));

    let Some((a_c2d_old, mut a_t21, mut a_t22)) = curve_on_surface_range(&a_e2, &a_f) else {
        return 1;
    };

    let mut a_c2d_old_c = a_c2d_old.clone_dyn();
    let b_is_to_reverse = is_split_to_reverse_edge(&a_e1, &a_e2, ctx).unwrap_or(false);
    if b_is_to_reverse {
        a_c2d_old_c.reverse();
        let a_t21r = a_c2d_old_c.reversed_parameter(a_t21);
        let a_t22r = a_c2d_old_c.reversed_parameter(a_t22);
        a_t21 = a_t22r;
        a_t22 = a_t21r;
    }
    let a_c2d_old_c: Arc<dyn Curve2d> = Arc::from(a_c2d_old_c);

    let mut a_c2dt: Arc<dyn Curve2d> =
        Arc::new(Geom2dTrimmedCurve::new(a_c2d_old_c, a_t21, a_t22));

    let Some(a_ce1) = BRepTool::edge_curve(&a_e1) else {
        return 2;
    };
    let (a_t11, a_t12) = BRepTool::edge_parameters(&a_e1);

    match same_range_2d(PCONFUSION, &a_c2dt, a_t21, a_t22, a_t11, a_t12) {
        Ok(pc) => a_c2dt = pc,
        Err(_) => return 2,
    }

    let Some(a_sf) = BRepTool::face_surface(&a_f) else {
        return 3;
    };
    let Some((a_tol_sp, _a_t_max)) =
        compute_tolerance(&a_ce1, &a_c2dt, &a_sf, a_t11, a_t12)
    else {
        return 3;
    };

    let a_tol = BRepTool::edge_tolerance(&a_e1);
    if (a_tol_sp > 10.0 * a_tol) && a_tol_sp > 0.1 {
        return 4;
    }

    let builder = TopoBuilder::new();
    let a_e1t = builder.make_edge(a_ce1.clone(), a_t11, a_t12);
    let reg = GeometryRegistry::global();
    reg.set_edge_tolerance(&a_e1t.0, a_tol);
    reg.set_ranges_all(&a_e1t.0, a_t11, a_t12);
    reg.set_same_range(&a_e1t.0, false);
    reg.set_same_parameter(&a_e1t.0, false);

    let face_key = GeometryRegistry::shape_key(&a_f.0);
    reg.set_edge_pcurve(&a_e1t.0, face_key, a_c2dt);
    reg.raise_edge_tolerance(&a_e1t.0, a_tol);

    // `cxx:128-136`: `BRepLib::SameParameter` / `SameRange` catch `Standard_Failure`.
    // The port has no throw path; still gate the calls so a future panic maps to 6.
    let sp_ok = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        same_parameter_edge_inplace(&a_e1t, CONFUSION);
        same_range(&a_e1t.0, CONFUSION);
    }));
    if sp_ok.is_err() {
        return 6;
    }

    if is_closed_for_attach(&a_e2, &a_f) {
        let i_ret = update_closed_pcurve(&a_e2, &a_e1t, &a_f, ctx);
        if i_ret != 0 {
            return 5;
        }
    }

    // `cxx:151`: `aBB.Transfert(aE1T, aE1)`.
    transfert_edge(&a_e1t, &a_e1);

    let a_new_tol = BRepTool::edge_tolerance(&a_e1t);
    for child in cumulated_children(&a_e1.0) {
        if child.is_vertex() {
            reg.raise_vertex_tolerance(&child, a_new_tol);
        }
    }
    0
}

/// `UpdateClosedPCurve` (`BOPTools_AlgoTools2D_1.cxx:164-285`).
fn update_closed_pcurve(
    a_e_old: &Edge,
    a_e_new: &Edge,
    a_f: &Face,
    ctx: &IntToolsContext,
) -> i32 {
    let a_tol = BRepTool::edge_tolerance(a_e_new);

    let Some((a_c2d_old_ct, a_t1, a_t2)) = curve_on_surface_range(a_e_new, a_f) else {
        return 1;
    };
    // `aC2Dold` is read in OCCT but unused after; kept for control-flow parity.
    let _a_c2d_old = curve_on_surface_range(a_e_old, a_f);

    let a_es_fwd = Edge(a_e_old.0.oriented(Orientation::Forward));
    let a_es_rev = Edge(a_e_old.0.oriented(Orientation::Reversed));
    let Some((a_c2ds1, a_ts1, a_ts2)) = curve_on_surface_range(&a_es_fwd, a_f) else {
        return 1;
    };
    // `cxx:200-203`: REVERSED edge selects PCurve2 via `BRep_Tool::CurveOnSurface`.
    let Some((a_c2ds2, _, _)) = curve_on_surface_range(&a_es_rev, a_f) else {
        return 1;
    };

    let a_ts = intermediate_point(a_ts1, a_ts2);
    let (a_p2ds1, a_v2ds1) = a_c2ds1.d1(a_ts);
    let (a_p2ds2, _a_v2ds2) = a_c2ds2.d1(a_ts);

    let a_v2ds12 = GpVec2d::new(a_p2ds2.x() - a_p2ds1.x(), a_p2ds2.y() - a_p2ds1.y());
    let Ok(a_d2ds12) = GpDir2d::from_vec2d(&a_v2ds12) else {
        return 1;
    };
    let a_d2dx = GpDir2d::new(1.0, 0.0).unwrap_or_default();
    let a_sc_pr = a_d2ds12.dot(&a_d2dx);
    let mut b_u_closed = true;
    if a_sc_pr.abs() < a_tol {
        b_u_closed = !b_u_closed;
    }

    let (a_us1, a_vs1) = (a_p2ds1.x(), a_p2ds1.y());
    let (a_us2, a_vs2) = (a_p2ds2.x(), a_p2ds2.y());

    let Some(a_s) = BRepTool::face_surface(a_f) else {
        return 1;
    };
    let a_p = a_s.d0(a_us1, a_vs1);

    let Some(a_t) = ctx.project_point_on_edge(a_e_new, &a_p) else {
        return 1;
    };

    let a_c2d = a_c2d_old_ct;
    let (a_p2d, a_v2d) = a_c2d.d1(a_t);
    let (a_u, a_v) = (a_p2d.x(), a_p2d.y());

    let a_c2d_new: Arc<dyn Curve2d> = Arc::from(a_c2d.clone_dyn());
    let mut a_c2d_tnew = Geom2dTrimmedCurve::new(a_c2d_new, a_t1, a_t2);

    let mut a_v2dt = a_v2ds12;
    if !b_u_closed {
        if (a_v - a_vs2).abs() < a_tol {
            a_v2dt = a_v2dt.reversed();
        }
    } else if (a_u - a_us2).abs() < a_tol {
        a_v2dt = a_v2dt.reversed();
    }

    let mut trsf = GpTrsf2d::identity();
    trsf.set_translation_vec(&a_v2dt);
    a_c2d_tnew.transform(&trsf);
    let a_c2d_tnew: Arc<dyn Curve2d> = Arc::new(a_c2d_tnew);

    let mut b_rev_order = false;
    let a_sc_pr = a_v2d.dot(&a_v2ds1);
    if a_sc_pr < 0.0 {
        b_rev_order = !b_rev_order;
    }

    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&a_f.0);
    if !b_rev_order {
        reg.set_edge_pcurves(
            &a_e_new.0,
            face_key,
            vec![a_c2d.clone(), a_c2d_tnew],
        );
    } else {
        reg.set_edge_pcurves(
            &a_e_new.0,
            face_key,
            vec![a_c2d_tnew, a_c2d],
        );
    }
    reg.raise_edge_tolerance(&a_e_new.0, a_tol);
    0
}

/// File-local `IsClosed` (`BOPTools_AlgoTools2D_1.cxx:287-311`):
/// `BRep_Tool::IsClosed(E, F)` and the edge appears twice under the face.
fn is_closed_for_attach(a_e: &Edge, a_f: &Face) -> bool {
    if !BRepTool::is_closed_edge_face(a_e, a_f) {
        return false;
    }
    let mut i_cnt = 0;
    for a_ex in explorer_shapes(&a_f.0, ShapeType::Edge) {
        if is_same(&a_ex, &a_e.0) {
            i_cnt += 1;
        }
    }
    i_cnt == 2
}
