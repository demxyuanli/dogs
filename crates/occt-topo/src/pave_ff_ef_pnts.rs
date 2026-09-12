//! `BOPAlgo_PaveFiller::GetEFPnts` (`_6.cxx:2608-2688`).
//!
//! For every E/F interference whose new vertex is a sub-shape of either face
//! of an F/F pair, evaluate the 3D point at `IntTools_CommonPrt::VertexParameter1`
//! (`IntTools_Tools::VertexParameter` when the stored parameter is inside
//! `Range1`) and project it onto the opposite face. When the edge already
//! carries a p-curve on the non-opposite face, that UV is paired with the
//! opposite-face projection; otherwise both faces are projected.

use std::collections::HashSet;

use occt_core::gp::GpPnt;
use occt_geom2d::Curve2d;

use crate::bopds::BopdsDS;
use crate::boptools_2d::curve_on_surface;
use crate::brep_tool::BRepTool;
use crate::int_tools_curve_box::vertex_parameter;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_data::IntRange;
use crate::pave_ff_misc::get_full_shape_map;
use crate::pave_ff_paves::PntOn2s;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Face};

fn full_shape_map(ds: &BopdsDS, n_f: usize, mi: &mut HashSet<usize>) {
    get_full_shape_map(ds, n_f, mi);
}

fn vertex_parameter1(common_first: f64, common_last: f64) -> f64 {
    let lo = common_first.min(common_last);
    let hi = common_first.max(common_last);
    vertex_parameter(IntRange::new_unchecked(lo, hi), 0.5 * (lo + hi))
}

fn push_pnt(
    n_f: usize,
    n_f1: usize,
    uv_on_f: (f64, f64),
    uv_on_opp: (f64, f64),
    out: &mut Vec<PntOn2s>,
) {
    let mut a_pnt = PntOn2s {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
    };
    if n_f == n_f1 {
        a_pnt.set_value(uv_on_f.0, uv_on_f.1, uv_on_opp.0, uv_on_opp.1);
    } else {
        a_pnt.set_value(uv_on_opp.0, uv_on_opp.1, uv_on_f.0, uv_on_f.1);
    }
    out.push(a_pnt);
}

fn project_with_pcurve(
    ctx: &IntToolsContext,
    a_par: f64,
    a_pcurve: &dyn Curve2d,
    a_point: &GpPnt,
    a_f_opp: &Face,
    n_f: usize,
    n_f1: usize,
    out: &mut Vec<PntOn2s>,
) {
    let a_p2d = a_pcurve.d0(a_par);
    if let Ok((u1, v1)) = ctx.project_point_on_face(a_f_opp, a_point) {
        push_pnt(n_f, n_f1, (a_p2d.x(), a_p2d.y()), (u1, v1), out);
    }
}

fn project_both_faces(
    ctx: &IntToolsContext,
    a_point: &GpPnt,
    a_f: &Face,
    a_f_opp: &Face,
    n_f: usize,
    n_f1: usize,
    out: &mut Vec<PntOn2s>,
) {
    if let (Ok((u1, v1)), Ok((u2, v2))) = (
        ctx.project_point_on_face(a_f, a_point),
        ctx.project_point_on_face(a_f_opp, a_point),
    ) {
        push_pnt(n_f, n_f1, (u1, v1), (u2, v2), out);
    }
}

/// One E/F interference that produced a new vertex (`HasIndexNew`).
struct EfNewVertex {
    n_e: usize,
    n_f_opp: usize,
    a_par: f64,
}

fn collect_ef_new_vertices(f: &PaveFiller, mi: &HashSet<usize>) -> Vec<EfNewVertex> {
    let mut out = Vec::new();
    for it in f.ds().interf_ef() {
        if it.get_index_new().is_none() {
            continue;
        }
        let n_e = it.index1();
        let n_f_opp = it.index2();
        if !mi.contains(&n_e) || !mi.contains(&n_f_opp) {
            continue;
        }
        let a_par = vertex_parameter1(it.common_first, it.common_last);
        out.push(EfNewVertex {
            n_e,
            n_f_opp,
            a_par,
        });
    }
    out
}

fn append_ef_point(
    f: &PaveFiller,
    n_f1: usize,
    n_f2: usize,
    rec: &EfNewVertex,
    out: &mut Vec<PntOn2s>,
) {
    let Some(e_shape) = f.ds().shape(rec.n_e).cloned() else {
        return;
    };
    let Some(f_opp_shape) = f.ds().shape(rec.n_f_opp).cloned() else {
        return;
    };
    let a_e = Edge(e_shape);
    let a_f_opp = Face(f_opp_shape);
    let Some(a_curve) = BRepTool::edge_curve(&a_e) else {
        return;
    };
    let n_f = if rec.n_f_opp == n_f1 { n_f2 } else { n_f1 };
    let Some(f_shape) = f.ds().shape(n_f).cloned() else {
        return;
    };
    let a_f = Face(f_shape);
    let a_point = a_curve.d0(rec.a_par);
    let ctx = f.context();
    if let Some(a_pcurve) = curve_on_surface(&a_e, &a_f) {
        project_with_pcurve(
            ctx,
            rec.a_par,
            a_pcurve.as_ref(),
            &a_point,
            &a_f_opp,
            n_f,
            n_f1,
            out,
        );
    } else {
        project_both_faces(ctx, &a_point, &a_f, &a_f_opp, n_f, n_f1, out);
    }
}

/// `BOPAlgo_PaveFiller::GetEFPnts` (`_6.cxx:2608`).
pub fn get_ef_pnts(f: &PaveFiller, n_f1: usize, n_f2: usize) -> Vec<PntOn2s> {
    let mut mi = HashSet::new();
    full_shape_map(f.ds(), n_f1, &mut mi);
    full_shape_map(f.ds(), n_f2, &mut mi);
    let recs = collect_ef_new_vertices(f, &mi);
    let mut out = Vec::new();
    for rec in &recs {
        append_ef_point(f, n_f1, n_f2, rec, &mut out);
    }
    out
}
