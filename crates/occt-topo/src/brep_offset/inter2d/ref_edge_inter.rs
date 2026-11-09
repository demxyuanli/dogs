//! `RefEdgeInter` (`BRepOffset_Inter2d.cxx:703-1070`): the reference-vertex
//! variant of `EdgeInter`, used while intersecting the edges of a vertex
//! neighbourhood (`ConnexIntByIntInVert`).
//!
//! Compared with `EdgeInter` the reference
//! * returns early for null edges (`cxx:726-729`) and marks `theCoincide` when
//!   the two pcurves are line-coincident (`cxx:779-793`) or when the only
//!   intersection is a line-line segment (`cxx:798-803`);
//! * overrides the computed orientations with `theOr1` / `theOr2` unless they
//!   are `TopAbs_EXTERNAL` (`cxx:907-916`);
//! * stores `Tol` instead of the vertex-pair tolerance in the endpoint test
//!   (`cxx:959-961`) and in the double purge (`cxx:992`);
//! * keeps only the solution closest to `theVref` when there is more than one
//!   (`cxx:1018-1047`) and records every kept vertex as an image of
//!   `theVref` in `theImageVV` (`cxx:1049-1062`);
//! * raises the store tolerance to the line-line `TolLL` (`cxx:1065-1069`).
//!
//! This function is driven by `ConnexIntByInt` and `ConnexIntByIntInVert`
//! (`cxx:1830-2094`, `:2100-2308`, paid as
//! [`connex_int_by_int`](super::connex_int_by_int) /
//! [`connex_int_by_int_in_vert`](super::connex_int_by_int_in_vert)).

use crate::abs::Orientation;
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_tool::BRepTool;
use crate::int_tools_curve_box::compute_int_range;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{is_same, vertex_position};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use occt_geom2d::geom2d_int::Geom2dIntGInter;

use super::edge_inter::{
    build_curve3d_defaults, collect_intersection_points, degenerate_point, make_intersection_vertex,
    prepend_endpoint_vertices, purge_doubles, TOL_DUB,
};
use super::store_vertices::store_on_two_edges;
use super::Dmvv;

/// `RefEdgeInter(F, BAsurf, E1, E2, theOr1, theOr2, AsDes, Tol, WithOri,
/// theVref, theImageVV, aDMVV, theCoincide)` (`cxx:703-1070`).
///
/// `theOr1` / `theOr2` are the `TopAbs_Orientation` overrides; pass
/// [`Orientation::External`] for the "keep the computed orientation" arm of
/// `cxx:907-916`.
///
/// Port note: the reference's `E1.IsNull() || E2.IsNull()` guard (`cxx:726-729`)
/// has no counterpart here - a [`Edge`] cannot hold a null shape, so the caller
/// has no null edge to pass.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ref_edge_inter(
    f: &Face,
    basurf: &SurfaceAdaptor,
    e1: &Edge,
    e2: &Edge,
    or1: Orientation,
    or2: Orientation,
    as_des: &mut AsDes,
    tol: f64,
    with_ori: bool,
    v_ref: &Vertex,
    image_vv: &mut BRepAlgoImage,
    dmvv: &mut Dmvv,
    coincide: &mut bool,
) {
    // `cxx:719-721`.
    *coincide = false;

    // `cxx:723-729`.
    if is_same(&e1.0, &e2.0) {
        return;
    }

    // `cxx:731-741`: the two pcurves; a null handle returns.
    let (Some((pc1, f1, l1)), Some((pc2, f2, l2))) = (
        crate::boptools_2d::curve_on_surface_range(e1, f),
        crate::boptools_2d::curve_on_surface_range(e2, f),
    ) else {
        return;
    };

    // `cxx:749-751`.
    build_curve3d_defaults(e1);
    build_curve3d_defaults(e2);

    // `cxx:753-775`.
    let deg_point = degenerate_point(basurf, f, e1, e2);

    // `cxx:777-793`: two lines. A parallel (or anti-parallel) pair is
    // "coincident" and returns without intersecting; otherwise the line-line
    // range is taken into the store tolerance.
    let mut tol_ll = 0.0;
    if pc1.is_geom2d_line() && pc2.is_geom2d_line() {
        let (Some(lin1), Some(lin2)) = (pc1.gp_lin2d(), pc2.gp_lin2d()) else {
            return;
        };
        let angle = lin1.direction().angle(&lin2.direction()).abs();
        if angle <= 1.0e-8 || std::f64::consts::PI - angle <= 1.0e-8 {
            *coincide = true;
            return;
        }
        // `TolLL = IntTools_Tools::ComputeIntRange(TolDub, TolDub, anAngle);`
        // `TolLL = std::min(TolLL, 1.e-5);`
        tol_ll = compute_int_range(TOL_DUB, TOL_DUB, angle).min(1.0e-5);
    }

    // `cxx:796-803`.
    let mut inter = Geom2dIntGInter::new();
    let d1 = Geom2dIntGInter::domain_over_range(pc1.as_ref(), f1, l1, TOL_DUB);
    let d2 = Geom2dIntGInter::domain_over_range(pc2.as_ref(), f2, l2, TOL_DUB);
    inter.perform(pc1.as_ref(), &d1, pc2.as_ref(), &d2, TOL_DUB, TOL_DUB);
    if !inter.is_done() || inter.nb_points() == 0 {
        *coincide = inter.nb_segments() > 0 && pc1.is_geom2d_line() && pc2.is_geom2d_line();
        return;
    }

    // `cxx:805-820`.
    let (res_points, res_params_e1, res_params_e2) =
        collect_intersection_points(basurf, &inter, deg_point);

    let mut lv1: Vec<TopoShape> = Vec::new();
    let mut lv2: Vec<TopoShape> = Vec::new();

    // `cxx:822-918`.
    for i in 0..res_points.len() {
        let t1 = res_params_e1[i];
        let t2 = res_params_e2[i];
        if !t1.is_finite() || !t2.is_finite() {
            continue;
        }
        let Some(nv) = make_intersection_vertex(
            f,
            basurf,
            e1,
            e2,
            pc1.as_ref(),
            pc2.as_ref(),
            t1,
            t2,
            &res_points[i],
            tol,
            with_ori,
        ) else {
            continue;
        };
        // `cxx:907-916`.
        let oo1 = if or1 != Orientation::External { or1 } else { nv.or1 };
        let oo2 = if or2 != Orientation::External { or2 } else { nv.or2 };
        lv1.push(nv.vertex.oriented(oo1));
        lv2.push(nv.vertex.oriented(oo2));
    }

    // `cxx:925-967`: `Tol` in both `UpdateVertex` calls, not the vertex pair's.
    prepend_endpoint_vertices(
        e1, e2, f1, l1, f2, l2, tol, false, as_des, &mut lv1, &mut lv2,
    );

    if lv1.is_empty() {
        return;
    }

    // `cxx:970-1013`.
    purge_doubles(&mut lv1, &mut lv2, tol, false);

    // `cxx:1018-1047`: keep the solution closest to the reference vertex.
    if lv1.len() > 1 {
        let pref = vertex_position(v_ref);
        let mut dmin = f64::MAX;
        let mut vmin: Option<TopoShape> = None;
        for v in &lv1 {
            let p = vertex_position(&Vertex(v.clone()));
            let d = p.square_distance(&pref);
            if d < dmin {
                dmin = d;
                vmin = Some(v.clone());
            }
        }
        if let Some(vmin) = vmin {
            // `if (!Vmin.IsSame(it1LV1.Value())) { LV1.Remove(it1LV1);
            // LV2.Remove(it1LV2); }` in list order.
            let mut i = 0;
            while i < lv1.len() {
                if !is_same(&lv1[i], &vmin) {
                    lv1.remove(i);
                    lv2.remove(i);
                } else {
                    i += 1;
                }
            }
        }
    }

    // `cxx:1049-1062`.
    for v in &lv1 {
        let mut new_vertex = v.clone();
        new_vertex.set_orientation(Orientation::Forward);
        let key = v_ref.0.oriented(Orientation::Forward);
        if image_vv.has_image(&key) {
            let _ = image_vv.add(&key, &new_vertex);
        } else {
            let _ = image_vv.bind(&key, &new_vertex);
        }
    }

    // `cxx:1064-1069`.
    let tol_store = (BRepTool::edge_tolerance(e1) + BRepTool::edge_tolerance(e2))
        .max(tol)
        .max(tol_ll);
    store_on_two_edges(&e1.0, &e2.0, &lv1, &lv2, tol_store, as_des, dmvv);
}


