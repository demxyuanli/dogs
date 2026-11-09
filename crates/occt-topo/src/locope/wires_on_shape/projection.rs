//! The three file-local `Project` helpers of `LocOpe_WiresOnShape.cxx`.
//!
//! * `Project(V, gp_Pnt2d, F, Edge&, Param&)` (`cxx:491-676`) - pick the face
//!   edge the 2d vertex parameter belongs to, then widen the vertex tolerance
//!   when the 2d gap does not fit inside the surface resolution.
//! * `Project(V, Edge)` (`cxx:679-698`) - `GeomAPI_ProjectPointOnCurve` over the
//!   edge range.
//! * `Project(V, gp_Pnt2d, Edge, Face)` (`cxx:702-718`) -
//!   `Geom2dAPI_ProjectPointOnCurve` over the p-curve range.
//!
//! All three answer `Precision::Infinite()` for "no projection", exactly as the
//! reference does.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{self, CONFUSION, INFINITE, PCONFUSION};
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::{edges_of, explorer_shapes};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `Project(V, p2d, F, theEdge, param)` (`LocOpe_WiresOnShape.cxx:491-676`).
///
/// `the_edge` is OCCT's `TopoDS_Edge&`: a null `TopoDS_Edge` is `None`. On entry
/// with a null edge the helper searches every edge of `F` for the closest one;
/// with an edge already chosen it only fills in `param` when that is infinite.
pub(super) fn project_vertex_on_face(
    v: &Vertex,
    p2d: &GpPnt2d,
    f: &Face,
    the_edge: &mut Option<Edge>,
    param: &mut f64,
) -> bool {
    let a_tol_v = BRepTool::vertex_tolerance(v);
    // `cxx:496-497`: `dmin = theEdge.IsNull() ? RealLast() : aTolV * aTolV`.
    let mut dmin = match the_edge {
        None => f64::MAX,
        Some(_) => a_tol_v * a_tol_v,
    };
    let mut valret = false;

    let Some(a_surf) = BRepTool::face_surface_world(f) else {
        return false;
    };

    if the_edge.is_none() {
        let toproj = BRepTool::vertex_point_world(v);
        // `cxx:506`: `TopExp_Explorer(F.Oriented(TopAbs_FORWARD), TopAbs_EDGE)`.
        for edg in edges_of(&f.0.clone().oriented(Orientation::Forward)) {
            let (f_par, l_par) = BRepTool::edge_parameters(&edg);
            let curve = BRepTool::edge_curve_world(&edg);
            let mut a_cur_dist = INFINITE;
            let mut a_cur_par = INFINITE;
            if let Some(c) = curve.as_ref() {
                // `cxx:515-523`.
                a_cur_par = project_vertex_on_edge(v, &edg);
                if precision::Precision::is_infinite(a_cur_par) {
                    continue;
                }
                let a_cur_p_bound = c.d0(a_cur_par);
                a_cur_dist = a_cur_p_bound.square_distance(&toproj);
            } else if !precision::Precision::is_infinite(p2d.x()) {
                // `cxx:526-547`: no 3D curve, project in the face's 2d space.
                let Some(_c2d) = boptools_2d::curve_on_surface(&edg, f) else {
                    continue;
                };
                a_cur_par = project_vertex_2d_on_edge(p2d, &edg, f);
                if precision::Precision::is_infinite(a_cur_par) {
                    continue;
                }
                let Some((pc, _cf, _cl)) = boptools_2d::curve_on_surface_range(&edg, f) else {
                    continue;
                };
                let a_p_proj = pc.d0(a_cur_par);
                let a_cur_p_bound = a_surf.d0(a_p_proj.x(), a_p_proj.y());
                a_cur_dist = a_cur_p_bound.square_distance(&toproj);
            }
            // `f_par` / `l_par` are read by the reference and never used past
            // this point (`cxx:510`); the range lookup itself has no effect.
            let _ = (f_par, l_par);

            // `cxx:549-554`.
            if a_cur_dist < dmin {
                *the_edge = Some(edg.clone());
                dmin = a_cur_dist;
                *param = a_cur_par;
            }
        }
        if the_edge.is_none() {
            return false;
        }
    } else if precision::Precision::is_infinite(*param) {
        // `cxx:562-566`.
        let edge = the_edge.as_ref().expect("edge checked above");
        let curve = BRepTool::edge_curve_world(edge);
        *param = match curve {
            Some(_) => project_vertex_on_edge(v, edge),
            None => project_vertex_2d_on_edge(p2d, edge, f),
        };
    }

    let edge = the_edge.as_ref().expect("edge checked above");
    let ttol = a_tol_v + BRepTool::edge_tolerance(edge);
    // `cxx:568-569`.
    if dmin <= ttol * ttol {
        valret = true;
        // `cxx:572-573`: `GeomAdaptor_Surface adSurf(aSurf)` over the surface's
        // own bounds.
        let ad_surf = SurfaceAdaptor::of_surface(a_surf.clone());
        let an_u_resolution = ad_surf.u_resolution(1.0);
        // OCCT reads `UResolution(1.)` for this one too (`cxx:576`); the copy is
        // in the reference and is kept here.
        let a_v_resolution = ad_surf.u_resolution(1.0);

        // `cxx:578-579`.
        if let Some((a_crv_bound, _cf, _cl)) = boptools_2d::curve_on_surface_range(edge, f) {
            let mut a_p_bound2d = a_crv_bound.d0(*param);

            // `cxx:586-587`.
            let mut dumax =
                0.01 * (ad_surf.last_u_parameter() - ad_surf.first_u_parameter());
            let mut dvmax =
                0.01 * (ad_surf.last_v_parameter() - ad_surf.first_v_parameter());

            let mut dumin = (p2d.x() - a_p_bound2d.x()).abs();
            let mut dvmin = (p2d.y() - a_p_bound2d.y()).abs();

            // `cxx:592-610`: the U period shift.
            if dumin > dumax && a_surf.is_u_periodic() {
                let mut x1 = a_p_bound2d.x();
                let a_shift = crate::pcurve_full::adjust_to_period(
                    x1,
                    ad_surf.first_u_parameter(),
                    ad_surf.last_u_parameter(),
                );
                x1 += a_shift;
                a_p_bound2d.set_x(x1);
                let mut x2 = p2d.x();
                let a_shift = crate::pcurve_full::adjust_to_period(
                    x2,
                    ad_surf.first_u_parameter(),
                    ad_surf.last_u_parameter(),
                );
                x2 += a_shift;
                dumin = (x2 - x1).abs();
                if dumin > dumax
                    && (dumin - a_surf.u_period()).abs() < PCONFUSION
                {
                    // `cxx:605-608`: `aX2 = aX1; dumin = 0.;`
                    dumin = 0.0;
                }
            }

            // `cxx:612-631`: the V period shift. The period comparison uses
            // `Precision::Confusion()` here (not `PConfusion`), as in the
            // reference (`cxx:625`).
            if dvmin > dvmax && a_surf.is_v_periodic() {
                let mut y1 = a_p_bound2d.y();
                let a_shift = crate::pcurve_full::adjust_to_period(
                    y1,
                    ad_surf.first_v_parameter(),
                    ad_surf.last_v_parameter(),
                );
                y1 += a_shift;
                a_p_bound2d.set_y(y1);
                let mut y2 = p2d.y();
                let a_shift = crate::pcurve_full::adjust_to_period(
                    y2,
                    ad_surf.first_v_parameter(),
                    ad_surf.last_v_parameter(),
                );
                y2 += a_shift;
                dvmin = (y1 - y2).abs();
                if dvmin > dvmax && (dvmin - a_surf.v_period()).abs() < CONFUSION {
                    // `cxx:626-629`.
                    dvmin = 0.0;
                }
            }

            // `cxx:632-648`.
            let mut a_dist3d = a_tol_v;
            if dumin > dumax || dvmin > dvmax {
                dumax = ad_surf.u_resolution(a_tol_v);
                dvmax = ad_surf.v_resolution(a_tol_v);
                let a_tol2d = 2.0 * dumax.max(dvmax);
                let a_dist2d = dumin.max(dvmin);
                if a_dist2d > a_tol2d {
                    let a_dist3d1 =
                        a_dist2d / an_u_resolution.max(a_v_resolution);
                    if a_dist3d1 > a_dist3d {
                        a_dist3d = a_dist3d1;
                    }
                }
            }

            // `cxx:651-663`: the 3d check, the same one
            // `BRepCheck_Wire::SelfIntersect` makes.
            let a_p_bound = a_surf.d0(a_p_bound2d.x(), a_p_bound2d.y());
            let a_p_v2d = a_surf.d0(p2d.x(), p2d.y());
            let a_dist_points_3d = a_p_v2d.square_distance(&a_p_bound);
            let a_max_dist = a_dist_points_3d.max(a_dist3d * a_dist3d);
            if a_tol_v * a_tol_v < a_max_dist {
                // `cxx:659-662`: `B.UpdateVertex(V, sqrt(aMaxDist))`.
                v.set_tolerance(a_max_dist.sqrt());
            }
        }
    }
    valret
}

/// `Project(V, theEdge)` (`LocPoe_WiresOnShape.cxx:679-698`).
///
/// `GeomAPI_ProjectPointOnCurve::Init(P, C, f, l)` followed by
/// `NbPoints() > 0 ? LowerDistanceParameter() : Precision::Infinite()`.
pub(super) fn project_vertex_on_edge(v: &Vertex, edge: &Edge) -> f64 {
    let toproj = BRepTool::vertex_point_world(v);
    let Some(c) = BRepTool::edge_curve_world(edge) else {
        return INFINITE;
    };
    let (first, last) = BRepTool::edge_parameters(edge);
    // `cxx:689-691`: the location is folded into the curve, which
    // `edge_curve_world` already did.
    match occt_geom::extrema_pc::point_curve_extrema_range(c.as_ref(), &toproj, first, last) {
        Some(pair) => pair.u1,
        None => INFINITE,
    }
}

/// `Project(V, p2d, theEdge, theFace)` (`cxx:702-718`).
pub(super) fn project_vertex_2d_on_edge(p2d: &GpPnt2d, edge: &Edge, face: &Face) -> f64 {
    // `cxx:714`: `PC = BRep_Tool::CurveOnSurface(theEdge, theFace, f, l)`.
    let Some((pc, first, last)) = boptools_2d::curve_on_surface_range(edge, face) else {
        return INFINITE;
    };
    // `cxx:716`: `proj.Init(p2d, PC, f, l)`.
    let proj = ProjectPointOnCurve2d::with_curve_range(p2d, &pc, first, last);
    // `cxx:718`: `proj.NbPoints() > 0 ? proj.LowerDistanceParameter() : Infinite`.
    proj.lower_distance_parameter().unwrap_or(INFINITE)
}

/// Unused import guard: `GpVec2d` / `explorer_shapes` / `ShapeType` are named
/// because the transposed reference reads them; keep the module warning-free.
#[allow(dead_code)]
fn _reference_markers() {
    let _ = GpVec2d::zero();
    let _ = explorer_shapes as fn(&crate::shape::TopoShape, ShapeType) -> Vec<crate::shape::TopoShape>;
}
