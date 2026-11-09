//! The 2D direction helpers of `LocOpe_SplitShape`:
//! `SameUV` (`cxx:75-96`), `GetDirection` (`cxx:1595-1629`), `ChoixUV`
//! (`cxx:1631-1689`), `ChooseDirection` (`cxx:1691-1776`) and
//! `checkOverlapping` (`cxx:735-777`).
//!
//! All five work in the face's parametric space: they read the p-curve of an
//! edge on a face, take its tangent at one end, and compare angles so the wire
//! walk continues along the branch that does not fold back over itself.

use std::f64::consts::PI;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_approx::surface_tool::SurfaceTool;
use crate::brep_tool::BRepTool;
use crate::locope::shape_maps::ShapeIndexedMap;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `SameUV(P1, P2, theBAS)` (`cxx:75-96`): equality modulo whole periods.
///
/// The reference starts from `isSame = true` and only ever *ands* the periodic
/// axes in, so a surface periodic in neither axis answers `true`.
pub fn same_uv(p1: &GpPnt2d, p2: &GpPnt2d, bas: &SurfaceAdaptor) -> bool {
    let mut is_same = true;
    // `cxx:81-84`.
    if SurfaceTool::is_u_periodic(bas) {
        is_same = (p1.x() - p2.x()).abs() < SurfaceTool::u_period(bas) * 0.5;
    }
    // `cxx:85-88`.
    if SurfaceTool::is_v_periodic(bas) {
        is_same = is_same && (p1.y() - p2.y()).abs() < SurfaceTool::v_period(bas) * 0.5;
    }
    is_same
}

/// `GetDirection(theEdge, theFace, thePnt, theDir, isFirstEnd)`
/// (`cxx:1595-1629`).
///
/// Takes the p-curve tangent at the requested end of the edge, one percent of
/// the range inwards, and reverses it for a `REVERSED` edge.
pub fn get_direction(e: &Edge, f: &Face, is_first_end: bool) -> Option<(GpPnt2d, GpVec2d)> {
    let (c2d, a_first, a_last) = boptools_2d::curve_on_surface_range(e, f)?;
    let an_or = e.0.orientation();
    // `cxx:1604-1607`.
    let take_first = (an_or == Orientation::Forward && is_first_end)
        || (an_or == Orientation::Reversed && !is_first_end);
    let dpar = (a_last - a_first) * 0.01;

    // `cxx:1609-1626`.
    let (the_pnt, mut the_dir) = if take_first {
        let the_pnt = c2d.d0(a_first);
        let a_next_pnt = c2d.d0(a_first + dpar);
        (the_pnt, vec_between(&the_pnt, &a_next_pnt))
    } else {
        let the_pnt = c2d.d0(a_last);
        let a_prev_pnt = c2d.d0(a_last - dpar);
        (the_pnt, vec_between(&a_prev_pnt, &the_pnt))
    };
    if an_or == Orientation::Reversed {
        the_dir.reverse();
    }
    Some((the_pnt, the_dir))
}

/// `gp_Vec2d(P1, P2)`: the vector `P2 - P1`.
fn vec_between(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}

/// `ChoixUV(Last, F, Poss, theResEdge, plst, dlst)` (`cxx:1631-1689`).
///
/// Among the candidate edges `poss` starting at `plst`, takes the one whose
/// initial tangent is the largest signed angle away from `dlst` - the reference
/// measures `gp_Dir2d::Angle`, so the choice is the extreme left turn.
///
/// `plst` / `dlst` are updated in place to the chosen edge's other end
/// (`cxx:1680-1684`). Returns `false` when no candidate starts at `plst`.
pub fn choix_uv(
    last: &Edge,
    f: &Face,
    poss: &ShapeIndexedMap<TopoShape>,
    pgp_plst: &mut GpPnt2d,
    pgp_dlst: &mut GpVec2d,
) -> Option<(Edge, GpPnt2d, GpVec2d)> {
    let surf = SurfaceAdaptor::of_face(f);
    let ref2d = GpDir2d::from_vec2d(pgp_dlst).ok()?;

    // `cxx:1649-1678`.
    let mut index = 1usize;
    let mut imin = 0usize;
    let mut angmax = -PI;
    while index <= poss.extent() {
        let an_edge = Edge::wrap(poss.find_key(index).clone())?;
        let (p2d, v2d) = get_direction(&an_edge, f, true)?;
        // `cxx:1658-1661`.
        if !same_uv(pgp_plst, &p2d, &surf) {
            index += 1;
            continue;
        }
        // `cxx:1662-1663`: `surf.D0(...)` (the 3D point is unused).
        let _ = surf.d0(p2d.x(), p2d.y());

        // `cxx:1665-1672`.
        let ang = if !is_same(&last.0, &an_edge.0) {
            ref2d.angle(&GpDir2d::from_vec2d(&v2d).ok()?)
        } else {
            -PI
        };
        // `cxx:1674-1677`.
        if ang > angmax {
            imin = index;
            angmax = ang;
        }
        index += 1;
    }

    // `cxx:1680-1687`.
    if imin == 0 {
        return None;
    }
    let the_res_edge = Edge::wrap(poss.find_key(imin).clone())?;
    let (plst, dlst) = get_direction(&the_res_edge, f, false)?;
    *pgp_plst = plst;
    *pgp_dlst = dlst;
    Some((the_res_edge, plst, dlst))
}

/// `ChooseDirection(RefDir, RefVertex, theFace, Ldirs)` (`cxx:1691-1776`).
///
/// Picks the shape of `ldirs` whose tangent at `ref_vertex` leaves `ref_dir`
/// with the smallest counter-clockwise angle. `ref_dir` is a shape (an edge,
/// usually a wire) whose *first* matching edge is the reference.
pub fn choose_direction(
    ref_dir: &TopoShape,
    ref_vertex: &crate::shape::Vertex,
    the_face: &Face,
    ldirs: &[TopoShape],
) -> Option<TopoShape> {
    // `cxx:1697-1712`: the reference edge of `RefDir` and the end it meets
    // `RefVertex` at.
    let mut ref_edge: Option<Edge> = None;
    let mut an_or = Orientation::Forward;
    for sub in crate::topo_tools::explorer_shapes(ref_dir, ShapeType::Edge) {
        let Some(e) = Edge::wrap(sub) else {
            continue;
        };
        let (v1, v2) = crate::topo_tools::edge_vertices(&e);
        if let Some(v1) = &v1 {
            if is_same(&v1.0, &ref_vertex.0) {
                an_or = Orientation::Reversed;
                ref_edge = Some(e);
                break;
            }
        }
        if let Some(v2) = &v2 {
            if is_same(&v2.0, &ref_vertex.0) {
                an_or = Orientation::Forward;
                ref_edge = Some(e);
                break;
            }
        }
        ref_edge = Some(e);
    }
    let ref_edge = ref_edge?;

    // `cxx:1715-1728`.
    let (ref_curve, ref_first, ref_last) = boptools_2d::curve_on_surface_range(&ref_edge, the_face)?;
    let ref_par = if an_or == Orientation::Forward { ref_last } else { ref_first };
    let (_, mut ref_vec) = ref_curve.d1(ref_par);
    if an_or == Orientation::Forward {
        ref_vec.reverse();
    }

    // `cxx:1730-1773`.
    let mut min_angle = f64::MAX; // `RealLast()`
    let mut target_dir: Option<TopoShape> = None;
    for a_shape in ldirs {
        let mut an_edge: Option<Edge> = None;
        for sub in crate::topo_tools::explorer_shapes(a_shape, ShapeType::Edge) {
            let Some(e) = Edge::wrap(sub) else {
                continue;
            };
            let (v1, v2) = crate::topo_tools::edge_vertices(&e);
            if let Some(v1) = &v1 {
                if is_same(&v1.0, &ref_vertex.0) {
                    an_or = Orientation::Forward;
                    an_edge = Some(e);
                    break;
                }
            }
            if let Some(v2) = &v2 {
                if is_same(&v2.0, &ref_vertex.0) {
                    an_or = Orientation::Reversed;
                    an_edge = Some(e);
                    break;
                }
            }
            an_edge = Some(e);
        }
        let Some(an_edge) = an_edge else {
            continue;
        };
        let Some((a_curve, a_first, a_last)) = boptools_2d::curve_on_surface_range(&an_edge, the_face)
        else {
            continue;
        };
        let a_par = if an_or == Orientation::Forward { a_first } else { a_last };
        let (_, mut a_vec) = a_curve.d1(a_par);
        if an_or == Orientation::Reversed {
            a_vec.reverse();
        }
        let mut an_angle = a_vec.angle(&ref_vec);
        if an_angle < 0.0 {
            an_angle += 2.0 * PI;
        }
        if an_angle < min_angle {
            min_angle = an_angle;
            target_dir = Some(a_shape.clone());
        }
    }
    target_dir
}

/// `checkOverlapping(theEdge1, theEdge2, theFace)` (`cxx:735-777`).
///
/// Samples the second edge at three interior points and demands that each one
/// projects onto the first within the 2D resolution of the edge tolerances and
/// with a parallel tangent. `true` therefore means "the two edges cover the
/// same stretch of p-curve", which makes the split degenerate.
pub fn check_overlapping(e1: &Edge, e2: &Edge, f: &Face) -> bool {
    let an_ad_s = SurfaceAdaptor::of_face(f);
    // `cxx:742`.
    let max_tol = BRepTool::edge_tolerance(e1) + BRepTool::edge_tolerance(e2);
    // `cxx:744-745`.
    let a_max_tol2d = an_ad_s.u_resolution(max_tol).max(an_ad_s.v_resolution(max_tol));
    let a_tol_ang = PI / 180.0;

    // `cxx:749-751`.
    let (Some((crv1, a_f1, a_l1)), Some((crv2, a_f2, a_l2))) = (
        boptools_2d::curve_on_surface_range(e1, f),
        boptools_2d::curve_on_surface_range(e2, f),
    ) else {
        return false;
    };

    // `cxx:752-773`.
    let nb_p = 4usize;
    let a_step = (a_l2 - a_f2) / nb_p as f64;
    for j in 1..nb_p {
        let par2 = a_f2 + a_step * j as f64;
        let (a_p2d, a_v2) = crv2.d1(par2);

        let mut proj = ProjectPointOnCurve2d::new();
        proj.init_range(&a_p2d, &crv1, a_f1, a_l1);
        // `cxx:762-765`.
        if proj.nb_points() == 0 || proj.lower_distance().unwrap_or(f64::MAX) > a_max_tol2d {
            return false;
        }
        // `cxx:766-771`.
        let Some(par1) = proj.lower_distance_parameter() else {
            return false;
        };
        let (_, a_v1) = crv1.d1(par1);
        if !a_v1.is_parallel(&a_v2, a_tol_ang) {
            return false;
        }
    }
    true
}
