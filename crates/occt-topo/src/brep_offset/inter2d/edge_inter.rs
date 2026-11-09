//! `EdgeInter` (`BRepOffset_Inter2d.cxx:416-699`): the 2D intersection of the
//! pcurves of two edges that lie on one face, producing the intersection
//! vertices in `theAsDes`.
//!
//! Control flow, in the reference's order:
//! 1. `E1.IsSame(E2)` returns at once (`cxx:429-432`).
//! 2. The ranges start as `BRep_Tool::Range(E1/E2)` (`cxx:437-438`) and are
//!    overwritten by `CurveOnSurface` (`cxx:481-482`) *only* on the no-common-
//!    vertex arm; the endpoint tests below therefore see the COS range exactly
//!    when the pcurves were intersected.
//! 3. `BRepLib::BuildCurve3d` on both edges, then `Geom2dInt_GInter` on the two
//!    pcurves with `TolDub = 1.e-7`; a degenerated edge contributes one point
//!    (the degenerated vertex's) instead of the pcurve point lifted to the
//!    surface (`cxx:459-500`).
//! 4. Each solution becomes a vertex at `P` with `INTERNAL` orientation, its
//!    parameter written on both edges (`cxx:516-518`), its tolerance raised to
//!    the largest of `|P1-P|`, `|P2-P|`, `|P1-P2|` (`cxx:519-527`), and an
//!    orientation - `REVERSED` unless the pcurve tangents say otherwise
//!    (`cxx:551-588`).
//! 5. The four endpoint pairs are tested for coincidence (`cxx:594-641`): a
//!    coincident pair prepends a vertex carrying the *edge* vertex orientation.
//! 6. Coincident solutions are purged repeatedly (`cxx:650-690`), then the two
//!    vertex lists are stored (`cxx:696-698`).
//!
//! Port guards (no OCCT counterpart): a missing pcurve stops the function
//! (OCCT builds a `Geom2dAdaptor_Curve` on a null handle and faults), and a
//! pcurve that cannot be evaluated skips the one solution or the one
//! orientation.

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::Curve2d;

use crate::abs::Orientation;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_algo::AsDes;
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, is_same, vertex_position};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::store_vertices::store_on_two_edges;
use super::{pcurve_value_on_face, Dmvv};

/// `Geom2dInt_GInter`'s `TolDub` (`cxx:434`, `:732`).
pub(crate) const TOL_DUB: f64 = 1.0e-7;

/// `BRepLib::BuildCurve3d(E)`'s defaults (`BRepLib.hxx:...`): tolerance
/// `1.0e-5`, continuity `GeomAbs_C1`, max degree 14, max segment 0.
pub(crate) fn build_curve3d_defaults(e: &Edge) {
    let _ = build_curve3d(e, 1.0e-5, Shape::C1, 14, 0);
}

/// The `DegPoint` half of the `WithDegen` block (`cxx:459-477`, `:757-775`):
/// `Some(P)` when one of the edges is degenerated, `None` otherwise.
pub(crate) fn degenerate_point(
    basurf: &SurfaceAdaptor,
    f: &Face,
    e1: &Edge,
    e2: &Edge,
) -> Option<GpPnt> {
    let d1 = BRepTool::is_degenerated(e1);
    let d2 = BRepTool::is_degenerated(e2);
    if !d1 && !d2 {
        return None;
    }
    let deg = if d1 { e1 } else { e2 };
    // `TopoDS_Iterator iter(EI[ideg]); if (iter.More()) DegPoint = Pnt(vdeg);`
    let (v1, v2) = edge_vertices(deg);
    if let Some(v) = v1.or(v2) {
        return Some(vertex_position(&v));
    }
    // `BRepAdaptor_Curve CEdeg(EI[ideg], F); Value(FirstParameter())`.
    let (pc, fpar, _lpar) = crate::boptools_2d::curve_on_surface_range(deg, f)?;
    pcurve_value_on_face(basurf, pc.as_ref(), fpar, deg)
}

/// `BRepAdaptor_Curve(E1, F)` / `BRepAdaptor_Curve(E2, F)` at one parameter,
/// through the pcurve on the face (`cxx:519-520`, `:839-840`).
pub(crate) fn curve_on_face_value(
    basurf: &SurfaceAdaptor,
    pc: &dyn Curve2d,
    u: f64,
    e: &Edge,
) -> Option<GpPnt> {
    pcurve_value_on_face(basurf, pc, u, e)
}

/// The `ResPoints` / `ResParamsOnE1` / `ResParamsOnE2` fill of
/// `cxx:486-500` and `cxx:805-820`: one entry per intersection point, the 3D
/// point being the degenerated vertex's when the edge is degenerated and the
/// surface value of the pcurve point otherwise.
pub(crate) fn collect_intersection_points(
    basurf: &SurfaceAdaptor,
    inter: &Geom2dIntGInter,
    deg_point: Option<GpPnt>,
) -> (Vec<GpPnt>, Vec<f64>, Vec<f64>) {
    let mut points = Vec::new();
    let mut params_e1 = Vec::new();
    let mut params_e2 = Vec::new();
    for i in 1..=inter.nb_points() {
        let ip = inter.point(i);
        let p3d = match deg_point {
            Some(p) => p,
            None => {
                let p2d = ip.value();
                // Port guard: `BAsurf.Value` cannot fail in OCCT; a surface
                // that answers nothing drops this solution.
                let Some(p) = basurf.value(p2d.x(), p2d.y()) else {
                    continue;
                };
                p
            }
        };
        points.push(p3d);
        params_e1.push(ip.param_on_first());
        params_e2.push(ip.param_on_second());
    }
    (points, params_e1, params_e2)
}

/// One solution of `cxx:504-588` / `cxx:822-918`: the `INTERNAL` vertex with
/// its parameters on both edges, its raised tolerance and the two orientations
/// the caller stores it with.
pub(crate) struct NewIntersectionVertex {
    pub vertex: TopoShape,
    pub or1: Orientation,
    pub or2: Orientation,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn make_intersection_vertex(
    f: &Face,
    basurf: &SurfaceAdaptor,
    e1: &Edge,
    e2: &Edge,
    pc1: &dyn Curve2d,
    pc2: &dyn Curve2d,
    t1: f64,
    t2: f64,
    p: &GpPnt,
    tol: f64,
    with_ori: bool,
) -> Option<NewIntersectionVertex> {
    // `BRepLib_MakeVertex(P)` (`BRepLib_MakeVertex.cxx:26-31`) uses
    // `BRepLib::Precision()` = `Precision::Confusion()`.
    let mut v = TopoBuilder::new().make_vertex(*p, CONFUSION);
    v.0.set_orientation(Orientation::Internal);

    // `cxx:517-518`.
    let _ = update_vertex_on_edge(&v.0, t1, &e1.0, tol);
    let _ = update_vertex_on_edge(&v.0, t2, &e2.0, tol);

    // `cxx:519-527`: `dist1 = max(|P1-P|, |P2-P|, |P1-P2|)` then
    // `B.UpdateVertex(aNewVertex, dist1)` (raise only).
    let p1 = curve_on_face_value(basurf, pc1, t1, e1);
    let p2 = curve_on_face_value(basurf, pc2, t2, e2);
    if let (Some(p1), Some(p2)) = (p1, p2) {
        let dist = p1
            .distance(p)
            .max(p2.distance(p))
            .max(p1.distance(&p2));
        GeometryRegistry::global().raise_vertex_tolerance(&v.0, dist);
    }

    // `cxx:551-588`.
    let mut or1 = Orientation::Reversed;
    let mut or2 = Orientation::Reversed;
    if with_ori {
        let pce1 = Curve2dAdaptor::of_edge_face(e1, f);
        let pce2 = Curve2dAdaptor::of_edge_face(e2, f);
        // Port guard: `PCE1.D1` faults on a pcurve-less edge.
        let (_, v1) = pce1.d1(t1)?;
        let (_, v2) = pce2.d1(t2)?;
        let mut v1or = v1;
        let mut v2or = v2;
        if e1.0.orientation() == Orientation::Reversed {
            v1or.reverse();
        }
        if e2.0.orientation() == Orientation::Reversed {
            v2or.reverse();
        }
        if v2or.crossed(&v1) > 0.0 {
            or1 = Orientation::Forward;
        }
        if v1or.crossed(&v2) > 0.0 {
            or2 = Orientation::Forward;
        }
    }

    Some(NewIntersectionVertex {
        vertex: v.0,
        or1,
        or2,
    })
}

/// The endpoint-coincidence block of `cxx:594-641` and `cxx:925-967`.
///
/// `aTol` in `EdgeInter` is `max(VertexTolerance(V1[j]), VertexTolerance(V2[k]))`
/// while `RefEdgeInter` passes `Tol`; `use_vertex_tol` selects the arm.
///
/// A pair is skipped when the two vertices are `IsSame` *and* the vertex
/// already has an ascendant in `theAsDes` (`cxx:617-621`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn prepend_endpoint_vertices(
    e1: &Edge,
    e2: &Edge,
    f1: f64,
    l1: f64,
    f2: f64,
    l2: f64,
    tol: f64,
    use_vertex_tol: bool,
    as_des: &AsDes,
    lv1: &mut Vec<TopoShape>,
    lv2: &mut Vec<TopoShape>,
) {
    // `TopExp::Vertices(E1, V1[0], V1[1])` / `(E2, ...)` (`cxx:597-598`,
    // `:927-928`) use the default `CumOri = false`, i.e. the stored vertex
    // order, *not* the orientation-cumulated `CommonVertex` order.
    let (v1a, v1b) = edge_vertices(e1);
    let (v2a, v2b) = edge_vertices(e2);
    let v1s = [v1a, v1b];
    let v2s = [v2a, v2b];
    for j in 0..2 {
        let Some(v1j) = v1s[j].as_ref() else {
            continue;
        };
        for k in 0..2 {
            let Some(v2k) = v2s[k].as_ref() else {
                continue;
            };
            if is_same(&v1j.0, &v2k.0) && as_des.has_ascendant(&v1j.0) {
                continue;
            }
            let p1 = vertex_position(v1j);
            let p2 = vertex_position(v2k);
            let dist = p1.distance(&p2);
            if dist >= tol {
                continue;
            }
            let a_tol = if use_vertex_tol {
                BRepTool::vertex_tolerance(v1j).max(BRepTool::vertex_tolerance(v2k))
            } else {
                tol
            };
            let mut v = TopoBuilder::new().make_vertex(p1, CONFUSION);
            let u1 = if j == 0 { f1 } else { l1 };
            let u2 = if k == 0 { f2 } else { l2 };
            v.0.set_orientation(Orientation::Internal);
            // `B.UpdateVertex(aLocalShape, U1, E1, aTol)`.
            let _ = update_vertex_on_edge(&v.0, u1, &e1.0, a_tol);
            let _ = update_vertex_on_edge(&v.0, u2, &e2.0, a_tol);
            // `LV1.Prepend(V.Oriented(V1[j].Orientation()))`.
            lv1.insert(0, v.0.oriented(v1j.0.orientation()));
            lv2.insert(0, v.0.oriented(v2k.0.orientation()));
        }
    }
}

/// The "remove all doubles" loop of `cxx:650-690` (`EdgeInter`) and
/// `cxx:970-1013` (`RefEdgeInter`).
///
/// The reference rescans the whole list from the start after every removal
/// (`while (Purge)`), which the port reproduces with the `'rescan` loop. The
/// tolerance is `max(tol1, tol2)` in `EdgeInter` and the flat `Tol` in
/// `RefEdgeInter`; `use_vertex_tol` selects the arm.
pub(crate) fn purge_doubles(
    lv1: &mut Vec<TopoShape>,
    lv2: &mut Vec<TopoShape>,
    tol: f64,
    use_vertex_tol: bool,
) {
    loop {
        let mut purged = false;
        let mut i = 0;
        'rescan: while i < lv1.len() {
            let p1 = vertex_position(&Vertex(lv1[i].clone()));
            let a_tol1 = BRepTool::vertex_tolerance(&Vertex(lv1[i].clone()));
            let mut j = 0;
            while j < i {
                let p2 = vertex_position(&Vertex(lv1[j].clone()));
                let a_tol = if use_vertex_tol {
                    a_tol1.max(BRepTool::vertex_tolerance(&Vertex(lv1[j].clone())))
                } else {
                    tol
                };
                if p1.distance(&p2) <= a_tol {
                    lv1.remove(i);
                    lv2.remove(i);
                    purged = true;
                    break 'rescan;
                }
                j += 1;
            }
            i += 1;
        }
        if !purged {
            break;
        }
    }
}

/// `EdgeInter(F, BAsurf, E1, E2, AsDes, Tol, WithOri, aDMVV)` (`cxx:416-699`).
pub(crate) fn edge_inter(
    f: &Face,
    basurf: &SurfaceAdaptor,
    e1: &Edge,
    e2: &Edge,
    as_des: &mut AsDes,
    tol: f64,
    with_ori: bool,
    dmvv: &mut Dmvv,
) {
    // `cxx:429-432`.
    if is_same(&e1.0, &e2.0) {
        return;
    }

    // `cxx:437-438`: `BRep_Tool::Range`.
    let (rf1, rl1) = BRepTool::edge_parameters(e1);
    let (rf2, rl2) = BRepTool::edge_parameters(e2);
    let mut fpar = [0.0, rf1, rf2];
    let mut lpar = [0.0, rl1, rl2];

    let mut lv1: Vec<TopoShape> = Vec::new();
    let mut lv2: Vec<TopoShape> = Vec::new();

    // `cxx:449-452`.
    if super::orient_edges::common_vertex(e1, e2).is_none() {
        // `cxx:453-454`.
        build_curve3d_defaults(e1);
        build_curve3d_defaults(e2);

        // `cxx:456-457`: `TolSum` is computed and never read, so it is unused
        // in the port as well.
        let deg_point = degenerate_point(basurf, f, e1, e2);

        // `cxx:481-482`.
        let (Some((pc1, cf1, cl1)), Some((pc2, cf2, cl2))) = (
            crate::boptools_2d::curve_on_surface_range(e1, f),
            crate::boptools_2d::curve_on_surface_range(e2, f),
        ) else {
            // Port guard: `Geom2dAdaptor_Curve` on a null pcurve raises.
            return;
        };
        fpar[1] = cf1;
        lpar[1] = cl1;
        fpar[2] = cf2;
        lpar[2] = cl2;

        // `cxx:483-486`: `Geom2dInt_GInter(GAC1, GAC2, TolDub, TolDub)`, the
        // domains being each adaptor's own `[f, l]`.
        let mut inter = Geom2dIntGInter::new();
        let d1 = Geom2dIntGInter::domain_over_range(pc1.as_ref(), fpar[1], lpar[1], TOL_DUB);
        let d2 = Geom2dIntGInter::domain_over_range(pc2.as_ref(), fpar[2], lpar[2], TOL_DUB);
        inter.perform(pc1.as_ref(), &d1, pc2.as_ref(), &d2, TOL_DUB, TOL_DUB);

        // `cxx:487-500`.
        let (res_points, res_params_e1, res_params_e2) =
            collect_intersection_points(basurf, &inter, deg_point);

        // `cxx:502-589`.
        for i in 0..res_points.len() {
            let t1 = res_params_e1[i];
            let t2 = res_params_e2[i];
            // `Precision::IsInfinite(aT1) || Precision::IsInfinite(aT2)`.
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
            lv1.push(nv.vertex.oriented(nv.or1));
            lv2.push(nv.vertex.oriented(nv.or2));
        }
    }

    // `cxx:594-641`.
    prepend_endpoint_vertices(
        e1,
        e2,
        fpar[1],
        lpar[1],
        fpar[2],
        lpar[2],
        tol,
        true,
        as_des,
        &mut lv1,
        &mut lv2,
    );

    // `cxx:643-698`.
    if !lv1.is_empty() {
        // `cxx:650-690`.
        purge_doubles(&mut lv1, &mut lv2, tol, true);
        // `cxx:694-698`: `TolStore = Tolerance(E1) + Tolerance(E2); TolStore =
        // std::max(TolStore, Tol)`.
        let tol_store = (BRepTool::edge_tolerance(e1) + BRepTool::edge_tolerance(e2)).max(tol);
        store_on_two_edges(&e1.0, &e2.0, &lv1, &lv2, tol_store, as_des, dmvv);
    }
}

