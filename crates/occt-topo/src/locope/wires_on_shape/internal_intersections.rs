//! File-local `FindInternalIntersections(theEdge, theFace, Splits, isOverlapped)`
//! (`LocOpe_WiresOnShape.cxx:1291-1511`).
//!
//! Finds the parameters where the splitting edge meets another edge of the face
//! in the middle of both (i.e. away from the two endpoint tolerance discs), and
//! records the pieces of `theEdge` those intersections split it into.

use occt_core::bnd::BndBox2d;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::Orientation;
use crate::bnd_lib_add2d::add_adaptor_range;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeIndexedMap;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, edges_of};

/// `TrimmedSquareDistances` of `Extrema_ExtCC` (`Extrema_ExtCC.cxx:190-245`,
/// `:377-395`): the squared distances between the four trimmed endpoints, in the
/// `(P1f-P2f, P1f-P2l, P1l-P2f, P1l-P2l)` order.
fn trimmed_square_distances(
    c1: &dyn occt_geom::curve::Curve,
    par1: (f64, f64),
    c2: &dyn occt_geom::curve::Curve,
    par2: (f64, f64),
) -> [f64; 4] {
    let p1f = c1.d0(par1.0);
    let p1l = c1.d0(par1.1);
    let p2f = c2.d0(par2.0);
    let p2l = c2.d0(par2.1);
    [
        p1f.square_distance(&p2f),
        p1f.square_distance(&p2l),
        p1l.square_distance(&p2f),
        p1l.square_distance(&p2l),
    ]
}

/// `FindInternalIntersections` (`LocOpe_WiresOnShape.cxx:1291-1511`).
pub(super) fn find_internal_intersections(
    the_edge: &Edge,
    the_face: &Face,
    splits: &mut ShapeIndexedMap<Vec<TopoShape>>,
    is_overlapped: &mut bool,
) {
    const TOL_EXT: f64 = PCONFUSION;
    // `cxx:1303`: `BRepAdaptor_Surface anAdSurf(theFace, false)` is built and
    // never read in the reference body, so it is not reproduced here.
    let mut split_pars: Vec<f64> = Vec::new();

    // `cxx:1306-1312`.
    let (v0, v1) = edge_vertices(the_edge);
    let (Some(v0), Some(v1)) = (v0, v1) else {
        return;
    };
    let the_pnt = [
        BRepTool::vertex_point_world(&v0),
        BRepTool::vertex_point_world(&v1),
    ];
    let a_tol_v = [
        BRepTool::vertex_tolerance(&v0),
        BRepTool::vertex_tolerance(&v1),
    ];
    // `cxx:1314-1318`: "= 4 * 4 - to avoid creating microedges, area around
    // vertices is increased up to 4 vertex tolerance".
    const EXT: f64 = 16.0;
    let a_tol_v_ext = [
        EXT * a_tol_v[0] * a_tol_v[0],
        EXT * a_tol_v[1] * a_tol_v[1],
    ];

    // `cxx:1319-1321`: `BndLib_Add2dCurve::Add(thePCurve, Tolerance(theEdge), theBox)`.
    let the_pcurve = Curve2dAdaptor::of_edge_face(the_edge, the_face);
    let mut the_box = BndBox2d::new();
    if let Some(c) = the_pcurve.curve() {
        add_adaptor_range(
            c.as_ref(),
            the_pcurve.first_parameter(),
            the_pcurve.last_parameter(),
            BRepTool::edge_tolerance(the_edge),
            &mut the_box,
        );
    }

    // `cxx:1323-1330`.
    let the_par = BRepTool::edge_parameters(the_edge);
    let Some(the_curve) = BRepTool::edge_curve(the_edge) else {
        return;
    };
    let a_tol_v2d = [
        the_curve.resolution(a_tol_v[0]).max(PCONFUSION),
        the_curve.resolution(a_tol_v[1]).max(PCONFUSION),
    ];
    let mut a_dist_max = CONFUSION * CONFUSION;

    // `cxx:1331-1419`: every edge of the face.
    for an_edge in edges_of(&the_face.0) {
        // `cxx:1335-1341`.
        let a_pcurve = Curve2dAdaptor::of_edge_face(&an_edge, the_face);
        let mut a_box = BndBox2d::new();
        if let Some(c) = a_pcurve.curve() {
            add_adaptor_range(
                c.as_ref(),
                a_pcurve.first_parameter(),
                a_pcurve.last_parameter(),
                BRepTool::edge_tolerance(&an_edge),
                &mut a_box,
            );
        }
        if the_box.is_out_box(&a_box) {
            continue;
        }

        // `cxx:1343-1347`.
        let (a_fpar, a_lpar) = BRepTool::edge_parameters(&an_edge);
        let Some(a_curve) = BRepTool::edge_curve(&an_edge) else {
            continue;
        };

        // `Extrema_ExtCC(theGAcurve, aGAcurve, TolExt, TolExt)` (`cxx:1347`).
        //
        // `Extrema_ExtCC::Perform` (`Extrema_ExtCC.cxx:176-316`) first tries the
        // analytic `Extrema_ExtElC` arms (line/elementary, two circles) and falls
        // back to `Extrema_ECC` = `Extrema_GGenExtCC`. The repo ports the general
        // engine only (`occt_geom::extrema_cc::GGenExtCC`), so the analytic arm is
        // UNPORTED and every pair here goes through the general engine, which is
        // the same engine `Extrema_ExtCC` uses for its fallback and the same
        // `IsParallel` / `NbExt` / `SquareDistance` / `Points` fields
        // `Extrema_ExtCC::PrepareResults(myECC, ...)` copies
        // (`Extrema_ExtCC.cxx:527-586`).
        let mut an_extrema = occt_geom::extrema_cc::GGenExtCC::new(
            the_curve.as_ref(),
            a_curve.as_ref(),
            the_par,
            (a_fpar, a_lpar),
        );
        // `cxx:180`: `myECC.SetTolerance(std::min(myTol[0], myTol[1]))`.
        an_extrema.set_tolerance(TOL_EXT.min(TOL_EXT));
        an_extrema.perform().ok();

        if !an_extrema.is_done() || an_extrema.nb_ext() == 0 {
            continue;
        }

        let a_max_tol = BRepTool::edge_tolerance(&an_edge);
        let a_max_tol2 = a_max_tol * a_max_tol;
        // `cxx:1355-1359`.
        if an_extrema.is_parallel() && an_extrema.square_distance(1) <= a_max_tol2 {
            *is_overlapped = true;
            return;
        }

        // `cxx:1360-1379`: an intersection at one of the four trimmed endpoints
        // is not "internal".
        let dists = trimmed_square_distances(
            the_curve.as_ref(),
            the_par,
            a_curve.as_ref(),
            (a_fpar, a_lpar),
        );
        for (i, d) in dists.iter().enumerate() {
            let j = if i < 2 { 0 } else { 1 };
            if *d < a_tol_v_ext[j] / EXT {
                return;
            }
        }

        // `cxx:1381-1418`.
        for n in 1..=an_extrema.nb_ext() {
            let a_dist = an_extrema.square_distance(n);
            if a_dist > a_max_tol2 {
                continue;
            }
            let (the_int_par, _, an_int_par, _) = an_extrema.points(n);

            // `cxx:1394-1399`: skip an intersection at one of `theEdge`'s two
            // endpoints.
            let mut at_endpoint = false;
            for j in 0..2 {
                let par_j = if j == 0 { the_par.0 } else { the_par.1 };
                if (the_int_par - par_j).abs() <= a_tol_v2d[j] {
                    at_endpoint = true;
                    break;
                }
            }
            if at_endpoint {
                continue;
            }

            // `cxx:1402-1415`.
            let a_point = a_curve.d0(an_int_par);
            let a_point_int = the_curve.d0(the_int_par);
            if a_point_int.square_distance(&the_pnt[0]) > a_tol_v_ext[0]
                && a_point_int.square_distance(&the_pnt[1]) > a_tol_v_ext[1]
                && a_point.square_distance(&the_pnt[0]) > a_tol_v_ext[0]
                && a_point.square_distance(&the_pnt[1]) > a_tol_v_ext[1]
            {
                split_pars.push(the_int_par);
                if a_dist > a_dist_max {
                    a_dist_max = a_dist;
                }
            }
        }
    }

    // `cxx:1421-1424`.
    if split_pars.is_empty() {
        return;
    }

    // `cxx:1426-1437`: ascending sort.
    let n = split_pars.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if split_pars[i] > split_pars[j] {
                split_pars.swap(i, j);
            }
        }
    }

    // `cxx:1439-1452`: drop repeats (`Sequence::Remove(i + 1)`, 1-based).
    let mut i = 0;
    while i + 1 < split_pars.len() {
        let pnt1 = the_curve.d0(split_pars[i]);
        let pnt2 = the_curve.d0(split_pars[i + 1]);
        if pnt1.square_distance(&pnt2) <= CONFUSION * CONFUSION {
            split_pars.remove(i + 1);
        } else {
            i += 1;
        }
    }

    // `cxx:1454-1511`: build the pieces.
    let mut new_edges: Vec<TopoShape> = Vec::new();
    let builder = TopoBuilder::new();
    let an_orient = the_edge.0.orientation();
    // The pieces keep the edge's representations (`BRep_TEdge::EmptyCopy`,
    // `BRep_TEdge.cxx:105-133`), which for this helper is the 3D curve and the
    // p-curve on `the_face`; both are copied onto each piece below.
    let the_cos = crate::boptools_2d::curve_on_surface_range(the_edge, the_face);

    let mut first_vertex = v0.clone();
    let mut first_par = the_par.0;
    for i in 0..=split_pars.len() {
        first_vertex.0.set_orientation(Orientation::Forward);
        let (last_par, mut last_vertex) = if i < split_pars.len() {
            let last_par = split_pars[i];
            let last_point = the_curve.d0(last_par);
            let lv = builder.make_vertex(last_point, 0.0);
            // `cxx:1469-1470`: `aB.UpdateVertex(LastVertex, sqrt(aDistMax))`.
            lv.set_tolerance(a_dist_max.sqrt());
            (last_par, lv)
        } else {
            (the_par.1, v1.clone())
        };
        last_vertex.0.set_orientation(Orientation::Reversed);

        // `cxx:1478-1486`: `EmptyCopied`, then `Range`, then the two vertices.
        let mut new_edge = builder.make_edge(the_curve.clone(), first_par, last_par);
        if let Some((pc, _, _)) = the_cos.as_ref() {
            // `BRep_Builder::Range` also re-ranges the COS reps (`BRep_Builder.cxx:1091-1119`),
            // so the copied p-curve keeps its own parameter range.
            crate::brep_offset::update_edge_pcurve(
                &new_edge,
                pc,
                the_face,
                BRepTool::edge_tolerance(&new_edge),
            );
        }
        builder.add(&mut new_edge.0, &first_vertex.0);
        builder.add(&mut new_edge.0, &last_vertex.0);
        new_edge.0.set_orientation(an_orient);

        // `cxx:1487-1493`: `ShapeAnalysis_Edge::CheckSameParameter(NewEdge,
        // theFace, amaxdev)` then `UpdateEdge(NewEdge, amaxdev)` when it answers
        // true. The repo's port (`shhealing::wire_fix::same_parameter`) answers
        // the deviation and not the `ShapeExtend_DONE` flag, so the update runs
        // whenever the deviation is finite - the flag is false only for a
        // degenerate edge or one without a 3D curve, and neither reaches here.
        let amaxdev = crate::shhealing::check_same_parameter(&new_edge, the_face);
        if amaxdev.is_finite() {
            GeometryRegistry::global().raise_edge_tolerance(&new_edge.0, amaxdev);
        }

        // `cxx:1495-1502`: FORWARD appends, REVERSED prepends.
        if an_orient == Orientation::Forward {
            new_edges.push(new_edge.0);
        } else {
            new_edges.insert(0, new_edge.0);
        }

        first_vertex = last_vertex;
        first_par = last_par;
    }

    // `cxx:1505-1509`.
    if !new_edges.is_empty() {
        splits.add(&the_edge.0, new_edges);
    }
}
