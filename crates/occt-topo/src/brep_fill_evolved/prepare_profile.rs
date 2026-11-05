//! `BRepFill_Evolved::PrepareProfile` (`BRepFill_Evolved.cxx:1553-1711`) and
//! its helper `CutEdgeProf` (`cxx:3122-3266`).
//!
//! The profile is projected on the `yOz` plane, cut at the extrema of the
//! distance to the `Oz` axis and at the intersections with `Oy`
//! (`CutEdgeProf`), then regrouped so that every work wire stays on one side
//! of `Oz` and holds only vertical, only planar or only other edges
//! (`PrepareProfile`).

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir, GpDir2d, GpLin2d, GpPln, GpPnt2d, GpTrsf};
use occt_core::precision::{Precision, REAL_EPSILON};
use occt_geom::curve::Curve;
use occt_geom::geom_proj_lib;
use occt_geom::plane::GeomPlane;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::Surface;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::extrema2d::ExtremaExtCC2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::bnd_lib_add2d;
use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, edges_of};

use super::shape_map::ShapeMap;
use super::{is_planar, is_vertical, FILL_CONFUSION};

/// `RealFirst()` / `RealLast()` (`Standard_Real.hxx:120-135`): the values
/// `Bnd_Box2d::Get` answers for a void box (`Bnd_Box2d.cxx`, `Get`).
const REAL_FIRST: f64 = -f64::MAX;
const REAL_LAST: f64 = f64::MAX;

/// `BRepLib_MakeEdge(C, V1, V2)` (`BRepLib_MakeEdge.cxx:217-230`): the range
/// comes from projecting both vertices on `C`.
fn make_edge_between(curve: &Arc<dyn Curve>, v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d_projecting(curve, Some(v1), Some(v2))?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}

/// `BRepFill_Confusion()` (`cxx:99-103`).
const CONFUSION: f64 = FILL_CONFUSION;

/// `CutEdgeProf(E, Plane, Line, Cuts, MapVerRefMoved)`
/// (`BRepFill_Evolved.cxx:3122-3266`).
///
/// Answers the cut edges in wire order: `FORWARD` pieces are appended,
/// `REVERSED` pieces prepended (`cxx:3243-3251`).
pub fn cut_edge_prof(
    e: &Edge,
    plane: &GeomPlane,
    line: &Geom2dLine,
    map_ver_ref_moved: &mut ShapeMap<Vertex>,
) -> Vec<TopoShape> {
    let mut cuts: Vec<TopoShape> = Vec::new();
    let builder = TopoBuilder::new();

    // `cxx:3144-3146`: `C = BRep_Tool::Curve(E, L, f, l)`,
    // `CT = new Geom_TrimmedCurve(C, f, l)`, `CT->Transform(L.Transformation())`.
    let (f, l) = BRepTool::edge_parameters(e);
    let Some(curve0) = BRepTool::edge_curve(e) else {
        // `BRep_Tool::Curve` on an edge without a 3D curve raises in OCCT.
        return cuts;
    };
    let Ok(mut ct) = GeomTrimmedCurveBasis::with_sense(curve0, f, l, true) else {
        return cuts;
    };
    let loc = e.0.location();
    if !loc.is_identity() {
        let t: GpTrsf = loc.transformation();
        ct.transform(&t);
    }
    let ct: Arc<dyn Curve> = Arc::new(ct);

    // `cxx:3147-3149`.
    let pln: GpPln = plane.pln().clone();
    let axis = pln.axis();
    let normal: GpDir = *axis.direction();
    let c = geom_proj_lib::project_on_plane(&ct, &pln, &normal);
    let plane_surface: Arc<dyn Surface> = Arc::new(plane.clone());
    let Some(c2d) = geom_proj_lib::curve_2d(&c, &plane_surface) else {
        return cuts;
    };

    // `cxx:3152-3164`: the box of the projection and `B.Get(xmin, U1, xmax, U2)`.
    let c2d_adaptor = Geom2dAdaptorCurve::from_curve(c2d);
    let mut b = BndBox2d::new();
    bnd_lib_add2d::add_adaptor(
        c2d_adaptor.curve().expect("loaded curve").as_ref(),
        CONFUSION,
        &mut b,
    );
    let (u1, u2) = match b.get() {
        // `Bnd_Box2d::Get` on a void box answers `RealFirst` / `RealLast`.
        None => (REAL_FIRST, REAL_LAST),
        Some((_, ymin, _, ymax)) => (ymin, ymax),
    };
    let (cf, cl) = (c2d_adaptor.first_parameter(), c2d_adaptor.last_parameter());

    let mut seq: Vec<f64> = Vec::new();

    // `cxx:3166-3178`: the extrema with the line, no sense for a point.
    if u1 != u2 {
        let c2d_curve = c2d_adaptor.curve().expect("loaded curve");
        let ext = ExtremaExtCC2d::new_range(line, c2d_curve.as_ref(), u1 - 1.0, u2 + 1.0, cf, cl);
        for i in 1..=ext.nb_ext() {
            let (_, _, param_on_second, _) = ext.points(i);
            seq.push(param_on_second);
        }
    }

    // `cxx:3180-3192`: the intersections with `Oy`.
    let mut intersector = Geom2dIntGInter::new();
    {
        let line_adaptor = Geom2dAdaptorCurve::from_curve(Arc::new(line.clone()));
        let domain_line = Geom2dIntGInter::domain_over_range(
            line_adaptor.curve().expect("loaded line").as_ref(),
            line_adaptor.first_parameter(),
            line_adaptor.last_parameter(),
            Precision::INTERSECTION,
        );
        let domain_c = Geom2dIntGInter::domain_over_range(
            c2d_adaptor.curve().expect("loaded curve").as_ref(),
            cf,
            cl,
            Precision::INTERSECTION,
        );
        intersector.perform(
            line_adaptor.curve().expect("loaded line").as_ref(),
            &domain_line,
            c2d_adaptor.curve().expect("loaded curve").as_ref(),
            &domain_c,
            0.0,
            Precision::INTERSECTION,
        );
    }
    for i in 1..=intersector.nb_points() {
        seq.push(intersector.point(i).param_on_second());
    }

    // `cxx:3195-3218`: the endpoint vertices of the projected curve.
    let (vrf, vrl) = edge_vertices(e);
    let (Some(vrf), Some(vrl)) = (vrf, vrl) else {
        return cuts;
    };
    let make_moved = |map: &mut ShapeMap<Vertex>, reference: &Vertex, param: f64| -> Vertex {
        let key = reference.0.clone();
        if let Some(v) = map.find(&key) {
            return v.clone();
        }
        let v = builder.make_vertex(
            c.d0(param),
            BRepTool::vertex_tolerance(reference),
        );
        map.bind(&key, v.clone());
        v
    };
    let mut vf = make_moved(map_ver_ref_moved, &vrf, f);
    let vl = make_moved(map_ver_ref_moved, &vrl, l);

    let push_cut = |cuts: &mut Vec<TopoShape>, mut ee: Edge| {
        ee.0.set_orientation(e.0.orientation());
        if ee.0.orientation() == Orientation::Forward {
            cuts.push(ee.0);
        } else {
            cuts.insert(0, ee.0);
        }
    };

    // `cxx:3220-3256`.
    if !seq.is_empty() {
        super::bubble(&mut seq);
        let mut cur_param = f;
        let mut empty = false;
        while !empty {
            let param = seq.remove(0);
            empty = seq.is_empty();
            if (param - cur_param).abs() > CONFUSION && (param - l).abs() > CONFUSION {
                // `cxx:3238`: `VV = BRepLib_MakeVertex(C->Value(Param))`.
                let vv = builder.make_vertex(c.d0(param), REAL_EPSILON);
                if let Some(ee) = make_edge_between(&c, &vf, &vv) {
                    push_cut(&mut cuts, ee);
                }
                cur_param = param;
                vf = vv;
            }
        }
    }

    // `cxx:3258-3265`: the last piece, up to the original last vertex.
    if let Some(ee) = make_edge_between(&c, &vf, &vl) {
        push_cut(&mut cuts, ee);
    }
    cuts
}

/// `BRepFill_Evolved::PrepareProfile(WorkProf, MapProf)`
/// (`BRepFill_Evolved.cxx:1553-1711`).
pub fn prepare_profile(
    my_profile: &Wire,
    work_prof: &mut Vec<TopoShape>,
    map_prof: &mut ShapeMap<TopoShape>,
) {
    let builder = TopoBuilder::new();
    // `cxx:1565-1566`.
    let plane = GeomPlane::new(super::yoz_plane());
    let line = Geom2dLine::from_pnt_dir(
        GpPnt2d::new(0.0, 0.0),
        GpDir2d::new(0.0, 1.0).expect("gp::OY2d"),
    );

    // `cxx:1568-1576`: `MapVerRefMoved` and the first work wire.
    let mut map_ver_ref_moved: ShapeMap<Vertex> = ShapeMap::new();
    let mut w = builder.make_wire(&[]);
    let mut wp: Vec<TopoShape> = vec![w.0.clone()];

    // `cxx:1578-1641`.
    let mut exp = WireExplorer::of_wire(my_profile);
    while exp.more() {
        let Some(e) = exp.current_edge() else {
            break;
        };
        let mut new_wire = false;

        let mut cuts = cut_edge_prof(&e, &plane, &line, &mut map_ver_ref_moved);

        if cuts.is_empty() {
            // `cxx:1590-1596`: neither extrema nor intersections.
            builder.add(&mut w.0, &e.0);
            map_prof.bind(&e.0, e.0.clone());
        } else {
            while !cuts.is_empty() {
                let ne = cuts.remove(0);
                map_prof.bind(&ne, e.0.clone());
                let (v1, v2) = edge_vertices(&Edge(ne.clone()));
                if let Some(v) = v1.as_ref() {
                    if !map_prof.is_bound(&v.0) {
                        map_prof.bind(&v.0, e.0.clone());
                    }
                }
                if let Some(v) = v2.as_ref() {
                    if !map_prof.is_bound(&v.0) {
                        map_prof.bind(&v.0, e.0.clone());
                    }
                }
                builder.add(&mut w.0, &ne);

                let d1 = v1.as_ref().map(super::distance_to_oz).unwrap_or(0.0);
                let d2 = v2.as_ref().map(super::distance_to_oz).unwrap_or(0.0);
                if d2 < CONFUSION && d1 > CONFUSION {
                    // `cxx:1613-1631`: the piece ends on `Oz`.
                    if cuts.is_empty() {
                        new_wire = true;
                    } else {
                        w = builder.make_wire(&[]);
                        wp.push(w.0.clone());
                    }
                }
            }
        }

        exp.next();
        if exp.more() && new_wire {
            // `cxx:1636-1640`.
            w = builder.make_wire(&[]);
            wp.push(w.0.clone());
        }
    }

    // `cxx:1643-1703`: split the wires that mix vertical / planar / other
    // pieces.
    for cur_w in wp {
        if cur_w.shape_type() != ShapeType::Wire {
            work_prof.push(cur_w);
            continue;
        }
        let cur_w = Wire(cur_w);
        let mut ya_modif = false;
        for ee in edges_of(&cur_w.0) {
            if is_vertical(&ee) || is_planar(&ee) {
                ya_modif = true;
                break;
            }
        }
        if ya_modif {
            // `cxx:1665-1698`.
            let mut status = 0;
            let mut nw = builder.make_wire(&[]);
            for ee in edges_of(&cur_w.0) {
                if is_vertical(&ee) {
                    if status != 3 {
                        nw = builder.make_wire(&[]);
                        work_prof.push(nw.0.clone());
                        status = 3;
                    }
                } else if is_planar(&ee) {
                    if status != 2 {
                        nw = builder.make_wire(&[]);
                        work_prof.push(nw.0.clone());
                        status = 2;
                    }
                } else if status != 1 {
                    nw = builder.make_wire(&[]);
                    work_prof.push(nw.0.clone());
                    status = 1;
                }
                builder.add(&mut nw.0, &ee.0);
            }
        } else {
            work_prof.push(cur_w.0);
        }
    }

    // `cxx:1705-1710`: bind every moved vertex as the image of the reference
    // one it replaces.
    let moved: Vec<(TopoShape, Vertex)> = map_ver_ref_moved
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    for (reference, moved_vertex) in moved {
        map_prof.bind(&moved_vertex.0, reference);
    }
}

/// The 2D line of `PrepareProfile` (`cxx:1566`): `gp::OY2d()`.
pub fn oy2d_line() -> GpLin2d {
    GpLin2d::from_pnt_dir(
        GpPnt2d::new(0.0, 0.0),
        GpDir2d::new(0.0, 1.0).expect("gp::OY2d"),
    )
}
