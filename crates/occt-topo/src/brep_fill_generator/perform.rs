//! `BRepFill_Generator::Perform` (`BRepFill_Generator.cxx:613-1186`).
//!
//! The two section wires of a pair are walked simultaneously
//! (`cxx:667-671`). A punctual (all-degenerated) wire does not advance
//! (`cxx:726-731`, `:1123-1133`), so a loft to a point pairs one edge with many.
//! Every edge pair yields one ruled face; the connector edges between two
//! consecutive patches are shared through the per-pair vertex map `Map`
//! (`cxx:884-972`), and the pcurves of the four boundary edges are the iso
//! lines of the patch surface (`cxx:1008-1077`).

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom::{GeomBezierCurve, GeomTrimmedCurve};
use occt_geom2d::{Curve2d, Geom2dBezierCurve, Geom2dLine};

use crate::abs::Orientation;
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Shell, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, vertices_of, ShapeMapKey};

use super::{create_kpart, detect_kpart, Generator, ThruSectionErrorStatus};

/// `Perform()` (`cxx:613-1186`).
pub(crate) fn perform(g: &mut Generator) {
    // `cxx:615`.
    g.status = ThruSectionErrorStatus::Done;

    let builder = TopoBuilder::new();
    let mut shell = Shell::new();
    let nb = g.wires.len();
    let mut degen_first = false;
    let mut degen_last = false;
    let mut u_closed = false;

    for i in 0..nb.saturating_sub(1) {
        let wire1 = g.wires[i].clone();
        let wire2 = g.wires[i + 1].clone();
        let outcome = rule_pair(
            g,
            &builder,
            &mut shell,
            &wire1,
            &wire2,
            i,
            nb,
            &mut degen_first,
            &mut degen_last,
            &mut u_closed,
        );
        if let Err(status) = outcome {
            g.status = status;
            return;
        }
    }

    // `cxx:1136-1139`: all vertices of the shell already belong to the section
    // wires, so `BRepLib::SameParameter(theShape)` (no reshaper) suffices.
    if g.mutable_input {
        same_parameter_shape(&shell.0, CONFUSION, false);
    }

    // `cxx:1159-1162`.
    if u_closed && degen_first && degen_last {
        shell.0.set_closed(true);
    }
    g.shell = Some(shell);
}

#[allow(clippy::too_many_arguments)]
fn rule_pair(
    g: &mut Generator,
    builder: &TopoBuilder,
    shell: &mut Shell,
    wire1: &Wire,
    wire2: &Wire,
    i: usize,
    nb: usize,
    degen_first: &mut bool,
    degen_last: &mut bool,
    u_closed: &mut bool,
) -> Result<(), ThruSectionErrorStatus> {
    // `cxx:640-654`.
    let mut w_point1 = false;
    if i == 0 {
        w_point1 = true;
        let mut ex1 = WireExplorer::of_wire(wire1);
        while ex1.more() {
            if let Some(e) = ex1.current() {
                w_point1 = w_point1 && BRepTool::is_degenerated(e);
            }
            ex1.next();
        }
        *degen_first = w_point1;

        // `TopExp::Vertices(Wire1, V1, V2); uClosed = V1.IsSame(V2)`.
        let verts = vertices_of(&wire1.0);
        *u_closed = match (verts.first(), verts.last()) {
            (Some(a), Some(b)) => is_same(a, b),
            _ => false,
        };
    }

    // `cxx:656-664`.
    let mut w_point2 = false;
    if i == nb - 2 {
        w_point2 = true;
        let mut ex2 = WireExplorer::of_wire(wire2);
        while ex2.more() {
            if let Some(e) = ex2.current() {
                w_point2 = w_point2 && BRepTool::is_degenerated(e);
            }
            ex2.next();
        }
        *degen_last = w_point2;
    }

    // `cxx:667-668`.
    let mut ex1 = WireExplorer::of_wire(wire1);
    let mut ex2 = WireExplorer::of_wire(wire2);

    // `cxx:669`: the per-pair vertex map.
    let mut map: Vec<(ShapeMapKey, Edge)> = Vec::new();

    let mut tant_que = ex1.more() && ex2.more();
    while tant_que {
        // `cxx:676-677`.
        let (Some(an_or_edge1), Some(an_or_edge2)) = (ex1.current_edge(), ex2.current_edge())
        else {
            break;
        };
        let degen1 = BRepTool::is_degenerated(&an_or_edge1);
        let degen2 = BRepTool::is_degenerated(&an_or_edge2);

        // `cxx:682-700`.
        let mut edge1 = if degen1 {
            let mut e = empty_copied_edge_with_geometry(&an_or_edge1);
            let (v1f, v1l) = edge_vertices(&an_or_edge1);
            let (Some(v1f), Some(v1l)) = (v1f, v1l) else {
                return Err(ThruSectionErrorStatus::Failed);
            };
            add_oriented_vertex(builder, &mut e, &v1f, Orientation::Forward);
            add_oriented_vertex(builder, &mut e, &v1l, Orientation::Reversed);
            GeometryRegistry::global().set_edge_range(&e.0, 0.0, 1.0);
            g.bind_old_new(&an_or_edge1.0, e.0.clone());
            e
        } else {
            an_or_edge1.clone()
        };

        // `cxx:702-718`.
        let mut edge2 = if degen2 {
            let mut e = empty_copied_edge_with_geometry(&an_or_edge2);
            let (v2f, v2l) = edge_vertices(&an_or_edge2);
            let (Some(v2f), Some(v2l)) = (v2f, v2l) else {
                return Err(ThruSectionErrorStatus::Failed);
            };
            add_oriented_vertex(builder, &mut e, &v2f, Orientation::Forward);
            add_oriented_vertex(builder, &mut e, &v2l, Orientation::Reversed);
            GeometryRegistry::global().set_edge_range(&e.0, 0.0, 1.0);
            g.bind_old_new(&an_or_edge2.0, e.0.clone());
            e
        } else {
            an_or_edge2.clone()
        };

        // `cxx:720-721`.
        let periodic = (edge_curve_is_closed(&edge1) || degen1)
            && (edge_curve_is_closed(&edge2) || degen2);

        // `cxx:726-731`.
        if !w_point1 {
            ex1.next();
        }
        if !w_point2 {
            ex2.next();
        }

        // `cxx:734-762`.
        let (v1a, v1b) = edge_vertices(&edge1);
        let (v1f, v1l) = if edge1.0.orientation() == Orientation::Reversed {
            (v1b, v1a)
        } else {
            (v1a, v1b)
        };
        let (v2a, v2b) = edge_vertices(&edge2);
        let (v2f, v2l) = if edge2.0.orientation() == Orientation::Reversed {
            (v2b, v2a)
        } else {
            (v2a, v2b)
        };
        let (Some(v1f), Some(v1l), Some(v2f), Some(v2l)) = (v1f, v1l, v2f, v2l) else {
            return Err(ThruSectionErrorStatus::Failed);
        };
        let (vf_to_map, vl_to_map) = if degen1 {
            (v2f.clone(), v2l.clone())
        } else {
            (v1f.clone(), v1l.clone())
        };

        // `cxx:764-769`.
        let itype = detect_kpart(&edge1, &edge2);
        if itype == -1 {
            return Err(ThruSectionErrorStatus::Null3DCurve);
        }

        // `cxx:771-868`: the patch surface.
        let surf: Arc<dyn Surface> = if itype == 0 {
            let mut c1 = section_curve(&edge1, degen1, &v1f, &v1l)?;
            if edge1.0.orientation() == Orientation::Reversed {
                c1.reverse();
            }
            let mut c2 = section_curve(&edge2, degen2, &v2f, &v2l)?;
            if edge2.0.orientation() == Orientation::Reversed {
                c2.reverse();
            }
            occt_geom::geom_fill::generator_surface(&[Arc::from(c1), Arc::from(c2)])
                .map_err(|_| ThruSectionErrorStatus::Failed)?
        } else {
            let (done, surface) = create_kpart(&edge1, &edge2, itype);
            if !done {
                return Err(ThruSectionErrorStatus::Null3DCurve);
            }
            surface.ok_or(ThruSectionErrorStatus::Failed)?
        };
        let mut face = builder.make_face_with_tolerance(surf.clone(), &[], CONFUSION);

        // `cxx:871-882`: `Surf->Bounds(f1, l1, f2, l2)`.
        let (f1, l1) = surf.u_range();
        let (f2, l2) = surf.v_range();
        let (first, last) = if itype == 0 { (f2, l2) } else { (0.0, 1.0) };

        // `cxx:884-922`: Edge3 (iso U = f1).
        let mut edge3_from_map = false;
        let mut edge3 = if let Some(stored) = map_get(&map, &vf_to_map) {
            edge3_from_map = true;
            reversed_edge(&stored)
        } else {
            let mut e = if is_same(&v1f.0, &v2f.0) {
                // `cxx:892-896`: a punctual connector.
                let e = Edge::new();
                GeometryRegistry::global().set_degenerated(&e.0, true);
                GeometryRegistry::global().set_curve3d_null(&e.0, true);
                e
            } else {
                let cc: Arc<dyn Curve> = if itype == 0 {
                    surf
                        .u_iso_curve(f1)
                        .ok_or(ThruSectionErrorStatus::Failed)?
                } else {
                    Arc::from(bezier2(
                        BRepTool::vertex_point_world(&v1f),
                        BRepTool::vertex_point_world(&v2f),
                    ))
                };
                make_edge_with_tolerance(cc, CONFUSION)
            };
            add_oriented_vertex(builder, &mut e, &v1f, Orientation::Forward);
            add_oriented_vertex(builder, &mut e, &v2f, Orientation::Reversed);
            GeometryRegistry::global().set_edge_range(&e.0, first, last);
            e
        };
        if !edge3_from_map {
            // `cxx:921-922`.
            edge3.0.reverse();
            map_bind(&mut map, &vf_to_map, &edge3);
        }

        // `cxx:924-972`: Edge4 (iso U = l1).
        let mut common_edge = false;
        if let Some(stored) = map_get(&map, &vl_to_map) {
            let common_e = reversed_edge(&stored);
            let (cv1, cv2) = edge_vertices(&common_e);
            common_edge = match (cv1, cv2) {
                (Some(a), Some(b)) => is_same(&a.0, &v1l.0) && is_same(&b.0, &v2l.0),
                _ => false,
            };
        }
        let edge4 = if common_edge {
            let stored = map_get(&map, &vl_to_map).ok_or(ThruSectionErrorStatus::Failed)?;
            reversed_edge(&stored)
        } else {
            let mut e = if is_same(&v1l.0, &v2l.0) {
                let e = Edge::new();
                GeometryRegistry::global().set_degenerated(&e.0, true);
                GeometryRegistry::global().set_curve3d_null(&e.0, true);
                e
            } else {
                let cc: Arc<dyn Curve> = if itype == 0 {
                    surf
                        .u_iso_curve(l1)
                        .ok_or(ThruSectionErrorStatus::Failed)?
                } else {
                    Arc::from(bezier2(
                        BRepTool::vertex_point_world(&v1l),
                        BRepTool::vertex_point_world(&v2l),
                    ))
                };
                make_edge_with_tolerance(cc, CONFUSION)
            };
            add_oriented_vertex(builder, &mut e, &v1l, Orientation::Forward);
            add_oriented_vertex(builder, &mut e, &v2l, Orientation::Reversed);
            GeometryRegistry::global().set_edge_range(&e.0, first, last);
            map_bind(&mut map, &vl_to_map, &e);
            e
        };

        // UNPORTED: `cxx:974-1006` (`!myMutableInput`: `CreateNewEdge`) and
        // `cxx:1136-1157` (`BRepTools_ReShape`) — see the module docs.

        // `cxx:1008-1077`: the pcurves.
        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&face.0);
        if itype != 4 {
            if edge1.0.orientation() == Orientation::Reversed {
                reg.set_edge_pcurve(&edge1.0, face_key, uv_line(0.0, f2, -1.0, 0.0));
                reg.set_pcurve_range(&edge1.0, face_key, -l1, -f1);
            } else {
                reg.set_edge_pcurve(&edge1.0, face_key, uv_line(0.0, f2, 1.0, 0.0));
                reg.set_pcurve_range(&edge1.0, face_key, f1, l1);
            }
            if edge2.0.orientation() == Orientation::Reversed {
                reg.set_edge_pcurve(&edge2.0, face_key, uv_line(0.0, l2, -1.0, 0.0));
                reg.set_pcurve_range(&edge2.0, face_key, -l1, -f1);
            } else {
                reg.set_edge_pcurve(&edge2.0, face_key, uv_line(0.0, l2, 1.0, 0.0));
                reg.set_pcurve_range(&edge2.0, face_key, f1, l1);
            }
        }
        if itype == 0 {
            if periodic {
                reg.set_edge_pcurves(
                    &edge3.0,
                    face_key,
                    vec![uv_line(l1, 0.0, 0.0, 1.0), uv_line(f1, 0.0, 0.0, 1.0)],
                );
            } else {
                reg.set_edge_pcurve(&edge3.0, face_key, uv_line(f1, 0.0, 0.0, 1.0));
                reg.set_edge_pcurve(&edge4.0, face_key, uv_line(l1, 0.0, 0.0, 1.0));
            }
        } else if periodic {
            reg.set_edge_pcurves(
                &edge3.0,
                face_key,
                vec![
                    uv_bezier(&[(l1, f2), (l1, l2)]),
                    uv_bezier(&[(f1, f2), (f1, l2)]),
                ],
            );
        } else if itype != 4 {
            reg.set_edge_pcurve(&edge3.0, face_key, uv_bezier(&[(f1, f2), (f1, l2)]));
            reg.set_edge_pcurve(&edge4.0, face_key, uv_bezier(&[(l1, f2), (l1, l2)]));
        }

        // `cxx:1079-1087`.
        for e in [&edge1, &edge2, &edge3, &edge4] {
            reg.set_same_parameter(&e.0, false);
            reg.set_same_range(&e.0, false);
        }

        // `cxx:1090-1103`.
        let mut wire = Wire::new();
        if !degen1 || itype != 4 {
            builder.add_edge(&mut wire, &edge1);
        }
        builder.add_edge(&mut wire, &edge4);
        if !degen2 || itype != 4 {
            builder.add_edge(&mut wire, &reversed_edge(&edge2));
        }
        builder.add_edge(&mut wire, &edge3);
        builder.add_wire(&mut face, &wire);
        builder.add_face(shell, &face);

        // `cxx:1106-1121`: complete `myMap` for the first edge.
        if !degen1 || itype != 4 {
            let a_re = if degen1 {
                edge1.0.clone()
            } else {
                an_or_edge1.0.clone()
            };
            g.append_generated(&a_re, &face);
        }

        // `cxx:1123-1133`.
        tant_que = ex1.more() && ex2.more();
        if w_point1 {
            tant_que = ex2.more();
        }
        if w_point2 {
            tant_que = ex1.more();
        }
    }
    Ok(())
}

/// `cxx:779-849`: the section curve of a patch, trimmed to the edge range and
/// moved by the edge location. A degenerated edge becomes the 2-pole Bezier
/// through its two vertices (`cxx:779-786`, `:794-799`); note the second one
/// lists the vertices in reverse order.
fn section_curve(
    edge: &Edge,
    degen: bool,
    v_first: &Vertex,
    v_last: &Vertex,
) -> Result<Box<dyn Curve>, ThruSectionErrorStatus> {
    if degen {
        let mut c = bezier2(
            BRepTool::vertex_point_world(v_first),
            BRepTool::vertex_point_world(v_last),
        )
        .clone_dyn();
        // `cxx:810-818`: `L1` is the identity for a degenerated edge.
        c.transform(&edge.0.location().transformation());
        return Ok(c);
    }
    let curve = BRepTool::edge_curve(edge).ok_or(ThruSectionErrorStatus::Null3DCurve)?;
    let (first, last) = BRepTool::edge_parameters(edge);
    let trimmed = (first - curve.first_parameter()).abs() > PCONFUSION
        || (last - curve.last_parameter()).abs() > PCONFUSION;
    let mut owned: Box<dyn Curve> = if trimmed {
        Box::new(GeomTrimmedCurve::new(curve, first, last))
    } else {
        curve.clone_dyn()
    };
    // `cxx:825-828`: `if (!SameLoc) C1->Transform(L1.Transformation())`.
    owned.transform(&edge.0.location().transformation());
    Ok(owned)
}

/// `TopoDS_Edge::EmptyCopied()` keeps the `BRep_TEdge` representations
/// ("a new TShape with the same geometry and no sub-shapes",
/// `TopoDS_Shape.hxx:290-301`), so the degenerated flag and the null 3D curve of
/// `BRep_Builder::Degenerated` survive the copy.
fn empty_copied_edge_with_geometry(edge: &Edge) -> Edge {
    let copy = Edge::wrap(edge.0.empty_copied()).expect("TopoDS_Edge::EmptyCopied");
    if let Some(geom) = GeometryRegistry::global().edge_geom(&edge.0) {
        GeometryRegistry::global().set_edge(&copy.0, geom);
    }
    copy
}

/// `BRep_Builder::MakeEdge(E, C, Tol)` (`BRep_Builder.cxx:755-770`).
fn make_edge_with_tolerance(curve: Arc<dyn Curve>, tolerance: f64) -> Edge {
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    let e = TopoBuilder::new().make_edge(curve, first, last);
    GeometryRegistry::global().set_edge_tolerance(&e.0, tolerance);
    e
}

/// `new Geom_BezierCurve(Extremities)` over two points (`cxx:786`, `:799`).
fn bezier2(p1: GpPnt, p2: GpPnt) -> Arc<dyn Curve> {
    Arc::new(GeomBezierCurve::new(vec![p1, p2]).expect("Geom_BezierCurve(2 poles)"))
}

fn add_oriented_vertex(builder: &TopoBuilder, edge: &mut Edge, v: &Vertex, ori: Orientation) {
    let mut s = v.0.clone();
    s.set_orientation(ori);
    builder.add(&mut edge.0, &s);
}

/// `TopoDS_Shape::Reversed()`: a new occurrence with the flipped orientation.
fn reversed_edge(e: &Edge) -> Edge {
    let mut r = e.clone();
    r.0.reverse();
    r
}

/// `BRep_Tool::IsClosed(Edge)` (`BRep_Tool.cxx:1750-1755`): the geometric
/// closure of the 3D curve.
fn edge_curve_is_closed(e: &Edge) -> bool {
    BRepTool::edge_curve(e).map(|c| c.is_closed()).unwrap_or(false)
}

fn map_get(map: &[(ShapeMapKey, Edge)], v: &Vertex) -> Option<Edge> {
    let key = shape_map_key(&v.0);
    map.iter().find(|(k, _)| *k == key).map(|(_, e)| e.clone())
}

fn map_bind(map: &mut Vec<(ShapeMapKey, Edge)>, v: &TopoShape, e: &Edge) {
    let key = shape_map_key(v);
    if let Some(slot) = map.iter_mut().find(|(k, _)| *k == key) {
        slot.1 = e.clone();
    } else {
        map.push((key, e.clone()));
    }
}

fn uv_line(x: f64, y: f64, dx: f64, dy: f64) -> Arc<dyn Curve2d> {
    Arc::new(Geom2dLine::from_pnt_dir(
        GpPnt2d::new(x, y),
        GpDir2d::new(dx, dy).expect("gp_Dir2d(pcurve direction)"),
    ))
}

fn uv_bezier(pts: &[(f64, f64)]) -> Arc<dyn Curve2d> {
    let poles: Vec<GpPnt2d> = pts.iter().map(|(x, y)| GpPnt2d::new(*x, *y)).collect();
    Arc::new(Geom2dBezierCurve::new(poles).expect("Geom2d_BezierCurve(2 poles)"))
}
