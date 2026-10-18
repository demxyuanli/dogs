//! `BRepFill::Face` and `BRepFill::Shell`.
//!
//! Source: `BRepFill.cxx:211-422` (`Face`) and `BRepFill.cxx:426-667` (`Shell`).
//! Both build the ruled patch with `GeomFill_Generator` (`cxx:295-298`,
//! `cxx:524-527`).

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edge_vertices, edges_of_wire, faces_of, is_same, wires_of_face};

/// `BRepFill::Face` (`BRepFill.cxx:211-422`).
pub fn face(edge1: &Edge, edge2: &Edge) -> Result<Face, String> {
    let (same_loc, location) = shared_location(edge1, edge2);
    let (c1, v1f, v1l) = prepare_curve(edge1, same_loc)?;
    let (c2, v2f, v2l) = prepare_curve(edge2, same_loc)?;
    let closed = is_same(&v1f.0, &v1l.0) && is_same(&v2f.0, &v2l.0);
    let surf = occt_geom::geom_fill::generator_surface(&[c1, c2])?;
    let (u0, u1) = surf.u_range();
    let (v0, v1) = surf.v_range();
    let mut face = TopoBuilder::new().make_face(surf.clone(), &[]);
    GeometryRegistry::global().set_face_tolerance(&face.0, CONFUSION);

    let mut edge3 = iso_edge(surf.as_ref(), u0, v0, v1, &v1f, &v2f)?;
    let edge4 = if closed {
        edge3.clone()
    } else {
        iso_edge(surf.as_ref(), u1, v0, v1, &v1l, &v2l)?
    };
    edge3.0.reverse();

    let mut wire = Wire::new();
    let builder = TopoBuilder::new();
    builder.add_edge(&mut wire, edge1);
    builder.add_edge(&mut wire, &edge4);
    let mut edge2_rev = edge2.clone();
    edge2_rev.0.reverse();
    builder.add_edge(&mut wire, &edge2_rev);
    builder.add_edge(&mut wire, &edge3);
    wire.0.set_closed(true);
    builder.add_wire(&mut face, &wire);

    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    bind_profile_pcurve(&reg, edge1, face_key, 0.0, v0, u0, u1);
    bind_profile_pcurve(&reg, edge2, face_key, 0.0, v1, u0, u1);
    if closed {
        reg.set_edge_pcurves(
            &edge3.0,
            face_key,
            vec![uv_line(u1, 0.0, 0.0, 1.0), uv_line(u0, 0.0, 0.0, 1.0)],
        );
    } else {
        reg.set_edge_pcurve(&edge3.0, face_key, uv_line(u0, 0.0, 0.0, 1.0));
        reg.set_edge_pcurve(&edge4.0, face_key, uv_line(u1, 0.0, 0.0, 1.0));
    }
    clear_same_flags(edge1);
    clear_same_flags(edge2);
    clear_same_flags(&edge3);
    clear_same_flags(&edge4);
    same_parameter_shape(&face.0, 1.0e-5, false);
    if same_loc {
        face.0.move_location(&location);
    }
    Ok(face)
}

/// `BRepFill::Shell` (`BRepFill.cxx:426-667`).
pub fn shell(wire1: &Wire, wire2: &Wire) -> Result<Shell, String> {
    let mut shell = Shell::new();
    let edges1 = edges_of_wire(wire1);
    let edges2 = edges_of_wire(wire2);
    let closed = wire1.0.is_closed() && wire2.0.is_closed();
    let mut edge4: Option<Edge> = None;
    let mut couture: Option<Edge> = None;
    let n = edges1.len().min(edges2.len());
    for i in 0..n {
        let edge1 = &edges1[i];
        let edge2 = &edges2[i];
        let periodic = edge_is_closed(edge1) && edge_is_closed(edge2);
        let (same_loc, location) = shared_location(edge1, edge2);
        let (c1, v1f, v1l) = prepare_curve(edge1, same_loc)?;
        let (c2, v2f, v2l) = prepare_curve(edge2, same_loc)?;
        let surf = occt_geom::geom_fill::generator_surface(&[c1, c2])?;
        let (u0, u1) = surf.u_range();
        let (v0, v1) = surf.v_range();
        let mut face = TopoBuilder::new().make_face(surf.clone(), &[]);
        GeometryRegistry::global().set_face_tolerance(&face.0, CONFUSION);

        let edge3 = if i == 0 {
            let built = iso_edge(surf.as_ref(), u0, v0, v1, &v1f, &v2f)?;
            if closed {
                couture = Some(built.clone());
            }
            let mut rev = built;
            rev.0.reverse();
            rev
        } else {
            let mut rev = edge4.clone().ok_or_else(|| "BRepFill::Shell".to_string())?;
            rev.0.reverse();
            rev
        };
        let edge4_now = if closed && i + 1 == n {
            couture.clone().ok_or_else(|| "BRepFill::Shell".to_string())?
        } else {
            iso_edge(surf.as_ref(), u1, v0, v1, &v1l, &v2l)?
        };

        let mut wire = Wire::new();
        let builder = TopoBuilder::new();
        builder.add_edge(&mut wire, edge1);
        builder.add_edge(&mut wire, &edge4_now);
        let mut edge2_rev = edge2.clone();
        edge2_rev.0.reverse();
        builder.add_edge(&mut wire, &edge2_rev);
        builder.add_edge(&mut wire, &edge3);
        wire.0.set_closed(true);
        builder.add_wire(&mut face, &wire);
        if same_loc {
            face.0.move_location(&location);
        }
        builder.add_face(&mut shell, &face);

        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&face.0);
        bind_profile_pcurve(&reg, edge1, face_key, 0.0, v0, u0, u1);
        bind_profile_pcurve(&reg, edge2, face_key, 0.0, v1, u0, u1);
        if periodic {
            reg.set_edge_pcurves(
                &edge3.0,
                face_key,
                vec![uv_line(u1, 0.0, 0.0, 1.0), uv_line(u0, 0.0, 0.0, 1.0)],
            );
        } else {
            reg.set_edge_pcurve(&edge3.0, face_key, uv_line(u0, 0.0, 0.0, 1.0));
            reg.set_edge_pcurve(&edge4_now.0, face_key, uv_line(u1, 0.0, 0.0, 1.0));
        }
        clear_same_flags(edge1);
        clear_same_flags(edge2);
        clear_same_flags(&edge3);
        clear_same_flags(&edge4_now);
        edge4 = Some(edge4_now);
    }
    shell.0.set_closed(shell_is_closed(&shell.0));
    same_parameter_shape(&shell.0, 1.0e-5, false);
    Ok(shell)
}

fn shared_location(edge1: &Edge, edge2: &Edge) -> (bool, TopLocLocation) {
    let l1 = edge1.0.location().clone();
    let l2 = edge2.0.location().clone();
    if l1.is_equal(&l2) {
        (true, l1)
    } else {
        (false, TopLocLocation::identity())
    }
}

/// Trim or copy, then `Transform` when the edges do not share a location,
/// then `Reverse` when the edge is reversed (`BRepFill.cxx:241-290`).
fn prepare_curve(edge: &Edge, same_loc: bool) -> Result<(Arc<dyn Curve>, Vertex, Vertex), String> {
    let reg = GeometryRegistry::global();
    let curve = reg
        .edge_curve(&edge.0)
        .ok_or_else(|| "BRep_Tool::Curve".to_string())?;
    let (first, last) = reg.edge_parameters(&edge.0);
    let trimmed = (first - curve.first_parameter()).abs() > PCONFUSION
        || (last - curve.last_parameter()).abs() > PCONFUSION;
    let mut owned = if trimmed {
        Box::new(occt_geom::geom_fill::ParameterTrim::new(curve, first, last)) as Box<dyn Curve>
    } else {
        curve.clone_dyn()
    };
    if !same_loc {
        owned.transform(&edge.0.location().transformation());
    }
    let (stored_first, stored_last) = edge_vertices(edge);
    let (v_first, v_last) = if edge.0.orientation() == Orientation::Reversed {
        owned.reverse();
        (stored_last, stored_first)
    } else {
        (stored_first, stored_last)
    };
    let v_first = v_first.ok_or_else(|| "TopExp::Vertices".to_string())?;
    let v_last = v_last.ok_or_else(|| "TopExp::Vertices".to_string())?;
    Ok((Arc::from(owned), v_first, v_last))
}

fn iso_edge(
    surf: &dyn Surface,
    u: f64,
    v0: f64,
    v1: f64,
    v_first: &Vertex,
    v_second: &Vertex,
) -> Result<Edge, String> {
    let iso = surf
        .u_iso_curve(u)
        .ok_or_else(|| "Geom_Surface::UIso".to_string())?;
    let tol = BRepTool::vertex_tolerance(v_first).max(BRepTool::vertex_tolerance(v_second));
    let mut edge = if iso.d0(v0).distance(&iso.d0(v1)) > tol {
        TopoBuilder::new().make_edge(iso, v0, v1)
    } else {
        TopoBuilder::new().make_degenerated_edge(&BRepTool::vertex_point(v_first))
    };
    GeometryRegistry::global().raise_edge_tolerance(&edge.0, CONFUSION);
    TopoBuilder::new().add_edge_vertices(&mut edge, v_first, v_second);
    GeometryRegistry::global().set_ranges_all(&edge.0, v0, v1);
    Ok(edge)
}

fn bind_profile_pcurve(
    reg: &GeometryRegistry,
    edge: &Edge,
    face_key: usize,
    _origin_u: f64,
    v: f64,
    u0: f64,
    u1: f64,
) {
    // `BRepFill.cxx:369-389`: a reversed edge stores the negative-X line and
    // the negated U range.
    if edge.0.orientation() == Orientation::Reversed {
        reg.set_edge_pcurve(&edge.0, face_key, uv_line(0.0, v, -1.0, 0.0));
        reg.set_pcurve_range(&edge.0, face_key, -u1, -u0);
    } else {
        reg.set_edge_pcurve(&edge.0, face_key, uv_line(0.0, v, 1.0, 0.0));
        reg.set_pcurve_range(&edge.0, face_key, u0, u1);
    }
}

fn uv_line(x: f64, y: f64, dx: f64, dy: f64) -> Arc<dyn Curve2d> {
    let dir = GpDir2d::new(dx, dy).unwrap_or(GpDir2d { x: dx, y: dy });
    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(x, y), dir))
}

fn clear_same_flags(edge: &Edge) {
    let reg = GeometryRegistry::global();
    reg.set_same_parameter(&edge.0, false);
    reg.set_same_range(&edge.0, false);
}

/// `BRep_Tool::IsClosed(Edge)` (`BRep_Tool.cxx:1750-1755`).
fn edge_is_closed(edge: &Edge) -> bool {
    match edge_vertices(edge) {
        (Some(first), Some(last)) => is_same(&first.0, &last.0),
        _ => false,
    }
}

/// `BRep_Tool::IsClosed(Shell)` (`BRep_Tool.cxx:1709-1728`).
fn shell_is_closed(shell: &TopoShape) -> bool {
    let mut odd: Vec<TopoShape> = Vec::new();
    let mut has_bound = false;
    for face in faces_of(shell) {
        for wire in wires_of_face(&face) {
            for edge in edges_of_wire(&wire) {
                if BRepTool::is_degenerated(&edge)
                    || matches!(
                        edge.0.orientation(),
                        Orientation::Internal | Orientation::External
                    )
                {
                    continue;
                }
                has_bound = true;
                toggle_same(&mut odd, &edge.0);
            }
        }
    }
    has_bound && odd.is_empty()
}

/// `TopTools_Map` toggle of `BRepLib_Fill`'s edge pairing: `IsSame` already
/// covers `TShape` and location (`TopoDS_Shape.hxx:262-267`).
fn toggle_same(odd: &mut Vec<TopoShape>, shape: &TopoShape) {
    if let Some(i) = odd.iter().position(|kept| is_same(kept, shape)) {
        odd.remove(i);
    } else {
        odd.push(shape.clone());
    }
}
