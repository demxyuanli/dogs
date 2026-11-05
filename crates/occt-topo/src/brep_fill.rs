//! `BRepFill::Face`, `BRepFill::Shell`, `BRepFill::ComputeACR`,
//! `BRepFill::InsertACR` and `BRepFill::SearchOrigin`.
//!
//! Source: `BRepFill.cxx:211-422` (`Face`), `cxx:426-667` (`Shell`),
//! `cxx:1027-1063` (`ComputeACR`), `cxx:1067-1155` (`InsertACR`) and
//! `cxx:885-1023` (`SearchOrigin`), plus the
//! two file-local statics `MakeWire` (`cxx:81-117`) and `CutEdge`
//! (`cxx:119-142`). `TrimEdge` (`cxx:143-207`) is shared with
//! `BRepFill_CompatibleWires` and lives in
//! [`crate::brep_fill_compatible_wires::statics::trim_edge`].
//!
//! `Face` and `Shell` build the ruled patch with `GeomFill_Generator`
//! (`cxx:295-298`, `cxx:524-527`).
//!
//! Translation notes for [`search_origin`]:
//! - `BRepExtrema_DistShapeShape(V, W)` (`cxx:904`) is re-derived from the two
//!   solution families the class itself enumerates for a vertex / wire pair:
//!   every vertex of `W` (vertex-vertex extrema, exact point distance) and every
//!   edge of `W` (point-curve extrema through
//!   [`crate::brep_extrema::closest_point_on_edge`], the port's analytic
//!   `Extrema_ExtPC`). `SupportTypeShape2` (`cxx:918`) is `IsVertex` exactly
//!   when the vertex family wins; the loop at `cxx:909-916` only keeps the
//!   minimum, so the solution order does not matter except for an exact tie.
//! - `ParOnEdgeS2(isol, theparam)` (`cxx:944`) is the parameter
//!   [`crate::brep_extrema::closest_point_on_edge`] returns for the winning
//!   edge.
//! - `BRepAdaptor_Curve(E).D1(theparam, Pe, Ve)` (`cxx:1005-1008`) is
//!   [`BRepTool::edge_curve_world`] evaluated at `theparam`; `E.Orientation()`
//!   is the edge's stored orientation.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpDir2d, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_find_surface::FindSurface;
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edge_vertices, edges_of_wire, explorer_shapes, faces_of, is_same, map_shapes_and_ancestors,
    shape_map_key, wires_of_face,
};

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

/// `BRepFill::ComputeACR(wire, ACR)` (`BRepFill.cxx:1027-1063`): the reduced
/// curvilinear abscissae of the wire's vertices.
///
/// `ACR` follows OCCT's `NCollection_Array1<double>(0, NbEdge)` bounds: entry
/// `0` holds the total length before the normalisation and is the divisor, so
/// the answer is returned as `ACR[0..=nb_edges]` with `ACR[0]` overwritten by
/// the normalised value `1.0` (or left at the total for a punctual wire).
///
/// `GCPnts_AbscissaPoint::Length(anEcur)` is the three-argument form: the
/// parameters are `BRepAdaptor_Curve`'s own first / last and there is no
/// tolerance argument.
pub fn compute_acr(wire: &Wire) -> Result<Vec<f64>, String> {
    let mut nb_edges = 0usize;
    let mut acr: Vec<f64> = vec![0.0];

    let mut wexp = crate::brep_tools_wire_explorer::WireExplorer::new();
    wexp.init(wire, None);
    while wexp.more() {
        nb_edges += 1;
        let current = wexp.current_edge();
        let Some(edge) = current else {
            wexp.next();
            continue;
        };
        // `ACR(nbEdges) = ACR(nbEdges - 1)` (`cxx:1039`).
        acr.push(acr[nb_edges - 1]);
        if !BRepTool::is_degenerated(&edge) {
            // `BRepAdaptor_Curve anEcur(Ecur); ACR(nbEdges) +=
            //  GCPnts_AbscissaPoint::Length(anEcur)` (`cxx:1042-1043`).
            let curve = BRepTool::edge_curve_adaptor_world(&edge)
                .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
            acr[nb_edges] += occt_geom::gcpnts::curve_length_between(
                &*curve,
                curve.first_parameter(),
                curve.last_parameter(),
                None,
            );
        }
        wexp.next();
    }

    // `ACR(0) = ACR(nbEdges)` (`cxx:1048`).
    acr[0] = acr[nb_edges];

    // `cxx:1051-1062`.
    if acr[0] > CONFUSION {
        for i in 1..=nb_edges {
            acr[i] /= acr[0];
        }
    } else {
        // Punctual wire.
        acr[nb_edges] = 1.0;
    }
    Ok(acr)
}

/// `BRepFill::InsertACR(wire, ACRcuts, prec)` (`BRepFill.cxx:1067-1155`): cut
/// `wire` at every reduced abscissa of `acr_cuts` and rebuild it with
/// `BRepLib_MakeWire`.
///
/// `TrimEdge` is the `BRepFill.cxx:143-207` copy, i.e.
/// [`crate::brep_fill_compatible_wires::statics::trim_edge`] with
/// `reject_empty_piece = false`.
pub fn insert_acr(wire: &Wire, acr_cuts: &[f64], prec: f64) -> Result<Wire, String> {
    // `cxx:1072-1079`.
    let acr_wire = compute_acr(wire)?;

    let mut mw = crate::brep_lib_make_wire::MakeWire::new();
    let mut t1 = 0.0f64;
    let mut nb_edges = 0usize;

    // `cxx:1089-1147`: edge by edge.
    let mut wexp = crate::brep_tools_wire_explorer::WireExplorer::new();
    wexp.init(wire, None);
    while wexp.more() {
        nb_edges += 1;
        let t0 = t1;
        t1 = acr_wire[nb_edges];

        // `cxx:1097-1105`: the cuts strictly inside this edge's interval.
        let mut paradec: Vec<f64> = Vec::new();
        for &cut in acr_cuts {
            if t0 + prec < cut && cut < t1 - prec {
                paradec.push(cut);
            }
        }

        let e = match wexp.current_edge() {
            Some(e) => e,
            None => {
                wexp.next();
                continue;
            }
        };
        let v = wexp.current_vertex().cloned();

        if paradec.is_empty() || BRepTool::is_degenerated(&e) {
            // `cxx:1108-1112`: copy the edge.
            let _ = mw.add(&e);
        } else {
            // `cxx:1116-1140`.
            let so = match (&v, edge_vertices(&e).0) {
                (Some(v), Some(first)) => is_same(&v.0, &first.0),
                _ => false,
            };
            // "the wire is always FORWARD / it is necessary to modify the
            //  parameter of cut if the edge is REVERSED" (`cxx:1124-1138`).
            let sr: Vec<f64> = if e.orientation() == Orientation::Forward {
                paradec.clone()
            } else {
                paradec.iter().rev().map(|p| t0 + t1 - p).collect()
            };
            let se = crate::brep_fill_compatible_wires::statics::trim_edge(
                &e, &sr, t0, t1, so, false,
            );
            for piece in &se {
                let _ = mw.add(piece);
            }
        }
        wexp.next();
    }

    // `cxx:1148-1153`: the result keeps the input wire's orientation.
    let mut wres = mw.wire();
    wres.0.set_orientation(wire.0.orientation());
    Ok(wres)
}

/// `MakeWire(Edges, rangdeb, forward, newwire)` (`BRepFill.cxx:81-117`).
/// `rangdeb` is 1-based like OCCT's `NCollection_Array1<TopoDS_Shape>(1, n)`:
/// the `forward` arm appends `Edges(rangdeb..=n)` then `Edges(1..rangdeb)`, the
/// reverse arm prepends `Edges((1..=rangdeb).rev())` (already reversed edges)
/// then `Edges((rangdeb+1..=n).rev())`. Both arms finish with
/// `Orientation(FORWARD)` and `Closed(true)`.
fn make_wire_from(
    edges: &[TopoShape],
    rangdeb: usize,
    forward: bool,
    newwire: &mut Wire,
) -> Result<(), String> {
    let nb_edges = edges.len();
    if rangdeb == 0 || rangdeb > nb_edges {
        return Err("Standard_OutOfRange: BRepFill::SearchOrigin rangdeb".to_string());
    }
    let b = TopoBuilder::new();
    let mut ordered: Vec<Edge> = Vec::with_capacity(nb_edges);
    let edge_at = |i: usize| -> Edge { Edge(edges[i - 1].clone()) };
    if forward {
        // `for (rang = rangdeb; rang <= nbEdges; rang++)`
        for rang in rangdeb..=nb_edges {
            ordered.push(edge_at(rang));
        }
        // `for (rang = 1; rang < rangdeb; rang++)`
        for rang in 1..rangdeb {
            ordered.push(edge_at(rang));
        }
    } else {
        // `for (rang = rangdeb; rang >= 1; rang--) E.Reversed()`
        for rang in (1..=rangdeb).rev() {
            let mut s = edge_at(rang).0;
            s.reverse();
            ordered.push(Edge(s));
        }
        // `for (rang = nbEdges; rang > rangdeb; rang--) E.Reversed()`
        for rang in ((rangdeb + 1)..=nb_edges).rev() {
            let mut s = edge_at(rang).0;
            s.reverse();
            ordered.push(Edge(s));
        }
    }
    *newwire = b.make_wire(&ordered);
    // `newwire.Orientation(TopAbs_FORWARD); newwire.Closed(true);`
    newwire.0.set_orientation(Orientation::Forward);
    newwire.0.set_closed(true);
    Ok(())
}

/// `CutEdge(CurrentEdge, Param, E1, E2, VRef)` (`BRepFill.cxx:119-142`).
///
/// Splits `current_edge` at `param`; `E1` is the piece that starts at `VRef`
/// (i.e. the piece the wire already owns and is parsed first).
fn cut_edge_at(
    current_edge: &Edge,
    param: f64,
    v_ref: &Vertex,
) -> Result<(Edge, Edge), String> {
    let c = BRepTool::edge_curve_world(current_edge)
        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
    let (first, last) = BRepTool::edge_parameters(current_edge);
    // `B.MakeVertex(Vi, C->Value(Param), Precision::Confusion())` (`cxx:129`).
    let b = TopoBuilder::new();
    let vi = b.make_vertex(c.d0(param), CONFUSION);
    let (vf, vl) = edge_vertices(current_edge);
    let (Some(vf), Some(vl)) = (vf, vl) else {
        return Err("Standard_NullObject: TopExp::Vertices".to_string());
    };
    let piece = |v1: &Vertex, v2: &Vertex, p1: f64, p2: f64| -> Result<Edge, String> {
        let build = crate::brep_lib_make_edge::init_3d_vertices(&c, Some(v1), Some(v2), p1, p2)
            .ok_or_else(|| "BRepLib_EdgeDone".to_string())?;
        let mut e = b.make_edge(c.clone(), build.first, build.last);
        crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
        Ok(e)
    };
    if is_same(&v_ref.0, &vf.0) {
        // `E1 = MakeEdge(C, Vf, Vi, first, Param); E2 = MakeEdge(C, Vi, Vl, ...)`
        Ok((
            piece(&vf, &vi, first, param)?,
            piece(&vi, &vl, param, last)?,
        ))
    } else {
        Ok((
            piece(&vi, &vl, param, last)?,
            piece(&vf, &vi, first, param)?,
        ))
    }
}

/// One `BRepExtrema_DistShapeShape(V, W)` solution: the distance, whether the
/// support on the wire is a vertex, and the supporting sub-shape.
struct DistSolution {
    /// `P.Distance(DSS.PointOnShape2(isol))` (`cxx:908`): the solution's
    /// `PointOnShape2` is only ever read through this distance.
    distance: f64,
    is_vertex: bool,
    shape: TopoShape,
}

/// Re-derivation of `BRepExtrema_DistShapeShape(V, W)` for a vertex `V` and a
/// wire `W` (`BRepFill.cxx:904`). See the module notes: the two solution
/// families are the wire's vertices (vertex-vertex) and the wire's edges
/// (point-curve extrema). The minimum wins; a vertex solution wins a tie, which
/// is the order `BRepExtrema_ShapeList` builds for a wire (vertices before
/// edges).
fn dist_vertex_wire(p: &GpPnt, wire: &Wire) -> Option<DistSolution> {
    let mut best: Option<DistSolution> = None;
    for shape in crate::topo_tools::explorer_shapes(&wire.0, crate::abs::ShapeType::Vertex) {
        let v = Vertex(shape.clone());
        let point = BRepTool::vertex_point_world(&v);
        let distance = point.distance(p);
        let better = match &best {
            None => true,
            Some(b) => distance < b.distance,
        };
        if better {
            best = Some(DistSolution {
                distance,
                is_vertex: true,
                shape,
            });
        }
    }
    for shape in crate::topo_tools::explorer_shapes(&wire.0, crate::abs::ShapeType::Edge) {
        let e = Edge(shape.clone());
        let (u, point) = crate::brep_extrema::closest_point_on_edge(&e, p, 0);
        if !u.is_finite() {
            continue;
        }
        let distance = point.distance(p);
        let better = match &best {
            None => true,
            Some(b) => distance < b.distance,
        };
        if better {
            best = Some(DistSolution {
                distance,
                is_vertex: false,
                shape,
            });
        }
    }
    best
}

/// `BRepFill::SearchOrigin(W, P, Dir, Tol)` (`BRepFill.cxx:885-1023`): move the
/// seam of the closed wire `w` onto `p`, so the wire starts where the frame
/// search should begin, and orient the parse direction towards `dir`.
///
/// `w` must be closed; OCCT throws `Standard_NoSuchObject` otherwise. `p` is
/// projected onto the wire. When the projection lands inside an edge, that edge
/// is cut at the projected parameter (`CutEdge`), otherwise the vertex that
/// already sits at `p` becomes the seam.
pub fn search_origin(w: &mut Wire, p: &GpPnt, dir: &GpVec, tol: f64) -> Result<(), String> {
    if !w.0.is_closed() {
        // `throw Standard_NoSuchObject("BRepFill::SearchOrigin : the wire must
        //  be closed")` (`cxx:887-890`).
        return Err("Standard_NoSuchObject: BRepFill::SearchOrigin : the wire must be closed".into());
    }

    let mut new_vertex = false;
    let mut theparam = 1.0e101;
    let mut v: Option<Vertex> = None;
    let mut e_ref: Option<Edge> = None;

    // `W.Orientation(TopAbs_FORWARD)` (`cxx:900`).
    w.0.set_orientation(Orientation::Forward);

    // `B.MakeVertex(V, P, Tol); BRepExtrema_DistShapeShape DSS(V, W)`
    // (`cxx:903-904`).
    if let Some(sol) = dist_vertex_wire(p, w) {
        // `cxx:907-916`: the loop only keeps the minimum solution.
        let _dss = sol.distance;
        if sol.is_vertex {
            // `V = TopoDS::Vertex(supp)` (`cxx:918-921`).
            v = Some(Vertex(sol.shape));
        } else {
            // `cxx:923-946`.
            let e = Edge(sol.shape);
            let (vf, vl) = edge_vertices(&e);
            let (Some(vf), Some(vl)) = (vf, vl) else {
                return Err("Standard_NullObject: TopExp::Vertices".to_string());
            };
            let dist = p.distance(&BRepTool::vertex_point_world(&vf));
            if dist < tol {
                v = Some(vl.clone());
            }
            let d = p.distance(&BRepTool::vertex_point_world(&vl));
            if d < tol && d < dist {
                v = Some(vf.clone());
                // `dist = d;` -- OCCT keeps this local, it feeds `NewVertex`.
            }
            let dist = if d < tol && d < dist { d } else { dist };
            new_vertex = dist > tol;
            if new_vertex {
                // `DSS.ParOnEdgeS2(isol, theparam)` (`cxx:944`).
                let (u, _) = crate::brep_extrema::closest_point_on_edge(&e, p, 0);
                theparam = u;
            }
            e_ref = Some(e);
        }
    }

    // `cxx:955-968`: count the edges; a new vertex adds one (the cut).
    let mut wexp = crate::brep_tools_wire_explorer::WireExplorer::new();
    wexp.init(w, None);
    let mut nb_edges = 0usize;
    while wexp.more() {
        nb_edges += 1;
        wexp.next();
    }
    if new_vertex {
        nb_edges += 1;
    }

    // `cxx:971-992`: the edge table and `rangdeb`.
    let mut rangdeb = 0usize;
    let mut edges: Vec<TopoShape> = Vec::with_capacity(nb_edges);
    wexp.init(w, None);
    let mut ii = 1usize;
    while wexp.more() {
        let Some(e) = wexp.current_edge() else {
            wexp.next();
            continue;
        };
        let split = match &e_ref {
            Some(r) => new_vertex && is_same(&e.0, &r.0),
            None => false,
        };
        if split {
            let current_vertex = wexp.current_vertex().cloned();
            let Some(current_vertex) = current_vertex else {
                return Err("Standard_NullObject: BRepTools_WireExplorer::CurrentVertex".into());
            };
            let (e1, e2) = cut_edge_at(&e, theparam, &current_vertex)?;
            if ii > nb_edges {
                return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges".to_string());
            }
            edges.push(e1.0.clone());
            ii += 1;
            if ii > nb_edges {
                return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges".to_string());
            }
            edges.push(e2.0.clone());
            rangdeb = ii;
        } else {
            if ii > nb_edges {
                return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges".to_string());
            }
            edges.push(e.0.clone());
        }
        if !new_vertex {
            if let (Some(v), Some(cv)) = (&v, wexp.current_vertex()) {
                if is_same(&v.0, &cv.0) {
                    rangdeb = ii;
                }
            }
        }
        ii += 1;
        wexp.next();
    }
    if rangdeb == 0 {
        rangdeb = nb_edges;
    }

    // `cxx:998-1018`: the direction of parsing.
    if rangdeb == 0 || rangdeb > edges.len() {
        return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges(rangdeb)".to_string());
    }
    let e = Edge(edges[rangdeb - 1].clone());
    if !new_vertex {
        // `theparam = BRep_Tool::Parameter(V, E)` (`cxx:1003`).
        let Some(v) = &v else {
            return Err("Standard_NullObject: BRep_Tool::Parameter".to_string());
        };
        theparam = BRepTool::parameter_of_vertex(&e.0, &v.0);
    }
    // `BRepAdaptor_Curve AC(E); AC.D1(theparam, Pe, Ve)` (`cxx:1005-1008`).
    let curve = BRepTool::edge_curve_world(&e)
        .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
    let (_pe, mut ve) = curve.d1(theparam);
    if e.0.orientation() == Orientation::Reversed {
        ve = ve.multiplied_scalar(-1.0);
    }
    let mut angle = ve.angle(dir);
    if angle > std::f64::consts::PI {
        angle = 2.0 * std::f64::consts::PI - angle;
    }
    let forward = angle <= std::f64::consts::PI / 2.0;

    // `MakeWire(Edges, rangdeb, forward, W); W.Closed(true)` (`cxx:1020-1022`).
    make_wire_from(&edges, rangdeb, forward, w)?;
    Ok(())
}

/// `BRepFill::Axe` (`BRepFill.cxx:672-878`).
///
/// Computes the coordinate system `AxeProf` of an evolved sweep: the normal of
/// the spine plane, and the point / tangent of the spine closest to the
/// profile. `ProfOnSpine` reports whether that closest point is within `tol`
/// (`cxx:874`), which the caller uses to reject an inconsistent
/// `ProfOnSpine` request.
///
/// Translation notes:
/// - `BRep_Tool::Surface(Face, L)` + `S->Transformed(L)` (`cxx:690-719`) is
///   [`BRepTool::face_surface`] plus an explicit composition of the
///   `BRepLib_FindSurface` location with the face location, the same way
///   [`crate::brep_fill_evolved::perform::FindLocation`] does it.
/// - `BRepExtrema_ExtPC::Perform(V)` + `IsMin(i)` (`cxx:845-857`) is
///   [`crate::brep_extrema::closest_point_on_edge`], which returns the
///   minimal solution of the point/curve extrema.
pub fn axe(spine: &TopoShape, profile: &Wire, tol: f64) -> Result<(GpAx3, bool), String> {
    // `cxx:684-709`: the face and its plane.
    let mut a_face = Face(TopoShape::new(ShapeType::Face));
    let mut surface: Option<Arc<dyn Surface>> = None;
    let mut location = TopLocLocation::identity();
    if spine.shape_type() == ShapeType::Face {
        let face = Face(spine.clone());
        location = face.0.location().clone();
        surface = BRepTool::face_surface(&face);
        // `if (!S->IsKind(STANDARD_TYPE(Geom_Plane)))` (`cxx:691`): the exact
        // type test, so a trimmed surface of a plane takes the
        // `BRepLib_FindSurface` branch.
        let is_plain_plane = surface
            .as_ref()
            .is_some_and(|s| s.rectangular_trimmed_basis().is_none() && s.gp_pln().is_some());
        if !is_plain_plane {
            // `cxx:693-702`.
            let fs = FindSurface::new(&face.0, -1.0, true, false);
            if !fs.is_found() {
                return Err(
                    "Standard_NoSuchObject: BRepFill_Evolved : The Face is not planar".into(),
                );
            }
            surface = fs.surface();
            location = TopLocLocation::composed(
                &fs.location().transformation(),
                face.0.location().clone(),
            );
        }
        a_face = face;
    } else if spine.shape_type() == ShapeType::Wire {
        // `cxx:705-709`.
        a_face = make_face_from_wire(&Wire(spine.clone()), true)
            .ok_or_else(|| "BRepLib_MakeFace: the wire is not planar".to_string())?;
        location = a_face.0.location().clone();
        surface = BRepTool::face_surface(&a_face);
    }

    // `cxx:711-714`.
    let mut surface = match surface {
        Some(s) => s,
        None => return Err("Standard_DomainError: BRepFill_Evolved::Axe".into()),
    };

    // `cxx:716-719`.
    if !location.is_identity() {
        surface = Arc::from(surface.transformed(&location.transformation()));
    }

    // `cxx:721`: `down_cast<Geom_Plane>(S)->Pln().Axis().Direction()`.
    let normal = surface
        .gp_pln()
        .map(|p| {
            let ax1 = p.axis();
            *ax1.direction()
        })
        .ok_or_else(|| "Standard_DomainError: BRepFill_Evolved::Axe".to_string())?;

    // `cxx:724-727`.
    let mut dist_min = INFINITE;
    // `double Tol2 = 1.e-10;` (`cxx:727`).
    let tol2 = 1.0e-10;

    // `cxx:734-760`: is there a vertex/vertex contact between the spine and the
    // profile? Note the loop does not stop at the first spine vertex in OCCT:
    // `SE.Next()` runs on every outer iteration, and the break only fires once
    // a contact is found (`cxx:739-760`).
    let forward_face = a_face.0.oriented(Orientation::Forward);
    let spine_vertices = explorer_shapes(&forward_face, ShapeType::Vertex);
    let profile_vertices = explorer_shapes(&profile.0, ShapeType::Vertex);

    let mut on_vertex: Option<Vertex> = None;
    for sv in &spine_vertices {
        let p1 = BRepTool::vertex_point_world(&Vertex(sv.clone()));
        let mut found = false;
        for pv in &profile_vertices {
            let p2 = BRepTool::vertex_point_world(&Vertex(pv.clone()));
            if p1.square_distance(&p2) <= tol2 {
                found = true;
                break;
            }
        }
        if found {
            on_vertex = Some(Vertex(sv.clone()));
            break;
        }
    }

    // `cxx:730`: `Par` only feeds the edge branch.
    let mut par = 0.0;
    let mut loc = GpPnt::new(0.0, 0.0, 0.0);
    let mut tang = GpVec::from_pnts(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 0.0, 0.0));

    if let Some(vonf) = on_vertex {
        // `cxx:762-833`.
        let map = map_shapes_and_ancestors(&forward_face, ShapeType::Vertex, ShapeType::Edge);
        let list = map
            .get(&shape_map_key(&vonf.0))
            .cloned()
            .ok_or_else(|| "Standard_NoSuchObject: TopExp_AncestorMap".to_string())?;
        let e1 = list
            .first()
            .map(|s| Edge(s.clone()))
            .ok_or_else(|| "Standard_NoSuchObject: TopExp_AncestorMap".to_string())?;
        let e2 = list
            .last()
            .map(|s| Edge(s.clone()))
            .ok_or_else(|| "Standard_NoSuchObject: TopExp_AncestorMap".to_string())?;

        // `CE1 = BRep_Tool::Curve(E1, L, f, l)` + `D1` (`cxx:778-789`).
        let (loc1, tang1) = edge_point_and_tangent(&e1, &vonf, &a_face)?;
        // `cxx:791-802`.
        let (loc2, tang2) = edge_point_and_tangent(&e2, &vonf, &a_face)?;

        // `cxx:804-820`: sum the |dot| of both tangents over the profile edges.
        let mut sca1 = 0.0;
        let mut sca2 = 0.0;
        for pe in explorer_shapes(&profile.0, ShapeType::Edge) {
            let e = Edge(pe);
            let (v1, v2) = edge_vertices(&e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                continue;
            };
            let p1 = BRepTool::vertex_point_world(&v1);
            let p2 = BRepTool::vertex_point_world(&v2);
            let vec = GpVec::from_pnts(&p1, &p2);
            sca1 += tang1.dot(&vec).abs();
            sca2 += tang2.dot(&vec).abs();
        }

        // `cxx:822-831`.
        if sca1.abs() < sca2.abs() {
            loc = loc1;
            tang = tang1;
        } else {
            loc = loc2;
            tang = tang2;
        }
        dist_min = 0.0;
    } else {
        // `cxx:834-872`.
        for se in explorer_shapes(&forward_face, ShapeType::Edge) {
            let e = Edge(se);
            for pv in &profile_vertices {
                let v = Vertex(pv.clone());
                let p2 = BRepTool::vertex_point_world(&v);
                // `BE.Perform(V)` (`cxx:845`): the minimal solution.
                let (u, _pt) = crate::brep_extrema::closest_point_on_edge(&e, &p2, 0);
                let dist = p2.distance(&evaluate_edge(&e, u)?);
                // `cxx:860-869`.
                if dist < dist_min {
                    dist_min = dist;
                    par = u;
                    // `BRepAdaptor_Curve BAC(E); BAC.D1(Par, Loc, Tang)`.
                    let curve = BRepTool::edge_curve_adaptor_world(&e)
                        .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
                    let (p, mut t) = curve.d1(par);
                    if e.0.orientation() == Orientation::Reversed {
                        t.reverse();
                    }
                    loc = p;
                    tang = t;
                }
            }
        }
    }

    // `cxx:874`.
    let prof_on_spine = dist_min < tol;

    // `cxx:876-877`: `gp_Ax3 A3(Loc, Normal, Tang)`.
    let tang_dir = GpDir::from_vec(&tang)
        .map_err(|_| "Standard_ConstructionError: gp_Ax3".to_string())?;
    let axe_prof = GpAx3::new(loc, normal, &tang_dir)
        .map_err(|_| "Standard_ConstructionError: gp_Ax3".to_string())?;
    Ok((axe_prof, prof_on_spine))
}

/// Evaluates the 3D curve of `e` at `u` in world coordinates
/// (`BRepAdaptor_Curve::D0`).
fn evaluate_edge(e: &Edge, u: f64) -> Result<GpPnt, String> {
    let curve = BRepTool::edge_curve_world(e)
        .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
    Ok(curve.d0(u))
}

/// `CE1 = BRep_Tool::Curve(E, L, f, l); CE1->D1(Par, Loc, Tang)` with the
/// location applied and the tangent reversed for a `REVERSED` edge
/// (`BRepFill.cxx:778-789`, `:791-802`).
fn edge_point_and_tangent(
    e: &Edge,
    v: &Vertex,
    face: &Face,
) -> Result<(GpPnt, GpVec), String> {
    let curve = BRepTool::edge_curve_world(e)
        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
    let par = BRepTool::parameter_of_vertex_on_face(&v.0, e, face);
    let (p, mut t) = curve.d1(par);
    if e.0.orientation() == Orientation::Reversed {
        t.reverse();
    }
    // `Tang1.Normalize()` (`cxx:805-806`).
    t.normalize();
    Ok((p, t))
}

