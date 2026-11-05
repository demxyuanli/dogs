//! `TrimFace` (`BRepFill_Evolved.cxx:2647-2711`) and its static
//! `AddDegeneratedEdge` (`cxx:2555-2645`).
//!
//! `TrimFace` chains the trimmed edges into wires with `BRepLib_MakeWire` and
//! turns every wire into a face over a copy of the trimmed surface.
//! `AddDegeneratedEdge` then inserts the degenerate edges that `TrimEdge` left
//! out: with the `yOz` working surface every wire must be closed on the
//! surface, so any UV gap between consecutive edges becomes a degenerate edge.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::boptools_2d;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::is_same;

use super::side::edge_vertices_oriented;

/// `TolConf` of `AddDegeneratedEdge` (`cxx:2573`).
const TOL_CONF: f64 = 1.0e-4;

/// `BRepLib_MakeEdge(C2d, S)` (`BRepLib_MakeEdge.cxx:217-262`): an edge whose
/// 3D curve is the p-curve lifted to the surface (`Adaptor3d_CurveOnSurface`).
///
/// The reference builds it over the *whole* range of `C2d` and re-trims with
/// `B.Range`; the port passes the final range straight to `Init` because the
/// port's `init_2d` rejects infinite parameters, which a `Geom2d_Line` always
/// has. The resulting curve, p-curve and range are the same.
fn make_edge_on_surface(
    pc: &Arc<dyn Curve2d>,
    surface: &Arc<dyn Surface>,
    face: &Face,
    first: f64,
    last: f64,
) -> Option<Edge> {
    let build = brep_lib_make_edge::init_2d(pc, surface, first, last)?;
    let curve3d: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(pc.clone(), surface.clone()));
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(curve3d, build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    let _ = boptools_2d::attach_existing_pcurve(&edge, face, pc.clone());
    Some(edge)
}

/// The degenerate edge of `cxx:2605-2618` / `cxx:2626-2640`: the segment
/// `a`..`b` in UV space, degenerated onto the vertex `v`.
fn add_degenerate_segment(
    face: &Face,
    surface: &Arc<dyn Surface>,
    wire: &mut TopoShape,
    a: &GpPnt2d,
    b: &GpPnt2d,
    v: &Vertex,
) {
    let vec = GpVec2d::new(b.x() - a.x(), b.y() - a.y());
    let Ok(dir) = GpDir2d::new(vec.x(), vec.y()) else {
        // `gp_Dir2d(V)` with a null vector raises in OCCT.
        return;
    };
    let line2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(*a, dir));
    let f = 0.0;
    let l = a.distance(b);
    let Some(edge) = make_edge_on_surface(&line2d, surface, face, f, l) else {
        return;
    };

    let builder = TopoBuilder::new();
    // `B.Degenerated(NE, true)`.
    GeometryRegistry::global().set_degenerated(&edge.0, true);
    let mut ne = edge.0.clone();
    // `B.Add(NE, V1.Oriented(TopAbs_FORWARD))` then the same vertex REVERSED.
    builder.add(&mut ne, &v.0.oriented(Orientation::Forward));
    builder.add(&mut ne, &v.0.oriented(Orientation::Reversed));
    GeometryRegistry::global().set_ranges_all(&ne, f, l);
    builder.add(wire, &ne);
}

/// `AddDegeneratedEdge(F, W)` (`BRepFill_Evolved.cxx:2555-2645`).
pub fn add_degenerated_edge(face: &Face, wire: &mut TopoShape, surface: &Arc<dyn Surface>) {
    // `cxx:2558-2572`: a plane (possibly behind a rectangular trim) needs no
    // degenerate edges.
    if let Some(basis) = surface.rectangular_trimmed_basis() {
        if basis.gp_pln().is_some() {
            return;
        }
    }
    if surface.gp_pln().is_some() {
        return;
    }

    let mut change = true;
    while change {
        change = false;
        let wire_view = Wire(wire.clone());
        let mut explorer = WireExplorer::of_wire_on_face(&wire_view, face);

        let mut vf: Option<Vertex> = None;
        let mut pf = GpPnt2d::new(0.0, 0.0);
        let mut prev_p = GpPnt2d::new(0.0, 0.0);
        let mut last_v2: Option<Vertex> = None;
        let mut last_p2 = GpPnt2d::new(0.0, 0.0);

        while explorer.more() {
            let Some(ce) = explorer.current_edge() else {
                break;
            };
            let (v1, v2) = edge_vertices_oriented(&ce);
            let (p_general_1, p_general_2) = BRepTool::uv_points(&ce, face);
            let (p1, p2) = if ce.0.orientation() == Orientation::Reversed {
                (p_general_2, p_general_1)
            } else {
                (p_general_1, p_general_2)
            };

            match vf.clone() {
                None => {
                    vf = v1.clone();
                    pf = p1;
                }
                Some(_) => {
                    if !p1.is_equal(&prev_p, TOL_CONF) {
                        // `cxx:2603-2618`: insert the missing degenerate edge.
                        change = true;
                        if let Some(v) = v1.as_ref() {
                            add_degenerate_segment(face, surface, wire, &prev_p, &p1, v);
                        }
                        break;
                    }
                }
            }
            prev_p = p2;
            last_v2 = v2;
            last_p2 = p2;
            explorer.next();
        }

        // `cxx:2626-2641`: the closed wire case, from the last UV point back
        // to the first one.
        if !change {
            if let (Some(vf), Some(v2)) = (vf.as_ref(), last_v2.as_ref()) {
                if is_same(&vf.0, &v2.0) && !pf.is_equal(&last_p2, TOL_CONF) {
                    change = true;
                    add_degenerate_segment(face, surface, wire, &last_p2, &pf, vf);
                }
            }
        }
    }
}

/// `TrimFace(Face, TheEdges, S)` (`BRepFill_Evolved.cxx:2647-2711`).
///
/// `TheEdges` is consumed; `S` receives one face per reconstructed wire.
pub fn trim_face(face: &Face, the_edges: &mut Vec<TopoShape>) -> Vec<TopoShape> {
    let builder = TopoBuilder::new();
    let mut out: Vec<TopoShape> = Vec::new();
    let Some(surface) = BRepTool::face_surface(face) else {
        return out;
    };
    let tol_f = BRepTool::face_tolerance(face);

    // `cxx:2663`: `while (!TheEdges.IsEmpty())`.
    while !the_edges.is_empty() {
        // `cxx:2666-2669`: `BRepLib_MakeWire MWire(first)`, `GoodWire = ...`.
        let first = Edge(the_edges[0].clone());
        let mut mwire = MakeWire::new();
        let _ = mwire.add(&first);
        let mut good_wire = mwire.wire().0;
        the_edges.remove(0);

        let mut nb_edges = the_edges.len();
        let mut new_wire = false;

        // `cxx:2671-2701`: grow the wire while an edge connects.
        while !new_wire {
            let mut add_edge = false;
            let mut i = 0;
            while i < nb_edges && !add_edge {
                let e = Edge(the_edges[i].clone());
                if BRepTool::is_degenerated(&e) {
                    // `cxx:2679-2685`.
                    the_edges.remove(i);
                    add_edge = true;
                    nb_edges = the_edges.len();
                    good_wire = mwire.wire().0;
                } else {
                    // `cxx:2687-2697`.
                    // `cxx:2687-2697`: `Add` returns `void` and only sets
                    // `myError`, which the `Error()` test below reads.
                    let _ = mwire.add(&e);
                    if mwire.error() == crate::brep_lib_make_wire::WireError::WireDone {
                        the_edges.remove(i);
                        add_edge = true;
                        nb_edges = the_edges.len();
                        good_wire = mwire.wire().0;
                    }
                }
                i += 1;
            }
            new_wire = !add_edge;
        }

        // `cxx:2703-2709`: the face over the trimmed surface.
        let mut face_cut = builder.make_face_with_tolerance(surface.clone(), &[], tol_f);
        face_cut.0.set_location(face.0.location());
        face_cut.0.set_orientation(Orientation::Forward);
        // `BRepTools::Update(FaceCut)` (`BRepTools.cxx:...)` only refreshes the
        // cached UV bounds of `BRep_TFace`; the port computes them on demand,
        // so there is nothing to update.
        add_degenerated_edge(&face_cut, &mut good_wire, &surface);
        builder.add(&mut face_cut.0, &good_wire);
        face_cut.0.set_orientation(face.0.orientation());
        out.push(face_cut.0);
    }

    out
}
