//! `BRepFill_CompatibleWires::SearchOrigin` (`cxx:2453-2629`).
//!
//! The open-wire branch: every section after the first is re-stacked either
//! forwards or backwards so that its first vertex lands on the side of the
//! previous section's `(Pdeb, Pfin)` chord that keeps the ruled surface
//! untwisted.
//!
//! The decision is a chord comparison for straight (non-planar) sections and a
//! chord *direction* comparison for planar ones; a section whose two boundary
//! vertices coincide on either side is nudged by `Precision::Confusion` along
//! its first curve before the directions are compared (`cxx:2572-2578`).
//!
//! Not tracked: `P0` (`cxx:2457`, assigned at `cxx:2617`) is only read by the
//! commented-out `cxx:2540-2571` block, so the plane of the previous section is
//! not carried.
//!
//! Reproduced as-is: `newwire.Oriented(TopAbs_FORWARD)` (`cxx:2604`) discards
//! the returned shape, so it does not change `newwire`; a `BRep_Builder` wire
//! is already `TopAbs_FORWARD`.

use occt_core::gp::GpVec;
use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::brep_fill_compatible_wires::statics::plane_of_wire_full;
use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_generator::ThruSectionErrorStatus;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Wire};
use crate::topo_tools::{is_same, wire_vertices};

/// `SearchOrigin()` (`cxx:2453-2629`).
pub(crate) fn search_origin(cw: &mut CompatibleWires) -> Result<(), String> {
    // `cxx:2465-2473`.
    let mut ideb = 1usize;
    let mut ifin = cw.my_work.len();
    if cw.my_degen1 {
        ideb += 1;
    }
    if cw.my_degen2 {
        ifin -= 1;
    }
    let v_closed = (!cw.my_degen1)
        && (!cw.my_degen2)
        && ideb <= ifin
        && is_same(&cw.my_work[ideb - 1], &cw.my_work[ifin - 1]);

    // `cxx:2478-2486`: every section must be open.
    let mut all_open = true;
    if ideb <= ifin {
        for i in ideb..=ifin {
            all_open = all_open && !cw.my_work[i - 1].is_closed();
        }
    }
    if !all_open {
        cw.my_status = ThruSectionErrorStatus::NotSameTopology;
        return Ok(());
    }

    // `cxx:2490-2497`.
    if ideb > ifin {
        return Ok(());
    }
    let mut wire1 = Wire(cw.my_work[ideb - 1].clone());
    wire1.0.set_orientation(Orientation::Forward);
    let (Some(vdeb0), Some(vfin0)) = wire_vertices(&wire1) else {
        cw.my_status = ThruSectionErrorStatus::Failed;
        return Ok(());
    };
    let mut pdeb = BRepTool::vertex_point_world(&vdeb0);
    let mut pfin = BRepTool::vertex_point_world(&vfin0);
    let (_, is_plane0) = plane_of_wire_full(&wire1);
    let mut isline0 = !is_plane0;
    cw.my_work[ideb - 1] = wire1.0.clone();
    // `cxx:2498-2499`.
    let mut e0: Option<Edge> = wire_first_edge(&wire1);

    for i in (ideb + 1)..=ifin {
        let mut wire = Wire(cw.my_work[i - 1].clone());
        wire.0.set_orientation(Orientation::Forward);

        // `cxx:2506-2516`: the edges in explorer order.
        let seq_edges = wire_edges(&wire);
        let (Some(vdeb), Some(vfin)) = wire_vertices(&wire) else {
            cw.my_status = ThruSectionErrorStatus::Failed;
            return Ok(());
        };
        let (_, is_plane) = plane_of_wire_full(&wire);
        let isline = !is_plane;

        let parcours = if isline0 || isline {
            // `cxx:2528-2534`: straight segments, compare the two chords.
            let p1 = BRepTool::vertex_point_world(&vdeb);
            let p2 = BRepTool::vertex_point_world(&vfin);
            let dist1 = pdeb.distance(&p1) + pfin.distance(&p2);
            let dist2 = pdeb.distance(&p2) + pfin.distance(&p1);
            dist2 >= dist1
        } else {
            // `cxx:2538-2579`.
            let p1 = BRepTool::vertex_point_world(&vdeb);
            let p1o = pdeb;
            let mut p2 = BRepTool::vertex_point_world(&vfin);
            let mut p2o = pfin;
            if p1.distance(&p2) <= CONFUSION || p1o.distance(&p2o) <= CONFUSION {
                // `Curve0.D0(Curve0.FirstParameter() + Confusion, P2o)`.
                if let Some(e) = e0.as_ref() {
                    if let Some(c0) = BRepTool::edge_curve_adaptor_world(e) {
                        let (first, _) = BRepTool::edge_parameters(e);
                        p2o = c0.d0(first + CONFUSION);
                    }
                }
                if let Some(e) = seq_edges.first() {
                    if let Some(c) = BRepTool::edge_curve_adaptor_world(e) {
                        let (first, _) = BRepTool::edge_parameters(e);
                        p2 = c.d0(first + CONFUSION);
                    }
                }
            }
            let v_deb_fin0 = GpVec::from_pnts(&p1o, &p2o);
            let v_deb_fin = GpVec::from_pnts(&p1, &p2);
            let a_straight = v_deb_fin0.angle(&v_deb_fin);
            a_straight < std::f64::consts::FRAC_PI_2
        };

        // `cxx:2582-2601`: re-stack the wire.
        let mut new_wire = Wire::new();
        let builder = TopoBuilder::new();
        if parcours {
            for e in &seq_edges {
                builder.add(&mut new_wire.0, &e.0);
            }
        } else {
            for e in seq_edges.iter().rev() {
                let mut er = e.clone();
                er.0.reverse();
                builder.add(&mut new_wire.0, &er.0);
            }
        }
        // `cxx:2604`: the result of `Oriented` is discarded by OCCT.
        cw.my_work[i - 1] = new_wire.0;

        // `cxx:2608-2619`.
        if parcours {
            pdeb = BRepTool::vertex_point_world(&vdeb);
            pfin = BRepTool::vertex_point_world(&vfin);
        } else {
            pfin = BRepTool::vertex_point_world(&vdeb);
            pdeb = BRepTool::vertex_point_world(&vfin);
        }
        isline0 = isline;
        e0 = seq_edges.first().cloned();
    }

    // `cxx:2622-2626`.
    if v_closed {
        let last = cw.my_work.len();
        cw.my_work[last - 1] = cw.my_work[0].clone();
    }

    Ok(())
}

/// `for (anExp.Init(w); anExp.More(); anExp.Next()) SeqEdges.Append(anExp.Current())`
/// (`cxx:2511-2515`).
fn wire_edges(w: &Wire) -> Vec<Edge> {
    let mut out = Vec::new();
    let mut exp = WireExplorer::new();
    exp.init(w, None);
    while exp.more() {
        if let Some(e) = exp.current_edge() {
            out.push(e);
        }
        exp.next();
    }
    out
}

/// `anExp.Init(w); anExp.Current()` (`cxx:2498-2499`).
fn wire_first_edge(w: &Wire) -> Option<Edge> {
    let mut exp = WireExplorer::new();
    exp.init(w, None);
    exp.current_edge()
}
