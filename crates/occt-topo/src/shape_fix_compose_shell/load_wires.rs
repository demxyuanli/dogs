//! Port of ShapeFix_ComposeShell::LoadWires (ShapeFix_ComposeShell.cxx:499-640).
//!
//! ShapeExtend_WireData is represented by an ordered Vec of Edge throughout
//! this module, so the OCCT wire-data operations below are plain vector work.

use crate::abs::Orientation;
use crate::shape::{Face, TopoShape, Vertex, Wire};

use super::wire_segment::WireSegment;

/// ShapeExtend_WireData(wire, chkseam, manifold) (ShapeExtend_WireData.cxx:80-120):
/// the stored edges keep the composed orientation of the wire, so a REVERSED
/// wire flips every edge.
fn wire_data_edges(wire: &Wire) -> Vec<crate::shape::Edge> {
    let edges = crate::topo_tools_full::edges_of_wire(wire);
    if wire.0.orientation() == Orientation::Reversed {
        edges
            .into_iter()
            .map(|mut e| {
                e.0.reverse();
                e
            })
            .collect()
    } else {
        edges
    }
}

/// ShapeFix_ComposeShell::LoadWires (cxx:499-640).
///
/// Context Apply (cxx:506) runs through the ComposeShell context (a
/// `MapReShape`; see `reshape.rs`). The ShapeFix_Wire::FixReorder block
/// (cxx:582-631) is ported through `shhealing::fix_reorder_wire*`.
pub fn load_wires(face: &Face) -> Vec<WireSegment> {
    let mut seqw: Vec<WireSegment> = Vec::new();
    let children: Vec<TopoShape> = {
        let ts = face.0.tshape.read().expect("poisoned TShape lock");
        ts.children.clone()
    };
    for child in children {
        if !child.is_wire() {
            if child.is_vertex() {
                let mut seg = WireSegment::new();
                seg.set_vertex(Some(Vertex(child.clone())));
                seg.set_orientation(child.orientation());
                seqw.push(seg);
            }
            continue;
        }
        let wire = Wire(child);
        let wo = wire.0.orientation();
        let is_non_manifold = wo != Orientation::Reversed && wo != Orientation::Forward;

        if is_non_manifold {
            let sbwd = wire_data_edges(&wire);
            if !sbwd.is_empty() {
                // non-manifold wires take INTERNAL orientation (cxx:541-543).
                seqw.push(WireSegment::with_edges(sbwd, Orientation::Internal));
            }
        } else {
            let mut sbwd_m: Vec<crate::shape::Edge> = Vec::new();
            let mut sbwd_nm: Vec<crate::shape::Edge> = Vec::new();
            for e in wire_data_edges(&wire) {
                let eo = e.0.orientation();
                if eo == Orientation::Forward || eo == Orientation::Reversed {
                    sbwd_m.push(e);
                } else {
                    sbwd_nm.push(e);
                }
            }
            if !sbwd_nm.is_empty() {
                seqw.push(WireSegment::with_edges(sbwd_nm, Orientation::Internal));
            }
            if !sbwd_m.is_empty() {
                // cxx:582-631: reorder the manifold part (FixReorder) and
                // correct its direction when the reorder flips IsOuterBound.
                let b = crate::builder::TopoBuilder::new();
                let wire_m = b.make_wire(&sbwd_m);
                let mut stat = 0i32;
                if let Some(surf) = crate::brep_tool::BRepTool::face_surface(face) {
                    if surf.is_u_periodic() && surf.is_v_periodic() {
                        // cxx:586-607: for torus-like shapes reorder in 2d first.
                        let face_fwd = Face(face.0.oriented(Orientation::Forward));
                        let mut sawo = crate::meshing::wire_order::WireOrder::new();
                        for e in &sbwd_m {
                            if let Some((c2d, f, l)) =
                                crate::boptools_2d::curve_on_surface_oriented(e, &face_fwd, true)
                            {
                                sawo.add_edge(c2d.d0(f), c2d.d0(l)); // cxx:601
                            }
                        }
                        sawo.perform();
                        stat = if matches!(
                            sawo.status(),
                            crate::meshing::wire_order::WireOrderStatus::Reversed
                        ) {
                            -1
                        } else {
                            1
                        };
                        let _ = crate::shhealing::fix_reorder_wire_with_order(&wire_m, &sawo);
                    }
                }
                // cxx:609: sfw->FixReorder() (3d).
                let (_ok, status3d) = crate::shhealing::fix_reorder_wire_3d(&wire_m);
                if matches!(
                    status3d,
                    crate::meshing::wire_order::WireOrderStatus::Reversed
                ) {
                    stat = -1; // cxx:610-613: StatusReorder(DONE3)
                }
                let mut edges_final: Vec<crate::shape::Edge> =
                    crate::topo_tools_full::edges_of_wire(&wire_m);
                if stat < 0 {
                    // cxx:615-631: reverse only when IsOuterBound flips.
                    if let Some(surf) = crate::brep_tool::BRepTool::face_surface(face) {
                        let before = crate::shape_analysis::is_outer_bound(
                            &b.make_face(surf.clone(), &[wire.clone()]),
                        );
                        let w = b.make_wire(&edges_final);
                        let after = crate::shape_analysis::is_outer_bound(&b.make_face(surf, &[w]));
                        if before != after {
                            super::wire_data::reverse_wire_data_on_face(&mut edges_final, face);
                        }
                    }
                }
                // cxx:634: ShapeFix_WireSegment seg(sbwdM, TopAbs_REVERSED).
                seqw.push(WireSegment::with_edges(edges_final, Orientation::Reversed));
            }
        }
    }
    seqw
}
