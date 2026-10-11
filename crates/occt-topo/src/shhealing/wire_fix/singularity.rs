//! Singularity split inside `ShapeFix_Wire::FixEdgeCurves`
//! (`ShapeFix_Wire.cxx:673-863`).
//!
//! The arm runs only when `FixAddPCurve` reports `ShapeExtend_DONE2`, which
//! is the projector's `ShapeExtend_DONE4` (`ShapeFix_Edge.cxx:536-538`): the
//! sampled pcurve jumped half a period over a singularity. A missing 3D
//! curve or a projection that is not `DONE4` leaves the edge whole.

use occt_core::precision::PCONFUSION;
use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::pcurve_full::{sa_is_u_closed, SurfaceSingularities};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::shhealing::shape_analysis_curve::project_adaptor_bounds;
use crate::shhealing::transfer_params::TransferParameters;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::edges_of_wire;

use super::pcurve_range::fix_add_pcurve;
use super::pcurve_range::is_seam_use;
use super::same_parameter::copy_replace_vertices_with;
use super::small::SHAPE_FIX_MIN_TOLERANCE;
use super::{first_vertex, last_vertex};

/// `ShapeFix_Wire.cxx:673-822`. Returns true when the edge was replaced, so
/// the caller repeats `FixAddPCurve` on the first new piece (`i--`).
pub(in crate::shhealing) fn split_edge_over_singularity(
    wire: &mut Wire,
    face: &Face,
    edge: &Edge,
    index: usize,
    ctx: &mut SharedReShape,
) -> bool {
    let reg = GeometryRegistry::global();
    let Some(curve) = reg.edge_curve(&edge.0) else {
        return false;
    };
    let (a0, b0) = reg.edge_parameters(&edge.0);
    if !a0.is_finite() || !b0.is_finite() {
        return false;
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let sings = SurfaceSingularities::compute(surf.as_ref());
    let mut seq: Vec<f64> = Vec::new();
    for (preci, p3d) in sings.points_within(SHAPE_FIX_MIN_TOLERANCE) {
        let proj = project_adaptor_bounds(
            curve.as_ref(),
            &p3d,
            SHAPE_FIX_MIN_TOLERANCE,
            a0,
            b0,
            true,
        );
        if proj.distance >= preci.max(SHAPE_FIX_MIN_TOLERANCE) {
            continue;
        }
        let split = proj.param;
        if split - a0 <= PCONFUSION || b0 - split <= PCONFUSION {
            continue;
        }
        let mut placed = false;
        for (k, existing) in seq.iter().enumerate() {
            if split < existing - PCONFUSION {
                seq.insert(k, split);
                placed = true;
                break;
            } else if split < existing + PCONFUSION {
                placed = true;
                break;
            }
        }
        if !placed {
            seq.push(split);
        }
    }
    if seq.is_empty() {
        return false;
    }

    let is_fwd = edge.0.orientation() == Orientation::Forward;
    let mut work = edge.clone();
    work.0.set_orientation(Orientation::Forward);
    let same_param = reg.same_parameter(&work.0);
    if same_param {
        reg.remove_pcurves_on_surface(&work.0, &face.0);
    } else if !reg.edge_pcurves(&work.0, GeometryRegistry::shape_key(&face.0)).is_empty() {
        if let Some((c2d, _, _)) = curve_on_surface_oriented(&work, face, true) {
            if c2d.trimmed_basis().is_none() {
                reg.remove_pcurves_on_surface(&work.0, &face.0);
            }
        }
    }

    let Some(mut v1) = first_vertex(&work) else {
        return false;
    };
    let Some(mut v2) = last_vertex(&work) else {
        return false;
    };
    let edge_tol = reg.edge_tolerance(&work.0);
    let builder = TopoBuilder::new();
    let mut a = a0;
    let mut b = b0;
    let mut pieces: Vec<Edge> = Vec::new();
    let mut k = 0usize;
    while k <= seq.len() {
        let split = if k < seq.len() { seq[k] } else { b };
        let v = if k < seq.len() {
            let v = builder.make_vertex(curve.d0(split), edge_tol);
            let d1 = BRepTool::vertex_point(&v1).distance(&BRepTool::vertex_point(&v));
            if d1 < BRepTool::vertex_tolerance(&v1) * 1.01 {
                v1.set_tolerance(d1.max(BRepTool::vertex_tolerance(&v1)));
                a = split;
                v1 = v;
                k += 1;
                continue;
            }
            let d2 = BRepTool::vertex_point(&v2).distance(&BRepTool::vertex_point(&v));
            if d2 < BRepTool::vertex_tolerance(&v2) * 1.01 {
                v.set_tolerance(d2.max(BRepTool::vertex_tolerance(&v2)));
                b = split;
                v2 = v;
                k += 1;
                continue;
            }
            v
        } else {
            v2.clone()
        };
        let mut piece = copy_replace_vertices_with(&work, Some(&v1), Some(&v));
        if same_param {
            reg.set_edge_range(&piece.0, a, split);
        } else {
            let transfer = TransferParameters::with_edge_face(&work, Some(face));
            transfer.transfer_range(&mut piece, a, split, false);
        }
        pieces.push(piece);
        a = split;
        v1 = v;
        k += 1;
    }
    if pieces.is_empty() {
        return false;
    }
    if !is_fwd {
        for piece in &mut pieces {
            piece.0.reverse();
        }
        pieces.reverse();
    }
    let split_wire = builder.make_wire(&pieces);
    let mut original = edge.clone();
    if !is_fwd {
        original.0.set_orientation(Orientation::Reversed);
    }
    ctx.replace(&original.0, &split_wire.0);
    let applied = ctx.apply(&original.0);
    splice_replacement(wire, index, &applied);
    true
}

/// `UpdateWire` (`ShapeFix_Wire.cxx:4293-4309`) for the one edge just replaced.
fn splice_replacement(wire: &mut Wire, index: usize, applied: &TopoShape) {
    let new_edges: Vec<Edge> = if applied.is_wire() {
        edges_of_wire(&Wire(applied.clone()))
    } else if applied.is_edge() {
        vec![Edge(applied.clone())]
    } else {
        return;
    };
    if new_edges.is_empty() {
        return;
    }
    let parent_rev = wire.0.orientation() == Orientation::Reversed;
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let mut seen = 0usize;
    let mut pos = None;
    for (p, child) in t.children.iter().enumerate() {
        if child.shape_type() == ShapeType::Edge {
            if seen == index {
                pos = Some(p);
                break;
            }
            seen += 1;
        }
    }
    let Some(pos) = pos else {
        return;
    };
    t.children.remove(pos);
    for (k, edge) in new_edges.iter().enumerate() {
        let mut stored = edge.0.clone();
        if parent_rev {
            stored.reverse();
        }
        t.children.insert(pos + k, stored);
    }
    t.set_modified(true);
}

/// `ShapeFix_Wire.cxx:829-857`: one edge was flagged `DONE2` but not split,
/// and the wire's pcurves travel one U period. Drop that pcurve and project
/// again with `AdjustOverDegenMode = false`.
pub(in crate::shhealing) fn adjust_overdegen_pcurve(
    wire: &Wire,
    face: &Face,
    index: usize,
    preci: f64,
    cache: &mut crate::pcurve_full::ProjectorCache,
) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    if !sa_is_u_closed(surf.as_ref(), preci) {
        return;
    }
    let (uf, ul) = surf.u_range();
    let u_range = (ul - uf).abs();
    let edges = edges_of_wire(wire);
    let mut dx = 0.0;
    for edge in &edges {
        let Some((c2d, cf, cl)) = curve_on_surface_oriented(edge, face, true) else {
            return;
        };
        let p1 = c2d.d0(cf);
        let p2 = c2d.d0(cl);
        dx += p2.x() - p1.x();
    }
    if (dx.abs() - u_range).abs() >= 0.1 * u_range {
        return;
    }
    let Some(edge) = edges.get(index) else {
        return;
    };
    GeometryRegistry::global().remove_pcurves_on_surface(&edge.0, &face.0);
    cache.adjust_over_degen_off = true;
    let is_seam = is_seam_use(&edges, index);
    let _ = fix_add_pcurve(edge, face, is_seam, preci, cache);
    cache.adjust_over_degen_off = false;
}
