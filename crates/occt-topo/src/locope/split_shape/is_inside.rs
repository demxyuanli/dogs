//! The two `IsInside` predicates of `LocOpe_SplitShape`
//! (`LocOpe_SplitShape.cxx:1459-1514` and `:1516-1593`).
//!
//! Both classify 2D points of one wire against a face built from another wire:
//! the splitter asks whether a section wire lies inside a face, and whether one
//! of the face's own wires falls on the side the section leaves behind.

use occt_core::precision::{INFINITE, PCONFUSION};

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};

/// `IsInside(F, W1, W2)` (`cxx:1459-1514`).
///
/// Builds a throwaway face over `w2` alone, then classifies the 2D start point
/// of `w1`'s first edge's first vertex against it. When `w2` encloses the
/// infinite point, the wire bounds a hole and the answer flips to "outside".
///
/// `false` is both the reference's `return false` for a missing p-curve
/// (`cxx:1494-1501`) and its answer when either explorer is empty (the
/// reference dereferences `exp.Current()` there, `cxx:1487-1490`).
pub fn is_inside_wire_pair(f: &Face, w1: &Wire, w2: &Wire) -> bool {
    // `cxx:1461-1471`.
    let Some(new_face) = face_with_single_wire(f, w2) else {
        return false;
    };
    // `cxx:1473-1483`.
    let classif = BrepTopAdaptorFClass2d::new(&new_face, PCONFUSION);
    let reversed = classif.perform_infinite_point() == FaceState::In;

    // `cxx:1486-1491`: the first edge of `W1` and its first vertex.
    let Some(edg) = first_edge(&w1.0) else {
        return false;
    };
    let Some(vtx) = first_vertex(&edg) else {
        return false;
    };

    // `cxx:1492-1495`.
    let prm = BRepTool::parameter_of_vertex(&edg.0, &vtx.0);
    // `cxx:1493-1501`.
    let Some((c2d, _f, _l)) = boptools_2d::curve_on_surface_range(&edg, f) else {
        return false;
    };
    let pt2d = c2d.d0(prm);

    // `cxx:1502-1512`.
    let stat = classif.perform(pt2d, true);
    if !reversed {
        stat == FaceState::In
    } else {
        stat == FaceState::Out
    }
}

/// `IsInside(F, W)` (`cxx:1516-1593`).
///
/// Every edge of `w` is sampled at its p-curve midpoint and classified against
/// `f`. A point `OUT` rejects the wire; a point `ON` is resolved by sampling
/// ten points across the edge, and the wire is rejected only when most of them
/// fall outside.
///
/// `false` also covers a missing p-curve (`C2d` is dereferenced unconditionally
/// in the reference, `cxx:1546`), which the port answers instead of crashing.
pub fn is_inside_wire(f: &Face, w: &Wire) -> bool {
    let classif = BrepTopAdaptorFClass2d::new(f, PCONFUSION);

    // `cxx:1519-1590`.
    for sub in crate::topo_tools::explorer_shapes(&w.0, ShapeType::Edge) {
        let Some(edg) = Edge::wrap(sub) else {
            continue;
        };
        let Some((c2d, f_par, l_par)) = boptools_2d::curve_on_surface_range(&edg, f) else {
            return false;
        };

        // `cxx:1527-1545`.
        let mut prm = midpoint_parameter(f_par, l_par);

        // `cxx:1547-1557`.
        let mut pt2d = c2d.d0(prm);
        let mut stat = classif.perform(pt2d, true);
        if stat == FaceState::Out {
            return false;
        }

        // `cxx:1559-1586`.
        if stat == FaceState::On {
            let nb_pnt = 10i32;
            let mut nb_out = 0i32;
            let mut nb_in = 0i32;
            let mut nb_on = 0i32;
            for j in 1..=nb_pnt {
                prm = f_par + (l_par - f_par) / nb_pnt as f64 * (j - 1) as f64;
                pt2d = c2d.d0(prm);
                stat = classif.perform(pt2d, true);
                if stat == FaceState::Out {
                    nb_out += 1;
                } else if stat == FaceState::In {
                    nb_in += 1;
                } else {
                    nb_on += 1;
                }
            }
            if nb_out > nb_in + nb_on {
                return false;
            }
        }
    }
    true
}

/// `cxx:1527-1545`: the interior sample parameter, with the infinite-range
/// fallbacks of `Precision::IsNegativeInfinite` / `IsPositiveInfinite`.
fn midpoint_parameter(f: f64, l: f64) -> f64 {
    let neg_inf = f < -INFINITE;
    let pos_inf = l > INFINITE;
    if !neg_inf && !pos_inf {
        (f + l) / 2.0
    } else if neg_inf && pos_inf {
        0.0
    } else if neg_inf {
        l - 1.0
    } else {
        f + 1.0
    }
}

/// `F.EmptyCopied()` + `Orientation(FORWARD)` + `B.Add(newFace, W)`
/// (`cxx:1462-1470`, `cxx:1463`).
fn face_with_single_wire(f: &Face, w: &Wire) -> Option<Face> {
    let mut new_face = Face::wrap(f.0.empty_copied())?;
    new_face.0.set_orientation(Orientation::Forward);
    TopoBuilder::new().add_wire(&mut new_face, w);
    Some(new_face)
}

/// `TopExp_Explorer(W, TopAbs_EDGE).Current()`.
fn first_edge(w: &TopoShape) -> Option<Edge> {
    crate::topo_tools::explorer_shapes(w, ShapeType::Edge)
        .into_iter()
        .next()
        .and_then(Edge::wrap)
}

/// `TopExp_Explorer(E, TopAbs_VERTEX).Current()`.
fn first_vertex(e: &Edge) -> Option<Vertex> {
    crate::topo_tools::explorer_shapes_avoiding(&e.0, ShapeType::Vertex, ShapeType::Edge)
        .into_iter()
        .next()
        .and_then(Vertex::wrap)
}
