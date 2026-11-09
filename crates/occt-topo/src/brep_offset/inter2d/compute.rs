//! `BRepOffset_Inter2d::Compute` (`BRepOffset_Inter2d.cxx:1725-1827`): the
//! pairwise 2D intersection of the edges the 3D phase attached to one face.
//!
//! The face's edge list is `AsDes->Descendant(F)` (`cxx:1754`), i.e. the edges
//! `BRepOffset_Tool::Inter3D` recorded, and every pair `(E1, E2)` with
//! `E2` earlier in the list is intersected by `EdgeInter`.
//!
//! A pair is skipped (`cxx:1772-1806`) when
//! * `theEdgeIntEdges` already ties one to the other (directly, or through a
//!   third edge bound to `E1`'s list), or
//! * both edges belong to the face itself *and* neither is a new edge - the
//!   face's own edges are left untouched (`cxx:1736-1744`).
//!
//! `Message_ProgressScope` (`cxx:1760-1764`) has no counterpart; the loop runs
//! to completion.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo::AsDes;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::edge_inter::edge_inter;
use super::{BRepOffsetInter2d, Dmvv};

impl BRepOffsetInter2d {
/// `BRepOffset_Inter2d::Compute(AsDes, F, NewEdges, Tol, theEdgeIntEdges,
/// theDMVV, theRange)` (`cxx:1725-1827`).
///
/// `new_edges` is the `NCollection_IndexedMap<TopoDS_Shape>` of `NewEdges`
/// (`IsSame` membership, hence [`ShapeMapKey`]) and `edge_int_edges` the
/// `NCollection_DataMap` of `theEdgeIntEdges`.
pub fn compute(
    as_des: &mut AsDes,
    f: &Face,
    new_edges: &HashSet<ShapeMapKey>,
    tol: f64,
    edge_int_edges: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    dmvv: &mut Dmvv,
) {
    // `cxx:1737-1744`: the face's own edges.
    let mut edges_of_face: HashSet<ShapeMapKey> = HashSet::new();
    for e in explorer_shapes(&f.0, ShapeType::Edge) {
        edges_of_face.insert(shape_map_key(&e));
    }

    // `cxx:1754`: the descendant list is read into `LE` before the loop;
    // `EdgeInter` only appends descendants of `E1`/`E2`, never of `F`, so the
    // copy is what the reference iterates as well.
    let le: Vec<TopoShape> = as_des.descendant(&f.0).to_vec();

    // `cxx:1758`: `BRepAdaptor_Surface BAsurf(F)`, built once for the face.
    let basurf = SurfaceAdaptor::of_face(f);
    // `cxx:1817`: `F.Oriented(TopAbs_FORWARD)` reaches `EdgeInter`.
    let f_forward = Face(f.0.oriented(Orientation::Forward));

    // `cxx:1765-1824`: `i` runs over the list, `j < i` over the earlier entries.
    for (i, e1_shape) in le.iter().enumerate() {
        let e1 = Edge(e1_shape.clone());
        for e2_shape in le.iter().take(i) {
            let e2 = Edge(e2_shape.clone());

            // `cxx:1774-1806`.
            let mut to_intersect = true;
            if let Some(e_list) = edge_int_edges.get(&shape_map_key(e1_shape)) {
                for it_edge in e_list {
                    if is_same(it_edge, e2_shape) {
                        to_intersect = false;
                    }
                }
                if to_intersect {
                    for edge in e_list {
                        if let Some(e_list2) = edge_int_edges.get(&shape_map_key(edge)) {
                            for it_edge2 in e_list2 {
                                if is_same(it_edge2, e2_shape) {
                                    to_intersect = false;
                                }
                            }
                        }
                    }
                }
            }

            // `cxx:1811-1820`.
            if to_intersect
                && (!edges_of_face.contains(&shape_map_key(e1_shape))
                    || !edges_of_face.contains(&shape_map_key(e2_shape)))
                && (new_edges.contains(&shape_map_key(e1_shape))
                    || new_edges.contains(&shape_map_key(e2_shape)))
            {
                edge_inter(&f_forward, &basurf, &e1, &e2, as_des, tol, true, dmvv);
            }
        }
    }
}
}
