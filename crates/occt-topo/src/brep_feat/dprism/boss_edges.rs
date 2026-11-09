//! `BRepFeat_MakeDPrism::BossEdges` / `TopEdges` / `LatEdges`
//! (`cxx:1015-1148`).
//!
//! `BossEdges` splits the edges of the result into the ones on the top of the
//! feature (`myTopEdges`) and the ones on its sides (`myLatEdges`): the top is
//! the edge set of the first or the last shape of the draft prism, and the
//! sides are the edges of the basis faces that touch a top edge without being
//! one of the top faces.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::locope::shape_maps::ShapeSet;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, edges_of, faces_of, is_same, explorer_shapes};

use super::MakeDPrism;

/// `V1.IsSame(VT1)` with a possibly null vertex on either side (`TopExp::
/// Vertices` leaves `V2` null on an open edge with one vertex).
fn same_vertex(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => is_same(&a.0, &b.0),
        _ => false,
    }
}

impl MakeDPrism {
    /// `BossEdges(signature)` (`cxx:1015-1125`).
    ///
    /// `signature` picks the top of the feature: `1` / `-1` the first shape of
    /// the draft prism, `2` / `-2` the last one. Any other value returns
    /// without touching the lists.
    pub fn boss_edges(&mut self, signature: i32) {
        // `cxx:1024-1037`.
        let the_last_shape: Vec<TopoShape> = if signature == 1 || signature == -1 {
            self.form.first_shape().to_vec()
        } else if signature == 2 || signature == -2 {
            self.form.last_shape().to_vec()
        } else {
            return;
        };

        // `cxx:1039-1049`: the edges of the top.
        for ls in &the_last_shape {
            for e in explorer_shapes(ls, ShapeType::Edge) {
                self.my_top_edges.push(e);
            }
        }

        // `cxx:1051-1054`.
        if signature < 0 {
            self.my_lat_edges = self.form.new_edges().to_vec();
        } else if signature > 0 {
            // `cxx:1056-1125`.
            let Some(shape) = self.form.my_shape.clone() else {
                return;
            };
            let mut map_e = ShapeSet::new();

            // `cxx:1063-1105`: once a face of the result touches a top edge
            // through one of its vertices, its edges are toggled in `MapE`.
            for ff in faces_of(&shape) {
                let mut found = false;
                for ls in &the_last_shape {
                    if is_same(&ff.0, ls) {
                        continue;
                    }
                    for e1 in edges_of(&ff.0) {
                        if found {
                            break;
                        }
                        // `cxx:1076-1077`.
                        let (v1, v2) = edge_vertices(&e1);
                        for e2 in &self.my_top_edges {
                            if found {
                                break;
                            }
                            // `cxx:1080-1085`.
                            let e2 = Edge(e2.clone());
                            let (vt1, vt2) = edge_vertices(&e2);
                            if same_vertex(&v1, &vt1)
                                || same_vertex(&v1, &vt2)
                                || same_vertex(&v2, &vt1)
                                || same_vertex(&v2, &vt2)
                            {
                                found = true;
                                // `cxx:1088-1101`.
                                for e3 in edges_of(&ff.0) {
                                    if map_e.contains(&e3.0) {
                                        map_e.remove(&e3.0);
                                    } else {
                                        map_e.add(&e3.0);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // `cxx:1108-1112`: the top edges themselves never belong to the
            // sides.
            for e in &self.my_top_edges {
                if map_e.contains(e) {
                    map_e.remove(e);
                }
            }

            // `cxx:1114-1124`.
            for k in map_e.iter() {
                if !BRepTool::is_degenerated(&Edge(k.clone())) {
                    self.my_lat_edges.push(k.clone());
                }
            }
        }
    }

    /// `TopEdges()` (`cxx:1133-1136`).
    pub fn top_edges(&self) -> &[TopoShape] {
        &self.my_top_edges
    }

    /// `LatEdges()` (`cxx:1142-1145`).
    pub fn lat_edges(&self) -> &[TopoShape] {
        &self.my_lat_edges
    }
}
