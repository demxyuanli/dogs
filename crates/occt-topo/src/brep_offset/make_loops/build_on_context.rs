//! `BRepOffset_MakeLoops::BuildOnContext` (`BRepOffset_MakeLoops.cxx:217-451`)
//! together with the file-local `IsBetweenCorks` (`:184-213`) it calls.

use std::collections::HashSet;

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_algo_loop::BRepAlgoLoop;
use crate::brep_offset::Analyse;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

use super::{substitute_vertices, BRepOffsetMakeLoops};

impl BRepOffsetMakeLoops {
    /// `BuildOnContext(LContext, Analyse, AsDes, Image, InSide, theRange)`
    /// (`cxx:217-451`): unwind the caps `LContext`, using the connectivity edges
    /// between caps (`AsDes` descendants whose every ascendant is a context
    /// face) as the only edges that may be cut.
    pub fn build_on_context(
        &mut self,
        l_context: &[TopoShape],
        analyse: &Analyse,
        as_des: &mut AsDes,
        image: &mut BRepAlgoImage,
        in_side: bool,
    ) {
        let mut loops = BRepAlgoLoop::new();
        // `cxx:229`: `Loops.VerticesForSubstitute(myVerVerMap)`.
        loops.set_vertices_for_substitute(self.my_ver_ver_map.clone());

        // `NCollection_Map<TopoDS_Shape> MapExtent` (`cxx:231`), `IsSame`
        // identity and insertion order irrelevant to the only use
        // (`MapExtent.Contains(E)`).
        let mut map_extent: HashSet<ShapeMapKey> = HashSet::new();

        for f in l_context {
            let face = Face(f.clone());

            // `NCollection_Map<TopoDS_Shape> MBound` (`cxx:241`).
            let mut m_bound: HashSet<ShapeMapKey> = HashSet::new();

            // `cxx:245-257`: `F` is reversed on the `InSide` branch because it
            // is later reversed again in the final result.
            if in_side {
                let reversed_f = Face(f.oriented(f.orientation().reversed()));
                loops.init(&reversed_f);
            } else {
                loops.init(&face);
            }

            // `cxx:260-279`: the edges of `F` that were not modified by
            // construction.
            let f_forward = face.0.oriented(Orientation::Forward);
            for ce_shape in explorer_shapes(&f_forward, ShapeType::Edge) {
                let ce = Edge(ce_shape);
                m_bound.insert(shape_map_key(&ce.0));
                if analyse.has_ancestor(&ce.0) {
                    if in_side {
                        loops.add_const_edge(&ce);
                    } else {
                        loops.add_const_edge(&Edge(ce.0.oriented(ce.0.orientation().reversed())));
                    }
                }
            }

            // `cxx:282-284`.
            let le: Vec<TopoShape> = as_des.descendant(&face.0).to_vec();
            let mut added_edges: Vec<TopoShape> = Vec::new();

            // `cxx:286-370`.
            for e_shape in &le {
                let e = Edge(e_shape.clone());
                if image.has_image(&e.0) {
                    // `cxx:287-320`: `E` was already cut on another face.
                    for ce_shape in image.image(&e.0) {
                        let mut ce = ce_shape.oriented(e.0.orientation());
                        if map_extent.contains(&shape_map_key(&e.0)) {
                            loops.add_const_edge(&Edge(ce));
                            continue;
                        }
                        if !m_bound.contains(&shape_map_key(&e.0)) {
                            ce.reverse();
                        }
                        if in_side {
                            loops.add_const_edge(&Edge(ce));
                        } else {
                            loops.add_const_edge(&Edge(
                                ce.oriented(ce.orientation().reversed()),
                            ));
                        }
                    }
                } else if self.is_between_corks(&e.0, as_des, l_context)
                    && as_des.has_descendant(&e.0)
                {
                    // `cxx:322-341`: the connection between two caps.
                    map_extent.insert(shape_map_key(&e.0));
                    if in_side {
                        let lv: Vec<TopoShape> = as_des
                            .descendant(&e.0)
                            .iter()
                            .map(|s| s.oriented(s.orientation().reversed()))
                            .collect();
                        loops.add_edge(&e, &lv);
                    } else {
                        let lv: Vec<TopoShape> = as_des.descendant(&e.0).to_vec();
                        loops.add_edge(&e, &lv);
                    }
                    added_edges.push(e.0.clone());
                } else if self.is_between_corks(&e.0, as_des, l_context) {
                    // `cxx:342-356`.
                    if in_side {
                        loops.add_const_edge(&e);
                    } else {
                        loops.add_const_edge(&Edge(e.0.oriented(e.0.orientation().reversed())));
                    }
                } else {
                    // `cxx:357-369`.
                    if in_side {
                        loops.add_const_edge(&Edge(e.0.oriented(e.0.orientation().reversed())));
                    } else {
                        loops.add_const_edge(&e);
                    }
                }
            }

            // `cxx:372-373`: unwind.
            loops.perform();
            loops.wires_to_faces();

            // `cxx:380-404`: `F => NewFaces`, then `E => NewEdges`.
            let nf: Vec<TopoShape> = loops.new_faces().to_vec();
            let _ = image.bind_list(&face.0, &nf);
            for e_shape in &added_edges {
                let e = Edge(e_shape.clone());
                let Some(new_edges) = loops.new_edges(&e).cloned() else {
                    continue;
                };
                if image.has_image(&e.0) {
                    let _ = image.add_list(&e.0, &new_edges);
                } else {
                    let _ = image.bind_list(&e.0, &new_edges);
                }
            }
        }

        // `cxx:405-408`.
        self.my_ver_ver_map = loops.vertices_for_substitute().clone();
        if self.my_ver_ver_map.is_empty() {
            return;
        }

        // `cxx:409-450`.
        let map = self.my_ver_ver_map.clone();
        substitute_vertices(l_context, image, as_des, &map);
    }
}
