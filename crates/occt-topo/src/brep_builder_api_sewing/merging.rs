//! `Merging` (`cxx:3788-4239`): merge the contiguous bounds (and their cutting
//! sections) through `MergedNearestEdges` + `SameParameterEdge`.

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};
use crate::shape_fix_compose_shell::MapReShape;

use super::shape_maps::{ShapeKeyedMap, ShapeSet};
use super::vertices::replace_edge;
use super::Sewing;

impl Sewing {
    /// `Merging(firstTime)` (`cxx:3788-4239`).
    pub(crate) fn merging(&mut self, _first_time: bool) {
        let builder = TopoBuilder::new();
        let bounds = self.bound_faces.keys();
        for bound in bounds {
            // If the bound was already merged - continue (`cxx:3414-3415`).
            if self.merged_edges.contains(&bound) {
                continue;
            }
            let list_faces = self.bound_faces.get(&bound).unwrap_or_default();
            if list_faces.is_empty() {
                // Merge a free edge - only the vertices (`cxx:3419-3442`).
                let (no1, no2) = match crate::topo_tools::edge_vertices(&Edge(bound.clone())) {
                    (Some(a), Some(b)) => (a.0, b.0),
                    _ => {
                        self.merged_edges.add(&bound);
                        continue;
                    }
                };
                let mut nno1 = self.vertex_node_free.get(&no1).unwrap_or_else(|| no1.clone());
                let mut nno2 = self.vertex_node_free.get(&no2).unwrap_or_else(|| no2.clone());
                if !crate::topo_tools::is_same(&no1, &nno1) {
                    nno1.set_orientation(no1.orientation());
                    self.reshape.replace(&no1, &nno1);
                }
                if !crate::topo_tools::is_same(&no2, &nno2) {
                    nno2.set_orientation(no2.orientation());
                    self.reshape.replace(&no2, &nno2);
                }
                self.merged_edges.add(&bound);
                continue;
            }

            // Check for a previous splitting, build the replacing wire
            // (`cxx:3444-3461`).
            let mut bound_wire: Option<TopoShape> = None;
            let mut is_prev_split = false;
            let has_cutting_sections = self.bound_sections.contains(&bound);
            if has_cutting_sections {
                let mut wire = TopoShape::new(ShapeType::Wire);
                wire.set_orientation(bound.orientation());
                if let Some(sections) = self.bound_sections.get(&bound) {
                    for section in &sections {
                        builder.add(&mut wire, section);
                        if self.merged_edges.contains(section) {
                            is_prev_split = true;
                        }
                    }
                }
                bound_wire = Some(wire);
            }

            // Merge with the bound (`cxx:3465-3565`).
            let mut merged_with_bound: ShapeKeyedMap<TopoShape> = ShapeKeyedMap::new();
            if !is_prev_split {
                let mut seq: Vec<TopoShape> = Vec::new();
                let mut seq_ori: Vec<bool> = Vec::new();
                if self.merged_nearest_edges(&bound, &mut seq, &mut seq_ori) {
                    merged_with_bound.add(&bound, bound.clone());
                    let mut ii = 1usize;
                    while ii <= seq.len() {
                        let iedge = seq[ii - 1].clone();
                        let mut is_rejected = self.merged_edges.contains(&iedge)
                            || merged_with_bound.contains(&iedge);
                        if !is_rejected {
                            if let Some(sections) = self.bound_sections.get(&iedge) {
                                for sec in &sections {
                                    if self.merged_edges.contains(sec)
                                        || merged_with_bound.contains(sec)
                                    {
                                        is_rejected = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if !is_rejected {
                            if let Some(bnd) = self.section_bound.get(&iedge) {
                                is_rejected = self.merged_edges.contains(&bnd)
                                    || merged_with_bound.contains(&bnd);
                            }
                        }
                        if is_rejected {
                            seq.remove(ii - 1);
                            seq_ori.remove(ii - 1);
                        } else {
                            merged_with_bound.add(&iedge, iedge.clone());
                            ii += 1;
                        }
                    }
                    let mut nb_merged = seq.len();
                    if nb_merged > 0 {
                        // Create the same-parameter edge (`cxx:3522-3529`).
                        let mut actually_merged = ShapeSet::new();
                        let mut reshape = std::mem::take(&mut self.reshape);
                        let merged_edge = self.same_parameter_edge_seq(
                            &bound,
                            &seq,
                            &seq_ori,
                            &mut ShapeSet::new(),
                            &mut actually_merged,
                            &mut reshape,
                            true,
                        );
                        self.reshape = reshape;
                        let is_forward = merged_edge
                            .as_ref()
                            .map(|e| e.0.orientation() == Orientation::Forward)
                            .unwrap_or(false);
                        // Process the actually merged edges (`cxx:3533-3558`).
                        let mut nb_actually_merged = 0usize;
                        for ii in 1..=nb_merged {
                            let iedge = seq[ii - 1].clone();
                            if actually_merged.contains(&iedge) {
                                nb_actually_merged += 1;
                                let mut orient = iedge.orientation();
                                if !is_forward {
                                    orient = orient.reversed();
                                }
                                if !seq_ori[ii - 1] {
                                    orient = orient.reversed();
                                }
                                if let Some(me) = &merged_edge {
                                    merged_with_bound.bind(&iedge, me.0.oriented(orient));
                                }
                            } else {
                                merged_with_bound.remove(&iedge);
                            }
                        }
                        if nb_actually_merged > 0 {
                            let mut orient = bound.orientation();
                            if !is_forward {
                                orient = orient.reversed();
                            }
                            if let Some(me) = &merged_edge {
                                merged_with_bound.bind(&bound, me.0.oriented(orient));
                            }
                        }
                        nb_merged = nb_actually_merged;
                    }
                    if nb_merged == 0 {
                        merged_with_bound.remove(&bound);
                    }
                }
            }
            let is_merged = !merged_with_bound.is_empty();

            // Merge with the cutting sections (`cxx:3568-3697`).
            let mut sections_reshape = MapReShape::new();
            let mut merged_with_sections: ShapeKeyedMap<TopoShape> = ShapeKeyedMap::new();
            if has_cutting_sections {
                let sections = self.bound_sections.get(&bound).unwrap_or_default();
                for section_in in sections {
                    let section = section_in;
                    if self.merged_edges.contains(&section) {
                        continue;
                    }
                    let mut seq: Vec<TopoShape> = Vec::new();
                    let mut seq_ori: Vec<bool> = Vec::new();
                    if self.merged_nearest_edges(&section, &mut seq, &mut seq_ori) {
                        merged_with_sections.add(&section, section.clone());
                        let mut ii = 1usize;
                        while ii <= seq.len() {
                            let iedge = seq[ii - 1].clone();
                            let mut is_rejected = self.merged_edges.contains(&iedge)
                                || merged_with_sections.contains(&iedge);
                            if !is_rejected {
                                if let Some(secs) = self.bound_sections.get(&iedge) {
                                    for sec in &secs {
                                        if self.merged_edges.contains(sec)
                                            || merged_with_sections.contains(sec)
                                        {
                                            is_rejected = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            if !is_rejected {
                                if let Some(bnd) = self.section_bound.get(&iedge) {
                                    is_rejected = self.merged_edges.contains(&bnd)
                                        || merged_with_sections.contains(&bnd);
                                }
                            }
                            if is_rejected {
                                seq.remove(ii - 1);
                                seq_ori.remove(ii - 1);
                            } else {
                                merged_with_sections.add(&iedge, iedge.clone());
                                ii += 1;
                            }
                        }
                        let mut nb_merged = seq.len();
                        if nb_merged > 0 {
                            let mut actually_merged = ShapeSet::new();
                            let reshape = std::mem::take(&mut self.reshape);
                            let merged_edge = self.same_parameter_edge_seq(
                                &section,
                                &seq,
                                &seq_ori,
                                &mut ShapeSet::new(),
                                &mut actually_merged,
                                &mut sections_reshape,
                                false,
                            );
                            self.reshape = reshape;
                            let is_forward = merged_edge
                                .as_ref()
                                .map(|e| e.0.orientation() == Orientation::Forward)
                                .unwrap_or(false);
                            let mut nb_actually_merged = 0usize;
                            for ii in 1..=nb_merged {
                                let iedge = seq[ii - 1].clone();
                                if actually_merged.contains(&iedge) {
                                    nb_actually_merged += 1;
                                    let mut orient = iedge.orientation();
                                    if !is_forward {
                                        orient = orient.reversed();
                                    }
                                    if !seq_ori[ii - 1] {
                                        orient = orient.reversed();
                                    }
                                    if let Some(me) = &merged_edge {
                                        let oedge = me.0.oriented(orient);
                                        merged_with_sections.bind(&iedge, oedge.clone());
                                        let applied = self
                                            .reshape
                                            .apply_option(&iedge)
                                            .unwrap_or_else(|| iedge.clone());
                                        replace_edge(&applied, &oedge, &mut sections_reshape);
                                    }
                                } else {
                                    merged_with_sections.remove(&iedge);
                                }
                            }
                            if nb_actually_merged > 0 {
                                let mut orient = section.orientation();
                                if !is_forward {
                                    orient = orient.reversed();
                                }
                                if let Some(me) = &merged_edge {
                                    let oedge = me.0.oriented(orient);
                                    merged_with_sections.bind(&section, oedge.clone());
                                    let applied = self
                                        .reshape
                                        .apply_option(&section)
                                        .unwrap_or_else(|| section.clone());
                                    replace_edge(&applied, &oedge, &mut sections_reshape);
                                }
                            }
                            nb_merged = nb_actually_merged;
                        }
                        if nb_merged == 0 {
                            merged_with_sections.remove(&section);
                        }
                    } else if is_merged {
                        // Reject the merging of the sections (`cxx:3686-3691`).
                        merged_with_sections.clear();
                        break;
                    }
                }
            }
            let is_merged_split = !merged_with_sections.is_empty();

            if !is_merged && !is_merged_split {
                // Nothing was merged in this iteration (`cxx:3694-3707`).
                if is_prev_split {
                    let a = self.reshape.apply_option(&bound).unwrap_or_else(|| bound.clone());
                    let b = bound_wire
                        .as_ref()
                        .and_then(|w| self.reshape.apply_option(w))
                        .unwrap_or_else(|| bound.clone());
                    self.reshape.replace(&a, &b);
                }
                continue;
            }

            // Set the splitting flag (`cxx:3710`).
            let mut is_splitted = (!is_merged && is_merged_split) || is_prev_split;

            // Choose between bound and sections merging (`cxx:3713-3759`).
            if is_merged && is_merged_split && !is_prev_split {
                let mut map_split_edges = ShapeSet::new();
                for (edge, _) in merged_with_sections.iter() {
                    map_split_edges.add(edge);
                }
                let bound_edges: Vec<TopoShape> = merged_with_bound.keys();
                for edge in bound_edges {
                    map_split_edges.remove(&edge);
                    if let Some(sections) = self.bound_sections.get(&edge) {
                        for sec in &sections {
                            map_split_edges.remove(sec);
                        }
                    }
                }
                // Calculate the section merging tolerance (`cxx:3745-3756`).
                let mut min_split_tol = f64::MAX;
                for ii in 1..=map_split_edges.extent() {
                    let key = map_split_edges.key_at(ii).unwrap();
                    if let Some(edge) = merged_with_sections.get(key) {
                        min_split_tol =
                            min_split_tol.min(BRepTool::edge_tolerance(&Edge(edge.clone())));
                    }
                }
                let bound_edge_tol = merged_with_bound
                    .get(&bound)
                    .map(|e| BRepTool::edge_tolerance(&Edge(e.clone())))
                    .unwrap_or(0.0);
                is_splitted = (min_split_tol < bound_edge_tol + self.min_tolerance)
                    || self.nonmanifold;
                is_splitted = !map_split_edges.is_empty() && is_splitted;
            }

            if is_splitted {
                // Merging of the cutting sections (`cxx:3763-3778`).
                let a = self.reshape.apply_option(&bound).unwrap_or_else(|| bound.clone());
                let b = bound_wire
                    .as_ref()
                    .and_then(|w| self.reshape.apply_option(w))
                    .unwrap_or_else(|| bound.clone());
                self.reshape.replace(&a, &b);
                let entries: Vec<(TopoShape, TopoShape)> = merged_with_sections
                    .iter()
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect();
                for (oldedge, newedge_value) in entries {
                    let newedge = sections_reshape
                        .apply_option(&newedge_value)
                        .unwrap_or_else(|| newedge_value.clone());
                    let applied = self
                        .reshape
                        .apply_option(&oldedge)
                        .unwrap_or_else(|| oldedge.clone());
                    let mut reshape = std::mem::take(&mut self.reshape);
                    replace_edge(&applied, &newedge, &mut reshape);
                    self.reshape = reshape;
                    self.merged_edges.add(&oldedge);
                    self.bound_sections.remove(&oldedge);
                }
            } else {
                // Merging of the initial bound (`cxx:3781-3796`).
                let entries: Vec<(TopoShape, TopoShape)> = merged_with_bound
                    .iter()
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect();
                for (oldedge, newedge) in entries {
                    let applied = self
                        .reshape
                        .apply_option(&oldedge)
                        .unwrap_or_else(|| oldedge.clone());
                    let mut reshape = std::mem::take(&mut self.reshape);
                    replace_edge(&applied, &newedge, &mut reshape);
                    self.reshape = reshape;
                    self.merged_edges.add(&oldedge);
                    self.bound_sections.remove(&oldedge);
                }
                self.bound_sections.remove(&bound);
                self.merged_edges.add(&bound);
            }
        }

        // `cxx:3798-3803`.
        self.nb_vertices = self.vertex_node.extent() + self.vertex_node_free.extent();
        self.node_sections.clear();
        self.vertex_node.clear();
        self.vertex_node_free.clear();
        self.cutting_node.clear();
        let _ = CONFUSION;
    }
}



