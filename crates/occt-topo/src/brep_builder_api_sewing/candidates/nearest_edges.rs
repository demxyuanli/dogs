//! `AnalysisNearestEdges` (cxx:1618-1782) and `MergedNearestEdges` (cxx:4242-4431).

use occt_core::precision::{CONFUSION, INFINITE};
use occt_geom::gcpnts::curve_length_range;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};

use super::super::shape_maps::{IntIndexedMap, IntSequence, ShapeSet};
use super::super::Sewing;


impl Sewing {
    /// `AnalysisNearestEdges(sequenceSec, seqIndCandidate, seqOrientations,
    /// evalDist)` (`cxx:1618-1782`): drop the candidates that belong to the
    /// same face as the worked edge, and (when `evalDist`) the ones that have a
    /// nearer candidate on an inner bound.
    pub(crate) fn analysis_nearest_edges(
        &self,
        sequence_sec: &[TopoShape],
        seq_ind_candidate: &mut IntSequence,
        seq_orientations: &mut IntSequence,
        eval_dist: bool,
    ) {
        let Some(work_index) = seq_ind_candidate.first() else {
            return;
        };
        let workedge = match sequence_sec.get((work_index - 1) as usize) {
            Some(s) => s.clone(),
            None => return,
        };
        let mut bnd = workedge.clone();
        if let Some(b) = self.section_bound.get(&workedge) {
            bnd = b;
        }
        let Some(workfaces) = self.bound_faces.get(&bnd) else {
            return;
        };
        if workfaces.is_empty() {
            return;
        }
        let mut map_faces = ShapeSet::new();
        for f in &workfaces {
            map_faces.add(f);
        }
        let mut seq_not_candidate = IntSequence::new();
        let _seq_new_forward = IntSequence::new();
        seq_not_candidate.append(work_index);
        let mut i = 1usize;
        while i <= seq_ind_candidate.len() {
            let index = seq_ind_candidate.value(i).unwrap();
            let mut is_remove = false;
            if index == work_index {
                seq_ind_candidate.remove(i);
                seq_orientations.remove(i);
                is_remove = true;
            }
            if !is_remove {
                let mut bnd2 = match sequence_sec.get((index - 1) as usize) {
                    Some(s) => s.clone(),
                    None => {
                        i += 1;
                        continue;
                    }
                };
                if let Some(b) = self.section_bound.get(&bnd2) {
                    bnd2 = b;
                }
                if let Some(listfaces) = self.bound_faces.get(&bnd2) {
                    let mut is_merged = true;
                    for f in &listfaces {
                        if !is_merged {
                            break;
                        }
                        if map_faces.contains(f) {
                            let face = Face(f.clone());
                            let is_merged_closed = match BRepTool::face_surface(&face) {
                                Some(surf) => {
                                    (self.is_u_closed_surface(&surf, &bnd2)
                                        || self.is_v_closed_surface(&surf, &bnd2))
                                        && self.is_merged_closed(
                                            &Edge(sequence_sec[(index - 1) as usize].clone()),
                                            &Edge(workedge.clone()),
                                            &face,
                                        )
                                }
                                None => false,
                            };
                            is_merged = is_merged_closed;
                        }
                    }
                    if !is_merged {
                        seq_not_candidate.append(index);
                        seq_ind_candidate.remove(i);
                        seq_orientations.remove(i);
                        is_remove = true;
                    }
                } else {
                    seq_ind_candidate.remove(i);
                    seq_orientations.remove(i);
                    is_remove = true;
                }
            }
            if !is_remove {
                i += 1;
            }
        }
        if seq_ind_candidate.len() == 0 || seq_not_candidate.len() == 1 {
            return;
        }
        if !eval_dist {
            return;
        }
        let nb_not = seq_not_candidate.len();
        let nb_cand = seq_ind_candidate.len();
        let mut tot_tab_dist = vec![vec![0.0f64; nb_cand + 1]; nb_not + 1];
        let mut map_index = IntIndexedMap::new();
        let mut seq_forward = IntSequence::new();
        for k in 1..=nb_not {
            let index1 = seq_not_candidate.value(k).unwrap();
            let Some(edge) = sequence_sec.get((index1 - 1) as usize) else {
                continue;
            };
            let mut tmp_seq: Vec<TopoShape> = vec![edge.clone()];
            for kk in 1..=nb_cand {
                let ic = seq_ind_candidate.value(kk).unwrap();
                tmp_seq.push(sequence_sec[(ic - 1) as usize].clone());
            }
            let leng_sec = tmp_seq.len();
            let mut tab_forward = vec![true; leng_sec];
            let mut tab_dist = vec![-1.0f64; leng_sec];
            let mut arr_len = vec![0.0f64; leng_sec];
            let mut tab_min_dist = vec![INFINITE; leng_sec];
            self.evaluate_distances(
                &tmp_seq,
                &mut tab_forward,
                &mut tab_dist,
                &mut arr_len,
                &mut tab_min_dist,
                1,
            );
            if k == 1 {
                for n in 1..leng_sec {
                    if tab_dist[n] == -1.0 || tab_dist[n] > self.tolerance {
                        map_index.add(n as i32);
                        continue;
                    }
                    tot_tab_dist[k][n] = tab_dist[n];
                    seq_forward.append(if tab_forward[n] { 1 } else { 0 });
                }
            } else {
                for n in 1..leng_sec {
                    if tab_dist[n - 1] == -1.0 || tab_dist[n - 1] > self.tolerance {
                        continue;
                    }
                    if tab_dist[n] < tot_tab_dist[1][n] {
                        map_index.add(n as i32);
                    }
                }
            }
        }
        let mut i2 = nb_cand as i32;
        while i2 >= 1 {
            if map_index.contains(i2) {
                seq_ind_candidate.remove(i2 as usize);
                seq_orientations.remove(i2 as usize);
            }
            i2 -= 1;
        }
    }
    /// `MergedNearestEdges(edge, SeqMergedEdge, SeqMergedOri)`
    /// (`cxx:4242-4431`): every bound connected to `edge` through its nodes,
    /// then the candidate search over that sequence.
    pub(crate) fn merged_nearest_edges(
        &self,
        edge: &TopoShape,
        seq_merged_edge: &mut Vec<TopoShape>,
        seq_merged_ori: &mut Vec<bool>,
    ) -> bool {
        // Retrieve the edge nodes (`cxx:3811-3818`).
        let (no1, no2) = match crate::topo_tools::edge_vertices(&Edge(edge.clone())) {
            (Some(a), Some(b)) => (a.0, b.0),
            _ => return false,
        };
        let is_node1 = self.vertex_node.contains(&no1);
        let is_node2 = self.vertex_node.contains(&no2);
        let nno1 = if is_node1 {
            self.vertex_node.get(&no1).unwrap_or(no1.clone())
        } else {
            no1.clone()
        };
        let nno2 = if is_node2 {
            self.vertex_node.get(&no2).unwrap_or(no2.clone())
        } else {
            no2.clone()
        };
        // Nodes connected to node #1 (`cxx:3821-3841`).
        let mut map_vert1 = ShapeSet::new();
        map_vert1.add(&nno1);
        if let Some(list) = self.cutting_node.get(&nno1) {
            for v1 in &list {
                map_vert1.add(v1);
                if !is_node1 {
                    if let Some(inner) = self.cutting_node.get(v1) {
                        for vn in &inner {
                            map_vert1.add(vn);
                        }
                    }
                }
            }
        }
        // Nodes connected to node #2 (`cxx:3843-3863`).
        let mut map_vert2 = ShapeSet::new();
        map_vert2.add(&nno2);
        if let Some(list) = self.cutting_node.get(&nno2) {
            for v1 in &list {
                map_vert2.add(v1);
                if !is_node2 {
                    if let Some(inner) = self.cutting_node.get(v1) {
                        for vn in &inner {
                            map_vert2.add(vn);
                        }
                    }
                }
            }
        }
        // Find all possible contiguous edges (`cxx:3865-3918`).
        let mut seq_edges: Vec<TopoShape> = vec![edge.clone()];
        let mut map_edges = ShapeSet::new();
        map_edges.add(edge);
        for i in 1..=map_vert1.extent() {
            let node1 = map_vert1.key_at(i).unwrap().clone();
            let Some(sections) = self.node_sections.get(&node1) else {
                continue;
            };
            for sec in &sections {
                if crate::topo_tools::is_same(sec, edge) {
                    continue;
                }
                let (vs1, vs2) = crate::topo_tools::edge_vertices(&Edge(sec.clone()));
                let (Some(vs1), Some(vs2)) = (vs1, vs2) else {
                    continue;
                };
                let vs1 = vs1.0;
                let vs2 = vs2.0;
                let vs1n = self.vertex_node.get(&vs1).unwrap_or(vs1.clone());
                let vs2n = self.vertex_node.get(&vs2).unwrap_or(vs2.clone());
                let connected = (map_vert1.contains(&vs1n) && map_vert2.contains(&vs2n))
                    || (map_vert1.contains(&vs2n) && map_vert2.contains(&vs1n));
                if !connected {
                    continue;
                }
                if !map_edges.add(sec) {
                    continue;
                }
                // Check for a rejected cutting (`cxx:3892-3914`).
                let mut is_rejected = self.merged_edges.contains(sec);
                if !is_rejected {
                    if let Some(its) = self.bound_sections.get(sec) {
                        for section in &its {
                            if self.merged_edges.contains(section) {
                                is_rejected = true;
                                break;
                            }
                        }
                    }
                }
                if !is_rejected {
                    if let Some(bnd) = self.section_bound.get(sec) {
                        is_rejected =
                            !self.bound_sections.contains(&bnd) || self.merged_edges.contains(&bnd);
                    }
                }
                if !is_rejected {
                    seq_edges.push(sec.clone());
                }
            }
        }
        map_edges.clear();

        let mut success = false;
        let nb_section = seq_edges.len();
        if nb_section > 1 {
            // Find the longest edge (`cxx:3925-3948`).
            let mut ind_ref = 1usize;
            if self.nonmanifold {
                let mut len_ref = 0.0;
                for i in 1..=nb_section {
                    let (f, l) = BRepTool::edge_parameters(&Edge(seq_edges[i - 1].clone()));
                    let len = match BRepTool::edge_curve_adaptor_world(&Edge(seq_edges[i - 1].clone()))
                    {
                        Some(c) => curve_length_range(c.as_ref(), f, l, CONFUSION),
                        None => 0.0,
                    };
                    if len > len_ref {
                        ind_ref = i;
                        len_ref = len;
                    }
                }
                if ind_ref != 1 {
                    seq_edges.swap(ind_ref - 1, 0);
                }
            }
            // Find the merging candidates (`cxx:3950-3971`).
            let mut seq_forward = IntSequence::new();
            let mut seq_candidates = IntSequence::new();
            let mut map_reference = IntIndexedMap::new();
            map_reference.add(ind_ref as i32);
            if self.find_candidates(
                &seq_edges,
                &mut map_reference,
                &mut seq_candidates,
                &mut seq_forward,
            ) {
                let nb_candidates = seq_candidates.len();
                for i in 1..=nb_candidates {
                    let ind = seq_candidates.value(i).unwrap();
                    let iedge = seq_edges[(ind - 1) as usize].clone();
                    let ori = seq_forward.value(i).unwrap_or(1) != 0;
                    seq_merged_edge.push(iedge);
                    seq_merged_ori.push(ori);
                    if !self.nonmanifold {
                        break;
                    }
                }
                success = nb_candidates != 0;
            }
        }
        success
    }
}
