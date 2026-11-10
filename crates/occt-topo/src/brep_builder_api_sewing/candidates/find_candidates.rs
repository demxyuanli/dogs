//! `FindCandidates` (BRepBuilderAPI_Sewing.cxx:1785-2086).

use occt_core::precision::{CONFUSION, INFINITE, REAL_SMALL};

use crate::shape::{Edge, Face, TopoShape};

use super::super::shape_maps::{IntIndexedMap, IntSequence, ShapeSet};
use super::super::Sewing;


impl Sewing {
    /// `FindCandidates(seqSections, mapReference, seqCandidates,
    /// seqOrientations)` (`cxx:1785-2086`).
    pub(crate) fn find_candidates(
        &self,
        seq_sections: &[TopoShape],
        map_reference: &mut IntIndexedMap,
        seq_candidates: &mut IntSequence,
        seq_orientations: &mut IntSequence,
    ) -> bool {
        let nb_sections = seq_sections.len();
        if nb_sections <= 1 {
            return false;
        }
        let Some(ind_reference) = map_reference.last() else {
            return false;
        };
        let mut nb_candidates = 0usize;
        let mut faces1 = ShapeSet::new();
        let edge1 = Edge(seq_sections[(ind_reference - 1) as usize].clone());
        {
            let mut bnd = edge1.0.clone();
            if let Some(b) = self.section_bound.get(&bnd) {
                bnd = b;
            }
            if let Some(faces) = self.bound_faces.get(&bnd) {
                for f in &faces {
                    faces1.add(f);
                }
            }
        }
        // Check the merging conditions for the candidates (`cxx:1626-1648`).
        let mut seq_sections_new: Vec<TopoShape> = Vec::with_capacity(nb_sections);
        let mut seq_candidates_new = IntSequence::new();
        for i in 1..=nb_sections {
            if i as i32 == ind_reference {
                seq_sections_new.insert(0, edge1.0.clone());
                seq_candidates_new.prepend(i as i32);
            } else {
                seq_sections_new.push(seq_sections[i - 1].clone());
                seq_candidates_new.append(i as i32);
            }
        }
        let nb_sections_new = seq_sections_new.len();
        if nb_sections_new > 1 {
            // Evaluate the distances (`cxx:1649-1656`).
            let mut arr_forward = vec![true; nb_sections_new];
            let mut arr_distance = vec![0.0f64; nb_sections_new];
            let mut arr_len = vec![0.0f64; nb_sections_new];
            let mut arr_min_dist = vec![INFINITE; nb_sections_new];
            self.evaluate_distances(
                &seq_sections_new,
                &mut arr_forward,
                &mut arr_distance,
                &mut arr_len,
                &mut arr_min_dist,
                1,
            );
            // Fill the sequence of candidate indices sorted by distance
            // (`cxx:1657-1687`).
            for i in 2..=nb_sections_new {
                let a_max_dist = arr_distance[i - 1];
                if a_max_dist >= 0.0
                    && a_max_dist <= self.tolerance
                    && arr_len[i - 1] > self.min_tolerance
                {
                    let mut is_inserted = false;
                    let ori = arr_forward[i - 1];
                    let mut j = 1usize;
                    while j <= seq_candidates.len() && !is_inserted {
                        let other = seq_candidates.value(j).unwrap();
                        let a_delta = arr_distance[i - 1] - arr_distance[(other - 1) as usize];
                        if a_delta < CONFUSION
                            && ((a_delta.abs() > REAL_SMALL)
                                || arr_min_dist[i - 1] < arr_min_dist[(other - 1) as usize])
                        {
                            seq_candidates.insert_before(j, i as i32);
                            seq_orientations.insert_before(j, if ori { 1 } else { 0 });
                            is_inserted = true;
                        }
                        j += 1;
                    }
                    if !is_inserted {
                        seq_candidates.append(i as i32);
                        seq_orientations.append(if ori { 1 } else { 0 });
                    }
                }
            }
            nb_candidates = seq_candidates.len();
            if nb_candidates == 0 {
                return false;
            }
            // Replace the candidate indices (`cxx:1691-1696`).
            for i in 1..=nb_candidates {
                let v = seq_candidates.value(i).unwrap();
                let nv = seq_candidates_new.value(v as usize).unwrap_or(v);
                seq_candidates.set(i, nv);
            }
        }
        if nb_candidates == 0 {
            return false;
        }
        if self.nonmanifold && nb_candidates > 1 {
            let mut seq_new_candidates = IntSequence::new();
            let mut seq_orientations_new = IntSequence::new();
            seq_candidates.prepend(1);
            seq_orientations.prepend(1);
            let mut k = 1usize;
            while k <= nb_sections && seq_candidates.len() > 1 {
                self.analysis_nearest_edges(
                    seq_sections,
                    seq_candidates,
                    seq_orientations,
                    k == 1,
                );
                if k == 1 && seq_candidates.len() == 0 {
                    return false;
                }
                if seq_candidates.len() > 0 {
                    seq_new_candidates.append(seq_candidates.first().unwrap());
                    seq_orientations_new.append(seq_orientations.first().unwrap());
                }
                k += 1;
            }
            seq_candidates.prepend_seq(&seq_new_candidates);
            seq_orientations.prepend_seq(&seq_orientations_new);
            return true;
        }
        // For the manifold case leave only one candidate from the
        // equidistant ones (`cxx:1724-1846`).
        while nb_candidates > 0 {
            let ind_candidate = seq_candidates.first().unwrap();
            if map_reference.contains(ind_candidate) {
                break;
            }
            map_reference.add(ind_candidate);
            let mut seq_candidates1 = IntSequence::new();
            let mut seq_orientations1 = IntSequence::new();
            let mut is_found = self.find_candidates(
                seq_sections,
                map_reference,
                &mut seq_candidates1,
                &mut seq_orientations1,
            );
            map_reference.remove_last();
            if is_found {
                is_found = seq_candidates1.len() > 0;
            }
            if is_found {
                let ind_candidate1 = seq_candidates1.first().unwrap();
                if ind_candidate1 == ind_reference {
                    break;
                }
                if map_reference.contains(ind_candidate1) {
                    seq_candidates.prepend(ind_candidate1);
                    nb_candidates += 1;
                    break;
                }
                is_found = false;
            }
            if !is_found {
                seq_candidates.remove(1);
                seq_orientations.remove(1);
                nb_candidates -= 1;
            }
        }
        // gka check (`cxx:1846-1871`).
        if nb_candidates > 0 {
            let an_ind = seq_candidates.value(1).unwrap();
            let edge2 = Edge(seq_sections[(an_ind - 1) as usize].clone());
            let mut bnd = edge2.0.clone();
            if let Some(b) = self.section_bound.get(&bnd) {
                bnd = b;
            }
            if let Some(faces) = self.bound_faces.get(&bnd) {
                let mut is_ok = true;
                for f in &faces {
                    if !is_ok {
                        break;
                    }
                    is_ok = !faces1.contains(f);
                    if !is_ok {
                        let face2 = Face(f.clone());
                        is_ok = self.is_merged_closed(&edge1, &edge2, &face2);
                    }
                }
                if !is_ok {
                    return false;
                }
            }
        }
        nb_candidates > 0
    }
}
