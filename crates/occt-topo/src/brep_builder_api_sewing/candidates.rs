//! Candidate search and distance evaluation of `BRepBuilderAPI_Sewing`:
//! `EvaluateAngulars` (`cxx:1189-1280`), `EvaluateDistances`
//! (`cxx:1287-1522`), `AnalysisNearestEdges` (`cxx:1618-1782`),
//! `FindCandidates` (`cxx:1785-2086`) and `MergedNearestEdges`
//! (`cxx:4242-4431`).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, INFINITE, REAL_SMALL, RESOLUTION};
use occt_geom::gcpnts::{abscissa_point_with_tolerance, curve_length_range};
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};

use super::cutting::project_points_on_curve;
use super::shape_maps::{IntIndexedMap, IntSequence, ShapeSet};
use super::Sewing;

/// `GCPnts_UniformAbscissa(adaptor, npt, first, last)`
/// (`GCPnts_UniformAbscissa.cxx:505-515`: `anAbscissa = L / (NbPoints - 1)`).
/// `None` when the arc-length inversion does not converge, which OCCT answers
/// with `IsDone() == false`.
#[allow(dead_code)] // OCCT calls it only from EvaluateAngulars, whose only caller is commented out.
pub(crate) fn uniform_abscissa_range(
    c: &dyn Curve,
    npt: usize,
    first: f64,
    last: f64,
) -> Option<Vec<f64>> {
    if npt < 2 {
        return Some(vec![first]);
    }
    let total = curve_length_range(c, first, last, CONFUSION);
    if total <= CONFUSION {
        return Some(vec![first; npt]);
    }
    let step = total / (npt - 1) as f64;
    let res = c.resolution(CONFUSION);
    let mut out = Vec::with_capacity(npt);
    out.push(first);
    let mut prev = first;
    for _ in 1..npt - 1 {
        let u = abscissa_point_with_tolerance(c, step, prev, res).ok()?;
        let u = u.min(last);
        out.push(u);
        prev = u;
    }
    out.push(last);
    Some(out)
}

impl Sewing {
    /// `EvaluateAngulars(sequenceSec, secForward, tabAng, indRef)`
    /// (`cxx:1189-1280`): the mean angle between the surface normal of the
    /// reference section and each following one, sampled at `4` uniform
    /// abscissa points.
    #[allow(dead_code)] // cxx:1769 is commented out: the only OCCT caller is dead source.
    pub(crate) fn evaluate_angulars(
        &self,
        sequence_sec: &[TopoShape],
        sec_forward: &[bool],
        tab_ang: &mut [f64],
        ind_ref: usize,
    ) {
        for a in tab_ang.iter_mut() {
            *a = -1.0;
        }
        let npt = 4usize;
        let leng_sec = sequence_sec.len();
        let mut norm_ref: Vec<GpVec> = vec![GpVec::zero(); npt];
        for i in ind_ref..=leng_sec {
            let edge = sequence_sec[i - 1].clone();
            let mut bnd = edge.clone();
            if let Some(b) = self.section_bound.get(&edge) {
                bnd = b;
            }
            let (face, c2d_opt) = if let Some(faces) = self.bound_faces.get(&bnd) {
                let face = Face(faces[0].clone());
                let c2d = crate::boptools_2d::curve_on_surface_range(&Edge(edge.clone()), &face)
                    .map(|(c, _, _)| c);
                (Some(face), c2d)
            } else if i == ind_ref {
                return;
            } else {
                continue;
            };
            let Some(face) = face else { continue };
            let Some(c2d) = c2d_opt else { continue };
            let Some(surf) = self.face_surface_located(&face) else {
                continue;
            };
            let Some(c3d) = BRepTool::edge_curve_adaptor_world(&Edge(edge.clone())) else {
                continue;
            };
            let (first, last) = BRepTool::edge_parameters(&Edge(edge));
            let Some(params) = uniform_abscissa_range(c3d.as_ref(), npt, first, last) else {
                continue;
            };
            let mut cumulate_angular = 0.0;
            let mut nb_computed = 0usize;
            for j in 1..=npt {
                let jj = if sec_forward.get(i - 1).copied().unwrap_or(false) || i == ind_ref {
                    j
                } else {
                    npt - j + 1
                };
                let p = c2d.d0(params[jj - 1]);
                let (_pnt, w1, w2) = surf.d1(p.x(), p.y());
                let n = w1.crossed(&w2);
                if i == ind_ref {
                    norm_ref[j - 1] = n;
                } else if n.magnitude() > RESOLUTION && norm_ref[j - 1].magnitude() > RESOLUTION {
                    nb_computed += 1;
                    let cos = n.dot(&norm_ref[j - 1]) / (n.magnitude() * norm_ref[j - 1].magnitude());
                    let mut angular = cos.clamp(-1.0, 1.0).acos();
                    if angular > std::f64::consts::FRAC_PI_2 {
                        angular = std::f64::consts::PI - angular;
                    }
                    cumulate_angular += angular;
                }
            }
            if nb_computed > 0 {
                tab_ang[i - 1] = cumulate_angular / nb_computed as f64;
            }
        }
    }

    /// `EvaluateDistances(sequenceSec, secForward, tabDst, arrLen, tabMinDist,
    /// indRef)` (`cxx:1287-1522`): mutual orientation and maximal deviation of
    /// every section against the reference one.
    ///
    /// The tail of the OCCT function (the commented-out `Project distant
    /// points` block, `cxx:1303-1361`) is dead source and not translated.
    pub(crate) fn evaluate_distances(
        &self,
        sequence_sec: &[TopoShape],
        sec_forward: &mut [bool],
        tab_dst: &mut [f64],
        arr_len: &mut [f64],
        tab_min_dist: &mut [f64],
        ind_ref: usize,
    ) {
        for v in sec_forward.iter_mut() {
            *v = true;
        }
        for v in tab_dst.iter_mut() {
            *v = -1.0;
        }
        for v in arr_len.iter_mut() {
            *v = 0.0;
        }
        for v in tab_min_dist.iter_mut() {
            *v = INFINITE;
        }
        let npt = 8usize;
        let leng_sec = sequence_sec.len();
        let mut pts_ref: Vec<GpPnt> = vec![GpPnt::zero(); npt];
        let mut pts_sec: Vec<GpPnt> = vec![GpPnt::zero(); npt];
        let mut c3d_ref: Option<std::sync::Arc<dyn Curve>> = None;
        let mut first_ref = 0.0f64;
        let mut last_ref = 0.0f64;

        for i in ind_ref..=leng_sec {
            let sec = Edge(sequence_sec[i - 1].clone());
            let (first, last) = BRepTool::edge_parameters(&sec);
            let Some(c3d) = BRepTool::edge_curve_adaptor_world(&sec) else {
                continue;
            };
            if i == ind_ref {
                c3d_ref = Some(c3d.clone());
                first_ref = first;
                last_ref = last;
            }
            let mut dist: f64;
            let mut dist_for = -1.0;
            let mut dist_rev = -1.0;
            let mut a_min_dist = INFINITE;
            let delta_t = (last - first) / (npt - 1) as f64;
            let mut a_len_sec2 = 0.0;
            let mut nb_found = 0usize;
            for j in 1..=npt {
                let t = if j == 1 {
                    first
                } else if j == npt {
                    last
                } else {
                    first + (j - 1) as f64 * delta_t
                };
                let pt = c3d.d0(t);
                if i == ind_ref {
                    pts_ref[j - 1] = pt;
                    if j > 1 {
                        a_len_sec2 += pt.square_distance(&pts_ref[j - 2]);
                    }
                } else {
                    pts_sec[j - 1] = pt;
                    if j > 1 {
                        a_len_sec2 += pt.square_distance(&pts_sec[j - 2]);
                    }
                    dist = pt.distance(&pts_ref[j - 1]);
                    if a_min_dist > dist {
                        a_min_dist = dist;
                    }
                    if dist_for < dist {
                        dist_for = dist;
                    }
                    dist = pt.distance(&pts_ref[npt - j]);
                    if a_min_dist > dist {
                        a_min_dist = dist;
                    }
                    if dist_rev < dist {
                        dist_rev = dist;
                    }
                    // Check that the point lies between the endpoints of the
                    // reference curve (`cxx:1245-1249`).
                    let p11 = pts_ref[0];
                    let p12 = pts_ref[npt - 1];
                    let a_vec1 = GpVec::from_pnts(&pt, &p11);
                    let a_vec2 = GpVec::from_pnts(&pt, &p12);
                    let a_vec_ref = GpVec::from_pnts(&p11, &p12);
                    if (a_vec_ref.dot(&a_vec1)) * (a_vec_ref.dot(&a_vec2)) < 0.0 {
                        nb_found += 1;
                    }
                }
            }
            let a_len_sec = a_len_sec2.sqrt();
            arr_len[i - 1] = a_len_sec;
            let is_forward = dist_for < dist_rev;
            sec_forward[i - 1] = is_forward;
            dist = if is_forward { dist_for } else { dist_rev };
            if i == ind_ref || (dist < self.tolerance && nb_found * 2 >= npt) {
                tab_dst[i - 1] = dist;
                tab_min_dist[i - 1] = a_min_dist;
            } else {
                nb_found = 0;
                a_min_dist = INFINITE;
                dist = -1.0;
                let mut arr_dist = vec![-1.0f64; npt];
                let mut arr_para = vec![0.0f64; npt];
                let mut arr_proj = vec![GpPnt::zero(); npt];
                if arr_len[ind_ref - 1] >= arr_len[i - 1] {
                    let c3d_ref = c3d_ref.as_ref().unwrap();
                    project_points_on_curve(
                        &pts_sec,
                        c3d_ref.as_ref(),
                        first_ref,
                        last_ref,
                        self.tolerance,
                        false,
                        &mut arr_dist,
                        &mut arr_para,
                        &mut arr_proj,
                    );
                } else {
                    project_points_on_curve(
                        &pts_ref,
                        c3d.as_ref(),
                        first,
                        last,
                        self.tolerance,
                        false,
                        &mut arr_dist,
                        &mut arr_para,
                        &mut arr_proj,
                    );
                }
                for j in 0..npt {
                    if arr_dist[j] < 0.0 {
                        continue;
                    }
                    if dist < arr_dist[j] {
                        dist = arr_dist[j];
                    }
                    if a_min_dist > arr_dist[j] {
                        a_min_dist = arr_dist[j];
                    }
                    nb_found += 1;
                }
                if nb_found > 1 {
                    tab_dst[i - 1] = dist;
                    tab_min_dist[i - 1] = a_min_dist;
                }
            }
        }
    }

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

/// `TopExp::Vertices(edge)` helper kept local so the module does not depend on
/// the `Shape` wrappers of the caller.
#[allow(dead_code)]
fn edge_end_points(edge: &TopoShape) -> Option<(TopoShape, TopoShape)> {
    match crate::topo_tools::edge_vertices(&Edge(edge.clone())) {
        (Some(a), Some(b)) => Some((a.0, b.0)),
        _ => None,
    }
}

#[allow(dead_code)]
fn face_of(shape: &TopoShape) -> Option<Face> {
    if shape.shape_type() != ShapeType::Face {
        return None;
    }
    Face::wrap(shape.clone())
}


