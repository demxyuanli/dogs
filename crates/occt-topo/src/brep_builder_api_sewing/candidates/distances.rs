//! `EvaluateDistances` (BRepBuilderAPI_Sewing.cxx:1287-1522).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::INFINITE;
use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

use super::super::cutting::project_points_on_curve;
use super::super::Sewing;


impl Sewing {
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
}
