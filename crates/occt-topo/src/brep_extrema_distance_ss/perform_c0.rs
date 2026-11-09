//! `TRI_SOLUTION`, `MIN_SOLUTION`, `isOnBoundary` and `PERFORM_C0`.

use super::prelude::*;


/// `TRI_SOLUTION` (`BRepExtrema_DistanceSS.cxx:60-75`): true when no already
/// recorded solution is within `Precision::Confusion()` of `pt`.
pub(super) fn tri_solution(seq_sol: &[SolutionElem], pt: &GpPnt) -> bool {
    for s in seq_sol {
        if s.point().distance(pt) <= Precision::CONFUSION {
            return false;
        }
    }
    true
}

/// `MIN_SOLUTION` (`cxx:79-97`): keep the pairs whose distance equals `dst_ref`
/// within `eps`.
pub(super) fn min_solution(
    seq_sol1: &[SolutionElem],
    seq_sol2: &[SolutionElem],
    dst_ref: f64,
    eps: f64,
) -> (Vec<SolutionElem>, Vec<SolutionElem>) {
    let mut out1 = Vec::new();
    let mut out2 = Vec::new();
    for (s1, s2) in seq_sol1.iter().zip(seq_sol2.iter()) {
        if (s1.dist() - dst_ref).abs() < eps {
            out1.push(s1.clone());
            out2.push(s2.clone());
        }
    }
    (out1, out2)
}

/// `isOnBoundary` (`cxx:637-653`): true when `(the_sol, the_param)` is the
/// edge's vertex at that parameter, within the vertex tolerance.
pub(super) fn is_on_boundary(the_edge: &Edge, the_sol: &GpPnt, the_param: f64, the_ptol: f64) -> bool {
    let kids = the_edge.0.tshape.read().unwrap().children.clone();
    for k in kids {
        if k.shape_type() != ShapeType::Vertex {
            continue;
        }
        let v = match Vertex::wrap(k.clone()) {
            Some(v) => v,
            None => continue,
        };
        let v_param = BRepTool::parameter_of_vertex(&the_edge.0, &k);
        if (v_param - the_param).abs() < the_ptol
            && BRepTool::vertex_point(&v).distance(the_sol) < BRepTool::vertex_tolerance(&v)
        {
            return true;
        }
    }
    false
}

/// `PERFORM_C0` (`cxx:488-630`).
#[allow(clippy::too_many_arguments)]
pub(super) fn perform_c0(
    s1: &Edge,
    s2: &Edge,
    seq_sol1: &mut Vec<SolutionElem>,
    seq_sol2: &mut Vec<SolutionElem>,
    dst_ref: f64,
    m_dst_ref: &mut f64,
    eps: f64,
) {
    // `cxx:495-498`.
    if BRepTool::is_degenerated(s1) || BRepTool::is_degenerated(s2) {
        return;
    }

    // `cxx:500-513`.
    for i_e in 0..2 {
        let (e, e_other) = if i_e == 0 { (s1, s2) } else { (s2, s1) };
        // `cxx:515-518`.
        let (a_first, a_last) = BRepTool::edge_parameters(e);
        let p_curv = match BRepTool::edge_curve(e) {
            Some(c) => c,
            None => continue,
        };
        let (a_f_other, a_l_other) = BRepTool::edge_parameters(e_other);
        let p_curv_other = match BRepTool::edge_curve(e_other) {
            Some(c) => c,
            None => continue,
        };

        // `cxx:520`.
        if p_curv.continuity() != 0 {
            continue;
        }

        // `cxx:522-536`.
        let eps_p = Precision::PCONFUSION;
        let arr_inter = extrema_pc::extrema_intervals(p_curv.as_ref(), 1);
        let arr_inter_other = extrema_pc::extrema_intervals(p_curv_other.as_ref(), 1);
        // `cxx:537-538`.
        let (u_deb, u_fin) = (a_f_other, a_l_other);

        // `cxx:541-628`.
        for &a_parameter in &arr_inter {
            let a_pnt = p_curv.d0(a_parameter);
            // `cxx:544`.
            let b = TopoBuilder::new();
            let v1 = b.make_vertex(a_pnt, Precision::CONFUSION);

            // `cxx:547-548`.
            let ext = ExtPC::new(&v1, e_other);
            let nb_extrema = if ext.is_done() { ext.nb_ext() } else { 0 };
            if nb_extrema > 0 {
                // `cxx:553-560`.
                let mut dstmin = ext.square_distance(1);
                for ii in 2..=nb_extrema {
                    let s_dst = ext.square_distance(ii);
                    if s_dst < dstmin {
                        dstmin = s_dst;
                    }
                }
                dstmin = dstmin.sqrt();

                // `cxx:563-604`.
                if dstmin < dst_ref - eps || (dstmin - dst_ref).abs() < eps {
                    for ii in 1..=nb_extrema {
                        if (dstmin - ext.square_distance(ii).sqrt()).abs() < eps {
                            let pt = ext.point(ii);
                            // `cxx:570-578`: the recorded point lies on
                            // `Eother`, so the already-known solutions are the
                            // mirror sequence.
                            let tri = if i_e == 0 {
                                tri_solution(seq_sol2, &pt)
                            } else {
                                tri_solution(seq_sol1, &pt)
                            };
                            if tri {
                                // `cxx:581-582`.
                                let t = ext.parameter(ii);
                                if (t - u_deb).abs() >= eps_p && (t - u_fin).abs() > eps_p {
                                    if *m_dst_ref > dstmin {
                                        *m_dst_ref = dstmin;
                                    }
                                    // `cxx:588-599`.
                                    let sol1 = SolutionElem::of_edge(
                                        dstmin,
                                        &a_pnt,
                                        SupportType::IsOnEdge,
                                        e,
                                        a_parameter,
                                    );
                                    let sol2 = SolutionElem::of_edge(
                                        dstmin,
                                        &pt,
                                        SupportType::IsOnEdge,
                                        e_other,
                                        t,
                                    );
                                    if i_e == 0 {
                                        seq_sol1.push(sol1);
                                        seq_sol2.push(sol2);
                                    } else {
                                        seq_sol1.push(sol2);
                                        seq_sol2.push(sol1);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // `cxx:606-627`.
            for &a_parameter_other in &arr_inter_other {
                let a_pnt_other = p_curv_other.d0(a_parameter_other);
                let dst = a_pnt.distance(&a_pnt_other);
                if dst < dst_ref - eps || (dst - dst_ref).abs() < eps {
                    if *m_dst_ref > dst {
                        *m_dst_ref = dst;
                    }
                    let sol1 = SolutionElem::of_edge(
                        dst,
                        &a_pnt,
                        SupportType::IsOnEdge,
                        e,
                        a_parameter,
                    );
                    let sol2 = SolutionElem::of_edge(
                        dst,
                        &a_pnt_other,
                        SupportType::IsOnEdge,
                        e_other,
                        a_parameter_other,
                    );
                    if i_e == 0 {
                        seq_sol1.push(sol1);
                        seq_sol2.push(sol2);
                    } else {
                        seq_sol1.push(sol2);
                        seq_sol2.push(sol1);
                    }
                }
            }
            let _ = a_first;
            let _ = a_last;
            let _ = a_f_other;
            let _ = a_l_other;
        }
    }
}
