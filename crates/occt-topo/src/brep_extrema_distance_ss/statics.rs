//! The typed pair solvers `vertex_*_impl` .. `face_face_impl`.

use super::prelude::*;

use super::perform_c0::*;

/// `Perform(V1, V2, ...)` (`cxx:780-802`).
pub(super) fn vertex_vertex_impl(
    the_s1: &Vertex,
    the_s2: &Vertex,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:787-788`.
    let p1 = BRepTool::vertex_point(the_s1);
    let p2 = BRepTool::vertex_point(the_s2);
    // `cxx:789-801`.
    let dst = p1.distance(&p2);
    if dst < *dst_ref - eps || (dst - *dst_ref).abs() < eps {
        if *dst_ref > dst {
            *dst_ref = dst;
        }
        *modif = true;
        seq1.push(SolutionElem::of_vertex(dst, &p1, SupportType::IsVertex, the_s1));
        seq2.push(SolutionElem::of_vertex(dst, &p2, SupportType::IsVertex, the_s2));
    }
}

/// `Perform(V, E, ...)` (`cxx:806-864`).
pub(super) fn vertex_edge_impl(
    the_s1: &Vertex,
    the_s2: &Edge,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:812-815`.
    if BRepTool::is_degenerated(the_s2) {
        return;
    }
    // `cxx:817-820`.
    let ext = ExtPC::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() { ext.nb_ext() } else { 0 };
    if nb_extrema == 0 {
        return;
    }
    // `cxx:822-831`.
    let mut dstmin = ext.square_distance(1);
    for i in 2..=nb_extrema {
        let s = ext.square_distance(i);
        if s < dstmin {
            dstmin = s;
        }
    }
    dstmin = dstmin.sqrt();
    // `cxx:832-861`.
    if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
        let p1 = BRepTool::vertex_point(the_s1);
        let eps_p = Precision::PCONFUSION;
        for i in 1..=nb_extrema {
            if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                let pt = ext.point(i);
                // `cxx:844`.
                if tri_solution(seq2, &pt) {
                    let t = ext.parameter(i);
                    // `cxx:847`.
                    if !is_on_boundary(the_s2, &pt, t, eps_p) {
                        if *dst_ref > dstmin {
                            *dst_ref = dstmin;
                        }
                        *modif = true;
                        seq1.push(SolutionElem::of_vertex(
                            dstmin,
                            &p1,
                            SupportType::IsVertex,
                            the_s1,
                        ));
                        seq2.push(SolutionElem::of_edge(
                            dstmin,
                            &pt,
                            SupportType::IsOnEdge,
                            the_s2,
                            t,
                        ));
                    }
                }
            }
        }
    }
}

/// `Perform(V, F, ...)` (`cxx:868-925`).
#[allow(clippy::too_many_arguments)]
pub(super) fn vertex_face_impl(
    the_s1: &Vertex,
    the_s2: &Face,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    flag: ExtremaExtFlag,
    algo: ExtremaExtAlgo,
    modif: &mut bool,
) {
    // `cxx:875-876`.
    let ext = ExtPF::new(the_s1, the_s2, flag, algo);
    let nb_extrema = if ext.is_done() { ext.nb_ext() } else { 0 };
    if nb_extrema == 0 {
        return;
    }
    // `cxx:878-888`.
    let mut dstmin = ext.square_distance(1);
    for i in 2..=nb_extrema {
        let s = ext.square_distance(i);
        if s < dstmin {
            dstmin = s;
        }
    }
    dstmin = dstmin.sqrt();
    // `cxx:890-922`.
    if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
        let p1 = BRepTool::vertex_point(the_s1);
        let mut classifier = FaceClassifier::new();
        let tol = BRepTool::face_tolerance(the_s2);
        for i in 1..=nb_extrema {
            if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                let pt = ext.point(i).2;
                // `cxx:902`.
                if tri_solution(seq2, &pt) {
                    let (u, v) = ext.parameter(i);
                    classifier.perform(the_s2, GpPnt2d::new(u, v), tol);
                    // `cxx:908`.
                    if classifier.state() == FaceState::In {
                        if *dst_ref > dstmin {
                            *dst_ref = dstmin;
                        }
                        *modif = true;
                        seq1.push(SolutionElem::of_vertex(
                            dstmin,
                            &p1,
                            SupportType::IsVertex,
                            the_s1,
                        ));
                        seq2.push(SolutionElem::of_face(
                            dstmin,
                            &pt,
                            SupportType::IsInFace,
                            the_s2,
                            u,
                            v,
                        ));
                    }
                }
            }
        }
    }
}

/// `Perform(E1, E2, ...)` (`cxx:929-1023`).
pub(super) fn edge_edge_impl(
    the_s1: &Edge,
    the_s2: &Edge,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:935-938`.
    if BRepTool::is_degenerated(the_s1) || BRepTool::is_degenerated(the_s2) {
        return;
    }
    // `cxx:940`: the reference captured before the extremum search.
    let dst_ref0 = *dst_ref;

    // `cxx:942-943`.
    let ext = ExtCC::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() {
        if ext.is_parallel() {
            0
        } else {
            ext.nb_ext()
        }
    } else {
        0
    };
    // `cxx:944-989`.
    if nb_extrema > 0 {
        let mut dstmin = ext.square_distance(1);
        for i in 2..=nb_extrema {
            let s = ext.square_distance(i);
            if s < dstmin {
                dstmin = s;
            }
        }
        dstmin = dstmin.sqrt();
        if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
            let eps_p = Precision::PCONFUSION;
            for i in 1..=nb_extrema {
                if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                    let pt1 = ext.point_on_e1(i);
                    let pt2 = ext.point_on_e2(i);
                    // `cxx:969`.
                    if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
                        let t1 = ext.parameter_on_e1(i);
                        let t2 = ext.parameter_on_e2(i);
                        // `cxx:973`.
                        if !is_on_boundary(the_s1, &pt1, t1, eps_p)
                            && !is_on_boundary(the_s2, &pt2, t2, eps_p)
                        {
                            if *dst_ref > dstmin {
                                *dst_ref = dstmin;
                            }
                            *modif = true;
                            seq1.push(SolutionElem::of_edge(
                                dstmin,
                                &pt1,
                                SupportType::IsOnEdge,
                                the_s1,
                                t1,
                            ));
                            seq2.push(SolutionElem::of_edge(
                                dstmin,
                                &pt2,
                                SupportType::IsOnEdge,
                                the_s2,
                                t2,
                            ));
                        }
                    }
                }
            }
        }
    }

    // `cxx:992-1022`.
    let mut seq_solution1 = Vec::new();
    let mut seq_solution2 = Vec::new();
    perform_c0(
        the_s1,
        the_s2,
        &mut seq_solution1,
        &mut seq_solution2,
        dst_ref0,
        dst_ref,
        eps,
    );
    let (sol1, sol2) = min_solution(&seq_solution1, &seq_solution2, *dst_ref, eps);
    for (s1, s2) in sol1.iter().zip(sol2.iter()) {
        let pt1 = *s1.point();
        let pt2 = *s2.point();
        // `cxx:1014`.
        if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
            seq1.push(s1.clone());
            seq2.push(s2.clone());
            *modif = true;
        }
    }
}

/// `Perform(E, F, ...)` (`cxx:1026-1187`).
pub(super) fn edge_face_impl(
    the_s1: &Edge,
    the_s2: &Face,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:1032-1035`.
    if BRepTool::is_degenerated(the_s1) {
        return;
    }

    let mut classifier = FaceClassifier::new();

    // `cxx:1039-1040`.
    let ext = ExtCF::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() {
        if ext.is_parallel() {
            0
        } else {
            ext.nb_ext()
        }
    } else {
        0
    };
    // `cxx:1041-1095`.
    if nb_extrema > 0 {
        let mut dstmin = ext.square_distance(1);
        for i in 2..=nb_extrema {
            let s = ext.square_distance(i);
            if s < dstmin {
                dstmin = s;
            }
        }
        dstmin = dstmin.sqrt();
        if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
            let tol = BRepTool::face_tolerance(the_s2);
            let eps_p = Precision::PCONFUSION;
            for i in 1..=nb_extrema {
                if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                    let pt1 = ext.point_on_edge(i);
                    let pt2 = ext.point_on_face(i);
                    // `cxx:1069`.
                    if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
                        let t1 = ext.parameter_on_edge(i);
                        // `cxx:1072`.
                        if !is_on_boundary(the_s1, &pt1, t1, eps_p) {
                            let (u, v) = ext.parameter_on_face(i);
                            classifier.perform(the_s2, GpPnt2d::new(u, v), tol);
                            // `cxx:1079`.
                            if classifier.state() == FaceState::In {
                                if *dst_ref > dstmin {
                                    *dst_ref = dstmin;
                                }
                                *modif = true;
                                seq1.push(SolutionElem::of_edge(
                                    dstmin,
                                    &pt1,
                                    SupportType::IsOnEdge,
                                    the_s1,
                                    t1,
                                ));
                                seq2.push(SolutionElem::of_face(
                                    dstmin,
                                    &pt2,
                                    SupportType::IsInFace,
                                    the_s2,
                                    u,
                                    v,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    // `cxx:1097-1186` (OCCT 8.0.1 #1382): null 3D curve aborts the C0 arm.
    let (a_first, a_last) = BRepTool::edge_parameters(the_s1);
    let p_curv = match BRepTool::edge_curve(the_s1) {
        Some(c) => c,
        None => return,
    };
    // `cxx:1107`.
    if p_curv.continuity() != 0 {
        return;
    }

    // `cxx:1100-1110`.
    let arr_inter = extrema_pc::extrema_intervals(p_curv.as_ref(), 1);
    let tol = BRepTool::face_tolerance(the_s2);

    let mut seq_solution1 = Vec::new();
    let mut seq_solution2 = Vec::new();
    let b = TopoBuilder::new();
    for &a_parameter in &arr_inter {
        // `cxx:1117-1119`.
        let a_pnt = p_curv.d0(a_parameter);
        let v1 = b.make_vertex(a_pnt, Precision::CONFUSION);
        // `cxx:1121`.
        let ext_pf = ExtPF::new(
            &v1,
            the_s2,
            ExtremaExtFlag::MinMax,
            ExtremaExtAlgo::Grad,
        );
        let nb_extrema_pf = if ext_pf.is_done() { ext_pf.nb_ext() } else { 0 };
        if nb_extrema_pf == 0 {
            continue;
        }
        // `cxx:1125-1134`.
        let mut dstmin = ext_pf.square_distance(1);
        for ii in 2..=nb_extrema_pf {
            let s = ext_pf.square_distance(ii);
            if s < dstmin {
                dstmin = s;
            }
        }
        dstmin = dstmin.sqrt();
        // `cxx:1136-1171`.
        if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
            for ii in 1..=nb_extrema_pf {
                if (dstmin - ext_pf.square_distance(ii).sqrt()).abs() < eps {
                    let (u, v) = ext_pf.parameter(ii);
                    classifier.perform(the_s2, GpPnt2d::new(u, v), tol);
                    // `cxx:1146`.
                    if classifier.state() == FaceState::In {
                        if *dst_ref > dstmin {
                            *dst_ref = dstmin;
                        }
                        *modif = true;
                        seq_solution1.push(SolutionElem::of_edge(
                            dstmin,
                            &a_pnt,
                            SupportType::IsOnEdge,
                            the_s1,
                            a_parameter,
                        ));
                        seq_solution2.push(SolutionElem::of_face(
                            dstmin,
                            &ext_pf.point(ii).2,
                            SupportType::IsInFace,
                            the_s2,
                            u,
                            v,
                        ));
                    }
                }
            }
        }
    }

    // `cxx:1173-1186`.
    let (sol1, sol2) = min_solution(&seq_solution1, &seq_solution2, *dst_ref, eps);
    for (s1, s2) in sol1.iter().zip(sol2.iter()) {
        seq1.push(s1.clone());
        seq2.push(s2.clone());
    }
    let _ = (a_first, a_last);
}

/// `Perform(F1, F2, ...)` (`cxx:1191-1206`).
pub(super) fn face_face_impl(
    the_s1: &Face,
    the_s2: &Face,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:1197-1198`.
    let ext = ExtFF::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() {
        if ext.is_parallel() {
            0
        } else {
            ext.nb_ext()
        }
    } else {
        0
    };
    if nb_extrema == 0 {
        return;
    }
    // `cxx:1201-1210`.
    let mut dstmin = ext.square_distance(1);
    for i in 2..=nb_extrema {
        let s = ext.square_distance(i);
        if s < dstmin {
            dstmin = s;
        }
    }
    dstmin = dstmin.sqrt();
    // `cxx:1211-1261`.
    if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
        let tol1 = BRepTool::face_tolerance(the_s1);
        let tol2 = BRepTool::face_tolerance(the_s2);
        let mut classifier = FaceClassifier::new();
        for i in 1..=nb_extrema {
            if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                let pt1 = ext.point_on_face1(i);
                let pt2 = ext.point_on_face2(i);
                // `cxx:1229`.
                if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
                    let (u1, v1) = ext.parameter_on_face1(i);
                    classifier.perform(the_s1, GpPnt2d::new(u1, v1), tol1);
                    // `cxx:1234`.
                    if classifier.state() == FaceState::In {
                        let (u2, v2) = ext.parameter_on_face2(i);
                        classifier.perform(the_s2, GpPnt2d::new(u2, v2), tol2);
                        // `cxx:1240`.
                        if classifier.state() == FaceState::In {
                            if *dst_ref > dstmin {
                                *dst_ref = dstmin;
                            }
                            *modif = true;
                            seq1.push(SolutionElem::of_face(
                                dstmin,
                                &pt1,
                                SupportType::IsInFace,
                                the_s1,
                                u1,
                                v1,
                            ));
                            seq2.push(SolutionElem::of_face(
                                dstmin,
                                &pt2,
                                SupportType::IsInFace,
                                the_s2,
                                u2,
                                v2,
                            ));
                        }
                    }
                }
            }
        }
    }
}
