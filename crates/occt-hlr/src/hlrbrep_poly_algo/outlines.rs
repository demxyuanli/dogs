//! Outline maintenance passes. Source: `HLRBRep_PolyAlgo::InsertOnOutLine`
//! (`cxx:2980-3161`), `CheckFrBackTriangles` (`cxx:3163-3511`) and
//! `UpdateOutLines` (`cxx:3998-4153`).

use occt_core::gp::GpXyz;

use super::{nmask, HlrBRepPolyAlgo, PolyInternalDataArray};
use crate::hlr_algo::{
    poly_mask, BiPointIndices, HlrAlgoBiPoint, HlrAlgoTriangleData,
};

/// `aTriangle.Flags |= HLRAlgo_PolyMask_EMskOutLin1/2/3` (`cxx:4071-4081`,
/// `:4100-4110`): mark which of the three edges of `tri` is the outline
/// segment `pd -> pf`. The node swap reproduces the `FMskOrBack` reordering
/// OCCT does inline before the three tests.
fn mark_out_line_edge(tri: &mut HlrAlgoTriangleData, pd: i32, pf: i32) {
    let mut tn1 = tri.node1;
    let tn2 = tri.node2;
    let mut tn3 = tri.node3;
    if tri.flags & poly_mask::FMSK_SIDE == 0 && tri.flags & poly_mask::FMSK_OR_BACK != 0 {
        std::mem::swap(&mut tn1, &mut tn3);
    }
    if (tn1 == pd && tn2 == pf) || (tn1 == pf && tn2 == pd) {
        tri.flags |= poly_mask::EMSK_OUT_LIN1;
    } else if (tn2 == pd && tn3 == pf) || (tn2 == pf && tn3 == pd) {
        tri.flags |= poly_mask::EMSK_OUT_LIN2;
    } else if (tn3 == pd && tn1 == pf) || (tn3 == pf && tn1 == pd) {
        tri.flags |= poly_mask::EMSK_OUT_LIN3;
    }
}

/// `isOutl` from the pair of triangles sharing a segment (`cxx:4033-4056`,
/// also `cxx:4204-4223`): the two `FMskSide` / `FMskFlat` / `FMskBack`
/// combinations that make the edge an outline.
fn is_out_line_edge(a: &HlrAlgoTriangleData, b: &HlrAlgoTriangleData) -> bool {
    let a_side = a.flags & poly_mask::FMSK_SIDE != 0;
    let b_side = b.flags & poly_mask::FMSK_SIDE != 0;
    if !a_side && !b_side {
        (a.flags & poly_mask::FMSK_BACK) != (b.flags & poly_mask::FMSK_BACK)
    } else if a_side && b_side {
        false
    } else if a_side {
        a.flags & poly_mask::FMSK_FLAT == 0 && b.flags & poly_mask::FMSK_BACK == 0
    } else {
        b.flags & poly_mask::FMSK_FLAT == 0 && a.flags & poly_mask::FMSK_BACK == 0
    }
}

impl HlrBRepPolyAlgo {
    /// `InsertOnOutLine(thePID)` (`cxx:2980-3161`): split every segment whose
    /// two nodes straddle the outline, moving a node onto the crossing when the
    /// crossing is close to one end and inserting a middle node otherwise.
    pub fn insert_on_out_line(&mut self, the_pid: &mut PolyInternalDataArray) {
        let a_nb_faces = self.f_map.extent();
        let tol_angular = self.tol_angular;
        let tol_sta = self.tol_sta;
        let tol_end = self.tol_end;

        // `cxx:2988-2989`: these out-parameters live outside the face loop.
        let mut m_p3_p1 = false;
        let mut a_coef3 = 0.0;
        let (mut x3, mut y3, mut z3) = (0.0, 0.0, 0.0);

        for a_face_iter in 1..=a_nb_faces {
            if the_pid[a_face_iter].is_none() {
                continue;
            }
            let Some(face) = self.f_map.find_key(a_face_iter).cloned() else {
                continue;
            };
            // `cxx:3007-3028`: `myBSurf.Initialize(aFace, false)` only feeds
            // `Normal` in OCCT; this function's sole surface use is
            // `myGSurf->Value(aU3, aV3)` below. `myGSurf` is
            // `BRep_Tool::Surface(aFace, aLoc)`, i.e. the *local* surface, which
            // is exactly what the registry stores.
            self.read_face_matrices(&face.location().transformation());
            let surf = occt_topo::tgeometry::GeometryRegistry::global().face_surface(&face);

            let pid = the_pid[a_face_iter].as_mut().unwrap();
            let mut is_int_out_l = false;
            let a_nb_segs = pid.nb_pi_seg();
            for a_seg_iter in 1..=a_nb_segs {
                let a_seg = pid.pi_seg()[a_seg_iter as usize];
                if a_seg.conex1 == 0 || a_seg.conex2 == 0 {
                    continue;
                }
                let ip1 = a_seg.lst_sg1;
                let ip2 = a_seg.lst_sg2;
                let nod1 = *pid.pi_nod()[ip1 as usize].data();
                let nod2 = *pid.pi_nod()[ip2 as usize].data();
                let flag1 = pid.pi_nod()[ip1 as usize].indices().flag;
                let flag2 = pid.pi_nod()[ip2 as usize].indices().flag;
                let straddles = (nod1.scal >= tol_angular && nod2.scal <= -tol_angular)
                    || (nod2.scal >= tol_angular && nod1.scal <= -tol_angular);
                if flag1 & nmask::OUT_L != 0 && flag2 & nmask::OUT_L != 0 {
                    is_int_out_l = true;
                } else if straddles {
                    is_int_out_l = true;
                    let mut is_ins_p3 = false;
                    if let Some((coef, move_p1)) =
                        Self::new_node(nod1.scal, nod2.scal, tol_angular)
                    {
                        a_coef3 = coef;
                        m_p3_p1 = move_p1;
                        is_ins_p3 = true;
                    }
                    if is_ins_p3 {
                        let (u3, v3) = Self::uv_node(
                            (nod1.uv.x(), nod1.uv.y()),
                            (nod2.uv.x(), nod2.uv.y()),
                            a_coef3,
                        );
                        match surf.as_ref() {
                            Some(s) => {
                                let p = s.d0(u3, v3);
                                x3 = p.x();
                                y3 = p.y();
                                z3 = p.z();
                            }
                            None => {
                                // `cxx:3074-3081`: simple averaging, the OCCT
                                // fallback when the face has no surface.
                                let a_coef2 = 1.0 - a_coef3;
                                x3 = a_coef2 * nod1.point.x() + a_coef3 * nod2.point.x();
                                y3 = a_coef2 * nod1.point.y() + a_coef3 * nod2.point.y();
                                z3 = a_coef2 * nod1.point.z() + a_coef3 * nod2.point.z();
                            }
                        }
                        self.ttmultiply(&mut x3, &mut y3, &mut z3, false);
                    }

                    if is_ins_p3 && m_p3_p1 {
                        // P1 ---> P3 (`cxx:3085-3105`).
                        if pid.pi_nod()[ip1 as usize].indices().flag & nmask::EDGE == 0
                            && a_coef3 < tol_sta
                        {
                            is_ins_p3 = false;
                            self.change_node(pid, ip1, ip2, a_coef3, x3, y3, z3, true);
                            pid.change_pi_nod()[ip1 as usize].change_data().scal = 0.0;
                            pid.change_pi_nod()[ip1 as usize].change_indices().flag |=
                                nmask::OUT_L;
                        }
                    }
                    if is_ins_p3 && !m_p3_p1 {
                        // P2 ---> P3 (`cxx:3107-3127`).
                        if pid.pi_nod()[ip2 as usize].indices().flag & nmask::EDGE == 0
                            && a_coef3 > tol_end
                        {
                            is_ins_p3 = false;
                            self.change_node(pid, ip1, ip2, a_coef3, x3, y3, z3, false);
                            pid.change_pi_nod()[ip2 as usize].change_data().scal = 0.0;
                            pid.change_pi_nod()[ip2 as usize].change_indices().flag |=
                                nmask::OUT_L;
                        }
                    }
                    if is_ins_p3 {
                        // p1 ip3 p2 (`cxx:3130-3141`).
                        let n1 = *pid.pi_nod()[ip1 as usize].data();
                        let n2 = *pid.pi_nod()[ip2 as usize].data();
                        let ip3 = pid.add_node(&n1, &n2, a_coef3, x3, y3, z3);
                        pid.update_links_split(ip1, ip2, ip3);
                        self.update_around_node(pid, ip3);
                        pid.change_pi_nod()[ip3 as usize].change_data().scal = 0.0;
                        pid.change_pi_nod()[ip3 as usize].change_indices().flag |= nmask::OUT_L;
                    }
                }
            }
            if is_int_out_l {
                pid.set_int_out_l(true);
            }
        }
    }

    /// `CheckFrBackTriangles(theList, thePID)` (`cxx:3163-3511`): nudge the
    /// longest edge of every front/back triangle outwards in XY (at most four
    /// passes), then publish the moved nodes into the bi-point list.
    pub fn check_fr_back_triangles(
        &mut self,
        the_list: &mut [HlrAlgoBiPoint],
        the_pid: &mut PolyInternalDataArray,
    ) {
        let a_nb_faces = self.f_map.extent();
        // `cxx:3166-3167`: these persist across triangles and faces.
        let (mut x1, mut y1, mut x2, mut y2) = (0.0, 0.0, 0.0, 0.0);
        let (mut x3, mut y3) = (0.0, 0.0);
        let mut is_modif = true;
        let mut i_loop = 0;
        while is_modif && i_loop < 4 {
            i_loop += 1;
            is_modif = false;
            let mut is_fr_back_in_list = false;
            for a_face_iter in 1..=a_nb_faces {
                if the_pid[a_face_iter].is_none() {
                    continue;
                }
                let pid = the_pid[a_face_iter].as_mut().unwrap();
                let a_nb_tris = pid.nb_t_data();
                for a_tri_iter in 1..=a_nb_tris {
                    let tdata = pid.t_data()[a_tri_iter as usize];
                    if tdata.flags & poly_mask::FMSK_SIDE != 0
                        || tdata.flags & poly_mask::FMSK_FR_BACK == 0
                    {
                        continue;
                    }
                    is_modif = true;
                    let n1 = *pid.pi_nod()[tdata.node1 as usize].data();
                    let n2 = *pid.pi_nod()[tdata.node2 as usize].data();
                    let n3 = *pid.pi_nod()[tdata.node3 as usize].data();
                    let f1 = pid.pi_nod()[tdata.node1 as usize].indices().flag;
                    let f2 = pid.pi_nod()[tdata.node2 as usize].indices().flag;
                    let f3 = pid.pi_nod()[tdata.node3 as usize].indices().flag;
                    let mut d1 = 0.0;
                    let mut d2 = 0.0;
                    let mut d3 = 0.0;

                    // `cxx:3227-3247`: node 1, opposite edge `node2 -> node3`.
                    if (f1 & nmask::EDGE == 0 || i_loop > 1)
                        && (f1 & nmask::OUT_L == 0 || i_loop > 1)
                        && (f1 & nmask::VERT == 0)
                    {
                        let v = nudge(&n2.point, &n3.point, &n1.point);
                        d1 = v.0;
                        x1 = v.1;
                        y1 = v.2;
                    }
                    // `cxx:3249-3271`: node 2, opposite edge `node3 -> node1`.
                    if (f2 & nmask::EDGE == 0 || i_loop > 1)
                        && (f2 & nmask::OUT_L == 0 || i_loop > 1)
                        && (f2 & nmask::VERT == 0)
                    {
                        let v = nudge(&n3.point, &n1.point, &n2.point);
                        d2 = v.0;
                        x2 = v.1;
                        y2 = v.2;
                    }
                    // `cxx:3272-3295`: node 3, opposite edge `node1 -> node2`.
                    if (f3 & nmask::EDGE == 0 || i_loop > 1)
                        && (f3 & nmask::OUT_L == 0 || i_loop > 1)
                        && (f3 & nmask::VERT == 0)
                    {
                        let v = nudge(&n1.point, &n2.point, &n3.point);
                        d3 = v.0;
                        x3 = v.1;
                        y3 = v.2;
                    }

                    let moved = if d1 > d2 && d1 > d3 {
                        Some((tdata.node1, x1, y1))
                    } else if d2 > d3 && d2 > d1 {
                        Some((tdata.node2, x2, y2))
                    } else if d3 > d1 && d3 > d2 {
                        Some((tdata.node3, x3, y3))
                    } else {
                        None
                    };
                    if let Some((node, dx, dy)) = moved {
                        let p = pid.pi_nod()[node as usize].data().point;
                        pid.change_pi_nod()[node as usize].change_data().point =
                            GpXyz::new(p.x() + dx, p.y() + dy, p.z());
                        pid.change_pi_nod()[node as usize].change_indices().flag |= nmask::MOVE;
                        self.update_around_node(pid, node);
                        is_fr_back_in_list = true;
                    }
                }
            }

            if is_fr_back_in_list {
                for bp in the_list.iter_mut() {
                    let ind: BiPointIndices = *bp.indices();
                    // `cxx:3380-3406`: face 1's moved nodes, propagated into
                    // face 2's twin node when the segment is shared.
                    if ind.face_conex1 != 0 {
                        let moved1 = the_pid[ind.face_conex1 as usize]
                            .as_ref()
                            .map(|p| {
                                (
                                    *p.pi_nod()[ind.face1_pt1 as usize].data(),
                                    p.pi_nod()[ind.face1_pt1 as usize].indices().flag,
                                )
                            });
                        if let Some((r, flag)) = moved1 {
                            if flag & nmask::MOVE != 0 {
                                let p = &mut bp.change_points();
                                p.pnt1 = r.point;
                                p.pntp1 = r.point;
                                self.timultiply_xyz(&mut p.pnt1, false);
                                if ind.face_conex2 != 0 {
                                    if let Some(pid2) =
                                        the_pid[ind.face_conex2 as usize].as_mut()
                                    {
                                        let q = pid2.pi_nod()[ind.face2_pt1 as usize]
                                            .data()
                                            .point;
                                        pid2.change_pi_nod()[ind.face2_pt1 as usize]
                                            .change_data()
                                            .point = GpXyz::new(r.point.x(), r.point.y(), q.z());
                                        self.update_around_node(pid2, ind.face2_pt1);
                                    }
                                }
                            }
                        }
                        let moved2 = the_pid[ind.face_conex1 as usize]
                            .as_ref()
                            .map(|p| {
                                (
                                    *p.pi_nod()[ind.face1_pt2 as usize].data(),
                                    p.pi_nod()[ind.face1_pt2 as usize].indices().flag,
                                )
                            });
                        if let Some((r, flag)) = moved2 {
                            if flag & nmask::MOVE != 0 {
                                let p = &mut bp.change_points();
                                p.pnt2 = r.point;
                                p.pntp2 = r.point;
                                self.timultiply_xyz(&mut p.pnt2, false);
                                if ind.face_conex2 != 0 {
                                    if let Some(pid2) =
                                        the_pid[ind.face_conex2 as usize].as_mut()
                                    {
                                        let q = pid2.pi_nod()[ind.face2_pt2 as usize]
                                            .data()
                                            .point;
                                        pid2.change_pi_nod()[ind.face2_pt2 as usize]
                                            .change_data()
                                            .point = GpXyz::new(r.point.x(), r.point.y(), q.z());
                                        self.update_around_node(pid2, ind.face2_pt2);
                                    }
                                }
                            }
                        }
                    }
                    // `cxx:3435-3489`: face 2's moved nodes.
                    if ind.face_conex2 != 0 {
                        let moved1 = the_pid[ind.face_conex2 as usize]
                            .as_ref()
                            .map(|p| {
                                (
                                    *p.pi_nod()[ind.face2_pt1 as usize].data(),
                                    p.pi_nod()[ind.face2_pt1 as usize].indices().flag,
                                )
                            });
                        if let Some((r, flag)) = moved1 {
                            if flag & nmask::MOVE != 0 {
                                let p = &mut bp.change_points();
                                p.pnt1 = r.point;
                                p.pntp1 = r.point;
                                self.timultiply_xyz(&mut p.pnt1, false);
                                if ind.face_conex1 != 0 {
                                    if let Some(pid1) =
                                        the_pid[ind.face_conex1 as usize].as_mut()
                                    {
                                        let q = pid1.pi_nod()[ind.face1_pt1 as usize]
                                            .data()
                                            .point;
                                        pid1.change_pi_nod()[ind.face1_pt1 as usize]
                                            .change_data()
                                            .point = GpXyz::new(r.point.x(), r.point.y(), q.z());
                                        self.update_around_node(pid1, ind.face1_pt1);
                                    }
                                }
                            }
                        }
                        let moved2 = the_pid[ind.face_conex2 as usize]
                            .as_ref()
                            .map(|p| {
                                (
                                    *p.pi_nod()[ind.face2_pt2 as usize].data(),
                                    p.pi_nod()[ind.face2_pt2 as usize].indices().flag,
                                )
                            });
                        if let Some((r, flag)) = moved2 {
                            if flag & nmask::MOVE != 0 {
                                let p = &mut bp.change_points();
                                p.pnt2 = r.point;
                                p.pntp2 = r.point;
                                self.timultiply_xyz(&mut p.pnt2, false);
                                if ind.face_conex1 != 0 {
                                    if let Some(pid1) =
                                        the_pid[ind.face_conex1 as usize].as_mut()
                                    {
                                        let q = pid1.pi_nod()[ind.face1_pt2 as usize]
                                            .data()
                                            .point;
                                        pid1.change_pi_nod()[ind.face1_pt2 as usize]
                                            .change_data()
                                            .point = GpXyz::new(r.point.x(), r.point.y(), q.z());
                                        self.update_around_node(pid1, ind.face1_pt2);
                                    }
                                }
                            }
                        }
                    }
                }

                // `cxx:3492-3507`: clear `NMsk_Move` on every node.
                for a_face_iter in 1..=a_nb_faces {
                    if let Some(pid) = the_pid[a_face_iter].as_mut() {
                        let a_nb_nodes = pid.nb_pi_nod();
                        for a_node_iter in 1..=a_nb_nodes {
                            pid.change_pi_nod()[a_node_iter as usize].change_indices().flag &=
                                !nmask::MOVE;
                        }
                    }
                }
            }
        }
    }

    /// `UpdateOutLines(theList, thePID)` (`cxx:3998-4153`): for every segment
    /// of a face that carries an internal outline, mark the outline edge on
    /// both adjacent triangles and append the projected two-face bi-point.
    pub fn update_out_lines(
        &mut self,
        the_list: &mut Vec<HlrAlgoBiPoint>,
        the_pid: &mut PolyInternalDataArray,
    ) {
        let a_nb_faces = self.f_map.extent();
        for a_face_iter in 1..=a_nb_faces {
            if the_pid[a_face_iter].is_none() {
                continue;
            }
            let pid = the_pid[a_face_iter].as_mut().unwrap();
            if !pid.int_out_l() {
                continue;
            }
            let a_nb_segs = pid.nb_pi_seg();
            for a_seg_iter in 1..=a_nb_segs {
                let psg = pid.pi_seg()[a_seg_iter as usize];
                let it1 = psg.conex1;
                let it2 = psg.conex2;
                if it1 == 0 || it2 == 0 || it1 == it2 {
                    continue;
                }
                let t1 = pid.t_data()[it1 as usize];
                let t2 = pid.t_data()[it2 as usize];
                if !is_out_line_edge(&t1, &t2) {
                    continue;
                }
                let pd = psg.lst_sg1;
                let pf = psg.lst_sg2;
                mark_out_line_edge(&mut pid.change_t_data()[it1 as usize], pd, pf);
                mark_out_line_edge(&mut pid.change_t_data()[it2 as usize], pd, pf);

                let n1 = pid.pi_nod()[pd as usize].data().point;
                let n2 = pid.pi_nod()[pf as usize].data().point;
                let (x1, y1, z1) = (n1.x(), n1.y(), n1.z());
                let (x2, y2, z2) = (n2.x(), n2.y(), n2.z());
                let (mut xt1, mut yt1, mut zt1) = (x1, y1, z1);
                let (mut xt2, mut yt2, mut zt2) = (x2, y2, z2);
                self.timultiply(&mut xt1, &mut yt1, &mut zt1, false);
                self.timultiply(&mut xt2, &mut yt2, &mut zt2, false);
                let f = a_face_iter as i32;
                the_list.push(HlrAlgoBiPoint::from_two_faces_with_flag(
                    xt1, yt1, zt1, xt2, yt2, zt2, x1, y1, z1, x2, y2, z2, f, f, pd, pf, f, pd,
                    pf, 12,
                ));
            }
        }
    }
}

/// `cxx:3228-3247` (and its two rotations): the outward nudge of the node
/// opposite the edge `p_from -> p_to`. Returns `(edge_length, dx, dy)`; the
/// `dd` sign test adds the `1%` push OCCT uses to keep the node off the edge.
fn nudge(p_from: &GpXyz, p_to: &GpXyz, p_node: &GpXyz) -> (f64, f64, f64) {
    let dx = p_to.x() - p_from.x();
    let dy = p_to.y() - p_from.y();
    let d = (dx * dx + dy * dy).sqrt();
    let nx = -dy / d;
    let ny = dx / d;
    let ax = p_node.x() - p_from.x();
    let ay = p_node.y() - p_from.y();
    let mut dd = -(ax * nx + ay * ny);
    if dd < 0.0 {
        dd -= d * 0.01;
    } else {
        dd += d * 0.01;
    }
    (d, nx * dd, ny * dd)
}
