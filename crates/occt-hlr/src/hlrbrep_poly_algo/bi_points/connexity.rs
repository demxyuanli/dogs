//! `InitBiPointsWithConnexity` and its `aNbConnex == 2` arm.
//! Source: `HLRBRep_PolyAlgo.cxx:878-1377`.

use super::*;

impl HlrBRepPolyAlgo {
    /// `InitBiPointsWithConnexity(theIEdge, theEdge, theList, thePID, theLS,
    /// theIsConnex)` (`cxx:878-1377`).
    ///
    /// `the_ls` is the `theLS` face list of `TopExp::MapShapesAndAncestors`
    /// (`cxx:582`): the faces incident to the edge, empty when the edge is
    /// isolated (`cxx:588-590`). `the_i_edge` is `theIEdge`, the index the
    /// bi-points carry as `ShapeIndex`. `the_rg` is
    /// `BRep_Tool::Continuity(theEdge, aF1, aF2)` (`cxx:1080`), the
    /// `GeomAbs_Shape` value the two-face arm turns into the bi-point flag
    /// (`cxx:1512-1519`).
    pub fn init_bi_points_with_connexity(
        &mut self,
        the_i_edge: usize,
        the_edge: &TopoShape,
        the_list: &mut Vec<HlrAlgoBiPoint>,
        the_pid: &mut PolyInternalDataArray,
        the_ls: &[TopoShape],
        the_is_connex: bool,
        the_rg: u8,
    ) {
        let edge = Edge(the_edge.clone());
        // `myBCurv.Initialize(theEdge)` (`cxx:891`).
        let bcurv = occt_topo::brep_tool::BRepTool::edge_curve(&edge);
        let edge_reversed = the_edge.orientation().is_reversed();
        let mut seg = Segment::default();

        if the_is_connex {
            let nb_connex = the_ls.len();
            if nb_connex == 1 {
                let f1 = &the_ls[0];
                let i1 = self.f_map.find_index(f1);
                if i1 == 0 {
                    // `thePID.Value(i1)` with `i1 == 0` would be out of range in
                    // OCCT (`cxx:903`); a face that was never mapped cannot be
                    // one of `theLS`'s entries, so this is unreachable.
                    return;
                }
                let Some(tr1) = face_triangulation(f1) else {
                    return;
                };
                let Some(hpol) = occt_topo::brep_tool::BRepTool::polygon_on_triangulation(
                    &edge,
                    &tr1,
                    edge_reversed,
                ) else {
                    return;
                };
                let Some(par) = hpol.parameters.as_ref() else {
                    return;
                };
                if i1 >= the_pid.len() || the_pid[i1].is_none() {
                    return;
                }
                let a_pol1: Vec<i32> = hpol.nodes.iter().map(|&n| n as i32).collect();
                let a_nb_pol = a_pol1.len();
                if a_nb_pol == 0 {
                    return;
                }
                // `myPC.Initialize(theEdge, aF1)` (`cxx:905`).
                let pc = pcurve_of(the_edge, f1, edge_reversed);

                {
                    let data = the_pid[i1].as_mut().unwrap();
                    // `cxx:917-919`: only the two ends get `NMsk_Vert`; for a
                    // closed polygon both writes hit the same node.
                    data.change_pi_nod()[a_pol1[0] as usize]
                        .change_indices()
                        .flag |= nmask::VERT;
                    data.change_pi_nod()[a_pol1[a_nb_pol - 1] as usize]
                        .change_indices()
                        .flag |= nmask::VERT;
                    // `cxx:921-935`.
                    for i_pol in 0..a_nb_pol {
                        let node = a_pol1[i_pol] as usize;
                        let par_i = par[i_pol];
                        let (edg1, _) = {
                            let ind = data.pi_nod()[node].indices();
                            (ind.edg1, ind.edg2)
                        };
                        if edg1 == 0 || edg1 == the_i_edge as i32 {
                            data.change_pi_nod()[node].change_indices().edg1 = the_i_edge as i32;
                            data.change_pi_nod()[node].change_data().pcu1 = par_i;
                        } else {
                            data.change_pi_nod()[node].change_indices().edg2 = the_i_edge as i32;
                            data.change_pi_nod()[node].change_data().pcu2 = par_i;
                        }
                    }
                }

                let mut i1p2 = a_pol1[0];
                // `cxx:938-940`: `aNode12Indices = aNode11Indices` and the same
                // for the values, i.e. node 12 is node 11 = `aPol1(1)`.
                let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
                let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
                seg.x2 = n12.point.x();
                seg.y2 = n12.point.y();
                seg.z2 = n12.point.z();
                seg.xt2 = seg.x2;
                seg.yt2 = seg.y2;
                seg.zt2 = seg.z2;
                if ind12.edg1 == the_i_edge as i32 {
                    seg.u2 = n12.pcu1;
                } else if ind12.edg2 == the_i_edge as i32 {
                    seg.u2 = n12.pcu2;
                }
                // `cxx:959-1027`.
                the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
                    .change_indices()
                    .flag |= nmask::EDGE;
                self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                if a_pol1[0] == a_pol1[a_nb_pol - 1] {
                    if let Some(pc) = &pc {
                        if pc.is_periodic() {
                            // `cxx:961-964`.
                            seg.u2 -= pc.period();
                        }
                    }
                }

                if a_nb_pol == 2 && occt_topo::brep_tool::BRepTool::is_degenerated(&edge) {
                    // `cxx:966-971`.
                    let data = the_pid[i1].as_mut().unwrap();
                    Self::check_degenerated_segment(
                        data,
                        self.tol_angular,
                        a_pol1[0],
                        a_pol1[a_nb_pol - 1],
                    );
                    self.update_around_node(data, a_pol1[0]);
                    self.update_around_node(data, a_pol1[a_nb_pol - 1]);
                } else {
                    for i_pol in 1..a_nb_pol {
                        let i1p1 = i1p2;
                        i1p2 = a_pol1[i_pol];
                        // `cxx:999-1008`.
                        seg.x1 = seg.x2;
                        seg.y1 = seg.y2;
                        seg.z1 = seg.z2;
                        seg.xt1 = seg.xt2;
                        seg.yt1 = seg.yt2;
                        seg.zt1 = seg.zt2;
                        seg.u1 = seg.u2;
                        let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
                        let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
                        seg.x2 = n12.point.x();
                        seg.y2 = n12.point.y();
                        seg.z2 = n12.point.z();
                        seg.xt2 = seg.x2;
                        seg.yt2 = seg.y2;
                        seg.zt2 = seg.z2;
                        if ind12.edg1 == the_i_edge as i32 {
                            seg.u2 = n12.pcu1;
                        } else if ind12.edg2 == the_i_edge as i32 {
                            seg.u2 = n12.pcu2;
                        }
                        the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
                            .change_indices()
                            .flag |= nmask::EDGE;
                        self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                        let side = FaceSide {
                            face: i1,
                            p1: i1p1,
                            p2: i1p2,
                        };
                        self.interpolation_one(&bcurv, the_list, &mut seg, the_i_edge, side, the_pid);
                    }
                }
            } else if nb_connex == 2 {
                self.init_bi_points_two_faces(
                    the_i_edge,
                    &edge,
                    the_edge,
                    the_list,
                    the_pid,
                    the_ls,
                    &bcurv,
                    &mut seg,
                    the_rg,
                );
            }
        } else {
            // `cxx:1311-1374`: no connexity, so the edge contributes its
            // `Poly_Polygon3D` segments directly.
            let Some(a_polyg) = occt_topo::brep_tool::BRepTool::edge_polygon3d(&edge) else {
                return;
            };
            // `cxx:1315-1335`: `aTT = aLoc.Transformation();
            // aTT.PreMultiply(aProjTrsf)` read into `TTMa`/`TTLo`. The port's
            // edge location is the identity, so this is the projection alone.
            self.read_face_matrices(&occt_core::gp::GpTrsf::identity());
            let a_pol = &a_polyg.nodes;
            if a_pol.is_empty() {
                return;
            }
            let a_nb_pol1 = a_pol.len();
            let p1 = a_pol[0];
            seg.x2 = p1.x();
            seg.y2 = p1.y();
            seg.z2 = p1.z();
            self.ttmultiply(&mut seg.x2, &mut seg.y2, &mut seg.z2, false);
            seg.xt2 = seg.x2;
            seg.yt2 = seg.y2;
            seg.zt2 = seg.z2;
            self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
            for j_pol in 1..a_nb_pol1 {
                seg.x1 = seg.x2;
                seg.y1 = seg.y2;
                seg.z1 = seg.z2;
                seg.xt1 = seg.xt2;
                seg.yt1 = seg.yt2;
                seg.zt1 = seg.zt2;
                let p2 = a_pol[j_pol];
                seg.x2 = p2.x();
                seg.y2 = p2.y();
                seg.z2 = p2.z();
                self.ttmultiply(&mut seg.x2, &mut seg.y2, &mut seg.z2, false);
                seg.xt2 = seg.x2;
                seg.yt2 = seg.y2;
                seg.zt2 = seg.z2;
                self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                // `cxx:1363-1364`: `Pnt*` are the 3D points (`XTI*`), `PntP*`
                // the projected ones (`X*`).
                the_list.insert(
                    0,
                    HlrAlgoBiPoint::with_flag(
                        seg.xt1, seg.yt1, seg.zt1, seg.xt2, seg.yt2, seg.zt2, seg.x1, seg.y1,
                        seg.z1, seg.x2, seg.y2, seg.z2, the_i_edge as i32, 0,
                    ),
                );
            }
        }
    }

    /// The `aNbConnex == 2` arm of `InitBiPointsWithConnexity`
    /// (`cxx:1064-1309`).
    #[allow(clippy::too_many_arguments)]
    fn init_bi_points_two_faces(
        &mut self,
        the_i_edge: usize,
        edge: &Edge,
        the_edge: &TopoShape,
        the_list: &mut Vec<HlrAlgoBiPoint>,
        the_pid: &mut PolyInternalDataArray,
        the_ls: &[TopoShape],
        bcurv: &Option<Arc<dyn Curve>>,
        seg: &mut Segment,
        the_rg: u8,
    ) {
        let f1 = &the_ls[0];
        let i1 = self.f_map.find_index(f1);
        let orig_reversed = the_edge.orientation().is_reversed();
        let f2 = &the_ls[1];
        let i2 = self.f_map.find_index(f2);
        // `cxx:1075-1078`: the same face twice means the edge is its own seam:
        // OCCT reverses the local edge so the second lookup answers the other
        // polygon of the closed triangulation.
        let reversed2 = if i1 == i2 {
            !orig_reversed
        } else {
            orig_reversed
        };
        if i1 == 0 || i2 == 0 {
            return;
        }
        let Some(tr1) = face_triangulation(f1) else {
            return;
        };
        let Some(tr2) = face_triangulation(f2) else {
            return;
        };
        let pol1 = occt_topo::brep_tool::BRepTool::polygon_on_triangulation(edge, &tr1, orig_reversed);
        let pol2 = occt_topo::brep_tool::BRepTool::polygon_on_triangulation(edge, &tr2, reversed2);
        let (Some(hpol1), Some(hpol2)) = (pol1, pol2) else {
            return;
        };
        let Some(par) = hpol1.parameters.as_ref() else {
            return;
        };
        if i1 >= the_pid.len()
            || i2 >= the_pid.len()
            || the_pid[i1].is_none()
            || the_pid[i2].is_none()
        {
            return;
        }
        let a_pol1: Vec<i32> = hpol1.nodes.iter().map(|&n| n as i32).collect();
        let a_pol2: Vec<i32> = hpol2.nodes.iter().map(|&n| n as i32).collect();
        let a_nb_pol1 = a_pol1.len();
        if a_nb_pol1 == 0 || a_pol2.len() < a_nb_pol1 {
            return;
        }
        // `cxx:1085`: `myPC.Initialize(theEdge, aF1)` — after the reverse, so
        // the reversed pcurve of a seam edge is the one that is consulted.
        let pc = pcurve_of(the_edge, f1, orig_reversed);

        // `cxx:1109-1112`: four end nodes get `NMsk_Vert`.
        for (face, pol) in [(i1, &a_pol1), (i2, &a_pol2)] {
            let data = the_pid[face].as_mut().unwrap();
            data.change_pi_nod()[pol[0] as usize]
                .change_indices()
                .flag |= nmask::VERT;
            data.change_pi_nod()[pol[a_nb_pol1 - 1] as usize]
                .change_indices()
                .flag |= nmask::VERT;
        }
        // `cxx:1114-1143`: both polygons' nodes take the edge parameter.
        for i_pol in 0..a_nb_pol1 {
            let par_i = par[i_pol];
            for (face, pol) in [(i1, &a_pol1), (i2, &a_pol2)] {
                let data = the_pid[face].as_mut().unwrap();
                let node = pol[i_pol] as usize;
                let edg1 = data.pi_nod()[node].indices().edg1;
                if edg1 == 0 || edg1 == the_i_edge as i32 {
                    data.change_pi_nod()[node].change_indices().edg1 = the_i_edge as i32;
                    data.change_pi_nod()[node].change_data().pcu1 = par_i;
                } else {
                    data.change_pi_nod()[node].change_indices().edg2 = the_i_edge as i32;
                    data.change_pi_nod()[node].change_data().pcu2 = par_i;
                }
            }
        }

        let mut i1p2 = a_pol1[0];
        let mut i2p2 = a_pol2[0];
        // `cxx:1151-1159`: node 12 is node 11, node 22 is node 21.
        let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
        let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
        seg.x2 = n12.point.x();
        seg.y2 = n12.point.y();
        seg.z2 = n12.point.z();
        seg.xt2 = seg.x2;
        seg.yt2 = seg.y2;
        seg.zt2 = seg.z2;
        if ind12.edg1 == the_i_edge as i32 {
            seg.u2 = n12.pcu1;
        } else if ind12.edg2 == the_i_edge as i32 {
            seg.u2 = n12.pcu2;
        }
        // `cxx:1169-1171`.
        the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
            .change_indices()
            .flag |= nmask::EDGE;
        the_pid[i2].as_mut().unwrap().change_pi_nod()[i2p2 as usize]
            .change_indices()
            .flag |= nmask::EDGE;
        self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
        if a_pol1[0] == a_pol1[a_nb_pol1 - 1] {
            if let Some(pc) = &pc {
                if pc.is_periodic() {
                    // `cxx:1172-1175`.
                    seg.u2 -= pc.period();
                }
            }
        }

        if a_nb_pol1 == 2 && occt_topo::brep_tool::BRepTool::is_degenerated(edge) {
            // `cxx:1177-1186`.
            {
                let data = the_pid[i1].as_mut().unwrap();
                Self::check_degenerated_segment(data, self.tol_angular, a_pol1[0], a_pol1[1]);
                self.update_around_node(data, a_pol1[0]);
                self.update_around_node(data, a_pol1[1]);
            }
            {
                let data = the_pid[i2].as_mut().unwrap();
                Self::check_degenerated_segment(data, self.tol_angular, a_pol2[0], a_pol2[1]);
                self.update_around_node(data, a_pol2[0]);
                self.update_around_node(data, a_pol2[1]);
            }
        } else {
            for i_pol in 1..a_nb_pol1 {
                let i1p1 = i1p2;
                let i2p1 = i2p2;
                i1p2 = a_pol1[i_pol];
                i2p2 = a_pol2[i_pol];
                // `cxx:1231-1239`.
                seg.x1 = seg.x2;
                seg.y1 = seg.y2;
                seg.z1 = seg.z2;
                seg.xt1 = seg.xt2;
                seg.yt1 = seg.yt2;
                seg.zt1 = seg.zt2;
                seg.u1 = seg.u2;
                let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
                let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
                seg.x2 = n12.point.x();
                seg.y2 = n12.point.y();
                seg.z2 = n12.point.z();
                seg.xt2 = seg.x2;
                seg.yt2 = seg.y2;
                seg.zt2 = seg.z2;
                if ind12.edg1 == the_i_edge as i32 {
                    seg.u2 = n12.pcu1;
                } else if ind12.edg2 == the_i_edge as i32 {
                    seg.u2 = n12.pcu2;
                }
                // `cxx:1256-1257`: the flag is OR-ed twice in OCCT (a duplicated
                // line); one OR is the same state.
                the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
                    .change_indices()
                    .flag |= nmask::EDGE;
                self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                let side1 = FaceSide {
                    face: i1,
                    p1: i1p1,
                    p2: i1p2,
                };
                let side2 = FaceSide {
                    face: i2,
                    p1: i2p1,
                    p2: i2p2,
                };
                self.interpolation_two(
                    bcurv, the_list, seg, the_i_edge, side1, side2, the_pid, the_rg,
                );
            }
        }
    }
}
