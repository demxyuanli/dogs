//! Two-face `Interpolation` / `MoveOrInsertPoint` overloads.
//! Source: `HLRBRep_PolyAlgo.cxx:1466-1853`, `:2116-2449`, `:2452-2978`.

use super::*;

impl HlrBRepPolyAlgo {
    /// `Interpolation(..., theRg, Nod11, Nod12, ..., Nod21, Nod22, ...)`
    /// (`cxx:1466-1853`): the two-face overload. One interpolated point per
    /// face, then the crossing nearest `P1` (or both) is inserted.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn interpolation_two(
        &mut self,
        bcurv: &Option<Arc<dyn Curve>>,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: usize,
        f1: FaceSide,
        f2: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        the_rg: u8,
    ) {
        // `cxx:1512-1519`: `flag += 1` for G1, `+= 2` for G2. The port's
        // `GeomAbs_Shape` values are these integers (`GeomAbs_Shape.hxx:26-29`).
        let mut flag = 0;
        if the_rg >= GEOMABS_G1 {
            flag += 1;
        }
        if the_rg >= GEOMABS_G2 {
            flag += 2;
        }
        let n11 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p1);
        let n12 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p2);
        let p3 = self.interpolation_point(bcurv, seg.u1, seg.u2, &n11, &n12);
        let n21 = node_data(the_pid[f2.face].as_ref().unwrap(), f2.p1);
        let n22 = node_data(the_pid[f2.face].as_ref().unwrap(), f2.p2);
        let p4 = self.interpolation_point(bcurv, seg.u1, seg.u2, &n21, &n22);
        if p3.ins || p4.ins {
            if !p4.ins {
                // p1 i1p3 p2
                self.move_or_insert_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f1,
                    f2,
                    the_pid,
                    &p3,
                    flag,
                );
            } else if !p3.ins {
                // p1 i2p4 p2 — the sides swap, so `f2` plays the "first" role.
                self.move_or_insert_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f2,
                    f1,
                    the_pid,
                    &p4,
                    flag,
                );
            } else if (p4.coef - p3.coef).abs() < self.tol_sta {
                // `cxx:1654`: the two crossings coincide -> keep the face-2 one.
                self.move_or_insert_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f2,
                    f1,
                    the_pid,
                    &p4,
                    flag,
                );
            } else if p4.coef < p3.coef {
                // p1 i2p4 i1p3 p2
                self.move_or_insert_two_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f2,
                    f1,
                    the_pid,
                    &p4,
                    &p3,
                    flag,
                );
            } else {
                // p1 i1p3 i2p4 p2
                self.move_or_insert_two_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f1,
                    f2,
                    the_pid,
                    &p3,
                    &p4,
                    flag,
                );
            }
        } else {
            // `cxx:1830-1851`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge as i32,
                    f1.face as i32,
                    f1.p1,
                    f1.p2,
                    f2.face as i32,
                    f2.p1,
                    f2.p2,
                    flag,
                ),
            );
        }
    }

    /// The two-face, one-point `MoveOrInsertPoint` (`cxx:2116-2449`).
    ///
    /// `f1`/`f2` are `(theNod1*, theI1*)` and `(theNod2*, theI2*)`; the point is
    /// `theP3` with its bundle.
    #[allow(clippy::too_many_arguments)]
    fn move_or_insert_two(
        &mut self,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: i32,
        f1: FaceSide,
        f2: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        p3: &Interp,
        the_flag: i32,
    ) {
        let mut an_ins3 = p3.ins;
        if an_ins3 && p3.m_p3_p1 {
            // `cxx:2168-2249`: P1 ---> P3 on both faces.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p1).flag & nmask::VERT != 0;
            if !is_vert && p3.coef < self.tol_sta {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f1.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p1);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f2.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f2.p1);
                }
                seg.x1 = p3.x;
                seg.y1 = p3.y;
                seg.z1 = p3.z;
                seg.xt1 = p3.xt;
                seg.yt1 = p3.yt;
                seg.zt1 = p3.zt;
                seg.u1 = p3.u;
                if let Some(front) = list.first_mut() {
                    let pts = front.change_points();
                    pts.pntp2 = GpXyz::new(p3.x, p3.y, p3.z);
                    pts.pnt2 = GpXyz::new(p3.xt, p3.yt, p3.zt);
                }
            }
        }
        if an_ins3 && !p3.m_p3_p1 {
            // `cxx:2251-2329`: P2 ---> P3 on both faces.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p2).flag & nmask::VERT != 0;
            if !is_vert && p3.coef > self.tol_end {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, false);
                    move_node_onto_point(data, f1.p2, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p2);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p3.coef, p3.x, p3.y, p3.z, false);
                    move_node_onto_point(data, f2.p2, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f2.p2);
                }
                seg.x2 = p3.x;
                seg.y2 = p3.y;
                seg.z2 = p3.z;
                seg.xt2 = p3.xt;
                seg.yt2 = p3.yt;
                seg.zt2 = p3.zt;
                seg.u2 = p3.u;
            }
        }
        if an_ins3 {
            // `cxx:2331-2423`: p1 i1p3 p2 with a node on each face.
            let (an_i1p3, an_i2p3) = {
                let data1 = the_pid[f1.face].as_mut().unwrap();
                let n11 = node_data(data1, f1.p1);
                let n12 = node_data(data1, f1.p2);
                let new1 = data1.add_node(&n11, &n12, p3.coef, p3.x, p3.y, p3.z);
                data1.change_pi_nod()[new1 as usize].change_indices().edg1 = the_i_edge;
                data1.change_pi_nod()[new1 as usize].change_indices().flag |=
                    nmask::OUT_L | nmask::EDGE;
                data1.change_pi_nod()[new1 as usize].change_data().pcu1 = p3.u;
                data1.change_pi_nod()[new1 as usize].change_data().scal = 0.0;
                let data2 = the_pid[f2.face].as_mut().unwrap();
                let n21 = node_data(data2, f2.p1);
                let n22 = node_data(data2, f2.p2);
                let new2 = data2.add_node(&n21, &n22, p3.coef, p3.x, p3.y, p3.z);
                data2.change_pi_nod()[new2 as usize].change_indices().edg1 = the_i_edge;
                data2.change_pi_nod()[new2 as usize].change_indices().flag |=
                    nmask::OUT_L | nmask::EDGE;
                data2.change_pi_nod()[new2 as usize].change_data().pcu1 = p3.u;
                data2.change_pi_nod()[new2 as usize].change_data().scal = 0.0;
                (new1, new2)
            };
            {
                let data = the_pid[f1.face].as_mut().unwrap();
                data.update_links_split(f1.p1, f1.p2, an_i1p3);
                self.update_around_node(data, an_i1p3);
            }
            {
                let data = the_pid[f2.face].as_mut().unwrap();
                data.update_links_split(f2.p1, f2.p2, an_i2p3);
                self.update_around_node(data, an_i2p3);
            }
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    p3.xt,
                    p3.yt,
                    p3.zt,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    p3.x,
                    p3.y,
                    p3.z,
                    the_i_edge,
                    f1.face as i32,
                    f1.p1,
                    an_i1p3,
                    f2.face as i32,
                    f2.p1,
                    an_i2p3,
                    the_flag,
                ),
            );
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    p3.xt,
                    p3.yt,
                    p3.zt,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    p3.x,
                    p3.y,
                    p3.z,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge,
                    f1.face as i32,
                    an_i1p3,
                    f1.p2,
                    f2.face as i32,
                    an_i2p3,
                    f2.p2,
                    the_flag,
                ),
            );
        } else {
            // `cxx:2426-2448`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge,
                    f1.face as i32,
                    f1.p1,
                    f1.p2,
                    f2.face as i32,
                    f2.p1,
                    f2.p2,
                    the_flag,
                ),
            );
        }
    }

    /// The two-face, two-point `MoveOrInsertPoint` (`cxx:2453-2976`).
    ///
    /// `f1`/`f2` pair with `p3`/`p4`; the snapping arms touch `f1`+`f2` for
    /// `P1 -> P3`, `f1`+`f2` for `P2 -> P4`, and then fall back to the
    /// one-point overloads when only one of the two points survived.
    #[allow(clippy::too_many_arguments)]
    fn move_or_insert_two_two(
        &mut self,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: i32,
        f1: FaceSide,
        f2: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        p3: &Interp,
        p4: &Interp,
        the_flag: i32,
    ) {
        let mut an_ins3 = p3.ins;
        let mut an_ins4 = p4.ins;
        if an_ins3 && p3.m_p3_p1 {
            // `cxx:2516-2597`: P1 ---> P3.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p1).flag & nmask::VERT != 0;
            if !is_vert && p3.coef < self.tol_sta {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f1.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p1);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f2.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f2.p1);
                }
                seg.x1 = p3.x;
                seg.y1 = p3.y;
                seg.z1 = p3.z;
                seg.xt1 = p3.xt;
                seg.yt1 = p3.yt;
                seg.zt1 = p3.zt;
                seg.u1 = p3.u;
                if let Some(front) = list.first_mut() {
                    let pts = front.change_points();
                    pts.pntp2 = GpXyz::new(p3.x, p3.y, p3.z);
                    pts.pnt2 = GpXyz::new(p3.xt, p3.yt, p3.zt);
                }
            }
        }
        if an_ins4 && !p4.m_p3_p1 {
            // `cxx:2599-2678`: P2 ---> P4.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p2).flag & nmask::VERT != 0;
            if !is_vert && p4.coef > self.tol_end {
                an_ins4 = false;
                // `cxx:2603-2630`: `ChangeNode` on the face-2 side first, then on
                // the face-1 side; each already re-orients its own node.
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p4.coef, p4.x, p4.y, p4.z, false);
                }
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p4.coef, p4.x, p4.y, p4.z, false);
                }
                // `cxx:2639-2665`: the parameter writes of the two sides.
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    move_node_onto_point(data, f1.p2, the_i_edge, p4.x, p4.y, p4.z, p4.u);
                    self.update_around_node(data, f1.p2);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    move_node_onto_point(data, f2.p2, the_i_edge, p4.x, p4.y, p4.z, p4.u);
                    self.update_around_node(data, f2.p2);
                }
                seg.x2 = p4.x;
                seg.y2 = p4.y;
                seg.z2 = p4.z;
                seg.xt2 = p4.xt;
                seg.yt2 = p4.yt;
                seg.zt2 = p4.zt;
                seg.u2 = p4.u;
            }
        }
        if an_ins3 || an_ins4 {
            if !an_ins4 {
                // `cxx:2681-2733`: p1 i1p3 p2.
                self.move_or_insert_two(
                    list, seg, the_i_edge, f1, f2, the_pid, p3, the_flag,
                );
            } else if !an_ins3 {
                // `cxx:2734-2786`: p1 i2p4 p2.
                self.move_or_insert_two(
                    list, seg, the_i_edge, f2, f1, the_pid, p4, the_flag,
                );
            } else {
                // `cxx:2787-2952`: p1 i1p3 i2p4 p2, a node on each face for each
                // of the two points. `AddNode` runs face-1/face-2/face-1/face-2 in
                // that order (`cxx:2787-2815`); the two nodes of one face keep
                // their relative numbering, which is all that matters.
                let an_i1p3 = {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    let n11 = node_data(data, f1.p1);
                    let n12 = node_data(data, f1.p2);
                    data.add_node(&n11, &n12, p3.coef, p3.x, p3.y, p3.z)
                };
                let an_i2p3 = {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    let n21 = node_data(data, f2.p1);
                    let n22 = node_data(data, f2.p2);
                    data.add_node(&n21, &n22, p3.coef, p3.x, p3.y, p3.z)
                };
                let an_i1p4 = {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    let n11 = node_data(data, f1.p1);
                    let n12 = node_data(data, f1.p2);
                    data.add_node(&n11, &n12, p4.coef, p4.x, p4.y, p4.z)
                };
                let an_i2p4 = {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    let n21 = node_data(data, f2.p1);
                    let n22 = node_data(data, f2.p2);
                    data.add_node(&n21, &n22, p4.coef, p4.x, p4.y, p4.z)
                };
                // `cxx:2831-2851`: `Edg1`/`PCu1`/`Scal = 0`/`OutL`/`Edge` per new
                // node, face 1 first.
                for (face, node, u) in [
                    (f1.face, an_i1p3, p3.u),
                    (f2.face, an_i2p3, p3.u),
                    (f1.face, an_i1p4, p4.u),
                    (f2.face, an_i2p4, p4.u),
                ] {
                    let data = the_pid[face].as_mut().unwrap();
                    data.change_pi_nod()[node as usize].change_indices().edg1 = the_i_edge;
                    data.change_pi_nod()[node as usize].change_indices().flag |=
                        nmask::OUT_L | nmask::EDGE;
                    data.change_pi_nod()[node as usize].change_data().pcu1 = u;
                    data.change_pi_nod()[node as usize].change_data().scal = 0.0;
                }
                // `cxx:2853-2888`: the four splits, in OCCT's order.
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    data.update_links_split(f1.p1, f1.p2, an_i1p3);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    data.update_links_split(f2.p1, f2.p2, an_i2p3);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    data.update_links_split(an_i2p3, f2.p2, an_i2p4);
                }
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    data.update_links_split(an_i1p3, f1.p2, an_i1p4);
                }
                // `cxx:2888-2890`.
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.update_around_node(data, an_i1p3);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.update_around_node(data, an_i2p3);
                }
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.update_around_node(data, an_i1p4);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.update_around_node(data, an_i2p4);
                }
                list.insert(
                    0,
                    HlrAlgoBiPoint::from_two_faces_with_flag(
                        seg.xt1,
                        seg.yt1,
                        seg.zt1,
                        p3.xt,
                        p3.yt,
                        p3.zt,
                        seg.x1,
                        seg.y1,
                        seg.z1,
                        p3.x,
                        p3.y,
                        p3.z,
                        the_i_edge,
                        f1.face as i32,
                        f1.p1,
                        an_i1p3,
                        f2.face as i32,
                        f2.p1,
                        an_i2p3,
                        the_flag,
                    ),
                );
                list.insert(
                    0,
                    HlrAlgoBiPoint::from_two_faces_with_flag(
                        p3.xt,
                        p3.yt,
                        p3.zt,
                        p4.xt,
                        p4.yt,
                        p4.zt,
                        p3.x,
                        p3.y,
                        p3.z,
                        p4.x,
                        p4.y,
                        p4.z,
                        the_i_edge,
                        f1.face as i32,
                        an_i1p3,
                        an_i1p4,
                        f2.face as i32,
                        an_i2p3,
                        an_i2p4,
                        the_flag,
                    ),
                );
                list.insert(
                    0,
                    HlrAlgoBiPoint::from_two_faces_with_flag(
                        p4.xt,
                        p4.yt,
                        p4.zt,
                        seg.xt2,
                        seg.yt2,
                        seg.zt2,
                        p4.x,
                        p4.y,
                        p4.z,
                        seg.x2,
                        seg.y2,
                        seg.z2,
                        the_i_edge,
                        f1.face as i32,
                        an_i1p4,
                        f1.p2,
                        f2.face as i32,
                        an_i2p4,
                        f2.p2,
                        the_flag,
                    ),
                );
            }
        } else {
            // `cxx:2956-2976`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge,
                    f1.face as i32,
                    f1.p1,
                    f1.p2,
                    f2.face as i32,
                    f2.p1,
                    f2.p2,
                    the_flag,
                ),
            );
        }
    }
}
