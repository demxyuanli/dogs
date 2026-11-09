//! One-face `Interpolation` / `MoveOrInsertPoint` and the shared crossing.
//! Source: `HLRBRep_PolyAlgo.cxx:1379-1463`, `:1857-1884`, `:1888-2112`.

use super::*;

impl HlrBRepPolyAlgo {
    /// `Interpolation(theU1, theU2, theNod1RValues, theNod2RValues, theP3,
    /// theX3..theZTI3, theCoef3, theU3, theIsNull, mP3P1)` (`cxx:1857-1884`):
    /// the crossing of the two nodes' `Scal` values, evaluated on the edge's 3D
    /// curve and projected with `TMultiply`.
    pub(super) fn interpolation_point(
        &self,
        bcurv: &Option<Arc<dyn Curve>>,
        u1: f64,
        u2: f64,
        n1: &PolyInternalNodeData,
        n2: &PolyInternalNodeData,
    ) -> Interp {
        if let Some((coef, m_p3_p1)) = Self::new_node(n1.scal, n2.scal, self.tol_angular) {
            if let Some(c) = bcurv {
                let u3 = u1 + (u2 - u1) * coef;
                let p = c.value(u3);
                let (mut x, mut y, mut z) = (p.x(), p.y(), p.z());
                let (xt, yt, zt) = (x, y, z);
                // `cxx:1878-1883`: `TMultiply` on the 3D point, `XTI3` keeps the
                // untransformed one.
                self.tmultiply(&mut x, &mut y, &mut z, false);
                return Interp {
                    x,
                    y,
                    z,
                    xt,
                    yt,
                    zt,
                    coef,
                    u: u3,
                    ins: true,
                    m_p3_p1,
                };
            }
        }
        Interp::default()
    }

    /// `HLRBRep_PolyAlgo::Interpolation(list, X1..ZTI2, theIEdge, U1, U2,
    /// Nod11, Nod12, I1p1, I1p2, I1, Pid1, ...)` (`cxx:1379-1463`), the
    /// one-face overload.
    pub(super) fn interpolation_one(
        &mut self,
        bcurv: &Option<Arc<dyn Curve>>,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: usize,
        f1: FaceSide,
        the_pid: &mut PolyInternalDataArray,
    ) {
        let n11 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p1);
        let n12 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p2);
        let p3 = self.interpolation_point(bcurv, seg.u1, seg.u2, &n11, &n12);
        self.move_or_insert_one(list, seg, the_i_edge as i32, f1, the_pid, &p3, 0);
    }

    /// The one-face `Interpolation` overload of `MoveOrInsertPoint`
    /// (`cxx:1888-2112`).
    #[allow(clippy::too_many_arguments)]
    fn move_or_insert_one(
        &mut self,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: i32,
        f1: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        p3: &Interp,
        the_flag: i32,
    ) {
        let mut an_ins3 = p3.ins;
        if an_ins3 && p3.m_p3_p1 {
            // `cxx:1935-1980`: P1 ---> P3.
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
                seg.x1 = p3.x;
                seg.y1 = p3.y;
                seg.z1 = p3.z;
                seg.xt1 = p3.xt;
                seg.yt1 = p3.yt;
                seg.zt1 = p3.zt;
                seg.u1 = p3.u;
                if let Some(front) = list.first_mut() {
                    // `cxx:1978-1979`.
                    let pts = front.change_points();
                    pts.pntp2 = GpXyz::new(p3.x, p3.y, p3.z);
                    pts.pnt2 = GpXyz::new(p3.xt, p3.yt, p3.zt);
                }
            }
        }
        if an_ins3 && !p3.m_p3_p1 {
            // `cxx:1981-2028`: P2 ---> P3.
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
            // `cxx:2029-2090`: p1 i1p3 p2.
            let an_i1p3 = {
                let data = the_pid[f1.face].as_mut().unwrap();
                let n11 = node_data(data, f1.p1);
                let n12 = node_data(data, f1.p2);
                let new_node = data.add_node(&n11, &n12, p3.coef, p3.x, p3.y, p3.z);
                let ind = data.change_pi_nod()[new_node as usize].change_indices();
                ind.edg1 = the_i_edge;
                ind.flag |= nmask::OUT_L | nmask::EDGE;
                let nd = data.change_pi_nod()[new_node as usize].change_data();
                nd.pcu1 = p3.u;
                nd.scal = 0.0;
                data.update_links_split(f1.p1, f1.p2, new_node);
                self.update_around_node(data, new_node);
                new_node
            };
            list.insert(
                0,
                HlrAlgoBiPoint::from_face1_with_flag(
                    seg.xt1, seg.yt1, seg.zt1, p3.xt, p3.yt, p3.zt, seg.x1, seg.y1, seg.z1, p3.x,
                    p3.y, p3.z, the_i_edge, f1.face as i32, f1.p1, an_i1p3, the_flag,
                ),
            );
            list.insert(
                0,
                HlrAlgoBiPoint::from_face1_with_flag(
                    p3.xt, p3.yt, p3.zt, seg.xt2, seg.yt2, seg.zt2, p3.x, p3.y, p3.z, seg.x2, seg.y2,
                    seg.z2, the_i_edge, f1.face as i32, an_i1p3, f1.p2, the_flag,
                ),
            );
        } else {
            // `cxx:2092-2111`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_face1_with_flag(
                    seg.xt1, seg.yt1, seg.zt1, seg.xt2, seg.yt2, seg.zt2, seg.x1, seg.y1, seg.z1,
                    seg.x2, seg.y2, seg.z2, the_i_edge, f1.face as i32, f1.p1, f1.p2, the_flag,
                ),
            );
        }
    }
}
