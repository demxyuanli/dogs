//! Node and adjacency layer of the outline update. Source:
//! `HLRBRep_PolyAlgo.cxx:3513-3996` (`FindEdgeOnTriangle`, `ChangeNode`,
//! `UpdateAroundNode`, `OrientTriangle`, `Triangles`, `NewNode`, `UVNode`,
//! `CheckDegeneratedSegment`) and `:719-875` (`AverageNormal`,
//! `AddNormalOnTriangle`).

use occt_core::gp::GpXyz;

use super::{
    blended_normal, nmask, xy_scaled_sum, HlrAlgoPolyInternalData, HlrAlgoPolyInternalSegment,
    HlrBRepPolyAlgo,
};
use crate::hlr_algo::poly_mask;

impl HlrBRepPolyAlgo {
    /// `FindEdgeOnTriangle(theTriangle, ip1, ip2, jtrouv, isDirect)`
    /// (`cxx:3513-3558`): which of the three triangle edges `ip1 -> ip2` is,
    /// and whether it follows the triangle winding. `None` when the edge is not
    /// one of the three.
    pub fn find_edge_on_triangle(
        the_triangle: &crate::hlr_algo::HlrAlgoTriangleData,
        ip1: i32,
        ip2: i32,
    ) -> Option<(i32, bool)> {
        let (n1, n2, n3) = (the_triangle.node1, the_triangle.node2, the_triangle.node3);
        if ip1 == n1 && ip2 == n2 {
            Some((0, true))
        } else if ip2 == n1 && ip1 == n2 {
            Some((0, false))
        } else if ip1 == n2 && ip2 == n3 {
            Some((1, true))
        } else if ip2 == n2 && ip1 == n3 {
            Some((1, false))
        } else if ip1 == n3 && ip2 == n1 {
            Some((2, true))
        } else if ip2 == n3 && ip1 == n1 {
            Some((2, false))
        } else {
            None
        }
    }

    /// `AddNormalOnTriangle(theITri, theINode, theJNode, ..., theX, theY, theZ,
    /// theIsOK)` (`cxx:802-875`).
    ///
    /// Returns `(jnode, normal)` where `jnode` is the neighbour the degenerate
    /// arm recorded (`theJNode`, `None` when no arm wrote it) and `normal` is
    /// the unit triangle normal to accumulate (`Some` iff `theIsOK` was set,
    /// i.e. all three edges are longer than `1e-10`).
    fn add_normal_on_triangle(
        data: &HlrAlgoPolyInternalData,
        itri: i32,
        inode: i32,
    ) -> (Option<i32>, Option<GpXyz>) {
        let tri = data.t_data()[itri as usize];
        let p1 = data.pi_nod()[tri.node1 as usize].data().point;
        let p2 = data.pi_nod()[tri.node2 as usize].data().point;
        let p3 = data.pi_nod()[tri.node3 as usize].data().point;
        let d1 = p2.subtracted(&p1);
        let d1_norm = d1.modulus();
        if d1_norm < 1e-10 {
            let j = if tri.node1 == inode {
                Some(tri.node2)
            } else if tri.node2 == inode {
                Some(tri.node1)
            } else {
                None
            };
            return (j, None);
        }
        let d2 = p3.subtracted(&p2);
        let d2_norm = d2.modulus();
        if d2_norm < 1e-10 {
            let j = if tri.node2 == inode {
                Some(tri.node3)
            } else if tri.node3 == inode {
                Some(tri.node2)
            } else {
                None
            };
            return (j, None);
        }
        let d3 = p1.subtracted(&p3);
        let d3_norm = d3.modulus();
        if d3_norm < 1e-10 {
            let j = if tri.node3 == inode {
                Some(tri.node1)
            } else if tri.node1 == inode {
                Some(tri.node3)
            } else {
                None
            };
            return (j, None);
        }
        let dn = d1.crossed(&d2).multiplied(1.0 / (d1_norm * d2_norm));
        if dn.modulus() > 1e-10 {
            (None, Some(dn))
        } else {
            (None, None)
        }
    }

    /// `AverageNormal(iNode, ...)` (`cxx:719-798`): the sum of the unit normals
    /// of every triangle around `iNode`. `None` where OCCT returns `false`
    /// (`theIsOK` never set, or the sum is shorter than `1e-10`).
    pub fn average_normal(data: &HlrAlgoPolyInternalData, inode: i32) -> Option<GpXyz> {
        let mut is_ok = false;
        let mut j_node = 0i32;
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        let mut iiii = data.pi_nod()[inode as usize].indices().nd_sg;
        while iiii != 0 && !is_ok {
            let seg: HlrAlgoPolyInternalSegment = data.pi_seg()[iiii as usize];
            let itri1 = seg.conex1;
            let itri2 = seg.conex2;
            if itri1 != 0 {
                let (j, n) = Self::add_normal_on_triangle(data, itri1, inode);
                if let Some(j) = j {
                    j_node = j;
                }
                if let Some(n) = n {
                    is_ok = true;
                    x += n.x();
                    y += n.y();
                    z += n.z();
                }
            }
            if itri2 != 0 {
                let (j, n) = Self::add_normal_on_triangle(data, itri2, inode);
                if let Some(j) = j {
                    j_node = j;
                }
                if let Some(n) = n {
                    is_ok = true;
                    x += n.x();
                    y += n.y();
                    z += n.z();
                }
            }
            iiii = if seg.lst_sg1 == inode {
                seg.nxt_sg1
            } else {
                seg.nxt_sg2
            };
        }

        if j_node != 0 {
            // `cxx:760-783`: the second pass walks `iNode`'s segments again but
            // accumulates around `jNode`. `kNode` is an uninitialised
            // out-parameter in OCCT (`cxx:730`).
            let mut iiii = data.pi_nod()[inode as usize].indices().nd_sg;
            while iiii != 0 && !is_ok {
                let seg: HlrAlgoPolyInternalSegment = data.pi_seg()[iiii as usize];
                let itri1 = seg.conex1;
                let itri2 = seg.conex2;
                if itri1 != 0 {
                    let (_, n) = Self::add_normal_on_triangle(data, itri1, j_node);
                    if let Some(n) = n {
                        is_ok = true;
                        x += n.x();
                        y += n.y();
                        z += n.z();
                    }
                }
                if itri2 != 0 {
                    let (_, n) = Self::add_normal_on_triangle(data, itri2, j_node);
                    if let Some(n) = n {
                        is_ok = true;
                        x += n.x();
                        y += n.y();
                        z += n.z();
                    }
                }
                iiii = if seg.lst_sg1 == j_node {
                    seg.nxt_sg1
                } else {
                    seg.nxt_sg2
                };
            }
        }

        let d = (x * x + y * y + z * z).sqrt();
        if is_ok && d < 1e-10 {
            return None;
        }
        if is_ok {
            Some(GpXyz::new(x, y, z))
        } else {
            None
        }
    }

    /// `ChangeNode(ip1, ip2, ..., theCoef1, theX3, theY3, theZ3, theIsFirst,
    /// ...)` (`cxx:3562-3627`): move node `ip1` (or `ip2` when `!is_first`) onto
    /// the interpolated point and re-orient its triangles.
    pub fn change_node(
        &mut self,
        data: &mut HlrAlgoPolyInternalData,
        ip1: i32,
        ip2: i32,
        coef1: f64,
        x3: f64,
        y3: f64,
        z3: f64,
        is_first: bool,
    ) {
        let coef2 = 1.0 - coef1;
        let (n1, n2) = Self::two_nodes(data, ip1, ip2);
        let target = if is_first { ip1 } else { ip2 };
        let node = data.change_pi_nod()[target as usize].change_data();
        node.point = GpXyz::new(x3, y3, z3);
        node.uv = xy_scaled_sum(&n1.uv, coef2, &n2.uv, coef1);
        node.scal = n1.scal * coef2 + n2.scal * coef1;
        node.normal = blended_normal(&n1.normal, coef2, &n2.normal, coef1);
        self.update_around_node(data, target);
    }

    pub(crate) fn two_nodes(
        data: &HlrAlgoPolyInternalData,
        i: i32,
        j: i32,
    ) -> (
        crate::hlr_algo::PolyInternalNodeData,
        crate::hlr_algo::PolyInternalNodeData,
    ) {
        (
            *data.pi_nod()[i as usize].data(),
            *data.pi_nod()[j as usize].data(),
        )
    }

    /// `UpdateAroundNode(theINode, ...)` (`cxx:3631-3695`): re-orient both
    /// triangles of every segment that starts at `theINode`.
    pub fn update_around_node(&mut self, data: &mut HlrAlgoPolyInternalData, inode: i32) {
        let mut iiii = data.pi_nod()[inode as usize].indices().nd_sg;
        while iiii != 0 {
            let seg: HlrAlgoPolyInternalSegment = data.pi_seg()[iiii as usize];
            if seg.conex1 != 0 {
                orient_triangle_in(data, seg.conex1, &self.proj);
            }
            if seg.conex2 != 0 {
                orient_triangle_in(data, seg.conex2, &self.proj);
            }
            iiii = if seg.lst_sg1 == inode {
                seg.nxt_sg1
            } else {
                seg.nxt_sg2
            };
        }
    }

    /// `Triangles(ip1, ip2, Nod1Indices, PISeg, iTri1, iTri2)`
    /// (`cxx:3890-3938`): the two triangles that share the edge `ip1 -> ip2`.
    pub fn triangles(data: &HlrAlgoPolyInternalData, ip1: i32, ip2: i32) -> Option<(i32, i32)> {
        let mut iiii = data.pi_nod()[ip1 as usize].indices().nd_sg;
        while iiii != 0 {
            let seg: HlrAlgoPolyInternalSegment = data.pi_seg()[iiii as usize];
            if seg.lst_sg1 == ip1 {
                if seg.lst_sg2 == ip2 {
                    return Some((seg.conex1, seg.conex2));
                }
                iiii = seg.nxt_sg1;
            } else if seg.lst_sg1 == ip2 {
                return Some((seg.conex1, seg.conex2));
            } else {
                iiii = seg.nxt_sg2;
            }
        }
        None
    }

    /// `NewNode(theNod1RValues, theNod2RValues, theCoef1, theToMoveP1)`
    /// (`cxx:3942-3960`): the parameter where the two scalars cross, when the
    /// pair straddles `±TolAngular/2`. `tol_angular` is the enclosing algorithm's
    /// `myTolAngular` (`cxx:3946`), which the OCCT member reads directly.
    pub fn new_node(nod1_scal: f64, nod2_scal: f64, tol_angular: f64) -> Option<(f64, bool)> {
        let a_tol_ang = tol_angular * 0.5;
        if (nod1_scal >= a_tol_ang && nod2_scal <= -a_tol_ang)
            || (nod2_scal >= a_tol_ang && nod1_scal <= -a_tol_ang)
        {
            let mut coef1 = nod1_scal / (nod2_scal - nod1_scal);
            if coef1 < 0.0 {
                coef1 = -coef1;
            }
            let move_p1 = coef1 < 0.5;
            Some((coef1, move_p1))
        } else {
            None
        }
    }

    /// `UVNode(theNod1RValues, theNod2RValues, theCoef1, theU3, theV3)`
    /// (`cxx:3964-3974`).
    pub fn uv_node(nod1_uv: (f64, f64), nod2_uv: (f64, f64), coef1: f64) -> (f64, f64) {
        let coef2 = 1.0 - coef1;
        (
            coef2 * nod1_uv.0 + coef1 * nod2_uv.0,
            coef2 * nod1_uv.1 + coef1 * nod2_uv.1,
        )
    }

    /// `CheckDegeneratedSegment(...)` (`cxx:3978-3994`): mark both nodes `Fuck`
    /// and, when the pair straddles `±TolAngular`, zero their scalars and mark
    /// them `OutL`.
    pub fn check_degenerated_segment(
        data: &mut HlrAlgoPolyInternalData,
        tol_angular: f64,
        ip1: i32,
        ip2: i32,
    ) {
        let (s1, s2) = {
            let n1 = data.pi_nod()[ip1 as usize].data().scal;
            let n2 = data.pi_nod()[ip2 as usize].data().scal;
            (n1, n2)
        };
        data.change_pi_nod()[ip1 as usize].change_indices().flag |= nmask::FUCK;
        data.change_pi_nod()[ip2 as usize].change_indices().flag |= nmask::FUCK;
        if (s1 >= tol_angular && s2 <= -tol_angular) || (s2 >= tol_angular && s1 <= -tol_angular)
        {
            data.change_pi_nod()[ip1 as usize].change_data().scal = 0.0;
            data.change_pi_nod()[ip1 as usize].change_indices().flag |= nmask::OUT_L;
            data.change_pi_nod()[ip2 as usize].change_data().scal = 0.0;
            data.change_pi_nod()[ip2 as usize].change_indices().flag |= nmask::OUT_L;
        }
    }
}

/// `OrientTriangle(iTri, ...)` (`cxx:3699-3886`) applied to triangle `itri` of
/// `data`. Free function because it needs three node slots plus the triangle
/// record of the same owner. `a_s = as3` (`cxx:3761`) is a dead store in OCCT
/// too; it is kept so the branch matches the source.
#[allow(unused_assignments)]
fn orient_triangle_in(
    data: &mut HlrAlgoPolyInternalData,
    itri: i32,
    proj: &crate::hlr_algo::HlrAlgoProjector,
) {
    let tri = data.t_data()[itri as usize];
    let (i1, i2, i3) = (tri.node1, tri.node2, tri.node3);
    let nd1 = *data.pi_nod()[i1 as usize].data();
    let nd2 = *data.pi_nod()[i2 as usize].data();
    let nd3 = *data.pi_nod()[i3 as usize].data();
    let o1 = data.pi_nod()[i1 as usize].indices().flag & nmask::OUT_L != 0;
    let o2 = data.pi_nod()[i2 as usize].indices().flag & nmask::OUT_L != 0;
    let o3 = data.pi_nod()[i3 as usize].indices().flag & nmask::OUT_L != 0;

    let mut flags = tri.flags;
    flags &= !poly_mask::FMSK_FLAT;
    flags &= !poly_mask::FMSK_ON_OUT_L;

    if o1 && o2 && o3 {
        flags |= poly_mask::FMSK_SIDE;
        flags &= !poly_mask::FMSK_BACK;
        flags |= poly_mask::FMSK_ON_OUT_L;
    } else {
        let (s1, s2, s3) = (nd1.scal, nd2.scal, nd3.scal);
        let mut as1 = s1;
        let mut as2 = s2;
        let mut as3 = s3;
        if s1 < 0.0 {
            as1 = -s1;
        }
        if s2 < 0.0 {
            as2 = -s2;
        }
        if s3 < 0.0 {
            as3 = -s3;
        }
        let mut s = 0.0;
        let mut a_s = 0.0;
        if !o1 {
            s = s1;
            a_s = as1;
        }
        if !o2 && a_s < as2 {
            s = s2;
            a_s = as2;
        }
        if !o3 && a_s < as3 {
            s = s3;
            a_s = as3;
        }
        if s > 0.0 {
            flags &= !poly_mask::FMSK_SIDE;
            flags |= poly_mask::FMSK_BACK;
        } else {
            flags &= !poly_mask::FMSK_SIDE;
            flags &= !poly_mask::FMSK_BACK;
        }
        let mut d12 = nd2.point.subtracted(&nd1.point);
        let d12_norm = d12.modulus();
        if d12_norm <= 1e-10 {
            flags |= poly_mask::FMSK_FLAT;
            flags |= poly_mask::FMSK_SIDE;
            flags &= !poly_mask::FMSK_BACK;
        } else {
            let mut d23 = nd3.point.subtracted(&nd2.point);
            let d23_norm = d23.modulus();
            if d23_norm < 1e-10 {
                flags |= poly_mask::FMSK_FLAT;
                flags |= poly_mask::FMSK_SIDE;
                flags &= !poly_mask::FMSK_BACK;
            } else {
                let d31 = nd1.point.subtracted(&nd3.point);
                let d31_norm = d31.modulus();
                if d31_norm < 1e-10 {
                    flags |= poly_mask::FMSK_FLAT;
                    flags |= poly_mask::FMSK_SIDE;
                    flags &= !poly_mask::FMSK_BACK;
                } else {
                    d12 = d12.multiplied(1.0 / d12_norm);
                    d23 = d23.multiplied(1.0 / d23_norm);
                    let mut d = d12.crossed(&d23);
                    let d_norm = d.modulus();
                    if d_norm < 1e-5 {
                        flags |= poly_mask::FMSK_FLAT;
                        flags |= poly_mask::FMSK_SIDE;
                        flags &= !poly_mask::FMSK_BACK;
                    } else {
                        let mut o;
                        if proj.perspective() {
                            d = d.multiplied(1.0 / d_norm);
                            // `myProj.Focus()` raises when the projection is not
                            // perspective (`Projector.cxx`); this arm is guarded
                            // by `Perspective() == true`, so the `Err` is
                            // unreachable.
                            let focus = proj.focus().unwrap_or(0.0);
                            o = d.z() * focus - d.dot(&nd1.point);
                        } else {
                            o = d.z() / d_norm;
                        }
                        if o < 0.0 {
                            flags |= poly_mask::FMSK_OR_BACK;
                            o = -o;
                        } else {
                            flags &= !poly_mask::FMSK_OR_BACK;
                        }
                        if o < 1e-10 {
                            flags |= poly_mask::FMSK_SIDE;
                            flags &= !poly_mask::FMSK_BACK;
                        }
                    }
                }
            }
        }
    }

    if (flags & poly_mask::FMSK_BACK == 0 && flags & poly_mask::FMSK_OR_BACK != 0)
        || (flags & poly_mask::FMSK_BACK != 0 && flags & poly_mask::FMSK_OR_BACK == 0)
    {
        flags |= poly_mask::FMSK_FR_BACK;
    } else {
        flags &= !poly_mask::FMSK_FR_BACK;
    }

    data.change_t_data()[itri as usize].flags = flags;
}
