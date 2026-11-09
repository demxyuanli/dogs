//! Outline update data of one polygon. Source: `HLRAlgo_PolyInternalData`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_PolyInternalData.hxx:30-116`,
//! `.lxx:28-117`, `.cxx:32-1039`).
//!
//! Storage convention: OCCT allocates `NCollection_HArray1(0, Mx)` and uses
//! indices `1..Nb`; the port keeps `Vec`s of length `Mx + 1` whose slot `0` is
//! the unused padding, so `t_data[i]` is `TData(i)`.
//!
//! The OCCT overloads hand out `NCollection_Array1*&` out-parameters purely so
//! that callers notice the reallocation performed by `IncTData` / `IncPISeg` /
//! `IncPINod`. The arrays belong to this object, so the port drops those
//! out-parameters and re-reads the members after every `inc_*`.
//!
//! UNPORTED: `Dump()` (`cxx:877-919`) — an `OCCT_DEBUG`-only `std::cout`
//! transcript of the three arrays, with no effect on the geometry.

use occt_core::gp::{GpXY, GpXyz};

use super::poly_internal_records::{
    HlrAlgoPolyInternalNode, HlrAlgoPolyInternalSegment, PolyInternalNodeData,
};
use super::poly_records::HlrAlgoTriangleData;

/// Outline update data of one polygon. Source: `HLRAlgo_PolyInternalData`.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrAlgoPolyInternalData {
    nb_t_data: i32,
    nb_pi_seg: i32,
    nb_pi_nod: i32,
    mx_t_data: i32,
    mx_pi_seg: i32,
    mx_pi_nod: i32,
    int_out_l: bool,
    planar: bool,
    t_data: Vec<HlrAlgoTriangleData>,
    pi_seg: Vec<HlrAlgoPolyInternalSegment>,
    pi_nod: Vec<HlrAlgoPolyInternalNode>,
}

impl HlrAlgoPolyInternalData {
    /// `HLRAlgo_PolyInternalData(nbNod, nbTri)` (`cxx:32-57`).
    pub fn new(nb_nod: i32, nb_tri: i32) -> Self {
        let mx_pi_seg = 2 + (3 * nb_tri + nb_nod) / 2;
        let mut pi_nod = vec![HlrAlgoPolyInternalNode::new(); mx_pi_nod_len(nb_nod)];
        for i in 1..=nb_nod as usize {
            pi_nod[i] = HlrAlgoPolyInternalNode::new();
        }
        Self {
            nb_t_data: nb_tri,
            nb_pi_seg: 0,
            nb_pi_nod: nb_nod,
            mx_t_data: nb_tri,
            mx_pi_seg: mx_pi_seg,
            mx_pi_nod: nb_nod,
            int_out_l: false,
            planar: false,
            t_data: vec![HlrAlgoTriangleData::default(); mx_t_data_len(nb_tri)],
            pi_seg: vec![HlrAlgoPolyInternalSegment::default(); mx_pi_seg_len(mx_pi_seg)],
            pi_nod,
        }
    }

    /// `DecTData()` (`lxx:28-31`).
    pub fn dec_t_data(&mut self) {
        self.nb_t_data -= 1;
    }

    /// `DecPISeg()` (`lxx:34-37`).
    pub fn dec_pi_seg(&mut self) {
        self.nb_pi_seg -= 1;
    }

    /// `DecPINod()` (`lxx:41-44`).
    pub fn dec_pi_nod(&mut self) {
        self.nb_pi_nod -= 1;
    }

    /// `NbTData()` (`lxx:51-54`).
    pub fn nb_t_data(&self) -> i32 {
        self.nb_t_data
    }

    /// `NbPISeg()` (`lxx:57-60`).
    pub fn nb_pi_seg(&self) -> i32 {
        self.nb_pi_seg
    }

    /// `NbPINod()` (`lxx:64-67`).
    pub fn nb_pi_nod(&self) -> i32 {
        self.nb_pi_nod
    }

    /// `Planar()` (`lxx:72-75`).
    pub fn planar(&self) -> bool {
        self.planar
    }

    /// `Planar(B)` (`lxx:78-81`).
    pub fn set_planar(&mut self, b: bool) {
        self.planar = b;
    }

    /// `IntOutL()` (`lxx:86-89`).
    pub fn int_out_l(&self) -> bool {
        self.int_out_l
    }

    /// `IntOutL(B)` (`lxx:92-95`).
    pub fn set_int_out_l(&mut self, b: bool) {
        self.int_out_l = b;
    }

    /// `TData()` (`lxx:98-101`); `t_data()[i]` is `TData(i)`.
    pub fn t_data(&self) -> &[HlrAlgoTriangleData] {
        &self.t_data
    }

    /// `TData()` mutable form.
    pub fn change_t_data(&mut self) -> &mut Vec<HlrAlgoTriangleData> {
        &mut self.t_data
    }

    /// `PINod()` (`lxx:104-110`); `pi_nod()[i]` is `PINod(i)`.
    pub fn pi_nod(&self) -> &[HlrAlgoPolyInternalNode] {
        &self.pi_nod
    }

    /// `PINod()` mutable form.
    pub fn change_pi_nod(&mut self) -> &mut Vec<HlrAlgoPolyInternalNode> {
        &mut self.pi_nod
    }

    /// `PISeg()` (`lxx:113-117`); `pi_seg()[i]` is `PISeg(i)`.
    pub fn pi_seg(&self) -> &[HlrAlgoPolyInternalSegment] {
        &self.pi_seg
    }

    /// `PISeg()` mutable form.
    pub fn change_pi_seg(&mut self) -> &mut Vec<HlrAlgoPolyInternalSegment> {
        &mut self.pi_seg
    }

    /// `UpdateLinks(theTData, thePISeg, thePINod)` (`cxx:60-429`): rebuild the
    /// node / segment adjacency of the first `NbTData()` triangles.
    pub fn update_links(&mut self) {
        let mut icsv = 0i32;
        let nb_t = self.nb_t_data;
        for i in 1..=nb_t {
            let td = self.t_data[i as usize];
            self.link_edge(i, td.node1, td.node2, &mut icsv);
            let td = self.t_data[i as usize];
            self.link_edge(i, td.node2, td.node3, &mut icsv);
            let td = self.t_data[i as usize];
            self.link_edge(i, td.node3, td.node1, &mut icsv);
        }
    }

    /// One of the three identical blocks of `UpdateLinks` (`cxx:77-193`,
    /// `:196-309`, `:312-427`): create or find the segment `n1 -> n2` of
    /// triangle `tri`.
    ///
    /// `icsv` is declared once in the OCCT function (`cxx:68`) and therefore
    /// carries over between blocks; it is threaded through unchanged.
    fn link_edge(&mut self, tri: i32, n1: i32, n2: i32, icsv: &mut i32) {
        let mut new_seg = false;
        let nd1 = self.pi_nod[n1 as usize].indices().nd_sg;
        let nd2 = self.pi_nod[n2 as usize].indices().nd_sg;
        if nd1 == 0 && nd2 == 0 {
            new_seg = true;
            self.nb_pi_seg += 1;
            self.pi_nod[n1 as usize].change_indices().nd_sg = self.nb_pi_seg;
            self.pi_nod[n2 as usize].change_indices().nd_sg = self.nb_pi_seg;
        } else {
            let mut iiii = nd1;
            if iiii != 0 {
                let mut find = 0i32;
                // `aSegIndices` of OCCT keeps pointing at the last segment
                // walked once `iiii` drops to 0 (`cxx:129`, `:246`, `:362`).
                let mut last = 0i32;
                while iiii != 0 && find == 0 {
                    last = iiii;
                    let seg = self.pi_seg[iiii as usize];
                    if seg.lst_sg1 == n1 {
                        if seg.lst_sg2 == n2 {
                            find = iiii;
                        } else {
                            iiii = seg.nxt_sg1;
                        }
                        *icsv = 1;
                    } else {
                        if seg.lst_sg1 == n2 {
                            find = iiii;
                        } else {
                            iiii = seg.nxt_sg2;
                        }
                        *icsv = 2;
                    }
                }
                if find == 0 {
                    new_seg = true;
                    self.nb_pi_seg += 1;
                    if *icsv == 1 {
                        self.pi_seg[last as usize].nxt_sg1 = self.nb_pi_seg;
                    } else {
                        self.pi_seg[last as usize].nxt_sg2 = self.nb_pi_seg;
                    }
                } else {
                    self.pi_seg[find as usize].conex2 = tri;
                }
            } else {
                new_seg = true;
                self.nb_pi_seg += 1;
                self.pi_nod[n1 as usize].change_indices().nd_sg = self.nb_pi_seg;
            }
            if new_seg {
                let mut iiii = self.pi_nod[n2 as usize].indices().nd_sg;
                if iiii != 0 {
                    let mut last = 0i32;
                    while iiii != 0 {
                        last = iiii;
                        let seg = self.pi_seg[iiii as usize];
                        if seg.lst_sg1 == n2 {
                            *icsv = 1;
                            iiii = seg.nxt_sg1;
                        } else {
                            *icsv = 2;
                            iiii = seg.nxt_sg2;
                        }
                    }
                    if *icsv == 1 {
                        self.pi_seg[last as usize].nxt_sg1 = self.nb_pi_seg;
                    } else {
                        self.pi_seg[last as usize].nxt_sg2 = self.nb_pi_seg;
                    }
                } else {
                    self.pi_nod[n2 as usize].change_indices().nd_sg = self.nb_pi_seg;
                }
            }
        }
        if new_seg {
            let seg = &mut self.pi_seg[self.nb_pi_seg as usize];
            seg.lst_sg1 = n1;
            seg.lst_sg2 = n2;
            seg.conex1 = tri;
            seg.conex2 = 0;
            seg.nxt_sg1 = 0;
            seg.nxt_sg2 = 0;
        }
    }

    /// `AddNode(theNod1RValues, theNod2RValues, ..., theCoef1, theX3, theY3,
    /// theZ3)` (`cxx:435-472`): interpolate a new node on the edge
    /// `Nod1 -> Nod2` and return its 1-based index.
    pub fn add_node(
        &mut self,
        nod1: &PolyInternalNodeData,
        nod2: &PolyInternalNodeData,
        coef1: f64,
        x3: f64,
        y3: f64,
        z3: f64,
    ) -> i32 {
        let coef2 = 1.0 - coef1;
        self.inc_pi_nod();
        let ip3 = self.nb_pi_nod;
        self.pi_nod[ip3 as usize] = HlrAlgoPolyInternalNode::new();
        {
            let node = &mut self.pi_nod[ip3 as usize];
            node.change_indices().nd_sg = 0;
            node.change_indices().flag = 0;
            let n = node.change_data();
            n.point = GpXyz::new(x3, y3, z3);
            n.uv = GpXY::new(
                coef2 * nod1.uv.x() + coef1 * nod2.uv.x(),
                coef2 * nod1.uv.y() + coef1 * nod2.uv.y(),
            );
            n.scal = nod1.scal * coef2 + nod2.scal * coef1;
            let a_xyz = GpXyz::new(
                coef2 * nod1.normal.x() + coef1 * nod2.normal.x(),
                coef2 * nod1.normal.y() + coef1 * nod2.normal.y(),
                coef2 * nod1.normal.z() + coef1 * nod2.normal.z(),
            );
            let a_norm = a_xyz.modulus();
            n.normal = if a_norm > 0.0 {
                GpXyz::new(a_xyz.x() / a_norm, a_xyz.y() / a_norm, a_xyz.z() / a_norm)
            } else {
                GpXyz::new(1.0, 0.0, 0.0)
            };
        }
        ip3
    }

    /// `UpdateLinks(ip1, ip2, ip3, ...)` (`cxx:477-872`): split the edge
    /// `ip1 -> ip2` at the new node `ip3`, splitting every incident triangle.
    pub fn update_links_split(&mut self, ip1: i32, ip2: i32, ip3: i32) {
        let mut find = 0i32;
        let mut iisv = 0i32;
        let mut icsv = 0i32;
        self.inc_pi_seg();
        self.inc_pi_seg();
        self.inc_pi_seg();
        self.nb_pi_seg -= 2;
        self.inc_t_data();
        self.inc_t_data();
        self.nb_t_data -= 2;

        let mut iip2 = 0i32;
        let mut cnx1 = 0i32;
        let mut cnx2 = 0i32;
        let mut iiii = self.pi_nod[ip2 as usize].indices().nd_sg;
        // OCCT keeps `aSegIndices` = the segment of the previous iteration
        // (`cxx:569`), which is the one whose `NxtSg*` link is rewritten.
        let mut prev = 0i32;

        while iiii != 0 && find == 0 {
            let seg2 = self.pi_seg[iiii as usize];
            if seg2.lst_sg1 == ip2 {
                if seg2.lst_sg2 == ip1 {
                    find = iiii;
                    cnx1 = seg2.conex1;
                    cnx2 = seg2.conex2;
                    self.pi_seg[iiii as usize].lst_sg1 = ip3;
                    iip2 = seg2.nxt_sg1;
                    self.pi_seg[iiii as usize].nxt_sg1 = self.nb_pi_seg;
                    if iisv == 0 {
                        self.pi_nod[ip2 as usize].change_indices().nd_sg = self.nb_pi_seg;
                    } else if icsv == 1 {
                        self.pi_seg[prev as usize].nxt_sg1 = self.nb_pi_seg;
                    } else {
                        self.pi_seg[prev as usize].nxt_sg2 = self.nb_pi_seg;
                    }
                } else {
                    iisv = iiii;
                    icsv = 1;
                }
                iiii = seg2.nxt_sg1;
            } else if seg2.lst_sg1 == ip1 {
                find = iiii;
                cnx1 = seg2.conex1;
                cnx2 = seg2.conex2;
                self.pi_seg[iiii as usize].lst_sg2 = ip3;
                iip2 = seg2.nxt_sg2;
                self.pi_seg[iiii as usize].nxt_sg2 = self.nb_pi_seg;
                if iisv == 0 {
                    self.pi_nod[ip2 as usize].change_indices().nd_sg = self.nb_pi_seg;
                } else if icsv == 1 {
                    self.pi_seg[prev as usize].nxt_sg1 = self.nb_pi_seg;
                } else {
                    self.pi_seg[prev as usize].nxt_sg2 = self.nb_pi_seg;
                }
            } else {
                iisv = iiii;
                icsv = 2;
            }
            prev = iiii;
            iiii = if seg2.lst_sg1 == ip2 {
                seg2.nxt_sg1
            } else {
                seg2.nxt_sg2
            };
        }

        if find == 0 {
            self.nb_pi_seg -= 1;
            return;
        }

        {
            let seg = &mut self.pi_seg[self.nb_pi_seg as usize];
            seg.nxt_sg1 = 0;
            seg.nxt_sg2 = iip2;
            seg.lst_sg1 = ip3;
            seg.lst_sg2 = ip2;
            seg.conex1 = cnx1;
            seg.conex2 = cnx2;
        }
        self.pi_nod[ip3 as usize].change_indices().nd_sg = find;

        // `cxx:593-607`.
        let mut new = [0i32; 4];
        new[0] = cnx1;
        new[2] = self.nb_t_data + 1;
        if cnx2 == 0 {
            new[1] = 0;
            new[3] = 0;
        } else {
            new[1] = cnx2;
            new[3] = self.nb_t_data + 2;
        }

        let mut itpk = [0i32; 2];
        for skip in 0..=1usize {
            let i_old = new[skip];
            let i_new = new[skip + 2];
            if i_old == 0 {
                continue;
            }
            let a_triangle = self.t_data[i_old as usize];
            let (n1, n2, n3) = (a_triangle.node1, a_triangle.node2, a_triangle.node3);
            let mut n_old = [n1, n2, n3];
            let mut n_new = [n1, n2, n3];
            let mut found = false;
            if n1 == ip1 && n2 == ip2 {
                found = true;
                n_old[1] = ip3;
                n_new[0] = ip3;
                itpk[skip] = n3;
            } else if n1 == ip2 && n2 == ip1 {
                found = true;
                n_old[0] = ip3;
                n_new[1] = ip3;
                itpk[skip] = n3;
            } else if n2 == ip1 && n3 == ip2 {
                found = true;
                n_old[2] = ip3;
                n_new[1] = ip3;
                itpk[skip] = n1;
            } else if n2 == ip2 && n3 == ip1 {
                found = true;
                n_old[1] = ip3;
                n_new[2] = ip3;
                itpk[skip] = n1;
            } else if n3 == ip1 && n1 == ip2 {
                found = true;
                n_old[0] = ip3;
                n_new[2] = ip3;
                itpk[skip] = n2;
            } else if n3 == ip2 && n1 == ip1 {
                found = true;
                n_old[2] = ip3;
                n_new[0] = ip3;
                itpk[skip] = n2;
            }
            if !found {
                continue;
            }
            self.nb_t_data += 1;
            let ip4 = itpk[skip];
            {
                let tri = &mut self.t_data[i_old as usize];
                tri.node1 = n_old[0];
                tri.node2 = n_old[1];
                tri.node3 = n_old[2];
            }
            {
                let flags = self.t_data[i_old as usize].flags;
                let tri2 = &mut self.t_data[i_new as usize];
                tri2.node1 = n_new[0];
                tri2.node2 = n_new[1];
                tri2.node3 = n_new[2];
                tri2.flags = flags;
            }
            self.nb_pi_seg += 1;
            let nd3 = self.pi_nod[ip3 as usize].indices().nd_sg;
            let nd4 = self.pi_nod[ip4 as usize].indices().nd_sg;
            {
                let seg = &mut self.pi_seg[self.nb_pi_seg as usize];
                seg.lst_sg1 = ip3;
                seg.lst_sg2 = ip4;
                seg.nxt_sg1 = nd3;
                seg.nxt_sg2 = nd4;
                seg.conex1 = i_old;
                seg.conex2 = i_new;
            }
            self.pi_nod[ip3 as usize].change_indices().nd_sg = self.nb_pi_seg;
            self.pi_nod[ip4 as usize].change_indices().nd_sg = self.nb_pi_seg;
        }

        // `cxx:701-871`.
        for i_tr in 0..=3usize {
            let i_new = new[i_tr];
            let skip = if i_tr < 2 { i_tr } else { i_tr - 2 };
            let i_old = new[skip];
            let ip4 = itpk[skip];
            if i_new == 0 {
                continue;
            }
            let a_triangle = self.t_data[i_new as usize];
            let (n1, n2, n3) = (a_triangle.node1, a_triangle.node2, a_triangle.node3);
            if !((n1 == ip3 && n2 == ip4) || (n2 == ip3 && n1 == ip4)) {
                fix_conex(&mut self.pi_seg, &self.pi_nod, n1, n2, i_old, i_new);
            }
            if !((n2 == ip3 && n3 == ip4) || (n3 == ip3 && n2 == ip4)) {
                fix_conex(&mut self.pi_seg, &self.pi_nod, n2, n3, i_old, i_new);
            }
            if !((n3 == ip3 && n1 == ip4) || (n1 == ip3 && n3 == ip4)) {
                fix_conex(&mut self.pi_seg, &self.pi_nod, n3, n1, i_old, i_new);
            }
        }
    }

    /// `IncTData(...)` (`cxx:923-958`): double the array, then `NbTData++`.
    pub fn inc_t_data(&mut self) {
        if self.nb_t_data >= self.mx_t_data {
            let k = 2 * self.mx_t_data;
            let mut nw = vec![HlrAlgoTriangleData::default(); mx_t_data_len(k)];
            for i in 1..=self.mx_t_data as usize {
                nw[i] = self.t_data[i];
            }
            self.mx_t_data = k;
            self.t_data = nw;
        }
        self.nb_t_data += 1;
    }

    /// `IncPISeg(...)` (`cxx:961-995`): double the array, then `NbPISeg++`.
    pub fn inc_pi_seg(&mut self) {
        if self.nb_pi_seg >= self.mx_pi_seg {
            let k = 2 * self.mx_pi_seg;
            let mut nw = vec![HlrAlgoPolyInternalSegment::default(); mx_pi_seg_len(k)];
            for i in 1..=self.mx_pi_seg as usize {
                nw[i] = self.pi_seg[i];
            }
            self.mx_pi_seg = k;
            self.pi_seg = nw;
        }
        self.nb_pi_seg += 1;
    }

    /// `IncPINod(...)` (`cxx:1000-1038`): double the array, then `NbPINod++`.
    pub fn inc_pi_nod(&mut self) {
        if self.nb_pi_nod >= self.mx_pi_nod {
            let k = 2 * self.mx_pi_nod;
            let mut nw = vec![HlrAlgoPolyInternalNode::new(); mx_pi_nod_len(k)];
            for i in 1..=self.mx_pi_nod as usize {
                nw[i] = self.pi_nod[i];
            }
            self.mx_pi_nod = k;
            self.pi_nod = nw;
        }
        self.nb_pi_nod += 1;
    }
}

/// `cxx:722-768`: rewrite the `Conex1` / `Conex2` of the segment `n1 -> n2`
/// from triangle `i_old` to triangle `i_new`.
fn fix_conex(
    pi_seg: &mut [HlrAlgoPolyInternalSegment],
    pi_nod: &[HlrAlgoPolyInternalNode],
    n1: i32,
    n2: i32,
    i_old: i32,
    i_new: i32,
) {
    let mut found = false;
    let mut iiii = pi_nod[n1 as usize].indices().nd_sg;
    while iiii != 0 && !found {
        let seg = pi_seg[iiii as usize];
        if seg.lst_sg1 == n1 {
            if seg.lst_sg2 == n2 {
                found = true;
                if seg.conex1 == i_old {
                    pi_seg[iiii as usize].conex1 = i_new;
                } else if seg.conex2 == i_old {
                    pi_seg[iiii as usize].conex2 = i_new;
                }
            } else {
                iiii = seg.nxt_sg1;
            }
        } else if seg.lst_sg1 == n2 {
            found = true;
            if seg.conex1 == i_old {
                pi_seg[iiii as usize].conex1 = i_new;
            } else if seg.conex2 == i_old {
                pi_seg[iiii as usize].conex2 = i_new;
            }
        } else {
            iiii = seg.nxt_sg2;
        }
    }
}

/// `NCollection_Array1(0, mx)`: `mx + 1` slots.
fn mx_t_data_len(mx: i32) -> usize {
    mx as usize + 1
}

fn mx_pi_seg_len(mx: i32) -> usize {
    mx as usize + 1
}

fn mx_pi_nod_len(mx: i32) -> usize {
    mx as usize + 1
}
