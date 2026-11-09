//! Final assembly of the poly data consumed by `HLRAlgo_PolyAlgo`. Source:
//! `HLRBRep_PolyAlgo::UpdateEdgesBiPoints` (`cxx:4155-4239`) and
//! `UpdatePolyData` (`cxx:4241-4329`).

use occt_core::gp::GpXyz;

use super::{HlrBRepPolyAlgo, PolyInternalDataArray};
use crate::hlr_algo::{
    poly_mask, BiPointIndices, HlrAlgoBiPoint, HlrAlgoPolyData, HlrAlgoPolyHidingData,
    HlrAlgoTriangleData,
};

impl HlrBRepPolyAlgo {
    /// `UpdateEdgesBiPoints(theList, thePID, theIsClosed)` (`cxx:4155-4239`):
    /// resolve the `OutLine` flag of every two-face bi-point from the flags of
    /// the two triangles that carry the edge.
    pub fn update_edges_bi_points(
        &self,
        the_list: &mut [HlrAlgoBiPoint],
        the_pid: &[Option<crate::hlr_algo::HlrAlgoPolyInternalData>],
        the_is_closed: bool,
    ) {
        for bp in the_list.iter_mut() {
            let ind: BiPointIndices = *bp.indices();
            if ind.face_conex1 == 0 || ind.face_conex2 == 0 {
                continue;
            }
            let pid1 = the_pid[ind.face_conex1 as usize].as_ref().unwrap();
            let pid2 = the_pid[ind.face_conex2 as usize].as_ref().unwrap();
            // `cxx:4177-4178`: `Triangles(...)` reports 0 when the edge is not
            // shared; `tbid` is an unused out-parameter in OCCT.
            let t1 = Self::triangles(pid1, ind.face1_pt1, ind.face1_pt2);
            let t2 = Self::triangles(pid2, ind.face2_pt1, ind.face2_pt2);
            let (Some((itri1, _)), Some((itri2, _))) = (t1, t2) else {
                continue;
            };
            if itri1 == 0 || itri2 == 0 {
                continue;
            }
            if ind.face_conex1 == ind.face_conex2 && itri1 == itri2 {
                continue;
            }
            let tri1 = pid1.t_data()[itri1 as usize];
            let tri2 = pid2.t_data()[itri2 as usize];
            if the_is_closed
                && (((tri1.flags & poly_mask::FMSK_BACK != 0)
                    && (tri2.flags & poly_mask::FMSK_BACK != 0))
                    || ((tri1.flags & poly_mask::FMSK_SIDE != 0)
                        && (tri2.flags & poly_mask::FMSK_SIDE != 0))
                    || ((tri1.flags & poly_mask::FMSK_BACK != 0)
                        && (tri2.flags & poly_mask::FMSK_SIDE != 0))
                    || ((tri1.flags & poly_mask::FMSK_SIDE != 0)
                        && (tri2.flags & poly_mask::FMSK_BACK != 0)))
            {
                bp.set_hidden(true);
            }
            let is_outl = if tri1.flags & poly_mask::FMSK_SIDE == 0
                && tri2.flags & poly_mask::FMSK_SIDE == 0
            {
                (tri1.flags & poly_mask::FMSK_BACK) != (tri2.flags & poly_mask::FMSK_BACK)
            } else if tri1.flags & poly_mask::FMSK_SIDE != 0
                && tri2.flags & poly_mask::FMSK_SIDE != 0
            {
                false
            } else if tri1.flags & poly_mask::FMSK_SIDE != 0 {
                tri1.flags & poly_mask::FMSK_FLAT == 0
                    && tri2.flags & poly_mask::FMSK_BACK == 0
            } else {
                tri2.flags & poly_mask::FMSK_FLAT == 0
                    && tri1.flags & poly_mask::FMSK_BACK == 0
            };
            bp.set_out_line(is_outl);
        }
    }

    /// `UpdatePolyData(thePD, thePID, theIsClosed)` (`cxx:4241-4329`): publish
    /// the internal nodes and triangles of every face into `thePD`, reversing
    /// the winding of the triangles marked `FMskOrBack` and flagging the hiding
    /// ones.
    ///
    /// OCCT aliases the `thePD` handles into `HLRAlgo_PolyShellData::PolyData()`
    /// (`cxx:499-508`); the port copies them into the shell (see
    /// `StoreShell`), so nothing here reads the shell array.
    pub fn update_poly_data(
        &self,
        the_pd: &mut [HlrAlgoPolyData],
        the_pid: &mut PolyInternalDataArray,
        the_is_closed: bool,
    ) {
        let a_nb_faces = self.f_map.extent();
        for a_face_iter in 1..=a_nb_faces {
            let Some(pid) = the_pid[a_face_iter].as_mut() else {
                continue;
            };
            let a_nb_nodes = pid.nb_pi_nod();
            let a_nb_tris = pid.nb_t_data();
            let mut nodes: Vec<GpXyz> = Vec::with_capacity(a_nb_nodes as usize);
            for a_node_iter in 1..=a_nb_nodes {
                nodes.push(pid.pi_nod()[a_node_iter as usize].data().point);
            }
            let mut trian: Vec<HlrAlgoTriangleData> = Vec::with_capacity(a_nb_tris.max(0) as usize);
            let mut a_nb_hidden = 0usize;
            for a_tri_iter in 1..=a_nb_tris {
                let mut ot = pid.t_data()[a_tri_iter as usize];
                if ot.flags & poly_mask::FMSK_SIDE == 0 {
                    if ot.flags & poly_mask::FMSK_OR_BACK != 0 {
                        let j = ot.node1;
                        ot.node1 = ot.node3;
                        ot.node3 = j;
                        ot.flags |= poly_mask::FMSK_BACK;
                    } else {
                        ot.flags &= !poly_mask::FMSK_BACK;
                        // Tri1Flags |= HLRAlgo_PolyMask_FMskBack;//OCC349
                    }
                }
                pid.change_t_data()[a_tri_iter as usize] = ot;
                let mut nt = HlrAlgoTriangleData::default();
                nt.node1 = ot.node1;
                nt.node2 = ot.node2;
                nt.node3 = ot.node3;
                nt.flags = ot.flags;
                if nt.flags & poly_mask::FMSK_SIDE == 0
                    && (nt.flags & poly_mask::FMSK_BACK == 0 || !the_is_closed)
                {
                    nt.flags |= poly_mask::FMSK_HIDING;
                    a_nb_hidden += 1;
                } else {
                    nt.flags &= !poly_mask::FMSK_HIDING;
                }
                trian.push(nt);
            }

            let ph_dat: Option<Vec<HlrAlgoPolyHidingData>> = if a_nb_hidden > 0 {
                Some(vec![HlrAlgoPolyHidingData::new(); a_nb_hidden])
            } else {
                None
            };
            let pd = &mut the_pd[a_face_iter];
            pd.set_h_nodes(nodes);
            pd.set_h_t_data(trian);
            pd.set_h_ph_dat(ph_dat);
            pd.set_face_index(a_face_iter as i32);
        }
    }
}
