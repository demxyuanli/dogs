//! `HLRBRep_Data::Write` (`HLRBRep_Data.cxx:548-608`): the splice that
//! merges one loaded shape data structure into an aggregate one.

use crate::hlrbrep_data::HlrBRepData;

impl HlrBRepData {
    /// `Write(DS, dv, de, df)` (`cxx:548-608`): appends the edges and faces of
    /// `ds` into this structure at the offsets `de` / `df`, shifting the vertex
    /// indices by `dv` and the wire edge indices by `de`.
    pub fn write(&mut self, ds: &HlrBRepData, dv: i32, de: i32, df: i32) {
        let n1edge = ds.nb_edges();
        let n1face = ds.nb_faces();

        // `HLRBRep_EdgeData* ed = &(myEData.ChangeValue(de)); ed++;`
        // (`cxx:557-559`): the first written slot is `de + 1`.
        let mut ed_idx = de + 1;
        let mut e1_idx = 1; // `ChangeValue(0)` then `e1++` (`cxx:560-561`)

        for _iedge in 1..=n1edge {
            if ed_idx < 0 || ed_idx as usize >= self.my_e_data.len() {
                break;
            }
            if (e1_idx as usize) >= ds.my_e_data.len() {
                break;
            }
            // `*ed = *e1;` (`cxx:569`).
            self.my_e_data[ed_idx as usize] = ds.my_e_data[e1_idx as usize].clone();

            if dv != 0 {
                // `ed->VSta(ed->VSta() + dv); ed->VEnd(ed->VEnd() + dv);`
                // (`cxx:571-575`).
                let ed = &mut self.my_e_data[ed_idx as usize];
                let vsta = ed.v_sta() + dv;
                ed.set_v_sta(vsta);
                let vend = ed.v_end() + dv;
                ed.set_v_end(vend);
            }

            // `myEMap.Add(DS->EdgeMap().FindKey(iedge));` (`cxx:577`).
            if let Some(key) = ds.my_e_map.get((e1_idx as usize).saturating_sub(1)) {
                self.my_e_map.push(key.clone());
            }

            ed_idx += 1;
            e1_idx += 1;
        }

        let mut fd_idx = df + 1;
        let mut f1_idx = 1;

        for _iface in 1..=n1face {
            if fd_idx < 0 || fd_idx as usize >= self.my_f_data.len() {
                break;
            }
            if (f1_idx as usize) >= ds.my_f_data.len() {
                break;
            }
            // `*fd = *f1;` (`cxx:585`).
            self.my_f_data[fd_idx as usize] = ds.my_f_data[f1_idx as usize].clone();

            if de != 0 {
                // `cxx:588-601`: every wire's every edge index is shifted.
                let fd = &mut self.my_f_data[fd_idx as usize];
                if let Some(wires) = fd.wires_mut() {
                    let nw = wires.nb_wires();
                    for iw in 1..=nw {
                        if let Some(eb) = wires.wire_mut(iw) {
                            let ne = eb.nb_edges();
                            for ie in 1..=ne {
                                let e = eb.edge(ie) + de;
                                eb.set_edge(ie, e);
                            }
                        }
                    }
                }
            }

            // `myFMap.Add(DS->FaceMap().FindKey(iface));` (`cxx:603`).
            if let Some(key) = ds.my_f_map.get((f1_idx as usize).saturating_sub(1)) {
                self.my_f_map.push(key.clone());
            }

            fd_idx += 1;
            f1_idx += 1;
        }
    }
}
