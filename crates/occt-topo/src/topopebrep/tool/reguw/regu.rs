//! `TopOpeBRepTool_REGUW::REGU`.
//!
//! Sources: `TopOpeBRepTool_REGUW.cxx:927-1132` (`REGU(istep, Scur, Splits)`)
//! and `:1136-1197` (`REGU()`).

use crate::abs::Orientation;
use crate::shape::TopoShape;
use crate::topo_tools::shape_map_key;

use super::super::closed::INTERNAL;
use super::super::topology::fun_tool_make_wire;
use super::Reguw;

impl Reguw {
    /// `REGU(istep, Scur, Splits)` (`cxx:927-1132`): one regularization step.
    ///
    /// Walks the UV-connexed edges of the wire into blocks; a block is closed
    /// when the current 2-D point comes back to the starting one. Every closed
    /// block becomes a new wire in `Splits`.
    pub fn regu(&mut self, istep: i32, splits: &mut Vec<TopoShape>) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        self.i_step = istep;

        splits.clear();

        if self.my_s.is_none() {
            return false;
        }
        // updated by the callers of chkp2dFORinE (`cxx:977-983`).
        self.my_tol2d = 1.0e-5;
        self.is_init0 = true;
        let mut lo_e_cur: Vec<TopoShape> = Vec::new();
        let mut lo_w: Vec<TopoShape> = Vec::new();

        let mut nite = 0;
        // `cxx:987-995`: INTERNAL/EXTERNAL edges do not take part in the count.
        let n_e = self
            .my_corr_iso
            .eds()
            .iter()
            .filter(|e| {
                !matches!(
                    e.orientation(),
                    Orientation::Internal | Orientation::External
                )
            })
            .count() as i32;

        while nite <= n_e {
            let mut fini = nite == n_e;

            // Getting starts elements if not defined (`cxx:1003-1020`)
            if self.is_init0 && !fini {
                if !self.init_block() {
                    return false;
                }
                if let Some(myed) = self.my_ed.clone() {
                    lo_e_cur.push(myed);
                }
                nite += 1;
                self.is_init0 = false;
                continue;
            }

            // Checking the wire is not already closed (`cxx:1024-1027`)
            let wire_is_closed = self.my_p2d.is_equal(&self.my_p2d0, self.my_tol2d);
            if wire_is_closed {
                // Adding INTERNAL/EXTERNAL edges to wires (`cxx:1031-1062`).
                // `cxx:1033` iterates `loEcur` while appending to it, so the
                // appended edges are visited as well.
                let mut i = 0;
                while i < lo_e_cur.len() {
                    let e = lo_e_cur[i].clone();
                    for v in crate::topo_tools::explorer_shapes(&e, crate::abs::ShapeType::Vertex)
                    {
                        let Some(co) = self.my_map_v_eds.get_mut(&shape_map_key(&v)) else {
                            continue;
                        };
                        let Some(le) = co.change_item(INTERNAL) else {
                            continue;
                        };
                        let mut j = 0;
                        while j < le.len() {
                            let ee = le[j].clone();
                            if matches!(
                                ee.orientation(),
                                Orientation::Internal | Orientation::External
                            ) {
                                lo_e_cur.push(ee);
                                le.remove(j);
                            } else {
                                j += 1;
                            }
                        }
                    }
                    i += 1;
                }
                fini = nite == n_e;

                // if Scur = <currentW> gives only one new wire, and has no new
                // split edges, <currentW> is valid and unchanged
                // (`cxx:1065-1071`)
                let one_w_ok = fini && lo_w.is_empty() && !self.hasnew_splits;
                if one_w_ok {
                    return true;
                }

                let mut new_w: Option<TopoShape> = None;
                let wi_ok = fun_tool_make_wire(&lo_e_cur, &mut new_w);
                if wi_ok {
                    if let Some(w) = new_w {
                        lo_w.push(w);
                    }
                } else {
                    return false;
                }
                self.is_init0 = true;
                lo_e_cur.clear();

                if fini {
                    splits.extend(lo_w.iter().cloned());
                    return true;
                }
                continue;
            }

            // Iteration on UV connexed edges (`cxx:1106-1124`)
            if !self.next_in_block() {
                return false;
            }
            if let Some(myed) = self.my_ed.clone() {
                lo_e_cur.push(myed);
            }
            nite += 1;
        }
        true
    }

    /// `REGU()` (`cxx:1136-1197`): step 1 on `S()`, then step 2 on every wire
    /// step 1 produced.
    ///
    /// `myOwNw` is bound even when nothing is regularized: `S()` maps to an
    /// empty split list, which is how `GetSplits` reports "unchanged".
    pub fn regu_all(&mut self) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        let Some(s) = self.my_s.clone() else {
            return false;
        };

        let to_regu = !self.my_list_v_multiple.is_empty() || self.hasnew_splits;
        let mut splits: Vec<TopoShape> = Vec::new();
        if !to_regu {
            self.my_ow_nw.insert(shape_map_key(&s), splits);
            return true;
        }

        // iStep = 1
        let mut lo_s: Vec<TopoShape> = Vec::new();
        if !self.regu(1, &mut lo_s) {
            return false;
        }

        // iStep = 2
        if lo_s.is_empty() {
            lo_s.push(s.clone()); // no new shape
        }

        for scur in lo_s {
            self.init_step(&scur);
            let _ = self.map_s();
            let to_regu1 = !self.my_list_v_multiple.is_empty();
            if !to_regu1 {
                splits.push(scur.clone());
                continue;
            }

            let mut sp: Vec<TopoShape> = Vec::new();
            if !self.regu(2, &mut sp) {
                return false;
            }
            if sp.is_empty() {
                sp.push(scur.clone()); // no new shape
            }
            splits.extend(sp);
        }
        self.my_ow_nw.insert(shape_map_key(&s), splits);
        true
    }
}
