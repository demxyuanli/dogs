//! `TopOpeBRepTool_REGUS::REGU` (`REGUS.cxx:520-666`) and the `FUN_update`
//! helper (`cxx:473-517`).
//!
//! `REGU` partitions `myS` into blocks: a block is a set of faces whose every
//! boundary edge has at most two ancestor faces *inside the block*. Each block
//! becomes a new shell, so a solid whose shell was glued through multi-ancestor
//! edges comes back as several closed shells.

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key};
use crate::topopebrep::tool::closed::is_closing_e;
use crate::topopebrep::tool::mk_shell::{mk_shell, remove_from_list};

use super::{EdgeSet, Regus};

/// `FUN_update(fcur, edstoconnect)` (`cxx:473-517`): toggle the block-boundary
/// edges of `f_cur` in `eds`.
///
/// An `INTERNAL`/`EXTERNAL` face (`cxx:486-489`) and an `INTERNAL`/`EXTERNAL`
/// or closing edge (`cxx:494-506`) are skipped; otherwise the edge is removed
/// when already bound (it now has two faces in the current block) and added
/// when not.
fn fun_update(f_cur: &TopoShape, eds: &mut EdgeSet) {
    let of_cur = f_cur.orientation();
    if of_cur == Orientation::Internal || of_cur == Orientation::External {
        return;
    }
    for e in explorer_shapes(f_cur, ShapeType::Edge) {
        let oe = e.orientation();
        if oe == Orientation::Internal || oe == Orientation::External {
            continue;
        }
        if is_closing_e(&e, f_cur) {
            continue;
        }
        if eds.contains(&e) {
            eds.remove(&e);
        } else {
            eds.add(&e);
        }
    }
}

impl Regus {
    /// `REGU()` (`cxx:520-666`): build one shell per block and bind it in
    /// `myOshNsh`.
    ///
    /// `cxx:529-533` answers `false` when nothing needs regularizing; the same
    /// `false` comes out of the `unchanged` arm (`cxx:579-582`), which leaves
    /// `myOshNsh` untouched — that is how `RegularizeShells` reports "nothing
    /// done" for a shell that was already valid.
    pub fn regu(&mut self) -> bool {
        let mut splits: Vec<TopoShape> = Vec::new();
        let to_regu = !self.my_mape_mult.is_empty() || self.my_n_f != self.my_old_n_f;
        if !to_regu {
            return false;
        }

        self.my_l_fin_block.clear();
        let mut nite = 0i32;
        while nite <= self.my_n_f {
            let mut start_block = self.my_l_fin_block.is_empty();
            let end_block = self.my_eds_to_connect.is_empty() && !start_block;

            // `cxx:572-602`: close the block.
            if end_block {
                let n_f_cur = self.my_l_fin_block.len() as i32;
                let unchanged = n_f_cur == self.my_old_n_f && self.my_n_f == self.my_old_n_f;
                if unchanged {
                    return false; // `cxx:582`: nyi analysis if we should raise
                }
                splits.push(mk_shell(&self.my_l_fin_block));
                self.my_l_fin_block.clear();
                start_block = true;
            }

            // `cxx:607-611`.
            if nite == self.my_n_f {
                break;
            }

            // `cxx:613-630`: first face of a block, or next one.
            let advance = if start_block || end_block {
                if !self.init_block() {
                    return false;
                }
                true
            } else {
                self.next_in_block()
            };

            // `cxx:632-633`.
            let Some(my_f) = self.my_f.clone() else {
                return false;
            };
            fun_update(&my_f, &mut self.my_eds_to_connect);

            // `cxx:635-646`.
            if !advance {
                let end_block = self.my_eds_to_connect.is_empty() && !start_block;
                if !end_block {
                    return false;
                }
                continue;
            }

            // `cxx:648-659`: the chosen face leaves `mymapeFs`.
            for e in explorer_shapes(&my_f, ShapeType::Edge) {
                if !self.is_bound(&e) {
                    continue; // ancestors faces of <e> are stored in Blocks
                }
                if let Some(lof) = self.change_find(&e) {
                    remove_from_list(lof, &my_f);
                }
            }

            self.my_l_fin_block.push(my_f);
            nite += 1;
        }

        self.my_osh_nsh
            .insert(shape_map_key(&self.my_s.clone().unwrap()), splits);
        true
    }
}
