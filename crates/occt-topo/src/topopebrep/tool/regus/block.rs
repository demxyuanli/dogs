//! `TopOpeBRepTool_REGUS` block traversal.
//!
//! Source (OCCT 8.0.0, `TopOpeBRepTool_REGUS.cxx`):
//! - `InitBlock()` (`cxx:670-706`)
//! - `NextinBlock()` (`cxx:709-851`)

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

use super::Regus;

impl Regus {
    /// `InitBlock()` (`cxx:670-706`): take the first face of the first edge
    /// that still has an untouched ancestor. Edges whose ancestor list ran dry
    /// are dropped from `mymapeFs` on the way.
    pub fn init_block(&mut self) -> bool {
        if self.my_eds_to_connect.len() != 0 {
            return false; // should be empty
        }

        for key in self.edge_keys() {
            let Some(lof) = self.my_mape_fs.get(&key).cloned() else {
                continue;
            };
            if lof.is_empty() {
                self.my_mape_fs.remove(&key);
                continue;
            }
            self.my_f = Some(lof[0].clone());
            return true;
        }
        false
    }

    /// `NextinBlock()` (`cxx:709-851`): the next face of the current block.
    ///
    /// The first loop walks the boundary edges already collected in
    /// `myedstoconnect`; when the current face shares none of them, OCCT falls
    /// back to the whole set (`cxx:726-850`) and picks the face that keeps the
    /// matter of the block continuous ([`Regus::nearest_f`]).
    pub fn next_in_block(&mut self) -> bool {
        // `cxx:716-725`: boundary edges of `myf`.
        let Some(my_f) = self.my_f.clone() else {
            return false;
        };
        let mut eds: Vec<TopoShape> = Vec::new();
        for e in explorer_shapes(&my_f, ShapeType::Edge) {
            if self.my_eds_to_connect.contains(&e) {
                eds.push(e);
            }
        }

        // `cxx:726-801`: `myf` shares no boundary edge.
        if eds.is_empty() {
            let candidates = self.my_eds_to_connect.items();
            for e in candidates {
                if !self.is_bound(&e) {
                    // all ancestor faces of <e> have been stored
                    self.my_eds_to_connect.remove(&e);
                    continue;
                }
                let Some(lof) = self.find(&e).cloned() else {
                    continue;
                };
                let nf = lof.len();
                if nf == 0 {
                    self.my_eds_to_connect.remove(&e);
                    self.unbind(&e);
                    continue;
                }

                if lof.len() == 1 {
                    self.my_f = Some(lof[0].clone());
                } else {
                    // `cxx:759-795`: a face of the block already connexed to
                    // `e`, then the nearest of the untouched ones.
                    let mapf: HashSet<ShapeMapKey> =
                        self.my_l_fin_block.iter().map(shape_map_key).collect();
                    let Some(lofc) = self.find_static(&e).cloned() else {
                        return false;
                    };
                    let fref = lofc
                        .iter()
                        .find(|fc| mapf.contains(&shape_map_key(fc)))
                        .cloned();
                    let Some(fref) = fref else {
                        return false; // `cxx:782-784`: !!!!!!!!!! a revoir 130499
                    };
                    self.my_f = Some(fref);
                    let Some(ffound) = self.nearest_f(&e, &lof) else {
                        return false;
                    };
                    self.my_f = Some(ffound);
                }
                return true;
            }
            return false;
        }

        // `cxx:803-848`: the first boundary edge of `myf` that leads somewhere.
        for e in eds {
            if !self.is_bound(&e) {
                self.my_eds_to_connect.remove(&e);
                continue;
            }
            let Some(lof) = self.find(&e).cloned() else {
                continue;
            };
            let nf = lof.len();
            if nf == 0 {
                self.my_eds_to_connect.remove(&e);
                self.unbind(&e);
                continue;
            }
            if nf == 1 {
                self.my_f = Some(lof[0].clone());
            } else {
                let Some(ffound) = self.nearest_f(&e, &lof) else {
                    return false;
                };
                self.my_f = Some(ffound);
            }
            return true;
        }
        false
    }
}
