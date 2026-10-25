//! `TopOpeBRepBuild_Builder::End` block 3 and `UpdateSplitAndMerged`.
//!
//! Sources, OCCT 8.0.0:
//! * `TopOpeBRepBuild_End.cxx:366-410` — the `TopOpeBRepTool_FuseEdges` block
//!   that fuses the edges the boolean operation created inside every merged
//!   shape of `myShape1`.
//! * `TopOpeBRepBuild_End.cxx:413-520` — `UpdateSplitAndMerged`.
//!
//! Control flow of the block (`cxx:368-409`):
//!
//! 1. `makeFE` (`cxx:368-371`). `TopOpeBRepBuild_GetcontextNOFE` only exists
//!    under `OCCT_DEBUG`, so the block always runs in a release build.
//! 2. `ls = ChangeMerged(myShape1, myState1)`; each `SFE` is replaced in place
//!    by its fused version (`cxx:375-396`).
//! 3. The edge sets of **both** input shapes go into `mapOldEdges`
//!    (`cxx:380-387`) and become `FE.AvoidEdges`: the fusion may only touch the
//!    edges created by the operation, never the original ones.
//! 4. `FE.Edges(mle)`; `nle == 0` skips the shape (`cxx:388-391`).
//! 5. `FE.Perform()`; `SFE = FE.Shape()`; `FE.ResultEdges(mre)`;
//!    `FE.Faces(mlf)` (`cxx:393-400`).
//! 6. `UpdateSplitAndMerged(mle, mre, mlf, state)` for `IN`, `OUT`, `ON`
//!    (`cxx:404-406`).
//!
//! `UpdateSplitAndMerged` walks `MSplit(state)` (`cxx:426-429`). For an edge
//! key `e` it takes `ChangeSplit(e, state)` (`cxx:438`) and, per split piece
//! `edgecur`, searches the packets of `mle` for one that `IsSame`s it
//! (`cxx:449-463`). On a hit (`cxx:465-482`) the piece is removed from the
//! split list and `ChangeMerged(e, stateMerged)` becomes the one-element list
//! `{ mre(iLst) }`, with `stateMerged` being `myState1` when `ShapeRank(e) == 1`
//! and `myState2` otherwise. `NCollection_List::Remove(iterator)` leaves the
//! iterator on the element that followed, so the loop moves on without
//! advancing. For a face key (`cxx:497-517`) a split face that `mlf` binds to a
//! rebuilt face is replaced in place: `InsertBefore(mlf(facecur), it)` then
//! `Remove(it)`, which is an in-place value substitution in the port.
//!
//! The port copies the keys of `MSplit(state)` out before mutating: OCCT
//! mutates `ChangeSplit`/`ChangeMerged` of the same maps while holding the
//! `MSplit` iterator, which the borrow checker does not allow. `e` is always a
//! key of the map being changed, so no entry is created; the iteration order is
//! `NCollection_DataMap`'s, reproduced by the `BTreeMap`-backed
//! [`DShapeOn1StateMap`](crate::topopebrep::ds::tool::DShapeOn1StateMap).
//!
//! `ls` itself is also moved out of `my_merged_*` for the duration of the loop
//! (`cxx:375` aliases it) and written back at the end.

use std::collections::BTreeMap;

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::{is_same_occurrence, map_shapes};
use crate::topopebrep::tool::fuse_edges::FuseEdges;

use super::Builder;

impl Builder {
    /// `End.cxx:366-410`: fuse the operation's edges in the merged results of
    /// `myShape1`.
    pub(super) fn end_fuse_edges(&mut self) {
        // `cxx:375`: `ChangeMerged(myShape1, myState1)`. `myShape1` is null only
        // when `End` was reached without a first shape, in which case `ls` is
        // empty in OCCT too.
        let Some(shape1) = self.my_shape1.clone() else {
            return;
        };
        let state1 = self.my_state1;
        let state2 = self.my_state2;

        // `cxx:381-387`: `mapOldEdges` = the edges of `myShape1` and `myShape2`.
        let mut old_edges: Vec<TopoShape> = Vec::new();
        for s in [self.my_shape1.clone(), self.my_shape2.clone()].into_iter().flatten() {
            if let Some(v) = map_shapes(&s, &[ShapeType::Edge]).remove(&ShapeType::Edge) {
                old_edges.extend(v);
            }
        }

        let merged = std::mem::take(self.change_merged(&shape1, state1));
        let mut result: Vec<TopoShape> = Vec::with_capacity(merged.len());

        for sfe_orig in merged {
            // `cxx:378`: `TopOpeBRepTool_FuseEdges FE(SFE)`.
            let mut fe = FuseEdges::new(&sfe_orig);
            // `cxx:380-387`
            fe.avoid_edges(old_edges.iter().cloned());

            // `cxx:388-392`
            let mle = fe.edges();
            if mle.is_empty() {
                result.push(sfe_orig);
                continue;
            }

            // `cxx:393-394`
            let sfe = match fe.shape() {
                Ok(s) => s,
                Err(_) => {
                    // `Perform` raises `Standard_ConstructionError`
                    // (`cxx:187-206` propagates it) which would abort `End`; the
                    // port keeps the unfused shape and drops the (partial) map.
                    result.push(sfe_orig);
                    continue;
                }
            };

            // `cxx:396-400`
            let mre: BTreeMap<i32, TopoShape> = match fe.result_edges() {
                Ok(m) => m,
                Err(_) => {
                    result.push(sfe);
                    continue;
                }
            };
            let mlf = match fe.faces() {
                Ok(m) => m,
                Err(_) => {
                    result.push(sfe);
                    continue;
                }
            };

            // `cxx:404-406`
            for state in [State::In, State::Out, State::On] {
                self.update_split_and_merged(&mle, &mre, &mlf, state, state1, state2);
            }

            // `cxx:395`: `SFE = FE.Shape()` writes the merged element in place.
            result.push(sfe);
        }

        *self.change_merged(&shape1, state1) = result;
    }

    /// `UpdateSplitAndMerged(mle, mre, mlf, state)` (`End.cxx:413-520`).
    fn update_split_and_merged(
        &mut self,
        mle: &[(i32, Vec<TopoShape>)],
        mre: &BTreeMap<i32, TopoShape>,
        mlf: &[(TopoShape, TopoShape)],
        state: State,
        state1: State,
        state2: State,
    ) {
        // `cxx:426-429`: the keys of `MSplit(state)`, copied for the borrow of
        // `ChangeSplit` below.
        let keys: Vec<_> = self.m_split(state).keys().cloned().collect();

        for key in keys {
            let Some(e) = self.shape_of_key(&key) else {
                continue;
            };

            // `cxx:430-493`: edges.
            if e.shape_type() == ShapeType::Edge {
                // `cxx:468-474`: `StateRank` does not change inside the loop.
                let state_merged = if self.shape_rank(&e) == 1 { state1 } else { state2 };
                // `cxx:479-482` writes `ChangeMerged(e, stateMerged)` on every hit;
                // consecutive writes to one key leave the last one, so the value
                // is recorded here and written once the split list is released.
                let mut last_packet: Option<i32> = None;
                {
                    let lst_split = self.change_split(&e, state);
                    let mut i = 0;
                    while i < lst_split.len() {
                        let edgecur = lst_split[i].clone();
                        let mut found = false;

                        // `cxx:449-487`: for each packet of edges to be fused.
                        for (i_lst, lmap_edg) in mle {
                            if found {
                                break;
                            }
                            for edgefuse in lmap_edg {
                                // `cxx:461-462`
                                if is_same_occurrence(&edgecur, edgefuse) {
                                    found = true;
                                    // `cxx:465`: `LstSplit.Remove(itSplitEdg)`;
                                    // the OCCT iterator then sits on the next
                                    // element, so no advance is made here either.
                                    lst_split.remove(i);
                                    last_packet = Some(*i_lst);
                                    break;
                                }
                            }
                        }

                        // `cxx:489-492`
                        if !found {
                            i += 1;
                        }
                    }
                }
                if let Some(i_lst) = last_packet {
                    let lst_merged: Vec<TopoShape> =
                        mre.get(&i_lst).cloned().into_iter().collect();
                    *self.change_merged(&e, state_merged) = lst_merged;
                }
            }
            // `cxx:495-518`: faces.
            else if e.shape_type() == ShapeType::Face {
                let lst_split = self.change_split(&e, state);
                let mut i = 0;
                while i < lst_split.len() {
                    let facecur = lst_split[i].clone();
                    // `cxx:509-515`: `mlf.IsBound(facecur)` uses
                    // `TopTools_ShapeMapHasher`, i.e. `IsSame`.
                    match mlf.iter().find(|(f, _)| is_same_occurrence(f, &facecur)) {
                        Some((_, rebuilt)) => {
                            // `InsertBefore(mlf(facecur), it)` + `Remove(it)` is
                            // an in-place substitution; the iterator lands on the
                            // following element either way.
                            lst_split[i] = rebuilt.clone();
                            i += 1;
                        }
                        None => i += 1,
                    }
                }
            }
        }
    }
}
