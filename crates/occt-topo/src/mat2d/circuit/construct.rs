//! `MAT2d_Circuit::ConstructCircuit` (`MAT2d_Circuit.cxx:360-480`),
//! `SubSequence` (`cxx:346-356`) and `SortRefToEqui` (`cxx:736-755`).
//!
//! The circuit starts on the items of the reference line up to the first
//! connexion, then follows the mini path: the portion between two consecutive
//! connexions of the same line, or the rest of the line plus its head when the
//! connexion jumps over the last item (`PassByLast`). `UpDateLink` records, for
//! every `(line, item)` of the figure, the items of the circuit it produced.

use super::Mat2dCircuit;
use crate::mat2d::bi_int::Mat2dBiInt;
use crate::mat2d::item::{Mat2dFigure, Mat2dItem};
use crate::mat2d::mini_path::Mat2dMiniPath;

impl Mat2dCircuit {
    /// `ConstructCircuit(FigItem, IndRefLine, Road)` (`cxx:360-480`).
    pub(crate) fn construct_circuit(
        &mut self,
        fig_item: &Mat2dFigure,
        ind_ref_line: i32,
        road: &Mat2dMiniPath,
    ) {
        let nb_connexions = road.path().len() as i32;
        let mut set_of_item: Vec<Mat2dItem> = Vec::new();

        // `cxx:373-380`: the circuit starts on the head of the reference line.
        let mut prev_c = road.path()[0].clone();
        let prev_first = prev_c.borrow().index_item_on_first();
        sub_sequence(
            &fig_item[(ind_ref_line - 1) as usize],
            1,
            prev_first,
            &mut self.geom_elements,
        );
        self.update_link(1, ind_ref_line, 1, prev_first);
        let mut i_last_item = self.geom_elements.len() as i32;
        self.bind_connexion(i_last_item + 1, prev_c.clone());

        // `cxx:385-434`: the portions of line delimited by two consecutive
        // connexions.
        for i in 2..=nb_connexions {
            let cur_c = road.path()[(i - 1) as usize].clone();
            let cur_line = cur_c.borrow().index_first_line();

            if self.pass_by_last(&prev_c.borrow(), &cur_c.borrow()) {
                // `cxx:388-414`: the portion goes through the last item.
                let ind_last = fig_item[(cur_line - 1) as usize].len() as i32;
                let prev_second = prev_c.borrow().index_item_on_second();
                sub_sequence(
                    &fig_item[(cur_line - 1) as usize],
                    prev_second,
                    ind_last,
                    &mut set_of_item,
                );
                self.update_link(i_last_item + 1, cur_line, prev_second, ind_last);
                self.geom_elements.extend_from_slice(&set_of_item);
                i_last_item = self.geom_elements.len() as i32;

                if fig_item[(cur_line - 1) as usize].len() > 1 {
                    let cur_first = cur_c.borrow().index_item_on_first();
                    sub_sequence(
                        &fig_item[(cur_line - 1) as usize],
                        1,
                        cur_first,
                        &mut set_of_item,
                    );
                    self.update_link(i_last_item + 1, cur_line, 1, cur_first);
                    self.geom_elements.extend_from_slice(&set_of_item);
                    i_last_item = self.geom_elements.len() as i32;
                }
                self.bind_connexion(i_last_item + 1, cur_c.clone());
            } else {
                // `cxx:415-432`: the portion stays inside the line.
                let prev_second = prev_c.borrow().index_item_on_second();
                let cur_first = cur_c.borrow().index_item_on_first();
                sub_sequence(
                    &fig_item[(cur_line - 1) as usize],
                    prev_second,
                    cur_first,
                    &mut set_of_item,
                );
                self.update_link(i_last_item + 1, cur_line, prev_second, cur_first);
                self.geom_elements.extend_from_slice(&set_of_item);
                i_last_item = self.geom_elements.len() as i32;
                self.bind_connexion(i_last_item + 1, cur_c.clone());
            }
            prev_c = cur_c;
        }

        // `cxx:436-452`: closure from the last connexion to the last item of the
        // reference line.
        let ind_last = fig_item[(ind_ref_line - 1) as usize].len() as i32;
        if ind_last == 1 {
            self.bind_connexion(1, prev_c.clone());
            self.unbind_connexion(i_last_item + 1);
        } else {
            let prev_second = prev_c.borrow().index_item_on_second();
            sub_sequence(
                &fig_item[(ind_ref_line - 1) as usize],
                prev_second,
                ind_last,
                &mut set_of_item,
            );
            self.update_link(i_last_item + 1, ind_ref_line, prev_second, ind_last);
            self.geom_elements.extend_from_slice(&set_of_item);
        }

        // `cxx:454-465`: the lists of the items reached through a connexion come
        // first.
        for key in self.link_ref_equi_keys() {
            if self.link_ref_equi_len(&key) > 1 {
                self.sort_ref_to_equi(&key);
            }
        }
    }

    /// `SortRefToEqui(BiRef)` (`cxx:736-755`): the items before the first
    /// connexion are moved to the end of the list.
    pub(crate) fn sort_ref_to_equi(&mut self, key: &Mat2dBiInt) {
        let Some(seq) = self.link_ref_equi.get(key) else {
            return;
        };
        let Some(first) = seq.iter().position(|&index| !self.connexion_on(index)) else {
            return;
        };
        // 1-based `i` of `cxx:742-748`.
        let i = first + 1;
        if i > 1 && i <= seq.len() {
            let seq = self.link_ref_equi_mut(key);
            let head: Vec<i32> = seq[..first].to_vec();
            let tail: Vec<i32> = seq[first..].to_vec();
            seq.clear();
            seq.extend_from_slice(&tail);
            seq.extend_from_slice(&head);
        }
    }
}

/// `SubSequence(S1, IF, IL, S2)` (`cxx:346-356`): `S2` gets the items `IF` to
/// `IL` of `S1`, `S2` being cleared first. `IL < IF` leaves `S2` empty.
///
/// `MAT2d_MiniPath::MinimumL1L2` leaves `IndexItemOnFirst` /
/// `IndexItemOnSecond` at their `0` initial value whenever `Extrema_ExtCC2d`
/// finds no extremum inside the trims (`MAT2d_MiniPath.cxx:301`, `:427-434`),
/// so the degenerate circuits built with `GeomAbs_Tangent` /
/// `GeomAbs_Intersection` do reach this call with `0`. OCCT's
/// `NCollection_Sequence::Value` raises `Standard_OutOfRange` there (verified
/// against the exception-enabled OCCT 7.9.3 build: every `Tangent` /
/// `Intersection` case of a multi-contour face answers
/// `Standard_OutOfRange NCollection_Sequence::Value` out of
/// `BRepMAT2d_BisectingLocus::Compute`, and `BRepFill_OffsetWire::Perform`
/// turns it into `myIsDone = false`), so the port raises too.
fn sub_sequence(s1: &[Mat2dItem], if_: i32, il: i32, s2: &mut Vec<Mat2dItem>) {
    s2.clear();
    let mut i = if_;
    while i <= il {
        assert!(
            i >= 1 && (i as usize) <= s1.len(),
            "NCollection_Sequence::Value: index {i} out of range 1..={}",
            s1.len()
        );
        s2.push(s1[(i - 1) as usize].clone());
        i += 1;
    }
}
