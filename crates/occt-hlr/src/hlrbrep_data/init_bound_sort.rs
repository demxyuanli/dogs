//! `HLRBRep_Data::InitBoundSort` (`HLRBRep_Data.cxx:984-1019`): the edge
//! sort used by the hiding pass.

use crate::hlr_algo::MinMaxIndices;
use crate::hlrbrep_data::{HlrBRepData, MINMAX_SIGN_MASK};

impl HlrBRepData {
    /// `InitBoundSort(MinMaxTot, e1, e2)` (`cxx:984-1019`): collects the
    /// selected-range edges still not fully hidden whose 16-bit min/max box
    /// overlaps `MinMaxTot`, in ascending edge order.
    pub fn init_bound_sort(&mut self, min_max_tot: &MinMaxIndices, e1: i32, e2: i32) {
        self.my_nbr_sort_ed = 0;
        let min_max_shap = min_max_tot;

        for e in e1..=e2 {
            let ed = match self.my_e_data.get(e as usize) {
                Some(ed) => ed,
                None => continue,
            };
            if ed.status().all_hidden() {
                continue;
            }
            let lem = ed.min_max();
            // `cxx:996-1011`: the eight `& 0x80008000` pair tests, both ways.
            let overlap = (0..8).all(|k| {
                (min_max_shap.max[k] - lem.min[k]) & MINMAX_SIGN_MASK == 0
                    && (lem.max[k] - min_max_shap.min[k]) & MINMAX_SIGN_MASK == 0
            });
            if overlap {
                // `myNbrSortEd++; myEdgeIndices(myNbrSortEd) = e;` (`cxx:1014-1015`).
                self.my_nbr_sort_ed += 1;
                let idx = self.my_nbr_sort_ed as usize;
                if idx < self.my_edge_indices.len() {
                    self.my_edge_indices[idx] = e;
                }
            }
        }
    }

}
