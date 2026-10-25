//! `TopOpeBRepBuild_Builder::GFillFacesWESK` — one of the three fill passes
//! of a same-domain face pair into a `WireEdgeSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ffwesk.cxx`
//! (OCCT 8.0.0, 259 lines):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GFillFacesWESK(LS1, LS2, G1, WES, Kfill)` | `47-258` |
//!
//! The same body runs twice: first with `LS1` against `LS2` under `G1`, then
//! with the two lists swapped under `G1.CopyPermuted()` (`cxx:178-180`). In
//! both halves a face whose geometric orientation differs from the reference
//! face's flips the `GTopo` (`cxx:137-149`, `:220-233`); the fill body then
//! depends on `Kfill`:
//!
//! - `1` — `GFillFaceWES` (`cxx:152-155` / `:237-240`),
//! - `2` — `GFillCurveTopologyWES` (`cxx:156-162` / `:241-244`),
//! - `3` — `GFillONPartsWES` (`cxx:163-172` / `:245-252`).
//!
//! UNPORTED: the `OCCT_DEBUG` block (`cxx:66-80`, `:92-107`, `:182-198`,
//! `:247-250`) and the `Modified by Sergey KHROMOV` comment blocks
//! (`cxx:118-135`, `:141-149`, `:211-233`) — both are dead code in OCCT.

use crate::abs::{Orientation, State};
use crate::shape::TopoShape;
use crate::topopebrep::ds::kind::Config;

use super::super::g_topo::GTopo;
use super::super::wire_edge_set::WireEdgeSet;
use super::Builder;

impl Builder {
    /// `GFillFacesWESK(LS1, LS2, G1, WES, Kfill)` (`cxx:47-258`).
    pub fn g_fill_faces_wes_k(
        &mut self,
        ls1: &[TopoShape],
        ls2: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
        kfill: i32,
    ) {
        // `cxx:55-58`
        if ls1.is_empty() {
            return;
        }
        let (tb1, tb2) = g1.states_on();

        // `cxx:62-69`
        let Some(f1) = ls1.first().cloned() else {
            return;
        };
        self.my_face_reference = Some(f1.clone());
        let (oref, _iref) = match self.my_data_structure.as_ref() {
            Some(hds) => {
                let iref = hds.same_domain_reference(&f1);
                let oref = hds
                    .shape(iref, true)
                    .map(|s| s.orientation())
                    .unwrap_or(Orientation::Forward);
                (oref, iref)
            }
            None => (Orientation::Forward, 0),
        };

        // `cxx:85-89`: half 1/2.
        let gistoreverse1 = g1.is_to_reverse1();
        let mut g = g1.clone();
        let tb = tb1;
        self.fill_faces_wes_k_half(ls1, ls2, &mut g, wes, kfill, tb, oref, gistoreverse1);

        // `cxx:178-181`: half 2/1.
        let mut g = g1.copy_permuted();
        let gistoreverse1 = g.is_to_reverse1();
        let tb = tb2;
        self.fill_faces_wes_k_half(ls2, ls1, &mut g, wes, kfill, tb, oref, gistoreverse1);
    }

    /// The `for` loop of `GFillFacesWESK` (`cxx:109-175`, `:200-256`): one
    /// half of the pair, with the two lists in the roles `(Lclass, Lother)`.
    #[allow(clippy::too_many_arguments)]
    fn fill_faces_wes_k_half(
        &mut self,
        lclass: &[TopoShape],
        lother: &[TopoShape],
        g: &mut GTopo,
        wes: &mut WireEdgeSet,
        kfill: i32,
        tb: State,
        oref: Orientation,
        gistoreverse1: bool,
    ) {
        for s in lclass {
            // `cxx:111-116`
            let tomerge = !self.is_merged(s, tb);
            if !tomerge {
                continue;
            }

            // `cxx:137-149`: `DIFFORIENTED` with the same orientation as the
            // reference, or `SAMEORIENTED` with a different one, flips it.
            let (os, conf) = match self.my_data_structure.as_ref() {
                Some(hds) => {
                    let is = hds.shape_index(s, true);
                    let os = hds
                        .shape(is, true)
                        .map(|x| x.orientation())
                        .unwrap_or(Orientation::Forward);
                    (os, hds.same_domain_orientation(s))
                }
                None => (Orientation::Forward, Config::UnshGeometry),
            };
            let b1 = conf == Config::DiffOriented && os == oref;
            let b2 = conf == Config::SameOriented && os != oref;
            if b1 || b2 {
                g.set_reverse(!gistoreverse1);
            }

            // `cxx:152-174`
            match kfill {
                1 => self.g_fill_face_wes(s, lother, g, wes),
                2 => self.g_fill_curve_topology_wes(s, g, wes),
                3 => self.g_fill_on_parts_wes(s, g, lother, wes),
                _ => {}
            }
        }
    }
}
