//! `TopOpeBRepDS_EIR` — the edge-interference reducer.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.hxx:30-46`,
//! `TopOpeBRepDS_EIR.cxx:1-1221` (OCCT 8.0.0).
//!
//! `ProcessEdgeInterferences()` walks every edge of the DS (`cxx:809-825`);
//! `ProcessEdgeInterferences(EIX)` groups the interferences of one edge by
//! geometry (`cxx:962-1160`) and hands each group to
//! [`fun_process_edge_interferences`], which is the reduction proper:
//!
//! ```text
//! LI -> (lF + lFE) + lE
//! lF  = {(T(face), G, S)}          lFE = {(T(face), G, S=edge)}
//! lE  = {(T(edge), G, S)}
//! lFE -> l3dFE + l2dFE + lFEresi   lF -> l3dF
//! lE  -> l1dE + l1dEsd
//! ```
//!
//! Translation notes:
//! - OCCT mutates `BDS.ChangeShapeInterferences(EIX)` in place through a
//!   reference while the helpers read `BDS.ShapeInterferences(EIX)` (the same
//!   list, at `FUN_select2dI`). The port therefore works on a snapshot of the
//!   handle list and writes the result back at the end of
//!   [`Eir::process_edge_interferences`], which is what the
//!   `LI.Clear(); LI.Append(LInew)` at `cxx:1049-1050` does.
//! - `performPNC` is the literal `false` at `cxx:1052`, so the block at
//!   `cxx:1060-1160` (the `FDSCNX_EdgeConnexitySameShape` point pruning) is
//!   unreachable and is left UNPORTED, as OCCT never runs it.
//! - The `FDS_repvg` / `FDS_reduceONFACEinterferences` / `FUN_reclSE` helpers
//!   are the ones at their own `.cxx` sites; see the submodules.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::Edge;

use super::data_structure::DataStructure;
use super::h_data_structure::HDataStructure;
use super::interference_iterator::InterferenceList;
use super::kind::Kind;
use super::process_interferences_tool::make_epv_interference_with_support;
use super::tki::TKI;

mod on_face;
mod process;
mod pure2d;
mod recl_se;
mod reducer_edge;
mod reducer_edge3d;

pub use on_face::fds_reduce_onface_interferences;
pub use process::{check_interference_is_valid, fun_process_edge_interferences};
pub use pure2d::{fun_keep_l3d_f, fun_reduce_pure2d_i, fun_reduce_pure2d_i0};
pub use recl_se::{fun_recl_se, fun_recl_se2};
pub use reducer_edge::{fun_reducer_edge, fun_reducer_sd_edge, fun_unkeep_evi_on_gb1};
pub use reducer_edge3d::fun_reducer_edge3d;

/// `TopOpeBRepDS_EIR` (`hxx:30-46`): the reducer, holding the HDS it rewrites.
pub struct Eir<'a> {
    /// `myHDS`
    hds: &'a mut HDataStructure,
}

impl<'a> Eir<'a> {
    /// `TopOpeBRepDS_EIR(HDS)` (`cxx:802-805`).
    pub fn new(hds: &'a mut HDataStructure) -> Self {
        Self { hds }
    }

    /// `ProcessEdgeInterferences()` (`cxx:809-825`): every edge shape of the
    /// DS, in index order.
    pub fn process_all_edge_interferences(&mut self) -> Result<(), String> {
        let nshape = self.hds.ds().nb_shapes();
        for i in 1..=nshape {
            let Some(s) = self.hds.ds().shape(i, true) else {
                continue;
            };
            if s.shape_type() == ShapeType::Edge {
                self.process_edge_interferences(i)?;
            }
        }
        Ok(())
    }

    /// `ProcessEdgeInterferences(EIX)` (`cxx:962-1160`).
    pub fn process_edge_interferences(&mut self, eix: i32) -> Result<(), String> {
        let bds: &mut DataStructure = self.hds.change_ds();

        // E is the edge, LI is list of interferences to compact
        let e = bds
            .shape(eix, true)
            .ok_or_else(|| "TopOpeBRepDS_EIR: BDS.Shape(EIX)".to_string())?;
        if BRepTool::is_degenerated(&Edge(e)) {
            return Ok(());
        }

        let snapshot = bds.shape_interferences_at(eix, true).to_vec();
        let mut newtki = TKI::new();
        newtki.fill_on_geometry(&snapshot);
        let mut tki = TKI::new();
        tki.fill_on_geometry(&snapshot);

        // `cxx:979-1037`: merge the interferences of a vertex with the ones of
        // its same-domain vertex, replacing `vG` by `OOvG`.
        let mut order: Vec<(Kind, i32)> = Vec::new();
        tki.init();
        while tki.more() {
            if let Some((k, g, _)) = tki.value() {
                order.push((k, g));
            }
            tki.next();
        }

        for (k, g) in order {
            let loi: InterferenceList = tki.interferences(k, g).to_vec();
            if k == Kind::Point {
                continue;
            }

            let vg = bds
                .shape(g, true)
                .ok_or_else(|| "TopOpeBRepDS_EIR: BDS.Shape(G)".to_string())?;
            let Some(oovg) = super::export::fun_ds_getoov(bds, &vg) else {
                continue;
            };
            let oog = bds.shape_index(&oovg, true);
            if oog == 0 {
                continue; // NYIRaise
            }

            let isb = newtki.is_bound(k, oog);

            // xpu201098 : cto904F6, e10,v6
            let isbound = loi
                .first()
                .map(|h| h.read().expect("interference lock").g_bound().unwrap_or(false))
                .unwrap_or(false);
            if isbound {
                // replacing vG with oovG
                let mut newloi: InterferenceList = Vec::new();
                for h in loi.iter() {
                    let (t, s, par, st) = {
                        let i = h.read().expect("interference lock");
                        let (_gt, _g1, st, s) = super::export::fds_data(&i);
                        let par = super::export::fds_parameter(&i)
                            .ok_or_else(|| "TopOpeBRepDS_EIR: FDS_Parameter".to_string())?;
                        (i.transition(), s, par, st)
                    };
                    match make_epv_interference_with_support(&t, s, oog, par, k, st, false) {
                        Some(newi) => newloi.push(newi),
                        // `MakeEPVInterference` leaves a null handle for any
                        // other geometry kind (`ProcessInterferencesTool.cxx:73-93`);
                        // a null handle cannot be carried by the port.
                        None => {}
                    }
                }
                if let Some(l) = newtki.change_interferences(k, g) {
                    l.clear();
                }
                if !isb {
                    newtki.add(k, oog)?;
                }
                if let Some(l) = newtki.change_interferences(k, oog) {
                    l.extend(newloi);
                }
                continue;
            }

            if !isb {
                continue;
            }
            let li: InterferenceList = newtki
                .change_interferences(k, oog)
                .map(|l| l.clone())
                .unwrap_or_default();
            if let Some(l) = newtki.change_interferences(k, g) {
                l.extend(li);
            }
        } // tki

        // `cxx:1039-1050`
        let mut li_new: InterferenceList = Vec::new();
        let mut order: Vec<(Kind, i32)> = Vec::new();
        newtki.init();
        while newtki.more() {
            if let Some((k, g, _)) = newtki.value() {
                order.push((k, g));
            }
            newtki.next();
        }
        for (k, g) in order {
            let mut loi: InterferenceList = newtki.interferences(k, g).to_vec();
            fun_process_edge_interferences(eix, k, g, bds, &mut loi)?;
            li_new.extend(loi);
        }
        let target = bds.change_shape_interferences_at(eix);
        target.clear();
        target.extend(li_new);

        // `cxx:1052-1057`: `performPNC` is `false` ("JYL 28/09/98 :
        // temporaire"), so the pruning at `cxx:1060-1160` never runs.
        Ok(())
    }
}


