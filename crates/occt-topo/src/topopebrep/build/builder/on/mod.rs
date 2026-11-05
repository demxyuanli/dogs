//! `TopOpeBRepBuild_BuilderON` — the part of `TopOpeBRepBuild_Builder` that
//! turns the ON parts of a face into start elements of a `WireEdgeSet`.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepBuild/`):
//! - `TopOpeBRepBuild_on.cxx:28-49` — `Builder::GFillONPartsWES`, the wrapper
//!   that picks `BuilderON::Perform2d` when `GLOBAL_faces2d` is set and
//!   `BuilderON::Perform` otherwise.
//! - `TopOpeBRepBuild_BuilderON.cxx:158-188` — `Perform`.
//! - `TopOpeBRepBuild_BuilderON.cxx:190-249` — `GFillONCheckI`.
//! - `TopOpeBRepBuild_BuilderON.cxx:251-290` — `GFillONPartsWES1`.
//! - `TopOpeBRepBuild_BuilderON.cxx:554-3123` — `GFillONPartsWES2`.
//!
//! `BuilderON` is a short-lived helper object holding the builder pointer, the
//! face, the `GTopo`, the class list and the `WireEdgeSet`
//! (`BuilderON.hxx:69-75`). The port has no such object: the same state is
//! passed explicitly to the methods, and `myPB` is the receiving `Builder`.
//!
//! `myPLSclass` is stored by both `Perform` overloads (`cxx:167`, `:75` of
//! `BuilderON2d.cxx`) and never read again — it is a dead field of the C++
//! class. The port keeps the parameter (callers pass it) and ignores it, which
//! is what OCCT does.
//!
//! `Perform2d` (`BuilderON2d.cxx:66-119`) and `GFillONParts2dWES2`
//! (`BuilderON2d.cxx:121-323`) - the `GLOBAL_faces2d` branch of the wrapper -
//! live in `perform_2d.rs`, together with `GLOBAL_DS2d` itself. That branch is
//! reached: `Builder::GMergeFaces` (`GridFF.cxx:409-413`) sets
//! `GLOBAL_faces2d = true` around its `GFillFacesWESK` double pass, and
//! `Builder1::SplitEdge` (`Builder1_1.cxx:670-674`) does the same around
//! `SplitEdge`, so the flag is not a debug switch.
//!
//! `GFillONPartsWES2` is complete: every branch of the `yapXX` chain is in
//! `wes2_yap00.rs` / `wes2_yapc.rs` / `wes2_fcx.rs` / `wes2_yap2.rs` /
//! `wes2_yap6.rs` / `wes2_yap3.rs` / `wes2_yap5.rs` / `wes2_yap4.rs` /
//! `wes2_yap7.rs`.
//!
//! `GFillONPartsWES` is reached from `GFillFacesWESK` (`ffwesk.cxx:163-172`,
//! `:245-252` → `builder/ffwesk.rs`, called at `grid_ff.rs:90-91` and
//! `:366-373`) and from `GSplitFace` (`GridSS.cxx:1381-1382` →
//! `builder/grid_ss_split.rs:46`, called at `ffsfs.rs:256`). The remaining
//! break in the chain is upstream: `Builder1::GFillShellSFS`/`GFillSolidSFS`
//! (the callers of `GFillFaceSFS`) are not invoked by `Builder::MergeShapes`
//! yet, so the ON parts do not reach a built face in the current gates (see
//! `_tasks.md` T-271 (5)).

mod helpers;
mod perform_2d;
mod wes2;

pub(crate) use perform_2d::{init_global_ds2d, with_global_ds2d_mut};
mod wes2_fcx;
mod wes2_yap0;
mod wes2_yap2;
mod wes2_yap3;
mod wes2_yap4;
mod wes2_yap5;
mod wes2_yap6;
mod wes2_yap7;
mod wes2_yapc;

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::interference::Interference;
use crate::topopebrep::ds::kind::Kind;

use super::super::g_topo::GTopo;
use super::super::wire_edge_set::WireEdgeSet;
use super::Builder;

impl Builder {
    /// `GFillONPartsWES(FOR, G, LSclass, WES)` (`on.cxx:28-49`).
    ///
    /// `BuilderON BON; if (GLOBAL_faces2d) BON.Perform2d(...); else
    /// BON.Perform(...)` (`cxx:33-45`): the 2D branch is
    /// [`Builder::perform_2d`] (`on/perform_2d.rs`), the other one is the walk
    /// below.
    pub fn g_fill_on_parts_wes(
        &self,
        for_: &TopoShape,
        g: &GTopo,
        _ls_class: &[TopoShape],
        wes: &mut WireEdgeSet,
    ) {
        // `on.cxx:34-42`: the `GLOBAL_faces2d` dispatch.
        if super::grid_ee::global_faces2d() {
            self.perform_2d(for_, g, wes);
            return;
        }
        // `Perform` (`BuilderON.cxx:158-188`).
        let Some(hds) = self.data_structure() else {
            return;
        };
        let li = hds.ds().shape_interferences(for_, true).to_vec();
        for i in li {
            let guard = i.read().expect("interference lock");
            self.g_fill_on_parts_wes1(for_, g, &guard, wes);
        }
    }

    /// `GFillONCheckI(I)` (`cxx:190-249`): `true` when `I` is a section edge's
    /// `EDGE`/`FACE` shape-shape interference whose transition is
    /// `(FACE, FACE)`, the edge is split ON with a non-empty split list, and
    /// both faces have a rank.
    pub fn g_fill_on_check_i(&self, for_: &TopoShape, i: &Interference) -> bool {
        let Some(hds) = self.data_structure() else {
            return false;
        };
        let bds = hds.ds();
        if !i.is_shape_shape() {
            return false;
        }
        let (gt, gi, st, si) = fds_data(i);
        if gt != Kind::Edge || st != Kind::Face {
            return false;
        }
        let Some(eg) = bds.shape(gi, false) else {
            return false;
        };
        let Some(fs) = bds.shape(si, false) else {
            return false;
        };

        let tfe = i.transition();
        let shab = tfe.shape_before();
        let shaa = tfe.shape_after();
        if shaa != ShapeType::Face || shab != ShapeType::Face {
            return false;
        }

        let is_rest = bds.is_section_edge(&eg, false);
        let is_split = self.is_split(&eg, State::On);
        if !(is_rest && is_split) {
            return false;
        }
        if self.splits(&eg, State::On).is_empty() {
            return false;
        }

        let rank_fs = self.shape_rank(&fs);
        let rank_for = self.shape_rank(for_);
        rank_fs != 0 && rank_for != 0
    }

    /// `GFillONPartsWES1(I)` (`cxx:251-290`): the guard of
    /// [`Builder::g_fill_on_check_i`], then every `Split(EG, ON)` piece is fed
    /// to `GFillONPartsWES2`.
    pub fn g_fill_on_parts_wes1(
        &self,
        for_: &TopoShape,
        g: &GTopo,
        i: &Interference,
        wes: &mut WireEdgeSet,
    ) {
        // `cxx:260`: `FDS_data(I, GT, GI, ST, SI)`.
        let (_gt, gi, _st, _si) = fds_data(i);

        if !self.g_fill_on_check_i(for_, i) {
            return;
        }

        let Some(hds) = self.data_structure() else {
            return;
        };
        let Some(eg) = hds.ds().shape(gi, true) else {
            return;
        };
        let l_esp_on = self.splits(&eg, State::On).to_vec();
        for esp_on in l_esp_on {
            self.g_fill_on_parts_wes2(for_, g, i, &esp_on, wes);
        }
    }
}
