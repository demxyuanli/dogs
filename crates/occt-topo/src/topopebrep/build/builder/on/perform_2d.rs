//! `TopOpeBRepBuild_BuilderON::Perform2d` / `GFillONParts2dWES2` - the
//! `GLOBAL_faces2d` branch of `GFillONPartsWES`.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepBuild/`):
//! - `TopOpeBRepBuild_BuilderON2d.cxx:66-119` - `Perform2d`, the walk of
//!   `GLOBAL_DS2d`'s interferences of the face `FOR`,
//! - `TopOpeBRepBuild_BuilderON2d.cxx:121-323` - `GFillONParts2dWES2`, the
//!   `yap6` decision that adds one `EspON` (an ON part of the edge `EG`) as a
//!   start element of the `WireEdgeSet`,
//! - `TopOpeBRepBuild_GridFF.cxx:331-410` - `GLOBAL_DS2d` and its fill, which
//!   [`crate::topopebrep::build::builder::grid_ff::Builder::g_merge_faces`]
//!   runs before it sets `GLOBAL_faces2d = true` (`cxx:409-413`) around the
//!   `GFillFacesWESK` double pass.
//!
//! `GLOBAL_DS2d` is a process-wide `TopOpeBRepDS_PDataStructure` created on
//! first use (`BuilderON2d.cxx:79-82`, `GridFF.cxx:343-347`) and re-`Init`ed by
//! every `GMergeFaces`. The port keeps it in a thread-local cell under the same
//! rules it already applies to `GLOBAL_faces2d` / `GLOBAL_issp`
//! (`build/builder/grid_ee.rs`): single-threaded visibility, no cross-thread
//! leakage.
//!
//! The interferences read here are *2D* interferences, i.e. the ones
//! `FUN_computeLIFfaces2d` (`GridFF.cxx:144-327`) writes into that separate
//! structure; they are not the interferences of the main DS that `Perform`
//! (`BuilderON.cxx:158-188`) walks.
//!
//! Divergences that are not heuristics:
//! - `TopoDS::Edge` / `TopoDS::Face` down-casts of a DS shape that answers
//!   another type raise `Standard_TypeMismatch` in OCCT; the port skips the
//!   entry, which is what the surrounding `FindShape` lookups already do for a
//!   missing shape,
//! - `SSI->GBound()` is `Standard_Boolean` with `true` as its constructor
//!   default (`TopOpeBRepDS_ShapeShapeInterference.cxx:26`); the port's
//!   `Interference::g_bound()` is an `Option<bool>` and the `None` arm answers
//!   `true`,
//! - `FUN_tool_orientEinFFORWARD` / `FUN_tool_orientEinF` answer `bool` plus an
//!   out-parameter in OCCT; the port answers `Option<Orientation>` and the
//!   unreached arm keeps `FORWARD`, which is the initial value of both OCCT
//!   locals,
//! - the `OCCT_DEBUG` blocks (`cxx:105-113`, `:138-166`, `:168-176`, `:179-190`,
//!   `:198-200`, `:214-218`, `:229-231`, `:252-255`, `:316-319`) are trace and
//!   dump only.

use std::cell::RefCell;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::complete_ds::fun_ds_share_g;
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::interference::{Interference, InterferenceHandle};
use crate::topopebrep::tool::{
    curves_same_oriented_point, fun_tool_bounds, fun_tool_orient_ein_f,
    fun_tool_orient_ein_f_forward,
};

use super::helpers::fun_build_tb;
use super::Builder;

thread_local! {
    /// `TopOpeBRepDS_PDataStructure GLOBAL_DS2d` (`GridFF.cxx:331`,
    /// `BuilderON2d.cxx:62`). `None` until the first `GMergeFaces`.
    static GLOBAL_DS2D: RefCell<Option<DataStructure>> = const { RefCell::new(None) };
}

/// The `GLOBAL_DS2d = new TopOpeBRepDS_DataStructure()` of
/// `GridFF.cxx:343-347` followed by `GLOBAL_DS2d->Init()` (`cxx:348`): a fresh
/// structure, which is what `Init()` leaves behind.
pub(crate) fn init_global_ds2d() {
    GLOBAL_DS2D.with(|cell| {
        *cell.borrow_mut() = Some(DataStructure::new());
    });
}

/// `GLOBAL_DS2d` as a mutable structure, or `None` when no `GMergeFaces` has
/// created it yet.
pub(crate) fn with_global_ds2d_mut<R>(f: impl FnOnce(&mut DataStructure) -> R) -> Option<R> {
    GLOBAL_DS2D.with(|cell| cell.borrow_mut().as_mut().map(f))
}

/// `GLOBAL_DS2d->ShapeInterferences(FOR)` (`BuilderON2d.cxx:84-85`).
fn global_ds2d_interferences(for_: &TopoShape) -> Vec<InterferenceHandle> {
    GLOBAL_DS2D.with(|cell| {
        cell.borrow()
            .as_ref()
            .map(|ds| ds.shape_interferences(for_, true).to_vec())
            .unwrap_or_default()
    })
}

impl Builder {
    /// `GFillONPartsWES(FOR, G, LSclass, WES)` (`TopOpeBRepBuild_on.cxx:28-49`)
    /// with `GLOBAL_faces2d` set: `BuilderON::Perform2d`.
    pub(super) fn perform_2d(&self, for_: &TopoShape, g: &GTopo, wes: &mut WireEdgeSet) {
        // `cxx:80-85`: the 2D structure is created on first use; `GMergeFaces`
        // is the only caller that creates it, and it does so before setting
        // `GLOBAL_faces2d`.
        let Some(hds) = self.data_structure() else {
            return;
        };
        let bds = hds.ds();
        // `cxx:73-77`: `myPB` / `myFace` / `myPG` / `myPLSclass` / `myPWES`
        // are the `BuilderON` object's fields; the port passes them.
        for handle in global_ds2d_interferences(for_) {
            let guard = handle.read().expect("interference lock");
            let (_gt, gi, _st, _si) = fds_data(&guard);
            // `cxx:105`: `BDS.Shape(GI)` down-cast to an edge.
            let Some(eg) = bds.shape(gi, false) else {
                continue;
            };
            if eg.shape_type() != ShapeType::Edge {
                continue;
            }
            // `cxx:108`: `myPB->Splits(EG, TopAbs_ON)`.
            let l_esp_on = self.splits(&eg, State::On).to_vec();
            for esp_on in l_esp_on {
                self.g_fill_on_parts_2d_wes2(&guard, &esp_on, for_, g, wes);
            }
        }
    }

    /// `GFillONParts2dWES2(I, EspON)` (`BuilderON2d.cxx:121-323`).
    fn g_fill_on_parts_2d_wes2(
        &self,
        i: &Interference,
        esp_on: &TopoShape,
        for_: &TopoShape,
        g: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let bds = hds.ds();
        // `cxx:129`: `myPG->StatesON(TB1, TB2)` is read and never used - the
        // state that matters below comes from `FUN_build_TB` (`cxx:156`).
        let _states_on = g.states_on();
        let (_gt, gi, _st, si) = fds_data(i);
        // `cxx:133-136`.
        let egbound_for = i.g_bound().unwrap_or(true);
        let ifor = bds.shape_index(for_, false);
        let Some(eg) = bds.shape(gi, false) else {
            return;
        };
        let Some(fs) = bds.shape(si, false) else {
            return;
        };
        // `cxx:152`: `rankEG = myPB->GShapeRank(EG)`.
        let rank_eg = self.shape_rank(&eg);
        // `cxx:156`: `TBEG = FUN_build_TB(myPB, rankEG)`.
        let tb_eg = fun_build_tb(self, rank_eg);
        // `cxx:164`: `eghassd = HDS->HasSameDomain(EG)`.
        let eg_has_sd = hds.has_same_domain(&eg, false);

        // `cxx:177-182`: `iFCX = SI; FFinSDSO = true; FFinSD = true; FCX = FS`.
        let ifcx = si;
        let fcx = fs;
        let ffin_sd = true;
        let ffin_sdso = true;

        // `cxx:194-201`: the orientation of EG inside FCX, and the second query
        // whose answer nothing reads.
        let oeg_fcx_f = fun_tool_orient_ein_f_forward(&eg, &fcx);
        let egbound_fcx = oeg_fcx_f.is_some();
        let oeg_fcx_f = oeg_fcx_f.unwrap_or(Orientation::Forward);
        let _ = fun_tool_orient_ein_f(&eg, &fcx);

        // `cxx:203-205`: the operation flags.
        let ope_fus = self.opefus();
        let ope_cut = self.opec12() || self.opec21();
        let ope_com = self.opecom();

        // `cxx:207-212`: `yap6`.
        let yap6 = ffin_sd && (egbound_for || egbound_fcx) && eg_has_sd;
        if !yap6 {
            return;
        }

        // `cxx:220`: `neworiE = TopAbs_FORWARD`. OCCT writes it before the
        // `rk1` test and both arms of `cxx:279-289` overwrite it, so the port
        // declares it where it is first given a live value.

        // `cxx:233-237`: only a rank-1 edge is kept.
        if rank_eg != 1 {
            return;
        }

        // `cxx:220` / `cxx:279-289`: the new orientation of the ON part.

        // `cxx:240-251`: `shareG` from the face that owns EG.
        let share_g = if egbound_fcx {
            fun_ds_share_g(hds, ifor, ifcx, gi, esp_on)
        } else if egbound_for {
            fun_ds_share_g(hds, ifcx, ifor, gi, esp_on)
        } else {
            None
        };
        let Some(share_g) = share_g else {
            return;
        };

        // `cxx:258-275`: keep the start element for this operation and state.
        let mut b = false;
        if ffin_sdso {
            if ope_fus {
                b = share_g;
            } else if ope_cut {
                b = !share_g && tb_eg == State::Out;
            } else if ope_com {
                b = share_g;
            }
        }
        if !b {
            return;
        }

        // `cxx:279-289`: the new orientation of the ON part.
        let mut new_ori_e;
        if !egbound_for {
            // `cxx:281-283` queries the orientation again and assigns the value
            // of the first query; the second result is dead.
            let _ = fun_tool_orient_ein_f_forward(&eg, &fcx);
            new_ori_e = oeg_fcx_f;
        } else {
            let oeg_for = fun_tool_orient_ein_f_forward(&eg, for_).unwrap_or(Orientation::Forward);
            new_ori_e = oeg_for;
        }

        // `cxx:291-303`: does `EspON` run the same way as `EG`?
        let (f, l) = fun_tool_bounds(esp_on);
        let par_on = (f + l) * 0.4352;
        let Some(eso) = curves_same_oriented_point(esp_on, par_on, &eg) else {
            return;
        };
        if !eso {
            new_ori_e = new_ori_e.complement();
        }

        // `cxx:305-311`: the two face orientations decide the rest.
        let o_for = bds
            .shape(ifor, false)
            .map(|s| s.orientation())
            .unwrap_or(Orientation::Forward);
        let o_fcx = bds
            .shape(ifcx, false)
            .map(|s| s.orientation())
            .unwrap_or(Orientation::Forward);
        if o_for != o_fcx {
            new_ori_e = new_ori_e.complement();
        }

        // `cxx:313-318`: `myPWES->AddStartElement(newE)`.
        let mut new_e = esp_on.clone();
        new_e.set_orientation(new_ori_e);
        wes.add_start_element(&new_e);
    }
}
